/* Copyright (c) Colorado School of Mines, 1995.*/
/* All rights reserved.                       */

/* SUSORT: $Revision: 1.1 $ ; $Date: 1997/02/26 09:55:59 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>
#include <errno.h>

#define ISOFF 200

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPSORT - memory based sort on any segy header keywords                        ",
"                                                                               ",
" susort2 <stdin >stdout [[+-]key1 [+-]key2 ...]                                ",
"                                                                               ",
" This is an adapted version of susort that both input and output are allowed to",
" be via the pipe. If data is read via a pipe a temporary file is created on the",
" fly in system memory. This speeds up the sorting but limits the size of input ",
" file to be sorted. To prevent excessive memory consumption, the max. size of  ",
" input file is fixed to 512 MB.                                                ",
"                                                                               ",
" For more info see selfdoc of susort                                           ",
"                                                                               ",
" Modified by Sanyu Ye from susort2 of DELPHI                                   ",
"                                                                               ",
    NULL
};

/* Credits:
 *	SEP: Einar, Stew
 *	CWP: Shuki, Jack
 *	Delphi: Alexander
 *	RWS: Ye, Sanyu
 *
 */
/**************** end self doc ***********************************/


#define NTRSTEP	1024	/* realloc() increment measured in traces */
#define MAXSIZE	512	/* max. size in MB of input data via pipe  */

segy tr;
static int nkey; /* number of keys to sort on	*/
static cwp_String type; /* header key types		*/

/* Prototypes */
Value negval(cwp_String type, Value val); /* reverse sign of value	*/
int cmp_list(); /* qsort comparison function	*/
/* extract maximum traces per ensemble*/
extern int get_maxtr_ens(Value *, long, int);
extern void shutdown(void);

/* set global pointers */
FILE *datafp = NULL; /* fp for data storage file		*/
char nametmp[32]; /* filename of temporary file		*/

int main(int argc, char **argv) {
    static Value *val_list; /* a list of the key values for each    */
    /* trace with each group headed by the	*/
    /* trace number of that trace		*/
    static int *index; /* header key indices			*/
    static cwp_Bool *up; /* sort direction (+ = up = ascending)	*/
    register Value *vptr; /* location pointer for val_list	*/
    int ngroup; /* size of unit in val_list (nkey + 1)	*/
    int nv; /* number of groups in val_list		*/
    int nvsize; /* size of group in val_list		*/
    Value val; /* value of a keyword			*/
    long ntr; /* number of traces from gettr		*/
    long itr; /* index of trace (0,1,2,...)		*/
    FileType ftypein; /* filetype of stdin			*/
    int nsegy; /* number of bytes on the trace		*/
    cwp_Bool ispipe; /* ftypein == PIPE ?			*/
    int iskey; /* sortkey code				*/
    register int i;


    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* set signal handler */
    signal(SIGHUP, (void (*) (int)) shutdown);
    signal(SIGINT, (void (*) (int)) shutdown);
    signal(SIGQUIT, (void (*) (int)) shutdown);
    signal(SIGTERM, (void (*) (int)) shutdown);

    /* The algorithm requires stdin to be rewound and */
    /* random access to either stdin or stdout.    */
    ftypein = filestat(STDIN);
    ispipe = (ftypein == PIPE || ftypein == FIFO) ? cwp_true : cwp_false;

    /* If pipe, prepare temporary file to hold data */
    if (ispipe) {
        sprintf(nametmp, "/dev/shm/spsort%d", (int) getpid());
        datafp = fopen(nametmp, "wb+");
        if (datafp <= 0) {
            err("Failed to open memory file %s: %s", nametmp, strerror(errno));
        }
    } else datafp = stdin;

    /* Set number of sort keys */
    if (argc == 1) nkey = 1; /* no explicit keys: default cdp key */
    else nkey = argc - 1; /* one or more explicit keys */


    /* Allocate lists for key indices, sort directions and key types */
    index = ealloc1int(nkey);
    up = (cwp_Bool *) ealloc1(nkey, sizeof (cwp_Bool));
    type = (char *) ealloc1(nkey, sizeof (char));


    /* Initialize index, type and up */
    if (argc == 1) {
        index[0] = getindex("cdp");
        up[0] = cwp_true;
        type[0] = 'l';
    } else {
        register int i;
        for (i = 0; i < nkey; ++i) {
            switch (**++argv) { /* sign char of next arg */
                case '+':
                    up[i] = cwp_true;
                    ++*argv; /* discard sign char in arg */
                    break;
                case '-':
                    up[i] = cwp_false;
                    ++*argv; /* discard sign char in arg */
                    break;
                default:
                    up[i] = cwp_true;
                    break;
            }
            index[i] = getindex(*argv);
            type[i] = hdtype(*argv)[0]; /* need single char */
        }
    }
    iskey = index[0] + 1;
    iskey += nkey > 1 ? ISOFF * (1 + index[1]) : 0;

    /* Allocate list of trace numbers + key values */
    ngroup = nkey + 1;
    nvsize = ngroup * sizeof (Value);
    nv = NTRSTEP; /* guess at required number */
    val_list = (Value *) ealloc1(nv, nvsize);
    vptr = val_list;


    /* Run through traces once to collect header values */
    ntr = 0;
    if (!(nsegy = gettr(&tr))) err("can't get first trace");

    do {
        itr = ntr++;
        /* realloc if out of room */
        if (0 == (ntr % NTRSTEP)) {
            if (ntr*tr.ns > MAXSIZE*250000 && ispipe) {
                err("input file too big, use susort2 or spdbwrite/read");
            }
            nv += NTRSTEP;
            val_list = (Value *) erealloc(val_list, nv * nvsize);
            vptr = val_list + itr * ngroup;
        }

        /* enter trace index in list and then key values */
        vptr++->l = itr; /* insert number and advance */
        {
            register int i;
            for (i = 0; i < nkey; ++i) {
                gethval(&tr, index[i], &val);
                *vptr++ = up[i] ? val : negval(type + i, val);
                /* negative values give reverse sort */
            }
        }
        if (ispipe) fputtr(datafp, &tr);
    } while (gettr(&tr));

    rewind(datafp);

    /* Sort the headers values */
    qsort(val_list, ntr, nvsize, cmp_list);

    /* run through sorted list and write output sequentially */
    for (i = 0; i < ntr; ++i) {
        itr = val_list[i * ngroup].l;
        fgettra(datafp, &tr, itr);
        tr.tracl = i + 1;
        puttr(&tr);
    }

    shutdown();
    return EXIT_SUCCESS;
}

/* Comparison routine for qsort */
int cmp_list(register Value *a, register Value *b) {
    register int i;
    Value va, vb;
    int compare;

    /* Can order as soon as any components are unequal */
    for (i = 0; i < nkey; ++i) {
        va = *++a;
        vb = *++b; /* advance and dereference */
        if ((compare = valcmp(type + i, va, vb)))
            return compare;
    }
    return 0;
}

/* Reverse sign of value */
Value negval(cwp_String type, Value val) {
    switch (*type) {
        case 'h':
            val.h = -val.h;
            break;
        case 'u':
            val.u = (short) -val.u;
            break;
        case 'l':
            val.l = -val.l;
            break;
        case 'v':
            val.u = (long) -val.u;
            break;
        case 'i':
            val.i = -val.i;
            break;
        case 'p':
            val.p = (int) -val.p;
            break;
        case 'f':
            val.f = -val.f;
            break;
        case 'd':
            val.d = -val.d;
            break;
        default: err("%d: mysterious type %s", __LINE__, type);
    }

    return val;
}

void shutdown() {
    /* remove temporary file if it exists */
    if (datafp) {
#ifdef DEBUG
        fprintf(stderr, "Removing temporary file \n");
#endif
        efclose(datafp);
        unlink(nametmp);
    }
}
