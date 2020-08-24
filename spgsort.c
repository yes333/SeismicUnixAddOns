/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                      */

/* SUSETGATHER: $Revision: 1.0 $ ; $Date: 2008/02/07 $        */

#include <su.h>
#include <segy.h>
#include <header.h>
#include <segyhdr.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPGSORT- Sort within single gather upon header key value change               ",
"                                                                               ",
"   spgsort <infile >outfile  [optional parameters]                             ",
"                                                                               ",
" Optional Parameters:                                                          ",
"   key=cdpt        gather sorting key                                          ",
"                                                                               ",
"   sort=           header keywords for sorting                                 ",
"                   a minus sign (-) preceding keyword for descending order     ",
"                                                                               ",
"   nmax=1024       max. number of traces of gathers                            ",
"                                                                               ",
"   verbose=0       >0 echo information                                         ",
"                                                                               ",
" Examples:                                                                     ",
"   1. read 3C, do taup transform, sort 4C trace by trace for wavefield         ",
"      decomposition and sort back according to wave mode                       ",
"   spradon key=duse < PZX.su |\\",
"   spgsort key=cdpt sort=sx,duse nmax=4000 |\\",
"   spwfs nc=3 modes=-2,2 |\\",
"   spgsort key=cdpt sort=-duse,sx nmax=4000 |\\",
"   spmovie key=cdpt names=duse desc=down-up-P-on-Z ",
"                                                                               ",
" Version 1.0 last modified Dec. 2011 by Sanyu Ye                               ",
NULL};

/*
* Credits: READ Well Service, Sanyu Ye, Dec. 2011.
*
*/
/**************** end self doc ********************************/

static int nkey; /* number of keys to sort on	*/
static cwp_String type; /* header key types		*/

/* Prototypes */
Value negval(cwp_String type, Value val); /* reverse sign of value	*/
int cmp_list(); /* qsort comparison function	*/

/** main **/
int main(int argc, char **argv)
{
    cwp_String *key, gkey;    /* header key word from segy.h */
    cwp_String gtype;   /* type of key                 */
    int *index, gindex;         /* index of key                */
    Value val, valnew, kval; /* value of key                */
    int nsegy;          /* Bytes read in for a trace   */
    int ns;             /* number of time samples of trace		*/
    int verbose;        /* flag for echoing information */
    int ntr, ntotal;    /* trace counters */
    int ngather;        /* gather counters */
    int nmax;           /* max. number of traces expected in gather */
    int eof = 0;        /* END OF File flag			*/
    int i, itr;
    segy tr, outtr;
    segyhdr *hdrs = NULL; /* placeholder for incoming trace headers	*/
    float **trdata = NULL; /* temporary array for trace data		*/

    /* hook up getpar to handle the parameters */
    initargs(argc,argv);
    requestdoc(1);

    if (!getparint("verbose", &verbose)) verbose= 0;
    if (!getparint("nmax", &nmax)) nmax = 1024;

    /* Get sort keys  */
    if ((nkey=countparval("sort"))!=0) {
        key = (cwp_String *) ealloc1(nkey, sizeof (cwp_String));
        getparstringarray("sort", key);
    } else {
        err("keywords for sorting must be given (sort=...)");
    }

    /* Allocate lists for key indices, sort directions and key types */
    index = ealloc1int(nkey);
    cwp_Bool *up = (cwp_Bool *) ealloc1(nkey, sizeof (cwp_Bool));
    type = (char *) ealloc1(nkey, sizeof (char));


    /* Initialize index, type and up */
    for (i = 0; i < nkey; ++i) {
        switch (key[i][0]) { /* sign char of next arg */
            case '+':
                up[i] = cwp_true;
                ++key[i]; /* discard sign char in arg */
                break;
            case '-':
                up[i] = cwp_false;
                ++key[i]; /* discard sign char in arg */
                break;
            default:
                up[i] = cwp_true;
                break;
        }
        index[i] = getindex(key[i]);
        type[i] = hdtype(key[i])[0]; /* need single char */
    }

    /* Allocate list of trace numbers + key values */
    int ngroup = nkey + 1;
    int nvsize = ngroup * sizeof (Value);
    Value* val_list = (Value *) ealloc1(nmax, nvsize);
    Value* vptr = val_list;

    /* Get info from first trace */
    if (!(nsegy = gettr(&tr))) err("can't read first trace");
    ns = (int) tr.ns;

    /* get gather sorting key */
    if (!getparstring("key", &gkey)) gkey = "cdpt";
    gtype = hdtype(gkey);
    gindex = getindex(gkey);
    gethval(&tr, gindex, &val);

    // allocate and reset memory for input and output traces
    hdrs = (segyhdr*) ealloc1(nmax, HDRBYTES);
    trdata = ealloc2float(ns, nmax);
    memset(hdrs, 0, nmax * HDRBYTES);
    memset(*trdata, 0, nmax * ns * FSIZE);

    /* Read headers and data while getting a count */
    ngather = 0;
    ntr = 0;
    ntotal = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, gindex, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(gtype, val, valnew)) { /* same key and more data*/
            if (ntr > nmax - 1) err("\nNumber of traces exceeding nmax=%d\n", nmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);
            memcpy(trdata[ntr], tr.data, FSIZE * ns);

            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (verbose > 1) warn("  processing %d traces of %d-th gather with %s=%d",
                                     ntr, ngather, key, vtoi(type, val));

            ntotal += ntr;

            /* Run through traces once to collect header values */

            for (itr = 0; itr < ntr; ++itr) {
                /* enter trace index in list and then key values */
                vptr++->l = itr; /* insert number and advance */
                {
                    for (i = 0; i < nkey; ++i) {
                        gethval((segy*) &hdrs[itr], index[i], &kval);
                        *vptr++ = up[i] ? kval : negval(type + i, kval);
                        /* negative values give reverse sort */
                    }
                }
            }

            /* Sort the headers values */
            qsort(val_list, ntr, nvsize, cmp_list);

            /* run through sorted list and write output sequentially */
            for (i = 0; i < ntr; ++i) {
                itr = val_list[i * ngroup].l;
                memcpy(&outtr, &hdrs[itr], HDRBYTES);
                memcpy(outtr.data, trdata[itr], FSIZE * ns);
                outtr.tracl = i + 1;
                puttr(&outtr);
            }

            val = valnew;
            // reset output data
            memset(hdrs, 0, nmax * HDRBYTES);
            memset(*trdata, 0, nmax * ns * FSIZE);
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose > 0) warn("  Totally %d gathers with %d traces processed", ngather, ntotal);

    return(CWP_Exit());
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

