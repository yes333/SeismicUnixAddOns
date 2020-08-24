/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* SUSTACK: $Revision: 1.25 $ ; $Date: 2006/10/31 22:03:52 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "sprinthlp.c"

/*********************** self documentation **********************/
char *sdoc[] = {
    "                                                                           ",
    " SUSTACK - stack adjacent traces having the same key header word           ",
    "                                                                           ",
    "     sustack <stdin >stdout [Optional parameters]                          ",
    "                                                                           ",
    " Required parameters:                                                      ",
    " 	none                                                                    ",
    "                                                                           ",
    " Optional parameters:                                                      ",
    " 	key=cdp		header key word to stack on                             ",
    " 	normpow=1.0	each sample is divided by the                           ",
    "			normpow'th number of non-zero values                    ",
    "			stacked (normpow=0 selects no division)                 ",
    "	average=key1,...averages values of specified keys, no default           ",
    "	fold=counit	stores fold in keyword, default counit                  ",
    "	repeat=0	=1 repeats the stack trace nrepeat times                ",
    "	nrepeat=10	repeats stack trace nrepeat times in                    ",
    "	          	output file                                             ",
    " 	verbose=0	verbose = 1 echos information                           ",
    "                                                                           ",
    " Notes:                                                                    ",
    " ------                                                                    ",
    " Header values of first input trace of each gather are overtaken           ",
    " for stacked trace output. Average header values can be computed           ",
    " if the corresponding keywords are specified with option                   ",
    "                                                                           ",
    "     average=key1,key2,...                                                 ",
    " 							        ",
    " Sushw can be used afterwards if this is not acceptable.	",
    "								",
    " For VSP users:						",
    " The stack trace appears ten times in the output file when	",
    " setting repeat=1 and nrepeat=10. Corridor stacking can be	",
    " achieved by properly muting the upgoing data with SUMUTE	",
    " before stacking.						",
    "								",
    NULL
};

/* Credits:
 *	SEP: Einar Kjartansson
 *	CWP: Jack K. Cohen, Dave Hale
 *	CENPET: Werner M. Heigl - added repeat trace functionality
 *	RWS: Sanyu Ye - added average key values functionality
 *
 * Note:
 *	The "valxxx" subroutines are in su/lib/valpkge.c.  In particular,
 *      "valcmp" shares the annoying attribute of "strcmp" that
 *		if (valcmp(type, val, valnew) {
 *			...
 *		}
 *	will be performed when val and valnew are different.
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: nhs, tracl, offset
 */
/**************** end self doc ***********************************/


int main(int argc, char **argv) {
    cwp_String key; /* header key word from segy.h		*/
    cwp_String type; /* ... its type				*/
    int indx; /* ... its index			*/

    cwp_String keyf; /* header key word with fold value		*/
    cwp_String typef; /* ... its type				*/
    int indxf; /* ... its index			*/
    Value valf; /* its value */

    cwp_String keya[SU_NKEYS]; /* array of keywords to average		*/
    cwp_String typea[SU_NKEYS]; /* array of name of type keywords to average */
    int indxa[SU_NKEYS]; /* index	of getparred key to match	*/
    Value vala[SU_NKEYS]; /* values of matching keys from trace */
    double dv[SU_NKEYS]; /* variable array that hold sum of  */
    int nkeys; /* number of keywords whose value should be averaged */
    int ikey; /* index of loop */

    int nt; /* number of data points on trace	*/
    int nsegy; /* number of bytes in the segy		*/
    Value val; /* value of key in current gather	*/
    Value valnew; /* value of key in trace being treated	*/
    int fold; /* number of traces stacked		*/
    int *nnz; /* number of non-zero values stacked	*/
    float normpow; /* divide by nnz[i]^normpow		*/
    int newtracl; /* tracl for stacked traces		*/
    int repeat; /* flag for stack trace repeating	*/
    int nrepeat; /* no. of times stack trace is repeated	*/
    int verbose; /* verbose flag				*/

    segy intrace, outtrace;


    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);


    /* Set parameters */
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparfloat("normpow", &normpow)) normpow = 1.0;
    if (!getparstring("key", &key)) key = "cdp";
    if (!getparint("repeat", &repeat)) repeat = 0;
    if (!getparint("nrepeat", &nrepeat)) nrepeat = 10;

    //if (STREQ(key,"offset"))		is_offset=cwp_true;

    type = hdtype(key);
    indx = getindex(key);

    if (!getparstring("fold", &keyf)) keyf = "counit";
    typef = hdtype(keyf);
    indxf = getindex(keyf);

    /* Get "average" keywords whose values are averaged */
    if ((nkeys = countparval("average")) != 0) {
        getparstringarray("average", keya);
        /* get types and indexes corresponding to the keys */
        for (ikey = 0; ikey < nkeys; ++ikey) {
            typea[ikey] = hdtype(keya[ikey]);
            indxa[ikey] = getindex(keya[ikey]);
        }
    }


    /* Set up for first trace (must compare new key field each time) */
    if ((nsegy = gettr(&intrace)) == 0) err("can't get first trace");
    nt = intrace.ns;
    memcpy((void *) &outtrace, (const void *) &intrace, nsegy);
    nnz = ealloc1int(nt);
    {
        register int i;
        for (i = 0; i < nt; i++) {
            if (intrace.data[i] != 0.0) nnz[i] = 1;
            else nnz[i] = 0;
        }
    }

    if (nkeys > 0) { /* loop over key fields to be averaged */
        for (ikey = 0; ikey < nkeys; ++ikey) {
            gethval(&intrace, indxa[ikey], &vala[ikey]);
            dv[ikey] = vtod(typea[ikey], vala[ikey]);
        }
    }

    /*gethval(&intrace, indxf, &valf);
    if ( (fold = vtoi(typef, valf)) <= 0) fold = 1; */
    fold = 1;

    /* Loop over traces */
    newtracl = 1;
    gethval(&intrace, indx, &val);
    while (nsegy) { /* While previous trace non-empty */
        nsegy = gettr(&intrace);
        gethval(&intrace, indx, &valnew);
        if (valcmp(type, val, valnew) || !nsegy) {
            /* Either val and valnew differ, indicating a  */
            /* new gather or nsegy is zero, indicating the */
            /* end of the traces.                          */
            if (verbose) {
                fprintf(stderr, "val=");
                fprintfval(stderr, type, val);
                fprintf(stderr, "\tfold=%d\n", fold);
            }

            /* Add header info and output stack */
            outtrace.tracl = newtracl++;
            //if(!is_offset) outtrace.offset = 0;
            valf.h = fold;
            puthval(&outtrace, indxf, &valf);

            if (nkeys > 0) { /* loop over key fields to be averaged */
                for (ikey = 0; ikey < nkeys; ++ikey) {
                    setval(typea[ikey], &vala[ikey], dv[ikey] / fold);
                    puthval(&outtrace, indxa[ikey], &vala[ikey]);
                }
            }

            if (normpow && fold != 1) {
                register int i;
                for (i = 0; i < nt; ++i) {
                    float nnzi = nnz[i];
                    if (nnzi)
                        outtrace.data[i] /= pow(nnzi, normpow);
                }
            }
            if (repeat) {
                register int i;
                for (i = 0; i < nrepeat; i++)
                    puttr(&outtrace);
            } else puttr(&outtrace);

            /* Set up for next gather */
            memcpy((void *) &outtrace,
                    (const void *) &intrace, nsegy);
            {
                register int i;
                for (i = 0; i < nt; i++) {
                    if (intrace.data[i] != 0.0) nnz[i] = 1;
                    else nnz[i] = 0;
                }
            }

            fold = 1;

            if (nkeys > 0) { /* loop over key fields to be averaged */
                for (ikey = 0; ikey < nkeys; ++ikey) {
                    gethval(&intrace, indxa[ikey], &vala[ikey]);
                    dv[ikey] = vtod(typea[ikey], vala[ikey]);
                }
            }

            val = valnew;
        } else { /* still in same gather */
            register int i;
            for (i = 0; i < nt; ++i) {
                float datum = intrace.data[i];
                if (!(datum == 0.0)) ++nnz[i];
                outtrace.data[i] += datum;
            }
            ++fold;
            if (nkeys > 0) { /* loop over key fields to be averaged */
                for (ikey = 0; ikey < nkeys; ++ikey) {
                    gethval(&intrace, indxa[ikey], &vala[ikey]);
                    dv[ikey] += vtod(typea[ikey], vala[ikey]);
                }
            }
        }
    }

    free1int(nnz);
    return (CWP_Exit());
}
