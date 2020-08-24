/* Copyright (c) READ Well Services, 2009.*/
/* All rights reserved.                       */

/* SPFTB: $Revision: 1.1 $ ; $Date: 2010/03/17 22:58:42 $	*/

#include "su.h"
#include "segy.h"
#include <segyhdr.h>
#include <cwp.h>
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPFTB - Determine First Time Break                                            ",
"                                                                               ",
" spftb <stdin >stdout [optional parameters]                                    ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
"   key=cdp             filter is applied gather-wise with constant key value   ",
"                                                                               ",
"   nxmax=1024          max. number of traces expected in a gather              ",
"                                                                               ",
"   mode=1              =0 amplitude as trigger mode                            ",
"                       =1 gradient, amplitude difference between two samples   ",
"                       =2 peak, i.e. local maximum                             ",
"                       =3 trough, i.e. local minimum (with negative amplitude) ",
"                                                                               ",
"   npoints=3           odd number (3 or 5) of samples to determine local peak/trough",
"                                                                               ",
"   f=2.0               threshold relative to rms amplitude of time window      ",
"                                                                               ",
" Parameters controlling reference FTB around which the search window is defined",
"                                                                               ",
"   option=0            =2 use previous FTB                                     ",
"                       =1 use velocity-depth function to predict (not implemented)",
"                       =0 use time window given below                          ",
"   search=-0.2,0.2 [s] search window around reference FTB for option=1/2       ",
"                       if only one value given, the window is symmetric        ",
"   vd=v1,d1,...        velocity-depth function in [m/s] and [m]                ",
"     =1800,0,2500,1000,3000,2000,4000,3000                                     ",
"                                                                               ",
" Search window parameters for option=0 when no reference FTB is given          ",
"   tmin=0.0            start time of the window                                ",
"   tmax=(from header)  end time of the window                                  ",
"   itmin=0             min time sample of the window                           ",
"   itmax=(from header) max time sample of the window                           ",
"   nt=itmax-itmin+1    number of time samples of the window                    ",
"                                                                               ",
"   maxdiff=20 [ms]     max. diviation between determined and median FTB        ",
"                                                                               ",
"   nm=11               odd number of traces as median filter width             ",
"                       =0  no second pass redetermination                      ",
"                                                                               ",
"   skey=laga           key to store FTB [1/10 ms]                              ",
"                                                                               ",
"   verbose=0           >0  echoes information                                  ",
"                                                                               ",
" FTBs are determined in a two-pass procedure:                                  ",
" In the first pass a gather is read in and rms amplitude within the search window",
" is computed for each trace. The search of the FTB also starts within the      ",
" window. When a value for the specified mode exceeds the threshold of rms      ",
" amplitude, the corresponding sample (time) is defined as FTB.                 ",
" In the second pass each FTB is checked against its neighbors for consistancy. ",
" If its deviation from the expected/median value is larger than the maxdiff,   ",
" the FTB is redetermined within a small window centered at the median FTB      ",
" with a half width of maxdiff.                                                 ",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
" spftb < h.su > h-ftb.su tmin=0.05 tmax=0.15 f=2                               ",
"                                                                               ",
" spftb < z.su > z-ftb.su tmin=0.2 tmax=0.6 mode=2                              ",
"                                                                               ",
" Version 2.1.0 last modified Oct. 2011 by Sanyu Ye                             ",
NULL};

/* Credits:
 *
 * RWS: Sanyu Ye, sanyu.ye@readgroup.com
 *
 * Trace header fields accessed: ns, dt, key=key
 * Trace header fields modified: skey
 *
 * Version 1.1  Last updated 2009.10.13
 */
/**************** end self doc ***********************************/

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float quick_select(float *arr, int n);
static void setval( cwp_String type, Value *valp, float fval);
int median(int nm, int ntr, float* indata, float* outdata);
int findFTB(int mode, int npoints, float threshold, int itmin, int itmax, float* tr);

int main(int argc, char **argv)
{
    char *key;              /* header key word from segy.h		*/
    char *type;             /* ... its type				*/
    int index;              /* ... its index			*/
    Value val, valnew;      /* ... its value			*/

    char *skey;              /* header key word for storing FTB		*/
    char *stype;             /* ... its type				*/
    int sindex;              /* ... its index			*/
    Value sval;             /* ... its value			*/

    int i, it;              /* sample counter			*/
    int itr;                /* trace counter			*/
    int ns;                 /* number of time samples 		*/
    int nsegy;              /* number of byte read 		*/
    int ntr, total, ngather;/* number of traces and gathers			*/
    int nxmax;              /* max. number of traces expected in gather */
    int eof=0;              /* END OF File flag			*/

    int npoints;             /* number of time samples for local trough/peak */
    int nt;                 /* number of time samples of the time window */
    int itmin;              /* smallest time sample (zero-based)    */
    int itmax;              /* largest time sample (zero-based)     */
    float tmin, tmax;       /* minimum/maximum time to calculate        	*/
    float tref, search[2];  /* reference FTB and search window around		*/
    float threshold;        /* threshold factor		*/
    float maxdiff;          /* max diviation in millisecond	*/

    float dt;               /* time sampling interval		*/

    int verbose;            /* flag for printing information	*/

    int mode, option;       /* trigger/search mode	*/
    int nm;                 /* no. of traces to median filter	*/

    segy tr, outtr;

    segyhdr *hdrs=NULL;     /* placeholder for incoming trace headers	*/
    float **indata=NULL;    /* temporary array			*/
    float *rmsamp = NULL;   /* array caching rms amplitude	*/
    float *medamp = NULL;   /* array caching median rms amplitude	*/
    float *medftb = NULL;   /* array caching median ftbs	*/
    float *ftb = NULL;      /* array caching tentative ftbs	*/

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("nxmax",&nxmax))	nxmax = 1024;
    if (!getparint("mode",&mode))	mode = 1;
    if (!getparint("option",&option))	option = 0;
    if (!getparint("nm",&nm) )          nm = 11;
    if (nm < 0 || (nm > 0 && !(nm%2) )) err("Filter width (nm=%d) must be a positive odd number", nm);
    if (!getparint("npoints",&npoints) )  npoints = 3;
    else if (npoints != 3 && npoints != 5) {
        npoints = 3;
        warn(" Number of samples is reset to 3 to dtermine local peak/trough");
    }

    if (!getparfloat("f", &threshold))		threshold = 2.0;
    if (!getparfloat("maxdiff", &maxdiff) )      maxdiff = 20.0;

    /* Get remaining parameters */
    if (!getparint("verbose", &verbose))	verbose = 0;

    /* Get info from first trace */
    if (!(nsegy = gettr(&tr))) err("can't read first trace");
    if (!tr.dt) err("dt header field must be set");
    dt   = ((double) tr.dt)/1000000.0;
    ns = (int) tr.ns;

    /* Time gating parameters */
    if (option == 2) {
        int nv = countparval("search");
        if (nv == 2) {
            getparfloat("search", search);
        } else if (nv == 1) {
            getparfloat("search", search);
            search[0] = -ABS(search[0]);
            search[1] = -search[0];
        } else {
            search[0] = -0.2;
            search[1] =  0.2;
        }
    } else if (option == 1) {
        err(" option=1 not implemented yet");
    } else {
        if (getparint("itmin", &itmin)) {
                tmin = itmin*dt;
        } else if (getparfloat("tmin", &tmin)) {
                itmin = NINT(tmin / dt);
        } else {
                itmin = 0;
                tmin = 0;
        }
        if (getparint("itmax", &itmax)) {
                tmax = itmax*dt;
                nt = itmax - itmin + 1;
        } else if (getparfloat("tmax", &tmax)) {
                itmax = NINT(tmax / dt);
                nt = itmax - itmin + 1;
        } else if (getparint("nt", &nt)) {
                itmax = itmin + nt - 1;
                tmax = itmax*dt;
        } else {
                itmax = tr.ns - 1;
                tmax = itmax*dt;
                nt = itmax - itmin + 1;
        }

        /* Check time gating values */
        if (itmin < 0)
                err("itmin=%d should be positive", itmin);
        if (nt > SU_NFLTS)
                err("nt=%d exceeds SU_NFLTS=%d", nt, SU_NFLTS);
        if (itmin > itmax)
                err("itmin=%d, itmax=%d conflict", itmin, itmax);
        if (tr.ns <= itmax)
                err("tr.ns=%d, itmax=%d window cannot extend over the trace length", tr.ns, itmax);
    }
    if (verbose > 1 && option > 0) warn("Search window %.3f ~ %.3f", search[0], search[1]);

    if (verbose > 10) warn(" tmin=%.3f  tmax=%.3f  itmin=%d itmax=%d nt=%d", tmin, tmax, itmin, itmax, nt);


    // allocate and reset memory for input and output traces
    hdrs = (segyhdr*) ealloc1(nxmax, HDRBYTES);
    indata = ealloc2float(ns, nxmax);
    rmsamp = ealloc1float(nxmax);
    medamp = ealloc1float(nxmax);
    medftb = ealloc1float(nxmax);
    ftb = ealloc1float(nxmax);
    if ( indata == NULL || hdrs == NULL)
        err("memory allocation error for input or/and output data");
    /* zero out data memory */
    memset(hdrs, 0, nxmax*HDRBYTES);
    memset(*indata, 0, nxmax*ns*FSIZE);

    /* Get key to store ftb */
    if (!getparstring("skey", &skey)) skey = "laga";
    stype = hdtype(skey);
    sindex = getindex(skey);

    /* Get sorting/gather key  */
    if (!getparstring("key", &key)) key = "cdp";
    type = hdtype(key);
    index = getindex(key);

    gethval(&tr, index, &val);

    /* Read headers and data while getting a count */
    ngather =0;
    ntr = 0;
    total = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nxmax-1) err("\nNumber of traces exceeding nxmax=%d\n", nxmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);
            memcpy(indata[ntr], tr.data, FSIZE*ns);
            val = valnew;
            ++ntr;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (verbose > 1) warn("  processing %d traces processed for %d-th gather (%s=%d)", ntr, ngather, key, vtoi(type, val));

            total += ntr;
            val = valnew;
            // first pass, reset caches
            memset(rmsamp, 0, nxmax*FSIZE);
            memset(medamp, 0, nxmax*FSIZE);
            memset(medftb, 0, nxmax*FSIZE);
            memset(ftb, 0, nxmax*FSIZE);
            for (itr=0; itr<ntr; ++itr) {  // loop over traces
                if (option == 2) {
                    gethval((segy*) &hdrs[itr], sindex, &sval);
                    tref = (float) vtoi(stype, sval)*0.0001;
                    if (tref < 10.0*dt) warn("  possibably invalid reference FTP=%.5f [s] for %d-th trace of %d-th gather (%s=%d)",
                            itr+1, ngather, key, vtoi(type, val));
                    itmin = MAX(0, NINT((tref + search[0])/dt));
                    itmax = MIN(NINT((tref + search[1])/dt), tr.ns - 1);
                }

                // compute rms amplitude
                for (i=0, it=itmin; it<=itmax; ++it) { // loop over samples
                    if ( indata[itr][it] != 0.0) { //skip hard zero
                        rmsamp[itr] += indata[itr][it]*indata[itr][it];
                        ++i;
                    }
                }
                if (i>0) rmsamp[itr] = sqrt(rmsamp[itr]/i);

                i = findFTB(mode, npoints, rmsamp[itr]*threshold, itmin, itmax, indata[itr]);
                if (i>0) { // find FTB
                    float amp = (mode == 1)? (indata[itr][i] - indata[itr][i-1]) : indata[itr][i];
                    ftb[itr] = dt * i;
                    if(verbose > 10) warn("FTB=%6.2f [ms] & trigger amp=%g ratio=%3.1f for %d-th trace of %d-th gather (%s=%d)",
                         ftb[itr]*1000.0, amp, amp/rmsamp[itr], itr+1, ngather, key, vtoi(type, val));
                } else {
                    if (option > 1) {
                        ftb[itr] = tref;
                        if (verbose) warn("keep old FTB=%6.4f [s] for %d-th trace of %d-th gather (%s=%d)", ftb[itr], itr+1, ngather, key, vtoi(type, val));
                    }
                    else if (verbose) warn("no FTB found for %d-th trace of %d-th gather (%s=%d)",  itr+1, ngather, key, vtoi(type, val));
                }
            }

            // second pass, check for ftb consistency
            if (nm > 0) {
                median(nm, ntr, ftb, medftb);
                for (itr=0; itr<ntr; ++itr) {  // loop over traces
                    if (fabs(ftb[itr] - medftb[itr])*1000.0 > maxdiff ) { // difference too big
                        //redetermine within smaller windows
                        itmin = NINT((medftb[itr] - (float) maxdiff / 1000.0) / dt);
                        itmax = NINT((medftb[itr] + (float) maxdiff / 1000.0) / dt);

                        i = findFTB(mode, npoints, rmsamp[itr] * threshold, itmin, itmax, indata[itr]);
                        if (i > 0) { // new ftb found
                            if (verbose > 1) warn("FTB old=%6.2f & new=%6.2f [ms] for %d-th trace of %d-th gather (%s=%d)",
                                    ftb[itr]*1000.0, i * 1000.0 * dt, itr + 1, ngather, key, vtoi(type, val));
                            ftb[itr] = i * dt;
                        } else if (ftb[itr] == 0.0) {
                            if (verbose) warn("Failed to find a new but set FTB to median=%6.2f [ms] for %d-th trace of %d-th gather (%s=%d)",
                                    medftb[itr]*1000.0, itr + 1, ngather, key, vtoi(type, val));
                            ftb[itr] = medftb[itr];
                        } else {
                            if (verbose) warn("Failed to correct FTB old=%6.2f within tolerance of median=%6.2f [ms] for %d-th trace of %d-th gather (%s=%d)",
                                    ftb[itr]*1000.0, medftb[itr]*1000.0, itr + 1, ngather, key, vtoi(type, val));
                        }
                    }
                }
            }

            // output traces
            for(itr=0; itr<ntr; ++itr) {
                memcpy(&outtr, &hdrs[itr], HDRBYTES);
                memcpy(outtr.data, indata[itr], FSIZE*ns);
                setval(stype, &sval, ftb[itr] * 10000.0);  // convert to 10-th millisecond
                puthval(&outtr, sindex, &sval);
                puttr(&outtr);
            }
            // reset counter for next gather
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose > 0) warn("  Totally %d gathers with %d traces processed", ngather, total);

    /* Clean up */
    //if (hdrs) free1(hdrs);
    //if (indata) free2float(indata);
    //if (rmsamp) free1float(rmsamp);
    //if (medamp) free1float(medamp);
    //if (medftb) free1float(medftb);
    //if (ftb) free1float(ftb);

    return(CWP_Exit());
}

//determine ftb, return it as number of samples
int findFTB(int mode, int npoints, float threshold, int itmin, int itmax, float* tr) {
    int i = 0;  // sample of ftb
    int it;

    for (it=itmin; it<=itmax; ++it) { // loop over samples
        switch ( mode ) {
            case 0: // amplitude
                if (tr[it] >= threshold ) i = it;
                break;
            case 1: // gradient
                if ( it > 0 && (tr[it] - tr[it-1]) >= threshold ) i = it;
                break;
            case 2: // peak, find first local maximum and check its amplitude agianst threshold, five point check
                if ( npoints == 3 && it > 1 && it < itmax-1 && tr[it] >= threshold
                        && (tr[it] - tr[it-1]) >= 0.0 && (tr[it] - tr[it+1]) >= 0.0 ) i = it;
                else if ( npoints == 5 && it > 2 && it < itmax-2 && tr[it] >= threshold
                        && (tr[it-1] - tr[it-2]) >= 0.0 && (tr[it+1] - tr[it+2]) >= 0.0
                        && (tr[it] - tr[it-1]) >= 0.0 && (tr[it] - tr[it+1]) >= 0.0 ) i = it;
                break;
            case 3: // trough, find first local minimum and check its amplitude agianst threshold
                if ( npoints == 3 && it > 2 && it < itmax-2 && tr[it] <= -threshold
                        && (tr[it] - tr[it-1]) <= 0.0 && (tr[it] - tr[it+1]) <= 0.0 ) i = it;
                else if ( npoints == 5 && it > 2 && it < itmax-2 && tr[it] <= -threshold
                        && (tr[it-1] - tr[it-2]) <= 0.0 && (tr[it+1] - tr[it+2]) <= 0.0
                        && (tr[it] - tr[it-1]) <= 0.0 && (tr[it] - tr[it+1]) <= 0.0 ) i = it;
                break;
        }
        if (i>0) return i;
    }

    return 0;
}

void setval( cwp_String type, Value *valp, float fval) {
    switch (*type) {
        case 's':
            err("can't set char header word");
            break;
        case 'h':
            valp->h = (short) fval;
            break;
        case 'u':
            valp->u = (unsigned short) fval;
            break;
        case 'l':
            valp->l = (long) fval;
            break;
        case 'v':
            valp->v = (unsigned long) fval;
            break;
        case 'i':
            valp->i = (int) fval;
            break;
        case 'p':
            valp->p = fval;
            break;
        case 'f':
            valp->f = fval;
            break;
        case 'd':
            valp->d = (double) fval;
        default:
            err("unknown type %s", type);
            break;
    }
    return;
}

// determine the median values of a 1-D array
int median(int nm, int ntr, float* indata, float* outdata) {
    int i, j=0, itr;
    int nhalf = (nm - 1)/2;  //half width of filter

    /* Allocate float array to hold samples to be filtered */
    float *dtemp1 = ealloc1float(nm);
    float *dtemp2 = ealloc1float(nm);

    for (itr=0; itr<ntr; ++itr) {  // loop over traces
        memset(dtemp1, 0, nm*FSIZE);
        memset(dtemp2, 0, nm*FSIZE);
        if (itr <= nhalf) {  //starting edge
            for (i=nhalf - itr; i < nm; ++i) {
                dtemp1[i] = indata[itr - nhalf + i];
            }
        } else if (ntr - itr < nhalf + 1) { //ending edge
            for (i=0; i < nhalf + ntr - itr; ++i) {
                dtemp1[i] = indata[itr - nhalf + i];
            }
        }
        else { // within traces
            for (i=0, j=0; (i - j) < nm && (itr - nhalf + i) < ntr; ++i) {
                if (indata[itr - nhalf + i] > 0.0) { // if valid ftb
                    dtemp1[i-j] = indata[itr - nhalf + i];
                } else { // skip and go further to next trace
                    ++j;
                }
            }
        }
        for (i=0, j=0; i < nm; ++i) {
            if (dtemp1[i] != 0.0 ) { // not hard zero, copy
                dtemp2[j++] = dtemp1[i];
            }
        }
        if (j > 0) outdata[itr] = quick_select(dtemp2, j);
    }

    free1float(dtemp1);
    free1float(dtemp2);
    return j;
}


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

float quick_select(float *arr, int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
		if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
		    do ll++; while (arr[low] > arr[ll]) ;
		    do hh--; while (arr[hh]  > arr[low]) ;

		    if (hh < ll)
		    break;

		    ELEM_SWAP(arr[ll], arr[hh]) ;
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]) ;

		/* Re-set active partition */
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}

    return (n%2)? arr[median] : 0.5*(arr[median] + arr[median + 1]);
}

#undef ELEM_SWAP

