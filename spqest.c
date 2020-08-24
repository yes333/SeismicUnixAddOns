/* Copyright (c) RWS, Read Well Services, Oslo, 2011.*/
/* All rights reserved.                       */

/* SPQEST: $Revision: 1.0 $ ; $Date: 2011/04/09 18:32:28 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPQEST - Q factor ESTimation using a panel of VSP data                        ",
"                                                                               ",
" spqest <stdin >stdout  [optional parameters]                                  ",
"                                                                               ",
" Optional parameters:                                                          ",
" nmax=256              max. number of traces input                             ",
"                                                                               ",
" fmin=6.0      [Hz]    lower limit of frequency band use for Q estimation      ",
" fstart=40     [Hz]    starting upper limit of frequency band                  ",
" fend=80       [Hz]    ending upper limit of frequency band                    ",
" df=5          [Hz]    frequency increment of upper limit                      ",
"                                                                               ",
" ftb=gwdep             key holding ftb                                         ",
" scale=0.1             scaling factor to convert ftb key value to ms           ",
"                                                                               ",
" qmin=0                min. value to be reset to if a Q estimate is lower      ",
" qmax=200              max. value to be reset to if a Q estimate is higher     ",
"                       median value is determined within qmin ~ qmax           ",
"                                                                               ",
" verbose=0             >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
"  Input is usually a common shot gather or raw stack of VSP sorted by receiver ",
"  depth. The first input trace is set to be the reference trace, which requires",
"  high S/N ratio with clean FTB pulse. Estimated Q values are output as seismic",
"  trace with first one containing the median values over all traces (receivers)",
"  which have a value within the range between qmin and qmax.                   ",
"                                                                               ",
" Typical flow example looks like:                                              ",
"                                                                               ",
" ......  |\\",
" suchw key1=tstat key2=gwdep a=-1000 b=0.1 |\\",
" sustatic hdrs=1 |\\",
" sumute mode=0 xmute=0 tmute=0.95 ntaper=10 type=4 |\\",
" sumute mode=1 xmute=0 tmute=1.25 ntaper=20 type=4 |\\",
" suwind tmin=0.9 tmax=1.9 |\\",
" spfk dim=1 |\\",
" suamp mode=amp |\\",
" spqest ftb=gwdep scale=0.1 fmin=4 fstart=50 fend=100 df=5 verbose=1 |\\",
" # plot graphs of Q of every f range for quick overview",
" suflip flip=-1 | suflip flip=3 |\\",
" suxgraph f1=1 d1=1 grid2=solid &",
" # to print out Q values for frequency range up to 55 Hz",
" suqest ...   |\\",
" suwind tmin=55 tmax=55 |\\",
" sustrip |\\",
" b2a n1=1 |\\",
" awk '{printf(\"%3.0f \", $1)}'",
" #or to print out median Q values for all frequency ranges",
" suqest ...  |\\",
" suwind key=qdel max=1 |\\",
" sustrip |\\",
" b2a n1=1 |\\",
" awk '{printf(\"%3.0f\\n\", $1)}' ",
"                                                                               ",
" assuming here key gwdep holding ftb in unit tenth of ms                       ",
"                                                                               ",
" Background info about the algorithm used:                                     ",
"  amplitude decay caused by Q for i-th trace/receiver can be approximated as   ",
"                                                                               ",
"   A(T(i), f) = A0(T0, f) * C(i) * exp(-PI*f*T(i)/Q)     C(i) correction factor",
"                                                                               ",
"  divided by j-th trace and take logarithm                                     ",
"                                                                               ",
"   log(A(T(i), f)/(T(j), f)) = log(C(i)/C(j)) + (PI*(T(j) - T(i))/Q) * f       ",
"                                                                               ",
"  the unknown coefficients a & b of this linear equation    c = b + a * f      ",
"                                                                               ",
"  cab be obtained by least square fitting over frequency range  (fmin ~ fmax)  ",
"                                                                               ",
"  thus get        Q = PI * (T(j) - T(i) / a                                    ",
"                                                                               ",
"                                                                               ",
" Version 1.0.0 last modified March. 2011 by Sanyu Ye                           ",
"                                                                               ",
 NULL};

/* Credits:
 *      RWS: Sanyu Ye, March. 2011
 *
 *
 * Notes:
 *
 */
/**************** end self doc *******************************************/

// forward declare prototype
float quick_select(float *arr, int n);
float CalcQ(int itr, int ifmin, int ifmax, float df, float* ftb, float** trdata);

int main(int argc, char **argv)
{
    cwp_String key, key5;       /* header key word from segy.h */
    cwp_String type, type5;     /* type of key	*/
    int index, index5;          /* index of key	*/
    Value val, valnew, vala;    /* value of key				*/
    int nsegy;                  /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;    /* number of total input traces and gathers */
    int itr, iq; /* counters				*/

    int verbose;
    int nmax; /* max dimensions of input/output trace array */

    float scale, qmin, qmax;
    float df, dfmax, fmin, f1, f2;

    segy tr;
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* Get parameters  */
    if (!getparint("verbose", &verbose))    verbose = 0;
    if (!getparint("nmax", &nmax))    nmax = 256;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    int nf = tr.ns;   // number of frequencies
    df = (tr.d1 > 0.0 && tr.d1 < 2.0) ? tr.d1 : tr.dz;


    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "gdel";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* get SU keys holding window parameters */
    if (!getparstring("ftb", &key5)) key5 = "gwdep";
    type5 = hdtype(key5);
    index5 = getindex(key5);

    if (!getparfloat("fmin",  &fmin)) fmin = 6.0;
    if (!getparfloat("fstart", &f1))    f1 = 40.0;
    if (!getparfloat("fend",   &f2))    f2 = 80.0;
    if (!getparfloat("df", &dfmax))  dfmax =  5.0;
    if (!getparfloat("scale", &scale)) scale = 0.1;
    if (!getparfloat("qmin", &qmin)) qmin = 0.0;
    if (!getparfloat("qmax", &qmax)) qmax = 200.0;

    int ndf = NINT(dfmax/df);           // number of frequency increments
    int nq = NINT((f2 - f1)/df/ndf) + 1; // number of Q/F-windows to be estimated/used along single f trace
    int ifmin = NINT(fmin/df);          // starting sample of lower limit of estimation window in f
    int if1 = NINT(f1/df);              // starting sample of upper limit of estimation window in f

    /* Allocate memory */
    float* ftb = ealloc1float(nmax);  // linear taper to merge overlapping windows
    float* qtmp = ealloc1float(nmax);
    float** q = ealloc2float(nq, nmax);
    float** trdata = ealloc2float(nf, nmax);
    segyhdr* hdrs1d = (segyhdr*) ealloc1(nmax, HDRBYTES);
    /* zero out data memory */
    memset(ftb, 0, nmax*FSIZE);
    memset(*q, 0, nmax*nq*FSIZE);
    memset(*trdata, 0, nmax*nf*FSIZE);
    memset(hdrs1d, 0, nmax*HDRBYTES);

    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = ntotal = 0;
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES ) { ///&& !valcmp(type, val, valnew)) { /* same key and more data */
            if (ntr > nmax - 1) err("  array dimension too small (%d < %d) traces input", nmax, ntr+1);
            memcpy(&hdrs1d[ntr], &tr, HDRBYTES);
            memcpy(trdata[ntr], tr.data, FSIZE*tr.ns);

            gethval(&tr, index5, &vala);
            ftb[ntr] = 0.001*scale*vtof(type5, vala);

            ++ntr;
        } else { // new gather or END_OF_FILE
            ++ngather;
            ntotal += ntr;
            if (verbose > 1) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, val));

            // loop over traces and calc Q
            for (itr=1; itr<ntr; ++itr) {
                for (iq=0; iq<nq; ++iq) {
                    float qest = CalcQ(itr, ifmin, if1 + iq*ndf, df, ftb, trdata);
                    if (qest < qmin) q[itr][iq] = qmin;
                    else if (qest > qmax) q[itr][iq] = qmax;
                    else  q[itr][iq] = qest;
                }
            }

            // find out median Q for every frequency band used
            for (iq=0; iq<nq; ++iq) {
                memset(qtmp, 0, nmax*FSIZE);
                int nq = 0;
                for (itr=1; itr<ntr; ++itr) {
                    if (q[itr-1][iq] > qmin && q[itr-1][iq] < qmax) qtmp[nq++] = q[itr-1][iq];
                }
                q[0][iq] = quick_select(qtmp, nq);
            }

            //  output
            for (itr=0; itr<ntr; ++itr) {
                   memcpy(&tr, &hdrs1d[itr], HDRBYTES); // copy back header to output trace
                   memcpy(tr.data, q[itr], nq*FSIZE);
                   tr.ns = nq;
                   tr.f1 = tr.fz = f1;
                   tr.d1 = tr.dz = dfmax;
                   //tr.f2 = itr + 1;
                   //tr.d2 = 1;
                   //tr.dt = 0;
                   tr.delrt = NINT(if1 * df) ;
                   tr.trid = 130;  // depth data
                   puttr(&tr);
            }

            /* zero out data memory */
            memset(ftb, 0, nmax*FSIZE);
            memset(*q, 0, nmax*nq*FSIZE);
            memset(*trdata, 0, nmax*nf*FSIZE);
            memset(hdrs1d, 0, nmax*HDRBYTES);

            val = valnew;
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces of %d gathers are processed", ntotal, ngather);

    return (CWP_Exit());
}

float CalcQ(int itr, int ifmin, int ifmax, float df, float* ftb, float** trdata)
{
    int i;

    int N = ifmax - ifmin + 1;

    float F = 0.0, F2 = 0.0, D = 0.0, FD = 0.0;
    for (i=ifmin; i<=ifmax; ++i) {
        float f = i*df;
        F  += f;
        F2 += f*f;
        D  +=     logf(trdata[0][i]/trdata[itr][i]);
        FD += f * logf(trdata[0][i]/trdata[itr][i]);
    }

    float a = (N * FD - D * F)/(N * F2 - F * F);

    return PI * (ftb[itr] - ftb[0]) / a;
}

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

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

