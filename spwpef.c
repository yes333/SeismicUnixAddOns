/* Copyright (c) Read Well Services, Oslo, 2011.*/
/* All rights reserved.                       */

/* SPWPEF: $Revision: 1.0.0 $ ; $Date: 2011/03/09 18:32:28 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPWPEF - Wiener predictive error filtering for 2/3D gather with variable windows",
"                                                                               ",
" spwpef <stdin >stdout  [optional parameters]                                  ",
"                                                                               ",
" Optional parameters:                                                          ",
" key=cdpt              gather sorting key                                      ",
" output=5              =5 output filtered traces                               ",
"                       =4 output prediction error filters                      ",
"                       =3 output Wiener filters                                ",
"                       =2 output weighted autocorrelations                     ",
"                       =1 output autocorrelations                              ",
"                                                                               ",
" nxmax=tr.sfs          max. number of trace in x direction                     ",
" nymax=tr.sfe          max. number of trace in y direction                     ",
"                                                                               ",
" trid=trid             key holding trace id                                    ",
" minlag=laga           key holding first/prediction lag of filter in milisecond",
" maxlag=lagb           key holding second lag, its difference to minlag is also",
"                           the operator length of Wiener filter                ",
" mincorr=stas          start of autocorrelation window in milisecond           ",
" maxcorr=stae          end of autocorrelation window                           ",
" nw=1,1                trace window in x/y direction for mixing autocorrelation",
" mix=1                 =1 boxcar; =2 Gausian; 3 cosine                         ",
"			  weighting function for averaging autocorrelations	",
" pnoise=0.1     [%]    percentage of relative additive noise level             ",
" verbose=0             >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
"                                                                               ",
" This is an extented version of SUPEF, capable of handling 3D gather and varying",
" parameters for both Wiener filter and autocorrelation window. These parameters",
" should be set/loaded using SUSHW. Input gathers are not required to be regularly",
" sampled, but sorted by shotline and numbered in X and Y position using SUSETGATHER",
" with index stored in keywords nhs and nvs. To improve S/N ratio, autocorrelations",
" of adjacent traces are weighted and summed together. The trace window in X/Y  ",
" is defined by parameter nw=nwx,nwy.                                           ",
" If gathers are not regularly sampled with varying number of traces, nxmax and ",
" nymax must be set big enough to hold the largest one. Only traces with non-zero",
" positive trace id value will be processed and output. If autocorrelation of a ",
" trace becomes zero, in case the correlation window is muted, the trace id is  ",
" negated and only processed if trace/autocorrlation mixing is enabled.         ",
"                                                                               ",
" Typical flow example looks like:                                              ",
"                                                                               ",
" spdbread select=\"cdpt+|ep+|sx+...\" |\\",
" ......  |\\",
" susetgather key=cdpt xkey=sx |\\",
" sushw match=offset key=laga,lagb,stas,stae infile=windows.tbl |\\",
" spwpef key=cdpt nw=11,5 |\\",
" su...  ",
"                                                                               ",
" assuming here cdpt hold receiver no, ep shotline number                       ",
"                                                                               ",
" More about parameters min and maxlag see selfdoc of SUPEF          		",
"                                                                               ",
" Version 1.0.0 last modified March. 2011 by Sanyu Ye                             ",
"                                                                               ",
 NULL};

/* Credits:
 *	CWP: Shuki Ronen, Jack K. Cohen, Ken Larner
 *      RWS: Sanyu Ye, Feb. 2011
 *
 *      Technical Reference:
 *	A. Ziolkowski, "Deconvolution", for value of maxlag default:
 *		page 91: imaxlag < nt/10.  I took nt/20.
 *
 * Notes:
 *	The prediction error filter is 1,0,0...,0,-wiener[0], ...,
 *	so no point in explicitly forming it.
 *
 *	If imaxlag < 2*iminlag - 1, then we don't need to compute the
 *	autocorrelation for lags:
 *		imaxlag-iminlag+1, ..., iminlag-1
 *	It doesn't seem worth the duplicated code to implement this.
 *
 */
/**************** end self doc *******************************************/


// forward declare prototype

#include "wienerhlp.c"


int main(int argc, char **argv)
{
    cwp_String key, key1, key2, key3, key4, key5;     /* header key word from segy.h */
    cwp_String type, type1, type2, type3, type4, type5;    /* type of key	*/
    int index, index1, index2, index3, index4, index5;     /* index of key	*/
    Value val, valnew, vala;  /* value of key				*/
    int nsegy;              /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;    /* number of total input traces and gathers */
    int nt;                 /* number of points on trace		*/
    int ix, iy, itr, it; /* counters				*/
    int output, mix;
    int nlag;       /* length of Wiener filter in samples	*/
    int ncorr;      /* length of corr window in samples	*/
    int lcorr;      /* length of autocorr in samples	*/

    int verbose;
    int nc, nw[2];   /* mixing/weighting windows			*/
    int nx, ny, nxmax, nymax; /* max dimensions of input/output trace array */

    float dt;       /* time sample interval (sec)	*/
    float pnoise;   /* pef additive noise level		*/
    float minlag, maxlag, mincorr, maxcorr, samplims;
    float* crosscorr;
    segy tr, tro;

    const float MIN_FILTER_LENGTH = 60.0;
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* Get parameters  */
    if (!getparint("verbose", &verbose))    verbose = 0;
    if (!getparint("mix", &mix))            mix = 1;
    if (!getparint("output", &output))      output = 5;
    if (!getparfloat("pnoise", &pnoise))    pnoise = 0.1;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    nt = tr.ns;
    dt = ((double) tr.dt) / 1000000.0;
    samplims = 1000.0*dt; // sampling interval in ms

    nx = tr.sfs;
    ny = tr.sfe;
    if(!getparint("nxmax", &nxmax)) nxmax = nx;
    if(!getparint("nymax", &nymax)) nymax = ny;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* get SU keys holding window parameters */
    if (!getparstring("minlag", &key1)) key1 = "laga";
    type1 = hdtype(key1);
    index1 = getindex(key1);
    if (!getparstring("maxlag", &key2)) key2 = "lagb";
    type2 = hdtype(key2);
    index2 = getindex(key2);
    if (!getparstring("mincorr", &key3)) key3 = "stas";
    type3 = hdtype(key3);
    index3 = getindex(key3);
    if (!getparstring("maxcorr", &key4)) key4 = "stae";
    type4 = hdtype(key4);
    index4 = getindex(key4);
    if (!getparstring("trid", &key5)) key5 = "trid";
    type5 = hdtype(key5);
    index5 = getindex(key5);

    // mixing window for autocorrelation
    nc = countparval("nw");
    if (nc > 0) {
        getparint("nw", nw);
        if (nc == 1) nw[1] = 1;
        nw[0] = 2*(nw[0]/2) + 1;
        nw[1] = 2*(nw[1]/2) + 1;
    } else {
        nw[0] = nw[1] = 1;
    }

    /* Allocate memory */
    int** trid = ealloc2int(nxmax, nymax);
    int** iminlag = ealloc2int(nxmax, nymax);
    int** imaxlag = ealloc2int(nxmax, nymax);
    int** imincorr = ealloc2int(nxmax, nymax);
    int** imaxcorr = ealloc2int(nxmax, nymax);
    float* AvgAC = ealloc1float(nt);   // averaged/mixed autocorr
    float* wiener = ealloc1float(nt);
    float* spiker = ealloc1float(nt);
    float** weight = ealloc2float(nw[0], nw[1]);
    float*** autocorr = ealloc3float(nt, nxmax, nymax);
    // allocate memory input/output data traces
    float*** trdata = ealloc3float(nt, nxmax, nymax);
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nxmax, nymax, HDRBYTES);
    /* zero out data memory */
    memset(*trid, 0, nxmax*nymax*ISIZE);
    memset(*iminlag, 0, nxmax*nymax*ISIZE);
    memset(*imaxlag, 0, nxmax*nymax*ISIZE);
    memset(*imincorr, 0, nxmax*nymax*ISIZE);
    memset(*imaxcorr, 0, nxmax*nymax*ISIZE);
    memset(*weight, 0, nw[0]*nw[1]*FSIZE);
    memset(**autocorr, 0, nymax*nxmax*nt*FSIZE);
    memset(**trdata, 0, nymax*nxmax*nt*FSIZE);
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

    // calculate weighting matrix
    CalcMix(mix, nw[0], nw[1], weight);
    
    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = ntotal = 0;
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            iy = tr.nvs - 1;
            ix = tr.nhs - 1;
            if (ntr > nxmax*nymax - 1) err("  array dimension too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nxmax*nymax, ntr+1, ngather+1, key, vtoi(type, val));
            if (ix > nxmax - 1) err("  array dimension nxmax too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nxmax, ix+1, ngather+1, key, vtoi(type, val));
            if (iy > nymax - 1) err("  array dimension nymax too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nymax, iy+1, ngather+1, key, vtoi(type, val));
            memcpy(&hdrs2d[iy][ix], &tr, HDRBYTES);
            memcpy(trdata[iy][ix], tr.data, FSIZE*tr.ns);

            gethval(&tr, index1, &vala);
            minlag = vtof(type1, vala);
            gethval(&tr, index2, &vala);
            maxlag = vtof(type2, vala);
            gethval(&tr, index3, &vala);
            mincorr = vtof(type3, vala);
            gethval(&tr, index4, &vala);
            maxcorr = vtof(type4, vala);
            gethval(&tr, index5, &vala);
            trid[iy][ix] = vtoi(type5, vala);
            iminlag[iy][ix] = NINT(0.001*minlag/dt);
            imaxlag[iy][ix] = NINT(0.001*maxlag/dt);
            imincorr[iy][ix] = NINT(0.001*mincorr/dt);
            imaxcorr[iy][ix] = NINT(0.001*maxcorr/dt);
            // check for consistency
            if (minlag < samplims || maxlag - minlag < samplims || maxlag >= (nt-1)*samplims)
                err("  Inconsistent min/max lag (%1.0f < %1.0f < %1.0f <= %1.0f) for %d of %d-th gather (%s=%d)",
                        samplims, minlag, maxlag, (nt-1)*samplims, ntr+1, ngather+1, key, vtoi(type, val));
            if (maxlag - minlag < MIN_FILTER_LENGTH)
                if (verbose > 0) warn(" Filter may be too short ( %1.0f - %1.0f < %1.0f) for %d of %d-th gather (%s=%d)",
                         maxlag, minlag, MIN_FILTER_LENGTH, ntr+1, ngather+1, key, vtoi(type, val));
            if (mincorr < 0 || maxcorr - mincorr < samplims || maxcorr > (nt-1)*samplims)
                err("  Inconsistent autocorrelation window not conform to (0 <= %1.0f < %1.0f <= %1.0f) for %d of %d-th gather (%s=%d)",
                        mincorr, maxcorr, (nt-1)*samplims, ntr+1, ngather+1, key, vtoi(type, val));
            if (maxcorr - mincorr < MIN_FILTER_LENGTH)
                if (verbose > 0) warn("  Too short autocorrelation window (%1.0f - %1.0f < %1.0f) for %d of %d-th gather (%s=%d)",
                        maxcorr, mincorr, MIN_FILTER_LENGTH, ntr+1, ngather+1, key, vtoi(type, val));
            if (trid[iy][ix] < 0)
                warn("  Negative trace id=%d for %d of %d-th gather (%s=%d)",
                        trid[iy][ix], ntr+1, ngather+1, key, vtoi(type, val));

            ++ntr;
        } else { // new gather or END_OF_FILE
            ++ngather;
            ntotal += ntr;
            if (verbose > 1) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, val));

            // compute autocrrelation for all input traces
            itr = 0;
            for (iy=0; iy<nymax; ++iy) {
                for (ix=0; ix<nxmax; ++ix) {
                    ++itr;
                    if (trid[iy][ix] <= 0) continue; // skip zero trace
                    ncorr = imaxcorr[iy][ix] - imincorr[iy][ix] + 1;
                    lcorr = imaxlag[iy][ix] + 1;
                    xcor(ncorr, 0, &trdata[iy][ix][imincorr[iy][ix]], ncorr, 0, &trdata[iy][ix][imincorr[iy][ix]], lcorr, 0, autocorr[iy][ix]);
                    if (autocorr[iy][ix][0] == 0.0) {
                        trid[iy][ix] = -trid[iy][ix];  // negate id for later identification
                        if (verbose > 3) warn("  Zero autocorrelation for trace %d of %d-th gather (%s=%d)",
                                itr, ngather, key, vtoi(type, val));
                    }
                }
            }

            itr = 0;
            for (iy=0; iy<nymax; ++iy) {
                for (ix=0; ix<nxmax; ++ix) {
                    if (trid[iy][ix] == 0) continue; // skip zero trace
                    memcpy(&tro, &hdrs2d[iy][ix], HDRBYTES);  // copy trace head
                    if (output == 1) { // output auto correlation
                        memcpy(tro.data, autocorr[iy][ix], nt*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    /* compute filter sizes and correlation number */
                    /* zero out filter vectors */
                    memset(wiener, 0, nt*FSIZE);
                    memset(spiker, 0, nt*FSIZE);
                    memset(AvgAC, 0, nt*FSIZE);
                    nlag = imaxlag[iy][ix] - iminlag[iy][ix] + 1;
                    ncorr = imaxcorr[iy][ix] - imincorr[iy][ix] + 1;
                    lcorr = imaxlag[iy][ix] + 1;
                    // compute avaraging autocorrelation
                    memset(AvgAC, 0, nt*FSIZE);
                    if (nw[0] > 1 || nw[1] > 1) {
                        CalcAvgAC(nxmax, nymax, trid, autocorr, nw[0], nw[1], weight, ix, iy, lcorr, AvgAC);
                    } else {
                        memcpy(AvgAC, autocorr[iy][ix], lcorr*FSIZE);
                    }

                    if (output == 2) { // output weighted auto correlation
                        memcpy(tro.data, AvgAC, nt*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    AvgAC[0] *= 1.0 + 0.01*pnoise; /* additive noise/whitening */

                    /* Set pointer to "cross" correlation */
                    crosscorr = AvgAC + iminlag[iy][ix];

                    /* Get inverse filter by Wiener-Levinson */
                    if (AvgAC[0] != 0.0) stoepf(nlag, AvgAC, crosscorr, wiener, spiker);

                    if (output == 3) { // output wiener filter
                        memcpy(tro.data, wiener, nt*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    if (output == 4) { // output prediction error filter
                        memset(tro.data, 0, nt*FSIZE);
                        tro.data[0] = 1.0;
                        for(it=0; it<nlag; ++it) tro.data[iminlag[iy][ix] + it] = -wiener[it];
                        puttr(&tro);
                        continue;
                    }

                    /* Convolve pefilter with trace - don't do zero multiplies */
                    for (it = 0; it < nt; ++it) {
                        register int j;
                        register int n = MIN(it, imaxlag[iy][ix]);
                        register float sum = trdata[iy][ix][it];

                        for (j = iminlag[iy][ix]; j <= n; ++j)
                            sum -= wiener[j - iminlag[iy][ix]] * trdata[iy][ix][it - j];

                        tro.data[it] = sum;
                    }

                    puttr(&tro);
                }
            }

            /* zero out data memory */
            memset(*trid, 0, nxmax*nymax*ISIZE);
            memset(*iminlag, 0, nxmax*nymax*ISIZE);
            memset(*imaxlag, 0, nxmax*nymax*ISIZE);
            memset(*imincorr, 0, nxmax*nymax*ISIZE);
            memset(*imaxcorr, 0, nxmax*nymax*ISIZE);
            memset(**autocorr, 0, nymax*nxmax*nt*FSIZE);
            memset(**trdata, 0, nymax*nxmax*nt*FSIZE);
            memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

            val = valnew;
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces of %d gathers are processed", ntotal, ngather);

    return (CWP_Exit());
}
