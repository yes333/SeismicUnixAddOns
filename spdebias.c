/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPDEBIAS: $Revision: 1.0 $ ; $Date: 2010/05/17 22:58:42 $	*/

#include "su.h"
#include "segy.h"
#include <segyhdr.h>
#include <cwp.h>
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPDEBIAS - Remove bias and overscaling of optical OBC sensor                  ",
"                                                                               ",
" spdebias <stdin >stdout [optional parameters]                                 ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
"   key=ep              gather sorting key                                      ",
"   nxmax=512           max. number of traces expected in a gather              ",
"   wkey=offset         key to select trace window                              ",
"   scale=tr.scalco     scaling factor, default tr.scalco for default key=offset",
"        =1.0           default 1.0 for other than default key=offset           ",
"   max=1500    [m]     max value for key header word                           ",
"   min=-max            min value for key header word                           ",
"                                                                               ",
" Threshold regarding reference amplitude for detecting edges above and below   ",
"                                                                               ",
"   thresa1=8            absolute sample amplitude                              ",
"   thresb1=8                                                                   ",
"   thresa2=4            amplitude difference between two consecutive samples   ",
"   thresb2=4                                                                 ",
"   thresa3=8            rms amplitude of time window                           ",
"   thresb3=8                                                                   ",
"                                                                               ",
" Threshold regarding amplitude ratio of win3 agianst win4 for detecting edges below",
"                                                                               ",
"   thresb5=20           rms amplitude                                          ",
"   thresb4=20           maximum absolute amplitude                             ",
"   thresb6=100000       avaerage amplitude                                     ",
"                                                                               ",
" Parameters to determine reference amplitude: 					",
"                                                                               ",
"   t0=1000     [ms]    time horizon direct water wave arrivals are aligned at  ",
"                       used also as ",
"   win0=40     [ms]    time window to find reference max. amplitude            ",
"   ntrace=10           trace window width                                      ",
"                                                                               ",
" Parameters to determine reference amplitude: 					",
"                                                                               ",
"   t1=1030     [ms]    start of time window to forward check overscale         ",
"                       should be set possibly close to earliest onset of overscaling",
"   t2=2000     [ms]    end of time window to forward check overscale           ",
"                       used also as start of time window to backward check overscale",
"                                                                               ",
" Parameters for various windows: 					",
"                                                                               ",
"   win1=100    [ms]    window length above of forward check    ",
"   win2=40     [ms]    window length below of forward check    ",
"   win3=40     [ms]    window length above of backward check    ",
"   win4=200    [ms]    window length below of forward check    ",
"                                                                               ",
"   nsa=3               number of samples for tapering above overscaling    ",
"   nsb=5               number of samples for tapering below overscaling    ",
"                                                                               ",
"   verbose=0           >0  echoes information				",
" 									",
" Approach: 								",
" ---------								",
" There are typically two types of overscaling, tooth-like impulse or bias after",
" a sudden jump. This program uses approach typical for first break detection,  ",
" using two windows, one before and one after to detect the onset of the sudden ",
" amplitude change. However to be able to tackle the tooth-like overscaling, two",
" sets of windows are used: one forward just little bit after expected direct   ",
" water wave arrival (t0) but closest to the earliest onset of overscaling t1   ",
" and another backwards from time t2 ",
"									",
" Notes: 								",
" ------								",
"                                                                       ",
" Examples:                                                             ",
"									",
" ... |\\                                                               ",
" sureduce mode=2 rv=1.48 t0=1 |\\                                      ",
" spdebias key=duse |\\                                                 ",
" sureduce mode=2 rv=-1.48 t0=1 |\\                                     ",
" ... |\\                                                               ",
"                                                                       ",
"									",
NULL
};

/* Credits:
 *
 * RWS: Sanyu Ye, sanyu.ye@readgroup.com
 *
 * Trace header fields accessed: ns, dt, key=key,wkey
 * Trace header fields modified: nvs when stack=1
 *
 */
/**************** end self doc ***********************************/

// forward prototyp declaration
float quick_select(float *arr, int n);
float calcAvg(int itbeg, int nt, float* data);
float calcMax(int itbeg, int nt, float* data, float base);
float calcRms(int itbeg, int nt, float* data, float base);
float calcMAX(int itr, int ntrace, int itbeg, int nt, int** bad, float** data, float base);
float calcRMS(int itr, int ntrace, int itbeg, int nt, int** bad, float** data, float base);
int debias(int* itw, int ntrace, int ntr, int ns, int nsa, int nsb,
    int it0, int itmin, int itmax, int* win, int** bad, float* thresa, float* thresb, float** indata);


int verbose; /* flag for printing information	*/

int main(int argc, char **argv) {
    char *key; /* header key word from segy.h		*/
    char *type; /* ... its type				*/
    int index; /* ... its index			*/
    Value val, valnew; /* ... its value			*/

    char *wkey; /* header key word for windowing		*/
    char *wtype; /* ... its type				*/
    int windex; /* ... its index			*/
    Value wval; /* ... its value			*/
    int wmin; /* min value    */
    int wmax; /* max value     */
    int *win = NULL; /* array containing flag if a trace within select window     */
    int **bad = NULL; /* array containing sample numbers if a trace is affected by overscaling */

    int itr; /* trace counter			*/
    int ns, nsa, nsb; /* number of time samples 		*/
    int nsegy; /* number of byte read 		*/
    int ntr, total, ngather, badtotal; /* number of traces and gathers			*/
    int nxmax; /* max. number of traces expected in gather */
    int eof = 0; /* END OF File flag			*/

    int ntrace; /* number of traces for calc rms  */
    int itmin; /* smallest time sample (zero-based)    */
    int itmax; /* largest time sample (zero-based)     */
    int   it0, itw[5];  /* time windows in number of samples */
    float t0, t1, t2; /* time window to check        	*/
    float win0, win1, win2, win3, win4; /* time windows to calculate rms or average amplitudes		*/
    float scale; /* (offset) scaling factor		*/
    float dt; /* time sampling interval		*/
    float thresa[3], thresb[6], thres; /* amplitude ratio	*/

    segy tr, outtr;

    segyhdr *hdrs = NULL; /* placeholder for incoming trace headers	*/
    float **indata = NULL; /* temporary array			*/

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparint("nxmax", &nxmax)) nxmax = 512;
    if (!getparint("ntrace", &ntrace)) ntrace = 10;
    if (!getparint("nsa", &nsa)) nsa = 3;
    if ( nsa < 0 ) err("  nsa must >= 0");
    if (!getparint("nsb", &nsb)) nsb = 5;
    if ( nsb < 0 ) err("  nsb must >= 0");
    if (!getparfloat("t0", &t0)) t0 = 1000;
    if (!getparfloat("t1", &t1)) t1 = 1030;
    if (!getparfloat("t2", &t2)) t2 = 2000;
    if (!getparfloat("thresa1", &thres)) thres =  8.0;   thresa[0] = thres;
    if (!getparfloat("thresa2", &thres)) thres =  4.0;   thresa[1] = thres;
    if (!getparfloat("thresa3", &thres)) thres =  8.0;   thresa[2] = thres;
    if (!getparfloat("thresb1", &thres)) thres =  8.0;   thresb[0] = thres;
    if (!getparfloat("thresb2", &thres)) thres =  4.0;   thresb[1] = thres;
    if (!getparfloat("thresb3", &thres)) thres =  8.0;   thresb[2] = thres;
    if (!getparfloat("thresb4", &thres)) thres = 20.0;   thresb[3] = thres;
    if (!getparfloat("thresb5", &thres)) thres = 20.0;   thresb[4] = thres;
    if (!getparfloat("thresb6", &thres)) thres = 1000000.0;   thresb[5] = thres;
    if (!getparfloat("win0", &win0)) win0 = 40;
    if (!getparfloat("win1", &win1)) win1 = 100;
    if (!getparfloat("win2", &win2)) win2 = 20;
    if (!getparfloat("win3", &win3)) win3 = 20;
    if (!getparfloat("win4", &win4)) win4 = 200;
    if (verbose) warn("  Time window to check: tmin=%.3f  tmax=%.3f [s]", t1/1000.0, t2/1000.0);
    if (verbose) warn("  Trace and Time window length for rms amplitude ntrace=%d win0=%3.0f [ms] t0=%.3f [s]", ntrace, win0, t0/1000.0);
    if (verbose) warn("  Thresholds for detecting edges of overscale above %3.1f %3.1f %3.1f   below: %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f",
                          thresa[0], thresa[1], thresa[2], thresb[0], thresb[1], thresb[2], thresb[3], thresb[4], thresb[5]);


    /* Get info from first trace */
    if (!(nsegy = gettr(&tr))) err("can't read first trace");
    if (!tr.dt) err("dt header field must be set");
    dt = ((double) tr.dt) / 1000000.0;
    ns = (int) tr.ns;

    /* time line and window length in number of samples */
    it0   = NINT(0.001*t0 / dt);
    itmin = NINT(0.001*t1 / dt);
    itmax = NINT(0.001*t2 / dt);
    itw[0]  = NINT(0.001*win0 / dt);
    itw[1]  = NINT(0.001*win1 / dt);
    itw[2]  = NINT(0.001*win2 / dt);
    itw[3]  = NINT(0.001*win3 / dt);
    itw[4]  = NINT(0.001*win4 / dt);

    /* Check time gating values */
    if (itmin < 0)
        err("itmin=%d should be positive", itmin);
    if (itmin > itmax)
        err("itmin=%d, itmax=%d conflict", itmin, itmax);
    if (tr.ns <= itmax)
        err("tr.ns=%d, itmax=%d window cannot extend over the trace length", tr.ns, itmax);

    /* Get windowing key  */
    if (!getparstring("wkey", &wkey)) wkey = "offset";
    wtype = hdtype(wkey);
    windex = getindex(wkey);
    if (!getparint("max", &wmax)) wmax =  1500;
    if (!getparint("min", &wmin)) wmin = -wmax;
    if (!getparfloat("scale", &scale)) {
        if (!strcmp(wkey, "offset")) {
            scale = (tr.scalco < 0) ? -1.0 / tr.scalco : (tr.scalco > 0) ? tr.scalco : 1.0;
        } else {
            scale = 1.0;
        }
    }
    if (verbose) warn("  Scaling factor for %s: f=%.3f", wkey, scale);

    /* Get sorting/gather key  */
    if (!getparstring("key", &key)) key = "cdp";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    // allocate and reset memory for input and output traces
    hdrs = (segyhdr*) ealloc1(nxmax, HDRBYTES);
    indata = ealloc2float(ns, nxmax);
    win = ealloc1int(nxmax);
    bad = ealloc2int(2, nxmax);
    /* zero out data memory */
    memset(win, 0, nxmax * sizeof (int));
    memset(*bad, 0, 2*nxmax * sizeof (int));
    memset(hdrs, 0, nxmax * HDRBYTES);
    memset(*indata, 0, nxmax * ns * FSIZE);

    /* Read headers and data while getting a count */
    ngather = ntr = total = badtotal = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nxmax - 1) err("\nNumber of traces exceeding nxmax=%d\n", nxmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);
            memcpy(indata[ntr], tr.data, FSIZE * ns);
            gethval(&tr, windex, &wval);
            if (vtod(wtype, wval) * scale >= wmin && vtod(wtype, wval) * scale <= wmax) win[ntr] = 1;
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (verbose) warn("  processing %d traces %d-th gather (%s=%d)", ntr, ngather, key, vtoi(type, val));

            if (ntr < ntrace) warn("  Fewer traces (%d < %d) than window length for %d-th gather (%s=%d)", ntr, ntrace, ngather, key, vtoi(type, val));

            total += ntr;
            // do debias
            int nbadtr = debias(itw, ntrace, ntr, ns, nsa, nsb, it0, itmin, itmax, win, bad, thresa, thresb, indata);
            badtotal += nbadtr;
            if (verbose) warn(" %d overscaled traces found for %d-th gather (%s=%d)", nbadtr, ngather, key, vtoi(type, val));

            for (itr = 0; itr < ntr; ++itr) {
                memcpy(&outtr, &hdrs[itr], HDRBYTES);
                outtr.muts = bad[itr][0];
                outtr.mute = bad[itr][1];
                memcpy(outtr.data, indata[itr], FSIZE * ns);
                puttr(&outtr);
            }
            val = valnew;
            // reset output data
            memset(win, 0, nxmax * sizeof (int));
            memset(*bad, 0, 2*nxmax * sizeof (int));
            memset(*indata, 0, nxmax * ns * FSIZE);
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn("  %d overscaled traces processed of %d total traces within %d gathers", badtotal, total, ngather);

    /* Clean up */
    if (hdrs) free1(hdrs);
    if (indata) free2float(indata);
    if (bad) free2int(bad);
    if (win) free1int(win);

    return (CWP_Exit());
}

int debias(int itw[5], int ntrace, int ntr, int ns, int nsa, int nsb, int it0, int itmin,
        int itmax, int* win, int** bad, float* thresa, float* thresb, float** indata) {
    int   itr, it, itb;
    float avg1, avg2, avg3, avg4, avga, avgb;
    float rms0, rms1, rms2, rms3, rms4, rmsb;
    float max0, max1, max2, max3, max4, maxb;
    float diff, sign;

    int nbadtr = 0;
    for (itr = 0; itr < ntr; ++itr) { // loop over traces
        // skip if trace out of selection window, however still used for calculating base amplitude
        if (!win[itr]) continue;
        // calc rms of previous ntrace within window itw[0] (win0) starting itmin (t0)
        rms0 = calcRMS(itr, ntrace, it0, itw[0], bad, indata, 0);
        max0 = calcMAX(itr, ntrace, it0, itw[0], bad, indata, 0);
        avga = calcAvg(itmin - itw[4], itw[4], indata[itr]);
        avgb = calcAvg(itmax, itw[4], indata[itr]);
        rmsb = calcRms(itmax, itw[4], indata[itr], avgb);

        // start search backword for back edge of overscaling
        for ( itb = itmax; itb > itmin; --itb ) {
            avg3 = calcAvg(itb - itw[3] + 1, itw[3], indata[itr]);
            avg4 = calcAvg(itb + 1, itw[4], indata[itr]);
            rms3 = calcRms(itb - itw[3] + 1, itw[3], indata[itr], avgb);
            rms4 = calcRms(itb + 1, itw[4], indata[itr], avgb);
            max3 = calcMax(itb - itw[3] + 1, itw[3], indata[itr], avgb);
            max4 = calcMax(itb + 1, itw[4], indata[itr], avgb);

            float ratio = ( rms3 > rms4 )? rms3/rms4 : rms4/rms3;
            diff = ABS(indata[itr][itb] - indata[itr][itb + 1]);
            if (     ABS(indata[itr][itb] - avgb) > thresb[0]*max0
                  || ( diff > thresb[1]*max0 && ratio > thresb[1] )
                  || rms3 > thresb[1]*rms0
                  || max3 > thresb[2]*max0
                  || rms3 > thresb[3]*rms4
                  || max3 > thresb[4]*max4
                  || ( avg4 != 0.0 && ABS(avg3) > thresb[5]*ABS(avg4) )
               )
            {
                bad[itr][1] = itb;  // found front edge of overscaling
                if (verbose > 1) {
                    warn(" Back edge at sample %d on %d-th trace: amp=%f diff=%f\n max=%f %f %f rms=%f %f %f Avg=%f %f %f",
                         itb, itr, indata[itr][itb], diff, max3, max4, max0, rms3, rms4, rms0, avg3, avg4, avgb);
                    if ( ABS(indata[itr][itb] - avgb) > thresb[0]*max0 ) warn(" ABS(indata - avgb) > thresb1*max0 ");
                    if ( diff > thresb[1]*max0 && ratio > thresb[1] ) warn("diff > thresb2*max0 && ratio > thresb2");
                    if ( rms3 > thresb[1]*rms0 ) warn("rms3 > thresb2*rms0 ");
                    if ( max3 > thresb[2]*max0 ) warn("max3 > thresb3*max0 ");
                    if ( rms3 > thresb[3]*rms4 ) warn("rms3 > thresb4*rms4 ");
                    if ( max3 > thresb[4]*max4 ) warn("max3 > thresb5]*max4");
                    if ( avg4 != 0.0 && ABS(avg3) > thresb[5]*ABS(avg4) ) warn("avg4 != 0.0 && ABS(avg3) > thresb6*ABS(avg4) ");
                }

                break;
            }
        }

        // start search forwards
        for (it = itmin; it <= (bad[itr][1]? bad[itr][1] : itmax); ++it) {
            // calc average and rms of the first set of forward sliding two windows
            avg1 = calcAvg(it - itw[1], itw[1], indata[itr]);
            avg2 = calcAvg(it, itw[2], indata[itr]);
            rms1 = calcRms(it - itw[1], itw[1], indata[itr], 0);
            rms2 = calcRms(it, itw[2], indata[itr], 0);
            max1 = calcMax(it - itw[1], itw[1], indata[itr], 0);
            max2 = calcMax(it, itw[2], indata[itr], 0);

            // check if threshold is exceeded
            diff = ABS(indata[itr][it] - indata[itr][it - 1]);
            if (    ABS(indata[itr][it]) > thresa[0]*max0
                 || diff > thresa[1]*max0
                 || rms2 > thresa[2]*rms0
                 || max2 > thresa[2]*max0
               )
            {
                bad[itr][0] = it;  // found front edge of overscaling
                if (verbose > 1) {
                    warn("Front edge at sample %d on %d-th trace: amp=%f diff=%f\n max=%f %f %f rms=%f %f %f Avg=%f %f",
                         it, itr, indata[itr][it], diff, max2, max1, max0, rms2, rms1, rms0, avg2, avg1);
                    if ( ABS(indata[itr][it]) > thresa[0]*max0 ) warn("ABS(indata) > thresa1*max0");
                    if ( diff > thresa[1]*max0 ) warn("diff > thresa2*max0");
                    if ( rms2 > thresa[2]*rms0 ) warn("rms2 > thresa3*rms0");
                    if ( max2 > thresa[2]*max0 ) warn("max2 > thresa3*max0");
                }
                break;
            }
        }

       if (bad[itr][0] | bad[itr][1]) ++nbadtr;

        // repair the trace
        if (  bad[itr][1] && bad[itr][0] && bad[itr][1] < bad[itr][0] ) { // can happen
            warn("  Conflicting sample interval (%d > %d) of overscaling for %d-th trace",
                bad[itr][0], bad[itr][1], itr);
            bad[itr][0] = bad[itr][1] = 0; // reset
        }
        else if (bad[itr][0] > 0 || bad[itr][1] > 0) {
            if ( bad[itr][1] == 0 ) bad[itr][1] = bad[itr][0];

            if ( bad[itr][0] > 0 )  { //edge above found
                // trace back the starting sample of overscaling to the crossing over zero
                diff = ABS(indata[itr][bad[itr][0]]);
                for (it=1; it<=itw[2]; ++it) {
                    if ( ABS(indata[itr][bad[itr][0] - it]) < diff) {
                        diff = ABS(indata[itr][bad[itr][0] - it]);
                        continue;
                    } else {
                        bad[itr][0] -= it;
                        break;
                    }
                }
            } else { // trace back to find first zero crossing 
                //sign = (avg3 > 0.0)? 1.0 : -1.0;
                //for (it=bad[itr][1] - 1; it>itmin; --it) {
                //    if ( sign*indata[itr][it] > 0.0 ) continue;
                //    else break;
                //}
                bad[itr][0] = itmin;
            }

            // trace forward the ending sample of overscaling to minimize diff to potential bias
            diff = ABS(indata[itr][bad[itr][1]] - avgb);
            for (it=1; it<=itw[3]; ++it) {
                if ( ABS(indata[itr][bad[itr][1]+it] - avgb) < diff ) {
                    diff = ABS(indata[itr][bad[itr][1]+it] - avgb);
                    continue;
                } else {
                    bad[itr][1] += it;
                    break;
                }
            }
//            *********/
            // always try to remove bias
            for (it = bad[itr][1] + nsb; it < ns; ++it) {
                if (indata[itr][it] == 0.0 ) continue; //skip hard zero
                avgb = calcAvg(MIN(it, ns - ns/20), ns/20, indata[itr]);
                indata[itr][it] -= avgb;  
            }
            // mute overscaled impulse
            for (it = bad[itr][0]; it <= bad[itr][1]; ++it) {
                indata[itr][it] = 0.0;
            }
            // linearly taper above
            for (it = 0; it < nsa; ++it) {
                indata[itr][bad[itr][0] - it] = indata[itr][bad[itr][0] - nsa] * it / (nsa + 1);
            }
            // linearly taper below
            for (it = 0; it < nsb; ++it) {
                indata[itr][bad[itr][1] + it] = indata[itr][bad[itr][1] + nsb] * it / (nsb + 1);
            }
        }

        if (!(bad[itr][0] | bad[itr][1])) { // debias for rest of traces until first zero crossing
            sign = (calcAvg(ns/2, ns/20, indata[itr]) > 0)? 1.0 : -1.0;
            for (it=ns-1; it > itmin; --it) {
                if (indata[itr][it] == 0.0 ) continue;  //skip hard zero
                if ( sign*indata[itr][it] < 0.0 ) break;
            }
            int start = it;
            for (it=start; it < ns; ++it) {
                avgb = calcAvg(MIN(it, ns - ns/20), ns/20, indata[itr]);
                indata[itr][it] -= avgb;
            }
        }
    }

    return nbadtr;
}

float calcAvg(int itbeg, int nt, float* data)
{
    int it, ns = nt;
    double sum = 0.0;
    for (it = 0; it < nt; ++it )
        if ( data[itbeg + it] == 0.0 ) {// skip hard zero
            --ns;
        } else {
            sum += data[itbeg + it];
        }
    return (ns > 0)? sum/ns : 0.0;
}

float calcRms(int itbeg, int nt, float* data, float base)
{
    double sum = 0.0;
    int nc, it;
    for (it = 0, nc = 0; it < nt; ++it ) {
        if (data[itbeg + it] != 0.0 ) {
            ++nc;
            sum += (data[itbeg + it] - base)*(data[itbeg + it] - base);
        }
    }
    return sqrt(sum/nc);
}

float calcMax(int itbeg, int nt, float* data, float base)
{
    float max = 0.0;
    int it;
    for (it = 0; it < nt; ++it ) {
        if ( ABS(data[itbeg + it] - base) > max)
            max = ABS(data[itbeg + it] - base);
    }

    return max;
}
float calcRMS(int itr, int ntrace, int itbeg, int nt, int** bad, float** data, float base)
{
    double sum = 0.0;
    int it, i = 1, j = 0, nc = 0;
    do {
        //if (bad == NULL || bad[itr - i] == 0) {
            for (it = 0; it < nt; ++it ) {
                if (itr - i >= 0 && data[itr - i][itbeg + it] != 0.0 ) {
                    ++nc;
                    sum += (data[itr - i][itbeg + it] - base)*(data[itr - i][itbeg + it] - base);
                }
            }
            ++j;
        //}
        ++i;
    } while (j < ntrace );
    return sqrt(sum/nc);
}

float calcMAX(int itr, int ntrace, int itbeg, int nt, int** bad, float** data, float base)
{
    float* maxs = ealloc1float(ntrace);
    memset(maxs, 0, ntrace*FSIZE);
    int it, i = 1, j = 0;
    do {
        //if (bad == NULL || bad[itr - i] == 0) { //exclude bad trace
            for (it = 0; it < nt; ++it ) {
                if (itr - i >= 0 && data[itr - i][itbeg + it] != 0.0 ) {
                    if ( ABS(data[itr - i][itbeg + it] - base) > maxs[j])
                        maxs[j] = ABS(data[itr - i][itbeg + it] - base);
                }
            }
            ++j;
        //}
        ++i;
    } while (j < ntrace );

    float medianmax = (ntrace == 1)? maxs[0] : quick_select(maxs, ntrace);
    free1float(maxs);
    return medianmax;
}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

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
