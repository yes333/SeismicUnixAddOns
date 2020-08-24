/* Copyright (c) READ Well Services, 2009.*/
/* All rights reserved.                       */

/* SPMEDIAN: $Revision: 1.2 $ ; $Date: 2010/04/07 22:58:42 $	*/

#include "su.h"
#include "segy.h"
#include <segyhdr.h>
#include <cwp.h>
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPMEDIAN - MEDIAN filter across multiple traces and samples within a gather 	",
"            with dip scan and repair/interpolation along max. semblence        ",
"                                                                               ",
" spmedian <stdin >stdout [optional parameters]                                 ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
"   key=cdp             gather sorting key                                      ",
"   nxmax=256           max. number of traces expected in a gather              ",
"   stack=0             =1 perform median stack for gather                      ",
"                       in this case following parameters are ignored           ",
"   nm=11               odd number of traces as median filter width             ",
"   nmt=1               odd number of samples as 2nd (temporal) filter dimension",
"   median=1            for median filter                                       ",
"                       =0 for trace mixing                                     ",
"   mix=w1,w2,...       array of length nm for weighting factors                ",
"                       default Gaussian function (see below)                   ",
"                                                                               ",
"   diff=0              =0  output filter                                       ",
"                       =1  output filtered data (input - filter)               ",
"                       =2  output filter and input data (for semblance weighted decon)",
"                                                                               ",
"   exclude=0           exclude bad data of anomalously high amplitude          ",
"   threshold=8.0       ratio to median amplitude for exclusion                 ",
"                                                                               ",
"   scan=0              =1 dip scan for median filter along max. semblance      ",
"   nw=5                number of samples of temporal semblance window          ",
"   nv=51               number of semblance scan of various dips                ",
"   rvmin=-0.5 [ms/m]   min. dip (slowness) semblance scan starts               ",
"   rvmax=0.5  [ms/m]   max. dip (slowness) semblance scan ends                 ",
"                       0.1 s/km ~ a moveout of 100 ms over 1 km                ",
"   repair=0            a non-zero trid set to repair/interpolate corresponding ",
"                       traces with median values of max. dip scan semblence    ",
"                       scan=1 set automatically                                ",
"                                                                               ",
"   tmin=0.0            start time of the window                                ",
"   tmax=(from header)  end time of the window                                  ",
"   itmin=0             min time sample of the window                           ",
    "   itmax=(from header) max time sample of the window			",
"   nt=itmax-itmin+1    number of time samples of the window                    ",
"                                                                               ",
"   wkey=offset         key to select trace window and offset for dip scan      ",
"   scale=tr.scalco     scaling factor, default tr.scalco                       ",
"   min=DBL_MIN         min value for key header word                           ",
"   max=DBL_MAX         max value for key header word                           ",
"                                                                               ",
"   verbose=0           >0  echoes information                                  ",
"                                                                               ",
" Notes:                                                                        ",
" ------                                                                        ",
" Median filtering is a process for suppressing a particular moveout on         ",
" seismic sections. Its advantage over traditional dip filtering is that        ",
" events with arbitrary moveout may be suppressed. Median filtering is          ",
" commonly used in up/down wavefield separation of VSP data.                    ",
"                                                                               ",
" This new version assumes that that particular moveout is already corrected (by",
" using other modules e.g. sushw & sustatic). Events that should be filtered out",
" are already aligned horizontally. Otherwise dip scan may be enabled to find the",
" dip with max. semblence. To avoid aliasing and reduce number of scans (thus   ",
" computing time), data should be moveout corrected that the dipping events to be",
" filtered appear possibly subhorizontal.                                       ",
"                                                                               ",
" Default mixing weighting factors are Gaussian function                        ",
"   W(n) = e**{-[(nh - n)/(nh - 1)]**2} where n=1..nm, nh=(nm+1)/2              ",
" Or input nh numbers like mix=0.2,0.5,0.8,1.0 for nm=7, nh=4                   ",
"                                                                               ",
" The second (temporal) dimension of median filter is used for trace mixing as  ",
"  well. Also the bad data exclusion flag exclude=1 if specified                ",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
" 1. trace mixing                                                               ",
" spmedian median=0 diff=0 nm=5 mix=0.5,0.8,1.0,0.8,0.5 <in.su >out.su          ",
"                                                                               ",
" 2. median stack                                                               ",
" spmedian key=offset stack=1 <in.su > out.su                                   ",
"                                                                               ",
" 3. moveout correction and median filter                                       ",
" su... |\\",
" sushw match=offset key=tstat interp=1 infile=x-t-polygon  |\\",
" sustatic hdrs=1 |\\",
" spmedian key=cdp tmin=0.5 tmax=3 |\\",
" sustatic hdrs=1 sign=-1 |\\",
" spmovie key=cdp",
"                                                                               ",
" 4. linear moveout and filter out direct water wave using dip scan             ",
" su... |\\",
" sureduce v=1.4 mode=1  |\\",
" spmedian tmin=0.5 tmax=1.5 scan=1 min=500 rvmin=-0.1 rvmax=0.1 |\\",
" sureduce v=-1.4 mode=1 |\\",
" su...                                                                         ",
"                                                                               ",
" 5. repair/interpolate traces                                                  ",
" su... |\\",
" sushw match=cdpt key=trid values=3,301,9,301 |\\",
" suchw key1=tstat key2=laga b=0.1 |\\",
" sustatic hdrs=1 |\\",
" spmedian key=cdp repair=301 nm=7 nmt=1 nw=5 nv=51 rvmin=-0.5 rvmax=0.5 |\\",
" sustatic hdrs=1 sign=-1 |\\",
" su,,,                                                                         ",
"                                                                               ",
" Version 1.3  last modified August 2012 by Sanyu Ye                           ",
"                                                                               ",
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

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float quick_select(float** arr2d, int nm, int nmt);
float median_select(float** arr2d, int nm, int nmt, int exclude, float threshold, int twoPass);
int medianstack(int nmt, int ntr, int ns, float** indata, float** outdata, int exclude, float threshold);
int medianfilter(int nm, int nmt, int ntr, int ns, int itmin, int itmax, float* mix, float** indata, float** outdata, int* win, int exclude, float threshold);
int mfws(int nm, int nmt, int ntr, int itmin, int itmax, int nw, int nv, int ns, float dt, float vmin, float vmax,
        double* offset, float** indata, float** outdata, int* win, int exclude, float threshold);
int interpolate(int repair, int* trid, int nm, int nmt, int ntr, int itmin, int itmax, int nw, int nv, int ns, float dt, float vmin, float vmax,
        double* offset, float** indata, float** outdata, int* win, int exclude, float threshold);
float interp(int it0, int itr0, int itr, int ns, float dt, float rv, double* offset, float** indata);

int main(int argc, char **argv) {
    char *key; /* header key word from segy.h		*/
    char *type; /* ... its type				*/
    int index; /* ... its index			*/
    Value val, valnew; /* ... its value			*/

    char *wkey; /* header key word for windowing and offset		*/
    char *wtype; /* ... its type				*/
    int windex; /* ... its index			*/
    Value wval; /* ... its value			*/
    double wmin; /* min value    */
    double wmax; /* max value     */
    double scale; /* (offset) scaling factor		*/
    double *offset = NULL; /* array containng offsets for windowing and dip scan     */
    int *win = NULL; /* array contaning flag if a trace within select window     */

    int it; /* sample counter			*/
    int itr; /* trace counter			*/
    int ns; /* number of time samples 		*/
    int nsegy; /* number of byte read 		*/
    int ntr, total, ngather; /* number of traces and gathers			*/
    int nmix; /* number of mixing weights */
    int nxmax; /* max. number of traces expected in gather */
    int eof = 0; /* END OF File flag			*/

    int nt; /* number of time samples of the time window */
    int itmin; /* smallest time sample (zero-based)    */
    int itmax; /* largest time sample (zero-based)     */
    float tmin; /* minimum time to calculate        	*/
    float tmax; /* maximum time to calculate		*/

    float dt; /* time sampling interval		*/

    int verbose; /* flag for printing information	*/

    int stack; /* flag for median stack instead of filter	*/
    int repair; /* flag for repair/interpolation using dip scan	*/
    int median; /* flag for median filter		*/
    int nm, nmt; /* no. of traces and samples of median filter	*/
    int diff ; /* flag for subtracting filter from input data	*/

    int scan; /* flag for semblance scan	*/
    int nw, nv; /* number of samples and scans          */
    float rvmin, rvmax; /* reversed velocity range of scan	*/

    int exclude; /* flag for excluding bad data	*/
    float threshold; /* amplitude ratio	*/

    segy tr, outtr;

    segyhdr *hdrs = NULL; /* placeholder for incoming trace headers	*/
    float **indata = NULL; /* temporary array			*/
    float **outdata = NULL; /* temporary array			*/
    float *mix = NULL; /* array of mix weighting values	*/

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparint("nxmax", &nxmax)) nxmax = 256;
    if (!getparint("stack", &stack)) stack = 0;
    if (!getparint("median", &median)) median = 1;
    if (!getparint("nm", &nm)) nm = 11;
    if (!(nm % 2)) err("Filter width must be a odd number");
    if (!getparint("nmt", &nmt)) nmt = 1;
    if (!(nmt % 2)) err("  2nd filter width must be a odd number");
    if (!getparint("diff", &diff)) diff  = 0;
    if (!getparint("exclude", &exclude)) exclude = 0;
    if (!getparfloat("threshold", &threshold)) threshold = 8.0;
    if (exclude && verbose) warn("  Exclusion of bad data with amplitude ratio to median over threshold=%3.1f", threshold);

    /* Get dip scan parameters */
    if (!getparint("scan", &scan)) scan = 0;
    if (getparint("repair", &repair)) {
        scan = 1;
    } else {
        repair = 0;
    }
    if (!getparint("nw", &nw)) nw = 5;
    if (!(nw % 2)) err("Semblance window must be a odd number");
    if (!getparint("nv", &nv)) nv = 51;
    if (nv < 5) err("nv=%d must larger than 5", nv);
    if (!getparfloat("rvmin", &rvmin)) rvmin = -0.5;
    if (!getparfloat("rvmax", &rvmax)) rvmax =  0.5;
    if (scan && verbose) warn("  Dip scan:  nv=%d  rvmin=%.3f  rvmax=%.3f ", nv, rvmin, rvmax);
    if (scan && nw < nmt) err("  Dip scan window should be longer than 2nd filter dimension ");

    /* Get mix weighting values values */
    if (!median) { //mixing filter
        mix = ealloc1float(nm);
        int i;
        int nhalf = nm/2;
        nmix = countparval("mix");
        if (nmix == 0) {
            for (i=0; i<nm; ++i) {
                float pow = ((float) (i - nhalf))/((float) nhalf);
                mix[i] = exp(-pow*pow);
            }
        } else if (nmix > nhalf && nmix <= nm) {
            getparfloat("mix", mix);
            for (i=0; i<nm-nmix; ++i) mix[nm - 1 - i] = mix[i];
        } else if (nmix <= nhalf) {
            err("  Number of weighting factors (%d < %d) must exceed half filter width", nmix, nhalf + 1);
        } else {
            err("  Too many weighting factors (%d > %d)", nmix, nm);
        }
        if (verbose) {
            fprintf(stderr, "  %d weighting factors for mixing: ", nm);
            for (i=0; i<nm; ++i) fprintf(stderr, "%4.2f ", mix[i]);
            fprintf(stderr, "\n");
        }
    }

    /* Get info from first trace */
    if (!(nsegy = gettr(&tr))) err("can't read first trace");
    if (!tr.dt) err("dt header field must be set");
    dt = ((double) tr.dt) / 1000000.0;
    ns = (int) tr.ns;

    /* Time gating parameters */
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
    if (verbose) warn("  Time window: tmin=%.3f  tmax=%.3f  itmin=%d itmax=%d nt=%d", tmin, tmax, itmin, itmax, nt);

    /* Check time gating values */
    if (itmin < 0)
        err("itmin=%d should be positive", itmin);
    if (nt > SU_NFLTS)
        err("nt=%d exceeds SU_NFLTS=%d", nt, SU_NFLTS);
    if (itmin > itmax)
        err("itmin=%d, itmax=%d conflict", itmin, itmax);
    if (tr.ns <= itmax)
        err("tr.ns=%d, itmax=%d window cannot extend over the trace length", tr.ns, itmax);

    /* Get windowing key  */
    if (!getparstring("wkey", &wkey)) wkey = "offset";
    wtype = hdtype(wkey);
    windex = getindex(wkey);
    if (!getpardouble("min", &wmin)) wmin = -DBL_MAX;
    if (!getpardouble("max", &wmax)) wmax =  DBL_MAX;
    if (!getpardouble("scale", &scale)) {
        scale = (tr.scalco < 0) ? -1.0 / tr.scalco : (tr.scalco > 0) ? tr.scalco : 1.0;
        //if (!strcmp(wkey, "offset")) {
        //    scale = (tr.scalco < 0) ? -1.0 / tr.scalco : (tr.scalco > 0) ? tr.scalco : 1.0;
        //} else {
        //    scale = 1.0;
        //}
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
    outdata = ealloc2float(ns, nxmax);
    offset = ealloc1double(nxmax);
    win = ealloc1int(nxmax);
    int* trid = ealloc1int(nxmax);
    /* zero out data memory */
    memset(trid, 0, nxmax * ISIZE);
    memset(win, 0, nxmax * ISIZE);
    memset(offset, 0, nxmax * DSIZE);
    memset(hdrs, 0, nxmax * HDRBYTES);
    memset(*indata, 0, nxmax * ns * FSIZE);
    memset(*outdata, 0, nxmax * ns * FSIZE);

    /* Read headers and data while getting a count */
    ngather =  ntr = total = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nxmax - 1) err("\nNumber of traces exceeding nxmax=%d\n", nxmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);
            memcpy(indata[ntr], tr.data, FSIZE * ns);
            gethval(&tr, windex, &wval);
            offset[ntr] = vtod(wtype, wval) * scale;
            trid[ntr] = tr.trid;
            if ( offset[ntr] >= wmin && offset[ntr] <= wmax) win[ntr] = 1;
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (verbose) warn("  processing %d traces %d-th gather (%s=%d)", ntr, ngather, key, vtoi(type, val));

            if (ntr < nm) warn("  Fewer traces (%d < %d) than window length for %d-th gather (%s=%d)", ntr, nm, ngather, key, vtoi(type, val));

            total += ntr;
            // do median filter
            if (stack != 0) {
                medianstack(nmt, ntr, ns, indata, outdata, exclude, threshold);
                memcpy(&outtr, &hdrs[0], HDRBYTES);
                memcpy(outtr.data, outdata[0], FSIZE * ns);
                outtr.nvs = (short) ntr;
                puttr(&outtr);
                if (verbose > 5) warn(" %d traces stacked for %d-th gather (%s=%d)", ntr, ngather, key, vtoi(type, val));

            } else {
                if (repair) interpolate(repair, trid, nm, nmt, ntr, itmin, itmax, nw, nv, ns, dt, rvmin, rvmax, offset, indata, outdata, win, exclude, threshold);
                else if (scan) mfws(nm, nmt, ntr, itmin, itmax, nw, nv, ns, dt, rvmin, rvmax, offset, indata, outdata, win, exclude, threshold);
                else medianfilter(nm, nmt, ntr, ns, itmin, itmax, mix, indata, outdata, win, exclude, threshold);
                //output traces
                for (itr = 0; itr < ntr; ++itr) {
                    memcpy(&outtr, &hdrs[itr], HDRBYTES);
                    if (diff && !repair) {
                        if (diff == 1) for (it = 0; it < ns; ++it) outtr.data[it] = indata[itr][it] - outdata[itr][it];
                        else if (diff == 2) {
                            memcpy(outtr.data, outdata[itr], FSIZE * ns);   // copy filter
                            puttr(&outtr);                                  // output filter
                            memcpy(outtr.data, indata[itr], FSIZE * ns);    // copy input data
                            outtr.trid++;
                        }
                    } else {
                        memcpy(outtr.data, outdata[itr], FSIZE * ns);
                    }
                    puttr(&outtr);
                }
            }
            val = valnew;
            // reset output data
            memset(win, 0, nxmax * ISIZE);
            memset(offset, 0, nxmax * DSIZE);
            memset(hdrs, 0, nxmax * HDRBYTES);
            memset(*indata, 0, nxmax * ns * FSIZE);
            memset(*outdata, 0, nxmax * ns * FSIZE);
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn("  Totally %d gathers with %d traces processed", ngather, total);

    /* Clean up */
    if (hdrs) free1(hdrs);
    if (indata) free2float(indata);
    if (outdata) free2float(outdata);
    if (mix) free1float(mix);
    if (win) free1int(win);

    return (CWP_Exit());
}

int medianstack(int nmt, int ntr, int ns, float** indata, float** outdata, int exclude, float threshold) {
    int i, j, it;

    /* Allocate float array to hold samples to be filtered */
    float **dtemp1 = ealloc2float(ntr, nmt);

    for (it = 0; it < ns; ++it) { // loop over samples
        memset(*dtemp1, 0, ntr * nmt * FSIZE);
        for (j=0; j<nmt; ++j) {
            int itt = it - nmt/2 + j;
            int outside = (itt < 0 || itt > ns - 1)? 1 : 0;
            for (i = 0; i < ntr; ++i) { // loop over traces
                dtemp1[j][i] = (outside)? 0.0 : indata[i][itt];
            }
        }
        outdata[0][it] = median_select(dtemp1, ntr, nmt, exclude, threshold, 1);
    }

    free2float(dtemp1);
    return 0;
}

int medianfilter(int nm, int nmt, int ntr, int ns, int itmin, int itmax, float* mix, float** indata, float** outdata, int* win, int exclude, float threshold) {
    int i, itr, it, j = 0;
    int nhalf = (nm - 1) / 2; //half width of filter

    /* Allocate float array to hold samples to be filtered */
    float **dtemp1 = ealloc2float(nm, nmt);

    for (it = itmin; it <= itmax; ++it) { // loop over samples
        for (itr = 0; itr < ntr; ++itr) { // loop over traces
            // skip if trace out of selection window, however still used for filter of neighboring traces
            if (!win[itr]) continue;
            if (indata[itr][it] == 0.0 ) continue;  // skip hard zero sample, preserve mute
            memset(*dtemp1, 0, nm * nmt * FSIZE);
            for (j=0; j<nmt; ++j) {
                int itt = it - nmt/2 + j;
                int outside = (itt < 0 || itt > ns - 1)? 1 : 0;
                if (itr <= nhalf) { //starting edge
                    for (i = nhalf - itr; i < nm; ++i) {
                        dtemp1[j][i] = (outside)? 0.0 : indata[itr - nhalf + i][itt];
                    }
                } else if (ntr - itr < nhalf + 1) { //ending edge
                    for (i = 0; i < nhalf + ntr - itr; ++i) {
                        dtemp1[j][i] = (outside)? 0.0 : indata[itr - nhalf + i][itt];
                    }
                } else { // fully within panel
                    for (i = 0; i < nm; ++i) {
                        dtemp1[j][i] = (outside)? 0.0 : indata[itr - nhalf + i][itt];
                    }
                }
            }
            if (mix) { //trace mixing
                if (exclude) {  // exclude bad value before mixing
                    median_select(dtemp1, nm, nmt, exclude, threshold, 0);
                }
                float fmix = 0.0;
                for (j=0; j<nmt; ++j) {
                    for (i = 0; i < nm; ++i) {
                        outdata[itr][it] += dtemp1[j][i] * mix[i];
                        if (dtemp1[j][i] != 0.0) fmix += mix[i];
                    }
                }
                /* Divide by all mixing weight on non-hard-zero samples (normalization)*/
                if (fmix) outdata[itr][it] /= fmix;
            } else { //median filter
                outdata[itr][it] = median_select(dtemp1, nm, nmt, exclude, threshold, 1);
            }
        }
    }

    free2float(dtemp1);
    return 0;
}

// median filter with semblance scan

int mfws(int nm, int nmt, int ntr, int itmin, int itmax, int nw, int nv, int ns, float dt, float rvmin, float rvmax,
        double* offset, float** indata, float** outdata, int* win, int exclude, float threshold) {
    int i, itr, it, iv, imax, j, n;
    int mhalf = (nw - 1) / 2; //half width of temporal semblance window
    int nhalf = (nm - 1) / 2; //half width of filter
    float drv, rv, semb, semb1, semb2, maxsemb, stk;

    /* Allocate float array to hold samples to be filtered */
    float ***dtemp1 = ealloc3float(nm, nw, nv);
    float **dtemp2 = ealloc2float(nm, nmt);

    drv = (rvmax - rvmin) / (nv - 1);
    for (it = itmin; it <= itmax; ++it) { // loop over samples
        for (itr = 0; itr < ntr; ++itr) { // loop over traces
            // skip if trace out of selection window, however still used for filter of neighboring traces
            if (!win[itr]) continue;
            memset(**dtemp1, 0, nm * nv * nw * FSIZE);
            memset(*dtemp2, 0, nm * nmt * FSIZE);
            for (iv = 0, rv = rvmin; iv < nv; ++iv, rv += drv) {
                if (itr <= nhalf) { //starting edge
                    for (j = 0; j < nw; ++j) {
                        for (i = nhalf - itr; i < nm; ++i) {
                            dtemp1[iv][j][i] = interp(it + j - mhalf, itr, itr - nhalf + i, ns, dt, rv, offset, indata);
                        }
                    }
                } else if (ntr - itr < nhalf + 1) { //ending edge
                    for (j = 0; j < nw; ++j) {
                        for (i = 0; i < nhalf + ntr - itr; ++i) {
                            dtemp1[iv][j][i] = interp(it + j - mhalf, itr, itr - nhalf + i, ns, dt, rv, offset, indata);
                        }
                    }
                } else { // within gather
                    for (j = 0; j < nw; ++j) {
                        for (i = 0; i < nm; ++i) {
                            dtemp1[iv][j][i] = interp(it + j - mhalf, itr, itr - nhalf + i, ns, dt, rv, offset, indata);
                        }
                    }
                }
            }

            // calculate semblance and find max
            imax = 0;
            maxsemb = 0.0;
            for (iv = 0; iv < nv; ++iv) {
                n = 0;
                semb = semb1 = semb2 = 0.0;
                for (j = 0; j < nw; ++j) {
                    stk = 0.0;
                    for (i = 0; i < nm; ++i) {
                        if (dtemp1[iv][j][i] != 0.0) { //always skip hard zero
                            ++n;
                            stk += dtemp1[iv][j][i];
                            semb2 += dtemp1[iv][j][i] * dtemp1[iv][j][i];
                        }
                    }
                    semb1 += stk*stk;
                }
                if (semb2 > 0.0) semb = nm * semb1 / semb2;
                if (semb > maxsemb) {
                    imax = iv;
                    maxsemb = semb;
                }
            }

            // copy center sample values to buffer for median calculation
            for (j = 0; j < nmt; ++j) {
                for (i = 0; i < nm; ++i) {
                    dtemp2[j][i] = dtemp1[imax][mhalf - nmt/2 + j][i];
                }
            }
            outdata[itr][it] = median_select(dtemp2, nm, nmt, exclude, threshold, 1);
        }
    }

    free3float(dtemp1);
    free2float(dtemp2);
    return 0;
}

// median filter with semblance scan

int interpolate(int repair, int* trid, int nm, int nmt, int ntr, int itmin, int itmax, int nw, int nv, int ns, float dt, float rvmin, float rvmax,
        double* offset, float** indata, float** outdata, int* win, int exclude, float threshold) {
    int i, itr, it, iv, imax, j, n;
    int mhalf = (nw - 1) / 2; //half width of temporal semblance window
    int nhalf = (nm - 1) / 2; //half width of filter
    float drv, rv, semb, semb1, semb2, maxsemb, stk;

    /* Allocate float array to hold samples to be filtered */
    float ***dtemp1 = ealloc3float(nm, nw, nv);
    float **dtemp2 = ealloc2float(nm, nmt);

    drv = (rvmax - rvmin) / (nv - 1);
    for (it = itmin; it <= itmax; ++it) { // loop over samples
        for (itr = 0; itr < ntr; ++itr) { // loop over traces
            // skip if trace id not equale to repair, but copy from input data
            if (trid[itr] != repair) {
                outdata[itr][it] = indata[itr][it];
                continue;
            }
            memset(**dtemp1, 0, nm * nv * nw * FSIZE);
            memset(*dtemp2, 0, nm * nmt * FSIZE);
            for (iv = 0, rv = rvmin; iv < nv; ++iv, rv += drv) {
                if (itr <= nhalf) { //starting edge
                    for (j = 0; j < nw; ++j) {
                        for (i = nhalf - itr; i < nm; ++i) {
                            if (trid[i] == repair) continue;
                            dtemp1[iv][j][i] = interp(it + j - mhalf, itr, itr - nhalf + i, ns, dt, rv, offset, indata);
                        }
                    }
                } else if (ntr - itr < nhalf + 1) { //ending edge
                    for (j = 0; j < nw; ++j) {
                        for (i = 0; i < nhalf + ntr - itr; ++i) {
                            if (trid[i] == repair) continue;
                            dtemp1[iv][j][i] = interp(it + j - mhalf, itr, itr - nhalf + i, ns, dt, rv, offset, indata);
                        }
                    }
                } else { // within gather
                    for (j = 0; j < nw; ++j) {
                        for (i = 0; i < nm; ++i) {
                            if (trid[i] == repair) continue;
                            dtemp1[iv][j][i] = interp(it + j - mhalf, itr, itr - nhalf + i, ns, dt, rv, offset, indata);
                        }
                    }
                }
            }

            // calculate semblance and find max
            imax = 0;
            maxsemb = 0.0;
            for (iv = 0; iv < nv; ++iv) {
                n = 0;
                semb = semb1 = semb2 = 0.0;
                for (j = 0; j < nw; ++j) {
                    stk = 0.0;
                    for (i = 0; i < nm; ++i) {
                        if (dtemp1[iv][j][i] != 0.0) { //always skip hard zero
                            ++n;
                            stk += dtemp1[iv][j][i];
                            semb2 += dtemp1[iv][j][i] * dtemp1[iv][j][i];
                        }
                    }
                    semb1 += stk*stk;
                }
                if (semb2 > 0.0) semb = nm * semb1 / semb2;
                if (semb > maxsemb) {
                    imax = iv;
                    maxsemb = semb;
                }
            }

            // copy center sample values to buffer for median calculation
            for (j = 0; j < nmt; ++j) {
                for (i = 0; i < nm; ++i) {
                    dtemp2[j][i] = dtemp1[imax][mhalf - nmt/2 + j][i];
                }
            }
            outdata[itr][it] = median_select(dtemp2, nm, nmt, exclude, threshold, 1);
        }
    }

    free3float(dtemp1);
    free2float(dtemp2);
    return 0;
}

// return linearly interpolated sample value for trace itr at that sample
// which is extrapolated with rv from center trace itr0 and sample it0

float interp(int it0, int itr0, int itr, int ns, float dt, float rv, double* offset, float** indata)
{
    float dd = offset[itr] - offset[itr0];
    float itf = 0.001* dd * rv / dt + it0;  // distance in km
    int it = (int) itf;
    if (it < 0 || it > ns - 2) return 0.0; //outside trace data
    return indata[itr][it] + (itf - it)*(indata[itr][it + 1] - indata[itr][it]);
}

/*
 *  two pass median filter to enable bad data exclusion
 */
float median_select(float** arr2d, int nm, int nmt, int exclude, float threshold, int twoPass)
{
    int i, j;

    float median = quick_select(arr2d, nm, nmt);
    if (!exclude) return median;

    for (j=0; j<nmt; ++j) {
        for (i=0; i<nm; ++i) {
            if (ABS(arr2d[j][i]) > threshold*ABS(median)) arr2d[j][i] = 0.0;
        }
    }
    return (twoPass)? quick_select(arr2d, nm, nmt) : median;
}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

float quick_select(float **arr2d, int nm, int nmt)
{
    int low, high;
    int median;
    int middle, ll, hh, i, j;

    float* arr = *arr2d;  //reuse input 2D array

    // copy data to 1d array while skipping always hard zero
    int n = 0;
    for (j=0; j<nmt; ++j) {
        for (i=0; i<nm; ++i) {
            if (arr2d[j][i] != 0.0 ) arr[n++] = arr2d[j][i];
        }
    }

    low = 0;
    high = n - 1;
    median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median];
        if (high == low + 1) { /* Two elements only */
            return 0.5*(arr[low] + arr[high]);
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high]) ELEM_SWAP(arr[middle], arr[high]);
        if (arr[low] > arr[high]) ELEM_SWAP(arr[low], arr[high]);
        if (arr[middle] > arr[low]) ELEM_SWAP(arr[middle], arr[low]);

        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low + 1]);

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (arr[low] > arr[ll]);
            do hh--; while (arr[hh] > arr[low]);

            if (hh < ll)
                break;

            ELEM_SWAP(arr[ll], arr[hh]);
        }

        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]);

        /* Re-set active partition */
        if (hh <= median) low = ll;
        if (hh >= median) high = hh - 1;
    }

    return (n%2)? arr[median] : 0.5*(arr[median] + arr[median + 1]);
}

#undef ELEM_SWAP

