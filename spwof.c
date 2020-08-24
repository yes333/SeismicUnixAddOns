/* Copyright (c) Read Well Services, 2011.*/
/* All rights reserved.                       */

/* SUWOF: $Revision: 1.0.0 $ ; $Date: 2011/04/09  $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPWOF - Wiener Optimum Filter for wavelet shaping/deterministic deconvolution ",
"                                                                               ",
" spwof <stdin >stdout  [optional parameters]                                  ",
"                                                                               ",
" Optional parameters:                                                          ",
" key=cdpt              gather sorting key                                      ",
" output=6              =6 output filtered traces                               ",
"                       =5 output normalized Wiener filters                     ",
"                       =4 output Wiener filters                                ",
"                       =3 output cross correlation                             ",
"                       =2 output weighted autocorrelations                     ",
"                       =1 output autocorrelations                              ",
" input=1               =1 input contain only (downgoing) traces for filter design",
"                          this option is used for QC purpose                   ",
"                       =2 designed filter applied to next trace (upgoing)      ",
"                          requiring down- and upgoing input trace by trace     ",
"                                                                               ",
" nxmax=tr.sfs          max. number of traces in x direction                    ",
" nymax=tr.sfe          max. number of traces in y direction                    ",
"                                                                               ",
" wavelet=              desired output wavelet. If not given, a spike at delay= ",
"                                                                               ",
" delay=0        [s]    delay time of max. energy of desired output wavelet     ",
"                         relative to that of autocorrelation window            ",
" fl=0.6         [s]    filter length of the causal/memory part                 ",
" al=0.2         [s]    anticipation length. if >0 noncausal, two-sided filter  ",
"                         total length of wiener filter ( al + fl )             ",
" mincorr=stas          key holding start of autocorrelation window in milisecond",
" maxcorr=stae          key holding end of autocorrelation window               ",
" nw=1,1                trace window in x/y direction for mixing autocorrelation",
" mix=1                 =1 boxcar; =2 Gausian; 3 cosine                         ",
"			  weighting function for averaging autocorrelations	",
" pnoise=0.1     [%]    percentage of relative additive noise level             ",
"                                                                               ",
" trid=trid             key holding trace id                                    ",
" use=                  trids of traces used to design filter. Default: all     ",
" apply=                trids of traces shaping filter is applied to. Default: all",
"                                                                               ",
" verbose=0             >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
"                                                                               ",
" This is a general implementation of Wiener Optimum Filter by solving the standard",
" equation R*f=g, where R is Toerplitz matrix of autocorrelation, f Wiener filter",
" vector and g cross correlation vector. It has a wide range of applications like",
" wavelet shaping/matching, designature filter, determininistic deconvolution etc.",
" While total length of Wiener filter is set to be constant by anticipation (al=)",
" and causal window (fl=), both start and end times of autocorrelation window are",
" fetched from given keywords, their values can be variable as set/loaded using ",
" SUSHW. To achieve best performance of Wiener optimum filter, the start time of",
" autocorrelation should be set in constant lag/delay to the time where the max.",
" energy/amplitude appears in input data. The same delay should be used to place",
" max. amplitude of desired output wavelet. This in turn means that the maximum ",
" energy/amplitude of autocorrelation windows and desired output wavelet should ",
" align when placed trace by trace. If the wavelet has its max. amplitude at a  ",
" different time (usually earlier/shorter because of short wavelet) than that of",
" autocorrelation window, this time difference can be compensated by parameter  ",
" (delay=).                                                                     ",
" For deterministic decon (design filter by downgoing, apply to upgoing), filter",
" length (fl=) should be just as long as enough to include multiples you want to",
" collaps. Anticipation length (al=) should be set short, usually 0.1 ~ 0.2 s.  ",
" In general shorter filter lengths (for both al and fl) are preferred over longer",
" ones. The window of autocorrelation should be set in such way that the autocor-",
" relation best reflects the characteristics of the wavetrain to be compresssed.",
" Input gathers are not required to be regularly sampled, but sorted by shotline",
" and numbered in X and Y position using SUSETGATHER with index stored in keywords",
" nhs and nvs. To improve S/N ratio, autocorrelations of adjacent traces are    ",
" weighted and summed together. The trace window in X/Y is defined by parameter ",
" nw=nwx,nwy.                                           ",
" If gathers have varying geometry with different number of traces, nxmax and ",
" nymax must be set big enough to hold the largest one. Only traces with non-zero",
" positive trace id value will be processed and output. If autocorrelation of a ",
" trace becomes zero, in case the correlation window is muted, the trace id is  ",
" negated and only processed if trace/autocorrlation mixing is enabled.         ",
"                                                                               ",
" Typical flow example for determinsitic decon looks like:                      ",
"                                                                               ",
" 1. create a desired wavelet, adapt bandpass filter to spectra of input data   ",
"",
"suspike nt=201 ntr=1 dt=0.002 offset=0 nspk=1 ix1=1 it1=101 |\\",
"subfilt fstoplo=3 fpasslo=6 fpasshi=60 fstophi=90 |\\",
"suwind tmin=0.1 tmax=0.3 |\\",
"tee wavelet.su |\\",
"suxwigb grid1=solid windowtitle=\"Desired output wavelet, bandpass filtered spike\" &",
"",
" 2. test various control parameters and QC results",
"",
"out=6; length=0.6; anti=0.1; lag=0.0;",
"",
"susetgather key=gdel nmax=1000 < dngp-1-ali.su |\\",
"sushw key=stas,stae a=0,1000 |\\",
"spwof input=1 output=$out al=$anti fl=$length delay=$lag wavelet=wavelet.su  |\\",
"suximage perc=99 grid1=solid windowtitle=\"$out anti=$anti length=$length delay=$lag\" &",
"",
" 3. Interleave upgoing data and apply decon with final parameters",
"",
"sumerge files=dngp-1-ali.su,upgpp-1-ali.su |\\",
"susetgather key=gdel nmax=2000 |\\",
"sushw key=stas,stae a=0,1000 |\\",
"spwof input=2 output=6 al=$anti fl=$length  delay=$lag wavelet=wavelet.su |\\",
"tee upgpp-decon.su |\\",
"suximage perc=99 grid1=solid windowtitle=\"decon upgoing data\" &",
"",
"                                                                               ",
" Version 1.0.0 last modified March. 2011 by Sanyu Ye                           ",
"                                                                               ",
 NULL};

/* Credits:
 *      RWS: Sanyu Ye, Feb. 2011
 *
 */
/**************** end self doc *******************************************/


// forward declare prototype

#include "wienerhlp.c"

int main(int argc, char **argv)
{
    cwp_String key, key3, key4, key5;     /* header key word from segy.h */
    cwp_String type, type3, type4, type5;    /* type of key	*/
    int index, index3, index4, index5;     /* index of key	*/
    Value val, valnew, vala;  /* value of key				*/
    int nsegy;              /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;    /* number of total input traces and gathers */
    int nt;                 /* number of points on trace		*/
    int i, ix, iy, itr, it; /* counters				*/
    int input, output, mix, lib;
    int ncorr;      /* length of corr window in samples	*/

    int verbose;
    int nc, nw[2];   /* mixing/weighting windows			*/
    int nx, ny, nxmax, nymax; /* max dimensions of input/output trace array */
    int *id2use, *id2apply;

    float dt;       /* time sample interval (sec)	*/
    float pnoise;   /* additive noise level		*/
    float mincorr, maxcorr, samplims;
    float delay, lwiener, lanti, rmax;
    segy tr, tro;
    FILE* fp = NULL;           /* file pointer for input wavelet    */
    cwp_String  filename;   /* filenames of input data */
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* Get parameters  */
    if (!getparint("verbose", &verbose))    verbose = 0;
    if (!getparint("lib", &lib))            lib = 0;
    if (!getparint("mix", &mix))            mix = 1;
    if (!getparint("output", &output))      output = 6;
    if (!getparint("input", &input))         input = 1;
    if (input < 1 || input > 2) err(" Invalid parameter input=%d for input trace pattern", input);
    if (!getparfloat("pnoise", &pnoise))    pnoise = 0.1;
    if (!getparfloat("rmax", &rmax))    rmax = 0.1;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    nt = tr.ns;
    dt = ((double) tr.dt) / 1000000.0;
    samplims = 1000.0*dt; // sampling interval in ms

    nx = tr.sfs;
    ny = tr.sfe;
    if(!getparint("nxmax", &nxmax)) nxmax = nx;
    if(!getparint("nymax", &nymax)) nymax = ny;
    if(input == 2) nxmax /= 2;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* get SU keys holding window parameters */
    if (!getparstring("mincorr", &key3)) key3 = "stas";
    type3 = hdtype(key3);
    index3 = getindex(key3);
    if (!getparstring("maxcorr", &key4)) key4 = "stae";
    type4 = hdtype(key4);
    index4 = getindex(key4);
    if (!getparstring("trid", &key5)) key5 = "trid";
    type5 = hdtype(key5);
    index5 = getindex(key5);

    // trace id
    int ncuse = countparval("use");
    if (ncuse > 0) {
        id2use = ealloc1int(ncuse);
        getparint("use", id2use);
    }
    int ncapply = countparval("apply");
    if (ncapply > 0) {
        id2apply = ealloc1int(ncapply);
        getparint("apply", id2apply);
    }
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

    if (!getparfloat("delay", &delay))  delay = 0.0;
    if (!getparfloat("fl", &lwiener))   lwiener = 0.6;
    if (!getparfloat("al", &lanti))       lanti = 0.2;
    int idelay = NINT(delay/dt);            // delay in samples of max. energy/spike
    int nwiener = NINT(lwiener/dt) + 1; // length in samples of wiener filter
    int ixcorr = (lanti > 0.0)? NINT(lanti/dt) : 0;  // starting index/offset for cross correlation in case of non-causal filter
    if (ixcorr > 0)  nwiener += ixcorr; // nwiener extended

    /* Allocate memory */
    int** trid = ealloc2int(nxmax, nymax);
    int** apply = ealloc2int(nxmax, nymax);
    int** use4acor = ealloc2int(nxmax, nymax);
    int** imincorr = ealloc2int(nxmax, nymax);
    int** imaxcorr = ealloc2int(nxmax, nymax);
    float* wavelet = ealloc1float(nt);
    float* AvgAC = ealloc1float(nt);   // averaged/mixed acorr
    float* xcorr = ealloc1float(nt);
    float* wiener = ealloc1float(nt);
    float* spiker = ealloc1float(nt);
    float** weight = ealloc2float(nw[0], nw[1]);
    float** ac0 = ealloc2float(nxmax, nymax);
    float*** acorr = ealloc3float(nwiener, nxmax, nymax);
    // allocate memory input/output data traces
    float*** trdata = ealloc3float(nt, nxmax, nymax);
    float*** tr2app = (input == 2)? ealloc3float(nt, nxmax, nymax) : NULL;
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nxmax, nymax, HDRBYTES);
    /* zero out data memory */
    memset(*trid, 0, nxmax*nymax*ISIZE);
    memset(*apply, 0, nxmax*nymax*ISIZE);
    memset(*use4acor, 0, nxmax*nymax*ISIZE);
    memset(*imincorr, 0, nxmax*nymax*ISIZE);
    memset(*imaxcorr, 0, nxmax*nymax*ISIZE);
    memset(wavelet, 0, nt*FSIZE);
    memset(xcorr, 0, nt*FSIZE);
    memset(*weight, 0, nw[0]*nw[1]*FSIZE);
    memset(*ac0, 0, nymax*nxmax*FSIZE);
    memset(**acorr, 0, nymax*nxmax*nwiener*FSIZE);
    memset(**trdata, 0, nymax*nxmax*nt*FSIZE);
    if (tr2app) memset(**tr2app, 0, nymax*nxmax*nt*FSIZE);
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

    int ns = 0;  // input wavelet sample number
    if (getparstring("wavelet", &filename)) { // if wavelet file provided
        fp = efopen(filename, "r");
        if(fgettr(fp, &tro) < HDRBYTES) {
            err(" Failed to read input wavelet file %s", filename);
        };
        ns = tro.ns;
        if (ns > nt) err(" wavelet longer than data trace (%d > %d)", ns, nt);
        if (tro.dt != tr.dt) err(" wavelet must have the same sampling interval as data trace (tr.dt = %d != %d", tro.dt, tr.dt);
        memcpy(wavelet, tro.data, ns*FSIZE);
        efclose(fp);
    } else {
        ns = 1;
        wavelet[0] = 1;
    }

    // calculate weighting matrix
    CalcMix(mix, nw[0], nw[1], weight);

    float wac = 1;  // wavelet maximum amplitude
    if (ns > 1) for(wac=0.0, it=0; it<ns; ++it) wac = MAX(wac, ABS(wavelet[it]));
    
    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = ntotal = 0;
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (nymax == 1) {
                iy = 0;
            } else {
                if (tr.nvs < 1) err("key nvs must set positive and sequentially ascending!");
                iy = tr.nvs - 1;
            }
            if (nxmax == 1) {
                ix = 0;
            } else {
                if (tr.nhs < 1) err("key nhs must set positive and sequentially ascending!");
                ix = (input == 2)? (tr.nhs + 1)/input - 1 : tr.nhs - 1;
            }
            if (ntr/input > nxmax*nymax - 1) err("  array dimension too small (%d < %d) traces input", nxmax*nymax, ntr+1);
            if (ix > nxmax - 1) err("  array dimension nxmax too small (%d < %d) traces input", nxmax, ix+1);
            if (iy > nymax - 1) err("  array dimension nymax too small (%d < %d) traces input", nymax, iy+1);
            if ((ntr + 1)%input == 0) memcpy(&hdrs2d[iy][ix], &tr, HDRBYTES);  // if upgoing trace also present, only copy its header
            memcpy((ntr%input)? tr2app[iy][ix] : trdata[iy][ix], tr.data, FSIZE*tr.ns);

            gethval(&tr, index3, &vala);
            mincorr = vtof(type3, vala);
            gethval(&tr, index4, &vala);
            maxcorr = vtof(type4, vala);
            gethval(&tr, index5, &vala);
            trid[iy][ix] = vtoi(type5, vala);
            imincorr[iy][ix] = NINT(0.001*mincorr/dt);
            imaxcorr[iy][ix] = NINT(0.001*maxcorr/dt);
            if (mincorr < 0 || maxcorr - mincorr < samplims || maxcorr >= (nt-1)*samplims)
                err("  Inconsistent autocorrelation window (0 < %1.0f < %1.0f <= %1.0f) for %d of %d-th gather (%s=%d)",
                        mincorr, maxcorr, (nt-1)*samplims, ntr+1, ngather+1, key, vtoi(type, val));
            if (maxcorr - mincorr < ns*dt)
                warn("Warning: Autocorrelation window shorter than input wavelet (%1.0f - %1.0f < %1.0f) for %d of %d-th gather (%s=%d)",
                        maxcorr, mincorr, 1000*ns*dt, ntr+1, ngather+1, key, vtoi(type, val));
            //if (trid[iy][ix] < 0)
            //    warn("  Negative trace id=%d for %d of %d-th gather (%s=%d)",
            //            trid[iy][ix], ntr, ngather, key, vtoi(type, val));

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
                    if (ncapply > 0) {
                        for (i=0; i<ncapply; ++i) if (trid[iy][ix] == id2apply[i])    apply[iy][ix] = 1;
                    } else {
                        apply[iy][ix] = 1;
                    }
                    if (ncuse > 0) {
                        for (i=0; i<ncuse; ++i)   if (trid[iy][ix] == id2use[i]) use4acor[iy][ix] = 1;
                    } else {
                        use4acor[iy][ix] = 1;
                    }

                    if (use4acor[iy][ix] == 0) continue;
                    // compute autocorrelation
                    ncorr = imaxcorr[iy][ix] - imincorr[iy][ix] + 1;
                    xcor(ncorr, 0, &trdata[iy][ix][imincorr[iy][ix]],
                         ncorr, 0, &trdata[iy][ix][imincorr[iy][ix]],
                         nwiener, 0, acorr[iy][ix]);
                    ac0[iy][ix] = acorr[iy][ix][0];  // cache zero lag ac
                    if (ac0[iy][ix] == 0.0) {
                        trid[iy][ix] = -trid[iy][ix];  // negate id for later identification
                        if (verbose > 3) warn("  Zero autocorrelation for trace %d of %d-th gather (%s=%d)",
                                itr, ngather, key, vtoi(type, val));
                    } //else if (nw[0] > 1 || nw[1] > 1) { // normalize when averaging
                    //    for (it=0; it<nwiener; ++it) acorr[iy][ix][it] /= ac0[iy][ix];
                    //}
                }
            }

            itr = 0;
            for (iy=0; iy<nymax; ++iy) {
                for (ix=0; ix<nxmax; ++ix) {
                    if (trid[iy][ix] == 0) continue; // skip zero trace
                    memcpy(&tro, &hdrs2d[iy][ix], HDRBYTES);  // copy trace head
                    if (output == 1) { // output auto correlation
                        if (use4acor[iy][ix] == 0) continue;
                        tro.ns = nwiener;
                        memcpy(tro.data, acorr[iy][ix], nwiener*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    //float norm = (nw[0] > 1 || nw[1] > 1)? 1.0/(sqrtf(ac0[iy][ix])*wac) : sqrtf(ac0[iy][ix]) / wac;  // mormalizing factor
                    float norm = sqrtf(ac0[iy][ix]) / wac;  // mormalizing factor

                    memset(AvgAC, 0, nwiener*FSIZE);
                    // compute avaraging autocorrelation
                    if (nw[0] > 1 || nw[1] > 1) {
                        CalcAvgAC(nxmax, nymax, use4acor, acorr, nw[0], nw[1], weight, ix, iy, nwiener, AvgAC);
                    } else {
                        memcpy(AvgAC, acorr[iy][ix], nwiener*FSIZE);
                        //if (AvgAC[0] != 0.0) for(it = nwiener - 1; it>=0; --it) AvgAC[it] /= AvgAC[0];  // normalize
                    }

                    if (output == 2) { // output weighted and normalized auto correlation
                        if (use4acor[iy][ix] == 0) continue;
                        tro.ns = nwiener;
                        memcpy(tro.data, AvgAC, nwiener*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    // compute cross correlation
                    memset(xcorr, 0, nt*FSIZE);
                    xcor(ncorr, -ixcorr, &trdata[iy][ix][imincorr[iy][ix]], ns, idelay, wavelet, nwiener, 0, xcorr);

                    if (output == 3) { // output cross correlation
                        if (use4acor[iy][ix] == 0) continue;
                        tro.ns = nwiener;
                        memcpy(tro.data, xcorr, nwiener*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    /* zero out filter vectors */
                    memset(wiener, 0, nt*FSIZE);
                    memset(spiker, 0, nt*FSIZE);
                    /* Get inverse filter by Wiener-Levinson */
                    if (AvgAC[0]) {
                        AvgAC[0] *= 1.0 + 0.01*pnoise; /* pre-whitening */
                        if (verbose == 33) {
                            warn("iy=%d ix=%d norm=%9.5g AC0=%9.5g XC0=%9.5g XCN=%9.5g WMA=%9.5g r=%9.5g", 
                              iy, ix, norm, ac0[iy][ix], xcorr[0], xcorr[ixcorr], wac, xcorr[0]/ac0[iy][ix]);
                        }
                        if ( ABS(xcorr[0]/ac0[iy][ix]/wac) > rmax ) {  // potentially instable
                            AvgAC[0] = ABS(xcorr[0]/wac)/rmax;
                            float r = AvgAC[0]/ac0[iy][ix];
                            for (it=1; it<nwiener; ++it) AvgAC[it] *= r;
                            if (verbose == 33) warn("  added white noise level %3.1f fold", r);
                        }
                        stoepf(nwiener, AvgAC, xcorr, wiener, spiker);
                    }

                    if (output == 4) { // output wiener filter
                        if (use4acor[iy][ix] == 0) continue;
                        tro.ns = nwiener;
                        memcpy(tro.data, wiener, nwiener*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    // scale/normalize filter
                    for(it=0; it<nwiener; ++it) wiener[it] *= norm;
                    if (output == 5) { // output wiener filter
                        if (use4acor[iy][ix] == 0) continue;
                        tro.ns = nwiener;
                        memcpy(tro.data, wiener, nwiener*FSIZE);
                        puttr(&tro);
                        continue;
                    }

                    if (apply[iy][ix] == 0) continue;
                    /* Convolve filter with trace  */
                    convolve_cwp(nwiener, -ixcorr, wiener, nt, 0, (input == 2)? tr2app[iy][ix] : trdata[iy][ix], nt, 0, tro.data);
                    puttr(&tro);
                }
            }

            /* zero out data memory */
            memset(*trid, 0, nxmax*nymax*ISIZE);
            memset(*use4acor, 0, nxmax*nymax*ISIZE);
            memset(*apply, 0, nxmax*nymax*ISIZE);
            memset(*imincorr, 0, nxmax*nymax*ISIZE);
            memset(*imaxcorr, 0, nxmax*nymax*ISIZE);
            memset(*ac0, 0, nymax*nxmax*FSIZE);
            memset(**acorr, 0, nymax*nxmax*nwiener*FSIZE);
            memset(**trdata, 0, nymax*nxmax*nt*FSIZE);
            if (tr2app) memset(**tr2app, 0, nymax*nxmax*nt*FSIZE);
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

