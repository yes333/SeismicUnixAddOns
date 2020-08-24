/* Copyright (c) RWS, Read Well Services, Oslo, 2011.*/
/* All rights reserved.                       */

/* SPFXMEDIAN: $Revision: 1.0 $ ; $Date: 2011/04/09 18:32:28 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"
#include <fftw3.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPFXMEDIAN - Multi-channel FX median filtering for 2/3D gather                ",
"                                                                               ",
" spfxmedian <stdin >stdout  [optional parameters]                              ",
"                                                                               ",
" Optional parameters:                                                          ",
" key=cdpt              gather sorting key                                      ",
" nxmax=tr.sfs          max. number of trace in x direction                     ",
" nymax=tr.sfe          max. number of trace in y direction                     ",
"                                                                               ",
" twin=1.0       [s]    center time subwindow for calculate frequency spectra   ",
" overl=0.1*twin [s]    overlapping window for merging adjacent subwindows      ",
" taper=0.1*twin [s]    taper windows added on both edges of center window      ",
" start=0.0      [s]    start of filter window                                  ",
" end=(ns*dt)    [s]    end of filter window. Default: trace length             ",
"                                                                               ",
" fmin=0.0      [Hz]    lower limit of frequency band to be filtered            ",
" fmax=fNyquist [Hz]    upper limit of frequency band to be filtered            ",
" dfmax=1.0     [Hz]    upper limit of sampling interval in frequency domain    ",
" nwx=11                trace window in x direction to compute median value     ",
" nwy=1                 trace window in y direction                             ",
" threshold=3.0         threshold over which spectral amplitude is replaced     ",
" level=1.0             with level of median value                              ",
"                                                                               ",
" trid=trid             key holding trace id                                    ",
" use=                  trids of traces used to calculate median value. Default: all",
" apply=                trids of traces median filter is applied to. Default: all",
"                                                                               ",
" verbose=0             >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
"                                                                               ",
" The center time window gives the subwindow length by which a trace is divided ",
" between the parameters start= and end= and median filter is performed. To smooth",
" the transition between two adjacent subwindows after filtering, samples within",
" an overlapping window specified by overl= are linearly tapered and summed. As ",
" usaul tapering beyond both edges are made. All together a window of data with ",
" length (taper + twin + overl + taper) is copied into a cache and appended with",
" zero samples to certain total length in term to garantee a sufficiently small ",
" sampling interval in frequency domain (as specified by dfmax=). A larger center",
" window is used to supress long wavelength noises like swell noises while short",
" one for noise bursts with window length corresponding to the typical length of",
" the noise bursts. In general a longer window length is preferred than a short ",
" one if both achieve similar result, not just due to consideration of computational",
" speed.                                                                        ",
" Input gathers are not required to be regularly sampled, but sorted by shotline",
" and numbered in X and Y position using SUSETGATHER with index stored in keywords",
" nhs and nvs. If gathers are not regularly sampled with varying number of traces,",
" nxmax and nymax must be set big enough to hold the largest one. Only traces with",
" non-zero trace id value will be processed and output. In case only certain trids",
" are used for computing median value, trace window nwx/nwy should be enlargered",
" accordingly.                                                                  ",
"                                                                               ",
" Typical flow example looks like:                                              ",
"                                                                               ",
" spdbread select=\"cdpt+|ep+|sx+|duse+\" |\\ ",
" ......  |\\ ",
" susetgather key=cdpt xkey=sx |\\ ",
" spfxmedian key=cdpt nwx=40 nwy=20 trid=duse use=1,3 apply=2,4  |\\  ",
" su...  ",
"                                                                               ",
" assuming here cdpt hold receiver no, ep shotline number                       ",
"                                                                               ",
"                                                                               ",
" Version 1.1.2 last modified Dec. 2011 by Sanyu Ye                             ",
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
void ApplyTaper(int ntaper, int nlive, float* etaper, float* cdata);
int  Conv2AmpPhase(int direction, int nwfft, int ifs, int ife, float** cdata);
int  DoFFT(int fft, fftwf_plan plan, int nwfft, int ntfft, float fftscale, float **cdata);
int GetMedianAmp(const int nxmax, const int nymax, const int nwfft, const int nwx, const int nwy, const int ifs, const int ife,
        int** use4median, int** apply, float**** cdata, float**** medianAmp);
//int  GetMedianAmp(int nxmax, int nymax, int nwfft, int nwx, int nwy, int ifs, int ife,
//        int** use4median, int** apply, float**** cdata, float**** medianAmp);
int  ApplyMedianAmp(int nxmax, int nymax, int nwfft, int ifs, int ife,
        int** apply, float threshold, float level, float**** cdata, float**** medianAmp, int** modified);
void MergeWindows(int nt, int nlive, int nmix, int ntaper, int its, int ite,
        int nwfft, float* mtaper, float** cdata, float* idata, float* odata);
float quick_select(float *arr, int n);

int main(int argc, char **argv)
{
    cwp_String key, key5;       /* header key word from segy.h */
    cwp_String type, type5;     /* type of key	*/
    int index, index5;          /* index of key	*/
    Value val, valnew, vala;    /* value of key				*/
    int nsegy;                  /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int nused, napplied;      /* number of traces actually used and filtered */
    int ngather, ntotal;    /* number of total input traces and gathers */
    int nt;                 /* number of points on trace		*/
    int i, ix, iy, itr, iw; /* counters				*/
    int ncuse, ncapply, *id2use, *id2apply;
    int nwx, nwy;       /* length of trace window/aperture	*/
    int ifs, ife, its, ite;      /* length of time window in samples	*/
    int nlive, nmix, ntaper, ntw, ntfft;      /* length of time window in samples	*/
    int nf;      /* length of time window in samples	*/

    int verbose;
    int nx, ny, nxmax, nymax; /* max dimensions of input/output trace array */

    float dt;       /* time sample interval (sec)	*/
    float thres, level;
    float df, dfmax, fmin, fmax, fNyq;
    float ts, te, tw0, tw1, tw2;

    segy tr, tro;
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* Get parameters  */
    if (!getparint("verbose", &verbose))    verbose = 0;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    nt = tr.ns;
    dt = ((double) tr.dt) / 1000000.0;
    fNyq = 0.5/dt; // Nyquist frequency

    nx = tr.sfs;
    ny = tr.sfe;
    if(!getparint("nxmax", &nxmax)) nxmax = nx;
    if(!getparint("nymax", &nymax)) nymax = ny;
    if(!getparint("nwx", &nwx)) nwx = 11;
    if(!getparint("nwy", &nwy)) nwy = 1;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* get SU keys holding window parameters */
    if (!getparstring("trid", &key5)) key5 = "trid";
    type5 = hdtype(key5);
    index5 = getindex(key5);

    if (!getparfloat("twin",  &tw0)) tw0 = 1.0;
    if (!getparfloat("overl", &tw1)) tw1 = 0.1*tw0;
    if (!getparfloat("taper", &tw2)) tw2 = 0.1*tw0;
    if (!getparfloat("start", &ts))  ts  = 0.0;
    if (!getparfloat("end",   &te))  te  = (nt - 1)*dt;
    if (!getparfloat("fmin",  &fmin)) fmin = 0.0;
    if (!getparfloat("fmax",  &fmax)) fmax = fNyq;
    if (!getparfloat("dfmax",&dfmax)) dfmax = 1.0;
    if (!getparfloat("level",&level)) level = 1.0;
    if (!getparfloat("threshold",&thres)) thres = 3.0;

    // trace id
    ncuse = countparval("use");
    if (ncuse > 0) {
        id2use = ealloc1int(ncuse);
        getparint("use", id2use);
    }
    ncapply = countparval("apply");
    if (ncapply > 0) {
        id2apply = ealloc1int(ncapply);
        getparint("apply", id2apply);
    }

    nlive = NINT(tw0/dt);  // central data samples without tapering
    nmix = NINT(tw1/dt);   // edge samples of central data linearly tapered and summed where adjacent subwindows merges
    ntaper = NINT(tw2/dt); // tapered edge samples added on both side
    ntw = nlive + nmix + 2*ntaper; // number of data samples used for fft
    its = NINT(ts/dt);  // starting sample of filter window in time
    ite = NINT(te/dt);  // ending sample of filter window in time
    if (its < 0 || ite - its < 1 || ite > (nt-1))
        err("  Inconsistent filter window in number of samples (0 =< %d < %d < %d)", its, ite, nt);
    // calculate number of individual live subwindow for fft
    int nwfft = (ite - its + nmix) / nlive;
    if (nwfft * nlive + nmix < ite - its ) ++nwfft;

    // figure out samling interval and number in frequncy domain
    ntfft = NINT(2.0 * ( fNyq / dfmax + 1));
    if (ntw > ntfft) {
        ntfft = 2*(ntw/2 + 1);
    }
    nf = ntfft/2;
    df = fNyq/(ntfft/2 - 1);
    ifs = NINT(fmin/df); // starting sample of filter window in f
    ife = NINT(fmax/df); // end sample of filter window in f

    /* Allocate memory */
    int** trid = ealloc2int(nxmax, nymax);
    int** apply = ealloc2int(nxmax, nymax);
    int** use4median = ealloc2int(nxmax, nymax);
    int** modified = ealloc2int(nxmax, nymax);
    float* mtaper = ealloc1float(nmix);  // linear taper to merge overlapping windows
    float* etaper = ealloc1float(ntaper);  // edge taper on both side of
    float**** medianAmp = ealloc4float(ife - ifs + 1, nwfft, nxmax, nymax);
    float**** cdata = ealloc4float(ntfft, nwfft, nxmax, nymax);
    // allocate memory input/output data traces
    float*** trdata = ealloc3float(nt, nxmax, nymax);
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nxmax, nymax, HDRBYTES);
    /* zero out data memory */
    memset(*trid, 0, nxmax*nymax*ISIZE);
    memset(*apply, 0, nxmax*nymax*ISIZE);
    memset(*modified, 0, nxmax*nymax*ISIZE);
    memset(*use4median, 0, nxmax*nymax*ISIZE);
    memset(mtaper, 0, nmix*FSIZE);
    memset(etaper, 0, ntaper*FSIZE);
    memset(***medianAmp, 0, nymax*nxmax*nwfft*(ife - ifs + 1)*FSIZE);
    memset(***cdata, 0, nymax*nxmax*nwfft*ntfft*FSIZE);
    memset(**trdata, 0, nymax*nxmax*nt*FSIZE);
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

    /* creat fftw plan*/
    fftwf_plan planf = NULL, planb = NULL;
    fftwf_complex* ctr = (fftwf_complex*) &cdata[0][0][0][0];
    int ntt   = ntfft - 1;
    int nhalf = ntfft/2;
    int nreal = ntfft*nwfft;
    int ncmpl = nreal/2;
    planf = fftwf_plan_many_dft_r2c(1, &ntt, nwfft, &cdata[0][0][0][0], NULL, 1, ntfft, ctr, NULL, 1, nhalf, FFTW_MEASURE);
    planb = fftwf_plan_many_dft_c2r(1, &ntt, nwfft, ctr, &ncmpl, 1, nhalf, &cdata[0][0][0][0], &nreal, 1, ntfft, FFTW_MEASURE);
    float fftscale = 1.0/((float) (ntt));  // fft scaling factor

    // calc linear merging taper
    for(i=0; i<nmix; ++ i) mtaper[i] = (i + 1.0)/(nmix + 1.0);  //increasing with index
    // calculate edge taper, Gausian taper
    for (i = 0; i < ntaper; i++) {
        float f = (i + 1.0) / (ntaper + 1.0);
        float x = 2.0 * (1.0 - f);
        etaper[i] = expf(-(x * x)); //increasing with index
    }

    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = ntotal = 0;
    nused = napplied = 0;
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data */
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

            gethval(&tr, index5, &vala);
            trid[iy][ix] = vtoi(type5, vala);

            ++ntr;
        } else { // new gather or END_OF_FILE
            ++ngather;
            ntotal += ntr;
            if (verbose) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, val));

            // load data of time subwindow into caches, apply taper and fft to spectra
            itr = 0;
            for (iy=0; iy<nymax; ++iy) {
                for (ix=0; ix<nxmax; ++ix) {
                    ++itr;
                    if (trid[iy][ix] == 0) continue; // skip zero trace

                    if (ncapply > 0) {
                        for (i=0; i<ncapply; ++i) {
                            if (trid[iy][ix] == id2apply[i]) {
                                ++napplied;
                                apply[iy][ix] = 1;
                            }
                        }
                    } else {
                        ++napplied;
                        apply[iy][ix] = 1;
                    }
                    if (ncuse > 0) {
                        for (i=0; i<ncuse; ++i) {
                            if (trid[iy][ix] == id2use[i]) {
                                ++nused;
                                use4median[iy][ix] = 1;
                            }
                        }
                    } else {
                        ++nused;
                        use4median[iy][ix] = 1;
                    }
                    for (iw=0; iw<nwfft; ++iw) { // loop over subwindows along each trace
                        int itws = its + iw * nlive - ntaper;  // starting sample to copy
                        int isoffset = (itws < 0)? -itws : 0; // if less than start trace, get start offset
                        int itwe = its + (iw + 1) * nlive + nmix + ntaper; // ending sample to copy
                        int ieoffset = (itwe > nt - 1)? itwe - nt : 0; // if larger than trace length, get end offset
                        int nwlength = itwe - itws - isoffset - ieoffset; // length to copy
                        memcpy(&cdata[iy][ix][iw][isoffset], &trdata[iy][ix][MAX(0, itws)], nwlength*FSIZE);
                        ApplyTaper(ntaper, nlive + nmix, etaper, cdata[iy][ix][iw]);
                    }
                    DoFFT(FFTW_FORWARD, planf, nwfft, ntfft, fftscale, cdata[iy][ix]);
                    Conv2AmpPhase(1, nwfft, ifs, ife, cdata[iy][ix]);  // only partially between frequency sample ifs ~ ife
                }
            }

            if (!nused || !napplied) {
                if (!napplied)
                    warn("  Check parameter trid= because no trace to be applied median filter for %d traces of %d-th gather (%s=%d)...",
                            ntr, ngather, key, vtoi(type, val));
                if (!nused)
                    warn("  Check parameter trid= because no trace to be used for median filter for %d traces of %d-th gather (%s=%d)...",
                            ntr, ngather, key, vtoi(type, val));
            }
            
            GetMedianAmp(nxmax, nymax, nwfft, nwx, nwy, ifs, ife, use4median, apply, cdata, medianAmp);
            ApplyMedianAmp(nxmax, nymax, nwfft, ifs, ife, apply, thres, level, cdata, medianAmp, modified);

            // inverse fft, merge and output
            itr = 0;
            for (iy=0; iy<nymax; ++iy) {
                for (ix=0; ix<nxmax; ++ix) {
                    ++itr;
                    if (trid[iy][ix] == 0) continue; // skip zero trace
                    memcpy(&tro, &hdrs2d[iy][ix], HDRBYTES); // copy back header to output trace
                    if (apply[iy][ix] == 0 || modified[iy][ix] == 0) { // forward unfiltered trace
                        memcpy(tro.data, trdata[iy][ix], nt*FSIZE);
                        puttr(&tro);
                        continue;
                    } else {
                        Conv2AmpPhase(-1, nwfft, ifs, ife, cdata[iy][ix]);
                        DoFFT(FFTW_BACKWARD, planb, nwfft, ntfft, 1.0, cdata[iy][ix]);
                        MergeWindows(nt, nlive, nmix, ntaper, its, ite, nwfft, mtaper, cdata[iy][ix], trdata[iy][ix], tro.data);
                        puttr(&tro);
                    }
                }
            }

            /* zero out data memory */
            memset(*trid, 0, nxmax*nymax*ISIZE);
            memset(*apply, 0, nxmax*nymax*ISIZE);
            memset(*modified, 0, nxmax*nymax*ISIZE);
            memset(*use4median, 0, nxmax*nymax*ISIZE);
            memset(***medianAmp, 0, nymax*nxmax*nwfft*(ife - ifs + 1)*FSIZE);
            memset(***cdata, 0, nymax*nxmax*nwfft*ntfft*FSIZE);
            memset(**trdata, 0, nymax*nxmax*nt*FSIZE);
            memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

            val = valnew;
            ntr = nused = napplied = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces of %d gathers are processed", ntotal, ngather);

    return (CWP_Exit());
}

void  ApplyTaper(int ntaper, int nlive, float* etaper, float* cdata)
{
    int it;

    for(it=0; it<ntaper; ++it) {
        cdata[it] *= etaper[it];
        cdata[2*ntaper + nlive - 1 - it] *= etaper[it];
    }
}

int Conv2AmpPhase(int direction, int nwfft, int ifs, int ife, float** cdata)
// only partially between frequency sample ifs ~ ife
{
    int iw, ifr;
    float tmp;

    for (iw=0; iw<nwfft; ++iw) {
        for (ifr=ifs; ifr<=ife; ++ifr) {
            if (direction < 0) { // from amp/phase to real/imaginary
                tmp = cdata[iw][2*ifr] * cosf(cdata[iw][2*ifr + 1]);
                cdata[iw][2*ifr + 1] = cdata[iw][2*ifr] * sinf(cdata[iw][2*ifr + 1]);
                cdata[iw][2*ifr] = tmp;
            } else {
                tmp = sqrtf(cdata[iw][2*ifr]*cdata[iw][2*ifr] + cdata[iw][2*ifr + 1]*cdata[iw][2*ifr + 1]);
                cdata[iw][2*ifr + 1] = atan2(cdata[iw][2*ifr + 1], cdata[iw][2*ifr]);
                cdata[iw][2*ifr] = tmp;
            }
        }
    }
    return ife - ifs + 1;
}

int DoFFT(int fft, fftwf_plan plan, int nwfft, int ntfft, float fftscale, float **cdata)
{
    int iw, it;

    fftwf_complex* ctr = (fftwf_complex*) &cdata[0][0];
    if (fft == FFTW_FORWARD) {
        fftwf_execute_dft_r2c(plan, &cdata[0][0], ctr);
        for (iw=0; iw<nwfft; ++iw)
            for (it=0; it<ntfft; ++it)
               cdata[iw][it] *=  fftscale;
    } else {
        fftwf_execute_dft_c2r(plan, ctr, &cdata[0][0]);
    }

    return fft;
}


int GetMedianAmp(const int nxmax, const int nymax, const int nwfft, const int nwx, const int nwy, const int ifs, const int ife,
        int** use4median, int** apply, float**** cdata, float**** medianAmp)
{
    // get the median amplitude value for frequencies to be filtered
    int ix, iy, iw, ifr, ixx, iyy;

    int nbuf = (nwx%2 ? nwx : nwx + 1)*(nwy%2 ? nwy : nwy + 1);
    float* amps = ealloc1float(nbuf);

    int ntr = 0;
    for (iy=0; iy<nymax; ++iy) {
        for (ix=0; ix<nxmax; ++ix) {
            if (apply[iy][ix] == 0) continue; // skip trace not to be filtered
            ++ntr;
            for (iw=0; iw<nwfft; ++iw) { // loop over windows along each trace
                for (ifr=ifs; ifr<=ife; ++ifr) { // loop over frequency
                    int ix1 = ix - nwx/2;
                    ix1 = MAX(0, ix1);
                    int ix2 = ix + nwx/2;
                    ix2 = MIN(nxmax - 1, ix2);
                    int iy1 = iy - nwy/2;
                    iy1 = MAX(0, iy1);
                    int iy2 = iy + nwy/2;
                    iy2 = MIN(nymax - 1, iy2);
                    memset(amps, 0, nbuf*FSIZE);
                    int nfm = 0;
                    for (iyy=iy1; iyy<=iy2; ++iyy) {
                        for (ixx=ix1; ixx<=ix2; ++ixx) {
                            if (use4median[iyy][ixx] > 0) {
                                amps[nfm] = cdata[iyy][ixx][iw][2 * ifr];
                                ++nfm;
                            }
                        }
                    }
                    medianAmp[iy][ix][iw][ifr - ifs] = quick_select(amps, nfm);
                }
            }
        }
    }
    free1float(amps);
    return ntr;
}

int ApplyMedianAmp(int nxmax, int nymax, int nwfft, int ifs, int ife,
        int** apply, float threshold, float level, float**** cdata, float**** medianAmp, int** modified)
{
    // apply the median amplitude value for frequencies to be filtered
    int ix, iy, iw, ifr;

    int ntr = 0;
    for (iy=0; iy<nymax; ++iy) {
        for (ix=0; ix<nxmax; ++ix) {
            if (apply[iy][ix] == 0) continue; // skip trace not to be filtered
            ++ntr;
            for (iw=0; iw<nwfft; ++iw) { // loop over windows along each trace
                for (ifr=ifs; ifr<=ife; ++ifr) { // loop over frequency
                    if (cdata[iy][ix][iw][2*ifr] > threshold * medianAmp[iy][ix][iw][ifr - ifs]) {
                        cdata[iy][ix][iw][2*ifr] =     level * medianAmp[iy][ix][iw][ifr - ifs];
                        ++modified[iy][ix];
                    }
                }
            }
        }
    }
    return ntr;
}

void MergeWindows(int nt, int nlive, int nmix, int ntaper, int its, int ite,
        int nwfft, float* mtaper, float** cdata, float* idata, float* odata)
{
    int i, iw;

    memset(odata, 0, nt*FSIZE);  // zero out to better debug
    if (its > 0) memcpy(&odata[0], &idata[0], its*FSIZE);  // copy data before filter window
    if (its > 0) { // blend first front edge for first subwindow
        for (i=0; i<nmix; ++i)
            odata[its + i] = (1.0 - mtaper[i])*idata[its + i] + mtaper[i]*cdata[0][ntaper + i];
    } else { // simply copy from first filtered subwindow
        memcpy(&odata[0], &cdata[0][ntaper], nmix*FSIZE);
    }
    memcpy(&odata[its + nmix], &cdata[0][ntaper + nmix], nlive*FSIZE);  // should length (nlive - nmix)

    for (iw=1; iw<nwfft; ++iw) {
        for (i=0; i<nmix; ++i)  // blend front edge of subwindow with end of previous subwindow
            odata[its + iw*nlive + i]
                = (1.0 - mtaper[i])*cdata[iw - 1][ntaper + nlive + i] + mtaper[i]*cdata[iw][ntaper + i];
        // copy middle part of subwindow, overshoot nmix sample in case for last subwindow end exactly at ite
        memcpy(&odata[its + iw*nlive + nmix], &cdata[iw][ntaper + nmix], nlive*FSIZE); // should length (nlive - nmix)
    }

    if (ite < nt) {
        // blend front edge for unfiltered part
        for (i=0; i < nmix && ite + i < nt; ++i)
            odata[ite + i] = (1.0 - mtaper[i])*odata[ite + i] + mtaper[i]*idata[ite + i];
        // copy last subwindow
        if (nt - ite > nmix) memcpy(&odata[ite + nmix], &idata[ite + nmix], (nt - ite - nmix)*FSIZE);
    }
    return;
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

