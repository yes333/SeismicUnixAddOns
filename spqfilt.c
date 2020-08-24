/* Copyright (c) RWS, Read Well Services, Oslo, 2011.*/
/* All rights reserved.                       */

/* SPQFILT: $Revision: 1.0 $ ; $Date: 2011/04/09 18:32:28 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <fftw3.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPQFILT - forward and inverse Q FILTering                                     ",
"                                                                               ",
" spqfilt <stdin >stdout  [optional parameters]                                 ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
"   opt=0               =0 inverse Q filtering                                  ",
"                       =1 forward Q filtering                                  ",
"   mode=2              =0 phase only                                           ",
"                       =1 amplitude only                                       ",
"                       =2 phase + amplitude                                    ",
"                                                                               ",
"   tmin=0.0       [s]  start of time window for Q filtering                    ",
"   tmax=(header)  [s]  end of time window for Q filtering                      ",
"                                                                               ",
"   tref=tmin      [s]  reference zero time                                     ",
"   twin=0.2       [s]  time step Q filter applied successively                 ",
"   taper=0.1*twin [s]  linear taper window for merging two adjacent windows    ",
"                                                                               ",
"   Q=100               Q factor                                                ",
"                                                                               ",
"   f0=50         [Hz]  the dominant frequency                                  ",
"                                                                               ",
"   limit=1.0           max. value of exponent for amplitude upscaling          ",
"                                                                               ",
"   verbose=0           >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
" Inverse Q filter is expressed by following formula                            ",
"                                                                               ",
"      U(T + dT, f) = U(T, f) * exp(j*2*PI*f*dT*rr) * exp(rr*PI*f*dT/Q)         ",
"                                                                               ",
"          where rr = (f/f0)^(-r),   r=1/(PI*Q),   j imaginary unit             ",
"                                                                               ",
" first exponent term describes phase correction, while second amplitude upscaling",
"                                                                               ",
" for more background info and discussion, see                                  ",
"                                                                               ",
" Wang Y, 2002: A stable and efficient approach of inverse Q filtering, GEOPHYSICS",
"                                                                               ",
" Version 1.0.0 last modified April. 2011 by Sanyu Ye                           ",
"                                                                               ",
 NULL};

/* Credits:
 *      RWS: Sanyu Ye, April. 2011
 *
 *
 * Notes:
 *
 */
/**************** end self doc *******************************************/

// forward declare prototype
float Qfilter(int opt, int mode, float limit, float Q, float T, float f0, int nf, float df,
        complex* ctr, complex* cfilt);
void MergeTraces(int nt, int nlive, int nmix, int nw, int its, int ite,
        float* mtaper, float** cdata, float* idata, float* odata);

int main(int argc, char **argv)
{
    int nsegy;              /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int itmin, itmax, itref, ntwin, ntaper, nfilt;
    int mode, opt;
    int verbose;
    int i, it;

    float tref, tmin, tmax, twin, taper;
    float Q, f0, T, limit, a;

    segy tr;
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* Get parameters  */
    if (!getparint("verbose", &verbose))    verbose = 0;
    if (!getparint("mode", &mode))          mode = 2;
    if (!getparint("opt", &opt))             opt = 0;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    int ns = tr.ns;   // number of samples
    float dt =  (float) tr.dt / 1000000.0;  // sampling interval in s


    /* get parameters */
    if (!getparfloat("Q", &Q)) Q = 100.0;
    if (!getparfloat("f0", &f0)) f0 = 50.0;
    if (!getparfloat("limit", &limit)) limit = 1.0;
    if (!getparfloat("tmin", &tmin)) tmin = 0.0;
    if (!getparfloat("tmax", &tmax)) tmax = (ns - 1)*dt;
    if (!getparfloat("tref", &tref)) tref = tmin;
    if (!getparfloat("twin", &twin)) twin = 0.2;
    if (!getparfloat("taper", &taper)) taper = 0.1*twin;
    
    itmin = MAX(0,  NINT(tmin/dt));
    itmax = MIN(ns, NINT(tmax/dt));
    itref = NINT(tref/dt);
    ntwin = NINT(twin/dt);
    nfilt = itmax - itmin + 1;  // number of samples of filter window
    ntaper =  NINT(taper/dt);

    if ( (itref - itmin) <= ntwin/2  )
        err(" Reference zero time tref=%5.3f more than half window length (twin=%5.3f) over tmin=%5.3f", tref, twin, tmin);

    int nstep = NINT((tmax - tmin + 0.5*twin)/twin);     // number of step
    int ntfft = 2*(nfilt/2 + 1);

    int nf = ntfft/2;
    float fNyq = 0.5/dt; // Nyquist frequency
    float df = fNyq/(ntfft/2 - 1);

    /* Allocate memory */
    float* mtaper = ealloc1float(ntaper);  // linear taper to merge overlapping windows
    float* data = ealloc1float(ntfft);          // buffer to hold data for fft
    float* odata = ealloc1float(ns);            // buffer for output data
    float** trcs = ealloc2float(ntfft, nstep);     // buffer for filtered data in various steps
    complex* cbuf = ealloc1complex(ntfft/2);
    memset(data, 0, ntfft*FSIZE);
    memset(*trcs, 0, nstep*ntfft*FSIZE);

    /* creat fftw plan*/
    fftwf_plan planf = NULL, planb = NULL;
    fftwf_complex* ctr = (fftwf_complex*) data;
    int ntt   = ntfft - 1;
    planf = fftwf_plan_dft_r2c_1d(ntt, data, ctr, FFTW_MEASURE);
    planb = fftwf_plan_dft_c2r_1d(ntt, ctr, data, FFTW_MEASURE);
    float fftscale = 1.0/((float) (ntt));  // fft scaling factor

    // calc linear merging taper
    for(i=0; i<ntaper; ++ i) mtaper[i] = (i + 1.0)/(ntaper + 1.0);  //increasing with index

    /* Read headers and data while getting a count */
    ntr = 0;
    do { /* Main loop over traces/gather */
        ++ntr;

        memcpy(data, &tr.data[itmin], nfilt*FSIZE);
        fftwf_execute_dft_r2c(planf, data, ctr);

        T = MAX(dt, tmin - tref + 0.5*twin);
        for (i=0; i<nstep; ++i, T += twin) {
            memset(cbuf, 0, (ntfft/2)*FSIZE);
            a = Qfilter(opt, mode, limit, Q, T, f0, nf, df, (complex*) ctr, cbuf);
            fftwf_execute_dft_c2r(planb, (fftwf_complex*) &cbuf[0], &trcs[i][0]);
            for (it=0; it<ntfft; ++it) trcs[i][it] *= fftscale;
        }

        // merge filtered window by window
        MergeTraces(ns, ntwin, ntaper, nstep, itmin, itmax, mtaper, trcs, tr.data, odata);
        memcpy(tr.data, odata, ns*FSIZE);
        puttr(&tr);

        memset(data, 0, ntfft*FSIZE);
        memset(*trcs, 0, nstep*ntfft*FSIZE);
    } while ((nsegy = gettr(&tr)) > HDRBYTES);

    if (verbose) warn(" Totally %d traces are processed", ntr);

    return (CWP_Exit());
}

float Qfilter(int opt, int mode, float limit, float Q, float T, float f0, int nf, float df, complex* ctr, complex* cfilt)
{
    int i;
    float f, rr, a, pshift;
    complex  p;

    float r = 1.0 / (PI * Q);   // gamma
    float fq = limit * Q / (PI * T);

    cfilt[0] = ctr[0];
    if (opt == 1) {
        for (i=1, f=df; i < nf; ++i, f += df) {
            rr = powf(f/f0, -r);
            a = expf(rr * PI * f * T / Q);
            pshift = -2.0 * PI * f * T * (rr - 1);
            p = cmplx(cosf(pshift), sinf(pshift));
            if (mode != 1) cfilt[i] = cmul(ctr[i], p);
            if (mode != 0) cfilt[i] = crmul(cfilt[i], 1.0/a);
        }

    } else {
        for (i=1, f=df; f <= fq && i < nf; ++i, f += df) {
            rr = powf(f/f0, -r);
            a = expf(rr * PI * f * T / Q);
            pshift = 2.0 * PI * f * T * (rr - 1);
            p = cmplx(cosf(pshift), sinf(pshift));
            if (mode != 1) cfilt[i] = cmul(ctr[i], p);
            if (mode != 0) cfilt[i] = crmul(cfilt[i], a);
        }

        for (; i < nf; ++i, f += df) {
            rr = powf(f/f0, -r);
            pshift = 2.0 * PI * f * T * (rr - 1);
            p = cmplx(cosf(pshift), sinf(pshift));
            if (mode != 1) cfilt[i] = cmul(ctr[i], p);
            if (mode != 0) cfilt[i] = crmul(cfilt[i], a);
        }
    }

    return a;  // max. value of amplitude scaling factor
}

void MergeTraces(int nt, int nlive, int nmix, int nw, int its, int ite, 
        float* mtaper, float** cdata, float* idata, float* odata)
{
    int i, iw;

    memset(odata, 0, nt*FSIZE);  // zero out to better debug
    if (its > 0) memcpy(&odata[0], &idata[0], its*FSIZE);  // copy data before filter window
    if (its > 0) { // blend first front edge for first subwindow
        for (i=0; i<nmix; ++i)
            odata[its + i] = (1.0 - mtaper[i])*idata[its + i] + mtaper[i]*cdata[0][i];
    } else { // simply copy from first filtered subwindow
        memcpy(&odata[0], &cdata[0][0], nmix*FSIZE);
    }
    memcpy(&odata[its + nmix], &cdata[0][nmix], (nlive - nmix)*FSIZE);  // should length (nlive - nmix)

    for (iw=1; iw<nw; ++iw) {
        for (i=0; i<nmix; ++i)  // blend front edge of subwindow with end of previous subwindow
            odata[its + iw*nlive + i]
                = (1.0 - mtaper[i])*cdata[iw - 1][iw*nlive + i] + mtaper[i]*cdata[iw][iw*nlive + i];
        // copy middle part of subwindow, overshoot nmix sample in case for last subwindow end exactly at ite
        memcpy(&odata[its + iw*nlive + nmix], &cdata[iw][iw*nlive + nmix], (nlive - nmix)*FSIZE); // should length (nlive - nmix)
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
