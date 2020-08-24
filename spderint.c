/* Copyright (c) READ well Service, 2010.   */
/* Written by Sanyu Ye
/* All rights reserved.                       */

/* SPDERINT: $Revision: 1.0 $ ; $Date: 20010/08/17 $		*/


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPDERINT - perform time derivative and integral on F-X traces                 ",
"                                                                               ",
"    spderint <infile >outfile  [optional parameters]                           ",
"                                                                               ",
" Optional Parameters:                                                          ",
"                                                                               ",
"   mode=0                      =0  time derivative of complex frequency trace  ",
"                               =1  time integral                               ",
"                                                                               ",
"   fmin=1.             [Hz]    low cut-off frequency                           ",
"   flow=2.             [Hz]    low pass frequency                              ",
"   fhign=0.75*fNyquest [Hz]    high pass frequency                             ",
"   fmax=fNyquest       [Hz]    high cut-off frequency                          ",
"                                                                               ",
"   type=3                      taper type					",
"                               =1 linear (default)                             ",
"                               =2 cosine					",
"                               =3 sine                                         ",
"                               =4 sine^2					",
"                               =5 gaussian (+/-2.0)                            ",
"                               =6 gaussian (+/-3.8)                            ",
"                               =n (n>6) even stronger gaussian (+/-n)          ",
"                                                                               ",
"   verbose=0                  >0 echo information                              ",
"                                                                               ",
" Note:                                                                         ",
"   Derivative or integral is computed by multiplying or dividing circular      ",
"   frequency. To avoid abnormal boost of either lower or higher end of spectrum",
"   min. and max. frequencies are used to taper and cut off amplitude spectrum. ",
"   spfk dim=1 should be used before and after to do t-f transform.             ",
"                                                                               ",
" Version 1.0.1   Last updated Nov, 2012 by Sanyu Ye                            ",
"                                                                               ",
NULL};

/*
 * Credits: RWS, Sanyu Ye (sanyu.ye@readgroup.com), August. 2010.
 */
/**************** end self doc ********************************/

/* forward declaration of prototype functions */
float* MakeTaper(int n, int type, int verbose);

/** main **/

int main(int argc, char **argv) {
    int nsegy;              /* number of bytes read in the segy trace */
    int ntr, nf;            /* number of actual input traces and frequencies */
    int mode, type;
    int i, imin, imax, ilow, ihigh;
    float fmin, fmax, flow, fhigh, df, fNyquist;
    float *taper1, *taper2;
    complex *data, czero=cmplx(0.0, 0.0);

    int verbose;	 /* flag for echoing information */

    segy tr;

    /* hook up getpar to handle the parameters */
    initargs(argc, argv);
    requestdoc(0);

    if(!getparint("verbose", &verbose)) verbose = 0;

    if(!getparint("mode", &mode)) mode = 0;
    if(!getparint("type", &type)) type = 3;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");
    if(!(tr.trid == FUNPACKNYQ)) // input must be complex trace in f-x domain
        err("input not complex freq data, trid=%d (!=%d)", tr.trid, FUNPACKNYQ);

    nf = tr.ns / 2;
    df = tr.dz;
    fNyquist = nf*df;

    if(!getparfloat("fmin", &fmin))     fmin = 1.0;
    if(!getparfloat("flow", &flow))     flow = 2.0;
    if(!getparfloat("fhigh", &fhigh))   fhigh = 0.75*fNyquist;
    if(!getparfloat("fmax", &fmax))     fmax = fNyquist;

    if (verbose) {
        warn("%s", mode? "Trace time derivative" : "Trace time integral");
        warn("Frequebcies input (nf=%d  fmin=%3.1f flow=%3.1f fhigh=%4.1f fmax=%4.1f fNyquist=%4.1f",
                nf, fmin, flow, fhigh, fmax, fNyquist);
    }
    imin  = NINT(fmin/df);
    imax  = NINT(fmax/df);
    ilow  = NINT(flow/df);
    ihigh = NINT(fhigh/df);

    if (ilow - imin < 2) {
        warn("CAUTION: Too few samples for tapering in lower f domain. Increase trace length");
    }

    taper1 = MakeTaper(ilow - imin,  type, verbose);
    if (!taper1) err("no sample for lower tapering window");
    taper2 = MakeTaper(imax - ihigh, type, verbose);
    if (!taper2) err("no sample for higher tapering window");

    /* loop over traces */
    ntr = 0;
    do {
        ++ntr;
        data = (complex*) tr.data;
        for (i = 0; i < imin; ++i) data[i] = cmul(data[i], czero);
        for (i = imin; i < ilow; ++i) {
            data[i] = mode? cdiv(data[i], cmplx(0, 2.0*PI*i*df)) : cmul(data[i], cmplx(0, 2.0*PI*i*df));
            data[i] = cmul(data[i], cmplx(taper1[i-imin], 0.0));
        }
        for (i = ilow; i < ihigh; ++i) {
            data[i] = mode? cdiv(data[i], cmplx(0, 2.0*PI*i*df)) : cmul(data[i], cmplx(0, 2.0*PI*i*df));
        }
        for (i = ihigh; i < imax; ++i) {
            data[i] = mode? cdiv(data[i], cmplx(0, 2.0*PI*i*df)) : cmul(data[i], cmplx(0, 2.0*PI*i*df));
            data[i] = cmul(data[i], cmplx(taper2[imax - i - 1], 0.0));
        }
        for (i = imax; i < nf; ++i) data[i] = cmul(data[i], czero);

        puttr(&tr);
    } while ((nsegy = gettr(&tr)) > HDRBYTES);
    
    return(CWP_Exit());
}

float* MakeTaper(int n, int type, int verbose)
{
    int i;
    float env = 0.0, f, x;
    const float min = 0.0, max = 1.0;
    const float EPS = 3.8090232;    /* exp(-EPS*EPS) = 5e-7, "noise" level  */

    if (n < 1) return NULL;

    float* w = ealloc1float(n);

    for (i = 0; i < n; i++) {
        f = (float) (i) / n;
        switch (type) {
            case 1: env = min + (max - min) * f;
                break;
            case 2: env = 0.5 * (1.0 - cos(PI * f));
                break;
            case 3: env = sin(PI * f / 2.);
                break;
            case 4: env = sin(PI * f / 2.)*sin(PI * f / 2.);
                break;
            case 5: x = 2.0 * (1.0 - f);
                env = exp(-(x * x));
                break;
            case 6: x = EPS * (1.0 - f);
                env = exp(-(x * x));
                break;
            default: x = (float) (type) * (1.0 - f);
                env = exp(-(x * x));
        }
        w[i] = env;
    }

    if (verbose > 1) {
        fprintf(stderr, "  %d weighting factors for taper: ", n);
        for (i=0; i<n; ++i) {
            if ( !(i%20) ) fprintf(stderr, "\n");
            fprintf(stderr, "%.3f ", w[i]);
        }
        fprintf(stderr, "\n");
    }

    return w;
}
