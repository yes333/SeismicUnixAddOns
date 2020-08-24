/* Copyright (c) Read Well Services, 2011.*/
/* All rights reserved.                       */

/* SUSVDI: $Revision: 1.0 $ ; $Date: 20011/04/07 17:11:15 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"

/*********************** self documentation *****************************/
char *sdoc[] = {
"                                                                               ",
" SPSVD - Singular Value Decomposition for extraction of aligned downgoing waves",
"           from VSP data                                                       ",
"                                                                               ",
" spsvd  <stdin >stdout [optional parameters]                                   ",
"                                                                               ",
" Optional parameters:                                                          ",
"   key=cdpt            gather sorting key                                      ",
"   nmax=128            max. number of traces in a gather                       ",
"                                                                               ",
"   output=0            =0 output eigensections or input data subtracted by them",
"                       =1 output eigenwavelet with singular values             ",
"                       =2 output scaling vectors                               ",
"                                                                               ",
"   mode=0              =0 all eigensections with singular values               ",
"                       =n up to n-th eigensections                             ",
"                                                                               ",
"   subtract=0          =1 subtract input by up to n-th eigensections           ",
"                                                                               ",
"   t1=0.0          [s] start of extraction time window                         ",
"   t2=(from hdr)   [s] end of extraction time window                           ",
"                                                                               ",
"   verbose=0           > 0 echo additional information                         ",
"                                                                               ",
" Notes:                                                                        ",
"   According to theory of Singular Value Decomposition, a matrix M[nt][ntr]    ",
"   can be factorized into nt wavelet vectors U, ntr singular values SV (global ",
"   scaling vector) and ntr (trace) scaling vectors V:                          ",
"                M[nt][ntr] = U[nt][nt]S[nt][ntr]VT[ntr][ntr]                   ",
"   where SV is diagonal matrix with S[i][i] = sv[i] sorted in descending order.",
"   In other words a seismic section [ntr x nt] can be decomposed into ntr eigen",
"   sections   e[i][itr][it] = sv[i]u[it][itr]v[i][itr],   i = 1 ... ntr;       ",
"   the first and second eigen sections represent the most dominant waves, which,",
"   on aligned VSP data, correspond the downgoing waves.                        ",
"   Input VSP data must be perfectly aligned with FTB",
"                                                                               ",
" Version 0.9.0 last modified April, 2011 by Sanyu Ye                           ",
"                                                                               ",
NULL};

/* 
 * Author: Sanyu Ye,
 *         Read Well Services, Oslo, 2011.
 *         E-mail: sanyu.ye@readgroup.com
 *
 *
 * References:
 *    Franco, R. de, and Musacchio, G., 2000: Polarization Filter with
 *       Singular Value Decomposition, submitted to Geophysics and
 *       published electronically in Geophysics online (www.geo-online.org).
 *    Jurkevics, A., 1988: Polarization analysis of three-comomponent
 *       array data, Bulletin of the Seismological Society of America, 
 *       vol. 78, no. 5.
 *    Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P.
 *       1996: Numerical Recipes in C - The Art of Scientific Computing,
 *       Cambridge University Press, Cambridge.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: none
 */
/**************** end self doc *******************************************/


/* function prototypes */
float computeEigenImg(const int j, const int ntr, const int nt, float** u, float* sval, float** v, float** eigenImg);

int  main(int argc, char **argv)
{
    cwp_String key;     /* header key word from segy.h */
    cwp_String type;    /* type of key	*/
    int index;     /* index of key	*/
    Value val, valnew;  /* value of key				*/
    int nsegy;              /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;    /* number of total input traces and gathers */
    int nt;                 /* number of points on trace		*/
    int i, itr, it; /* counters				*/
    int mode, subtract, output;
    int verbose;
    int nmax; /* max dimensions of input/output trace array */
    float dt;       /* time sample interval (sec)	*/
    float t1, t2;   // start and end of extraction window

    segy tr, trout;

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);
    
    /* Get parameters  */
    if (!getparint("verbose", &verbose))    verbose = 0;
    if (!getparint("output", &output))      output = 0;
    if (output < 0 || output > 2 ) err(" Invalid parameter output=%d. Valid value [0, 1, 2]", output);
    if (!getparint("mode", &mode))          mode = 0;
    if (!getparint("subtract", &subtract))  subtract = 0;
    if(!getparint("nmax", &nmax)) nmax = 128;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    nt = tr.ns;
    dt = ((double) tr.dt) / 1000000.0;

    if (!getparfloat("t1", &t1))    t1 = 0.0;
    if (!getparfloat("t2", &t2))    t2 = (nt-1)*dt;
    int its = NINT(t1/dt);   // start sample
    int ite = NINT(t2/dt);   // end sample
    int ntw = ite - its + 1; // window length in samples
        
    /* allocate space and reset to zero*/
//    float** eigenImg = ealloc2float(nt, nmax);    /* input data */
    float** trdata = ealloc2float(nt, nmax);    /* input data */
    segyhdr* hdrs = (segyhdr*) ealloc1(nmax, HDRBYTES);
    memset(*trdata, 0, nmax*nt*FSIZE);
    memset(hdrs, 0, nmax*HDRBYTES);
    
    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = ntotal = 0;
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nmax - 1) err("\nNumber of traces exceeding nmax=%d\n", nmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);  // copy header
            memcpy(trdata[ntr], tr.data, FSIZE*nt);
       
            ++ntr;
        } else { // new gather or END_OF_FILE
            ++ngather;
            ntotal += ntr;
            if (verbose > 1) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, val));

            float** udata = ealloc2float(ntr, ntw);  /* SVD matrix U */
            float** vdata = ealloc2float(ntr, ntr);    /* SVD matrix V */
            float* sval = ealloc1float(ntr);       /* singular values */
            float** eigenImg = ealloc2float(ntw, ntr);       /* eigen sections */
            float**   sumImg = ealloc2float(ntw, ntr);       /* eigen sections */
            memset(sval, 0, ntr*FSIZE);
            memset(*udata, 0, ntr*ntw*FSIZE);
            memset(*vdata, 0, ntr*ntr*FSIZE);
            memset(*eigenImg, 0, ntr*ntw*FSIZE);
            memset(*sumImg, 0, ntr*ntw*FSIZE);

            /* copy data into matrix for SVD */
            for (itr=0; itr<ntr; itr++) {
                for (it=its; it<ite; it++)
                    udata[it - its][itr] = trdata[itr][it];
            }

            /* perform singular value decomposition (SVD) */
            compute_svd(udata, ntw, ntr, sval, vdata);

            /* sort singular values in descending order */
            svd_sort(udata, sval, vdata, ntr, ntw);


            // output
            if (output == 1) {
                for(itr=0; itr<ntr; ++itr) {
                    memcpy(&trout, &hdrs[itr], HDRBYTES);  // copy trace head
                    memset(trout.data, 0, nt*FSIZE);
                    trout.ungpow = sval[itr];
                    for (it=0; it<ntw; ++it) trout.data[its + it] = udata[it][itr];
                    puttr(&trout);
                }
            }
            else if (output == 2) {
                for(itr=0; itr<ntr; ++itr) {
                    memcpy(&trout, &hdrs[itr], HDRBYTES);  // copy trace head
                    memset(trout.data, 0, nt*FSIZE);
                    trout.ungpow = sval[itr];
                    for (it=0; it<ntr; ++it) trout.data[it] = vdata[it][itr];;
                    puttr(&trout);
                }
            } else {
                int nmode = (mode == 0)? ntr : MIN(mode, ntr);
                for (i = 0; i < nmode; ++i) {
                    computeEigenImg(i, ntr, ntw, udata, sval, vdata, eigenImg);
                    if (mode == 0) {
                        for(itr=0; itr<ntr; ++itr) {
                            memcpy(&trout, &hdrs[itr], HDRBYTES);  // copy trace head
                            memcpy(&trout.data[its], eigenImg[itr], ntw*FSIZE);
                            trout.nhs = itr + 1;
                            trout.nvs = i + 1;
                            trout.ungpow = sval[i];
                            puttr(&trout);
                        }
                    } else {
                        for(itr=0; itr<ntr; ++itr) {
                            for (it=0; it<ntw; ++it) sumImg[itr][it] += eigenImg[itr][it];
                        }
                    }
                }

                if (mode != 0) {
                    for(itr=0; itr<ntr; ++itr) {
                        memcpy(&trout, &hdrs[itr], HDRBYTES);  // copy trace head
                        if (subtract == 1) {
                            trout.data[its + it] -= sumImg[itr][it];
                        } else {
                            trout.data[its + it]  = sumImg[itr][it];
                        }
                        puttr(&trout);
                    }
                }
            }

            // free memory
            free2float(eigenImg);
            free2float(sumImg);
            free2float(udata);
            free2float(vdata);
            free1float(sval);
    
            /* zero out data memory */
            memset(*trdata, 0, nmax*nt*FSIZE);
            memset(hdrs, 0, nmax*HDRBYTES);

            val = valnew;
            ntr = 0;
            continue;
        }
        nsegy = gettr(&trout);
    } while (!eof);

    if (verbose) warn(" Totally %d traces of %d gathers are processed", ntotal, ngather);

    return (CWP_Exit());
}

float computeEigenImg(const int j, const int ntr, const int nt, float** u, float* sval, float** v, float** eigenImg)
{
    int i, it;

    /* construct eigenimages */
    for (i=0; i<ntr; i++) {
        for (it=0; it<nt; it++) {
            eigenImg[i][it] = u[it][j] * v[i][j] * sval[j];
        }
    }
    return sval[j];
}