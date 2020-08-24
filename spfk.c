/* Copyright (c) READ well Service, 2010.*/
/* All rights reserved.                       */

/* SPFK: $Revision: 0.1 $ ; $Date: 20010/03/07 $		*/


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>
#include <fftw3.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SPFK - forward and inverse T-X,Y <--> F-Kx,Ky Fourier transform	",
"        using FFTW routines                                            ",
"                                                                       ",
"    spfk <infile >outfile  [optional parameters]                 	",
"                                                                       ",
" Optional Parameters:                                                  ",
"   key=cdpt                gather sort key				",
"                                                                       ",
"   dim=3                   (=2 or =1) Dimension of dataset        	",
"                                                                       ",
"   fft=-1                  =-1 forward Fourier transform          	",
"                           =1  inverse Fourier transform       	",
"                                                                       ",
"   plan=0                  flag for finding out best FFTW algorithm    ",
"                           =0  ESTIMATE, more for single test use      ",
"                           =1  MEASURE, for repeated transform         ",
"                           =2  PATIENT, for mass production            ",
"                                                                       ",
"   phaseshift=0            =1 apply e^(-i*k*x) on k axis to correct    ",
"                              spatial offset (xmin or ymin != 0)       ",
"                                                                       ",
" Dimensions of place holder for input/output gather in forward transform:      ",
"   ntmax=(tr.ns)           max. number of t/f samples                  ",
"   nxmax=(tr.sfs)          max. number of inline traces                ",
"   nymax=(tr.sfe)          max. number of x-line traces                ",
"   overs=1.0,1.0,1.0       Oversampling factors in T, X, Y dimensions  ",
"                                                                       ",
" NOTE: if gathers have varying dimensions in X or Y, the corresponding ",
"  values should be set >= max to make output gathers having same dimension.",
"  Otherwise dimension of first gather will be used.                    ",
"  Here the array nymax*nxmax*ntmax is padded with zero traces          ",
"                                                                       ",
" For forward transform                                                 ",
"   dx=tr.dx      [m]       x spacing                                   ",
"   dy=tr.dy      [m]       y spacing                                   ",
"   xmin=tr.fx    [m]       starting x of the gather                    ",
"   ymin=tr.fy    [m]       starting y                                  ",
"                                                                       ",
" Note:                                                                 ",
"   resulting f,df,Kx,dKx,Ky,dKy in unit [Hz] and [1/km] are stored in  ",
"   fz,dz,fx,dx,fy,dy while t0, dt, x, dx, y, dy in unit [sec] and [m]  ",
"                                                                       ",
"   verbose=0                  >0 echo information			",
"                                                                       ",
" Phase shift must be applied for further transform into tau-p. Header  ",
" words tr.fx/fy should be properly set using SUCHW.                    ",
" Oversampling factors should be used only for forward transform. For   ",
" further transform into Tau-P, an oversampling factor 4 is needed in Kx/Ky ",
" to suppress dipping noises (artifact) regardless if there is aliasing ",
" in FK domain.                                                         ",
" For inverse transform, oversampling is better facilitated by program  ",
" SPTAPER using wrap around and tapering. For example:                  ",
" sptaper wrap=50 taper=30 type=6                                       ",
"                                                                       ",
" Output of forward transform is scaled down by 1/N but not backward one        ",
" 1D transform is just like usual T-F transform, trace by trace                 ",
" 2D transform is like usual TX-FK                                              ",
" For 2/3D transform, input data must be sorted ascendingly in key, Y, X        ",
" Number of samples (tr.ns) of input X-T data for forward transfor should       ",
" usually be set as odd number to retain the same when inverse transformed.     ",
" Also the trace length should not be manipulated after forward transform.      ",
"                                                                               ",
" Run susetgather first to set gather dimensions correctly                      ",
"                                                                               ",
" Version 1.1   Last modified June 2011 by Sanyu Ye                             ",
NULL};

/*
 * Credits: RWS, Sanyu Ye (sanyu.ye@readgroup.com), March. 2010.
 *
 * Reference:
 *    WWW.FFTW.ORG
 *
 */
/**************** end self doc ********************************/


/** main **/

int main(int argc, char **argv) {
    cwp_String key;     /* header key word from segy.h		*/
    cwp_String type;    /* type of key				*/
    int index;          /* index of key				*/
    Value val, valnew;  /* value of key				*/
    int nsegy;		/* number of bytes read in the segy trace */
    int ntr;            /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;   /* number of total input traces and gathers */
    int nt, nf;         /* number of time or frequency samples */
    int nx, ny;		/* number of horizontal samples */
    int fft, dim;       /* direction of fft and dimension of input dataset (gather) */
    float dt, df;       /* Time/Frequecy sample interval */
    float dx, dy;       /* horizontal sample interval */
    float xmin, xmax, ymin, ymax;   /* min. max. extents of original gather before transform, recorded on first trace  */
    float dkx, dky;             /* increments for kx, ky transform */
    float xshift, yshift, phase; /* minimum x, y extent and phase */
    float fmin, fmax;		/* minimum frequency of interest */
    float kxmin, kymin;         /* actual extent of input gather [m] */
    float scalco, fftscale;		/* coordinate scale */
    float overs[3];		/* oversampling factors */

    int     ntmax, nxmax, nymax;    /* max imensions of input/output trace array */
    int     ixzero, iyzero;         /* position index of zero k, when centered */
    int     eof, IsNewGather, flag, phaseshift, center;
    int     ix, iy, it, i, j;
    float   T, f0, fNyquist, kxNyquist, kyNyquist;
    complex** ps = NULL;  // phaseshift matrix

    fftwf_plan plan1 = NULL, plan2 = NULL, plan3 = NULL;

    int verbose;	 /* flag for echoing information */

    segy tr, trout;

    /* hook up getpar to handle the parameters */
    initargs(argc, argv);
    requestdoc(1);

    if(!getparint("verbose", &verbose)) verbose = 0;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    if(!getparint("dim", &dim)) dim = 3;
    if(!getparint("fft", &fft)) fft = FFTW_FORWARD;
    if(!getparint("phaseshift", &phaseshift)) phaseshift = 0;
    if(fft == FFTW_BACKWARD && tr.trid != FUNPACKNYQ)
        err("input not complex freq data for inverse FFT, trid=%d (!=%d)", tr.trid, FUNPACKNYQ);
    if(!getparint("center", &center)) center = 1;
    if(!getparint("plan", &flag)) flag = 0;
    if(flag == 2)       flag = FFTW_PATIENT;
    else if(flag == 1)  flag = FFTW_MEASURE;
    else                flag = FFTW_ESTIMATE;

    // apply oversampling if specified
    int nc = countparval("overs");
    if (nc == 0) {
        for (i=0; i<dim; ++i) {
            overs[i] = 1.0;
        }
    } else if (nc > dim) {
        err("  Too many oversampling factors (%d > %d)", nc, dim);
    } else {
        getparfloat("overs", overs);
        for (i=0; i<nc; ++i) {
            if (overs[i] < 1.0) err("  Oversampling factor must be larger than 1");
        }
        for (i=nc; i<dim; ++i) overs[i] = overs[nc - 1];  // replicate last one specified
        //if (dim > 1) nxmax = NINT(overs[1]*nxmax);
        //if (dim > 2) nymax = NINT(overs[2]*nymax);
        if (verbose) {
            fprintf(stderr, "  %d oversampling factors: ", nc);
            for (i=0; i<dim; ++i) fprintf(stderr, "%4.2f ", overs[i]);
            fprintf(stderr, "\n");
        }
    }

    dt = df = dx = dy = dkx = dky = 0.0;
    xmin = xmax = ymin = ymax = T = 0.0;

    nt = tr.ns;
    if (fft == FFTW_FORWARD) { /* forward transform */
        nt = NINT(overs[0]*nt);
        int ntfft = 2*(nt/2 + 1);  // default/intrinsic 1st fft dimension
        // figure out cache dimension of trace (ntmax)
        if(!getparint("ntmax", &ntmax)) {
            ntmax = ntfft ;
        } else {
            if (ntmax % 2) // odd number
                err("ntmax (=%d) must be even", ntmax);
            if (ntmax < ntfft) {
                warn("too small ntmax, reset (%d --> %d) to accommondate trace length (tr.ns=%d)", ntmax, ntfft, tr.ns);
                ntmax = ntfft;
            } else {    // aquivalent to padding trace
                nt = ntmax - 1;
                warn("reset nt (%d --> %d) to adapt input ntmax (=%d)", tr.ns, nt, ntmax);
            }
        }

        dt = ((float) tr.dt)*1.e-6;
    } else {
        ntmax = nt;  // for backward fftw, nt can't be manipulated after forward fftw

        kxmin = tr.fx;
        kymin = tr.fy;
        dkx = tr.dx;
        dky = tr.dy;
        f0 = tr.fz;
        df = tr.dz;
        nf = tr.ns/2;

        if(!getparfloat("xmin", &xmin)) xmin = tr.fx;
        if(!getparfloat("ymin", &ymin)) ymin = tr.fy;
    }

    // always reduce logical dimension to odd number, because for backward fftw no way to know original number
    nt = ntmax - 1;  

    if(!getparfloat("xmin", &xmin)) xmin = tr.fx;
    if(!getparfloat("ymin", &ymin)) ymin = tr.fy;
    if(!getparfloat("dx", &dx)) dx = tr.dx;
    if(!getparfloat("dy", &dy)) dy = tr.dy;
    // get actual x & y-dimension
    nx = (dim < 2)? 1 : tr.sfs;
    ny = (dim < 3)? 1 : tr.sfe;
    if(!getparint("nxmax", &nxmax)) nxmax = 0;
    if(!getparint("nymax", &nymax)) nymax = (dim > 2)? 0 : 1;
    if ( dim > 2 ) { // 3D
        if ( nymax < NINT(overs[2]*ny) )  nymax = NINT(overs[2]*ny);
        if (ny < nymax/4) {
            warn("Number of trace along Y may not be properly set (ny=%d --> %d)", ny, nymax);
            ny = nymax;
        }
    }
    if (dim > 1) {  // 2/3 D
        if ( nxmax < NINT(overs[1]*nx) )  nxmax = NINT(overs[1]*nx);
        if (nx < nxmax/4) {
            warn("Number of trace along X may not be properly set (nx=%d --> %d)", nx, nxmax);
            //nx = nxmax;
        }
    }
    if ( dim == 2 ) {
        nymax = ny = 1;
        ymin = dy = dky = 0.0;
    }
    if ( dim == 1 ) {
        nxmax = nx = 1;
        nymax = ny = 1;
        xmin = dx = dkx = 0.0;
        ymin = dy = dky = 0.0;
        phaseshift = 0;  // no ps for 1-d
    }

    // determine sampling intervals
    if (fft == FFTW_FORWARD) {
        fNyquist = 0.5/dt;
        df = fNyquist/(ntmax/2 - 1);
        if (dim > 1) {
            if ( dx == 0.0 ) err ("dx must be set");
            if ( nxmax <= 1 ) err ("nxmax (=%d) too small");
            kxNyquist = 0.5*1000.0/dx;
            dkx = 2.0*kxNyquist/(nxmax - 1);
        }
        if (dim > 2) {
            if ( dy == 0.0 ) err ("dy must be set");
            if ( nymax <= 1 ) err ("nymax (=%d) too small");
            kyNyquist = 0.5*1000.0/dy;
            dky = 2.0*kyNyquist/(nymax - 1);
        }
    } else {
        T = 1.0/df, dt = T/(nt - 1);
        if (dim > 1) {
            if ( dkx == 0.0 ) err ("dkx (tr.dx) must be set");
            if ( nxmax <= 1 ) err ("nxmax (=%d) too small");
            xmax = 0.5*1000.0/dkx;
            dx = 2.0*xmax/(nxmax - 1);
        }
        if (dim > 2) {
            if ( dky == 0.0 ) err ("dky (tr.dy) must be set");
            if ( nymax <= 1 ) err ("nymax (=%d) too small");
            ymax = 0.5*1000.0/dky;
            dy = 2.0*ymax/(nymax - 1);
        }
    }

    if (verbose && dim > 1) {
        if (fft == FFTW_FORWARD) {
            warn("%d-D forward transform: XMIN=%.1f DX=%.2f  NX=%d  NXMAX=%d ", dim, xmin, dx, nx, nxmax);
        } else {
            warn("%d-D backward transform: KXMIN=%.1f DKX=%.4f  NX=%d  NXMAX=%d ", dim, kxmin, dkx, nx, nxmax);
        }
        if (phaseshift) warn("  Phaseshift XMIN=%.2f ", xmin);
        if (dim > 2) {
            if (fft == FFTW_FORWARD) {
                warn("3-D forward transform: YMIN=%.1f DY=%.2f  NY=%d  NYMAX=%d ", ymin, dy, ny, nymax);
            } else {
                warn("3-D backward transform: KYMIN=%.1f DKY=%.4f  NY=%d  NYMAX=%d ", kymin, dky, ny, nymax);
            }
            if (phaseshift) warn("  Phaseshift YMIN=%.2f ", ymin);
        }

    }

    if (verbose) warn("array dimensions used (NT x NX x NY): %d x %d x %d", ntmax, nxmax, nymax);

    if(!getparfloat("fmin", &fmin)) fmin = 3.0;
    if(!getparfloat("fmax", &fmax)) fmax = 60.0;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* allocate memory input/output data traces*/
    float*** trdata = ealloc3float(ntmax, nxmax, nymax);
    fftwf_complex* ctr = (fftwf_complex*) &trdata[0][0][0];
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nxmax, nymax, HDRBYTES);
    /* zero out data memory */
    memset(**trdata, 0, nymax*nxmax*ntmax*FSIZE);
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

    if(phaseshift) ps = ealloc2complex(nxmax, nymax);

    int nhalf = ntmax/2;
    int nreal = ntmax*nxmax*nymax;
    int ncmpl = nreal/2;
    if ( fft == FFTW_FORWARD ) {
        plan1 = fftwf_plan_many_dft_r2c(1, &nt, nxmax*nymax, &trdata[0][0][0], NULL, 1, ntmax, ctr, NULL, 1, nhalf, flag);
        if ( dim > 1) plan2 = fftwf_plan_many_dft(1, &nxmax, nhalf, ctr, NULL, nhalf, 1, ctr, NULL, nhalf, 1, FFTW_BACKWARD, flag);
        if ( dim > 2) plan3 = fftwf_plan_many_dft(1, &nymax, nhalf*nxmax, ctr, NULL, nhalf*nxmax, 1, ctr, NULL, nhalf*nxmax, 1, FFTW_BACKWARD, flag);
    } else {
        plan1 = fftwf_plan_many_dft_c2r(1, &nt, nxmax*nymax, ctr, &ncmpl, 1, nhalf, &trdata[0][0][0], &nreal, 1, ntmax, flag);
        if ( dim > 1) plan2 = fftwf_plan_many_dft(1, &nxmax, nhalf, ctr, NULL, nhalf, 1, ctr, NULL, nhalf, 1, FFTW_FORWARD, flag);
        if ( dim > 2) plan3 = fftwf_plan_many_dft(1, &nymax, nhalf*nxmax, ctr, NULL, nhalf*nxmax, 1, ctr, NULL, nhalf*nxmax, 1, FFTW_FORWARD, flag);
    }
    //fftscale = 1.0/sqrt((float) (nymax*nxmax*nt));  // fft scaling factor
    fftscale = 1.0/((float) (nymax*nxmax*nt));  // fft scaling factor
    memset(**trdata, 0, nymax*nxmax*ntmax*FSIZE);  // reset to zero 

    /* loop over traces */
    ntr = 0, ntotal = 0, ngather = 0;
    iy = 0, ix = 0; eof = 0;
    do {
        if (nsegy > HDRBYTES) {
            gethval(&tr, index, &valnew);
            IsNewGather = valcmp(type, val, valnew);
        } else {
            eof = 1; //END_OF_FILE
        }

        if(ntr == 0 && dim > 1) { // first trace, check gather specific parameters
            xshift = xmin/1000.0;
            yshift = ymin/1000.0;
            if ( fft == FFTW_BACKWARD ) {
                xshift = -xshift;
                yshift = -yshift;
            }
            // calculate phaseshift matrix
            if (phaseshift) {
                for (j = 0; j < nymax; ++j) {
                    for (i = 0; i < nxmax; ++i) {
                        phase = 2.0*PI*xshift*i*dkx;
                        ps[j][i] = cmplx(cos(phase), sin(phase));
                    }
                    if (dim == 3) {
                        phase = 2.0*PI*yshift*j*dky;
                        ps[j][i] = cmul(ps[j][i], cmplx(cos(phase), sin(phase)));
                    }
                }
            }

            if ( tr.sfs != 0 && nx != tr.sfs ) {
                warn("  Diviate array dimension nx (%d != %d) at %d-th gather (%s=%d)", nx, tr.sfs, ngather+1, key, vtoi(type, valnew));
                nx = tr.sfs;
            }
            if ( dim == 3 && tr.sfe != 0 && ny != tr.sfe ) {
                warn("  Diviate array dimension ny (%d != %d) at %d-th gather (%s=%d)", ny, tr.sfe, ngather+1, key, vtoi(type, valnew));
                ny = tr.sfe;
            }

            ixzero = nx/2;
            iyzero = ny/2;
        }

        // cache in input data
        if ( (nsegy > HDRBYTES && !IsNewGather) /* same key and more data*/
                || (dim == 1 && !eof) ) { // 1-D but not eof always
            j = iy = ntr/nx;
            i = ix = ntr%nx;

            // apply phaseshift before back transform
            if (phaseshift && fft == FFTW_BACKWARD) {
                complex* cdata=(complex*) tr.data;
                    for(it=0; it<tr.ns/2; ++it, ++cdata) *cdata = cmul(*cdata, ps[iy][ix]);
            }

            if (ntr > nxmax*nymax - 1) err("  array dimension too small (%d < %d) traces input", nxmax*nymax, ntr+1);
            if (ix > nxmax - 1) err("  array dimension nxmax too small (%d < %d) traces input", nxmax, ix+1);
            if (iy > nymax - 1) err("  array dimension nymax too small (%d < %d) traces input", nymax, iy+1);
            // reposition trace according to transform type because of in-order placement of FFTW
            if (center && dim > 1 &&  fft == FFTW_BACKWARD ) {
                i -= ixzero; // shift back
                if (i < 0) i += nx; // if negative, shift to right
                if ( dim > 2) {
                    j -= iyzero; // shift back
                    if (j < 0) i += ny; // if negative, shift to right
                }

                if (verbose == 333) fprintf(stderr, "i=%4d ix=%4d kx=%7.3f   j=%4d iy=%4d kx=%7.3f\n", i, ix, tr.fx, j, iy, tr.fy);
            }
            memcpy(&hdrs2d[j][i], &tr, HDRBYTES);
            memcpy(&trdata[j][i][0], tr.data, FSIZE*tr.ns);


            ++ntr;
        }

        // do transform
        if ( (dim == 1 && !eof) // 1-D but not eof always do transform
                || (dim > 1 && (IsNewGather || eof) ) ) { // new gather or END_OF_FILE
            ++ngather;
            if (verbose > 1) warn("  Transforming %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, valnew));

            ntotal += ntr;

            // do fft
            if ( fft == FFTW_FORWARD ) {
                fftwf_execute(plan1);
                if ( dim > 1) { // loop over Y for T-X plane
                    for (j = 0; j < nymax; ++j ) {
                        ctr = (fftwf_complex*) trdata[j][0];
                        fftwf_execute_dft(plan2, ctr, ctr);
                    }
                }
                ctr = (fftwf_complex*) trdata[0][0];
                if ( dim > 2) fftwf_execute(plan3);
            } else {
                if ( dim > 2) fftwf_execute(plan3);
                if ( dim > 1) { // loop over Y for T-X plane
                    for (j = 0; j < nymax; ++j ) {
                        ctr = (fftwf_complex*) trdata[j][0];
                        fftwf_execute_dft(plan2, ctr, ctr);
                    }
                }
                ctr = (fftwf_complex*) trdata[0][0];
                fftwf_execute(plan1);
            }


            // output trace
            ixzero = (nxmax + 1)/2;
            iyzero = (nymax + 1)/2;
            for (j = 0; j < nymax; ++j) {
                for (i = 0; i < nxmax; ++i) {
                    ix = i, iy = j;

                    // reposition output f-k traces
                    if (center && dim > 1 &&  fft == FFTW_FORWARD) {
                        ix += ixzero;   // zero kx index
                        if (ix > nxmax - 1) ix -= nxmax;
                        if (dim > 2) {
                            iy += iyzero;   // zero kx index
                            if (iy > nymax - 1) iy -= nymax;
                        }
                    }

                    // use header of last real trace for traces beyond
                    int ixh = (ix < nx)? ix : nx-1;
                    int iyh = (iy < ny)? iy : ny-1;

                    // set grid parameters of trace header
                    memcpy(&trout, &hdrs2d[iyh][ixh], HDRBYTES);
                    if ( fft == FFTW_FORWARD ) {
                        if ( dkx != 0.0 ) trout.dx = trout.d2 = dkx;
                        if ( dky != 0.0 ) trout.dy = dky;
                        if ( dkx != 0.0 ) trout.fx = trout.f2 = ((center)? -kxNyquist : 0.0) + i*dkx;
                        if ( dky != 0.0 ) trout.fy = ((center)? -kyNyquist : 0.0) + j*dky;
                    } else {
                        if ( dx != 0.0 ) trout.dx = trout.d2 = dx;
                        if ( dy != 0.0 ) trout.dy = dy;
                        if ( dx != 0.0 ) trout.fx = trout.f2 = xmin + i*dx;
                        if ( dy != 0.0 ) trout.fy = ymin + j*dy;
                    }

                    if (verbose == 333) fprintf(stderr, "i=%4d ix=%4d kx=%7.3f  j=%4d iy=%4d ky=%7.3f\n", i, ix, tr.fx, j, iy, tr.fy);

                    // set common parameters
                    trout.dt = NINT(1000000.0*dt);
                    trout.counit = -fft*dim;  // set tranform type for next step use
                    if (dim >= 2) {
                        trout.shortpad = nxmax;
                        trout.sfs = nxmax;
                        trout.nhs = i + 1;
                    }
                    if (dim == 3) {
                        trout.ntr = nymax*nxmax;
                        trout.sfe = nymax;
                        trout.nvs = j + 1;
                    }
                    if ( fft == FFTW_FORWARD ) {
                        trout.trid = FUNPACKNYQ;
                        trout.fz = trout.f1 = 0.0;
                        trout.dz = trout.d1 = df;
                        trout.ns = ntmax;
                        if ( dx != 0.0 ) trout.afilf  = dx*100000.0;  // in centimeter, cache for later use
                        if ( dx != 0.0 ) trout.afils  = xshift*1000.0;  // in meter
                        if ( dy != 0.0 ) trout.nofilf = dy*100000.0;  // in centimeter, cache for later use
                        if ( dy != 0.0 ) trout.nofils = yshift*1000.0;  // in meter
                    } else {
                        trout.trid = 1;
                        trout.fz = trout.f1 = 0.0;
                        trout.dz = trout.d1 = dt;
                        trout.ns = nt;
                        if ( dkx != 0.0 ) trout.lcf = dkx*100000.0;  // in centimeter, cache for later use
                        if ( dkx != 0.0 ) trout.hcf = kxNyquist*100.0;  // in 10 m
                        if ( dky != 0.0 ) trout.lcs = dky*100000.0;  // in centimeter, cache for later use
                        if ( dky != 0.0 ) trout.hcs = kyNyquist*100.0;  // in 10 m
                    }
                    memcpy(trout.data, &trdata[iy][ix][0], FSIZE * trout.ns);
                    // apply phaseshift matrix for forward fft
                    if (phaseshift && fft == FFTW_FORWARD) {
                        complex* cdata=(complex*) trout.data;
                        for(it=0; it<trout.ns/2; ++it, ++cdata)
                            *cdata = cmul(*cdata, ps[j][i]);
                    }
                    if ( fft == FFTW_FORWARD ) for(it=0; it<trout.ns; ++it) trout.data[it] *= fftscale;
                    puttr(&trout);
                }
            }

            val = valnew;
            // reset data cache
            memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);
            memset(&trdata[0][0][0], 0, nymax*nxmax*ntmax*FSIZE);
            ntr = 0, iy = 0, ix = 0;
            if (dim > 1) continue; // skip reading next trace upon new gather
        }
        nsegy = gettr(&tr);
    } while (!eof);

    fftwf_destroy_plan(plan1);
    if (dim > 1) fftwf_destroy_plan(plan2);
    if (dim > 2) fftwf_destroy_plan(plan3);
    free2((void**)hdrs2d);
    free3float(trdata);
    if(ps) free2complex(ps);
    
    return(CWP_Exit());
}
