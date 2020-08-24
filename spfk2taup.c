/* 
 * File:   spfk2taup.c
 * Author: yes
 *
 * Created on 22 February 2010, 10:01
 */

#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>
#include <fftw3.h>

/*
 * 
 */

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SPFK2TAUP - forward and inverse F-Kx,Ky <--> F/Tau-Px,Py transform	",
"                                                    ",
"                                                                       ",
"    spfk2taup <infile >outfile  [optional parameters]                 	",
"                                                                       ",
" Optional Parameters:                                                  ",
"   key=cdpt                gather sort key				",
"                                                                       ",
"   dim=3                   (=2 or =1) Dimension of dataset        	",
"                                                                       ",
"   dir=1                   =1  forward transform          	",
"                           =-1 inverse transform       	",
"                                                                       ",
"   ftau=1                  flag for in/output traces in Tau-P or F-P domain  ",
"                           =1  tau-p,        ",
"                           =0  f-p          ",
"                                                                       ",
" Dimensions of place holder for input/output gather:                   ",
"   ntmax=(tr.ns)           max. number of t/f samples          ",
"   nxmax=(tr.sfs)          max. number of inline traces      ",
"   nymax=(tr.sfe)          max. number of x-line traces              ",
"                                                                       ",
" In general the input and output parameters ",
" f/t,df/dt,Kx/Px,dKx/dPx,Ky/Py,dKy/dPy are stored in fz,dz,fx,dx,fy,dy  ",
"                                                                       ",
" For forward transform                                                                  ",
"   npx=nxmax                number of slopes for Px transform	",
"   npy=nymax                number of slopes for Py transform	",
"   pxmax=1.0    [s/km]      maximum slope for Tau-P transform 	",
"   vxmax=1000   [m/s]       or alternatively using velocity	",
"   pxmin=-pxmax [s/km]      minimum slope for Tau-P transform (s/km)	",
"   pymax=pxmax  [s/km]      maximum slope for Tau-P transform (s/km)	",
"   vymax=1000   [m/s]       or alternatively using velocity	",
"   pymin=-pymax             minimum slope for Tau-P transform (s/km)	",
"                                                                   ",
" For inverse transform                                                                  ",
"   nkx=nxmax                number of samples of Kx	",
"   nky=nymax                number of samples of Ky 	",
"   kxmax=40     [1/km]      maximum wave number of Kx	",
"   dx=12.5      [m]         or using corresponding x sampling interval	",
"   kymax=40     [1/km]      maximum wave number of Ky	",
"   dy=12.5      [m]         or using corresponding y sampling interval	",
"								",
"   verbose=0                >0 echo information			",
"                                                                       ",
" NOTE: input and output gathers in K domain must be always symmetric   ",
"  with zero K centered ",
"                             ",
"",
" Output transform is scaled by 1/sqrt(NT) to be symmetric if fft applied",
" Input data must be sorted ascendingly in key, Y, X",
"                                                       	",
"                                                                       ",
NULL};

/*
 * Credits: RWS, Sanyu Ye (sanyu.ye@readgroup.com), April. 2010.
 *
 * Reference:
 *    WWW.FFTW.ORG
 *
 */
/**************** end self doc ********************************/

void fk2fp(int dir, int nf, float f0, float df, int nkx, float kxmin, float dkx, int nky, float kymin, float dky,
	  int npx, float pxmin, float dpx, int npy, float pymin, float dpy, float ***trdata);
void intp_wk2wp (int nw, float fw, float dw, int nk, float fk, float dk,
	 int np, float fp, float dp, complex **tr_k);
void intp_wp2wk (int nw, float fw, float dw, int nk, float fk, float dk,
	 int np, float fp, float dp, complex **tr_k);

/** main **/

int main(int argc, char **argv) {
    cwp_String key;     /* header key word from segy.h		*/
    cwp_String type;    /* type of key				*/
    int index;          /* index of key				*/
    Value val, valnew;  /* value of key				*/
    int nsegy;		/* number of bytes read in the segy trace */
    int ntr;            /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;/* number of total input traces and gathers */
    int nt, nf;         /* number of time or frequency samples */
    int nx, ny, nkx, nky, npx, npy, nxout, nyout; /* number of horizontal samples */
    int dir, dim;       /* direction of transform and dimension of input dataset (gather) */
    float dt, df;       /* Time/Frequecy sample interval */
    float dx, dy;       /* horizontal sample interval */
    float pxmin, pxmax, pymin, pymax;   /* min. max. extents of Px-Py  */
    float vxmin, vxmax, vymin, vymax;   /* equivalent of Vx, Vy for min. max. extents of Px-Py  */
    float xmin, xmax, ymin, ymax;   /* min. max. extents of original gather before transform, recorded on first trace  */
    float dkx, dky, dpx, dpy;             /* increments for kx, ky transform */
    float fmin, fmax;		/* minimum frequency of interest */
    float kxmin, kxmax, kymin, kymax;	/* actual extent of input gather [m] */
    float scalco, fftscale;		/* coordinate scale */

    int     ntmax, nxmax, nymax; /* max imensions of input/output trace array */
    int     eof, IsNewGather, ftau;
    int     ix, iy, it, i, j;
    float   T, f0, fNyquist, kxNyquist, kyNyquist;
    cwp_Bool IsComplex;  // complex trace data for f domain

    int verbose;	 /* mode for echoing information */

    segy tr, trout;
    segyhdr **hdrs2d;      /* input  header array */

    /* hook up getpar to handle the parameters */
    initargs(argc, argv);
    requestdoc(1);

    if(!getparint("verbose", &verbose)) verbose = 0;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    if(!getparint("dim", &dim)) dim = 3;
    if (dim < 2 || dim > 3) err("  Only 2- or 3-D transform supported");
    if(!getparint("dir", &dir)) dir = 1;
    if ( !(dir == -1 || dir == 1) ) err("dir (=%d) can only be -1 or 1", dir);
    if(!getparint("ftau", &ftau)) ftau = 1;
    if ( !(ftau == 0 || ftau == 1) ) err("ftau (=%d) can only be 0 or 1", ftau);

    IsComplex = tr.trid == FUNPACKNYQ;
    if(dir > 0 && !IsComplex ) // for forward transform input always f-k
        err("input not complex freq data, trid=%d (!=%d)", tr.trid, FUNPACKNYQ);
    if(dir < 0 && ftau != 1 && !IsComplex)
        err("input not complex freq data, trid=%d (!=%d)", tr.trid, FUNPACKNYQ);

    dt = df = dpx = dpy = dkx = dky = 0.0;
    kxmin = kxmax = kymin = kymax = 0.0;
    pxmin = pxmax = pymin = pymax = 0.0;
    xmin = xmax = ymin = ymax = 0.0;

    nt = tr.ns;
    if (dir > 0) { /* forward transform */
        // figure out cache dimension of trace (ntmax)
        if(!getparint("ntmax", &ntmax)) ntmax = nt;
        if (ntmax < nt) {
            warn("too small ntmax, reset (%d --> %d) to accommondate trace length (tr.ns=%d)", ntmax, nt, tr.ns);
            ntmax = nt;
        }
        if (ntmax % 2) { // odd number
            warn("odd ntmax, reset (%d --> %d) to accommondate trace length (tr.ns=%d)", ntmax, 2*(ntmax/2 + 1), tr.ns);
            ntmax = 2*(ntmax/2 + 1);
        }
        if (ftau == 1) { // output tau-p
            nt = ntmax - 1;
            if (verbose) warn("reset output nt (%d --> %d) to adapt input ntmax (=%d)", tr.ns, nt, ntmax);
        }

        if(!getparfloat("pxmax", &pxmax)) {
            if(getparfloat("vxmax", &vxmax)) pxmax = 0.001/vxmax;
            else pxmax = 1.0;
        }
        if(!getparfloat("pxmin", &pxmin)) {
            if(getparfloat("vxmin", &vxmin)) pxmin = 0.001/vxmin;
            else pxmin = -pxmax;
        }
        if(!getparfloat("pymax", &pymax)) {
            if(getparfloat("vymax", &vymax)) pymax = 0.001/vymax;
            else pymax = 1.0;
        }
        if(!getparfloat("pymin", &pymin)) {
            if(getparfloat("vymin", &vymin)) pymin = 0.001/vxmin;
            else pymin = -pymax;
        }

        dkx = tr.dx;
        dky = (dim > 2)? tr.dy : 0.0;
        kxmin = tr.fx;
        kymin = (dim > 2)? tr.fy : 0.0;
        f0 = tr.fz;
        df = tr.dz;
        nf = tr.ns/2;
    } else {
        // figure out cache dimension of trace (ntmax)
        if(!getparint("ntmax", &ntmax)) ntmax = (ftau == 1)? 2*(nt/2 +1) : nt;
        else {
            if (ntmax < nt) {
                warn("too small ntmax, reset (%d --> %d) to accommondate trace length (tr.ns=%d)", ntmax, 2*(nt/2 + 1), tr.ns);
                ntmax = (ftau == 1)? 2*(nt/2 +1) : nt;
            }
            if (ntmax % 2) { // odd number
                warn("odd ntmax, reset (%d --> %d) to accommondate trace length (tr.ns=%d)", ntmax, 2*(ntmax/2 + 1), tr.ns);
                ntmax = 2*(ntmax/2 + 1);
            }
            if (ftau == 1) { // input tau-p  
                nt = ntmax - 1;
                if (verbose) warn("reset output nt (%d --> %d) to adapt input ntmax (=%d)", tr.ns, nt, ntmax);
            }
        }
        if(!getparfloat("kxmax", &kxmax)) {
            if(getparfloat("dx", &dx)) kxmax = 500.0/dx;
            else kxmax = 40;
        }
        kxmin = -kxmax;
        if(!getparfloat("kymax", &kymax)) {
            if(getparfloat("dy", &dy)) kymax = 500.0/dy;
            else kymax = 40;
        }
        kymin = -kymax;

        pxmin = tr.fx;
        pymin = (dim > 2)? tr.fy : 0.0;
        dpx = tr.dx;
        dpy = (dim > 2)? tr.dy : 0.0;
        nf = tr.ns/2;
        if (IsComplex) {
            f0 = tr.fz;
            df = tr.dz;
        } else {
            dt = 0.000001*tr.dt;
        }
    }

    // get actual x & y-dimension
    nx = tr.sfs;
    ny = (dim > 2)? tr.sfe : 1;
    if ( nx < 10 ) err("  X dimension of gather (=%d) seems not set correctly", nx);
    if ( ny < 10 && dim == 3) err("  Y dimension of gather (=%d) seems not set correctly", ny);
    if(!getparint("nxmax", &nxmax)) nxmax = nx;
    if( dim > 2 ) { if (!getparint("nymax", &nymax)) nymax = ny; }
    else nymax = 1;
    if (dir > 0) {
        nkx = nx;
        nky = ny;
        kxmax = kxmin + dkx*nkx;
        kymax = kymin + dky*nky;
        if(!getparint("npx", &npx)) npx = nxmax;
        if(!getparint("npy", &npy)) npy = nymax;
        dpx = (pxmax -pxmin)/(npx - 1);
        dpy = (npy > 1)? (pymax -pymin)/(npy - 1) : 0.0;
    } else {
        npx = nx;
        npy = ny;
        if(!getparint("nkx", &nkx)) nkx = nxmax;
        if(!getparint("nky", &nky)) nky = nymax;
        dkx = 2.0*kxmax/nkx;
        dky = (nky > 1)? 2.0*kymax/nky : 0.0;
    }

    if (verbose) warn("Tau-P extension: PX: %6.3f - %6.3f  PY: %6.3f - %6.3f   %d x %d ", pxmin, pxmax, pymin, pymax, npx, npy);
    if (verbose) warn("KX-KY extension: KX: %6.3f - %6.3f  KY: %6.3f - %6.3f   %d x %d ", -kxmax, kxmax, -kymax, kymax, nkx, nky);

    /* determine trace holder max dimensions */
    if ( npx > nxmax ) nxmax = npx;
    if ( npy > nymax ) nymax = npy;
    if ( nkx > nxmax ) nxmax = nkx;
    if ( nky > nymax ) nymax = nky;

    if (verbose) warn("array dimensions used (NT x NX x NY): %d x %d x %d", ntmax, nxmax, nymax);

    if ( dir < 0 && ftau == 1 ) {
        fNyquist = 0.5/dt;
        nf = ntmax/2;
        df = fNyquist/(nf - 1);
    } else if ( dir > 0 && ftau == 1 ) {
        T = 1.0/df;
        dt = T/(nt - 1);
    }

    if(!getparfloat("fmin", &fmin)) fmin = 3.0;
    if(!getparfloat("fmax", &fmax)) fmax = 60.0;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* allocate memory input/output data traces*/
    float*** trdata = ealloc3float(ntmax, nxmax, nymax);
    fftwf_complex *ctrout;  ctrout = (fftwf_complex*) &trdata[0][0][0];
    hdrs2d = (segyhdr **) alloc2(nxmax, nymax, HDRBYTES);
    /* zero out data memory */
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);
    memset(**trdata, 0, nymax*nxmax*ntmax*FSIZE);

    // create plan for transform
    fftwf_plan plan = NULL;
    fftwf_complex *ctr = (fftwf_complex*) trdata[0][0];
    if ( dir < 0 && ftau == 1 ) {  //tau to f before fp2fk
        plan = fftwf_plan_dft_r2c_1d(nt, trdata[0][0], ctr, FFTW_MEASURE);
    } else if ( dir > 0 && ftau == 1 ) {  // f to tau after fk2fp
        plan = fftwf_plan_dft_c2r_1d(nt, ctr, trdata[0][0], FFTW_MEASURE);
    }
    fftscale = 1.0/sqrt((float) (nt));  // fft scaling factor

    /* loop over traces */
    ntr = 0, ntotal = 0, ngather = 0;
    iy = 0, ix = 0; eof = 0;
    do {
        if (nsegy > HDRBYTES) {
            gethval(&tr, index, &valnew);
            IsNewGather = valcmp(type, val, valnew);
        } else eof = 1; //END_OF_FILE

        // cache in input data
        if ( nsegy > HDRBYTES && !IsNewGather ) { /* same key and more data*/
            iy = ntr/nx;
            ix = ntr%nx;
            if (ntr > nxmax*nymax - 1) err("  array dimension too small (%d < %d) traces input", nxmax*nymax, ntr+1);
            if (ix > nxmax - 1) err("  array dimension nxmax too small (%d < %d) traces input", nxmax, ix+1);
            if (iy > nymax - 1) err("  array dimension nymax too small (%d < %d) traces input", nymax, iy+1);
            memcpy(&hdrs2d[iy][ix], &tr, HDRBYTES);
            memcpy(trdata[iy][ix], tr.data, FSIZE*tr.ns);

            if(ntr == 0 && dim > 1) { // first trace, check gather specific parameters
                if ( nx != tr.sfs ) {
                    warn("  Diviate array dimension nx (%d != %d) at %d-th gather (%s=%d)", nx, tr.sfs, ngather+1, key, vtoi(type, valnew));
                    nx = tr.sfs;
                }
                if ( dim == 3 && ny != tr.sfe ) {
                    warn("  Diviate array dimension ny (%d != %d) at %d-th gather (%s=%d)", ny, tr.sfe, ngather+1, key, vtoi(type, valnew));
                }
            }

            ++ntr;
        }

        // do transform
        if ( IsNewGather || eof ) { // new gather or END_OF_FILE
            ++ngather;
            if (verbose > 1) warn("  Transforming %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, valnew));

            ntotal += ntr;

            if ( dir < 0 && ftau == 1) { // convert tau to f
                for (j = 0; j < npy; ++j) {
                    for (i = 0; i < npx; ++i) {
                        ctr = (fftwf_complex*) trdata[j][i];
                        fftwf_execute_dft_r2c(plan, trdata[j][i], ctr);
                    }
                }
            }
            
            // do transform
            fk2fp(dir, nf, f0, df, nkx, kxmin, dkx, nky, kymin, dky, npx, pxmin, dpx, npy, pymin, dpy, trdata);

            // output trace
            nxout = (dir > 0)? npx : nkx;
            nyout = (dir > 0)? npy : nky;
            for (j = 0; j < nyout; ++j) {
                for (i = 0; i < nxout; ++i) {
                    ix = i, iy = j;

                    // set grid parameters of trace header
                    memcpy(&trout, &hdrs2d[iy][ix], HDRBYTES);
                    if ( dir < 0 ) {
                        if ( dkx != 0.0 ) trout.dx = trout.d2 = dkx;
                        if ( dky != 0.0 ) trout.dy = dky;
                        if ( dkx != 0.0 ) trout.fx = trout.f2 = -kxmax + i*dkx;
                        if ( dky != 0.0 ) trout.fy = -kymax + j*dky;
                    } else {
                        if ( dpx != 0.0 ) trout.dx = trout.d2 = dpx;
                        if ( dpy != 0.0 ) trout.dy = dpy;
                        if ( dpx != 0.0 ) trout.fx = trout.f2 = pxmin + i*dpx;
                        if ( dpy != 0.0 ) trout.fy = pymin + j*dpy;
                    }

                    // set common parameters
                    trout.counit = dir*dim;  // set tranform type for next step use
                    if (dim >= 2) {
                        trout.shortpad = nxout;
                        trout.sfs = nxout;
                        trout.nhs = i + 1;
                    }
                    if (dim == 3) {
                        trout.ntr = nyout*nxout;
                        trout.sfe = nyout;
                        trout.nvs = j + 1;
                    }
                    if ( dir < 0 ) {  // inverse, tau-p to f-k
                        trout.trid = FUNPACKNYQ;
                        trout.fz = trout.f1 = 0.0;
                        trout.dz = trout.d1 = df;
                        trout.ns = ntmax;
                    } else if (ftau == 1) { // output tau-p
                        trout.trid = 1;
                        trout.dt = NINT(1000000.0*dt);
                        trout.fz = trout.f1 = 0.0;
                        trout.dz = trout.d1 = dt;
                        trout.ns = nt;
                    } else { // output f-p
                        trout.trid = FUNPACKNYQ;
                        trout.dt = NINT(1000000.0*dt);
                        trout.fz = trout.f1 = 0.0;
                        trout.dz = trout.d1 = df;
                        trout.ns = ntmax;
                    }
                    if ( dir > 0 && ftau == 1) { // convert f to tau
                        ctr = (fftwf_complex*) trdata[iy][ix];
                        fftwf_execute_dft_c2r(plan, ctr, trdata[iy][ix]);
                    }
                    memcpy(trout.data, trdata[iy][ix], FSIZE * trout.ns);
                    if (ftau) for(it=0; it<trout.ns; ++it) trout.data[it] *= fftscale;
                    puttr(&trout);
                }
            }

            val = valnew;
            // reset data cache
            memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);
            memset(**trdata, 0, nymax*nxmax*ntmax*FSIZE);
            ntr = 0, iy = 0, ix = 0;
            continue; // skip reading next trace upon new gather
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (plan) fftwf_destroy_plan(plan);
    free2((void**)hdrs2d);
    free3float(trdata);

    return(CWP_Exit());
}

void fk2fp(int dir, int nf, float f0, float df, int nkx, float kxmin, float dkx, int nky, float kymin, float dky,
	  int npx, float pxmin, float dpx, int npy, float pymin, float dpy, float ***trdata)
{
    complex**   fkslice;
    int         nmax, ix, iy, nx, ny, nynew, nxmax, nymax;

    nx = (dir > 0)? nkx : npx;
    ny = (dir > 0)? nky : npy;
    nynew = (dir > 0)? npy : nky;
    nxmax = MAX(nkx, npx);
    nymax = MAX(nky, npy);

    nmax = MAX(MAX(nkx, nky), MAX(npx, npy));
    if ((fkslice=(complex**)malloc(nmax*sizeof(void*)))==NULL)
        err("  cannot allocate memory for fkslice");
    if (ny > 1) { // 3D
        for (ix=0; ix<nx; ++ix) { // loop over x
            // transform first in y direction
            for (iy=0; iy<nymax; ++iy) fkslice[iy] = (complex*) trdata[iy][ix];
            (dir > 0)? intp_wk2wp(nf, 2*PI*f0, 2*PI*df, nky, 2*PI*kymin, 2*PI*dky, npy, pymin, dpy, fkslice)
                     : intp_wp2wk(nf, 2*PI*f0, 2*PI*df, nky, 2*PI*kymin, 2*PI*dky, npy, pymin, dpy, fkslice);
        }
        ny = nynew;  // change dimension along y after k-p transform
    }
    // transform in x direction
    for (iy=0; iy<ny; ++iy) {
        for (ix=0; ix<nxmax; ++ix) fkslice[ix] = (complex*) trdata[iy][ix];
        (dir > 0)? intp_wk2wp(nf, 2*PI*f0, 2*PI*df, nkx, 2*PI*kxmin, 2*PI*dkx, npx, pxmin, dpx, fkslice)
                 : intp_wp2wk(nf, 2*PI*f0, 2*PI*df, nkx, 2*PI*kxmin, 2*PI*dkx, npx, pxmin, dpx, fkslice);
    }
    free(fkslice);
}
/******************************************************************************

        transform FK to FP via sinc8 interpolation

******************************************************************************/
void intp_wk2wp (int nw, float fw, float dw, int nk, float fk, float dk,
	 int np, float fp, float dp, complex **tr_k)
/******************************************************************************
Input:
nw              number of frequency samples
dw              ANGULAR frequency sampling interval
fw              start/minimum ANGULAR frequency
nk              number of horizontal samples (traces)
dk		sampling interval of k
fk		start/minimum ANGULAR k
np              number of slopes
dp		slope sampling interval
fp              start/minimum slope for tau-p transform

Input and Output:
tr_k          2-D array of input/output complex traces in f-k/f-p domain

modified from SU taup.c. wrap around is done in module SPTAPER. Horizontal phase
shift is implemented in SPFK
******************************************************************************/
{
    int iw, ik, ip; /* loop counters */
    float w, p; /* frequency, slope  */
    complex *tr_ka; /*  */
    complex *hp; /* slant stacked single trace */
    float *kp; /* K-transformed slant stacked single trace */
    complex czero; /* complex number zero */

    czero = cmplx(0.0, 0.0);

    /* allocate working space */
    tr_ka = alloc1complex(nk);
    hp = alloc1complex(np);
    kp = alloc1float(np);

    /* loop over w */
    for (iw = 0, w = fw; iw < nw; ++iw, w += dw) {

        // copy complex input trace along k axis
        for (ik = 0; ik < nk; ik++)
            tr_ka[ik] = tr_k[ik][iw];

        /* compute k values at which to interpolate tr_k */
        for (ip = 0, p = fp; ip < np; ip++, p += dp)
            kp[ip] = w*p;

        /* 8-point sinc interpolation of tr_k to obtain h(p) */
        ints8c(nk, dk, fk, tr_ka, czero, czero, np, kp, hp);

        // copy back to input array
        for (ip = 0; ip < np; ip++) {
            tr_k[ip][iw] = hp[ip];
        }
        // zero out rest of array if there are
        for (ip = np; ip < nk; ip++) {
            tr_k[ip][iw] = czero;
        }
    }

    /* clean up */
    free1complex(tr_ka);
    free1complex(hp);
    free1float(kp);
}

/******************************************************************************

        transform FP to FK via sinc8 interpolation

******************************************************************************/
void intp_wp2wk (int nw, float fw, float dw, int nk, float fk, float dk,
	 int np, float fp, float dp, complex **tr_k)
/******************************************************************************
Input:
nw              number of frequency samples
dw              ANGULAR frequency sampling interval
fw              start/minimum ANGULAR frequency
nk              number of horizontal samples (traces)
dk		sampling interval of k
fk		start/minimum ANGULAR k
np              number of slopes
dp		slope sampling interval
fp              start/minimum slope for tau-p transform

Input and Output:
tr_k          2-D array of input/output complex traces in f-k/f-p domain
******************************************************************************/
{
    int iw, ik, ip; /* loop counters */
    float w, k; /* frequency, k  */
    complex *tr_p; /*  */
    complex *hk; /* uniformly sampled k single trace */
    float *pk; /* p value array where k is evenly sampled */
    complex czero; /* complex number zero */

    czero = cmplx(0.0, 0.0);

    /* allocate working space */
    tr_p = ealloc1complex(np);
    hk = ealloc1complex(nk);
    pk = ealloc1float(nk);

    /* loop over w */
    for (iw = 0, w = fw; iw < nw; ++iw, w += dw) {
        // copy complex input trace along k axis
        for (ip = 0; ip < np; ip++)
            tr_p[ip] = tr_k[ip][iw];

        /* compute p values at which to interpolate k */
        for (ik = 0, k = fk; ik < nk; ik++, k += dk)
            pk[ik] = (w == 0.0) ? 0.0 : k / w;

        /* 8-point sinc interpolation of tr_k to obtain h(p) */
        ints8c(np, dp, fp, tr_p, czero, czero, nk, pk, hk);

        // copy back to input array
        for (ik = 0; ik < nk; ik++) {
            tr_k[ik][iw] = hk[ik];
        }
        // zero out rest of array if there are
        for (ik = nk; ik < np; ik++) {
            tr_k[ik][iw] = czero;
        }
    }

    /* clean up */
    free1complex(tr_p);
    free1complex(hk);
    free1float(pk);
}
