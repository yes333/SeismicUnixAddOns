/* Copyright (c) READ well Service, 2008.*/
/* All rights reserved.                       */

/* SPTAUP3D: $Revision: 0.1 $ ; $Date: 2008/01/07 $		*/

#include <time.h>
#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>
#include <taup.h>


/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SPTAUP3D - forwared and inverse T-X,Y <--> Tau-Px,Py slant stacks	",
"                                                                       ",
"    sptaup3d <infile >outfile  [optional parameters]                 	",
"                                                                       ",
" Optional Parameters:                                                  ",
" mode=1                     =1 for 2D forward T-X domain computation	",
"                            =2 for 2D inverse T-X domain computation	",
"                            =3 for 3D forward T-X domain computation	",
"                            =4 for 3D inverse T-X domain computation	",
"                            =5 for 3D forward F-K domain computation	",
"                            =6 for 3D inverse F-K domain computation	",
"                                                                   ",
" Dimensions of place holder for input/output traces                   ",
" ntmax=max(nt, ntau)        max. number of t/tau samples of input/output gather ",
" nxmax=401                  max. number of inline traces of input/output gather",
" nymax=1                    max. number of x-line traces of input/output gather",
"                                                                   ",
" For forward transform                                                                  ",
" ntau=tr.ns                 number of time samples in Tau-P domain  ",
" pxmax=1.0                  maximum slope for Tau-P transform (s/km)	",
" pxmin=-pxmax               minimum slope for Tau-P transform (s/km)	",
" npx=401                    number of slopes for Px transform	",
" pymax=pxmax                maximum slope for Tau-P transform (s/km)	",
" pymin=-pymax               minimum slope for Tau-P transform (s/km)	",
" npy=(npx-1)/2 + 1          number of slopes for Py transform	",
"                                                                   ",
" For inverse transform                                                                  ",
" nt=tr.ns                   number of time samples in X-T domain  ",
" dx=12.5                    x spacing	",
" dy=dx                      y spacing	",
" xmax=4500.0                maximum reference x (m)	",
" xmin=-xmax                 minimum reference x (m)	",
" ymax=2000.0                maximum reference y (m)	",
" ymin=-ymax                 minimum reference y (m)	",
"                                                                       ",
" scalco=tr.scalco           scaler for sx/px,sy/py,gx,gy (e.g. -100=0.01 )   ",
" taper=0.0,0.0,0.0,0.0      spatial edge taper, ratios to xmin,xmax,ymin,ymax",
" fmin=3                     minimum frequency of interest 	        ",
" fmax=60                    maximum frequency of interest 	        ",
" ofactor=4                  oversampling in k domain		",
" npoints=71                 number of points for rho filter		",
"								",
" verbose=0                  >0 echo information			",
" key=cdp                    gather sort key				",
" j=1                        process every j-th gather ...           ",
" s=0                        ... based at s  (if ((key - s)%j) == 0) ",
"								",
" for 2D transform, offset is used to compute tau       	",
"         tau = t + px * offset ;  npy = 1, py = 0		",
" for 3D  tau = t + px * (sx - gx) + py * (sy - gy)		",
"                                                       	",
" Notes for pseudo parallel processing:                         ",
" The cascade of multiple (j) modules with different s=(0,1,,,j-1) ",
" enables some kind of parallel processing. Normally j <= c the number of",
" physical CPU cores. ",
"                                                                       ",
NULL};

/*
 * Credits: RWS, Sanyu Ye, Jan. 2008.
 *
 * Reference:
 *    Levin, F., editor, 1991, Slant-Stack Processing, Geophysics Reprint
 *         Series #14, SEG Press, Tulsa.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: dt,d2,f2
 */
/**************** end self doc ********************************/

/** function prototype */
int inv_tx_sltstack(int mode, float dt, int nt, int ntau, int nx, int ny, float xmin, float ymin, float dx, float dy, int npx, float pxmin, float pxmax,
        int npy, float pymin, float pymax, float scalco, float **scale, int npoints, float* rho,  float ***traces, float ***outtraces);
int fwd_tx_sltstack(int mode, float dt, int nt, int ntau, int nxy, int npx, float pxmin, float pxmax,
        int npy, float pymin, float pymax, float scalco, float *taper, segyhdr *hdrs, float **traces, float ***outtraces);
int fwd_sstack(int mode, int nxy, int nt, int ntau, float dt, float px, float py, float scalco, segyhdr *hdrs, float **traces, float *ptrace);
int putgather(int mode, int nt, float dt, int nx, int ny, float xmin, float ymin, float dx, float dy, float x0, float y0, float ***traces, segy *ptr);
float** calc_edge_taper(int npx, int npy, float xtaper, float ytaper);

void fk_trans(float *ctr_p, complex *tr_xy, complex *tr_x, complex *tr_k, complex *tr_ka, complex *hp, complex *ole,
        float *kp, int ny, int nx, int nyfft, int nxfft, float dky, float dkx, float fkay, float fkax,
        int nky, int nkx, float dy, float dx, int npy, int npx, int nw, float fpy, float fpx,
        float dpy, float dpx, float fy, float fx, int iw, float w, int nkay, int nkax, int ntfft,
        float ymin, float xmin, int nxmax);
void fwd_taup_3d_(float *traces, int *ntt, int *nttmax, int *nxxmax, int *nyymax, int *nxx, int *nyy, float *dtt, float *dxx,
        float *dyy, float *xxmin, float *yymin, int *npxx, int *npyy, float *pxxmin, float *pyymin,
        //           float *maxfreqq, int *nxxfft, int *nyyfft, float *work)
        float *pxxmax, float *pyymax, float *maxfreqq, int *nxxfft, int *nyyfft, float *slice, float *work);
void inv_taup_3d_(float *traces, int *ntt, int *nttmax, int *nxxmax, int *nyymax, int *nxx, int *nyy, float *dtt, float *dxx,
        float *dyy, float *xxmin, float *yymin, int *npxx, int *npyy, float *pxxmin, float *pyymin,
        float *pxxmax, float *pyymax, float *maxfreqq, int *nxxfft, int *nyyfft, float *slice, float *work);

const short TRID_TP = 333;


/** main **/

int main(int argc, char **argv) {
    cwp_String key_type;/* type of key				*/
    int key_index;	/* index of key				*/
    cwp_String key;     /* header key word from segy.h		*/
    Value val, val_last;/* value of key				*/
    int nsegy;		/* number of bytes read in the segy trace */
    int ikey;		/* loop counters */
    int ntr;            /* number of actual input traces of the gather just read in*/
    int ngathers, ntotal;   /* number of total input traces and gathers */
    int nt, ntau;	/* number of time samples */
    int nx, ny;		/* number of horizontal samples */
    int npx, npy;	/* number of slopes for slant stack */
    int mode;		/* flag for requested operation */
    float dt;           /* Time sample interval */
    float dx=25.0, dy=25.0;       /* horizontal sample interval */
    float xmin, xmax, ymin, ymax;   /* min. max. extents of original gather before transform, recorded on first trace  */
    float pymin=0.0, pymax=0.0;             /* min. max. slope for Tau-P transform */
    float pxmin=0.0, pxmax=0.0;             /* min. max. slope for Tau-P transform */
    float dpx, dpy;             /* increments for Tau-P transform */
    float fmin, fmax;		/* minimum frequency of interest */
    float gxmin, gxmax, gymin, gymax, x0, y0;	/* actual extent of input gather [m] */
    float scalco;		/* coordinate scale */
    float *rho = NULL;		/* auxiliary array for rho filter */
    float **scale = NULL;		/* 2D edge taper scaler for tau-p traces  */

    int     npoints;		/* number of points for rho filter */
    int     type;  /* trace type marked by trid value */
    int     ntmax, nxmax, nymax; /* max imensions of input/output trace array */
    int     n, n11, n12, n21, n22, n1=0, n2=0, iline, xline;
    int     ntfft, nxfft, nyfft, ofactor;
    float   d1, d2, f1, f2;
    time_t  t0, t1, t2;
    float ***traces2d=NULL, **traces1d=NULL, ***outtraces=NULL;	/* array[nxmax][nymax][ntmax] of input/output traces */
    float *taper;	/* array of edge tapers */
    float *slice;       /* data buffer for t->f fft */
    float *work;        /* data buffer for x->k fft */

    int j;		/* take every jth trace ...		*/
    int s;		/* ... starting at the sth trace ...	*/
    int verbose;	/* flag for echoing information */
    int newGatherFound=0, endGatherFound=0;
    
    segy tr, trout;
    segyhdr **hdrs2d, *hdrs1d;      /* input  header array */

    /* hook up getpar to handle the parameters */
    initargs(argc, argv);
    requestdoc(1);

    t2 = t0 = time(NULL);

    if(!getparint("verbose", &verbose)) verbose = 0;

    if (!getparint("j", &j))		j = 1;
    if (!getparint("s", &s))		s = 0;

    /* get taper's */
    taper=ealloc1float(4);
    if ((n=countparval("taper"))!=0) {
        if (n<5) {
            getparfloat("taper", taper);
            if (n<4) {
                for (ikey=n; ikey<4; ++ikey) {
                    taper[ikey]=taper[n-1];
                }
            }
        }
    } else { /* set default */
        for (ikey=0; ikey<4; ++ikey)
            taper[ikey]=0.0;
    }

    if(!getparint("ofactor", &ofactor)) ofactor = 4;
    if(!getparint("npoints", &npoints)) npoints = 71;

    if(!getparfloat("dx", &dx)) dx = 12.5;
    if(!getparfloat("dy", &dy)) dy = dx;

    /* read first trace */
    nsegy = gettr(&tr);

    if(!getparfloat("scalco", &scalco)) scalco = 0.0;
    /* find out coordinate scaler out of trace header*/
    if ( scalco == 0.0 ) scalco = tr.scalco;
    if ( scalco == 0.0 ) scalco = 1.0;
    else if (scalco < 0.0) scalco = -1.0/scalco;

    if(!getparint("mode", &mode)) mode = 1;
    if(!getparint("nxmax", &nxmax)) nxmax = 201;
    if(!getparint("nymax", &nymax)) nymax = 1;
    if(!getparint("ntmax", &ntmax)) ntmax = tr.ns;
    if(!getparint("nt", &nt)) nt = tr.ns;
    if(!getparint("ntau", &ntau)) ntau = tr.ns;
    if (mode%2) { /* forward transform */
        if(!getparint("npx", &npx)) npx = 401;
        if(!getparint("npy", &npy)) npy = (npx-1)/2 + 1;
    } else {
        pxmin = (float) tr.grnors / 1000.0;
        pxmax = (float) tr.grnofr / 1000.0;
        pymin = (float) tr.grnlof / 1000.0;
        pymax = (float) tr.gaps / 1000.0;
        npx = tr.sfs;
        npy = tr.sfe;
    }

    if(!getparfloat("pxmax", &pxmax)) pxmax = 1.0;
    if(!getparfloat("pxmin", &pxmin)) pxmin = -pxmax;
    if(!getparfloat("pymax", &pymax)) pymax = pxmax;
    if(!getparfloat("pymin", &pymin)) pymin = -pymax;
    if (verbose>=2) warn("Tau-p extension: PX: %.3f - %.03f  PY: %.3f - %.3f   %d x %d ", pxmin, pxmax, pymin, pymax, npx, npy);

    if(!getparfloat("xmax", &xmax)) xmax = 4500.0;
    if(!getparfloat("xmin", &xmin)) xmin = -xmax;
    if(!getparfloat("ymax", &ymax)) ymax = 2000.0;
    if(!getparfloat("ymin", &ymin)) ymin = -ymax;
    nx = NINT((xmax - xmin)/dx) + 1;
    ny = NINT((ymax - ymin)/dy) + 1;
    if (verbose>=2) warn("reference extension: X: %.0f - %.0f  Y: %.0f - %.0f   %d x %d ", xmin, xmax, ymin, ymax, nx, ny);

    /* determine trace holder max dimensions */
    if ( nt > ntmax ) ntmax = nt;
    if ( ntau > ntmax) ntmax = ntau;
    if ( npx > nxmax ) nxmax = npx;
    if ( npy > nymax ) nymax = npy;
    if ( nx > nxmax ) nxmax = nx;
    if ( ny > nymax ) nymax = ny;
    if (verbose>=2) warn("dimensions used: %d x %d x %d", ntmax, nxmax, nymax);

    if(!getparfloat("fmin", &fmin)) fmin = 3.0;
    if(!getparfloat("fmax", &fmax)) fmax = 60.0;

    if(mode ==1 || mode == 2) {  /*2D*/
        nymax = 1;
        npy = 1;
        pymin=pymax=0.0;
    }

    dpx = (npx > 1)? (pxmax - pxmin)/(npx - 1) : 0;
    dpy = (npy > 1)? (pymax - pymin)/(npy - 1) : 0;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdp";
    key_type = hdtype(key);
    key_index = getindex(key);
    gethval(&tr, key_index, &val_last);

    /* allocate memory input/output data traces*/
    traces2d = ealloc3float(ntmax, nxmax, nymax);
    traces1d = *traces2d;
    hdrs2d = (segyhdr **) alloc2(nxmax, nymax, HDRBYTES);
    hdrs1d = (segyhdr*) *hdrs2d;
    if ( traces2d == NULL || hdrs2d==NULL)
        err("memory allocation error for input or/and output data");
    /* zero out data memory */
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);
    memset(**traces2d, 0, nymax*nxmax*ntmax*FSIZE);

    if ( mode == 1 || mode == 3) {
        outtraces = ealloc3float(ntau, npx, npy);
        if ( outtraces == NULL )
            warn("memory allocation error for output data buffer");
    }
    else if ( mode == 2 || mode == 4) {
        outtraces = ealloc3float(nt, nx, ny);
        if ( outtraces == NULL )
            warn("memory allocation error for output data buffer");
    }

    /* estimate buffer size for fft and allocate memory */
    ntfft = npfar(1.5*ntmax);
    nxfft = npfa(ofactor*(nxmax+16));
    nyfft = npfa(ofactor*(nymax+16));

    if (verbose>=2) warn("fft sampling used: ntfft*nxfft*nyfft = %d x %d x %d", ntfft, nxfft, nyfft);

    slice = alloc1float((mode==5)? 4*nxfft*nyfft : 2*ntfft);
    work = alloc1float(4*ntfft + 4*nxfft + 6*(nxfft+nyfft));
    if (slice == NULL || work==NULL)
        err("memory allocation error for data buffer");

    dt = ((float) tr.dt)*1.e-6;

    if ( mode == 2 || mode == 4 ) { //inverse t-x computation
	rho = ealloc1float(npoints);

	/* compute rho filter */
	rho_filter (npoints, ntau, dt, rho);
        
	/* compute scaler */
        scale = calc_edge_taper(npx, npy, taper[0], taper[2]);
    }


    /* loop over traces */
    ntr = 1;
    ntotal = 0;
    ngathers = 0;
    while (nsegy > HDRBYTES) {   /* While previous trace non-empty */
        cwp_Bool passthrough = cwp_false;
        type = tr.trid;
        if (type == TRID_TP && mode%2 ) passthrough = cwp_true;
        if (type != TRID_TP && (mode-1)%2 ) passthrough = cwp_true;
        else {
            int ival = 0;
            /* get key value */
            gethval(&tr, key_index, &val);
            ival = vtoi(key_type, val);
            if ( (ival - s)%j ) passthrough = cwp_true;
        }
        if (passthrough) {
            puttr(&tr);
            nsegy=gettr(&tr);
            continue;
        }

        iline = tr.stas;
        xline = tr.stae;
        if (verbose>=4) warn("ntr = %d  nt = %d  inline = %d  xline = %d", ntr, tr.ns, iline, xline);
        if ( mode < 4 ) {
            memcpy(&hdrs1d[ntr-1], &tr, HDRBYTES);
            memcpy(traces1d[ntr-1], tr.data, FSIZE*tr.ns);
            if (verbose>=4) warn("2 - cdp = %d  nt = %d dt = %d", hdrs1d[ntr-1].cdp, hdrs1d[ntr-1].ns, hdrs1d[ntr-1].dt);
        } else { /* both iline and xline must be set by previous transform */
            if (iline < 1 || xline < 1) // consistancy check
                err("Invalid input iline=%d or/and xline=%d of %d-th gather (tr.stas/tr.stae)!", iline, xline, ngathers+1);
            
            memcpy(&hdrs2d[xline-1][iline-1], &tr, HDRBYTES);
            memcpy( traces2d[xline-1][iline-1], tr.data, FSIZE*tr.ns);
        }
        if (verbose>=4) warn("ntr = %d  nt = %d  inline = %d  xline = %d", ntr, tr.ns, iline, xline);

        if ( ntr == 1 ) {  /* first trace of gather */
            // try to find out the real extent of gather before tau-p transform and use them for inverse
            gxmin = (float) tr.lcf;	/* real min-max in meters */
            gxmax = (float) tr.hcf;
            gymin = (float) tr.lcs;
            gymax = (float) tr.hcs;
            n11 = NINT(gxmin/dx);
            n12 = NINT(gxmax/dx);
            n21 = NINT(gymin/dy);
            n22 = NINT(gymax/dy);
            gxmin = n11*dx;     /* binned to interval */
            if ( gxmin < xmin ) gxmin = xmin;
            gxmax = n12*dx;
            if ( gxmax > xmax ) gxmax = xmax;
            gymin = n21*dy;
            if ( gymin < ymin ) gymin = ymin;
            gymax = n22*dy;
            if ( gymax > ymax ) gymax = ymax;
            n1 = NINT((gxmax - gxmin)/dx) + 1;
            n2 = NINT((gymax - gymin)/dy) + 1;
            if (verbose>=2) warn("extension of %d-th gather %s=%d: X: %.0f - %.0f  Y: %.0f - %.0f  n1 x n2: %d x %d ( %d - %d ) x ( %d - %d)",
                    ngathers+1, key, vtoi(key_type, val), gxmin, gxmax, gymin, gymax, n1, n2, n11, n12, n21, n22);
        }

        endGatherFound = tr.mark;
        if ( endGatherFound ) { // skip reading next trace
            newGatherFound = 0;
        } else {        /*read next trace */
            nsegy = gettr(&tr);
            gethval(&tr, key_index, &val);
            if ( nsegy > HDRBYTES && !valcmp(key_type, val, val_last) ) { /*  new trace read and same key */
                ++ntr;
                val_last = val;
                newGatherFound=0;
                endGatherFound=0;
                continue;
            } else {  
                newGatherFound=1;
            }
        }
        
        if (newGatherFound || endGatherFound) { /* new gather or end of gather, do the transform */
            ntotal += ntr;
            ++ngathers;
            t1 = time(NULL);
            if (verbose) warn("Elapsed time for %d-th gather of total %d traces I/O = %.3f", ngathers, ntr, difftime(t1,t2));

            x0 = (float) hdrs2d[0][0].afilf;
            y0 = (float) hdrs2d[0][0].afils;
            // copy first input header to output trace in order to preserve many common values
            memcpy(&trout, *hdrs2d, HDRBYTES);

            if ( mode == 1 || mode == 3 ) {
                fwd_tx_sltstack(mode, dt, nt, ntau, ntr, npx, pxmin, pxmax, npy, pymin, pymax, scalco, taper, hdrs1d, traces1d, outtraces);
                if (outtraces) memset(**outtraces, 0, npy*npx*ntau*FSIZE); //reset output data cache to zero
            } else if ( mode == 2 || mode == 4 ) {
                scalco = 1.0;
                if (verbose>=2) warn("nt=%d n1=%d n2=%d, npx=%d, npy=%d dt=%f dx=%.1f dy=%.1f xmin=%.1f ymin=%.1f px=%.3f - %.3f py=%.3f - %.3f ",
                        nt, n1, n2, npx, npy, dt, dx, dy, gxmin, gymin, pxmin, pxmax, pymin, pxmax);
                inv_tx_sltstack(mode, dt, ntau, nt, n1, n2, gxmin, gymin, dx, dy, npx, pxmin, pxmax, npy, pymin, pymax,
                        scalco, scale, npoints, rho, traces2d, outtraces);
                putgather(mode, nt, dt, n1, n2, gxmin, gymin, dx, dy, x0, y0, outtraces, &trout);
                memset(**outtraces, 0, ny*nx*nt*FSIZE); //reset output data cache to zero
            } else if ( mode == 5 ) {
                float pxxmin = pxmin;
                float pyymin = pymin;
                float pxxmax= (1-taper[0])*pxxmin;  //muting taper
                float pyymax= (1-taper[2])*pyymin;
                n1 = tr.sfs; //get actual t-x-y gather size 
                n2 = tr.sfe;
                if (verbose>=2) warn("nt=%d n1=%d n2=%d, npx=%d, npy=%d dt=%f dx=%.1f dy=%.1f xmin=%.1f ymin=%.1f px=%.3f - %.3f py=%.3f - %.3f ",
                        nt, n1, n2, npx, npy, dt, dx, dy, gxmin, gymin, pxmin, pxmax, pymin, pxmax);
                fwd_taup_3d_(**traces2d, &nt, &ntmax, &nxmax, &nymax, &nx, &ny, &dt, &dx,
                        &dy, &xmin, &ymin, &npx, &npy, &pxxmin, &pyymin, &pxxmax, &pyymax, &fmax, &nxfft, &nyfft, slice, work);
                if (verbose>=2) warn(" x0=%.1f y0=%.1f pxmin=%.0f pymin=%.0f dpx=%.1f dpy=%.1f npx=%d npy=%d",
                        x0, y0, gxmin, gymin, dx, dy, n1, n2);
                putgather(mode, ntau, dt, npx, npy, pxmin, pymin, dpx, dpy, x0, y0, traces2d, &trout);
            } else if ( mode == 6 ) {
                /* use previous gather dimension */
                float pxxmin = pxmin;
                float pyymin = pymin;
                float pxxmax= (1-taper[0])*pxxmin;  //muting taper
                float pyymax= (1-taper[2])*pyymin;
                if (verbose>=2) warn("nt=%d n1=%d n2=%d, npx=%d, npy=%d dt=%.3f dx=%.1f dy=%.1f xmin=%.1f ymin=%.1f px=%.3f - %.3f py=%.3f - %.3f ",
                        nt, n1, n2, npx, npy, dt, dx, dy, gxmin, gymin, pxmin, pxmax, pymin, pxmax);

                inv_taup_3d_(**traces2d, &nt, &ntmax, &nxmax, &nymax, &npx, &npy, &dt, &dx, &dy, &gxmin, &gymin,
                        &n1, &n2, &pxxmin, &pyymin, &pxxmax, &pyymax, &fmax, &nxfft, &nyfft, slice, work);

                if (verbose>=2) warn(" x0=%.1f y0=%.1f xmin=%.0f ymin=%.0f dx=%.1f dy=%.1f n1=%d n2=%d",
                        x0, y0, gxmin, gymin, dx, dy, n1, n2);
                putgather(mode, nt, dt, n1, n2, gxmin, gymin, dx, dy, x0, y0, traces2d, &trout);
            }

            t2 = time(NULL);
            if (verbose) warn("Elapsed time for %d-th gather transform = %.3f", ngathers, difftime(t2,t1));

            /* zero out data cache */
            memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);
            memset(**traces2d, 0, nymax*nxmax*ntmax*FSIZE);
            val_last = val; /* cache key value */
            ntr = 1;    /* reset gather trace counter */
            
            if (endGatherFound) {
                nsegy=gettr(&tr);  // read next trace
                gethval(&tr, key_index, &val_last);
            }
        }
    }

    t2 = time(NULL);
    if (verbose) warn("Total CPU-time for total %d traces of %d gathers = %.3f", ntotal, ngathers, difftime(t2,t0));

    if (rho)   free1float(rho);
    if (work)  free1float(work);
    if (slice) free1float(slice);
    if (scale) free2float(scale);
    if (outtraces)  free3float(outtraces);
    free1float(taper);
    free2((void**)hdrs2d);
    free3float(traces2d);
    return(CWP_Exit());
}

/******************************************************************************

    Subroutine to output t-x/tau-p trace gather after transform

 ******************************************************************************/
int putgather(int mode, int nt, float dt, int nx, int ny, float xmin, float ymin, float dx, float dy, float x0, float y0, float ***traces, segy *ptr)
/******************************************************************************
Input:
dt              time sampling interval
nt              number of time samples
nx, ny          number of x-y horizontal samples (input traces)
xmin, ymin      minimum x-y coordinates
traces          3-D array of input traces in t-x-y domain

ptr             segy trace with header template input
 ******************************************************************************/
{
    int ix, iy;	/* loop counters */
    float x, y;	/* auxilairy variables */
    float scalco;	/* auxilairy variables */

    /* prepare output t-x trace header, set common parameters */
    ptr->ns = nt; /* number of samples */
    ptr->sfs = nx; /* number of x samples */
    ptr->sfe = ny; /* number of y samples */
    ptr->afilf = NINT(x0);
    ptr->afils = NINT(y0);

    if ( mode % 2 ) { /* forward transform */
        ptr->trid = TRID_TP; /* set trid of tau-p domain */
        scalco = 0.0001;
        ptr->gx = 0;
        ptr->gy = 0;
        ptr->grnors = NINT(1000.0*xmin);
        ptr->grnofr = NINT(1000.0*(xmin+(nx-1)*dx));
        ptr->grnlof = NINT(1000.0*ymin);
        ptr->gaps   = NINT(1000.0*(ymin+(ny-1)*dy));
    } else { /* inverse */
        ptr->trid = 1; /* set trid of tau-p domain */
        scalco = 0.1;
        ptr->gx = NINT(x0/scalco);
        ptr->gy = NINT(y0/scalco);
        ptr->lcf = NINT(xmin);
        ptr->hcf = NINT(xmin+(nx-1)*dx);
        ptr->lcs = NINT(ymin);
        ptr->hcs = NINT(ymin+(ny-1)*dy);
    }
    ptr->scalco = NINT(scalco < 1.0 ? -1.0/scalco : scalco);

    ptr->stae = 0; /* reset x-line no */
    /* loop over traces to output them one by one */
    for (iy=0, y=ymin; iy<ny; y=++iy*dy + ymin) {
        ptr->stas = 0; /* reset inline no */
        ptr->stae++; /* x-line no */
        for (ix=0, x=xmin; ix<nx; x=++ix*dx + xmin) {
            /* set trace individual parameters */
            ptr->stas++; /* inline no */
            ptr->sx = x / scalco ;
            ptr->sy = y / scalco;
            ptr->offset = sqrt(x*x + y*y)/scalco;
            if ( mode % 2 ) { /* forward transform */
                ptr->sx = NINT(x/scalco);  // x=px
                ptr->sy = NINT(y/scalco);  // y=py
            } else { /* inverse */
                ptr->sx = NINT((x0 + x)/scalco);
                ptr->sy = NINT((y0 + y)/scalco);
            }
            memcpy(&(ptr->data[0]), traces[iy][ix], nt*FSIZE);
            puttr(ptr);
        }
    }
    return nx*ny;
}


/******************************************************************************

    Subroutine to compute a forward slant stack (taup transform)
                            in t-x domain

 ******************************************************************************/
int fwd_tx_sltstack(int mode, float dt, int nt, int ntau, int nxy, int npx, float pxmin, float pxmax,
        int npy, float pymin, float pymax, float scalco, float *taper, segyhdr *hdrs, float **traces, float ***outtraces)
        /******************************************************************************
Input:
dt              time sampling interval
nt              number of time samples
ntau            number of time samples in tau-p domain
nxy             number of x-y horizontal samples (input traces)
npx,pxmin,pxmax number, minimum and maximum of px slopes
npy,pymin,pymax number, minimum and maximum of py slopes
scalco          coordinate scale
traces          2-D array of input traces in t-x-y domain

Output:
ptrace          output trace data in tau-px-py domain

Credits:
written by Sanyu Ye, READ Well Service, 2008
         ******************************************************************************/
{
    int ipx, ipy, ixy, it;		/* loop counters */
    float px, py, dpx, dpy;	/* auxilairy variables */
    float xmin, xmax, ymin, ymax, dx, dy, dist, dmin, dmax;
    double dfrac;
    segy tr;

    /* find out range of x-y */
    xmin = xmax = hdrs[0].sx - hdrs[0].gx;
    ymin = ymax = hdrs[0].sy - hdrs[0].gy;
    dmin = dmax = hdrs[0].offset;
    for (ixy=0; ixy<nxy; ++ixy) {
        dx = hdrs[ixy].sx - hdrs[ixy].gx;
        dy = hdrs[ixy].sy - hdrs[ixy].gy;
        dist = hdrs[ixy].offset;
        if ( dx > xmax ) xmax = dx;
        else if ( dx < xmin ) xmin = dx;
        if ( dy > ymax ) ymax = dy;
        else if ( dy < ymin ) ymin = dy;
        if ( dist > dmax ) dmax = dist;
        else if ( dist < dmin ) dmin = dist;
    }
    warn(" minmax=%.0f - %.0f, %.0f - %.0f, %.0f - %.0f", xmin, xmax, ymin, ymax, dmin, dmax);

    for (ixy=0; ixy<nxy; ++ixy) {  /* loop over trace to set taper/weighting factors */
        dx = hdrs[ixy].sx - hdrs[ixy].gx;
        dy = hdrs[ixy].sy - hdrs[ixy].gy;
        dist = hdrs[ixy].offset;
        dfrac = 1.0;
        if (mode==1) { /* 2D */
            if ( taper[0] != 0.0 && dmin != 0.0 && fabs(dmin - dist) < fabs(taper[0]*dmin) ) dfrac = fabs((dmin - dist)/taper[0]*dmin);
            else if ( taper[1] != 0.0 && dmax != 0.0 && fabs(dmax - dist) < fabs(taper[1]*dmax) ) dfrac = fabs((dmax - dist)/taper[1]*dmax);
        }
        else if ( taper[0] !=0.0 && xmin != 0.0 && fabs(xmin - dx) < fabs(taper[0]*xmin) ) dfrac = fabs(xmin - dx)/fabs(taper[0]*xmin);
        else if ( taper[1] !=0.0 && xmax != 0.0 && fabs(xmax - dx) < fabs(taper[1]*xmax) ) dfrac = fabs(xmax - dx)/fabs(taper[1]*xmax);
        else if ( taper[2] !=0.0 && ymin != 0.0 && fabs(ymin - dx) < fabs(taper[2]*ymin) ) dfrac = fabs(ymin - dy)/fabs(taper[2]*ymin);
        else if ( taper[3] !=0.0 && ymax != 0.0 && fabs(ymax - dx) < fabs(taper[3]*ymax) ) dfrac = fabs(ymax - dy)/fabs(taper[3]*ymax);

        hdrs[ixy].unscale = (float) dfrac;
    }

    dpx = (npx == 1)? 0.5 * (pxmax - pxmin) : (pxmax - pxmin) / (npx - 1);
    if ( mode == 1 ) { /* 2D */
        pymin = 0.0;
        pymax = 0.0;
        dpy = 0.0;
        npy = 1;
    }
    else if (npy == 1) { pymin = 0.5 * (pymax - pymin); dpy = 0; }
    else dpy = (pymax - pymin) / (npy - 1);

    /* prepare output tau-p trace header, set common parameters */
    memcpy(&tr, &hdrs[0], HDRBYTES);
    tr.trid = TRID_TP; /* set trid of tau-p domain */
    /*hdr.ntr = npx*npy;*/
    tr.f1 = 0.0;
    tr.scalco = -10000;  /* fixed */
    tr.sfs = (short) npx;
    tr.sfe = (short) npy;
    tr.afilf = (short) (((float)tr.gx)*scalco);  /* x origin */
    tr.afils = (short) (((float)tr.gy)*scalco);  /* y origin */
    tr.stas = 0;
    tr.stae = 0;
    tr.nofilf = (short) (dmin*scalco);
    tr.nofils = (short) (dmax*scalco);
    tr.lcf = (short) (xmin*scalco);
    tr.hcf = (short) (xmax*scalco);
    tr.lcs = (short) (ymin*scalco);
    tr.hcs = (short) (ymax*scalco);

    tr.grnors = 1000.0*pxmin;
    tr.grnofr = 1000.0*pxmax;
    tr.grnlof = 1000.0*pymin;
    tr.gaps   = 1000.0*pymax;

    tr.nhs = 0; /* inline no */
    tr.nvs = 0; /* crossline no */
    tr.gx = 0; /* zero out */
    tr.gy = 0; /* receiver coordinates */


    warn(" mode=%d dt=%f nt=%d ntau=%d nxy=%d npx=%d npy=%d px=%f - %f dpx=%f py=%f -%f dpy=%f scalco=%f",
            mode, dt, nt, ntau, nxy, npx, npy, pxmin, pxmax, dpx, pymin, pymax, dpy, scalco);

    /* loop over slopes */
    for (ipy=0, py=pymin; ipy<npy; py=pymin+(++ipy)*dpy) {
        ++tr.stae;  /* x-line no */
        tr.stas = 0;  /* inline no reset */
        for (ipx=0, px=pxmin; ipx<npx; px=pxmin+(++ipx)*dpx) {
            /* loop over traces */
            fwd_sstack(mode, nxy, nt, ntau, dt, px, py, scalco, hdrs, traces, tr.data);

            /* output trace */
            if (outtraces) {
                memcpy(outtraces[ipy][ipx], tr.data, ntau*FSIZE);
            }
            else {
                /* set trace individual parameters */
                tr.sx = (int) (10000.0 * px);
                tr.sy = (int) (10000.0 * py);
                tr.offset = (int) (10000.0 * sqrt(px*px + py*py));
                ++tr.stas; /* inline no */
                puttr(&tr);
            }
        }
    }
    if ( !outtraces) return EXIT_SUCCESS;

    /* loop over slopes to out put trace from buffer */
    tr.stae = 0;  /* x-line no */
    for (ipy=0, py=pymin; ipy<npy; py=pymin+(++ipy)*dpy) {
        ++tr.stae;  /* x-line no */
        tr.stas = 0;  /* inline no reset */
        for (ipx=0, px=pxmin; ipx<npx; px=pxmin+(++ipx)*dpx) {
            /* set trace individual parameters */
            tr.sx = (int) (10000.0 * px);
            tr.sy = (int) (10000.0 * py);
            tr.offset = (int) (10000.0 * sqrt(px*px + py*py));
            ++tr.stas; /* inline no */
            memcpy(tr.data, outtraces[ipy][ipx], ntau*FSIZE);
            puttr(&tr);
        }
    }
    return EXIT_SUCCESS;
}
/******************************************************************************

    Subroutine to compute a inverse slant stack (taup transform)
                            in t-x domain

 ******************************************************************************/
int inv_tx_sltstack(int mode, float dt, int nt, int ntau, int nx, int ny, float xmin, float ymin, float dx, float dy, int npx, float pxmin, float pxmax,
        int npy, float pymin, float pymax, float scalco, float **scale, int npoints, float* rho, float ***traces, float ***outtraces)
        /******************************************************************************
Input:
dt              time sampling interval
nt              number of time samples of input traces
ntau            number of time samples of output traces
dx, dy          sampling interval in meter of x-y horizontal samples (output traces)
nx, ny          number of x-y horizontal samples (output traces)
npx,pxmin,pxmax number, minimum and maximum of px slopes
npy,pymin,pymax number, minimum and maximum of py slopes
scalco          coordinate scale
scale           2D array containing weighting factors due to edge tapering
traces          2-D array of input traces in tau-p domain

Output:
outtraces       output trace data in x-t domain

Credits:
written by Sanyu Ye, READ Well Service, 2008
         ******************************************************************************/
{
    int fit, lit;			/* first and last time samples */
    int ipx, ipy, ix, iy, it, id;	/* loop counters */
    int np2=npoints/2;	/* half number of points in rho filter */
    float px, py, dpx, dpy, x, y;	/* auxilairy variables */
    float dfrac, delay;
    float *rho_trace; 
    
    dpx = (npx == 1)? 0.0 : (pxmax - pxmin) / (npx - 1);
    dpy = (npy == 1)? 0.0 : (pymax - pymin) / (npy - 1);
    if ( mode == 2 ) { /* 2D */
        ymin = 0.0;
        dy = 0.0;
        ny = 1;
    }
  
    // conv input traces with rho filter
    rho_trace = ealloc1float(nt);
    for (ipy=0; ipy<npy; ++ipy) {
        for (ipx=0; ipx<npx; ++ipx) {
            convolve_cwp(nt,-np2,traces[ipy][ipx],npoints,0,rho,nt,0,rho_trace);
            memcpy(traces[ipy][ipx], rho_trace, nt*FSIZE); //copy back
        }
    }
    free1float(rho_trace);
    
    /* loop over x-y */
    for (iy=0, y=ymin; iy<ny; y=ymin+(++iy)*dy) {
        for (ix=0, x=xmin; ix<nx; x=xmin+(++ix)*dx) {
            /* loop over traces */
            for (ipy=0, py=pymin; ipy<npy; py=pymin+(++ipy)*dpy) {
                for (ipx=0, px=pxmin; ipx<npx; px=pxmin+(++ipx)*dpx) {
                    /* compute two point interpolator, note p in s/km hence factor 0.001 */
                    delay = -0.001*(px*x + py*y)*scalco/dt;
                    if (delay>=0) {
                        id = (int)delay;
                        fit = id+1;
                        lit = (ntau + id < nt)? ntau + id - 1 : nt-1;
                    } else {
                        id = (int)delay-1;
                        fit = 1;
                        lit = (ntau - id < nt)? ntau + id : nt+id;
                    }
                    dfrac = delay-id;

                    /* compute the actual slant stack */
                    for (it=fit; it<lit; it++) {
                        outtraces[iy][ix][it-id] += scale[ipy][ipx]*(traces[ipy][ipx][it] + dfrac*(traces[ipy][ipx][it+1]-traces[ipy][ipx][it]));
                    }
                }
            }
        }
    }
    return EXIT_SUCCESS;
}

/******************************************************************************/
int fwd_sstack(int mode, int nxy, int nt, int ntau, float dt, float px, float py, float scalco, segyhdr *hdrs, float **traces, float *ptrace)
/******************************************************************************/
{
    /******************************************************************************
    Input:
    dt              time sampling interval
    nt              number of time samples
    ntau            number of time samples in tau-p domain
    nxy             number of x-y horizontal samples (input traces)
    px 			    px slopes
    py              py slopes
    scalco          coordinate scale
    traces          2-D array of input traces in t-x-y domain

    Output:
    ptrace          output trace data in tau-px-py domain
     ******************************************************************************/
    int fit, lit;		/* first and last time samples */
    int id, ixy, it;                 /* auxiliary variables */
    float dfrac, delay, *trace;	/* more auxiliary variables */

    /* zero out output trace */
    memset(ptrace, 0, ntau*FSIZE);

    for (ixy=0; ixy<nxy; ++ixy) {
        /* compute two point interpolator, note p in s/km hence factor 0.001 */
        delay = (mode == 1)? 0.001*(px*((float)hdrs[ixy].offset))*scalco/dt : \
            0.001*(px*((float)(hdrs[ixy].sx - hdrs[ixy].gx)) + py*((float)(hdrs[ixy].sy - hdrs[ixy].gy)))*scalco/dt;
        if (delay>=0) {
            id = (int)delay;
            fit = id+1;
            lit = (ntau + id < nt)? ntau + id - 1 : nt-1;
        } else {
            id = (int)delay-1;
            fit = 1;
            lit = (ntau - id < nt)? ntau + id : nt+id;
        }
        dfrac = delay-id;

        trace = traces[ixy];
        /* compute the actual slant stack */
        for (it=fit; it<lit; it++) {
            /* TODO: areal weighting factor by triangulation */
            ptrace[it-id] += hdrs[ixy].unscale*(trace[it] + dfrac*(trace[it+1]-trace[it]));
        }
    }
    return EXIT_SUCCESS;
}

/*
 */
float** calc_edge_taper(int npx, int npy, float xtaper, float ytaper) {
    int ix, iy, nxtaper, nytaper;
    float dfrac;
    float** scale;

    scale = ealloc2float(npx, npy);
    // assumming symetric tau-p gather
    nxtaper = (npx == 1)? 0 : NINT(0.5 * xtaper * (float)npx);
    nytaper = (npy == 1)? 0 : NINT(0.5 * ytaper * (float)npy);

    for (iy=0; iy<npy; ++iy) {  /* loop over trace to set taper/weighting factors */
        for (ix=0; ix<npx; ++ix) {  /* loop over trace to set taper/weighting factors */
            dfrac = 1.0;
            if ( nxtaper > 0 && ix < nxtaper ) dfrac = (float)(ix)/(float)(nxtaper);
            else if ( nxtaper > 0 && (npx - ix + 1) < nxtaper ) dfrac = (float)(npx -ix + 1)/(float)(nxtaper);
            else if ( nytaper > 0 && iy < nytaper ) dfrac = (float)(iy)/(float)(nytaper);
            else if ( nytaper > 0 && (npy - iy + 1) < nytaper ) dfrac = (float)(npy -iy + 1)/(float)(nytaper);

            scale[iy][ix] = dfrac;
        }
    }
    return scale;
}

#include "fwd_taup_3d.c"
#include "inv_taup_3d.c""
