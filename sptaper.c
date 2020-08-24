/* Copyright (c) READ well Service, 2010.*/
/* All rights reserved.                       */

/* SPFK: $Revision: 0.1 $ ; $Date: 20010/03/07 $		*/


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SPTAPER - Wrap around and taper edges  (X - Y) of gather             	",
"                                                                       ",
"    sptaper <infile >outfile  [optional parameters]                 	",
"                                                                       ",
" Optional Parameters:                                                  ",
"   key=cdpt                gather sort key				",
"                                                                       ",
"   dim=3                   (=2) Dimension of gather                    ",
"                                                                       ",
"   wrap=0                  percentage of wrap around added to the edges",
"                           use 2 values for X and Y separately         ",
"                                                                       ",
"   pad=0                   number of times trace at edge is repeated   ",
"                           up to 4 values for all edges separately     ",
"                                                                       ",
"   ntaper=10               number of traces to be tapered at edges     ",
"                           up to 4 values for all edges separately     ",
"   taper=5                 percentages of trace window to be tapered   ",
"                           up to 4 values for all edges separately     ",
"                           ntaper has precedence over taper            ",
"   type=1                  taper type					",
"                           =1 linear (default)                         ",
"                           =2 cosine					",
"                           =3 sine					",
"                           =4 sine^2					",
"                           =5 gaussian (+/-2.0)			",
"                           =6 gaussian (+/-3.8)			",
"                           =n (n>6) even stronger gaussian (+/-n)      ",
"                                                                       ",
" Dimensions of place holder for input/output gather:                   ",
"   nxmax=(tr.sfs)          max. number of inline traces                ",
"   nymax=(tr.sfe)          max. number of x-line traces                ",
"                                                                       ",
" Following headers specify grid coodinates and updated accordingly after wrap around",
"   dx=tr.dx                x spacing                                   ",
"   dy=tr.dy                y spacing                                   ",
"   xmin=tr.fx              starting x of the gather                    ",
"   ymin=tr.fy              starting y                                  ",
"                                                                       ",
"   verbose=0               >0 echo information			",
"                                                                       ",
" NOTE: if gathers have varying dimensions in X or Y, the corresponding ",
"  values should be set >= max plus wrap around (1 + 2*0.01*wrap)*nmax. ",
"  Otherwise dimension of first gather will be used as nmax.            ",
"  Wrap around and taper windows (number of traces) remain constant once",
"  set according to first gather dimensions.                            ",
"                                                                       ",
" Input data must be sorted ascendingly in key, Y, X                    ",
"                                                                       ",
" Run susetgather first to set gather dimensions correctly              ",
"                                                                       ",
NULL};

/*
 * Credits: RWS, Sanyu Ye (sanyu.ye@readgroup.com), March. 2010.
 *
 * Reference:
 *    
 *
 */
/**************** end self doc ********************************/

/* forward declaration of prototype functions */
float* MakeTaper(int n, int type, int verbose);

/** main **/

int main(int argc, char **argv) {
    cwp_String key;     /* header key word from segy.h		*/
    cwp_String type;    /* type of key				*/
    int index;          /* index of key				*/
    Value val, valnew;  /* value of key				*/
    int nsegy;		/* number of bytes read in the segy trace */
    int ntr;            /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;   /* number of total input traces and gathers */
    int nt;             /* number of time or frequency samples */
    int nx, ny;		/* number of horizontal samples */
    int nx1, nx2, ny1, ny2; /* number of traces for taper windows */
    int nw[4];          /* number of traces of wrap around or padding windows */
    int dim, ttype;     /* dimension of input dataset (gather)  and typte of taper */
    int wrap=0, pad=0;  /* wrap and pad flag, exclude each other, not at same time */
    float dt, df;       /* Time/Frequecy sample interval */
    float dx, dy;       /* horizontal sample interval */
    float xmin, xmax, ymin, ymax;   /* min. max. extents of original gather before transform, recorded on first trace  */
    float *sx1, *sx2, *sy1, *sy2;   /* array of scaling factors */
    float *fwrap = NULL, *ftaper = NULL;  // array for holding input window percentages
    int     ntaper[4];    // array holding number of trace to be tapered at edges
    int     ntmax, nxmax, nymax; /* max imensions of input/output trace array */
    int     eof, IsNewGather;
    int     ix, iy, it, i, j;

    int verbose;	 /* flag for echoing information */

    segy tr, trout;

    /* hook up getpar to handle the parameters */
    initargs(argc, argv);
    requestdoc(0);

    if(!getparint("verbose", &verbose)) verbose = 0;

    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    if(!getparint("dim", &dim)) dim = 3;
    if(dim < 2) err( "Tapering is only meaningful for 2 or 3D");
    if(!getparint("type", &ttype)) ttype = 1;
    if(ttype < 1) err( " Unkonwn taper type (=%d)", ttype);

    dt = df = dx = dy = 0.0;
    xmin = xmax = ymin = ymax = 0.0;
    nx1 = nx2 = ny1 = ny2 = nx = ny = 0.0;
    sx1 = sx2 = sy1 = sy2 = NULL;

    ntmax = nt = tr.ns;

    if(!getparfloat("xmin", &xmin)) xmin = tr.fx;
    if(!getparfloat("dx", &dx)) dx = tr.dx;
    if (dim == 3) if(!getparfloat("ymin", &ymin)) ymin = tr.fy;
    if (dim == 3) if(!getparfloat("dy", &dy)) dy = tr.dy;

    // get actual x & y-dimension
    nx = tr.sfs;
    ny = (dim ==3)? tr.sfe : 1;
    if(!getparint("nxmax", &nxmax)) nxmax = 0;
    if(!getparint("nymax", &nymax)) nymax = 1;
    if (nx < nxmax/4) {
        warn("Number of trace along X is not properly set (nx=%d --> %d)", nx, nxmax);
        nx = nxmax;
    }
    if ( nx > nxmax ) {
        nxmax = nx;
        warn("%d-D transform: XMIN=%.1f DX=%.2f  NX=%d  NXMAX=%d ", dim, xmin, dx, nx, nxmax);
    }
    if ( dim == 3 ) { // 3D
        if (ny < nymax/4) {
            warn("Number of trace along Y is not properly set (ny=%d --> %d)", ny, nymax);
            ny = nymax;
        }
        if ( ny > nymax ) {
            nymax = ny;
            warn("3-D transform: YMIN=%.1f DY=%.2f  NY=%d  NYMAX=%d ", ymin, dy, ny, nymax);
        }
    }

    nw[0] = nw[1] = nw[2] = nw[3] = 0;
    // get wrap around parameters if specified
    fwrap = ealloc1float(2);
    memset(fwrap, 0, 2*FSIZE);
    wrap = countparval("wrap");
    if (wrap > dim - 1) err("  Too many wrap around windows (%d > %d)", wrap, dim - 1);
    if (wrap) {
        getparfloat("wrap", fwrap);
        for (i=wrap; i<dim - 1; ++i) fwrap[i] = fwrap[wrap - 1];  // replicate last one specified
        nw[0] = nw[1] = NINT(0.01*fwrap[0]*nxmax);
        if (dim > 2) {
            nw[2] = nw[3] = NINT(0.01*fwrap[1]*nymax);
        }
    }
    if (verbose && wrap) warn(" Wrap around window:  NXW=%d (%.0f%%)    NYW=%d (%.0f%%)", nw[0], fwrap[0], nw[2], fwrap[1]);

    // get padding parameters if specified
    pad = countparval("pad");
    if (pad && wrap) err("  Either pad or wrap around can be specified, but not both");
    if (pad > 2*(dim - 1)) err("  Too many pad windows (%d > %d)", pad, 2*(dim - 1));
    if (pad) {
        getparint("pad", nw);
        for (i=pad; i<2*(dim - 1); ++i) nw[i] = nw[pad - 1];  // replicate last one specified
    }
    if (verbose && pad) warn(" Pad window:  NWX=%d - %d    NWY=%d - %d", nw[0], nw[1], nw[2], nw[3]);

    nxmax += nw[0] + nw[1];
    nymax += nw[2] + nw[3];
    if (verbose) warn("array dimensions used (NT x NX x NY): %d x %d x %d", ntmax, nxmax, nymax);

    // get taper window parameters if specified
    int nc = countparval("ntaper");
    if (nc > 2*(dim - 1)) err("  Too many taper windows (%d > %d)", nc, 2*(dim - 1));
    if (nc > 0) {
        getparint("ntaper", ntaper);
        for (i=nc; i<2*(dim - 1); ++i) ntaper[i] = ntaper[nc - 1];  // replicate last one specified
        nx1 = ntaper[0];
        nx2 = ntaper[1];
        ny1 = ntaper[2];
        ny2 = ntaper[3];
        if (verbose) warn("Constant taper window width used (NX1=%d  NX2=%d   NY1=%d  NY2=%d",
                nx1, nx2, ny1, ny2);
        sx1 = MakeTaper(nx1, ttype, verbose);
        sx2 = MakeTaper(nx2, ttype, verbose);
        sy1 = MakeTaper(ny1, ttype, verbose);
        sy2 = MakeTaper(ny1, ttype, verbose);
    } else {
        ftaper = ealloc1float(4);
        memset(ftaper, 0, 4*FSIZE);
        nc = countparval("taper");
        if (nc == 0) {
            for (i=0; i<2*(dim - 1); ++i) ftaper[i] = 5.0;
        } else if (nc > 2*(dim - 1)) {
            err("  Too many taper windows (%d > %d)", nc, 2*(dim - 1));
        } else {
            getparfloat("taper", ftaper);
            for (i=0; i<nc; ++i) {
                if (ftaper[i] > 10 ) warn("  Probably too larger taper window (=%.0f%%)", ftaper[i]);
            }
            for (i=nc; i<2*(dim - 1); ++i) ftaper[i] = ftaper[nc - 1];  // replicate last one specified
        }
        if (verbose) warn("Percentage taper window width used (NX1=%.0f%%  NX2=%.0f%%    NY1=%.0f%%  NY2=%.0f%%",
                ftaper[0], ftaper[1], ftaper[2], ftaper[3]);
    }


    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* allocate memory input/output data traces*/
    float*** trdata = ealloc3float(ntmax, nxmax, nymax);
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nxmax, nymax, HDRBYTES);
    /* zero out data memory */
    memset(**trdata, 0, nymax*nxmax*ntmax*FSIZE);
    memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);

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
            if(ntr == 0) { // first trace, check gather specific parameters
                if ( xmin != tr.fx ) {
                    warn("  Diviate xmin (%.3f != %.3f) at %d-th gather (%s=%d)", xmin, tr.fx, ngather+1, key, vtoi(type, valnew));
                    xmin = tr.fx;
                }
                if ( nx != tr.sfs ) {
                    warn("  Diviate array dimension nx (%d != %d) at %d-th gather (%s=%d)",
                        nx, tr.sfs, ngather+1, key, vtoi(type, valnew));
                    nx = tr.sfs;
                }
                if (nx + nw[0] + nw[1] < nx1 + nx2)
                    warn("  array dimension nx smaller than taper width (%d < %d) at %d-th gather (%s=%d)",
                        nx + nw[0] + nw[1], nx1 + nx2, ngather+1, key, vtoi(type, valnew));

                if ( dim ==3 && ymin != tr.fy ) {
                    warn("  Diviate ymin (%.3f != %.3f) at %d-th gather (%s=%d)", ymin, tr.fy, ngather+1, key, vtoi(type, valnew));
                    ymin = tr.fy;
                }
                if ( dim == 3 && ny != tr.sfe ) {
                    warn("  Diviate array dimension ny (%d != %d) at %d-th gather (%s=%d)",
                        ny, tr.sfe, ngather+1, key, vtoi(type, valnew));
                    ny = tr.sfe;
                }
                if (dim == 3 && ny + nw[2] + nw[3] < ny1 + ny2)
                    warn("  array dimension ny smaller than taper width (%d < %d) at %d-th gather (%s=%d)",
                        ny + nw[2] + nw[3], ny1 + ny2, ngather+1, key, vtoi(type, valnew));

                if (ftaper) {  // percentage taper used
                    nx1 = NINT(0.01*ftaper[0]*nx);
                    nx2 = NINT(0.01*ftaper[1]*nx);
                    ny1 = (dim < 3)? 0 : NINT(0.01*ftaper[2]*ny);
                    ny2 = (dim < 3)? 0 : NINT(0.01*ftaper[3]*ny);

                    sx1 = MakeTaper(nx1, ttype, verbose);
                    sx2 = MakeTaper(nx2, ttype, verbose);
                    sy1 = MakeTaper(ny1, ttype, verbose);
                    sy2 = MakeTaper(ny1, ttype, verbose);
                }
            }

            iy = ntr/nx;
            ix = ntr%nx;
            if (ntr > nxmax*nymax - 1) err("  array dimension too small (%d < %d) traces input at %d-th gather (%s=%d)",
                      nxmax*nymax, ntr+1, ngather+1, key, vtoi(type, valnew));
            if (ix > nxmax - 1) err("  array dimension nxmax on X too small (%d < %d) traces input at %d-th gather (%s=%d)",
                     nxmax, ix+1, ngather+1, key, vtoi(type, valnew));
            if (iy > nymax - 1) err("  array dimension nymax on Y too small (%d < %d) traces input at %d-th gather (%s=%d)",
                     nymax, iy+1, ngather+1, key, vtoi(type, valnew));
            memcpy(&hdrs2d[iy + nw[2]][ix + nw[0]], &tr, HDRBYTES);
            memcpy(trdata[iy + nw[2]][ix + nw[0]], tr.data, FSIZE*tr.ns);

            ++ntr;
        }

        // do tapering
        if ( IsNewGather || eof ) { // new gather or END_OF_FILE
            ++ngather;
            ntotal += ntr;
            if (verbose > 1) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, valnew));

            //pad wrap around
            if (wrap) {
                for(j=nw[2]; j<ny+nw[2]; ++j) { // first along X
                    for(i=0; i<nw[0]; ++i) {
                        memcpy(trdata[j][i], trdata[j][nx + i], nt*FSIZE);
                        memcpy(trdata[j][i + nx + nw[0]], trdata[j][nw[0] + i], nt*FSIZE);
                    }
                }
                for(i=0; i<nx + nw[0] + nw[1]; ++i) { // then along Y
                    for(j=0; j<nw[2]; ++j) {
                        memcpy(trdata[j][i], trdata[ny + j][i], nt*FSIZE);
                        memcpy(trdata[j + ny + nw[2]][i], trdata[nw[2] + j][i], nt*FSIZE);
                    }
                }
            }

            //pad trace at edge
            if (pad) {
                for(j=nw[2]; j<ny+nw[2]; ++j) { // first along X
                    for(i=0; i<nw[0]; ++i) {
                        memcpy(trdata[j][i], trdata[j][nw[0]], nt*FSIZE);
                    }
                    for(i=0; i<nw[1]; ++i) {
                        memcpy(trdata[j][i + nx + nw[0]], trdata[j][nx + nw[0] - 1], nt*FSIZE);
                    }
                }
                for(i=0; i<nx + nw[0] + nw[1]; ++i) { // then along Y
                    for(j=0; j<nw[2]; ++j) {
                        memcpy(trdata[j][i], trdata[nw[2]][i], nt*FSIZE);
                    }
                    for(j=0; j<nw[3]; ++j) {
                        memcpy(trdata[j + ny + nw[2]][i], trdata[ny + nw[2] - 1][i], nt*FSIZE);
                    }
                }
            }

            // do tapering
            for (it=0; it<nt; ++it) { //loop over time sample
                for(j=0; j<ny+nw[2]+nw[3]; ++j) { // first along X
                    for(i=0; i<nx1; ++i) trdata[j][i][it] *= sx1[i];
                    for(i=0; i<nx2; ++i) trdata[j][nx + nw[0] + nw[1] - i - 1][it] *= sx2[i];
                }
                for(i=0; i<nx+nw[0]+nw[1]; ++i) { // then along Y
                    for(j=0; j<ny1; ++j) trdata[j][i][it] *= sy1[j];
                    for(j=0; j<ny2; ++j) trdata[ny + nw[2 ]+ nw[3] - j - 1][i][it] *= sy2[j];
                }
            }

            // output trace
            for (j = 0; j < ny+nw[2]+nw[3]; ++j) {
                for (i = 0; i < nx+nw[0]+nw[1]; ++i) {
                    // use header on the edge if traces are added by wrap around
                    ix = (i<nw[0])? nw[0] : ( (i>=nx+nw[0])? nx+nw[0]-1 : i );
                    iy = (j<nw[2])? nw[2] : ( (j>=ny+nw[2])? ny+nw[2]-1 : j );

                   // set grid parameters of trace header
                    memcpy(&trout, &hdrs2d[iy][ix], HDRBYTES);
                    trout.dx = trout.d2 = dx;
                    trout.fx = trout.f2 = xmin + (i - nw[0])*dx;
                    if (dim == 3) {
                        trout.dy = dy;
                        trout.fy = ymin + (j - nw[2])*dy;
                    }

                    // set common parameters
                    trout.shortpad = nx+nw[0]+nw[1];
                    trout.sfs = nx+nw[0]+nw[1];
                    trout.nhs = i + 1;
                    if (dim == 3) {
                        trout.ntr = (ny+nw[2]+nw[3])*(nx+nw[0]+nw[1]);
                        trout.sfe = ny+nw[2]+nw[3];
                        trout.nvs = j + 1;
                    }
                    memcpy(trout.data, trdata[j][i], FSIZE * trout.ns);
                    puttr(&trout);
                }
            }

            val = valnew;
            // reset data cache
            memset(*hdrs2d, 0, nxmax*nymax*HDRBYTES);
            memset(**trdata, 0, nymax*nxmax*ntmax*FSIZE);
            ntr = 0, iy = 0, ix = 0;
            if (ftaper) {  // percentage taper used, free memory
                if(sx1) free1float(sx1);
                if(sx2) free1float(sx2);
                if(sy1) free1float(sy1);
                if(sy2) free1float(sy2);
            }
            continue; // skip reading next trace upon new gather
        }
        nsegy = gettr(&tr);
    } while (!eof);

    free2((void**)hdrs2d);
    free3float(trdata);
    
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
