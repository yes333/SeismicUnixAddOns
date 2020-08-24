/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* SULINE2GRID: $Revision: 1.1 $ ; $Date: 2010/04/09 16:17:07 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SULINE2GRID - spread a 2D line of synthetic data to 3D grid	",
" 								",
" suline2grid <data1 >data2 [optional parameters]		",
" 								",
" Required parameters:						",
"   none							",
" 								",
" Optional parameters:						",
"   rotate=0            no rotation for hydrophe & vertical geophone	",
"                       =1 rotate (radial geophone) to output X component",
"                       =2 rotate (radial geophone) to output Y component",
"   dx=12.5             shot interval along x axis of output gather",
"   dy=dx               shot interval along y axis              	",
"   xmax=3000           max. offset along x axis (positive toward north)",
"   xmin=-xmax          min. offset along x axis                        ",
"   ymax=xmax           max. offset along y axis (positive toward east) ",
"   ymin=-ymax          min. offset along y axis                        ",
"   scalco=tr.scalco    scaler for offset,sx,sy,gx,gy (e.g. -100=0.01 )   ",
"   nmax=6000           max. number of input traces                     ",
" 								",
"   verbose=1           verbose = 0 not echo info		",
" 								",
" NOTE:  input data must be sorted by offset            	",
"        the grid takes the nearest trace without interpolation	",
" 								",
" 								",
" Caveat:  offset of input traces must be positive.		",
"   input data is assumed to be shot gather as output of MSEIS  ",
"   output data is in receiver gather with polarity flipped for ",
"   horizontal components                                       ",
NULL};

/* Credits:
 *	Created by: Sanyu Ye, READ Well Service, March, 2008
 *      Last modified:  April 20010
 */
/**************** end self doc ***********************************/


int
main(int argc, char **argv)
{
    int rotate;             /* flag dictating sense of rotation	*/
    int verbose;            /* flag for echoing info		*/
    int trdt;               /* time sample rate as integer	*/
    int nt;                 /* samples per trace on input	*/
    int ntr;                /* traces in input data		*/
    int i, itr, nx ,ny;     /* counter			 	*/
    int nmax;               /* max. number of input gather allowed */
    float scalco;           /* coordinate scale */
    float scale;            /* scaling factor for trace data */
    float x, y;             /* x y coordinate of out put traces */
    float offset;           /* offset of output gather */
    float dxx, x1, x2;      /* trace spacing, min and max offset of input gather  */
    float dx, xmin, xmax;   /* x (inline) dimension of output gather  */
    float dy, ymin, ymax;   /* y (crossline) dimension of output   */
    float **trdata;         /* buffer of trace data	*/
    segy tr;
	
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* Set parameters */
    if (!getparint("verbose", &verbose))	verbose = 1;
    if (!getparint("rotate",  &rotate))    	rotate = 0;
    if ( rotate != 0 && rotate != 1 && rotate != 2 )
            err("rotate = %d, flag must be 0, 1 or 2", rotate);

    if (!getparfloat("xmax", &xmax))	xmax = 3000;
    if (!getparfloat("xmin", &xmin))	xmin = -xmax;
    if (!getparfloat("ymax", &ymax))	ymax = xmax;
    if (!getparfloat("ymin", &ymin))	ymin = -ymax;
    if (!getparfloat("dx", &dx))	dx = 12.5;
    if (!getparfloat("dy", &dy))	dy = dx;

    /* Get info from first trace */
    if (!gettr(&tr))  err("can't get first trace");
    nt = tr.ns;
    trdt = tr.dt;

    if(!getparfloat("scalco", &scalco)) scalco = 0.0;
    /* find out coordinate scaler out of trace header*/
    if ( scalco == 0.0 ) scalco = tr.scalco;
    if ( scalco == 0.0 ) scalco = 1.0;
    else if (scalco < 0.0) scalco = -1.0/scalco;

    x1 = (float) tr.offset * scalco;  // first
    
    /* Allocate data matrices */
    if (!getparint("nmax",  &nmax))     nmax = 6000;
    trdata = ealloc2float(nt, nmax);
    
    ntr = 0;
    do {
        memcpy(trdata[ntr++], tr.data, nt*FSIZE);
        x2 = (float) tr.offset * scalco;  // last offset
    } while (gettr(&tr) && ntr < nmax);
    
    dxx = (x2 - x1)/(ntr - 1);  // input trace spacing
    
    if (verbose) warn(" Total %d traces read in, xmin=%.2f xmax=%.2f dx=%.3f", ntr, x1, x2, dxx);
    
    // set common parameters of output grid
    tr.tracr = 0;
    tr.gx = 0;
    tr.gy = 0;
    tr.dx = tr.d2 = dx;
    tr.dy = dy;
    tr.scalco = NINT( (scalco < 1.0)? -1.0/scalco : scalco );  
    tr.sfs = NINT((xmax - xmin)/dx + 1);
    tr.sfe = NINT((ymax - ymin)/dy + 1);
    tr.lcf = NINT(xmin);
    tr.hcf = NINT(xmax);
    tr.lcs = NINT(ymin);
    tr.hcs = NINT(ymax);
    
    if (verbose) warn("Grid info: %d x %d traces with cell size=%.1f x %.1f  x=%.1f~%.1f y=%.1f~%.1f ", 
            tr.sfs, tr.sfe, dx, dy, xmin, xmax, ymin, ymax);

    
    // loop over grid
    ny = 0;
    for(y=ymin; y<=ymax; y=y+dy) {
        tr.nvs = ++ny;
        tr.sy = NINT(y/scalco);
        nx = 0;
        for (x=xmin; x<=xmax; x=x+dx) {
            offset = sqrt(x*x + y*y);
            // find nearest trace
            itr = (offset < x1)? 0 : NINT(fabs(offset - x1)/dxx);
            if ( itr >= ntr ) {
                warn(" Max. offset exceeds the input one (%.2f > %.2f), program terminated abnormally", offset, x2);
                break;
            }
            
            // input data is actually a shot gather with positive polarity
            if (offset == 0.0) {
                scale = rotate? 0.0 : 1.0;
            }
            else if(rotate == 1) scale = -x / offset; // minus sign because output is receiver gather
            else if(rotate == 2) scale = -y / offset; // flip polarity for positive x or y source coordinate
            else scale = 1.0;
            for (i=0; i<nt; ++i) tr.data[i] = scale*trdata[itr][i];

            tr.fx = tr.f2 = x;
            tr.fy = y;
            tr.tracr++;
            tr.nhs = ++nx;
            tr.sx = NINT(x/scalco);
            tr.offset = NINT(offset/scalco);
            
            if (verbose>1) warn("ny=%d nx=%d  x=%.2f y=%.2f offset=%.2f scale=%.3f ntr=%d ", ny, nx, x, y, offset, scale, itr);

            puttr(&tr);
        }
    }
    
    free2float(trdata);

    return(CWP_Exit());
}

