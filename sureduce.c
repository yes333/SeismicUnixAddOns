/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* SUREDUCE: $Revision: 1.11 $ ; $Date: 2005/10/04 16:42:43 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SUREDUCE - shift traces uniformly in reduced time, linear or hyperbolic       ",
"                                                                               ",
" sureduce <stdin >stdout                                                       ",
"                                                                               ",
" Optional parameters:                                                          ",
"   rv=1.480            reduction velocity in km/sec                            ",
"                       <0 remove time shift                                    ",
"                                                                               ",
"   mode=2              =0 hyperbolic, nmo-like reduction                       ",
"                       =1 linear as for refraction seismic                     ",
"                       =2 direct water wave arrivals                           ",
"                       =3 first multiple of water wave arrivals                ",
"                       =4 tau-p of first water wave arrivals                   ",
"                                                                               ",
"   t0=1.0              time origin used to compute hyperbolic reduction        ",
"   scale=tr.scalco     offset scaler                                           ",
"                                                                               ",
"                                                                               ",
"                                                                               ",
" Note:                                                                         ",
" time shift is computed by:                                                    ",
"   mode=0  tshift = sqrt((scale*offset)**2/rv**2 + t0**2) - t0                 ",
"       =1  tshift = abs(scale*offset)/rv - t0                                  ",
"       =2  tshift = sqrt(offset**2 + gwdep**2)*scale/rv - t0                   ",
"       =3  tshift = 3*sqrt((offset/3)**2 + gwdep**2)*scale/rv - t0             ",
"       =4  tshift = gwdep**scale*sqrt((1/rv)**2 - p**2) - t0  for p <= 1/rv    ",
"                                                                               ",
" Trace header fields accessed: dt, ns, offset, scalco                          ",
" Trace header fields modified: none                                            ",
NULL};

/*
 * Author: UC Davis: Mike Begnaud  March 1995
 * extended by Sanyu Ye, READ Well Service, Jan. 2008
 *
 * Trace header fields accessed: ns, dt, offset
 */
/**************** end self doc ***********************************/
segy tr;

int
main(int argc, char **argv)
{
	int nt, mode, remove=0, i, rnt;
	float dt, rv, t0, offset, wdepth, bt, scale, px, py, p2;

	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Get parameters */
	if (!getparfloat("rv", &rv))	 	rv = 1.48;
	if (!getparfloat("t0", &t0))	 	t0 = 1.0;
	if (!getparfloat("scale", &scale))	scale = 0.0;
	if (!getparint("mode", &mode))	 	mode = 2;

	if (rv < 0 ) remove=1;
        
	/* Get info from first trace */
	if (!fgettr(stdin, &tr)) err("can't read first trace");
	if (!tr.dt) MUSTGETPARFLOAT("dt", &dt);
	nt = (int) tr.ns;
	dt = ((double) tr.dt)/1000000.0;
	if (scale == 0.0) scale = tr.scalco;
	if (scale == 0.0) scale = 1.0;
        else if(scale <  0.0) scale = -1.0/scale;

        float* data = ealloc1float(nt);
        memset(data, 0, nt*FSIZE);

	/* Loop over traces */
	do {
		offset = 0.001*scale*ABS(tr.offset);
                wdepth = 0.001*scale*ABS(tr.gwdep);
                if (mode == 0) {
                    bt = sqrt(offset*offset/(rv*rv) + t0*t0) - t0;
                } else if (mode == 1) {
                    bt = offset/ABS(rv) - t0;
                } else if (mode == 2) { // mode=2
                    bt = sqrt(offset*offset + wdepth*wdepth)/ABS(rv) - t0;
                } else if ( mode == 3) { // mode=3
                    bt  = 3.0*sqrt(offset*offset/9.0 + wdepth*wdepth)/ABS(rv) - t0;
                } else if ( mode == 4) {
                    px = tr.fx;
                    py = tr.fy;
                    p2 = px*px + py*py;
                    if ( p2 <= 1.0/(rv*rv) ) {
                        bt  = wdepth*sqrt(1.0/(rv*rv) - p2) - t0;
                    } else {
                        bt = - t0;
                    }
                } else {
                    err(" Invalid mode =%d)", mode);
                }
                rnt    = NINT(bt/dt);

                if ( !remove ) {  /*move upwards */
                    for (i = 0; i < nt; ++i) {
                            register int j = i + rnt;
                            data[i] = (0 <= j && j < nt) ?  tr.data[j] : 0.0;
                    }
                }
                else {  /* move downwards */
                    for (i = nt-1; i >= 0; --i) {
                            register int j = i - rnt;
                            data[i] = (0 <= j && j < nt) ?  tr.data[j] : 0.0;
                    }
                }
                memcpy(tr.data, data, nt*FSIZE);
		puttr(&tr);
                memset(data, 0, nt*FSIZE);
	} while (gettr(&tr));
	
	return(CWP_Exit());
}

