/***************************************************************************
 *   Copyright (C) 2009 by Sanyu Ye   *
 *   sanyu.ye@readgroup.com   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SUNMO3 - NMO for an arbitrary velocity function at a pseudo grid",
"									     ",
"   sunmo3 <stdin >stdout [optional parameters]			     ",
"									     ",
" There are three ways to input time-velocity function		     ",
"  Usually via ascii table similar to sushw			     ",
"									     ",
"   key=            =[[k3,]k2] is the list of header keys as they appear in vfile",
"   vfile=          ascii table of values for fields specified by	",
"                     [k3 [k2]] tnmo vnmo		",
"                     ...               ...				",
" 									",
"  If regularly sampled, t-v function can be also input as velocity cube      ",
"   key=            =[[k3,]k2]     keys for third/slowest and second/fast dimensions",
"                                  first/fastest dimension is always the time   ",
"   n=              =n1[,n2[,n3]] dimensions of velocity cube",
"   d=              =d1[,d2[,d3]] sampling intervals of velocity cube",
"   o=              =o1[,o2[,o3]] origins of velocity cube",
"   vfile=          binary file containing velocity cube in native floating point",
" 									",
" 									",
"  Alterntively t-v function can be input in the same way as sunmo (2D only)",
" 									",
"   key=            header word that denotes receiver number                     ",
"   cdp=            CDPs (receivers) for which vnmo & tnmo are specified",
"   tv=             arrays of trms/vrms pairs, a replacement for seperate  ",
"                     tnmo/vnmo arrays as below ",
"   tnmo=           times corresponding to velocities in vnmo	     ",
"   vnmo=           rms velocities corresponding to times in tnmo	     ",
"									     ",
" More optional Parameters:						     ",
" 									     ",
"   maxlines=512    max. number of lines in ascii vfile                       ",
"   max=2.0         max. scaling factor computed ",
"                      (2.0 corresponds incident angle of 60 degree) ",
"   scale=1.0,...   scaling factors for header values specified by key=k3,k2 ",
"   scalof=scalco   scaling factor for offset, default tr.scalco             ",
"                                                                            ",
"   smute=1.5       samples with NMO stretch exceeding smute are zeroed  ",
"   lmute=25        length (in samples) of linear ramp for stretch mute  ",
"   sscale=1        =1 to divide output samples by NMO stretch factor    ",
"   invert=0        =1 to perform (approximate) inverse NMO	             ",
"   upward=0        =1 to scan upward to find first sample to kill       ",
"                                                                            ",
" Notes:								     ",
"									     ",
" This is sunmo extended with 3D velocity input, however without anisotropy",
"									     ",
"  If t-v function is given via an ascii table, it must be sorted ascendingly ",
"  by k3, k2, trms, i.e. the slowest (e.g. crossline), fast (inline/cdp), and",
"  fastest (time) dimension. ",
"  If they are regularly sampled, it can be read as a binary file, which stores",
"  velocity values in the sequence of fastest, fast and slowest dimensions (as",
"  e.g. k3 k2 sorted su file with header-stripped, or seplib .H@ file).       ",
"  If parameter n is specified, vfile given is automatically regarded as binary",
"  floating point velocity file. ",
"  If t-v function is input via multiple arrays of tv or tnmo/vnmo, key values",
"  specified by cdp= must be sorted ascendingly, and accordingly the order of ",
"  appearance of tv or tnmo/vnmo arrays",
"",
" Caveat: ",
"",
"  If given in 3D, velocities at any point are first linearly interpolated at ",
"  the intersecting points along two neighboring inlines, then again linearly",
"  interpolated between these tow points. Velocity functions should be sampled",
"  densely enough along inlines to avoid strange behavior of the interpolation",
"",
"  Since float is used to represent all key values in vfile, they should not",
"  exceed a precision of 6-7 digits.",
"",
" Examples: ",
"",
" sunmo3 key=ep,cdp vfile=ep-cdp-t-v.table maxlines=2000 <in.su >out.su   ",
"",
" sunmo3 key=gy,gx scale=0.1,0.1 vfile=vcube.rsf         \\ ",
"        n=1001,20,10 d=0.008,500,500 o=0,3000,2000 <in.su >out.su   ",
"",
" sunmo3 tv=0.400,1500,1.100,1834,1.900,2140,2.567,2393,3.138,2630 \\ .  ",
"        scalof=1.0 <in.su >out.su   ",
"",
" sunmo3 key=cdp cdp=1,240 <in.su >out.su  \\",
"        tnmo=0.400,1.100,1.900,2.567,3.138  vnmo=1500,1834,2140,2393,2630 \\",
"        tnmo=0.600,1.300,2.100,2.867,3.438  vnmo=1500,1734,2040,2493,2730   ",
"            ",
"									     ",
NULL};

/* Credits:
    *	RWS: Sanyu Ye, 2009, sanyu.ye@readgroup.com
    *
    * Trace header fields accessed: ns, dt, delrt, offset, cdp, sy
*/
/**************** end self doc *******************************************/

int FindInterval(int ikey, float dv, float *dprev, float **keys, int *order,
                int nRows, int *nbegp,  int* n1p, int* n2p, float* d1p, float* d2p);

void interp1(int ikey, int n1, int n2, float d, float **keys, float **var, int nv, float* v);

void interp2(float d, float d1, float d2, float *v1, float *v2, int nv, float* v);

int interp(int nkeys, float *dv, float **keys, int *sort, int nRows, float **var, int nv, float* v);


float xtscale(float x, float t, float depth, float* vrms, float ft, float dt, float fmax);

segy tr;

int verbose = 0;

int
        main(int argc, char **argv)
{
    cwp_String key[3];  /* array of keywords			*/
    cwp_String type[3]; /* array of keywords			*/
    int index[3];       /* name of type	of getparred key	*/
    Value val[3], valast[3];      /* ... its value			*/
    float fkeys[3];     /* key values for which vnmo sought through interpolation */
    float fscale[3];     /* scaling factors for header keys */
    float **keys = NULL;  /* pointer containing values of matching keys */
    float **v = NULL;   /* pointer containing velocity values */
    float scalof;       /* scaling factors for header */

    int maxRows = 0;        /* max. number of lines of ascii infile */
    int nRow = 0;           /* actual number of lines of ascii infile read in*/

    int ikey;           /* key counter 				*/
    int nkeys;          /* number of header fields set		*/
    int nt;		/* number of time samples per trace */
    float dt;           /* time sampling interval */
    float ft, t;        /* time of first sample and later */
    int it;             /* time sample index */
    int ncdp;           /* number of cdps specified */
    float *cdp;         /* array[ncdp] of cdps */
    int icdp;           /* index into cdp array */
    int nnmo;           /* number of t-v arrays */
    int ntv, ntvnmo;    /* number of tvnmos specified */
    float *tvnmo;       /* array[ntvnmo] of tvnmos */
    int nvnmo;          /* number of vnmos specified */
    float *vnmo;        /* array[nvnmo] of vnmos */
    int ntnmo;          /* number of tnmos specified */
    float *tnmo;        /* array[ntnmo] of tnmos */
    int i, j;           /* index used in loop */
    int n[3];           /* dimensions of velocity cube */
    float o[3];            /* origin of velocity cube  */
    float d[3];            /* sampling intervals of velocity cube */
    char *vfile="";     /* name of input file of header values and time-velocity	*/
    FILE *vfp=NULL;     /* pointer to input vfile		*/
    int nn=0;           /*number of dimension of velocity cube */
    float smute;	/* zero samples with NMO stretch exceeding smute */
    float osmute;	/* 1/smute */
    int lmute;          /* length in samples of linear ramp for mute */
    int itmute=0;       /* zero samples with indices less than itmute */
    int sscale;         /* if non-zero, apply NMO stretch scaling */
    int invert;         /* if non-zero, do inverse NMO */
    int upward;         /* scans upward if it's nonzero. */
    float tn;           /* NMO time (time after NMO correction) */
    float *qtn;         /* NMO-corrected trace q(tn) */
    float *ttn;         /* time t(tn) for NMO */
    float *atn;         /* amplitude a(tn) for NMO */
    float *qt;          /* inverse NMO-corrected trace q(t) */
    float *tnt;         /* time tn(t) for inverse NMO */
    float *at;          /* amplitude a(t) for inverse NMO */
    float tsq;          /* temporary float */

    /* hook up getpar */
    initargs(argc, argv);
    requestdoc(1);

    /* get information from the first header */
    if (!gettr(&tr)) err("can't get first trace");
    nt = tr.ns;
    dt = ((float) tr.dt)/1000000.0;
    ft = tr.delrt/1000.0;

    if (!getparint("maxlines", &maxRows))	maxRows = 512;
    if (!getparint("verbose", &verbose))	verbose = 0;

    if ((nkeys=countparval("key")) > 0) {
        getparstringarray("key", key);
        /* get types and indexes corresponding to the keys */
        for (ikey=0; ikey<nkeys; ++ikey) {
            type[ikey]=hdtype(key[ikey]);
            index[ikey]=getindex(key[ikey]);
        }
    }

// first read in time-velocity function
    if (getparstring("vfile", &vfile)) { // velocity file specified
        if((vfp=efopen(vfile, "r"))==NULL) err("cannot open vfile=%s\n", vfile);
        if ( (nn = countparval("n")) > 0 ) { // dimension of v cube given, binary vfile
            if (nn != nkeys+1)
                err("Number of keys (=%d) must one less than that of dimension (=%d)", nkeys, nn);
            if ( countparval("d") != nn ) err("number of d values must be the same (n=%d)", nn);
            if ( countparval("o") != nn ) err("number of o values must be the same (n=%d)", nn);
            getparint("n", n);
            getparfloat("d", d);
            getparfloat("o", o);
            for(i=0, maxRows=1; i<nn; ++i) maxRows *= n[i];
            keys = ealloc2float(nn, maxRows);
            v = ealloc2float(1, maxRows);
            if (fread(&v[0][0], sizeof(float), maxRows, vfp) != maxRows)
                err("error reading vfile=%s\n", vfile);
            int ng1 = 1, ng2 = 1;
        // populate keys and time
            for(i = 1; i <= nn; ++i) {
                ng1 *= n[nn - i ];
                ng2 *= ( i > 1 )? n[nn - i + 1] : 1;
                for(nRow=0; nRow<maxRows; ++nRow) {
                    keys[nRow][nn-i] = o[nn-i] + ((int)((nRow)%ng1)/ng2)*d[nn-i];
                }
            }
        } else {/* reading from a acsii file/table or value array*/
            int bEnd = 0;

            keys = ealloc2float(nkeys + 1, maxRows);
            v = ealloc2float(1, maxRows);  // currently only one value, 2D array due to compatibility with sushw routings
            if (verbose) warn("number of keys including tnmo m=%d file=%s" , nkeys+1, vfile);
            // reading all data from ascii infile */
        //  sscanf returns: 0   : characters there, but no conversion (error)
        //		  EOF : eof before conversion
        //		  else: number of conversions
            //
            for (nRow=0; nRow < maxRows; ++nRow) {
                for (ikey = 0; ikey < nkeys; ++ikey ) {
                    if ( fscanf(vfp, "%f", &keys[nRow][ikey]) == EOF) {
                        bEnd = 1; break;  /* else everything is okay: get out of the loop */
                    }
                }
                if (bEnd) break;
                if ( fscanf(vfp, "%f", &keys[nRow][nkeys]) == EOF ) break;
                if ( fscanf(vfp, "%f", &v[nRow][0]) == EOF ) break;
            }
        }
    } else { //t-v function read in via tv or tnms/vnms arrays
        ncdp = countparval("cdp");
        nnmo = (ncdp == 0)? 1 : (nkeys > 1)? ncdp / nkeys : ncdp;  // number of t-v function
        ntv = countparname("tv");
        nvnmo = countparname("vnmo");
        ntnmo = countparname("tnmo");
        if (nnmo > 1) { // more than one t-v functions
            if (nkeys == 0) err("Header key word must be given by key=[[k3,]k2]");
            else if ((ncdp%nkeys)) err("number of key values (=%d) must be multiple of number of keys (=%d)", ncdp, nkeys);
            if (ntv > 0 && ntv != nnmo) err("a tv array must be specified for each cdp");
            if (ntnmo > 0 && ntnmo != nnmo) err("a tnmo array must be specified for each cdp");
            if (nvnmo > 0 && nvnmo != nnmo) err("a vnmo array must be specified for each cdp");

            cdp = ealloc1float(ncdp);
            getparfloat("cdp",cdp);
            if (nkeys < 2) // just one key, check for sorting
                for (i=1; i<ncdp; ++i)
                    if (cdp[i]<=cdp[i-1])
                        err("key (%s) values must increase monotonically", key[0]);

        } else {  // just one t-v function, keys are irrelevant
            nkeys = 0;  // reset, because it is meaningless
            if (ntv>1)  err("only one tv array must be specified");
            if (nvnmo>1) err("only one vnmo array must be specified");
            if (ntnmo>1) err("only one tnmo array must be specified");
        }
        if ( ntv == 0 && ntnmo == 0 && nvnmo == 0 ) {
            err("one tv or tnmo/vnmo array must be specified");
        }

        keys = ealloc2float(nkeys + 1, maxRows);
        v = ealloc2float(1, maxRows);  // currently only one value, 2D array due to compatibility with sushw routings

        for (icdp=0; icdp<nnmo; ++icdp) {
            if (ntv > 0) {  // t-v function via tv arrays
                ntvnmo = countnparval(icdp+1,"tv");
                if ( ntvnmo%2 ) // odd number of values
                    err("number of values of t-v pairs must be even for %d-th tv array", icdp+1);
                tvnmo = ealloc1float(ntvnmo);
                getnparfloat(icdp+1,"tv",tvnmo);
                for (i=2; i<ntvnmo; i += 2)
                    if (tvnmo[i]<=tvnmo[i-2])
                        err("tnmo values must increase monotonically for %d-th tv array", icdp+1);
                for (i=0; i<ntvnmo/2; ++i) {
                    if ( nRow > maxRows -1 )
                        err("vfile %s exceeds limit of %d lines! Increase maxRows", vfile, maxRows);
                    v[nRow][0] = tvnmo[2*i + 1];
                    keys[nRow][nkeys] = tvnmo[2*i];
                    if ( nkeys > 0 )
                        for (j=0; j<nkeys; ++j)
                            keys[nRow][j] = cdp[nkeys*icdp + j];
                    ++nRow;
                }
                free1float(tvnmo);
            } else {
                nvnmo = countnparval(icdp+1,"vnmo");
                ntnmo = countnparval(icdp+1,"tnmo");
                if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
                    err("number of vnmo and tnmo values must be equal");
                if (nvnmo==0) nvnmo = 1;
                if (ntnmo==0) ntnmo = nvnmo;
                /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
                vnmo = ealloc1float(nvnmo);
                tnmo = ealloc1float(nvnmo);
                if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 1500.0;
                if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
                for (it=1; it<ntnmo; ++it)
                    if (tnmo[it]<=tnmo[it-1])
                        err("tnmo values must increase monotonically");
                for (i=0; i<ntnmo; ++i) {
                    v[nRow][0] = vnmo[i];
                    keys[nRow][nkeys] = tnmo[i];
                    if ( nkeys > 0 )
                        for (j=0; j<nkeys; ++j)
                            keys[nRow][j] = cdp[nkeys*icdp + j];
                    ++nRow;
                }
                free1float(vnmo);
                free1float(tnmo);
            }
        }
    }

    if (verbose) warn("Header key words used=%d; totally %d lines read", nkeys, nRow);
    if ( verbose>30 ) { // print out info for visual check
        for (i=0; i<nRow; ++i) {
            for (ikey=0; ikey<nkeys; ++ikey){
                fprintf(stderr, "%5s=%7d ", key[ikey], (int) keys[i][ikey]);
            }
            fprintf(stderr, " %6.3f ", keys[i][nkeys]);
            fprintf(stderr, "%4.0f ", v[i][0]);
            fprintf(stderr, "\n");
        }
    }

    /* get other optional parameters */
    if (!getparfloat("scalof", &scalof)) {
        scalof = ( tr.scalco < 0 )? -1.0/tr.scalco : (tr.scalco > 0)? tr.scalco : 1.0;
        if ( verbose > 0) warn(" Offset scaling factor=%.4f", scalof);
    }
    if (countparval("scale") > 0) { //
        int nf = getparfloat("scale", fscale);
        if (nf < nkeys) for(i=nf; i<nkeys; ++i) fscale[i] = fscale[nf -1];
    } else {
        for(i=0; i<nkeys; ++i)  fscale[i] = 1.0;
    }
    /* get other optional parameters */
    if (!getparfloat("smute",&smute)) smute = 1.5;
    if (smute<=0.0) err("smute must be greater than 0.0");
    if (!getparint("lmute",&lmute)) lmute = 25;
    if (!getparint("sscale",&sscale)) sscale = 1;
    if (!getparint("invert",&invert)) invert = 0;
    if (!getparint("upward",&upward)) upward = 0;

    int sort[3] = {1, 1, 1};  // sorting is always ascending
    /* allocate workspace */
    vnmo = ealloc1float(nt);
    ttn = ealloc1float(nt);
    atn = ealloc1float(nt);
    qtn = ealloc1float(nt);
    tnt = ealloc1float(nt);
    at = ealloc1float(nt);
    qt = ealloc1float(nt);

    /* loop over traces */
    int ntr = 0, ntotal = 0;
    do {
        ++ ntotal;
        /* loop over matching key fields and get values */
        cwp_Bool isEqual = cwp_true;
        /* get header values */
        for (ikey=0; ikey<nkeys; ++ikey) {
            gethval(&tr, index[ikey], &val[ikey]);
            fkeys[ikey] = fscale[ikey]*vtof(type[ikey], val[ikey]);
            isEqual = isEqual && !valcmp(type[ikey], val[ikey], valast[ikey]);
            valast[ikey] = val[ikey];
        }

        // get time-velocity function
        if (!ntr || !isEqual) { // not first or the same key values as last trace
            for (it=0, t=ft; it<nt; ++it, t += dt) {
                fkeys[nkeys] = t;
                interp(nkeys + 1, fkeys, keys, sort, nRow, v, 1, &vnmo[it]);
            }
            if ( verbose>30 ) { // print out info for visual check
                for (ikey=0; ikey<nkeys; ++ikey){
                    fprintf(stderr, "%5s=%d ", key[ikey], (int) fkeys[ikey]);
                }
                for (it=0; it<nt ; ++it){
                    if (!((it)%25)) fprintf(stderr, "\nt=%4.2f ", ft+it*dt);
                    fprintf(stderr, "%4.0f ", vnmo[it]);
                }
                fprintf(stderr, "\n");
            }
        }

        ++ntr;
        float x = scalof*tr.offset;
        /* compute time t(tn) (normalized) */
        float temp = x*x/(dt*dt);
        for (it=0,tn=ft/dt; it<nt; ++it,tn+=1.0) {
            tsq = temp/(vnmo[it]*vnmo[it]);
            ttn[it] = sqrt (tn*tn + tsq);
        }
        /* compute inverse of stretch factor a(tn) */
        atn[0] = ttn[1]-ttn[0];
        for (it=1; it<nt; ++it)
            atn[it] = ttn[it]-ttn[it-1];

        /* determine index of first sample to survive mute */
        osmute = 1.0/smute;
        if( !upward ) {
            for (it=0; it<nt-1 && atn[it]<osmute; ++it)
                ;
        } else {
            /* scan samples from bottom to top */
            for (it=nt-1; it>0 && atn[it]>=osmute; --it)
                ;
        }
        itmute = it;

        /* if inverse NMO will be performed */
        if (invert) {

            /* compute tn(t) from t(tn) */
            yxtoxy(nt-itmute,1.0,ft/dt+itmute,&ttn[itmute],
                    nt-itmute,1.0,ft/dt+itmute,
                    ft/dt-nt,ft/dt+nt,&tnt[itmute]);

            /* adjust mute time */
            itmute = 1.0+ttn[itmute]-ft/dt;
            itmute = MIN(nt-2,itmute);

            /* compute a(t) */
            if (sscale) {
                for (it=itmute+1; it<nt; ++it)
                    at[it] = tnt[it]-tnt[it-1];
                at[itmute] = at[itmute+1];
            }
        }
		
        /* if forward (not inverse) nmo */
        if (!invert) {
	
            /* do nmo via 8-point sinc interpolation */
            ints8r(nt,1.0,ft/dt,tr.data,0.0,0.0,
                   nt-itmute,&ttn[itmute],&qtn[itmute]);
			
            /* apply mute */
            for (it=0; it<itmute; ++it)
                qtn[it] = 0.0;
			
            /* apply linear ramp */
            for (it=itmute; it<itmute+lmute && it<nt; ++it)
                qtn[it] *= (float)(it-itmute+1)/(float)lmute;
			
            /* if specified, scale by the NMO stretch factor */
            if (sscale)
                for (it=itmute; it<nt; ++it)
                    qtn[it] *= atn[it];
			
            /* copy NMO corrected trace to output trace */
            memcpy( (void *) tr.data,
                     (const void *) qtn, nt*sizeof(float));
		
            /* else inverse nmo */
        } else {
	
            /* do inverse nmo via 8-point sinc interpolation */
            ints8r(nt,1.0,ft/dt,tr.data,0.0,0.0,
                   nt-itmute,&tnt[itmute],&qt[itmute]);
			
            /* apply mute */
            for (it=0; it<itmute; ++it)
                qt[it] = 0.0;
			
            /* if specified, undo NMO stretch factor scaling */
            if (sscale)
                for (it=itmute; it<nt; ++it)
                    qt[it] *= at[it];
			
            /* copy inverse NMO corrected trace to output trace */
            memcpy( (void *) tr.data,
                     (const void *) qt,nt*sizeof(float));
        }

        puttr(&tr); //  after writing output trace */
    } while (gettr(&tr)); // get next trace

    if (verbose > 0) warn("  %d of total %d traces are processed", ntr, ntotal);

    free1float(vnmo);
    free1float(ttn);
    free1float(atn);
    free1float(qtn);
    free1float(tnt);
    free1float(at);
    free1float(qt);
    if(keys) free2float(keys);
    free2float(v);

    return(CWP_Exit());
}

float xtscale(float x,      // offset
            float t,      // time
            float depth,  // depth of receiver in meter
            float* vrms,  // velocity
            float ft,     // time of first sample
            float dt,     // sampling interval
            float fmax    // max. scaling factor
            )
{
    const float maxa = 0.866;
    const float sinmax = sqrt(1.0 - 1.0/(fmax*fmax));
    int it = NINT((t - ft) / dt);
    float a = (t < dt)? maxa : x/(t * vrms[it]);
    if ( a > maxa ) a = maxa;
    float t0 = t*sqrt(1 - a*a);  // time at zero offset (apex)
    if (t0 < ft) t0 = ft;
    it = NINT((t0 - ft) / dt);
    float sintheta = x / (t*vrms[it]*vrms[it]/vrms[0] - 3.0*depth*(vrms[it]*vrms[it]/(vrms[0]*vrms[0]) - 1)/sqrt(1 - a*a));
    if ( sintheta < 0 || sintheta > sinmax ) sintheta = sinmax;
    return (float) (1.0/sqrt(1 - sintheta*sintheta));
}

/* Copyright (c) READD Well Service, 2007.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
                        /******************************************************************************
* Functions with prototypes for functions used internally to interpolate
* header values based on up to 3 keywrods
* Written by Sanyu Ye, sanyu.ye@readgroup.com
******************************************************************************
*
******************************************************************************/
/**************** end self doc ********************************/

int FindInterval(
                int ikey,        /* key level/index to search (which match key) */
                float dv,       /* value to be searched */
                float *dprev,   /* values of previous keys that must be matched during the search */
                float **keys,  /* array containing the input table of match keys */
                int *sort,      /* sorting order of the input match keys, >0 ascending, <0 descending */
                int nRows,       /* Total number of the rows/lines of input match keys */
                int *nbegp,      /* input: starting search line; output: new starting line for next search */
                int* n1p, int* n2p, /* interval of line number between which the searched value located */
                float* d1p, float* d2p /* interval of key values between which the search value located */
                        /* return value =-2 beyond first row; =+-1 between two rows; =0 exact match; =2 beyond last row */
                )
{
    int startIntervalFound = 0, match = sort[ikey];
    float prevKeyValue = (float) INT_MIN;
    for (int n = *nbegp; n < nRows; ++n) {  // loop over rows of match keys
        int isMatchPrevKeys = 1;
    // in case of search for 2nd or third key, find first match in previous keys
        for (int i=0; i<ikey && isMatchPrevKeys; ++i) {
            isMatchPrevKeys = isMatchPrevKeys && dprev[i] == keys[n][i];
        }

        if (!isMatchPrevKeys ) continue; /* skip to first line that previous key values matched*/

        if ( keys[n][ikey] == prevKeyValue ) continue;  // skip line with same key value
        else { // key value change
            prevKeyValue = keys[n][ikey];

            if ( dv == keys[n][ikey] ) { // exact match found, break
                                    *nbegp = n;  /* new start point for next level search */
                                    *d1p = *d2p = keys[n][ikey];
                                    *n1p = *n2p = n;
                                    match = 0;
                                    break;
            } else if ( sort[ikey]*dv > sort[ikey]*keys[n][ikey] ) { // first value of interval found
                                    *nbegp = n;  /* new start point for next level search */
                                    *d1p = *d2p = keys[n][ikey];
                                    *n1p = *n2p = n;
                                    startIntervalFound = 1;
            } else if (sort[ikey]*dv < sort[ikey]*keys[n][ikey]) /*check if second value of interval found */ {
                *d2p = keys[n][ikey];
                *n2p = n;

                if (!startIntervalFound) { // even the start is not found
                                        *nbegp = n;  /* new start point for next level search */
                                        *d1p = *d2p;
                                        *n1p = *n2p;
                                        match = 2;
                }
                break;
            }
        }
    }
    return match;
}

void interp1(
            int     ikey,   /* index of match key */
            int     n1,     /* first row */
            int     n2,     /* last row between them the values are interpolated */
            float   d,      /* key value of trace to which the interpolation is done */
            float  **keys, /* 2-D array of input table values of match keys */
            float  **var, /* 2-D array of input table values */
            int     nv,     /* array size */
            float*  v   /* output values interpolated */
            )
{
    int i;
    for ( i=0; i<nv; ++i) {
        if (d <= keys[n1][ikey]) v[i] = var[n1][i];
        else if (d >= keys[n2][ikey]) v[i] = var[n2][i];
        else if( n1 == n2 || keys[n2][ikey] == keys[n1][ikey] ) v[i] = (var[n1][i] + var[n2][i])/2;
        else v[i] = var[n1][i] + (var[n2][i] - var[n1][i])*(d - keys[n1][ikey])/(keys[n2][ikey] - keys[n1][ikey]);
    }
    return;
}

void interp2(
            float  d,      /* key value of trace to which the interpolation is done */
            float  d1,      /* first key value  */
            float  d2,      /* second key value  */
            float* v1,  /* first 1-D array to interpolation between */
            float* v2,  /* second 1-D array to interpolation between */
            int    nv,    /* array size */
            float* v    /* output values interpolated */
            )
{
    int i;
    for ( i=0; i<nv; ++i) {
        if (d <= d1) v[i] = v1[i];
        else if (d >= d2) v[i] = v2[i];
        else if( d1 == d2 ) v[i] = (v1[i] + v2[i])/2;
        else v[i] = v1[i] + (v2[i] - v1[i])*(d - d1)/(d2 - d1);
    }
    return;
}

int interp(
        int nkeys,      /* Number of the match keys */
        float *d,      /* 1D array of key values to be searched */
        float **keys,  /* 2D array containing the input table of match keys */
        int *sort,      /* sorting order of the input match keys, >0 ascending, <0 descending */
        int nRows,      /* Total number of the rows/lines of input match keys */
        float **var,   /* 2D array containing the input table of values */
        int nv,         /* Total number of the keys to be set */
        float* v       /* interval of key values between which the search value located */
                /* return value =-2 beyond first row; =+-1 between two rows; =0 exact match; =2 beyond last row */
        )
{
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14;
    int nbeg1=0, nbeg3, nbeg5, nbeg7, nbeg9, nbeg11, nbeg13;
    float d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14;
    float *v3=NULL, *v5=NULL, *v7=NULL, *v9=NULL, *v11=NULL, *v13=NULL;
    float dprev[2] = {0, 0};
    int match;

    if ( nkeys > 3 ) err("interpolation not implemented for more than 3D");

    match = FindInterval(0, d[0], dprev, keys, sort, nRows, &nbeg1, &n1, &n2, &d1, &d2);
    if (verbose > 100)
        fprintf(stderr, "level=%d d0=%f dprev=%f k1=%f k2=%f nbeg=%d n1=%d n2=%d d1=%f d2=%f\n",
                0, d[0], dprev[0], keys[n1][0], keys[n2][0], nbeg1, n1, n2, d1, d2);
    if ( 1 == nkeys ) { /* just 1-D interpolation */
        interp1(0, n1, n2, d[0], keys, var, nv, v );
        return 1;
    } else { /* one level more */
        dprev[0] = d1; nbeg3 = n1;
        match = FindInterval(1, d[1], dprev, keys, sort, nRows, &nbeg3, &n3, &n4, &d3, &d4);
        if (verbose>100)
            fprintf(stderr, "level=%d d0=%f dprev=%f k1=%f k2=%f nbeg=%d n1=%d n2=%d d1=%f d2=%f\n",
                    1, d[1], dprev[0], keys[n3][1], keys[n4][1], nbeg3, n3, n4, d3, d4);
        if ( d1 != d2 ) {
            dprev[0] = d2; nbeg5 = n2;
            match = FindInterval(1, d[1], dprev, keys, sort, nRows, &nbeg5, &n5, &n6, &d5, &d6);
            if (verbose>100)
                fprintf(stderr, "level=%d d0=%f dprev=%f k1=%f k2=%f nbeg=%d n1=%d n2=%d d1=%f d2=%f\n",
                        1, d[1], dprev[0], keys[n5][1], keys[n6][1], nbeg5, n5, n6, d5, d6);
        } else {
            n5 = n3; n6 = n4; d5 = d3; d6 = d4;
        }
        if ( 2 == nkeys ) { /* 2D */
            v3 = ealloc1float(nv);
            v5 = ealloc1float(nv);
            interp1(1, n3, n4, d[1], keys, var, nv, v3 );
            interp1(1, n5, n6, d[1], keys, var, nv, v5 );
            interp2(d[0], d1, d2, v3, v5, nv, v );
            free1float(v3);
            free1float(v5);
            return 2;
        } else { /* 3D, one level more */
            dprev[0] = d1; dprev[1] = d3; nbeg7 = n3;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg7, &n7, &n8, &d7, &d8);
            dprev[0] = d1; dprev[1] = d4; nbeg9 = n4;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg7, &n9, &n10, &d9, &d10);
            dprev[0] = d2; dprev[1] = d5; nbeg11 = n5;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg11, &n11, &n12, &d11, &d12);
            dprev[0] = d2; dprev[1] = d6; nbeg13 = n6;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg13, &n13, &n14, &d13, &d14);
            v3 = ealloc1float(nv);
            v5 = ealloc1float(nv);
            v7 = ealloc1float(nv);
            v9 = ealloc1float(nv);
            v11 = ealloc1float(nv);
            v13 = ealloc1float(nv);
            interp1(2, n7, n8, d[2], keys, var, nv, v7 );
            interp1(2, n9, n10, d[2], keys, var, nv, v9 );
            interp1(2, n11, n12, d[2], keys, var, nv, v11 );
            interp1(2, n13, n14, d[2], keys, var, nv, v13 );
            interp2(d[1], d3, d4, v7, v9, nv, v3 );
            interp2(d[1], d5, d6, v11, v13, nv, v5 );
            interp2(d[0], d1, d2, v3, v5, nv, v );
            free1float(v3);
            free1float(v5);
            free1float(v7);
            free1float(v9);
            free1float(v11);
            free1float(v13);
            return 3;
        }
    }
}