/* Copyright (c) READ Well Services, Oslo, 2009.*/
/* All rights reserved.                       */

/* S: $Revision: 1.00 $ ; $Date: 2009/05/31 22:25:06 $		*/

#include "su.h"
#include "segy.h"

#ifndef VALUETYPE
#define VALUETYPE float
#define EALLOC1VALUETYPE(x) ealloc1float(x)
#define FREE1VALUETYPE(x) free1float(x)
#endif

/*********************** self documentation ******************************/
char *sdoc[] = {
"									     ",
" SPXTSCALE - x-t variantly scale Z component for better multiple suppression",
"									     ",
"   spxtscale <stdin >stdout [optional parameters]			     ",
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
"   tnmo=           times [s] corresponding to velocities in vnmo	     ",
"   vnmo=           rms velocities [m] corresponding to times in tnmo	     ",
"									     ",
" More optional Parameters:						     ",
" 									     ",
"   maxlines=512    max. number of lines in ascii vfile                      ",
"   scale=1.0,...   scaling factors for header values specified by key=k3,k2 ",
" 									     ",
"   power=1.0       power of scaling factor                                  ",
"   max=2.0         max. scaling factor applied                              ",
" 									     ",
"   scalof=scalco   scaling factor for offset, default tr.scalco             ",
"                                                                            ",
"   trid=trid       trace id key to its certain value the scaling is applied ",
"   apply=2         value of trace id key to which the scaling is applied    ",
"                   =0 apply to all traces                                   ",
"   depth=gwdep     header key word that denotes water depth at receiver     ",
"   scalde=scalco   scaling factor for depth, default tr.scalel              ",
"                                                                            ",
" Notes:								     ",
"									     ",
" Scaling factor is calculated by incident angle of first water bottom multiple",
"									     ",
"  f = 1/cos(i) = 1/sqrt(1-sin(i)^2) = 1+sin(i)^2/2+3*sin{i}^4/8+5*sin(i)^8/16  ",
"",
"  sin(i) ~= X/(T*V) < 1 and T > sqrt(X*X + D*D)/V1 time of direct water arrival",
"									     ",
"  here X=offset, T=time, V=rms velocity, V1=1480 water velocity             ",
"  D=tr.gwdep/tr.scalel water depth at receiver ",
"",
"  f^power is the total scaling factor applied, which is limited to max      ",
"",
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
"  Since float is used to represent all key values in vfile, they should not",
"  exceed a precision of 6-7 digits.",
"",
" Examples: ",
"",
" spxtscale key=ep,cdp vfile=ep-cdp-t-v.table maxlines=2000 <in.su >out.su   ",
"",
" spxtscale key=gy,gx scale=0.1,0.1 vfile=vcube.rsf trid=ep apply=6       \\ ",
"           n=1001,20,10 d=0.008,500,500 o=0,3000,2000 <in.su >out.su   ",
"",
" spxtscale tv=0.400,1500,1.100,1834,1.900,2140,2.567,2393,3.138,2630 \\ .  ",
"           depth=gelev scalde=-0.1 scalof=1.0 <in.su >out.su   ",
"",
" spxtscale key=cdp cdp=1,240 <in.su >out.su  \\",
"        tnmo=0.400,1.100,1.900,2.567,3.138  vnmo=1500,1834,2140,2393,2630 \\",
"        tnmo=0.600,1.300,2.100,2.867,3.438  vnmo=1500,1734,2040,2493,2730   ",
"            ",
"									     ",
NULL};

/* Credits:
            *	RWS: Sanyu Ye, 2009, sanyu.ye@readgroup.com
            *
            * Trace header fields accessed: ns, dt, delrt, offset ...
*/
            /**************** end self doc *******************************************/


float xtscale(float x, float t, float depth, float* vrms, float ft, float dt, float fmax);

segy tr;

int verbose = 0;

#include "interp_incl.c"

int main(int argc, char **argv)
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
    float fmax;         /* max. scaling factor for amplitude */
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
    float power;        /* power of scaling factor */
    float scalde;       /* depth scaling factor */
    cwp_String depth, dtype;   /* keyword for depth */
    int dindex;
    cwp_String trid, itype;   /* keyword for trace id */
    int iindex, apply;

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
    if ( verbose>10 ) { // print out info for visual check
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
    if (!getparfloat("max",&fmax))              fmax = 2;
    if (!getparfloat("power",&power))           power = 1.0;
    if (countparval("scale") > 0) { //
        int nf = getparfloat("scale", fscale);
        if (nf < nkeys) for(i=nf; i<nkeys; ++i) fscale[i] = fscale[nf -1];
    } else {
        for(i=0; i<nkeys; ++i)  fscale[i] = 1.0;
    }
    if (!getparfloat("scalof", &scalof)) {
        scalof = ( tr.scalco < 0 )? -1.0/tr.scalco : (tr.scalco > 0)? tr.scalco : 1.0;
        if ( verbose > 0) warn(" Offset scaling factor=%.4f", scalof);
    }
    if (!getparstring("trid", &trid))     trid = "trid";
    itype  = hdtype(trid);
    iindex = getindex(trid);
    if (!getparint("apply", &apply))      apply = 12;
    if (apply && verbose > 0)
        warn(" X-T scaling is applied only to traces with header value tr.%s=%d", trid, apply);

    if (!getparstring("depth", &depth))     depth = "gwdep";
    dtype  = hdtype(depth);
    dindex = getindex(depth);
    if (!getparfloat("scalde", &scalde)) {
        scalde = ( tr.scalel < 0 )? -1.0/tr.scalel : (tr.scalel > 0)? tr.scalel : 1.0;
        if ( verbose > 0) warn(" Depth scaling factor=%.4f", scalof);
    }

    int sort[3] = {1, 1, 1};  // sorting is always ascending
    vnmo = ealloc1float(nt);

    /* loop over traces */
    int ntr = 0, ntotal = 0;
    do {
        ++ ntotal;
        // check trid id
        gethval(&tr, iindex, &val[0]);
        if ( apply && (apply != vtoi(itype, val[0])) ) {
            puttr(&tr);
            continue;
        }
        /* loop over matching key fields and get values */
        cwp_Bool isEqual = cwp_true;
        /* get header values */
        for (ikey=0; ikey<nkeys; ++ikey) {
            gethval(&tr, index[ikey], &val[ikey]);
            fkeys[ikey] = fscale[ikey]*vtof(type[ikey], val[ikey]);
            isEqual = isEqual && !valcmp(type[ikey], val[ikey], valast[ikey]);
            valast[ikey] = val[ikey];
        }

        // get depth
        gethval(&tr, dindex, &val[0]);
        float fdepth = scalde*vtof(dtype, val[0]);
        
        // get time-velocity function
        if (!ntr || !isEqual) { // not first or the same key values as last trace
            for (it=0, t=ft; it<nt; ++it, t += dt) {
                fkeys[nkeys] = t;
                interp(nkeys + 1, fkeys, keys, sort, nRow, v, 1, &vnmo[it]);
            }
            if ( verbose>10 ) { // print out info for visual check
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
        cwp_Bool toPrint = cwp_true;
        float flast=1;
        for (it=nt-1, t=ft+(nt-1)*dt; it>=0; --it, t -= dt) { //loop from end of trace to start
            float x = scalof*tr.offset;
            float tw = sqrtf(x*x + fdepth*fdepth)/vnmo[0]; // direct water wave time
            float f = xtscale(x, t, fdepth, vnmo, ft, dt, fmax);
            if (verbose==33 && !(it%10)) { // print every scaling for visual check
                fprintf(stderr, "x=%4.0f t=%5.3f f=%5.3f\n", x, t, f);
            }
            // Normally f should increase monotonically with decreasing time
            // if f starts to decrease or t less than direct water wave arrival, stop, use last value
            if ( (f - flast) < -0.01 || t < (tw - 0.02) ) {
                f = flast;
                if ( verbose>10 && toPrint) { // print out info for visual check
                    for (ikey=0; ikey<nkeys; ++ikey){
                        fprintf(stderr, "%5s=%d ", key[ikey], (int) fkeys[ikey]);
                    }
                    fprintf(stderr, "x=%4.0f t=%5.3f fstop=%5.3f\n", x, t, flast);
                    toPrint = cwp_false;
                }
            }
            float s = powf(f, power);
            if ( s > fmax ) s = fmax;
            tr.data[it] *= s;
            flast = f;
        }
        puttr(&tr); //  after writing output trace */
    } while (gettr(&tr)); // get next trace

    free1float(vnmo);
    free2float(keys);
    free2float(v);

    if (verbose > 0) warn("  %d of total %d traces are x-t scaled", ntr, ntotal);

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
    const float maxa = 1.0;
    int it = NINT((t - ft) / dt);
    float a = (t < dt)? maxa : fabs(x)/(t * vrms[it]);
    if ( a > maxa ) a = maxa;
    a = a*a;
    return 1.0 + 0.5*a + 3.0*a*a/8.0 + 5.0*a*a*a/16 + 35.0*a*a*a*a/64;
}
