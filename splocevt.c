/* Copyright (c) Sanyu Ye <sanyu_ye@yahoo.com>, 2012.   */
/* All rights reserved.                                 */

/* SPLOCEVT written by Sanyu Ye <sanyu_ye@yahoo.com>     */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPLOCEVT - LOCate EVenT by least-square fitting of arrival times and estimating",
"            receiver to source polarization vector/azimuth for microseismic data",
"                                                                               ",
"   splocevt <stdin >stdout [options]                                            ",
"                                                                               ",
" Optional Parameters:                                                          ",
"                                                                               ",
"   key=ep          gather sort key                                             ",
"   nmax=12         max. number of receivers expected in a gather               ",
"   nc=4            number of components expected in a gather                   ",
"                                                                               ",
"   level=3         =0: project to P/SV/SH towards source/injection point       ",
"                   =1: estimate average velocity by mimimizing traveltime error",
"                   =2: estimate source lateral location by minimizing tt error ",
"                   =3: estimate azimuth/location by computing polarization vectors",
"                       in conjunction of results of level=2                    ",
"   tolt=4    [ms]  tolerance of tt standard deviation of level=2 location      ",
"                   beyond which the location attempt is regarded as failed     ",
"   tola=30  [deg]  tolerance of azimuth standard deviation of level=3 location ",
"                   beyond which the location attempt is regarded as failed     ",
"   bias=0.0        factor to favour P or S wave type in case the equal success ",
"                   level. should be -1 < bias < 1. if < 0 favour for P; > 0 for S",
"   threst=1.6      exclusion threshold a tt pick is discarded from location    ",
"                   if its tt/azimuth error exceeds threshold x standard deviation",
"   thresa=1.0      exclusion threshold an azimuth is discarded from location   ",
"                   if its azimuth error exceeds threshold x standard deviation ",
"   perc=75         minimum percentage of receivers needed for location         ",
"                                                                               ",
"   mute1=20  [ms]  time offset to predicted P arrivals for top mute            ",
"   mute2=40  [ms]  time offset to predicted S arrivals for bottom mute         ",
"                                                                               ",
"   ps=0            separate P from S waves by projecting P on receiver to source",
"                   vector while S on transverse vertical and horizontal (only  ",
"                   for level=3)                                                ",
"                                                                               ",
"   vp=      [m/s]  P velocity, if not given, fetch from header keyword  wevel  ",
"   vs=      [m/s]  S velocity, if not given, fetch from header keyword swevel  ",
"                                                                               ",
"   ixyz=x,y,z [m]  center coordinates of injection points                      ",
"                   defauft fetched from header keywords sx, sy, sdepth         ",
"                                                                               ",
"   dr=5       [m]  search step in radial direction                             ",
"   dz=2       [m]  z search step                                               ",
"                                                                               ",
"   rwin=1500  [m]  radius of search window in radial direction from monitor well",
"   zwin=300   [m]  half length of search window in z centered at source/injection point",
"                                                                               ",
"                                                                               ",
"   awin=90         search (half) windows in degree for azimuth                 ",
"                   centered around vector of receiver to injection point       ",
"   da=3,1   [deg]  increment of rotation angle for first and second run        ",
"                                                                               ",
"   tref=laga       keyword where S arrival time [1/10 ms] is saved or fetched  ",
"   twin=0,20       time windows in ms around S arrival for compute polarization",
"                                                                               ",
"   gdepth=gelev    keyword holding geophone depth, scaled by scalco as usual   ",
"                                                                               ",
"   verbose=0       >0 output info                                              ",
"                                                                               ",
" Note:                                                                         ",
"   Microseismic data are input in the strict order of ZNE, receiver by receiver,",
"   event by event.                                                             ",
"   First it is attempted to determine the depth and lateral distance of source ",
"   by least square fitting of arrival times picked automatically during level 3",
"   SPEVENT. In an iterative process the pickes with largest errors are excluded",
"   and remaining arrival time are refitted. In a second stage data around the  ",
"   predicted arrival time, defined by twin=beg,end are used to determine the   ",
"   polarization vector (incidence angle in azimuth/inclination). Median azimuth",
"   of all receivers is calculated, again those with largest errors are excluded",
"   in an iterative process.                                                    ",
"   Three successive performance levels of location process are defined:        ",
"     SUCCESS=1 if traveltime fitting yields a standard deviation less than the ",
"       given tolerance (tErr < tolt).                                          ",
"     SUCCESS=2 if in addition the standard deviation of azimuth is less than   ",
"       the given tolerance (azErr < tola).                                     ",
"     SUCCESS=3 if in addition the standard deviation of inclination is less than",
"       the given tolerance (inErr < tola).                                     ",
"   The above procedure is performed twice, first assuming P arrivals, then S.  ",
"   A higher success level will set the corresponding wave type and its results ",
"   valid as final. In case of both wave types achieving the same success level ",
"   (usually =1 for traveltime fitting), the product of time and angle errors   ",
"   are used to determine the wave type. the criteria is expressed as:          ",
"        tErrP*(azErrP + inErrP)  vs  (1 - bias) * tErrS*(azErrS + inErrS)      ",
"   if left hand side, i.e. of P is smaller than right hand side, i.e. of S, P wave",
"   is chosen, otherwise S. Usually a positive bias factor is used to favour S  ",
"   over P because in noisy recording, S waves with higher energy are more apparent",
"   than P and thus dominate recording.                                         ",
"   The resulting success level is stored in header corr, along with the wave type",
"   identified in counit (0 for S, 1 for P). The fitted/predicted traveltimes are",
"   in header lagb. Times for top muting above P and bottom muting below S are  ",
"   stored in muts and mute respectively. Source position and its errors are in ",
"   fx/fy/fz and dx/dy/dz. In addition the median azimuth and its standard deviation",
"   are stored in unscale/ungpow. These parameters are used for next stage location",
"   based on 3C migration SPMS3CMIG.",
"                                                                               ",
" Caveats:                                                                      ",
" Lateral distance estimated by S traveltime fitting ususally smaller while by P",
" larger compared with result based on migration location. If the source depth is",
" out of receiver array, i.e. no apex in traveltime curve, the error would become",
" larger, particularly when fluctuation in time picking of arrivals is bigger.  ",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
" 1: QC the location result by checking the mute after location                 ",
"                                                                               ",
" mirf2su file=$RecNo.rcd length=8.1 verbose=10 overlap=0.1 |\\",
" sushw key=styp a=$StageNo |\\",
" subfilt fstoplo=10 fpasslo=20 fpasshi=240 fstophi=360 |\\",
" sushw match=styp key=sx,sy,sdepth infile=injpoint.xyz |\\",
" sushw match=duse key=unscale values=1,-1,2,1,3,-1 |\\",
" suhtmath op=mult key=unscale |\\",
" sushw match=cdp key=otrav infile=rcv-orient.tbl scale=10 |\\",
" sp2crot key=cdp mode=0 base=0 vsp=-1 a=otrav scale=0.1 |\\",
" spgsort key=fldr sort=ep |\\",
" spevent key=ep level=2 eeout=1 verbose=1 wine=200 thres=3.0 |\\",
" spevent key=ep level=3 nc=4 wina=10 verbose=1 |\\",
" suwind key=corr min=1 |\\",
" sushw key=wevel,swevel a=4200,2400 |\\",
" splocevt nc=4 level=3  verbose=1 ps=1 tola=35 twin=-2,18 mute1=30 mute2=50 |\\",
" tee $StageNo-$RecNo-l3.su |\\",
" sugain pbal=1 |\\",
" sumute mode=0 tkey=muts |\\",
" sumute mode=1 tkey=mute |\\",
" spsort ep duse |\\",
" suxwigb perc=99.9 grid1=solid windowtitle=\"Mute according to predicted P & S arrivals\" & ",
"                                                                               ",
"                                                                               ",
" Version 1.1.0 last modified  Dec. 2012 by Sanyu Ye                            ",
NULL
};

/* Credits:
 *	Sanyu Ye  <sanyu_ye@yahoo.com>, first created May. 2012
 */
/**************** end self doc *******************************************/

// forward declaration
float calcVel(const int nrcv, float* xrcv, float* yrcv, float* zrcv,
         float* ixyz, float* t, float* t0, float* pv);
float seekLateralLoc(const int nrcv, int* notuse, const float v, float* xrcv, float* yrcv, float* zrcv,
         float* ixyz, const float radius, const float zwin, float* dxyz, float* Terr, float* t, float* t0, float* x0, float* z0);
float LateralLoc(const int nrcv, int* notuse, const float v, float* xrcv, float* yrcv, float* zrcv,
        float* ixyz, const float radius, const float zwin, float* dxyz, float* tref, float* Terr, float* t0, float* xs, float* zs, int* nused,
        const int verbose, int* iorder, const float tolt, const float perc, const float thres, float* T2sort);
float calcExtPol(const int mode, const int opt, const int nrcv, const int nt, const float dt, const float* t0,
        const float* twin, const float* awin, const float* da,
        const float azim0, float* incl0, float*** data, float* azim, float* incl);
float calcExtP(int mode, int itmin, int itmax, float azim0, float wazm, float dazm,
        float** zyx, float* azim, float* incl);
float seekExtP(int mode, int itmin, int itmax, float azim0, float incl0, float wazm, float dazm,
        float** zyx, float* azim, float* incl);
int EstAzInc(const int verbose, const int nrcv, int* notuse, int* iorder, const float perc, const float thresa,
        float tola, float* T2sort, float* azim, float* incl, float* am, float* aErr, float* im, float* iErr);
void Project3C2P(const int itmin, const int itmax, const float azim, const float incl, float** zyx, float* p);
float calcMedianAngle(int n, float* angles, int* notuse);
float stddev(float* a, int n, float avg, int* notuse);
float calcSrcLoc(const int nrcv, const int* notuse, const float* xr, const float* yr, const float* zr,
        const float* azim, const float* incl, const float azimuth, float* x, float* y, float* z, float* d);
float seekSrcLoc(const int nrcv, float* xr, float* yr, float* zr, float* azim, float* incl,
        const float xc, const float yc, const float zc, const float* cube, const float* dxyz, float* px, float* py, float* pz);
float calcDistP2L(const float x, const float y, const float z, const float x1,
        const float y1, const float z1, const float azim, const float incl);
void  calcDirection(const float x1, const float y1, const float z1, const float x2,
                    const float y2, const float z2, float* azim, float* incl);
float calcRMS(int nt, float* data);


#include "sprinthlp.c"

int main(int argc, char **argv)
{
    cwp_String key, tkey, key_g;    /* header key word from segy.h		*/
    cwp_String type, ttype, type_g; /* type of key				*/
    Value val, valnew, vala;        /* value of key                             */
    int index, tindex, index_g;     /* index of key				*/
    int nt, nc, nmax, perc;
    int nsegy;                  /* length bytes read for segy trace	*/
    int deviated, ps, opt;
    int verbose;		/* =0 no info, >0 output info  */
    int level;
    float v[2], threst, thresa, tolt, tola, bias, mute1, mute2;
    float azErr[2], inErr[2], tErr[2], dErr, zErr;
    float dt;                
    float xs, ys, zs[2], rs[2], xc, yc, zc, t0[2];
    float da[2], awin[2], twin[2];
    float azimuth[2], azimc, inclc, azims, incls;
    float azimrs, inclrs;
    float ixyz[3] = {0, 0, 0}, radius, zwin, dxyz[3];
    int i, j, itr, ntr, nrcv, ngather, total;
    segy tr, outtr;

    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparint("nmax", &nmax))       nmax = 12;
    if (!getparint("nc", &nc))             nc = 4;
    if (!getparint("level", &level))    level = 3;
    if (!getparint("ps", &ps)) ps = 0;
    if (!getparint("opt", &opt)) opt = 1;
    if (!getparint("perc", &perc)) perc = 75;

    if(!getparfloat("threst", &threst))   threst = 1.6;
    if(!getparfloat("thresa", &thresa))   thresa = 1.0;
    if(!getparfloat("tolt",   &tolt))     tolt  = 4.0;
    if(!getparfloat("tola",   &tola))     tola  = 30.0;
    if(!getparfloat("mute1", &mute1))     mute1 = 20.0;
    if(!getparfloat("mute2", &mute2))     mute2 = 40.0;
    if(!getparfloat("bias",  &bias))       bias = 0.0;
    if (bias <= -1.0 || bias >= 1.0) err("Invalid parameter bias (=%3.1f) must be (-1 < bias < 1)", bias);
    
    /* get gather sorting key */
    if (!getparstring("key", &key))	 key="ep";
    type = hdtype(key);
    index = getindex(key);
    /* get key holding S-arrival time in 1/10 ms*/
    if (!getparstring("tref", &tkey))	 tkey="laga";
    ttype = hdtype(tkey);
    tindex = getindex(tkey);
    /* get key holding geophone depth (z true depth)*/
    if (!getparstring("gdepth", &key_g)) key_g="gelev";
    type_g = hdtype(key_g);
    index_g = getindex(key_g);

    int n = countparval("da");
    if ( n > 0 && n <= 2 ) {
        getparfloat("da", da);
        for(i=n; i<2; ++i) da[i] = da[n-1]/3.0;
    } else {
        da[0] = 3.0; da[1] = 1.0;
    }

    n = countparval("awin");
    if ( n > 0 && n <= 2 ) {
        getparfloat("awin", awin);
        for(i=n; i<2; ++i) awin[i] = awin[n-1];
    } else {
        awin[0] = 90.0; awin[1] = 90.0;
    }

    n = countparval("twin");
    if ( n == 2 ) {
        getparfloat("twin", twin);
    } else if (n == 0) {
        twin[0] = 0.0; twin[1] = 20;
    } else {
        err("Time window (twin=begin,end) must be specified");
    }

    if(!getparfloat("dr",   &dxyz[0]))     dxyz[0] = 5.0;
    if(!getparfloat("dy",   &dxyz[1]))     dxyz[1] = dxyz[0];
    if(!getparfloat("dz",   &dxyz[2]))     dxyz[2] = 2.0;

    if(!getparfloat("rwin", &radius))  radius = 1500;
    if(!getparfloat("zwin", &zwin))      zwin =  300;

    /* get first trace */
    if ( (nsegy=gettr(&tr)) < HDRBYTES) err("can't get first trace");
    gethval(&tr, index, &val);

    if(!getparfloat("vp", &v[1])) v[1] = (float) tr.wevel;
    if (level == 2 && v[1] < 2000.0) err("Probably too small P velocity (=%3.0f)", v[1]);
    if(!getparfloat("vs", &v[0])) v[0] = (float) tr.swevel;
    if (level == 2 && v[0] < 1000.0) err("Probably too small S velocity (=%3.0f)", v[0]);
    if ( verbose > 1) warn("velocity Vp=%4.0f  Vs=%4.0f", v[1], v[0]);

    dt = ((double) tr.dt) / 1000000.0;
    nt = tr.ns;
    float scalco = ( tr.scalco < 0 )? -1.0/tr.scalco : (tr.scalco > 0)? tr.scalco : 1.0;
    if ( verbose > 1) warn(" Offset scaling factor=%.4f", scalco);

    n = countparval("ixyz");
    if ( n == 3 ) {
        getparfloat("ixyz", ixyz);
    } else {
        if (verbose > 9) warn("Source/Injection point corrdinates fetched from keywords sx, sy, sdepth");
    }

    segyhdr** hdrs = (segyhdr **) ealloc2(nc, nmax, HDRBYTES);
    float***  indata = ealloc3float(nt, nc, nmax);
    float***  psdata = ealloc3float(nt,  3, nmax);
    float** azim = ealloc2float(nmax, 2);
    float** incl = ealloc2float(nmax, 2);
    float* dist = ealloc1float(nmax);
    float* xrcv = ealloc1float(nmax);
    float* yrcv = ealloc1float(nmax);
    float* zrcv = ealloc1float(nmax);
    float* tref = ealloc1float(nmax);
    float* t    = ealloc1float(nmax);  // calculated arrivals for P or S
    float* buf = ealloc1float(nt);
    float**  Terr = ealloc2float(nmax, 2);
    float**    tt = ealloc2float(nmax, 2);  // calculated arrivals for P and S for muting purpose
    float** incll = ealloc2float(nmax, 2);  // calculated inclination from every receiver to lateral source position
    int* notuse = ealloc1int(nmax);
    int* exclud = ealloc1int(nmax);
    int* iorder = ealloc1int(nmax);

    // reset with zero
    memset(*hdrs, 0, nc * nmax * HDRBYTES);
    memset(**indata, 0, nc * nmax * nt * FSIZE);
    memset(**psdata, 0,  3 * nmax * nt * FSIZE);
    memset(*azim, 0, 2 * nmax * FSIZE);
    memset(*incl, 0, 2 * nmax * FSIZE);
    memset(dist, 0, nmax * FSIZE);
    memset(xrcv, 0, nmax * FSIZE);
    memset(yrcv, 0, nmax * FSIZE);
    memset(zrcv, 0, nmax * FSIZE);
    memset(tref, 0, nmax * FSIZE);
    memset(*tt, 0, 2*nmax * FSIZE);
    memset(buf, 0, nt * FSIZE);
    memset(notuse, 0, nmax * ISIZE);

    xc = yc = zc = 0.0; // init center coodinate of receiver array
    int eof = 0;
    ngather = ntr = total = 0;
    do {  // Read headers and data while getting a count
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nc*nmax - 1) err("\nNumber of traces exceeding nmax=%d\n", nmax);
            int ircv = ntr/nc;
            int icmp = ntr%nc;
            memcpy(&hdrs[ircv][icmp], &tr, HDRBYTES);
            memcpy(indata[ircv][icmp], tr.data, nt*FSIZE);
            if (!ntr && ixyz[2] == 0) {  // first trace
                // get injection point coordinates
                ixyz[0] = scalco*tr.sy;
                ixyz[1] = scalco*tr.sx;
                ixyz[2] = scalco*tr.sdepth;
            }
            if (!icmp) { // first component of every receiver
                yrcv[ircv] = scalco*tr.gx;  // internal coordinate x-N, y-E, z-Z downwards
                xrcv[ircv] = scalco*tr.gy;
                gethval(&tr, index_g, &vala);
                zrcv[ircv] = scalco*vtof(type_g, vala);
                gethval(&tr, tindex, &vala);
                tref[ircv] = 0.0001*vtof(ttype, vala);
                if (tref[ircv] < 10.0*dt) notuse[ircv] = 1;
            }
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            nrcv = ntr/nc;
            total += ntr;
            ++ngather;
            
            int n = 0;
            for (i=0; i<nrcv; ++i) n += notuse[i];
            if ( n > nrcv/4 ) {
                warn("Opps, should never happen: %d receivers of %d-th gather (%s=%d) have no time pick!", n, ngather, key, vtoi(type, val));
                continue;
            }

            if (verbose) warn(" Processing %d traces %d-th gather (%s=%d) ...", ntr, ngather, key, vtoi(type, val));
            if (ntr%nc) err(" number of traces not multiple of number of components for %d-th gather (%s=%d)",
                    ntr, ngather, key, vtoi(type, val));

            // do coordinate transform if needed

            azimc = 180/PI*atan2(ixyz[1], ixyz[0]);
            if ( azimc < 0.0 ) azimc += 360.0;
            if (verbose > 1 && ngather == 1) {
                fprintf(stderr, "Source/Injection point coordinates: E/X/y=%3.1f  N/Y/x=%3.1f   Z=%5.1f   Azimuth=%3.1f [deg]\n",
                        ixyz[1], ixyz[0], ixyz[2], azimc);
                fprintf(stderr, "Receiver coordinates E/N/Z:");
                for (i=0; i<nrcv; ++i) fprintf(stderr, " %3.1f/%3.1f/%5.1f", yrcv[i], xrcv[i], zrcv[i]);
                fprintf(stderr, "\n");
            }
            if (verbose > 1) {
                fprintf(stderr, "\nTime picks:");
                for (i=0; i<nrcv; ++i) fprintf(stderr, " %d/%5.3f", i + 1, tref[i]);
                fprintf(stderr, "\n");
            }
         
            int PorS = 0;  // default presumed S arrivals
            int successful[2] = {0, 0};  // default presumed not valid/successful location
            int nrcvused[2];
            if (level == 1) {
                PorS = 1;  // always using P
                tErr[1] = calcVel(nrcv, xrcv, yrcv, zrcv, ixyz, tref, &t0[1], &v[1]);
                if (verbose) warn("Best estimate: V=%4.0f [m/s]  T0=%5.3f  tErr=%5.3f", v[1], t0[1], tErr[1]);
            } else if (level >= 2) {
                for (j=1; j>=0; --j) {  // loop over P and S 
                    if (verbose > 1) warn(" Assuming %s arrivals with velocity %4.0f and trying to locate ...", j? "P" : "S", v[j]);
                    memcpy(exclud, notuse, nmax*ISIZE);
                    tErr[j] = LateralLoc(nrcv, exclud, v[j], xrcv, yrcv, zrcv, ixyz, radius, zwin, dxyz, tref, Terr[j], &t0[j], &rs[j], &zs[j], &nrcvused[j],
                            verbose, iorder, tolt, perc, threst, buf);

                    if ( 1000.0*tErr[j] < tolt ) {
                        successful[j] = 1;
                        if (verbose>1) warn("Successful location depth/dist/tErr %4.0f/%1.0f/%3.1f for %s arrivals",
                            zs[j], rs[j], 1000.0*tErr[j], j? "P" : "S");
                    } else {
                        if (verbose>1) warn(" tt error too big (%3.1f > %3.1f [ms]) to be successful for lateral location", tErr[j]*1000.0, tolt);
                    }
                    
                    if (level > 2 && successful[j]) {
                        // calc center coordinates of receiver array
                        if (!xc && !yc && !zc) {  // if not set
                            for (i=0; i<nrcv; ++i) {
                                xc += xrcv[i];
                                yc += yrcv[i];
                                zc += zrcv[i];
                            }
                            xc /= nrcv;
                            yc /= nrcv;
                            zc /= nrcv;

                            calcDirection(xc, yc, zc, ixyz[0], ixyz[1], ixyz[2], &azimc, &inclc);
                            if ( azimc < 0.0 )   azimc += 360.0;
                            if ( azimc > 360.0 ) azimc -= 360.0;
                        }
                        
                        for (i=0; i<nrcv; ++i) { // calc arrival times
                            float d = sqrtf(rs[j]*rs[j] + (zs[j] - zrcv[i])*(zs[j] - zrcv[i]));
                            t[i] = (d/v[j] - t0[j]);
                            // calc inclination of lateral source position to receivers
                            incll[j][i] = 180.0/PI*atan2(zs[j] - zrcv[i], rs[j]);
                        }
                        
                        if (verbose > 9) {
                            fprintf(stderr, " Best fitting arrival times (picked/calculated):\n");
                            for (i=0; i<nrcv; ++i) fprintf(stderr, " %d/%5.4f/%5.4f", i+1, tref[i], t[i]);
                            fprintf(stderr, "\n Expected inclination [deg] to source position estimated by traveltime fitting:\n");
                            for (i=0; i<nrcv; ++i) fprintf(stderr, " %d/%1.0f", i+1, incll[j][i]);
                            fprintf(stderr, "\n");
                        }
                        
                        // do calculate polarization angle
                        calcExtPol(j, opt, nrcv, nt, dt, t, twin, awin, da, azimc, incll[j], indata, azim[j], incl[j]);
                        //seekPol(nrcv, nt, dt, tref, twin, awin, da, azimc, inclc, indata, azim, incl);

                        
                        int nexcl = EstAzInc(verbose, nrcv, notuse, iorder, perc, thresa, tola, buf, incll[j], azim[j], incl[j], &azimuth[j], &azErr[j], &inErr[j]);

                        azimuth[j] += azimc;  // back
                        if ( azimuth[j] <   0.0 ) azimuth[j] += 360.0;
                        if ( azimuth[j] > 360.0 ) azimuth[j] -= 360.0;

                        if (verbose > 1) {
                            warn("Estimated Azimuth=%3.1f with standard deviation (azim/incl)=%3.1f/%3.1f by %d receivers", 
                                    azimuth[j], azErr[j], inErr[j], nrcv - nexcl);
                        }
                        
                        if (azErr[j] < tola && inErr[j] < tola) {
                            successful[j] = 3;
                        } else if (azErr[j] < tola ) {
                            successful[j] = 2;
                        //} else {
                        //    azimuth[j] = azimc;
                        }
                        if (successful[j] > 0) {
                            if (verbose > 1) {
                                warn("%s successful estimate: Azimuth=%3.1f with az/inErr=%3.1f/%3.1f based on %d receivers",
                                       (successful[j] > 2)? "Very" : (successful[j] > 1)? "Partially" : "Little",
                                        azimuth[j], azErr[j], inErr[j], nrcv - nexcl);
                            }
                        } /*else {
                            if (verbose > 1) {
                                warn(" Error too big for both azimuth/inclination(%3.1f/%3.1f > %3.1f [deg]) to be successful",
                                         azErr[j], inErr[j], tola);
                                warn(" Using injection point for azimuth=%3.1f", azimc);
                            }
                        }*/
                    }
                }

                if ( successful[1] > successful[0] && tErr[1]) {
                    PorS = 1;
                } else if ( successful[1] == successful[0]) {
                    if ( rs[1] > 1.1*dxyz[0] // lateral distance by P not the starting value of search window
                        && tErr[1]*(azErr[1] + inErr[1])/nrcvused[1] < (1.0 - bias) * tErr[0]*(azErr[0] + inErr[0])/nrcvused[0])  PorS = 1;
                }
                zErr = tErr[PorS]*v[PorS];

                for (i=0; i<nrcv; ++i) { // calc P & S arrival times
                    float d = sqrtf(rs[PorS]*rs[PorS] + (zs[PorS] - zrcv[i])*(zs[PorS] - zrcv[i]));
                    tt[0][i] = (d/v[0] - t0[PorS]);
                    tt[1][i] = (d/v[1] - t0[PorS]);
                }

                dErr = rs[PorS] * azErr[PorS] * PI / 180;
                ys   = rs[PorS] * sinf(azimuth[PorS]*PI/180.0);
                xs   = rs[PorS] * cosf(azimuth[PorS]*PI/180.0);
                if (successful[PorS] > 0) warn("Source position by %s/%d: Z=%1.0f/%3.1f  R=%1.0f  Azimuth=%3.1f  X/E=%1.0f  Y/N=%1.0f\n              tErr=%3.1f [ms] Err=%3.1f/%3.1f [deg]  AErr=%3.1f [m]\n",
                        PorS? "P" : "S", successful[PorS], zs[PorS], zErr, rs[PorS], azimuth[PorS], ys, xs, tErr[PorS]*1000.0, azErr[PorS], inErr[PorS], dErr);

                // calc source location
                //dErr = calcSrcLoc(nrcv, notuse, xrcv, yrcv, zrcv, azim[0], incl[0], azimuth[0], &xs, &ys, &zs, dist);
                // calculate center vector to source just located
                calcDirection(xc, yc, zc, xs, ys, zs[PorS], &azims, &incls);
            }

            if (level == 0 || (level > 2 && ps && successful[PorS] > 0)) {
                for (itr = 0; itr < nrcv; ++itr) {
                    if ( level == 0 ) {
                        xs = ixyz[0];
                        ys = ixyz[1];
                        zs[1] = ixyz[2];
                    }
                    // calc polarization vector to source from individual receiver
                    calcDirection(xrcv[itr], yrcv[itr], zrcv[itr], xs, ys, zs[PorS], &azimrs, &inclrs);
                    // project/decompose input data to receiver-source polarization vector (P wave)
                    // and transverse horizontal/verical (S-waves)
                    Project3C2P(0, nt - 1, azimrs - 90.0,    0.0, indata[itr], psdata[itr][2]);  // SH
                    Project3C2P(0, nt - 1, azimrs, inclrs - 90.0, indata[itr], psdata[itr][1]);  // SV
                    Project3C2P(0, nt - 1, azimrs, inclrs,        indata[itr], psdata[itr][0]);  // P
                }
            }
            for (itr = 0; itr < nrcv; ++itr) {
                for (i = 0; i < nc; ++i) {
                    memcpy(&outtr, &hdrs[itr][i], HDRBYTES);
                    memcpy(outtr.data, indata[itr][i], FSIZE * nt);
                    if(level == 1) {
                        outtr.wevel = NINT(v[1]);
                        outtr.lagb = NINT(10000.0*t0[1]);
                    }
                    if (level >= 2) {
                        outtr.dz = zErr;
                        outtr.dy = zErr;
                        outtr.fy = rs[PorS];
                        outtr.fx = 0.0;
                        outtr.fz = zs[PorS];
                        outtr.lagb = NINT(10000.0*tt[PorS][itr]);
                        outtr.corr = successful[PorS];
                        outtr.counit = PorS;
                        outtr.muts = MAX(0, NINT((tt[1][itr]*1000.0 - mute1)));
                        outtr.mute = NINT(MIN((nt - 1)*dt*1000.0, tt[0][itr]*1000.0 + mute2));
                    }
                    if (level > 2) {
                        outtr.lcf  = (short) (10.0*azimrs);
                        outtr.hcf  = (short) (10.0*inclrs);
                        outtr.lcs  = (short) (10.0*azim[PorS][itr]);
                        outtr.hcs  = (short) (10.0*incl[PorS][itr]);
                        outtr.fx = ys;
                        outtr.dx = dErr;
                        outtr.fy = xs;
                        outtr.unscale = azimuth[PorS];
                        outtr.ungpow = azErr[PorS];
                        outtr.duse = i + 1;
                        outtr.lcf  = (short) (10.0*azimrs);
                        outtr.hcf  = (short) (10.0*inclrs);
                        //if (successful[PorS] > 2) {  // reset source/injection point to new source location
                        //    outtr.sx = NINT(ys/scalco);
                        //    outtr.sy = NINT(xs/scalco);
                        //    outtr.sdepth = NINT(zs[PorS]/scalco);
                        //}
                    }
                    puttr(&outtr);
                }
                if (level == 0 || (level == 3 && ps && successful[PorS] > 0)) {
                    for (i = 0; i < 3; ++i) {
                        memcpy(outtr.data, psdata[itr][i], FSIZE * nt);
                        outtr.duse = i + nc + 1;
                        puttr(&outtr);
                    }
                }
            }
            // reset output data
            memset(**indata, 0, nc * nmax * nt * FSIZE);
            memset(*hdrs, 0, nc * nmax * HDRBYTES);
            memset(**psdata, 0,  3 * nmax * nt * FSIZE);
            memset(*azim, 0, 2 * nmax * FSIZE);
            memset(*incl, 0, 2 * nmax * FSIZE);
            memset(dist, 0, nmax * FSIZE);
            memset(tref, 0, nmax * FSIZE);
            memset(*Terr, 0, 2*nmax * FSIZE);
            memset(buf, 0, nt * FSIZE);
            memset(notuse, 0, nmax * ISIZE);
            memset(*tt, 0, 2 * nmax * FSIZE);
            
            val = valnew;
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces for each of %d compomnents of %d gathers are processed", total/nc, nc, ngather);
    
    return(CWP_Exit());
}

float calcVel(const int nrcv, float* xrcv, float* yrcv, float* zrcv,
         float* ixyz, float* t, float* t0, float* pv)
{
    // least square solving object function sum(t[i] + t0 - d[i]/v)^2 regarding t0 and v

    int i;

    float d2, sumd2 = 0.0, sumd = 0.0, sumt = 0.0, sumtd = 0.0;
    for(i=0; i<nrcv; ++i) {
        d2 = (ixyz[0] - xrcv[i]) * (ixyz[0] - xrcv[i])
           + (ixyz[1] - yrcv[i]) * (ixyz[1] - yrcv[i])
           + (ixyz[2] - zrcv[i]) * (ixyz[2] - zrcv[1]);
        sumd2 += d2;
        sumd  += sqrtf(d2);
        sumt  += t[i];
        sumtd += t[i]*sqrtf(d2);
    }
    *pv  = (sumd*sumd - nrcv*sumd2)/(sumd*sumt - nrcv*sumtd);
    *t0 = (sumd/(*pv) - sumt)/nrcv;

    float tErr = 0;
    for(i=0; i<nrcv; ++i) {
        d2 = (ixyz[0] - xrcv[i]) * (ixyz[0] - xrcv[i])
           + (ixyz[1] - yrcv[i]) * (ixyz[1] - yrcv[i])
           + (ixyz[2] - zrcv[i]) * (ixyz[2] - zrcv[1]);
        tErr += (t[i] - *t0 - sqrtf(d2)/(*pv))*(t[i] - *t0 - sqrtf(d2)/(*pv));
    }
    return sqrtf(tErr/nrcv);  // std deviation of tt
}

float seekLateralLoc(const int nrcv, int* notuse, const float v, float* xrcv, float* yrcv, float* zrcv,
         float* ixyz, const float radius, const float zwin, float* dxyz, float* t, float* Terr, float* toffset, float* x0, float* z0)
{
    // search for minimum sum{(t + t0 - d[i]/v)^2}
    // t0 = (sum{d[i]}/v - sum{t[i]})/nrcv
    // d[i] = sqrt(x^2 + (z - zrcv[i])^2)
    // regarding two variable x and z
    int i;
    float x, z, minErr = FLT_MAX;

    // horizontal distance
    float z1 = ixyz[2] - zwin;
    float z2 = ixyz[2] + zwin;
    int n = nrcv;
    for (i=0; i<nrcv; ++i) n -= notuse[i];
    for (x = dxyz[0]; x <= radius; x += dxyz[0]) {
        for (z = z1; z <= z2; z += dxyz[2]) {
            float sumt = 0.0, sumd = 0.0;
            for(i=0; i<nrcv; ++i) {
                if (notuse[i]) continue;
                sumt += t[i];
                float d = sqrtf(x * x + (z - zrcv[i]) * (z - zrcv[i]));
                sumd += d;
            }
            float t0 = (sumd/v - sumt)/n;
            float tErr = 0.0;
            for(i=0; i<nrcv; ++i) {
                if (notuse[i]) continue;
                float d = sqrtf(x * x + (z - zrcv[i]) * (z - zrcv[i]));
                tErr += (t[i] + t0 - d/v)*(t[i] + t0 - d/v);
            }
            if (tErr < minErr) {
                minErr = tErr;
                *toffset = t0;
                *x0 = x;
                *z0 = z;
            }
        }
    }
    // calc tErr of individual receivers
    for(i=0; i<nrcv; ++i) {
        //if (notuse[i]) continue;
        float d = sqrtf((*x0) * (*x0) + (*z0 - zrcv[i]) * (*z0 - zrcv[i]));
        Terr[i] = t[i] + (*toffset) - d/v;
    }
    return sqrtf(minErr/nrcv);
}

float calcExtPol(const int mode, const int opt, const int nrcv, const int nt, const float dt, const float* t0,
        const float* twin, const float* awin, const float* da,
        const float azim0, float* incl0, float*** data, float* azim, float* incl)
{
    int itr, itmin, itmax;
    float rms = 0.0;

    for (itr=0; itr<nrcv; ++itr) {
        itmin = NINT((t0[itr] + 0.001*twin[0])/dt);
        itmin = MAX(0, itmin);
        itmax = NINT((t0[itr] + 0.001*twin[1])/dt);
        itmax = MIN(nt - 1, itmax);

        if (!opt) {
            // first run with big window and coarse grid
            rms = calcExtP(mode, itmin, itmax, azim0, awin[0], da[0], data[itr], &azim[itr], &incl[itr]);
            // second run, much small window and grid
            rms = calcExtP(mode, itmin, itmax, azim[itr], 2*da[0], da[1], data[itr], &azim[itr], &incl[itr]);
        } else {
            // first run with big window and coarse grid
            rms = seekExtP(mode, itmin, itmax, azim0, 0.0, awin[0], da[0], data[itr], &azim[itr], &incl[itr]);
            // second run, much small window and grid
            rms = seekExtP(mode, itmin, itmax, azim[itr], incl[itr], 2*da[0], da[1], data[itr], &azim[itr], &incl[itr]);
        }

        azim[itr] -= azim0;  // refer to center azimuth
        incl[itr] -= incl0[itr];
        
        // considering ambiguity of inclination like a or +-a
        //float inc = incl[itr];
        //if ( inc * incl0[itr] < 0.0 )  inc *= -1.0;  // not same sign
        //if (ABS(incl[itr] - incl0[itr]) > ABS(inc - incl0[itr])) incl[itr] = inc - incl0[itr];
        //else incl[itr] -= incl0[itr];
    }

    return rms;
}

float calcExtP(int mode, int itmin, int itmax, float azim0, float wazm, float dazm,
        float** zyx, float* azim, float* incl)
{
    int it;
    float azm;
    float x2 = 0.0, y2 = 0.0, z2 = 0.0, xy = 0.0, xz = 0.0, yz = 0.0;
    float prmsmin = FLT_MAX, prmsmax = 0.0;
    for (it=itmin; it<=itmax; ++it) {
        x2 += zyx[2][it]*zyx[2][it];
        y2 += zyx[1][it]*zyx[1][it];
        z2 += zyx[0][it]*zyx[0][it];
        xy += zyx[2][it]*zyx[1][it];
        xz += zyx[2][it]*zyx[0][it];
        yz += zyx[1][it]*zyx[0][it];
    }
    for (azm = azim0-wazm; azm <= azim0+wazm; azm += dazm) {
        float cosa = cosf(azm*PI/180.0);
        float sina = sinf(azm*PI/180.0);
        float inclm = atan((y2 - x2)*cosa*sina + xy*(1.0 - 2.0*sina*sina)/(xz*sina - yz*cosa));
        float cosi = cosf(inclm);
        float sini = sinf(inclm);
        float prms = x2*cosa*cosa*cosi*cosi + y2*sina*sina*cosi*cosi + z2*sini*sini
                   + 2.0*(xy*cosa*sina*cosi*cosi + xz*cosa*cosi*sini + yz*sina*cosi*sini);
        if (mode == 0 && prms < prmsmin) {  // S
            *azim = azm;
            *incl = inclm*180.0/PI;
            prmsmin = prms;
        } else if (mode > 0 && prms > prmsmax) {  // P
            *azim = azm;
            *incl = inclm*180.0/PI;
            prmsmax = prms;
        }
    }
    return sqrtf( (mode == 0)? prmsmin : prmsmax );
}

float seekExtP(int mode, int itmin, int itmax, float azim0, float incl0, float wazm, float dazm,
        float** zyx, float* azim, float* incl)
{
    int it;
    float azm, inc;
    float x2 = 0.0, y2 = 0.0, z2 = 0.0, xy = 0.0, xz = 0.0, yz = 0.0;
    float prmsmin = FLT_MAX, prmsmax = 0.0;
    for (it=itmin; it<=itmax; ++it) {
        x2 += zyx[2][it]*zyx[2][it];
        y2 += zyx[1][it]*zyx[1][it];
        z2 += zyx[0][it]*zyx[0][it];
        xy += zyx[2][it]*zyx[1][it];
        xz += zyx[2][it]*zyx[0][it];
        yz += zyx[1][it]*zyx[0][it];
    }
    for (azm = azim0-wazm; azm <= azim0+wazm; azm += dazm) {
        for (inc = incl0-wazm; inc <= incl0+wazm; inc += dazm) {
            float cosa = cosf(azm*PI/180.0);
            float sina = sinf(azm*PI/180.0);
            float cosi = cosf(inc*PI/180.0);
            float sini = sinf(inc*PI/180.0);
            float prms = x2*cosa*cosa*cosi*cosi + y2*sina*sina*cosi*cosi + z2*sini*sini
                       + 2.0*(xy*cosa*sina*cosi*cosi + xz*cosa*cosi*sini + yz*sina*cosi*sini);
            if (mode == 0 && prms < prmsmin) {  // S
                *azim = azm;
                *incl = inc;
                prmsmin = prms;
            } else if (mode > 0 && prms > prmsmax) {  // P
                *azim = azm;
                *incl = inc;
                prmsmax = prms;
            }
        }
    }
    return sqrtf( (mode == 0)? prmsmin : prmsmax );
}

void Project3C2P(const int itmin, const int itmax, const float azim, const float incl, float** zyx, float* p)
{
    int it;

    float sina = sinf(azim*PI/180.0);
    float cosa = cosf(azim*PI/180.0);
    float sini = sinf(incl*PI/180.0);  // inclination positive downward from horizon (-90/up ~ 90/down)
    float cosi = cosf(incl*PI/180.0);
    for (it=itmin; it<=itmax; ++it) {
        p[it - itmin] = cosa*cosi*zyx[2][it] + sina*cosi*zyx[1][it] + sini*zyx[0][it];
    }
    return;
}

float calcMedianAngle(int n, float* angles, int* notuse)
{
    int i, m = 0;
    float a[n];
    for (i=0; i<n; ++i) {
        if ( notuse && notuse[i] ) ;
        else a[m++] = angles[i];
    }
    float am = quick_select(a, m);
    //free1float(a);
    return am;
}

float stddev(float* a, int n, float avg, int* notuse)
{
    int i, m = 0;
    float diff;
    float asum = 0.0;
    for (i=0; i<n; ++i) {
        if ( notuse && notuse[i] ) ;
        else {
            diff = a[i] - avg;
            if (ABS(diff) > 180.0) diff = 360.0 - ABS(diff);
            asum += diff*diff;
            ++m;
        }
    }
    return (m>1)? sqrtf(asum/m) : 0.0;
}

float calcSrcLoc(const int nrcv, const int* notuse, const float* xr, const float* yr, const float* zr,
   const float* azim, const float* incl, const float azimuth, float* x, float* y, float* z, float* d)
{
    // compute source location that has minimum of square distances to all rays through receiver
    
    int i, n=0;
    float d2sum = 0.0;;
    
    float ** A = ealloc2float(3, 3);
    float *  B = ealloc1float(3);
    float *  D = ealloc1float(3);
    int*   idx = ealloc1int(3);
    memset(*A, 0, 9*FSIZE);
    memset( B, 0, 3*FSIZE);
    
    for (i=0; i<nrcv; ++i) {
        if (notuse[i]) continue;
        //float a1 = cosf(azim[i]*PI/180.0)*cosf(incl[i]*PI/180.0);
        //float a2 = sinf(azim[i]*PI/180.0)*cosf(incl[i]*PI/180.0);
        float a1 = cosf(azimuth*PI/180.0)*cosf(incl[i]*PI/180.0);
        float a2 = sinf(azimuth*PI/180.0)*cosf(incl[i]*PI/180.0);
        float a3 = sinf(incl[i]*PI/180.0);
        A[0][0] +=  a2*a2 + a3*a3;
        A[0][1] += -a1*a2;
        A[0][2] += -a1*a3;
        B[0]    +=  a1*a2*yr[i] - (a2*a2 + a3*a3)*xr[i] + a1*a3*zr[i];
        A[1][0] += -a1*a2;
        A[1][1] += -a1*a1 + a3*a3;
        A[1][2] += -a2*a3;
        B[1]    +=  a1*a2*xr[i] - (a1*a1 + a3*a3)*yr[i] + a2*a3*zr[i];
        A[2][0] += -a1*a3;
        A[2][1] += -a2*a3;
        A[2][2] +=  a1*a1 + a2*a2;
        B[2]    +=  a1*a3 + a2*a3 - (a1*a1 + a2*a2)*zr[i];
    }
    
    LU_decomposition(3, A, idx, D);
    backward_substitution(3, A, idx, B);
        
    *x = B[0]; 
    *y = B[1]; 
    *z = B[2];
    
    // loop over receivers
    for (i = 0; i < nrcv; ++i) {
        if (notuse[i]) continue;
        ++n;
        d[i]   = sqrtf(calcDistP2L(*x, *y, *z, xr[i], yr[i], zr[i], azim[i], incl[i]));
        d2sum += d[i]*d[i];
    }
    return sqrtf(d2sum/n);
}

float seekSrcLoc(const int nrcv, float* xr, float* yr, float* zr, float* azim, float* incl,
        const float xc, const float yc, const float zc, const float* cube, const float* dxyz, float* px, float* py, float* pz)
{
    int i;
    float x, y, z, d2sum, d2min;
    // loop over the whole cube
    for(x = xc - 0.5*cube[0]; x <= xc + 0.5*cube[0]; x += dxyz[0]) {
        for(y = yc - 0.5*cube[1]; y <= yc + 0.5*cube[1]; y += dxyz[1]) {
            for(z = zc - 0.5*cube[2]; z <= zc + 0.5*cube[2]; z += dxyz[2]) {
                d2sum = 0.0;
                d2min = FLT_MAX;
                for (i = 0; i < nrcv; ++i) {
                    d2sum += calcDistP2L(x, y, z, xr[i], yr[i], zr[i], azim[i], incl[i]);
                }
                if (d2sum < d2min) {
                    d2min = d2sum;
                    *px = x;
                    *py = y;
                    *pz = z;
                }
            }
        }
    }
    return sqrtf(d2min/nrcv);
}

float calcDistP2L(const float x, const float y, const float z, const float x1,
        const float y1, const float z1, const float azim, const float incl)
{
    // d= |(X2-X1)x(X1-X)|/|X2-X1| = |(a1, a2, a3)x(x1-x,y1-y,z1-z)|
    //  =|(a2(z1-z) - a3(y1-y), a3(x1-x) - a1(z1-z), a1(y1-y) - a2(x1-x)|

    float a1 = cosf(azim*PI/180.0)*cosf(incl*PI/180.0);
    float a2 = sinf(azim*PI/180.0)*cosf(incl*PI/180.0);
    float a3 = sinf(incl*PI/180.0);

    float dist2 = (a2*(z1 - z) - a3*(y1 - y)) * (a2*(z1 - z) - a3*(y1 - y))
                + (a3*(x1 - x) - a1*(z1 - z)) * (a3*(x1 - x) - a1*(z1 - z))
                + (a1*(y1 - y) - a2*(x1 - x)) * (a1*(y1 - y) - a2*(x1 - x));

    return dist2;
}

void  calcDirection(const float x1, const float y1, const float z1, const float x2,
                    const float y2, const float z2, float* azim, float* incl)
{
    float dist = sqrtf((y2 - y1)*(y2 - y1) + (x2 - x1)*(x2 - x1));
    *azim = 180/PI*atan2(y2 - y1, x2 - x1);
    *incl = 180/PI*atan2(z2 - z1, dist);
    return;
}

float calcRMS(int nt, float* data)
{
    float sum = 0.0;
    int it, nc = 0;
    for (it = 0; it < nt; ++it ) {
        if (data[it] != 0.0 ) {
            ++nc;
            sum += data[it]*data[it];
        }
    }
    return (nc == 0)? 0.0 : sqrt(sum/nc);
}

float LateralLoc(const int nrcv, int* notuse, const float v, float* xrcv, float* yrcv, float* zrcv,
        float* ixyz, const float radius, const float zwin, float* dxyz, float* tref, float* Terr, float* t0, float* xs, float* zs, int* nused,
        const int verbose, int* iorder, const float tolt, const float perc, const float thres, float* T2sort)
{
    int i;
    int low  = NINT((float) nrcv * perc / 100.0);
    int high = nrcv - 1;
    float tErr = seekLateralLoc(nrcv, notuse, v, xrcv, yrcv, zrcv, ixyz, radius, zwin, dxyz, tref, Terr, t0, xs, zs);
    if (verbose > 9) fprintf(stderr, "\n Lateral location depth/dist/tErr %4.0f/%1.0f/%3.1f  T0=%3.1f for V=%4.0f\n",
            *zs, *xs, 1000.0*tErr, 1000.0*(*t0), v);

    int nexcl = 0;
    int itera = 0;
    do {
        // check time errors to exclude outliers
        for (i=0; i<nrcv; ++i) {
            iorder[i] = i;
            T2sort[i] = ABS(Terr[i]);
        }
        qkisort(nrcv, T2sort, iorder);

        nexcl = 0;
        for(i=high; i >= low; --i) {
            if (T2sort[iorder[i]] > thres*tErr) {
                ++nexcl;
                notuse[iorder[i]] = 1;
                if (verbose > 9) fprintf(stderr, "  tErr[%d]=%3.1f", iorder[i] + 1, 1000.0*Terr[iorder[i]]);
            }
        }
        ++itera;

        if (verbose > 9 && nexcl) fprintf(stderr, " --> %d receivers are excluded\n", nexcl);

        tErr = seekLateralLoc(nrcv, notuse, v, xrcv, yrcv, zrcv, ixyz, radius, zwin, dxyz, tref, Terr, t0, xs, zs);

        if (verbose>9) {
            fprintf(stderr, " Iteration %d: depth=%4.0f  r=%3.1f [m]   T0=%3.1f  tErr=%3.1f [ms] for V=%4.0f based on %d receivers\n",
                itera, *zs, *xs, 1000.0*(*t0), 1000.0*tErr, v, nrcv - nexcl);
        }
    } while (1000.0*tErr > tolt &&  nrcv - low > itera);

    *nused = nrcv - nexcl;
    
    return tErr;
}

int EstAzInc(const int verbose, const int nrcv, int* notuse, int* iorder, const float perc, const float thresa,
        float tola, float* T2sort, float* incll, float* azim, float* incl, float* am, float* aErr, float* iErr)
{
    int i;

    float azimuth = calcMedianAngle(nrcv, azim, NULL);
    float azErr = stddev(azim, nrcv, azimuth, NULL);
    float inErr = stddev(incl, nrcv, 0.0, NULL);

    if (verbose > 9) {
        fprintf(stderr, "\n DiffAzimth/Incl to InjPts/Source azimuth=%1.0f/%1.0f  medianInclination=%1.0f/%1.0f\n", azimuth, azErr, incll[nrcv/2], inErr);
        for (i=0; i<nrcv; ++i) fprintf(stderr, " %d/%1.0f/%1.0f ", i+1, azim[i], incl[i] + incll[i]);
        fprintf(stderr, "\n");
    }


    int nexcl = 0;
    int itera = 0;
    int low  = NINT((float) nrcv * perc / 100.0);
    int high = nrcv - 1;
    do {
        for (i=0; i<nrcv; ++i) {
            iorder[i] = i;
            T2sort[i] = ABS(azim[i] - azimuth);
        }
        qkisort(nrcv, T2sort, iorder);

        ++itera;
        nexcl = 0;
        memset(notuse, 0, nrcv*ISIZE);
        for(i=high; i >= low; --i) {
            if (T2sort[iorder[i]] > thresa*azErr) {
                ++nexcl;
                notuse[iorder[i]] = 1;
                if (verbose > 9) fprintf(stderr, "  az/inErr[%d]=%3.1f/%3.1f", iorder[i] + 1, azim[iorder[i]] - azimuth, incl[iorder[i]]);
            } else {
                notuse[iorder[i]] = 0;
            }
        }
        if (verbose > 9 && nexcl) fprintf(stderr, " --> %d receivers are excluded\n", nexcl);
        azimuth = calcMedianAngle(nrcv, azim, notuse);
        azErr = stddev(azim, nrcv, azimuth, notuse);
        inErr = stddev(incl, nrcv, 0.0, notuse);
        if (verbose>9) {
            fprintf(stderr, " Iteration %d: Revised Angle=%3.0f with Error(azim/incl)=%1.0f/%1.0f based on %d receivers\n",
                itera, azimuth, azErr, inErr, nrcv - nexcl);
        }
    } while (azErr > tola &&  nrcv - low > itera);

    *am   = azimuth;
    *aErr = azErr;
    *iErr = inErr;

    return nexcl;

}
