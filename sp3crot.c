/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPCOMPROT v. 1.0 written by Sanyu Ye */


#include <math.h>
#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SP3CROT - Compute and apply 3C orientation/projection for 4C OBC data         ",
"                                                                               ",
"   spcomprot <stdin >stdout [options]                                          ",
"                                                                               ",
" Optional Parameters:                                                          ",
"                                                                               ",
"   nc=4          number of input components in order [H,] Z, X and Y trace by trace",
"                 =3  Vertical(Z), Inline(X) and Crossline(Y)                   ",
"                                                                               ",
"   pol=1,1,1     polarities of Z, X, Y                                         ",
"                                                                               ",
"   mode=0        =0 apply projection with given orientation of 3C              ",
"                 =1 compute orientation/rotation angle                         ",
"                 =2 compute orientation and apply rotation/projection          ",
"                                                                               ",
"   step=1        =1 compute orientation/rotation angles using Oloffson method  ",
"                 =2 refine by maximizing energy using individual rotation angles",
"                 =3 refine by maximizing energy using median rotation angles   ",
"                                                                               ",
"   key=cdpt      receiver gather sort key, used to compute angle               ",
"   nxmax=512     max. number of traces for a component expected in a gather    ",
"                                                                               ",
"   gdepth=gwdep  keyword holding geophone depth, scaled by scalco as usual     ",
"   obcazim=hour  keyword holding cable azimuth, in 100th degrees (-180 ~ 180)  ",
"                                                                               ",
" Parameters used to determine orientation angles of 3 geophone component:      ",
"                                                                               ",
"   t0=1.0  [s]   time line peaks of first arrivals are aligned on              ",
"                 i.e. first positive peak for OBC data                         ",
"   vw=1.48       P wave velocity of water                                      ",
"   vp=vw         P wave velocity of seafloor                                   ",
"                                                                               ",
"   as=3,5,5      angles (size) of sector to compute median value of rotations  ",
"                 azimuth, tilt and roll                                        ",
"                                                                               ",
"   amp=0         =0 output original amplitude of polarization vector           ",
"                 =1 output amplitude multiplied by polarities (pol=z,x,y)      ",
"                                                                               ",
"   t1=0.980 [s]  start of time window to search for local max. total 3C amplitude",
"                 should be set possibly close to onset of first arrivals       ",
"   t2=1.020 [s]  end of time window                                            ",
"                                                                               ",
" Parameters used for further refining steps by maximizing total energy projected",
" on source-receiver vector within the time window given above                  ",
"                                                                               ",
"   aw=5,10,15    half width of search window for azimuth, tilt and roll        ",
"   da=0.5  [deg] increment of rotation angle                                   ",
"                                                                               ",
"   verbose=0     >0 output info                                                ",
"                                                                               ",
" Note:                                                                         ",
"   For details of algorithm see Olofsson et al. 2007 SEG extended abstract     ",
"   Determining 3C geophone orientation from a single shot.                     ",
"                                                                               ",
"   Inline (X), Crossline (Y), and Vertical (Z) build a right hand Cartesian    ",
"   coordinate system. With Z axis positive downwards, and Inline (X) normally  ",
"   coinciding with surface X axis (to right/east), the Crossline Y is positive ",
"   down/south, opposite of surface Y axis, which is normally positive up/north.",
"   3C geophones built-in within OBC may differ from above requirement. If Tilt ",
"   angle and its stardand deviation show unreasonable values (normally < 10/5),",
"   most likely the polarity of inline component is wrong and should be reversed.",
"   If standard deviation of computed Roll angle is very large, either Z or Y or",
"   both are reversed in polarity. The parameter pol= can be used to set polarities",
"   of Z/X/Y correctly.                                                         ",
"   In case of multiple sources in a gather, distribution of determined rotation",
"   angles is computed. The sector with highest concentration and its two neighbouring",
"   sectors are used to determine the median value.                             ",
"   Input data must be shifted in the way that the FTBs are aligned at specified",
"   time (t0, 1 sec by defauft). The FTB vector is used to determined tilt and roll ",
"   of the cable stored in keyword minute and sec respectively. Azimuths of the ",
"   cable (Inline component) must be given in keyword specified by obcazim=     ",
"   (default =hour). The source-to-receiver azimuth and transmission angle at   ",
"   geophone are stored in year and day. Median azimuth, tilt and roll are hold ",
"   in keywords timbas, grnors, grnlof and its respective std deviation in keywords",
"   trwf, grnofr, gaps when a gather has multiple sources. The orientation (azimuth",
"   and inclination) of individual component are stored in fx, dx, fy, dy and fz",
"   and dz. Note the sense of inclination is positive downwards. Amplitudes of  ",
"   polarization vector are stored in keyword ungpow for later QC plot. In case ",
"   of refining steps (step > 1) the rms amplitudes are extracted.              ",
"   Internally right-hand coordinate system is used, i.e. X points to east, Y south",
"   and Z downwards; so remain the same the sense of particle motion after      ",
"   rotation is applied. Parameter pol= should be used to set the sense of FTB of",
"   3C geophones all positive for shot located NW of receiver.                  ",
"   Note t1= and t2= are used to search local max. total amplitude of 3C geophones.",
"   This is necessary for optical sensors because of phase shift between hydrophone",
"   and geophones. In other words FTB (e.g. peak) of hydrophone does not correspond",
"   peak of total 3C geophone amplitude (in strict sense zero crossing since the",
"   geophones record acceleration, the time derivative of velocity which is     ",
"   recorded by hydrophone.                                                     ",
"   When used for applying 3C rotation, more correctly, projection, the orientations",
"   of 3C X/Y/Z, given in azimuth and inclination (positive downward), should be",
"   loaded in keywords fx/dx, fy/dy and fz/dz in unit of degrees. The projected ",
"   resultss are output in the same component as input with the sense of polarity",
"   as the coordinate system used, namely positive X east, Y south, Z downward. ",
"                                                                               ",
" Major steps for determining 3C orientation:                                   ",
" 1. Calulate cable azimuth                                                     ",
" 2. Use hydrophone to determine FTB                                            ",
" 3. Load FTB to geophones and aligned them, load azimuth                       ",
" 4. Calculate 3C orientation/rotation                                          ",
" 5  Refine calculation by maximizing energy (rms amplitude) on polarization vector",
"                                                                               ",
" sp3crot nc=3 key=cdpt < input.su > output.su                        ",
"                                                                               ",
"                                                                               ",
" Version 1.0.4 Last Modified by Sanyu Ye on 07.03.11                           ",
NULL};

/* Credits:
 *	Sanyu Ye: READ Well Services, July. 2010
 */
/**************** end self doc *******************************************/

// forward declaration
void rotate3C(int itmin, int itmax, float azim, float tilt, float roll, const float* pol, float* x, float* y, float* z, float* x1, float* y1, float* z1);
void projectons2r(int nt, float azim, float incident, float* x, float* y, float* z, float* s2r);
void Project3C(const int ntr, const int itmin, const int itmax, const float* azim, const float* incl, const float* pol, float** x, float** y, float** z);
void Project1C(const int itmin, const int itmax, const float azim, const float incl, const float polarity, const float* x, float* x1, float* y1, float* z1);
void calc3COrient(const float gamma, const float theta, const float phi, const float* pol, float* azim3c, float* incl3c);
void calc3C(const float azims2r, const float incident, const float azimx, float* tilt, float* roll, const float x, const float y, const float z);
void calcOrientations(int nsrc, int itmin, int itmax, int* itftb, float* azims2r, float* incident,
        const float azim, float* tilt, float* roll, float** iline, float** xline, float** z, const float* pol, float** amp);
void refineOrientations(int nsrc, int itmin, int itmax, float* azims2r, float* incident, float* aw, const float da,
        float* azim, float* tilt, float* roll, float** iline, float** xline, float** z, const float* pol, float** amp);
double calcRMS(int nt, float* data);
float getMedianAngle(int ntr, float* ax, float alimit, float* stddiv, int* n);
float stddev(float* a, int n, float avg);
float quick_select(float *arr, int n);

int verbose;		/* =0 no info, >0 output info  */

int main(int argc, char **argv) {
    cwp_String key, key_a, key_g;	/* header key word from segy.h		*/
    cwp_String type, type_a, type_g;    /* type of key				*/
    Value val, valnew, vala;                  /* value of key			*/
    int index, index_a, index_g;         /* index of key				*/
    int nt, *itftb, itmin, itmax;       /* nums amps as int			*/
    int xaxis, yaxis, zaxis;
    int nc, nxmax;                     /* number of components			*/
    int nsegy;                  /* length bytes read for segy trace	*/
    int mode, step, ampi;
    int nam, ntm, nrm;
    float azim3c[3], incl3c[3], Azim, Tilt, Roll;   /* rotation angle computed or applied   */
    float stddvAzim, stddvTilt, stddvRoll;          /* standard deviation   */
    float aw[3], da, as[3], pol[3];                                /* search window widths and increment   */
    float vp, vw;                                   /* P velocity for seafloor and sea water   */
    float t0, t1, t2, dt, scale;                    /* input trace px py header values */
    int i, itr, isrc, ntr, ngather, total;
    segy tr, outtr;

    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose=0;

    if (!getparstring("obcazim", &key_a)) key_a="hour";
    type_a = hdtype(key_a);
    index_a = getindex(key_a);
    if (!getparstring("gdepth", &key_g)) key_g="gwdep";
    type_g = hdtype(key_g);
    index_g = getindex(key_g);
    if (!getparstring("key", &key))	 key="cdpt";
    type = hdtype(key);
    index = getindex(key);

    if (!getparfloat("vw", &vw)) vw = 1.48;
    if (!getparfloat("vp", &vp)) vp = vw;
    if (!getparfloat("scale", &scale)) scale = 0.01;
    if (!getparfloat("da", &da))  {da = 0.5;}
    //if (!getparint("na", &na))      na = 180;
    if (!getparint("nxmax", &nxmax)) nxmax = 512;
    int n = countparval("aw");
    if ( n > 0 && n <= 3 ) {
        getparfloat("aw", aw);
        for (i = n; i < 3; ++i) aw[i] = aw[n-1];
    }
    else {
        {aw[0] = 5.0; aw[1] = 10.0; aw[2] = 15.0;}
    }
    n = countparval("as");
    if ( n > 0 && n <= 3 ) {
        getparfloat("as", as);
        for (i = n; i < 3; ++i) as[i] = as[n-1];
    }
    else {
        {as[0] = 3.0; as[1] = 5.0; as[2] = 5.0;}
    }
    n = countparval("pol");
    if ( n > 0 && n <= 3 ) {
        getparfloat("pol", pol);
        for (i = n; i < 3; ++i) pol[i] = pol[n-1];
    }
    else {
        {pol[0] = 1.0; pol[1] = 1.0; pol[2] = 1.0;}
    }

    if (!getparint("nc", &nc))        nc = 4;
    if (!getparint("mode", &mode))  mode = 0;
    if (!getparint("step", &step))  step = 1;
    if (!getparint("amp",  &ampi))  ampi = 0;
    // find out component index corresponding to x and y axis
    xaxis = (nc == 4)? 2 : 1;
    yaxis = (nc == 4)? 3 : 2;
    zaxis = (nc == 4)? 1 : 0;

    /* get first trace */
    if ( (nsegy=gettr(&tr)) < HDRBYTES) err("can't get first trace");
    gethval(&tr, index, &val);

    dt = ((double) tr.dt) / 1000000.0;
    nt = tr.ns;

    if (!getparfloat("t0", &t0)) t0 = 1.0;
    if (!getparfloat("t1", &t1)) t1 = 0.980;
    if (!getparfloat("t2", &t2)) t2 = 1.020;
    //itftb = NINT(t0 / dt);
    itmin = NINT(t1 / dt);
    itmax = NINT(t2 / dt);
    if (mode > 0) { /* Check time gating values */
        if (itmin < 0)
            err("itmin=%d should be positive", itmin);
        //if (itftb < itmin)
        //    err("t0 %5.3f should be larger than t1 %5.3f", t0, t1);
        if (itmin > itmax)
            err("itmin=%d, itmax=%d conflict", itmin, itmax);
        if (tr.ns <= itmax)
            err("tr.ns=%d, itmax=%d window cannot extend over the trace length", tr.ns, itmax);
    }
    segyhdr** hdrs = (segyhdr **) ealloc2(nxmax, nc, HDRBYTES);
    float***  indata = ealloc3float(nt, nxmax, nc);
    float** amp = ealloc2float(3, nxmax);
    float* azim = ealloc1float(nxmax);
    float* tilt = ealloc1float(nxmax);
    float* roll = ealloc1float(nxmax);
    float* azims2r = ealloc1float(nxmax);
    float* incident = ealloc1float(nxmax);
    itftb = ealloc1int(nxmax);
    
    memset(*hdrs, 0, nc * nxmax * HDRBYTES);
    memset(**indata, 0, nc * nxmax * nt * FSIZE);
    memset(*amp, 0, 3*nxmax * FSIZE);
    memset(azim, 0, nxmax * FSIZE);
    memset(tilt, 0, nxmax * FSIZE);
    memset(roll, 0, nxmax * FSIZE);
    memset(azims2r, 0, nxmax * FSIZE);
    memset(incident, 0, nxmax * FSIZE);
    memset(itftb, 0, nxmax * sizeof(int));

    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = total = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nc*nxmax - 1) err("\nNumber of traces exceeding nxmax=%d\n", nxmax);
            memcpy(&hdrs[ntr%nc][ntr/nc], &tr, HDRBYTES);
            memcpy(indata[ntr%nc][ntr/nc], tr.data, FSIZE * nt);
            if (!(ntr%nc)) {
                isrc = ntr/nc;
                float x = tr.gx - tr.sx;
                float y = tr.gy - tr.sy;
                float dist = sqrtf(x*x + y*y);
                gethval(&tr, index_g, &vala);
                float z = vtoi(type_g, vala) - tr.sdepth;
                if ( !x && !y ) {
                    warn(" no azimuth can be calculated (x=0, y=0) for %d-th source of %d-th gather(%s=%d)",
                    isrc + 1, ngather, key, vtoi(type, val));
                } else {
                    float inci = atan2f(dist, z);
                    if (fabs(sinf(inci)) > vw/vp) {
                        warn(" total reflection for %d-th source (offset=%d) of %d-th gather(%s=%d)",
                        isrc + 1, tr.offset/(!tr.scalco ? 1 : ABS(tr.scalco)), ngather, key, vtoi(type, val));
                    } else {
                        azims2r[isrc] = 180/PI*atan2f(x, y);
                        incident[isrc] = 180/PI*asinf(sin(inci) * vp/vw);
                        gethval(&tr, index_a, &vala);
                        azim[isrc] = 0.01*vtoi(type_a, vala);  // OBC cable azimuth
                    }
                }
                if (step == 2) { // refine calculation, retrieve rotation angles
                    azim[isrc] = 0.01*tr.hour;
                    tilt[isrc] = 0.01*tr.minute;
                    roll[isrc] = 0.01*tr.sec;
                } else if (step == 3) {
                    azim[isrc] = 0.01*tr.timbas;
                    tilt[isrc] = 0.01*tr.grnors;
                    roll[isrc] = 0.01*tr.grnlof;
                }
            }
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (ntr%nc) {
                warn("!!! number of traces  not multiple of number of components (%d/%d) for %d-th gather (%s=%d)", ntr, nc, ngather, key, vtoi(type, val));
                warn("!!! cdpt=%d cdp=%d fldr=%d tracf=%d ep=%d offset=%d sx=%d sy=%d",
                        hdrs[0][0].cdpt,  hdrs[0][0].cdp,  hdrs[0][0].fldr,  hdrs[0][0].tracf, hdrs[0][0].ep,  hdrs[0][0].offset,  hdrs[0][0].sx,  hdrs[0][0].sy);
            }

            total += ntr;

            int nsrc = ntr/nc;
            if (verbose) warn("  processing %d traces a la %d component for %d-th gather (%s=%d)", nsrc, nc, ngather, key, vtoi(type, val));
            if ( mode == 0 ) { // fetch angle
                azim3c[2] = hdrs[0][0].fz;
                incl3c[2] = hdrs[0][0].dz;
                azim3c[1] = hdrs[0][0].fy;
                incl3c[1] = hdrs[0][0].dy;
                azim3c[0] = hdrs[0][0].fx;
                incl3c[0] = hdrs[0][0].dx;
            } else { // calculate rotation angle
                if (step >= 2) {
                    refineOrientations(nsrc, itmin, itmax, azims2r, incident,  aw, da, azim, tilt, roll,
                            indata[xaxis], indata[yaxis], indata[zaxis], pol, amp);
                } else {
                    calcOrientations(nsrc, itmin, itmax, itftb, azims2r, incident, azim[0], tilt, roll,
                            indata[xaxis], indata[yaxis], indata[zaxis], pol, amp);
                }
                azim3c[0] = azim[0];
                Azim = azim[0];
                Tilt = tilt[0];
                Roll = roll[0];
                nam = ntm = nrm = 0;
                stddvAzim = stddvTilt = stddvRoll = 0.0;
                if (nsrc > 1) {
                    if (step > 1) Azim = getMedianAngle(nsrc, azim, as[0], &stddvAzim, &nam);
                    Tilt = getMedianAngle(nsrc, tilt, as[1], &stddvTilt, &ntm);
                    Roll = getMedianAngle(nsrc, roll, as[2], &stddvRoll, &nrm);
                    if (step > 1) stddvAzim = stddev(azim, nsrc, Azim);
                    stddvTilt = stddev(tilt, nsrc, Tilt);
                    stddvRoll = stddev(roll, nsrc, Roll);
                }
                calc3COrient(Azim*PI/180.0, Tilt*PI/180.0, Roll*PI/180.0, pol, azim3c, incl3c);
                if (verbose) warn("  Median/StdDev/N%%:  Azimuth=%3.1f/%3.1f/%1.0f  Tilt=%3.1f/%3.1f/%1.0f   Roll=%3.1f/%3.1f/%1.0f",
                        Azim, stddvAzim, 100.0*nam/nsrc, Tilt, stddvTilt, 100.0*ntm/nsrc, Roll, stddvRoll, 100.0*nrm/nsrc);
                if (verbose) warn("  Azimuth/Inclination Inline=%3.1f/%3.1f  Xline=%3.1f/%3.1f  Z=%3.1f/%3.1f\n",
                        azim3c[0], incl3c[0], azim3c[1], incl3c[1], azim3c[2], incl3c[2]);
            }
            if (mode==0 || mode == 2) { // apply rotation
                Project3C(nsrc, 0, nt-1, azim3c, incl3c, pol, indata[xaxis], indata[yaxis], indata[zaxis]);
            }

            int iz = nc - 3;  // index for z component
            for (itr = 0; itr < nsrc; ++itr) {
                for (i = 0; i < nc; ++i) {
                    memcpy(&outtr, &hdrs[i][itr], HDRBYTES);
                    memcpy(outtr.data, indata[i][itr], FSIZE * nt);
                    if (mode > 0) {
                        outtr.lcf  = nsrc;
                        outtr.hcf  = nam;
                        outtr.lcs  = ntm;
                        outtr.hcs  = nrm;
                        outtr.lagb = NINT(10000.0*dt*itftb[itr]);
                        outtr.year = NINT(100.0*azims2r[itr]);
                        outtr.day  = NINT(100.0*incident[itr]);
                        outtr.hour = NINT(100.0*azim[itr]);
                        outtr.minute = NINT(100.0*tilt[itr]);
                        outtr.sec  = NINT(100.0*roll[itr]);
                        outtr.timbas = NINT(100.0*Azim);
                        outtr.grnors = NINT(100.0*Tilt);
                        outtr.grnlof = NINT(100.0*Roll);
                        if(stddvAzim > 0.0) outtr.trwf    = NINT(100.0*stddvAzim);
                        if(stddvTilt > 0.0) outtr.grnofr  = NINT(100.0*stddvTilt);
                        if(stddvRoll > 0.0) outtr.gaps    = NINT(100.0*stddvRoll);
                        if(i >= iz) {
                            outtr.ungpow = (step < 2 && ampi)? amp[itr][i-iz]*pol[i-iz]
                                                             : amp[itr][i-iz];
                        }
                        outtr.fz = azim3c[2];
                        outtr.dz = incl3c[2];
                        outtr.fy = azim3c[1];
                        outtr.dy = incl3c[1];
                        outtr.fx = azim3c[0];
                        outtr.dx = incl3c[0];
                    }
                    puttr(&outtr);
                }
            }
            val = valnew;
            // reset cache
            memset(*hdrs, 0, nc * nxmax * HDRBYTES);
            memset( **indata, 0, nc * nxmax * nt * FSIZE);
            memset(*amp, 0, 3*nxmax * FSIZE);
            memset(azim, 0, nxmax * FSIZE);
            memset(tilt, 0, nxmax * FSIZE);
            memset(roll, 0, nxmax * FSIZE);
            memset(azims2r, 0, nxmax * FSIZE);
            memset(incident, 0, nxmax * FSIZE);
            memset(itftb, 0, nxmax * sizeof(int));
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces for each of %d compomnents are processed", total/nc, nc);
    
    return(CWP_Exit());
}

void refineOrientations(int nsrc, int itmin, int itmax, float* azims2r, float* incident, float* aw, const float da,
        float* azim, float* tilt, float* roll, float** iline, float** xline, float** z, const float* pol, float** amp)
{
    int itr, ns;
    float az, ro, ti, Azim, Tilt, Roll;
    double rms, maxrms;
    float *nord, *east, *vert, *s2r;

    ns = itmax - itmin + 1;
    nord = ealloc1float(ns);
    east = ealloc1float(ns);
    vert = ealloc1float(ns);
    s2r  = ealloc1float(ns);
    for (itr=0; itr<nsrc; ++itr) {
        if (!azims2r[itr] && !incident[itr] && !azim[itr]) continue;  // skip
        maxrms = 0.0;
        for (az = azim[itr] - aw[0]; az <= azim[itr] + aw[0]; az += da) {
            for (ti = tilt[itr] - aw[1]; ti <= tilt[itr] + aw[1]; ti += da) {
                for (ro = roll[itr] - aw[2]; ro <= roll[itr] + aw[2]; ro += da) {
                    rotate3C(itmin, itmax, az, ti, ro, pol, iline[itr], xline[itr], z[itr], nord, east, vert);
                    projectons2r(ns, azims2r[itr], incident[itr], nord, east, vert, s2r);
                    rms = calcRMS(ns, s2r);
                    if (rms > maxrms) {
                        Azim = az;
                        Tilt = ti;
                        Roll = ro;
                        maxrms = rms;
                    }
                }
            }
        }
        // calculate rms of each component
        amp[itr][0] = calcRMS(ns, &z[itr][itmin]);
        amp[itr][1] = calcRMS(ns, &iline[itr][itmin]);
        amp[itr][2] = calcRMS(ns, &xline[itr][itmin]);

        // copy back new refine values
        azim[itr] = Azim;
        tilt[itr] = Tilt;
        roll[itr] = Roll;
    }

    free1float(nord);
    free1float(east);
    free1float(vert);
    free1float(s2r);
    return;
}

void calcOrientations(int nsrc, int itmin, int itmax, int* itftb, float* azims2r, float* incident,
        const float azim, float* tilt, float* roll, float** iline, float** xline, float** z, const float* pol, float** amp)
{
    int itr, it;
    float ampl, maxampl;

    for (itr=0; itr<nsrc; ++itr) {
        if (!azims2r[itr] && !incident[itr] && !azim) continue;  // skip
        // search for max. amplitude of polarization vector because ftb of hydrophone may not good due to phase shift
        maxampl = 0.0;
        for (it = itmin; it <= itmax; ++it) {
            ampl = sqrtf(iline[itr][it]*iline[itr][it] + xline[itr][it]*xline[itr][it] + z[itr][it]*z[itr][it]);
            if ( ampl > maxampl ) {
                maxampl = ampl;
                itftb[itr] = it;
            }
        }
        // output 3C amplitudes for QC plot
        amp[itr][0] = z[itr][itftb[itr]];
        amp[itr][1] = iline[itr][itftb[itr]];
        amp[itr][2] = xline[itr][itftb[itr]];

        calc3C(azims2r[itr], incident[itr], azim, &tilt[itr], &roll[itr], 
                pol[1]*iline[itr][itftb[itr]], pol[2]*xline[itr][itftb[itr]], pol[0]*z[itr][itftb[itr]]);
    }

    return;
}

void calc3C(const float azims2r, const float incident, const float azimx, float* tilt, float* roll, const float x, const float y, const float z)
{
    // implemented according to Olofsson et al. Determining 3C geophone orientation from a single shot
    // SEG 2007 extended abstract

    //const float PI = 3.1416;
    
    float epsilon, kappa, theta, phi;

    float incid =  incident*PI/180;
    float sin_i =  sinf(incid);
    float cos_i =  cosf(incid);
    float gamma =  azimx*PI/180;
    float sin_g =  sinf(gamma);
    float cos_g =  cosf(gamma);
    // determine tilt theta
    //float za = sqrtf(x*x + y*y);
    //      za *= (incident==0.0)? 100.0 : 1.0/tanf(incid);
    float A = sqrtf(x*x + y*y + z*z);
    float a3 = x/A;
    float a4 = y/A;
    float a5 = z/A;
    float bi = sin_i*cosf(azims2r*PI/180);
    float bc = sin_i*sinf(azims2r*PI/180);
    float bv = cos_i;
    float bir =  bi*cos_g + bc*sin_g;
    float bcr = -bi*sin_g + bc*cos_g;
    float bvr = bv;
    float fsb = sqrtf(bir*bir + bvr*bvr);
    if ( fsb > ABS(a3) ) {   // it can happen fsb <= abs(a3)
        epsilon = acosf(a3/fsb);
    } else if (a3 >= 0.0 ) {
        epsilon = 0.0;
    } else { // a3 < 0 )
        epsilon = PI;
    }
    kappa = atan2f(bir, bvr);
    theta = kappa - 0.5*PI + epsilon;  // for OBC

    // nummerically proofed  lhs=a3=rhs
    //float rhs  = bir*cosf(theta) - bvr*sinf(theta);
    //float rhsn = bir*cosf(-theta) - bvr*sinf(-theta);
    //warn(" lhs=%7.6f =? %7.6f=rhs  =?  %7.6f=rhsn ", a3, rhs, rhsn);
    
    // determine roll phi
    float bvtr = bir*sinf(theta) + bvr*cosf(theta);
    phi = atan2f(bvtr*a4 - bcr*a5, bvtr*a5 + bcr*a4);
    // note: phi, the roll angle is defined here conventionally,
    // anti-clockwise from a4 to a5 (crossline to Z, reversed as defined
    // in paper, in which the formula is wrong

    *tilt = theta*180.0/PI;
    *roll = phi*180.0/PI;

    return;
}

void calc3COrient(const float gamma, const float theta, const float phi, const float* pol, float* azim3c, float* incl3c)
{
    float sinr = sinf(gamma);
    float cosr = cosf(gamma);
    float sint = sinf(theta);
    float cost = cosf(theta);
    float sinp = sinf(phi);
    float cosp = cosf(phi);

    azim3c[0] = gamma*180.0/PI;  // inline == cable, Y in internal coord. system           // x comp
    incl3c[0] = -theta*180.0/PI;    // inline inclination, equal tilt of cable
    //if (pol[1] < 0.0) {
    //    azim3c[0] = (azim3c[0] >= 0.0)? azim3c[0] - 180.0 : azim3c[0] + 180.0;
    //    incl3c[0] = -incl3c[0];
    //}
    azim3c[1] = atan2f(cosr*cosp + sinr*sint*sinp, -sinr*cosp + cosr*sint*sinp)*180.0/PI;  // y comp
    incl3c[1] = asinf(cost*sinp)*180.0/PI;
    //if (pol[2] < 0.0) {
    //    azim3c[1] = (azim3c[1] >= 0.0)? azim3c[1] - 180.0 : azim3c[1] + 180.0;
    //   incl3c[1] = -incl3c[1];
    //}
    azim3c[2] = atan2f(-cosr*sinp + sinr*sint*cosp, sinr*sinp + cosr*sint*cosp)*180.0/PI;  // z comp
    incl3c[2] = asinf(cost*cosp)*180.0/PI;
    //if (pol[0] < 0.0) {
    //    azim3c[2] = (azim3c[2] >= 0.0)? azim3c[2] - 180.0 : azim3c[2] + 180.0;
    //    incl3c[2] = -incl3c[2];
    //}
    return;
}

float stddev(float* a, int n, float avg)
{
    int i;
    float diff;
    float asum = 0.0;
    for (i=0; i<n; ++i) {
        diff = a[i] - avg;
        if (ABS(diff) > 180.0) diff = 360.0 - ABS(diff);
        asum += diff*diff;
    }
    return (n>1)? sqrtf(asum/n) : 0.0;
}

float getMedianAngle(int ntr, float* ax, float alimit, float* stddiv, int* n)
{
    int i, itr, na = 0, imaxsector = 0, maxcount = 0 ;
    int nsector = 360/alimit;
    int* count = ealloc1int(nsector);
    memset(count, 0, nsector*sizeof(int));
    for (itr=0; itr<ntr; ++itr) {
        if (ax[itr] == 0.0) continue;  // exclude hard zero, probably not determined value
        i = (int)((ax[itr] + 180.0)/alimit);
        if (i < 0 || i > nsector - 1) continue;
        //if (i < 0) i += nsector;
        //if (i > nsector -1) i -= nsector;
        count[i]++;
    }
    for (i = 0; i < nsector; ++i) {
        if (count[i] > maxcount) {
            imaxsector = i;
            maxcount = count[i];
        }
    }
    float amin = (imaxsector - 1) * alimit - 180.0;
    float amax = (imaxsector + 2) * alimit - 180.0;
    for (itr=0; itr<ntr; ++itr) {
        if (ax[itr] == 0.0) continue;  // exclude hard zero, probably not determined value
        if (imaxsector == nsector - 1 && ax[itr] < -180.0 + alimit) ax[itr] += 360.0;  // -180 ~ -180 + 6
        if (imaxsector == 0 && ax[itr] >  180.0 - alimit) ax[itr] -= 360.0;            // 180 - 6 ~ 180
        if ( ax[itr] < amax && ax[itr] > amin ) ax[na++] = ax[itr];
    }
    if (na == 0) {
        warn("no non-hardzero present, return 0");
        return 0.0;
    }

    float axm = quick_select(ax, na);
    if (stddiv) *stddiv = stddev(ax, na, axm);
    if (n) *n = na;

    if (axm >  180.0) axm -= 360.0;
    if (axm < -180.0) axm += 360.0;

    free1int(count);
    return axm;
}

void rotate3C(int itmin, int itmax, float azim, float tilt, float roll, const float* pol,
        float* x, float* y, float* z, float* x1, float* y1, float* z1)
{
    int it;  // loop counter
    float sina, cosa, sint, cost, sinr, cosr;

    sina = sinf(azim*PI/180.0);
    cosa = cosf(azim*PI/180.0);
    sint = sinf(tilt*PI/180.0);
    cost = cosf(tilt*PI/180.0);
    sinr = sinf(roll*PI/180.0);
    cosr = cosf(roll*PI/180.0);

    for (it=itmin; it<=itmax; ++it) { // loop over samples
        x1[it - itmin] = cosa*cost*x[it]*pol[1] + (cosa*sint*sinr - sina*cosr)*y[it]*pol[2] + (cosa*sint*cosr + sina*sinr)*z[it]*pol[0];
        y1[it - itmin] = sina*cost*x[it]*pol[1] + (sina*sint*sinr + cosa*cosr)*y[it]*pol[2] + (sina*sint*cosr - cosa*sinr)*z[it]*pol[0];
        z1[it - itmin] = -sint*x[it]*pol[1] + cost*sinr*y[it]*pol[2] + cost*cosr*z[it]*pol[0];
    }
    return;
}

void projectons2r(int nt, float azim, float incident, float* x, float* y, float* z, float* s2r)
{
    int it;
    float sina, cosa, sini, cosi;

    sina = sinf(azim*PI/180.0);
    cosa = cosf(azim*PI/180.0);
    sini = sinf(incident*PI/180.0);
    cosi = cosf(incident*PI/180.0);

    for (it=0; it<nt; ++it) { // loop over samples
        s2r[it] = cosa*sini*x[it] + sina*sini*y[it] + cosi*z[it];
    }
    return;
}

void Project3C(const int ntr, const int itmin, const int itmax, const float* azim, const float* incl, const float* pol,
        float** x, float** y, float** z)
{
    int it, itr;
    int nt = itmax - itmin + 1;
    float* x1 = ealloc1float(nt);
    float* y1 = ealloc1float(nt);
    float* z1 = ealloc1float(nt);
    for (itr=0; itr<ntr; ++itr) {
        memset(x1, 0, nt*FSIZE);
        memset(y1, 0, nt*FSIZE);
        memset(z1, 0, nt*FSIZE);
        Project1C(itmin, itmax, azim[0], incl[0], pol[1], x[itr], x1, y1, z1);
        Project1C(itmin, itmax, azim[1], incl[1], pol[2], y[itr], x1, y1, z1);
        Project1C(itmin, itmax, azim[2], incl[2], pol[0], z[itr], x1, y1, z1);
        for (it=itmin; it<=itmax; ++it) {
            x[itr][it] = x1[it - itmin];
            y[itr][it] = y1[it - itmin];
            z[itr][it] = z1[it - itmin];
        }
    }
    free1float(x1);
    free1float(y1);
    free1float(z1);
    return;
}

void Project1C(const int itmin, const int itmax, const float azim, const float incl, const float polarity, const float* x, float* x1, float* y1, float* z1)
{
    int it;  // loop counter
    float xtmp, ytmp, ztmp;

    float sina = sinf(azim*PI/180.0);
    float cosa = cosf(azim*PI/180.0);
    float sini = sinf(incl*PI/180.0);
    float cosi = cosf(incl*PI/180.0);

    for (it=itmin; it<=itmax; ++it) { // loop over samples
        if (x[it] == 0.0) continue;
        ytmp = cosa*cosi*x[it]*polarity;  // project to north, azim=0,  positive y on surface coordinate system
        xtmp = sina*cosi*x[it]*polarity;  // project to east,  azim=90  positive x on surface coordinate system
        ztmp = sini*x[it]*polarity;
        x1[it - itmin] += xtmp;
        y1[it - itmin] -= ytmp;           // reverse polarity of Y, pointing positive south (avoid later rotation problem)
        z1[it - itmin] += ztmp;
    }
    return;
}

double calcRMS(int nt, float* data)
{
    double sum = 0.0;
    int it, nc = 0;
    for (it = 0; it < nt; ++it ) {
        if (data[it] != 0.0 ) {
            ++nc;
            sum += data[it]*data[it];
        }
    }
    return (nc>0)? sqrt(sum/nc): 0.0;
}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float quick_select(float *arr, int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
		if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
		    do ll++; while (arr[low] > arr[ll]) ;
		    do hh--; while (arr[hh]  > arr[low]) ;

		    if (hh < ll)
		    break;

		    ELEM_SWAP(arr[ll], arr[hh]) ;
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]) ;

		/* Re-set active partition */
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}

    return (n%2)? arr[median] : 0.5*(arr[median] + arr[median + 1]);
}

#undef ELEM_SWAP
