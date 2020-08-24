/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPCOMPROT v. 1.0 written by Sanyu Ye */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SP2CROT  -  Compute and apply 2/3 component rotation/orientation for gimbaled ",
"               or fixed 3C geophone along diviated/horizontal well             ",
"                                                                               ",
"   sp2crot <stdin >stdout [options]                                            ",
"                                                                               ",
" Optional Parameters:                                                          ",
"                                                                               ",
"   key=cdpt      receiver gather sort key                                      ",
"   nxmax=512     max. number of shotpoints/traces for a component expected in a gather",
"                                                                               ",
"   nc=3          number of input components in order [H,] C3, C2 and C1, trace by trace",
"                 for VSP in order Z, H1/E, H2/N                                ",
"                 =4 for OBC Hydrophone(H), Vertical(Z), Inline(X) and Crossline(Y)",
"                                                                               ",
"   pol=1,1,1     polarities of C3/Z, C2/X/H1/E, C1/Y/H2/N                      ",
"                                                                               ",
"   axis=-3       rotation axis, count backward from last component             ",
"                 =-3 rotate along C3 (Z) from C1 (Y/H2/N) towards C2 (X/H1/E)  ",
"                 =-2 rotate along C2 (X/H1/E) from C3 (Z) towards C1 (Y/H2/N)  ",
"                 =-1 rotate along C1 (Y/H2/N) from C2 (X/H1/E) towards C3 (Z)  ",
"                                                                               ",
"   mode=0        =0 apply rotation, angle positive from C1 to C2 while along C3, and so further",
"                 =1 compute orientation/rotation angle                         ",
"                 =2 compute orientation and apply rotation                     ",
"                                                                               ",
"   base=0        =0  no base line, usually with mode=0 for just plain rotation ",
"                 =1  use azimuth receiver to source as base line (H1-H2 rotation)",
"                 =2  use inclination angle receiver to source,                 ",
"                     positive upward from vertical (Z-H1 rotation)",
"                 =3  use dihedral angle (non-gimballed H1-H2 rotation along diviated well path)",
"                                                                               ",
"   a=otrav       keyword rotation angle is saved or fetched                    ",
"   azim=igc      keyword holding azimuth of well path at receiver              ",
"   incl=igi      keyword holding inclination of well path at receiver          ",
"   scale=0.1     scaling factor to convert angle stored in header in degrees   ",
"                                                                               ",
"   gdepth=gelev  keyword holding geophone depth, scaled by scalco as usual     ",
"   sdepth=sdepth keyword holding source depth, scaled by scalco as usual       ",
"                                                                               ",
" Parameters used to control applying rotation:                                 ",
"                                                                               ",
"   vsp=0         =-1 apply back projection onto ZEN for microseismic           ",
"                 =0  apply plain rotation without base line                    ",
"                 =1  rotate 1st component C1/H2 with max. energy               ",
"                 =2  rotate 2nd component C2/H1 with max. energy               ",
"                 =3  rotate in TAU-P domain (Z-H1 rotation)                    ",
"   v=vevel [m/s] keyword containing velocity at receiver level for computing   ",
"                   rotation angle=asin(P*V) for vsp=3. P is fetched from tr.fx ",
"                                                                               ",
" Parameters used to determine orientation angle of first of the two components:",
"                                                                               ",
"   t0=1.0   [s]  time line peaks of first arrivals are aligned on              ",
"                 i.e. FTB for VSP, first positive peak for OBC data            ",
"   t1=0.98  [s]  start of time window to determine rotation angles             ",
"                 should be set possibly close to onset of first arrivals       ",
"   t2=1.04  [s]  end of time window to determine rotation angle                ",
"                                                                               ",
"   da=0.1  [deg] increment of rotation angle                                   ",
"   as=30.0 [deg] angle (size) of sector to compute median value of orientation ",
"                                                                               ",
"   verbose=0     =0 quiet mode; >0 output info                                 ",
"                                                                               ",
" Note:                                                                         ",
"   C1, C2, C3 or H2, H1, Z must build a right hand Cartesian coordinate system.",
"   For OBC input data must be sorted by hydrophone (H), if present, vertical (Z),",
"   inline (X) and crossline (Y) while for VSP in Z, H1/E, H2/N.                ",
"   Z (C3) axis is positive downwards, so the H1/E (inline) component that normally",
"   coincides with surface X axis is treated by this program as y (C2) internally",
"   while surface Y i.e. H2/N (crossline) as x (C1).                            ",
"   If only 2 components input, first trace is considered as C2.                ",
"   In case of multiple sources in a gather, distribution of determined orientation",
"   angles is computed. The sector with highest concentration and its two neighbouring",
"   sectors are used to determine the median value.                             ",
"   The program computes rms amplitude of the given window while rotate from first",
"   (x) to second axis (y). The orientation of x is determined separately by max",
"   energy on x, y and ratio y/x. Normally all three results should give same   ",
"   values.                                                                     ",
"   For gimbaled geophone/vertical well the orientation determined here (base=1)",
"   refers to azimuthal angles of component C1 (N/Y/H2 on surface). The values  ",
"   determined by individual shot are stored in keyword fx, while their median  ",
"   values in keyword specified by a=keyword.                                   ",
"   For deviated or even horizontal well with non-gimbaled geophone the C3/Z is ",
"   refered to the inline component along well path. The azimuth and inclination",
"   of well path must be loaded in the keywords specified by (azim= and incl=)  ",
"   for every receiver. Note the inclination is computed from vertical, positive",
"   upnward, in accord with diviation of well from vertical as usually given for",
"   well path.                                                                  ",
"   To enable 2C rotation analysis along Z for deviated well path, two coaxial  ",
"   planes are constructed: A vertical plane through the well path at receiver  ",
"   is used as reference plane, another one through polarization vector and the ",
"   same well path. The angle between these two planes is called dihedral angle ",
"   and used as base angle. The azimuthal angle determined by 2C rotation analysis",
"   is refered to the vertical reference plane instead the North. This angle is ",
"   transformed back to the orientations (azimuth/inclination) of C1/C2/C3, and ",
"   their values are stored in keyword fy/dy, fx/dx and fz/dz.                  ",
"   The parameter vsp= controls how rotation is applied. In plain rotation mode ",
"   (mode=0 vsp=0), a right-hand rotation by an angle fetched from keyword as set",
"   by parameter a=keyword is performed along the axis specified by the parameter",
"   axis=. In case of usual rotation of horizontal components, C1/C2 are rotated",
"   clockwise, viewd from top. To get result as if C2/H1 and C1/H2 are along East",
"   and North, as for microseismic study, back projection (vsp=-1) should be used.",
"   For deviated well (base=3) the back projection to ZEN (vsp=-1) is the only  ",
"   meaningful option. Here the orientation angles (azim/incl) are fetched from ",
"   keyword fy/dy, fx/dx and fz/dz, and the output traces C3/C2/C1 correspond to",
"   Z/E/N.                                                                      ",
"   Polarities of C3/C2/C1 can be corrected on the fly by parameter (pol=) before",
"   rotation analysis/application. Only in case of deviated well (base=3) the   ",
"   polarity of component Z is checked and warning given if it does not correspond",
"   its naturally recorded (positive downwards).                                ",
"   FTB must be strictly aligned at same level (t0=) before rotation analysis.  ",
"   Correct polarities of 3C geophone are essential for correct analysis result ",
"   and particularly for back projection in case of deviated well.              ",
"   3C geophones of Avalon VSP tool string that we are using is constructed in  ",
"   left hand system, therefore C1/H2/Y component should be flipped before rotation",
"   analysis (pol=1,1,-1). Also the FTB polarity of Z component of Avalon tools ",
"   is negative, and needs to be flipped when applying back projecttion/rotation",
"   for the deviated well (warning is given during analysis). Vertical well is a",
"   special case of deviated well (azimuth & inclination equal zero). This also ",
"   means that the polarities of Z components can be checked by false deviated  ",
"   well (base=3) after resetting to zero the header keywords specified by      ",
"   parameters azim= and incl=.                                                 ",
"                                                                               ",
" Examples:                                                                     ",
" 1. align FTB at 1 sec and compute orientation angle of H2 for VSP             ",
" suchw key1=tstat key2=gwdep b=0.1 a=-1000 < gimbaledZH1H2.su |\\",
" sustatic hdrs=1 |\\",
" sp2crot mode=2 nc=3 axis=-3 vsp=2 key=gdel > ZRT.su",
"                                                                               ",
" 2. above example run in mode deviated well and back project to ZEN            ",
" suchw key1=tstat key2=gwdep b=0.1 a=-1000 < gimbaledZH1H2.su |\\",
" sustatic hdrs=1 |\\",
" sushw key=igc,igi a=0,0 |\\",
" sp2crot mode=2 nc=3 axis=-3 base=3 vsp=-1 key=gdel azim=igc incl=igi > ZEN.su",
"                                                                               ",
" 3. load azimuth/inclination for well path before analysis/application         ",
"    note to limit accuracy of input azim/incl to 1/10 degree when using short 16bit keywords",
" !!!NOTE: it is very tricky to prepare azimuth to avoid interpolation error ",
"    when two adjacent values cross zero line (0 ~ 360) or 180 (-180 ~ 180)",
" sushw match=offset key=minute,sec infile=azim-incl.tbl scale=10 < ZXY.su interp=1 |\\",
" sp2crot mode=2 nc=3 axis=-3 key=offset base=3 vsp=-1 azim=minute incl=sec scale=0.1 > ZEN.su",
"                                                                               ",
"                                                                               ",
" Version 3.0.0 last modified  Aug. 2012 by Sanyu Ye                            ",
NULL};

/* Credits:
 *	Sanyu Ye: READ Well Services, May. 2010
 */
/**************** end self doc *******************************************/

#include "sprinthlp.c"

// forward declaration
float calcRot(int ntr, int itftb, int itmin, int itmax, float da,
               float* baseangle, float** x, float** y, float** angles);
void Rotate(const int ntr, const int nt, const int vsp, const float* baseangle, const float angle, float** x, float** y);
float calcRMS(int nt, float* data);
float getMedianAngle(int ntr, float* ax, float alimit, float* stddiv, int* n);
float stddev(float* a, int n, float avg);
float calcDihedral(const float x, const float y, const float z, const float azim, const float incl, float* zPol);
void calc3COrient(const float azim, const float incl, const float a3, float* azim3c, float* incl3c);
void Project3C(const int ntr, const int itmin, const int itmax, const float* azim, const float* incl, float** x, float** y, float** z);

int verbose;		/* =0 no info, >0 output info  */

int main(int argc, char **argv) {
    cwp_String key, key_a, key_g, key_s, key_az, key_in, key_v;           /* header key word from segy.h		*/
    cwp_String type, type_a, type_g, type_s, type_az, type_in, type_v;      /* type of key				*/
    Value val, valnew, vala;                                /* value of key			*/
    int index, index_a, index_g, index_s, index_az, index_in, index_v;        /* index of key				*/
    int nt, it, itftb, itmin, itmax;       /* nums amps as int			*/
    int axis, xaxis, yaxis, zaxis, base;
    int nc, nxmax;                     /* number of components			*/
    int nsegy;                  /* length bytes read for segy trace	*/
    int mode, vsp;
    float a, asector, da, azim, incl;                    /* rotation angle computed or applied   */
    float t0, t1, t2, dt, scale;                /* input trace px py header values */
    float am[3], stddv[3], dv[3];
    float azim3c[3], incl3c[3], pol[3];   /* orientation angle computed for x/y/z (N/E/Z-downward) in degrees (-180 - 180)   */
    float v, pv;
    int nam[3];
    int i, itr, ntr, nsrc, ngather, total;
    segy tr, outtr;

    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose=0;

    /* get keys */
    if (!getparstring("a", &key_a))	 key_a="otrav";
    type_a = hdtype(key_a);
    index_a = getindex(key_a);
    if (!getparstring("azim", &key_az))	 key_az="igc";
    type_az = hdtype(key_az);
    index_az = getindex(key_az);
    if (!getparstring("incl", &key_in))	 key_in="igi";
    type_in = hdtype(key_in);
    index_in = getindex(key_in);
    if (!getparstring("gdepth", &key_g)) key_g="gelev";
    type_g = hdtype(key_g);
    index_g = getindex(key_g);
    if (!getparstring("sdepth", &key_s)) key_s="sdepth";
    type_s = hdtype(key_s);
    index_s = getindex(key_s);
    if (!getparstring("v", &key_v)) key_v="wevel";
    type_v = hdtype(key_v);
    index_v = getindex(key_v);
    if (!getparstring("key", &key))	 key="cdpt";
    type = hdtype(key);
    index = getindex(key);

    if (!getparfloat("scale", &scale))  scale = 0.1;
    if (!getparfloat("as", &asector))   asector = 30.0;
    if (!getparfloat("da", &da))        da = 0.1;
    if (!getparint("nxmax", &nxmax))    nxmax = 512;

    if (!getparint("nc", &nc))        nc = 3;
    if (!getparint("mode", &mode))  mode = 0;
    if (!getparint("axis", &axis))  axis = -3;
    else {
        if (axis >=0 || axis < -3 || axis < -nc)
            err("invalid axis (=%d), should be between -1 and -3, count backwards", axis);
    }
    if (!getparint("base", &base)) {
        base = 0;
    } else {
        if ( base < 0 || base > 3 )
            err("invalid base line (=%d), should be none (=0) or azimuth (=1) or incident (=2) or dihedral angle (=3)", base);
    }
    if (!getparint("vsp", &vsp))     vsp = 0;
    if (vsp == 3) base = 0;  // reset base line angle to all zero
    if ( vsp < -1 || vsp > 3 )  err("invalid vsp option (-1 <= vsp=%d <= 3)", vsp);

    if ( (mode == 0 || mode == 2) && vsp == -1 && verbose) {
        if (base == 3) warn("Back projection (vsp=-1) read 3C orientation from fy/dy fx/dy fz/dz");
        else           warn("Back projection (vsp=-1) rotates back C1/H2 and C2/H1 to N and E");
    }

    int n = countparval("pol");
    if ( n > 0 && n <= 3 ) {
        getparfloat("pol", pol);
        for (i = n; i < 3; ++i) pol[i] = pol[n-1];
    }
    else {
        {pol[0] = 1.0; pol[1] = 1.0; pol[2] = 1.0;}
    }

    // find out component index corresponding to x and y axis
    if (nc == 2) {
        xaxis = 1;
        yaxis = 0;
    } else {
        xaxis = nc + axis - 1;
        yaxis = nc + axis - 2;
        zaxis = nc + axis - 3;
        if (xaxis < nc - 3) xaxis += nc;
        if (yaxis < nc - 3) yaxis += nc;
        if (zaxis < nc - 3) zaxis += nc;
    }

    /* get first trace */
    if ( (nsegy=gettr(&tr)) < HDRBYTES) err("can't get first trace");
    gethval(&tr, index, &val);

    dt = ((double) tr.dt) / 1000000.0;
    nt = tr.ns;

    if (!getparfloat("t0", &t0)) t0 = 1.0;
    if (!getparfloat("t1", &t1)) t1 = 0.98;
    if (!getparfloat("t2", &t2)) t2 = 1.04;
    itftb = NINT(t0 / dt);
    itmin = NINT(t1 / dt);
    itmax = NINT(t2 / dt);
    if (mode > 0) { // Check time gating values for analysis mode
        if (itmin < 0)
            err("itmin=%d should be positive", itmin);
        if (itmin >= itmax)
            err("itmin=%d, itmax=%d conflict", itmin, itmax);
        if (itftb < itmin || itftb > itmax)
            err("t0 %5.3f should be larger than t1 %5.3f but smaller than t2 %5.3f", t0, t1, t2);
        if (tr.ns <= itmax)
            err("tr.ns=%d, itmax=%d window cannot extend over the trace length", tr.ns, itmax);
    }
    
    segyhdr** hdrs = (segyhdr **) ealloc2(nxmax, nc, HDRBYTES);
    float***  indata = ealloc3float(nt, nxmax, nc);
    float*** outdata = ealloc3float(nt, nxmax, nc);
    float** baseangle = ealloc2float(nxmax, 4);
    float** angles = ealloc2float(nxmax, 3);
    float*  zPol = ealloc1float(nxmax);
    
    // reset with zero
    memset(*hdrs, 0, nc * nxmax * HDRBYTES);
    memset( **indata, 0, nc * nxmax * nt * FSIZE);
    memset(**outdata, 0, nc * nxmax * nt * FSIZE);
    memset(*baseangle, 0, 4 * nxmax * FSIZE);
    memset(*angles, 0, 3 * nxmax * FSIZE);
    memset(zPol, 0, nxmax * FSIZE);

    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = total = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nc*(nxmax - 1)) err("\nNumber of traces exceeding nxmax=%d\n", nxmax);
            int isrc = ntr/nc;
            int icmp = ntr%nc;
            memcpy(&hdrs[icmp][isrc], &tr, HDRBYTES);
            for (it = 0; it < nt; ++it) indata[icmp][isrc][it] = pol[icmp]*tr.data[it];
            if (!icmp) {
                float y = tr.sx - tr.gx;  // internal coordinate x-N, y-E, z-Z downwards
                float x = tr.sy - tr.gy;
                float dist = sqrtf(x*x + y*y);
                gethval(&tr, index_g, &vala);
                float zg = vtoi(type_g, vala);
                gethval(&tr, index_s, &vala);
                float zs = vtoi(type_s, vala);
                float z = zs - zg;
                baseangle[1][isrc] = 180/PI*atan2(y, x); // receiver to source azimuth
                baseangle[2][isrc] = 180/PI*atan2(z, dist); // // receiver to source inclination
                if (base == 3) {
                    gethval(&tr, index_az, &vala);
                    azim = scale*vtof(type_az, vala);
                    gethval(&tr, index_in, &vala);
                    incl = 90.0 - scale*vtof(type_in, vala);
                    baseangle[3][isrc] = calcDihedral(x, y, z, azim, incl, &zPol[isrc]);
                }
            }
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            int keyno = vtoi(type, val);

            ++ngather;
            if (verbose > 1) warn(" Processing %d traces of %d-th gather (%s=%d) ...", ntr, ngather, key, keyno);
            if (ntr%nc) err(" number of traces not multiple of number of components for %d-th gather (%s=%d)",
                    ntr, ngather, key, keyno);

            nsrc = ntr/nc;
            total += ntr;

            stddv[0] = stddv[1] = stddv[2] = 0.0;
            if ( mode == 0 ) { // fetch angle
                if (base == 3) {
                    azim3c[2] = hdrs[0][0].fz;
                    incl3c[2] = 90.0 - hdrs[0][0].dz;
                    azim3c[1] = hdrs[0][0].fy;
                    incl3c[1] = 90.0 - hdrs[0][0].dy;
                    azim3c[0] = hdrs[0][0].fx;
                    incl3c[0] = 90.0 - hdrs[0][0].dx;
                } else if (vsp == 3) { // compute angle from ray parameter p=tr.fx and velocity contained in key_v
                    gethval((segy*) hdrs[0], index_v, &vala);
                    v = 0.001*vtod(type_v, vala);  //  convert to km/s
                    pv = -tr.fx*v;  // negate to account for natural polarity of X (positive on -X or -P side)
                    if ( ABS(pv) > 1.0) {
                        if (verbose) warn (" rotation angle reset to 90 degree for P=%8.6f V=%4.0f of %d-th gather (%s=%d)",
                                tr.fx, 1000.0*v, ngather, key, keyno);
                        pv = (pv > 0)? 1.0 : -1.0;
                    }
                    a = asinf(pv) * 180.0 / PI;
                    vala.h = (short) NINT(a/scale);
                } else {
                    gethval((segy*) hdrs[0], index_a, &vala);
                    a = vtod(type_a, vala)*scale;
                }
            } else {
                // check z polarity
                if (nc > 2 && base == 3) {
                    int nFalsePol = 0;
                    for (i = 0; i < nsrc; ++i) {
                        if (indata[zaxis][i][itftb]*zPol[i] < 0.0) ++nFalsePol;
                    }
                    if (verbose) warn("  %d out of total %d traces of Z (inline to well) have wrong apparent polarity", nFalsePol, nsrc);
                }
                // do calculate rotation angle
                calcRot(nsrc, itftb, itmin, itmax, da, baseangle[base], indata[xaxis], indata[yaxis], angles);
                int iFlip = 0, nShotUsed = 0;
                //float maxdv = 360.0;
                for (i=0; i<3; ++i) {
                    am[i] = getMedianAngle(nsrc, angles[i], asector, &stddv[i], &nam[i]);
                    dv[i] = (nam[i] > 0) ? stddv[i]/nam[i] : 360.0;
                    if (nam[i] > nShotUsed) {
                        nShotUsed = nam[i];
                        iFlip = i;
                    }
                }

                a = am[iFlip];

                vala.h = (short) NINT(a/scale);
                //vala.i = (int) (a/scale);
                if (verbose > 1) {
                    fprintf(stderr,     "  Shot   Baselines: Azimuth/Inclination/Dihedral  Orientation for  C1  /   -C1   /   -C2\n");
                    for (i = 0; i < nsrc; ++i) {
                        fprintf(stderr, "  %3d              %7.2f /   %6.2f  / %7.2f               %7.2f / %7.2f / %7.2f\n", \
                            i+1, baseangle[1][i], 90.0 - baseangle[2][i], (base == 3)? baseangle[3][i] : 0.0, angles[0][i], angles[1][i], angles[2][i]);
                    }
                }

                if (verbose) fprintf(stderr, " Median/StdDev/nShots of C1, -C1 and -C2: %4.2f/%4.2f/%d %4.2f/%4.2f/%d %4.2f/%4.2f/%d\n",
                        am[0], stddv[0], nam[0], am[1], stddv[1], nam[1], am[2], stddv[2], nam[2]);
                if (verbose) {
                    fprintf(stderr, " Likely Orientation/Azimuth of C%d = %.2f (degrees) ", nc - xaxis, a);
                    if (iFlip>0) fprintf(stderr, "of %s flipped ", (iFlip == 1)? "C1/Y" : "C2/X");
                    fprintf(stderr, "for %d-th gather (%s=%d)\n",  ngather, key, keyno);
                }
                        
                
                if (base == 3) { // for dipping well
                    if (verbose) fprintf(stderr, "  relative to well track orientation Azimuth = %.2f , Inclination = %.2f (degrees)\n", azim, 90.0 - incl);
                    calc3COrient(azim, incl, a, azim3c, incl3c);
                    if (verbose) warn("  Azimuth/Inclination N/C1/Y=%3.1f/%3.1f  E/C2/X=%3.1f/%3.1f  C3/Z=%3.1f/%3.1f\n",
                            azim3c[1], 90.0 - incl3c[1], azim3c[0], 90.0 - incl3c[0], azim3c[2], 90.0 - incl3c[2]);
                }
            }

            if (mode == 0 || mode == 2) { // apply rotation
                if (base == 3) {
                    Project3C(nsrc, 0, nt - 1, azim3c, incl3c, indata[xaxis], indata[yaxis], indata[zaxis]);
                } else {
                    Rotate(nsrc, nt, vsp, baseangle[base],  a, indata[xaxis], indata[yaxis]);
                }
            }

            for (itr = 0; itr < nsrc; ++itr) {
                for (i = 0; i < nc; ++i) {
                    memcpy(&outtr, &hdrs[i][itr], HDRBYTES);
                    memcpy(outtr.data, indata[i][itr], FSIZE * nt);
                    if (mode > 0) {
                        outtr.gaps    = NINT((base != 2 ? baseangle[base][itr] : 90.0 - baseangle[base][itr])/scale);
                        outtr.grnors  = NINT(angles[0][itr]/scale);
                        outtr.grnofr  = NINT(am[0]/scale);
                        outtr.grnlof  = NINT(stddv[0]/scale);
                        if (base == 3) {
                            outtr.fz = azim3c[2];
                            outtr.dz = 90.0 - incl3c[2];
                            outtr.fy = azim3c[1];
                            outtr.dy = 90.0 - incl3c[1];
                            outtr.fx = azim3c[0];
                            outtr.dx = 90.0 - incl3c[0];
                        } else {
                            outtr.fx = am[0];
                            outtr.dx = stddv[0];
                            outtr.fy = angles[0][itr];
                            outtr.dy = nsrc;
                            outtr.fz = baseangle[1][itr];
                            outtr.dz = 90.0 - baseangle[2][itr];
                        }
                        puthval(&outtr, index_a, &vala);
                    } else if (vsp == 3) {
                        puthval(&outtr, index_a, &vala);
                    }
                    puttr(&outtr);
                }
            }
            
            val = valnew;
            // reset output data
            memset(**indata, 0, nc * nxmax * nt * FSIZE);
            memset(*hdrs, 0, nc * nxmax * HDRBYTES);
            memset(*baseangle, 0, 4 * nxmax * FSIZE);
            memset(*angles, 0, 3 * nxmax * FSIZE);
            memset(zPol, 0, nxmax * FSIZE);
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces for each of %d compomnents are processed", total/nc, nc);
    
    return(CWP_Exit());
}

float calcRot(int nsrc, int itftb, int itmin, int itmax, float da, float* baseangle, float** x, float** y, float** angles)
{
    int itr;
    float a, axftb, ax, a1, a2;
    float xrms, xrmsmax;

    int nt = itmax - itmin + 1;
    float* x1 = ealloc1float(nt);
    float* y1 = ealloc1float(nt);

    for (itr=0; itr<nsrc; ++itr) {
        xrmsmax = 0;
        for (a = 0; a <= 180.0; a += da) {
            Rotate2C(itmin, itmax, 1, a, x[itr], y[itr], x1, y1);
            xrms = calcRMS(nt, x1);
            if (xrms > xrmsmax ) {
                axftb = x1[itftb-itmin];
                ax = baseangle[itr] - a + 180;
                xrmsmax = xrms;
            }
        }
        
        //if (axftb < 0) ax -= 180.0;
        //if (ax >  180) ax -= 360.0;
        //if (ax < -180) ax += 360.0;
        if (axftb < 0) ax -= 180.0;
        if (ax >  180) ax -= 360.0;
        if (ax < -180) ax += 360.0;
       
        // calc angle as if C1 is flipped s2r - ax = 180 - (s2r - a1)  ; s2r = baseangle + 180
        a1 = -ax + 2.0*baseangle[itr] + 180;
        if (a1 >  180) a1 -= 360.0;
        if (a1 >  180) a1 -= 360.0;
        if (a1 < -180) a1 += 360.0;
        if (a1 < -180) a1 += 360.0;
        // calc angle as if C2 is flipped s2r - ax = -(s2r - a1)
        a2 = -ax + 2.0*baseangle[itr];
        //a2 = 2.0*ax - s2r;
        if (a2 >  180) a2 -= 360.0;
        if (a2 >  180) a2 -= 360.0;
        if (a2 < -180) a2 += 360.0;
        if (a2 < -180) a2 += 360.0;

        angles[0][itr] = ax;
        angles[1][itr] = a1;
        angles[2][itr] = a2;
    }

    free1float(x1);
    free1float(y1);

    return 0.0;
}


float getMedianAngle(int ntr, float* ai, float alimit, float* stddiv, int* n)
{
    int i, itr, na = 0, imaxsector = 0, maxcount = 0 ;
    int nsector = 360/alimit;
    int* count = ealloc1int(nsector);
    memset(count, 0, nsector*sizeof(int));
    float* ax = ealloc1float(ntr);
    memcpy(ax, ai, ntr*FSIZE);
    float* a = ealloc1float(ntr);
    memset(a, 0, ntr*FSIZE);

    for (itr=0; itr<ntr; ++itr) {
        //if (ax[itr] == 0.0) continue;  // exclude hard zero, probably not determined value
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
        if ( ax[itr] < amax && ax[itr] > amin ) a[na++] = ax[itr];
    }
    if (na == 0) {
        warn("no non-hardzero present, return 0");
        return 0.0;
    }

    float axm = quick_select(a, na);
    if (stddiv) *stddiv = stddev(a, na, axm);
    if (n) *n = na;

    if (axm >  180.0) axm -= 360.0;
    if (axm < -180.0) axm += 360.0;

    free1int(count);
    free1float(a);
    free1float(ax);
    return axm;
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


// rotate x y according to angle a
void Rotate(const int ntr, const int nt, const int vsp, const float* baseangle, const float angle, float** x, float** y)
{
    int itr;  // loop counter
    float a = angle;

    for (itr=0; itr<ntr; ++itr) { // loop over traces
        if (vsp == -1) {
            a = -angle;
        } else if (vsp == 1) {
            a = 180.0 + baseangle[itr] - angle;
        } else if (vsp == 2) {
            a =  90.0 + baseangle[itr] - angle;
        }
        Rotate2C(0, nt, 0, a, x[itr], y[itr], x[itr], y[itr]);
    }
}

float calcDihedral(const float x, const float y, const float z, const float azim, const float incl, float* zPol)
{
    float ax, ay, bx, by, bz, angle;
    float cx, cy, cz, dx, dy, dz, d;
    float cx2, cy2, cz2, cc, zz;
    
    // calc directional cos a for vertical plane along well path
    // ax*x + ay*y + az*z + d1 = 0

    ay =  cosf(azim*PI/180.0);
    ax = -sinf(azim*PI/180.0);
    // c1 = d1 = 0.0;

    // calc directional cos for well path (line)
    cx = cosf(azim*PI/180.0)*cosf(incl*PI/180.0);
    cy = sinf(azim*PI/180.0)*cosf(incl*PI/180.0);
    cz = sinf(incl*PI/180.0);

    // calc directional cos for receiver to source (line)
    d = sqrtf(x*x + y*y + z*z);
    dx = x/d;
    dy = y/d;
    dz = z/d;

    // calc dot product of reversed ray (r2s) and well track
    zz = dx*cx + dy*cy + dz*cz;
    *zPol = (zz > 0.0)? -1.0 : 1.0;

    // calc normal directional cos b for plane through source and well path
    // b = c x d = | (c2d3 - c3b2)  (c3d1 - c1d3)  (c1d2 - c2d1) |   cross product

    bx = cy*dz - cz*dy;
    by = cz*dx - cx*dz;
    bz = cx*dy - cy*dx;
    // should be d x c = - c x d
    //bx = -bx;
    //by = -by;
    //bz = -bz;

    // dehydral angle
    angle = acosf((ax*bx + ay*by)/sqrt(bx*bx + by*by + bz*bz));

    // the angle is smaller 180, use cross product to figure out rotation sense
    // c2 = a x b
    cx2 = ay*bz;
    cy2 = -ax*bz;
    cz2 = ax*by - ay*bx;

    cc = cx*cx2 + cy*cy2 + cz*cz2;

    if (cc < 0.0) angle = -angle;

    return angle*180.0/PI;
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

float calcMax(int nt, float* data)
{
    float max = 0;
    int it;
    for (it = 0; it < nt; ++it ) {
        if (data[it] != 0.0 ) {
            if ( ABS(data[it]) > max) max = ABS(data[it]);
        }
    }
    return max;
}

void calc3COrient(const float azim, const float incl, const float a3, float* azim3c, float* incl3c)
{
    // a3 last rotation angle along z (downwards but dip with an angle=a2 to vertical)
    // a2 second rotation along y (X for data) from z (vertical downwards)
    // a1 first rotation along z (vertical downwards) the azimuth of well track at receiver

    //  x   | cosa1  -sina3    0 | | cosa2   0   sina2 | | cosa3  -sina3   0  | | x3 |
    //  y = | sina1   cosa1    0 |*|   0     1     0   |*| sina3   cosa3   0  |*| y2 |
    //  z   |   0       0      1 | | -sina2  0   cosa2 | |   0       0     1  | | z1 |

    //  x   | cosa1cosi1   cosa2cosi2   cosa3cosi3 | | x' |
    //  y = | sina1cosi1   sina2cosi2   sina3cosi3 |*| y' |
    //  z   |    sini1        sini2        sini3   | | z' |


    // a1 = azim  well track azimuth, a2 = 90 - incl  (incl = well track inclination)
    float a1 = azim;
    float a2= 90.0 - incl;

    float sina1 = sinf(a1*PI/180.0);
    float cosa1 = cosf(a1*PI/180.0);
    float sina2 = sinf(a2*PI/180.0);
    float cosa2 = cosf(a2*PI/180.0);
    float sina3 = sinf(a3*PI/180.0);
    float cosa3 = cosf(a3*PI/180.0);

    azim3c[2] = azim;  // z comp, positive downwards
    incl3c[2] = incl;   // inline (to well) inclination, equal 90 - angle (diviation of well from vertical)

    azim3c[1] = atan2f(sina1*cosa2*cosa3 + cosa1*sina3 , cosa1*cosa2*cosa3 - sina1*sina3)*180.0/PI;    // N/X/y, X in internal but y on surface coord. system
    incl3c[1] = asinf(-cosa3*sina2)*180.0/PI;

    azim3c[0] = atan2f(-sina1*cosa2*sina3 + cosa1*cosa3 , -cosa1*cosa2*sina3 - sina1*cosa3)*180.0/PI;  // E/Y/x, Y in internal but x on surface coord. system
    incl3c[0] = asinf(sina2*sina3)*180.0/PI;    

    return;
}

void Project3C(const int ntr, const int itmin, const int itmax, const float* azim, const float* incl, float** x, float** y, float** z)
{
    //  x   | cosa1cosi1   cosa2cosi2   cosa3cosi3 | | x' |
    //  y = | sina1cosi1   sina2cosi2   sina3cosi3 |*| y' |
    //  z   |    sini1        sini2        sini3   | | z' |

    int it, itr;
    float xtmp, ytmp, ztmp;

    float sina1 = sinf(azim[1]*PI/180.0);
    float cosa1 = cosf(azim[1]*PI/180.0);
    float sini1 = sinf(incl[1]*PI/180.0);
    float cosi1 = cosf(incl[1]*PI/180.0);
    float sina2 = sinf(azim[0]*PI/180.0);
    float cosa2 = cosf(azim[0]*PI/180.0);
    float sini2 = sinf(incl[0]*PI/180.0);
    float cosi2 = cosf(incl[0]*PI/180.0);
    float sina3 = sinf(azim[2]*PI/180.0);
    float cosa3 = cosf(azim[2]*PI/180.0);
    float sini3 = sinf(incl[2]*PI/180.0);
    float cosi3 = cosf(incl[2]*PI/180.0);

    for (itr=0; itr<ntr; ++itr) {
        for (it=itmin; it<=itmax; ++it) {
            xtmp = cosa1*cosi1*x[itr][it] + cosa2*cosi2*y[itr][it] + cosa3*cosi3*z[itr][it];
            ytmp = sina1*cosi1*x[itr][it] + sina2*cosi2*y[itr][it] + sina3*cosi3*z[itr][it];
            ztmp = sini1*x[itr][it] + sini2*y[itr][it] + sini3*z[itr][it];
            x[itr][it] = xtmp;
            y[itr][it] = ytmp;
            z[itr][it] = ztmp;
        }
    }
    return;
}

