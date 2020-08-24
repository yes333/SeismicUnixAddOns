
/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPRELOBC v. 1.0 written by Sanyu Ye */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPRELOBC  -  RELocate OBC position using FTB of hydrophone component          ",
"                                                                               ",
"   sprelobc <stdin >stdout [options]                                          ",
"                                                                               ",
" Optional Parameters:                                                          ",
"                                                                               ",
"   key=cdpt      receiver gather sort key                                      ",
"   nxmax=1024    max. number of traces for a component expected in a gather    ",
"                                                                               ",
"   t0=1.0  [s]   time FTBs are aligned on. If set to zero, FTB is not corrected by ",
"                 nominal traveltime sqrt((sx-gx)^2+(sy-gy)^2+(sdepth-gwdep)^2)/v0",
"   td1=20  [ms]  min. time of airgun delay for search window                   ",
"   td2=60  [ms]  max. time of airgun delay for search window                   ",
"   dtd=tr.dt     search increment                                              ",
"                                                                               ",
"   ftb=laga      keyword FTB is fetched from                                   ",
"   scale=0.0001  scaling factor to convert FTB to second                       ",
"                                                                               ",
"   gdepth=gwdep  keyword holding geophone depth, scaled by scalco as usual     ",
"                                                                               ",
"   xy=25   [m]   max position correction of x and y                            ",
"   dxy=0.1 [m]   step of position shift for x and y                            ",
"                 if not given, calculated based on averaged delays             ",
"   z=0     [m]   max position correction of z (water depth)                    ",
"   dz=1.0  [m]   step of position shift for z                                  ",
"                                                                               ",
"   s=50    [m]   receiver spacing, used to constrain distance from last one    ",
"   smin=s - 2    min. distance allowed                                         ",
"   smax=s + 0.1  max. distance allowed                                         ",
"                                                                               ",
"   v0=1480 [m/s] water velocity used to calculate norminal traveltime          ",
"   vmin=v0       min velocity                                                  ",
"   vmax=v0       max velocity                                                  ",
"   dv=1          search increment                                              ",
"                                                                               ",
"   verbose=0     >0 output info                                                ",
"                                                                               ",
" Note:                                                                         ",
"   The sophistication of determining all possible parameters should be used with",
"   care bearing in mind the brute-force approach behind: it can take hours to run",
"   through millions of iterations when all paramters are actively sought. First",
"   the default setting with more coarse dxy and dtd should be used to find the ",
"   average airgun delay and then set td1 and td2 to the same value to search for",
"   optimal water velocity because these two parameters should be fixed for entire",
"   survey. In general coarse search steps should be used to narrow search windows",
"   before using fine steps.                                                    ",
"   The new posisition and its correction against the old are stored in keywords",
"   fx,dx,fy,dy,fz,dz. Furthermore the distance from last receiver and the cumulative",
"   cable length are stored in ungpow and unscale.                              ",
"                                                                               ",
" Examples:                                                                     ",
" 1. align water wave at 1 sec, find FTB and look for average airgun delay time ",
" sureduce mode=2 rv=1.48 t0=1.0 < input.su |\\                                 ",
" sushw key=laga a=10000 |\\                                                    ",
" spftb mode=2 option=2 search=-0.02,0.04 |\\                                   ",
" sprelobc dtb=4 dxy=0.2 > stage1.su                                            ",
"                                                                               ",
"                                                                               ",
"                                                                               ",
NULL};

/* Credits:
 *	Sanyu Ye: READ Well Services, May. 2010
 */
/**************** end self doc *******************************************/

float calcRMS(int nt, float* data);

int verbose;		/* =0 no info, >0 output info  */

int main(int argc, char **argv) {
    cwp_String key, key_a, key_g;	/* header key word from segy.h		*/
    cwp_String type, type_a, type_g;    /* type of key				*/
    Value val, valnew, vala;                  /* value of key			*/
    int index, index_a, index_g;         /* index of key				*/
    int nt;       /* nums amps as int			*/
    int nxmax;                     /* number of traces			*/
    int nsegy;                  /* length bytes read for segy trace	*/
    int itr, ntr, ngather, ntotal;
    float t0, td0, td1, td2, dtd, dt, dT;
    float v, v0, dv, dvv, vmin, vmax, wV;
    double xy, dxy, zz, dzz, s, smin, smax, dx, dy, dz;
    double gx, gy, gz, gX, gY, gZ, xp, yp, zp, sz, x, y, z;
    double scale, scalco, spacing, d=0.0, l=0.0;
    double dti, dtmin, dtsum;
    segy tr, outtr;

    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose=0;

    /* get keys */
    if (!getparstring("ftp", &key_a))	 key_a="laga";
    type_a = hdtype(key_a);
    index_a = getindex(key_a);
    if (!getparstring("gdepth", &key_g)) key_g="gwdep";
    type_g = hdtype(key_g);
    index_g = getindex(key_g);
    if (!getparstring("key", &key))	 key="cdpt";
    type = hdtype(key);
    index = getindex(key);

    if (!getparint("nxmax", &nxmax))    nxmax = 1024;

    if (!getpardouble("scale", &scale)) scale = 0.0001;
    if (!getpardouble("xy",   &xy))        xy = 25.0;
    if (!getpardouble("dxy", &dxy))       dxy = 0.1;
    if (!getpardouble("z",    &zz))        zz = 0.0;
    if (!getpardouble("dz",   &dzz))      dzz = 1.0;
    if (!getpardouble("s",    &s))         s = 50.0;
    if (!getpardouble("smin", &smin))   smin = s - 2.0;
    if (!getpardouble("smax", &smax))   smax = s + 0.1;
    if (!getparfloat("dv",   &dvv))      dvv = 1.0;
    if (!getparfloat("v0",   &v0))        v0 = 1480.0;
    if (!getparfloat("vmin", &vmin))    vmin = v0;
    if (!getparfloat("vmax", &vmax))    vmax = v0;
    if (!getparfloat("td1",  &td1))      td1 = 20;
    if (!getparfloat("td2",  &td2))      td2 = 60;
    if (!getparfloat("t0",   &t0))        t0 = 1.0;

    /* get first trace */
    if ( (nsegy=gettr(&tr)) < HDRBYTES) err("can't get first trace");
    gethval(&tr, index, &val);

    dt = ((double) tr.dt) / 1000000.0;
    if (!getparfloat("dtd",  &dtd))      dtd = 1000*dt;
    nt = tr.ns;
    scalco = tr.scalco;
    if ( scalco == 0.0 ) scalco = 1.0;
    else if (scalco < 0.0) scalco = -1.0/scalco;
    
    segyhdr* hdrs = (segyhdr *) ealloc1(nxmax, HDRBYTES);
    float** indata = ealloc2float(nt, nxmax);
    double*  sx   = ealloc1double(nxmax);
    double*  sy   = ealloc1double(nxmax);
    float*  ftb  = ealloc1float(nxmax);
    float*  tred = ealloc1float(nxmax);
    
    //decomposing wavefield
    memset(hdrs, 0,  nxmax * HDRBYTES);
    memset(*indata, 0, nxmax * nt * FSIZE);
    memset(sx,  0, nxmax * DSIZE);
    memset(sy,  0, nxmax * DSIZE);
    memset(ftb, 0, nxmax * FSIZE);
    memset(tred,0, nxmax * FSIZE);

    /* Read headers and data while getting a count */
    int eof = 0;
    ngather = ntr = ntotal = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nxmax - 1) err("\nNumber of traces exceeding nxmax=%d\n", nxmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);
            memcpy(indata[ntr], tr.data, FSIZE * nt);
            gethval(&tr, index_a, &vala);
            ftb[ntr] = scale*vtod(type_a, vala);
            sx[ntr] = scalco*tr.sx;
            sy[ntr] = scalco*tr.sy;
            if (!ntr) { // first trace
                sz = scalco*tr.sdepth;
                gethval(&tr, index_g, &vala);
                gz = scalco*vtod(type_g, vala);
                gx = scalco*tr.gx;
                gy = scalco*tr.gy;
            }
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (verbose) warn("  processing %d traces %d-th gather (%s=%d)", ntr, ngather, key, vtoi(type, val));

            ntotal += ntr;
            
            // calculate reduction time
            for (itr = 0; itr < ntr; ++itr) {
                if (!t0) tred[itr] = 0.0;
                else tred[itr] = sqrt((sx[itr]-gx)*(sx[itr]-gx) + (sy[itr]-gy)*(sy[itr]-gy) + gz*gz)/v0;
            }

            // find minimum
            dtmin = 1E20;
            for (x = gx - xy; x <= gx + xy + 0.1*dxy; x += dxy) {
                for (y = gy - xy; y <= gy + xy + 0.1*dxy; y += dxy) {
                    for (z = gz - zz; z <= gz + zz + 0.1*dzz; z += dzz) {
                        spacing = (ngather == 1)? s : sqrt((x-xp)*(x-xp) + (y-yp)*(y-yp) + (z-zp)*(z-zp));
                        if ( spacing <= smin || spacing >= smax) continue;
                        for (v = vmin; v <= vmax + 0.1*dvv; v += dvv) {
                            for(td0 = td1; td0 <= td2 + 0.1*dtd; td0 += dtd) {
                                dtsum = 0.0;
                                for (itr = 0; itr < ntr; ++itr) {
                                    dti = ftb[itr] - t0 - 0.001*td0 + tred[itr] -
                                        sqrt((sx[itr]-x)*(sx[itr]-x) + (sy[itr]-y)*(sy[itr]-y) + (z-sz)*(z-sz))/v;
                                    dtsum += dti*dti;
                                }
                                if (dtsum < dtmin) {
                                    dtmin = dtsum;
                                    gX = x;
                                    gY = y;
                                    gZ = z;
                                    wV = v;
                                    dT = td0;
                                    //if (verbose>3) warn("  new position found for %d-th gather (%s=%d)", ngather, key, vtoi(type, val));
                                }
                            }
                        }
                    }
                }
            }

            dx = gX - gx;
            dy = gY - gy;
            dz = gZ - gz;
            dv = wV - v0;
            if (ngather>1) d = sqrt((gX-xp)*(gX-xp) + (gY-yp)*(gY-yp) + (gZ-zp)*(gZ-zp));
            if (ngather>1) l += d;

            // write out
            for (itr = 0; itr < ntr; ++itr) {
                memcpy(&outtr, &hdrs[itr], HDRBYTES);
                memcpy(outtr.data, indata[itr], FSIZE * nt);
                outtr.fx = gX;
                outtr.fy = gY;
                outtr.fz = gZ;
                outtr.dx = dx;
                outtr.dy = dy;
                outtr.dz = dz;
                outtr.wevel = wV;
                outtr.lagb = dT*10;
                outtr.ungpow = d;
                outtr.unscale = l;
                puttr(&outtr);
            }

            if(verbose) {
                warn("Position correction: dt=%1.0f dx=%3.1f dy=%3.1f dz=%3.1f z=%3.0f dv=%1.0f d=%4.1f  l=%4.1f",
                      dT, dx, dy, dz, gZ, dv, d, l);
            }
            
            val = valnew;
            // reset output data
            memset(hdrs, 0,  nxmax * HDRBYTES);
            memset(*indata, 0, nxmax * nt * FSIZE);
            memset(sx,  0, nxmax * DSIZE);
            memset(sy,  0, nxmax * DSIZE);
            memset(ftb, 0, nxmax * FSIZE);
            memset(tred,0, nxmax * FSIZE);
            ntr = 0;
            xp = gX; yp = gY; zp = gZ;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces of %d gathers are processed", ntotal, ngather);
    
    return(CWP_Exit());
}

