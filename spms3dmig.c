/* Copyright (c) RWS, Read Well Services, Oslo, 2012.*/
/* All rights reserved.                       */

/* SPMS3DMIG: $Revision: 1.0 $ ; $Date: 2012/09/19 18:32:28 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"
#include <fftw3.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPMS3DMIG - 3D/3C MIGration/location of Micro-Seismic events                  ",
"                                                                               ",
" spmig3dms zmin= zmax= xmin= xmax= ymin= ymax= <stdin >stdout  [optional parameters]",
"                                                                               ",
" Parameters:                                                                   ",
" key=ep                gather sorting key                                      ",
" nc=3                  number of components/panels within a gather             ",
" nlevel=12             number of VSP receiver levels/stations                  ",
"                                                                               ",
" dx=10           [m]   x sampling interval of image cube                       ",
" dy=dx           [m]   y sampling interval of image cube                       ",
" dz=dx           [m]   z sampling interval of image cube                       ",
"                                                                               ",
" zmax=           [m]   max. z of image cube                                    ",
" zmin=           [m]   min. z of image cube                                    ",
" xmax=           [m]   max. x of image cube                                    ",
" xmin=           [m]   min. x of image cube                                    ",
" ymax=           [m]   max. y of image cube                                    ",
" ymin=           [m]   min. y of image cube                                    ",
"                       above parameters, if not given, are set equal to subcube",
"                                                                               ",
" xwin=300        [m]   half length in x of subcube centered at injection point ",
" ywin=xwin       [m]   half length in y of subcube centered at injection point ",
" zwin=xwin       [m]   half length in z of subcube centered at injection point ",
"                       within which migration is performed                     ",
"                                                                               ",
" fmin=20.0      [Hz]   lower limit of frequency band to be filtered            ",
" fmax=240.0     [Hz]   upper limit of frequency band to be filtered            ",
" dfmax=10.0     [Hz]   upper limit of sampling interval in frequency domain    ",
"                                                                               ",
" vfp=1.0               velocity factor of P                                    ",
" vfsv=1.0              velocity factor of SV                                   ",
" vfsh=1.0              velocity factor of SH                                   ",
"                                                                               ",
" nxy=                  horizontal/lateral dimension of traveltime/ray-angle table",
" ttp=                  file of travel time table for P                         ",
" ttsv=                 file of travel time table for SV                        ",
" ttsh=                 file of travel time table for SH                        ",
" rap=                  file of ray angle table for P                           ",
" rasv=                 file of ray angle table for SV                          ",
" rash=                 file of ray angle table for SH                          ",
"                                                                               ",
" vp=3000       [m/s]   velocity of P                                  ",
" vsv=vp/2      [m/s]   velocity of SV                                 ",
" vsh=vsv       [m/s]   velocity of SH                                 ",
"                                                                               ",
" output=0              =-2 output event location with data for next iteration  ",
"                       =-1 output input data rotated toward new source position",
"                       =0 output only slices of image subcube through event location",
"                       =1 output only SH image subcube                         ",
"                       =2 output only SV image subcube                         ",
"                       =3 output SH and SV image subcube                       ",
"                       =4 output only SH image cube                            ",
"                       =5 output only SV image cube                            ",
"                       =6 output SH and SV image cube                          ",
"                                                                               ",
" verbose=0             >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
"                                                                               ",
" This is a general implementation of 3D/3C migration/location of microseismic  ",
" events with semblence weighted Wiener decon filter as proposed by Jakob Haldorsen",
" (Haldorsen et at., : Locating microseismic sources using migration-based ",
" deconvolution",
" The program computes object function within a cube centered at injection point",
" (source) and determines the location as point with max. semblance/cross correlation ",
" between P and SH/SV.",
" Input data must be sorted by event and component in order of ZEN. Depths of",
" receivers and source are stored in header keyword gelev and sdepth respectively.",
" The computational time is proportional to the size of search grid within image",
" subcube. It can be dramatically reduced by cascade with coarse and big, then  ",
" fine and smaller search grid.                                                 ",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
" 1. Locate and examine the result ",
"                                                                               ",
" spmig3dms < ZEN.su output=0 dx=10 xwin=300 |\\",
" suximage perc=99 windowtitle=\"X Y Z cross sections of object function through source location\" & |\\ ",
"                                                                               ",
" 2. Cascade to speed up for fine resolution",
"                                                                               ",
" spmig3dms output=-1 dx=20 xwin=300 < ZEN.su |\\",
" spmig3dms output=0  dx=5  xwin=40  |\\",
" suximage perc=99 windowtitle=\"X Y Z cross sections of object function through source location\" & |\\ ",
"                                                                               ",
" 3. Locate and examine thw whole image subcube ",
"                                                                               ",
" spmig3dms < ZEN.su output=3 dx=10 xwin=300 |\\",
" suxmovie perc=99 n2=61 dframe=10 fframe=-320 windowtitle=\"Cross section at X = %3.0f\" & |\\ ",
"                                                                               ",
" Version 0.9.0 last modified Sept. 2012 by Sanyu Ye                             ",
"                                                                               ",
 NULL};

/* Credits:
 *      RWS: Sanyu Ye, Sept. 2012
 *
 *
 * Notes:
 *
 */
/**************** end self doc *******************************************/

#include "sprinthlp.c"

// forward declare prototype
int  Conv2AmpPhase(int direction, int nwfft, int ifs, int ife, float** cdata);
int  DoFFT(int fft, fftwf_plan plan, int nwfft, int ntfft, float fftscale, float **cdata);
int mig3d(int nrcv, float** rxyz, float* xi, float* yi, float* zi, float* vf,
   int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int ifs, int ife, float df,
   float*** cdata, float*** ttp, float*** ttsv, float*** ttsh, float*** rap, float*** rasv, float*** rash, float*** impsv, float*** imgpsh);
int ms3dmig(int nrcv, int nt, float dt, float** rxyz, float* sxyz, float xwin, float ywin, float zwin,
    int nx, float* xi, int ny, float* yi, int nz, float* zi, float* vf, float fmin, float fmax,  \
   float*** trdata, float*** ttp, float*** ttsv, float*** ttsh, float*** rap, float*** rasv, float*** rash, float*** imgpsv, float*** imgpsh);

int main(int argc, char **argv)
{
    cwp_String key;       /* header key word from segy.h */
    cwp_String type;     /* type of key	*/
    int index;          /* index of key	*/
    Value val, valnew;    /* value of key				*/
    int nsegy;                  /* number of bytes read in the segy trace */
    int ntr;                /* number of actual input traces of the gather just read in*/
    int ngather, ntotal;    /* number of total input traces and gathers */
    int nt;                 /* number of points on trace		*/
    int i, j, ix, iy, iz, itr; /* counters				*/

    int verbose;
    int output;
    int nc, nrcv;
    int nx, ny, nz, nxy; /* max dimensions of input/output trace array */
    int ixpsv, iypsv, izpsv, ixpsh, iypsh, izpsh; /* index of microseismic event/max. cross-correlation */

    float dt;       /* time sample interval (sec)	*/
    float dfmax, fmin, fmax, fNyq;
    float dx, dy, dz, xmin, xmax, ymin, ymax, zmin, zmax;
    float xwin, ywin, zwin;
    float vf[3], v[3], sxyz[3];

    char *file;                 /* name of time table file	*/

    FILE *filefp = NULL;        /* file pointer for file	*/

    segy tr, tro;
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* Get parameters  */
    if(!getparint("verbose", &verbose))    verbose = 0;
    if(!getparint("nc", &nc)) nc = 3;
    if(!getparint("nlevel", &nrcv)) nrcv = 12;
    if(!getparint("output", &output)) output = 0;


    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    float scale = (float) tr.scalco;
    if (scale == 0.0) scale = 1.0;
    else if (scale <  0.0) scale = -1.0/scale;

    // get source (injection point) receiver coordinates
    sxyz[0] = scale*tr.sx;
    sxyz[1] = scale*tr.sy;
    sxyz[2] = scale*tr.sdepth;

    nt = tr.ns;
    dt = ((double) tr.dt) / 1000000.0;
    fNyq = 0.5/dt; // Nyquist frequency
    
    if(!getparfloat("dx",   &dx))     dx = 10.0;
    if(!getparfloat("dy",   &dy))     dy = dx;
    if(!getparfloat("dz",   &dz))     dz = dx;

    if(!getparfloat("xwin", &xwin)) xwin = 300;
    if(!getparfloat("ywin", &ywin)) ywin = xwin;
    if(!getparfloat("zwin", &zwin)) zwin = xwin;

    if(!getparfloat("zmin", &zmin)) {
        zmin = sxyz[2] - zwin;
    }
    if(!getparfloat("zmax", &zmax)) {
        zmax = sxyz[2] + zwin;
    }
    if(!getparfloat("xmax", &xmax)) {
        xmax= sxyz[0] + xwin;
    }
    if(!getparfloat("xmin", &xmin)) {
        xmin = sxyz[0] - xwin;
    }
    if(!getparfloat("ymax", &ymax)) {
        ymax = sxyz[1] + ywin;
    }
    if(!getparfloat("ymin", &ymin)) {
        ymin = sxyz[1] - ywin;
    }

    nx = NINT((xmax - xmin)/dx) + 1;
    ny = NINT((ymax - ymin)/dy) + 1;
    nz = NINT((zmax - zmin)/dz) + 1;

    if(!getparint("nxy", &nxy)) nxy = nx;

    if(!getparfloat("vfp",  &vf[0]))  vf[0] = 1.0;
    if(!getparfloat("vfsv", &vf[1]))  vf[1] = 1.0;
    if(!getparfloat("vfsh", &vf[2]))  vf[2] = 1.0;
    if(!getparfloat("vp",  &v[0]))  v[0] = 3000;
    if(!getparfloat("vsv", &v[1]))  v[1] = 0.5*v[0];
    if(!getparfloat("vsh", &v[2]))  v[2] = v[1];
    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "ep";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    if (!getparfloat("fmin",  &fmin)) fmin = 20.0;
    if (!getparfloat("fmax",  &fmax)) fmax = 240;
    if ( fmax > fNyq) fmax = fNyq;
    if (!getparfloat("dfmax",&dfmax)) dfmax = 10.0;



    // allocate memory input data traces
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nrcv, nc, HDRBYTES);
    float*** trdata = ealloc3float(nt, nrcv, nc);
    float*** imgpsv= (output < 0)? NULL : ealloc3float(nz, nx, ny);
    float*** imgpsh= (output < 0)? NULL : ealloc3float(nz, nx, ny);
    float**  rxyz  = ealloc2float(3, nrcv);
    float*   xi    = ealloc1float(nx);
    float*   yi    = ealloc1float(ny);
    float*   zi    = ealloc1float(nz);

    /* zero out data memory */
    memset(**trdata, 0, nc*nrcv*nt*FSIZE);
    memset(*hdrs2d,  0, nc*nrcv*HDRBYTES);
    if (imgpsv) memset(**imgpsv, 0, ny*nx*nz*FSIZE);
    if (imgpsh) memset(**imgpsh, 0, ny*nx*nz*FSIZE);

    float*** ttp  = ealloc3float(nz, nxy, nrcv);
    float*** ttsv = ealloc3float(nz, nxy, nrcv);
    float*** ttsh = ealloc3float(nz, nxy, nrcv);
    float*** rap  = ealloc3float(nz, nxy, nrcv);
    float*** rasv = ealloc3float(nz, nxy, nrcv);
    float*** rash = ealloc3float(nz, nxy, nrcv);

    memset(**ttp,   0, nrcv*nxy*nz*FSIZE);
    memset(**ttsv,  0, nrcv*nxy*nz*FSIZE);
    memset(**ttsh,  0, nrcv*nxy*nz*FSIZE);
    memset(**rap,   0, nrcv*nxy*nz*FSIZE);
    memset(**rasv,  0, nrcv*nxy*nz*FSIZE);
    memset(**rash,  0, nrcv*nxy*nz*FSIZE);

    for (ix=0; ix<nx; ++ix) xi[ix] = xmin + ix*dx;
    for (iy=0; iy<ny; ++iy) yi[iy] = ymin + iy*dy;
    for (iz=0; iz<nz; ++iz) zi[iz] = zmin + iz*dz;

    // figure out indices for image subcube
    int ixmin = MAX(0,  NINT((sxyz[0] - xwin - xmin)/dx));
    int ixmax = MIN(nx, NINT((sxyz[0] + xwin - xmin)/dx));
    int iymin = MAX(0,  NINT((sxyz[1] - ywin - ymin)/dy));
    int iymax = MIN(ny, NINT((sxyz[1] + ywin - ymin)/dy));
    int izmin = MAX(0,  NINT((sxyz[2] - zwin - zmin)/dz));
    int izmax = MIN(nz, NINT((sxyz[2] + zwin - zmin)/dz));
    
    if (output == 0) {
        if ( ixmax - ixmin != iymax - iymin || izmax - izmin != iymax - iymin )
        {
            err("unequal subcube sample numbers, please check");
        }
    }

    // read in traveltime table files
    size_t npdata = nrcv*nxy*nz;
    size_t npread = 0;
    if (getparstring("ttp", &file)) {
       filefp = efopen(file, "r");
       npread = fread((char *) &ttp[0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
       fclose(filefp);
       if (npread != npdata) err("TT file %s has wrong size read in (%d != %d)", npdata, npread);

        if (getparstring("ttsv", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &ttsv[0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("TT file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (getparstring("ttsh", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &ttsh[0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("TT file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (getparstring("rap", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &rap[0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("Ray angle file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (getparstring("rasv", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &rasv[0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("Ray angle file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (getparstring("rash", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &rash[0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("Ray angle file %s has wrong size read in (%d != %d)", npdata, npread);
        }
    }


    /* Read headers and data while getting a count */
    int eof = 0;
    int icmp  = 0;
    ngather = ntr = ntotal = 0;
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data */
            icmp = ntr/nrcv;
            ix = ntr%nrcv;  // receiver level index
            if (ntr > nc*nrcv - 1) err("  array dimension too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nc*nrcv, ntr+1, ngather+1, key, vtoi(type, val));
            if (ix > nx - 1) err("  array dimension nx too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nx, ix+1, ngather+1, key, vtoi(type, val));
            memcpy(&hdrs2d[icmp][ix], &tr, HDRBYTES);
            memcpy(trdata[icmp][ix], tr.data, FSIZE*tr.ns);

            if (!icmp) { // get source (injection point) receiver coordinates
                rxyz[ix][0] = scale*tr.gx;
                rxyz[ix][1] = scale*tr.gy;
                rxyz[ix][2] = scale*tr.gelev;
            }

            ++ntr;
        } else { // new gather or END_OF_FILE
            if (!ttp[0][0][0] && !ttp[nrcv -1][nxy - 1][nz -1]) {
                // populate traveltime table and ray angle table
                for (i=0; i<nrcv; ++i) {
                    for (iz=0; iz<nz; ++iz) {
                        for (ix=0; ix<nxy; ++ix) {
                            //float x = xi[ix] - rxyz[ircv][0];
                            //float y = yi[iy] - rxyz[i][1];
                            float z = zi[iz] - rxyz[i][2];
                            float x = xmin + ix*dx;
                            float d = sqrtf(z*z + x*x);
                            ttp[i][ix][iz]  = d / v[0];
                            ttsv[i][ix][iz] = d / v[1];
                            ttsh[i][ix][iz] = d / v[2];
                            float angle = atan2f(z, -x);  // -x in accordance with Jakob's angle
                            rap[i][ix][iz]  = angle;
                            rasv[i][ix][iz] = angle;
                            rash[i][ix][iz] = angle;
                        }
                    }
                }
            }

            ++ngather;
            ntotal += ntr;
            if (verbose) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, vtoi(type, val));

            itr = 0;

            ms3dmig(nrcv, nt, dt, rxyz, sxyz, xwin, ywin, zwin, nx, xi, ny, yi, nz, zi, vf, fmin, fmax,
               trdata, ttp, ttsv, ttsh, rap, rasv, rash, imgpsv, imgpsh);
            float maxpsv = 0.0, maxpsh = 0.0;
            for (iy=iymin; iy<=iymax; ++iy) {
                for (ix=ixmin; ix<=ixmax; ++ix) {
                    for (iz=izmin; iz<=izmax; ++iz) {
                        if (imgpsv[iy][ix][iz] > maxpsv) {
                            ixpsv = ix;
                            iypsv = iy;
                            izpsv = iz;
                            maxpsv= imgpsv[iy][ix][iz];
                        }
                        if (imgpsh[iy][ix][iz] > maxpsh) {
                            ixpsh = ix;
                            iypsh = iy;
                            izpsh = iz;
                            maxpsh= imgpsh[iy][ix][iz];
                        }
                    }
                }
            }

            // compute azimuth and inclination of event to first receiver
            float azim = (180.0/PI)*atan2(xi[ixpsh] - rxyz[0][0], yi[iypsh] - rxyz[0][1]);
            if (azim < 0.0) azim += 360.0;
            float incl = (180.0/PI)*atan2(sqrtf((xi[ixpsh] - rxyz[0][0])*(xi[ixpsh] - rxyz[0][0]) + (yi[iypsh] - rxyz[0][1])*(yi[iypsh] - rxyz[0][1])), zi[izpsh] - rxyz[0][2]);

            //  output
            memcpy(&tro, &hdrs2d[0][0], HDRBYTES); // copy back header to output trace
            int ns = (output > 3)? nz : izmax - izmin + 1;
            tro.ns = ns;
            //tro.sx = NINT(xi[ixpsh]/scale);  // event location X from SH
            //tro.sy = NINT(yi[iypsh]/scale);  // event location Y from SH
            //tro.sdepth = NINT(zi[izpsh]/scale); // event location Z from SV
            tro.fx = xi[ixpsh];  // event location from SH
            tro.dx = xi[ixpsv];  // event location from SV
            tro.fy = yi[iypsh];
            tro.dy = yi[iypsv];
            tro.fz = zi[izpsh];
            tro.dz = zi[izpsv];
            tro.otrav = NINT(10.0*azim);
            tro.ungpow  = azim;
            tro.unscale = incl;
            tro.d1 = dz;
            tro.d2 = dx;
            tro.f1 = (output > 3)? zi[0] : zi[izmin];
            tro.dt = NINT(1000.0*dz);
            if (output == -2) {  // output original input data with estimated source for next run
                for (icmp=0; icmp<nc; ++icmp) {
                    for (i=0; i<nrcv; ++i) {
                        memcpy(&tro, &hdrs2d[icmp][i], HDRBYTES);
                        memcpy(tro.data, trdata[icmp][i], tro.ns*FSIZE);
                        tro.sx = NINT(xi[ixpsh]/scale);  // event location X from SH
                        tro.sy = NINT(yi[iypsh]/scale);  // event location Y from SH
                        tro.sdepth = NINT(zi[izpsh]/scale); // event location Z from SV
                        puttr(&tro);
                    }
                }
            } else if (output == -1) {  // output original input data projected as P/SV/SH towards estimated source for QC
                for (i=0; i<nrcv; ++i) {
                    float x = xi[ixpsh] - rxyz[i][0];
                    float y = yi[iypsh] - rxyz[i][1];
                    float z = zi[izpsh] - rxyz[i][2];
                    float ah = (180.0/PI)*atan2(y, x);
                    float av = (180.0/PI)*atan2(sqrtf(x*x + y*y), z);
                    Rotate2C(0, nt-1, 0, ah, trdata[1][i], trdata[2][i], trdata[1][i], trdata[2][i]);
                    Rotate2C(0, nt-1, 0, av, trdata[0][i], trdata[1][i], trdata[0][i], trdata[1][i]);
                }
                for (icmp=0; icmp<nc; ++icmp) {
                    for (i=0; i<nrcv; ++i) {
                        memcpy(&tro, &hdrs2d[icmp][i], HDRBYTES);
                        memcpy(tro.data, trdata[icmp][i], tro.ns*FSIZE);
                        tro.fx = xi[ixpsh];  // event location from SH
                        tro.dx = xi[ixpsv];  // event location from SV
                        tro.fy = yi[iypsh];
                        tro.dy = yi[iypsv];
                        tro.fz = zi[izpsh];
                        tro.dz = zi[izpsv];
                        tro.otrav = NINT(10.0*azim);
                        tro.ungpow  = imgpsh[iypsh][ixpsh][izpsh];
                        tro.unscale = imgpsv[iypsv][ixpsv][izpsv];
                        puttr(&tro);
                    }
                }
            } else if (output == 0) { // output slices
                tro.ungpow  = imgpsh[iypsh][ixpsh][izpsh];
                tro.unscale = imgpsv[iypsv][ixpsv][izpsv];
                tro.cdp = 1;
                tro.nvs = iypsh;
                tro.cdpt = 1;
                for (i=ixmin; i<=ixmax; ++i) {  // slice along X/SH
                    tro.nhs = i;
                    tro.f2 = xi[i];
                    memcpy(tro.data, &imgpsh[iypsh][i][izmin], ns*FSIZE);
                    puttr(&tro);
                }
                tro.cdp = 2;
                tro.nvs = iypsv;
                tro.cdpt = 1;
                for (i=ixmin; i<=ixmax; ++i) {  // slice along X/SV
                    tro.nhs = i;
                    tro.f2 = xi[i];
                    memcpy(tro.data, &imgpsv[iypsv][i][izmin], ns*FSIZE);
                    puttr(&tro);
                }
                tro.cdp = 1;
                tro.nhs = ixpsh;
                tro.cdpt = 2;
                for (i=iymin; i<=iymax; ++i) {  // slice along Y/SH
                    tro.nvs = i;
                    tro.f2 = yi[i];
                    memcpy(tro.data, &imgpsh[i][ixpsh][izmin], ns*FSIZE);
                    puttr(&tro);
                }
                tro.cdp = 2;
                tro.nhs = ixpsv;
                tro.cdpt = 2;
                for (i=iymin; i<=iymax; ++i) {  // slice along Y/SV
                    tro.nvs = i;
                    tro.f2 = yi[i];
                    memcpy(tro.data, &imgpsv[i][ixpsv][izmin], ns*FSIZE);
                    puttr(&tro);
                }
                tro.dt = NINT(1000.0*dy);
                tro.f1 = yi[iymax];
                tro.d1 = -dy;
                tro.d2 = dx;
                tro.cdp = 1;
                tro.nvs = izpsh;
                tro.cdpt = 3;
                for (i=ixmin; i<=ixmax; ++i) {  // time slice along X/SH
                    tro.nhs = i;
                    tro.f2 = xi[i];
                    for(j=0; j<ns; ++j) tro.data[j] = imgpsh[iymax - j][i][izpsh];
                    puttr(&tro);
                }
                tro.cdp = 2;
                tro.nvs = izpsv;
                tro.cdpt = 3;
                for (i=ixmin; i<=ixmax; ++i) {  // time slice along X/SH
                    tro.nhs = i;
                    tro.f2 = xi[i];
                    for(j=0; j<ns; ++j) tro.data[j] = imgpsv[iymax - j][i][izpsv];
                    puttr(&tro);
                }
            } else {
                tro.ungpow  = imgpsh[iypsh][ixpsh][izpsh];
                tro.unscale = imgpsv[iypsv][ixpsv][izpsv];
                iz = (output > 3)? 0 : izmin;
                if (output != 2 && output != 5 ) {
                    tro.cdp = 1;
                    for (iy = (output > 3)? 0 : iymin; iy <= (output > 3)? ny-1 : iymax; ++iy) {
                        tro.nvs = iy;
                        for (ix = (output > 3)? 0 : ixmin; ix <= (output > 3)? nx-1 : ixmax; ++ix) {
                            tro.f2 = xi[ix];
                            tro.nhs = ix;
                            memcpy(tro.data, &imgpsh[iy][ix][iz], ns*FSIZE);
                            puttr(&tro);
                        }
                    }
                }
                if (output != 1 && output != 4 ) {
                    tro.cdp = 2;
                    for (iy = (output > 3)? 0 : iymin; iy <= (output > 3)? ny-1 : iymax; ++iy) {
                        tro.nvs = iy;
                        for (ix = (output > 3)? 0 : ixmin; ix <= (output > 3)? nx-1 : ixmax; ++ix) {
                            tro.f2 = xi[ix];
                            tro.nhs = ix;
                            memcpy(tro.data, &imgpsv[iy][ix][iz], ns*FSIZE);
                            puttr(&tro);
                        }
                    }
                }
            }
            
            /* zero out data memory */
            if (imgpsv) memset(**imgpsv, 0, nx*ny*nz*FSIZE);
            if (imgpsh) memset(**imgpsh, 0, nx*ny*nz*FSIZE);
            memset(**trdata, 0, nc*nrcv*nt*FSIZE);
            memset(*hdrs2d,  0, nc*nrcv*HDRBYTES);

            val = valnew;
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) warn(" Totally %d traces of %d gathers are processed", ntotal, ngather);

    return (CWP_Exit());
}

int Conv2AmpPhase(int direction, int nwfft, int ifs, int ife, float** cdata)
// only partially between frequency sample ifs ~ ife
{
    int iw, ifr;
    float tmp;

    for (iw=0; iw<nwfft; ++iw) {
        for (ifr=ifs; ifr<=ife; ++ifr) {
            if (direction < 0) { // from amp/phase to real/imaginary
                tmp = cdata[iw][2*ifr] * cosf(cdata[iw][2*ifr + 1]);
                cdata[iw][2*ifr + 1] = cdata[iw][2*ifr] * sinf(cdata[iw][2*ifr + 1]);
                cdata[iw][2*ifr] = tmp;
            } else {
                tmp = sqrtf(cdata[iw][2*ifr]*cdata[iw][2*ifr] + cdata[iw][2*ifr + 1]*cdata[iw][2*ifr + 1]);
                cdata[iw][2*ifr + 1] = atan2(cdata[iw][2*ifr + 1], cdata[iw][2*ifr]);
                cdata[iw][2*ifr] = tmp;
            }
        }
    }
    return ife - ifs + 1;
}

int DoFFT(int fft, fftwf_plan plan, int nwfft, int ntfft, float fftscale, float **cdata)
{
    int iw, it;

    fftwf_complex* ctr = (fftwf_complex*) &cdata[0][0];
    if (fft == FFTW_FORWARD) {
        fftwf_execute_dft_r2c(plan, &cdata[0][0], ctr);
        for (iw=0; iw<nwfft; ++iw)
            for (it=0; it<ntfft; ++it)
               cdata[iw][it] *=  fftscale;
    } else {
        fftwf_execute_dft_c2r(plan, ctr, &cdata[0][0]);
    }

    return fft;
}

int ms3dmig(int nrcv, int nt, float dt, float** rxyz, float* sxyz, float xwin, float ywin, float zwin,
   int nx, float* xi, int ny, float* yi, int nz, float* zi, float* vf, float fmin, float fmax,
   float*** trdata, float*** ttp, float*** ttsv, float*** ttsh, float*** rap, float*** rasv, float*** rash, float*** impsv, float*** imgpsh)
{
    int i, icmp, ix, nc=3;

    // figure out indices for image subcube
    float xmin = xi[0];
    float   dx = xi[1] - xi[0];
    float ymin = yi[0];
    float   dy = yi[1] - yi[0];
    float zmin = zi[0];
    float   dz = zi[1] - zi[0];
    int ixmin = MAX(0,  NINT((sxyz[0] - xwin - xmin)/dx));
    int ixmax = MIN(nx, NINT((sxyz[0] + xwin - xmin)/dx));
    int iymin = MAX(0,  NINT((sxyz[1] - ywin - ymin)/dy));
    int iymax = MIN(ny, NINT((sxyz[1] + ywin - ymin)/dy));
    int izmin = MAX(0,  NINT((sxyz[2] - zwin - zmin)/dz));
    int izmax = MIN(nz, NINT((sxyz[2] + zwin - zmin)/dz));

    float dfmax=10.0; // max. sampling interval in f doamin, corresponding 0.1 s trace length

    // figure out samling interval and number in frequncy domain
    float fNyq = 0.5/dt; // Nyquist frequency
    int ntfft = NINT(2.0 * ( fNyq / dfmax + 1));
    if (nt > ntfft) {
        ntfft = 2*(nt/2 + 1);
    }
    //int nf = ntfft/2;
    float df = fNyq/(ntfft/2 - 1);
    int ifs = NINT(fmin/df); // starting sample of filter window in f
    int ife = NINT(fmax/df); // end sample of filter window in f

    float*** cdata = ealloc3float(ntfft, nrcv, nc);
    /* zero out data memory */
    memset(**cdata,  0, nc*nrcv*ntfft*FSIZE);

    /* creat fftw plan*/
    fftwf_plan planf = NULL;
    fftwf_complex* ctr = (fftwf_complex*) &cdata[0][0][0];
    int ntt   = ntfft - 1;
    int nhalf = ntfft/2;
    planf = fftwf_plan_many_dft_r2c(1, &ntt, nrcv, &cdata[0][0][0], NULL, 1, ntfft, ctr, NULL, 1, nhalf, FFTW_MEASURE);
    float fftscale = 1.0/((float) (ntt));  // fft scaling factor

    for (icmp=0; icmp<nc; ++icmp) {
        for (ix=0; ix<nrcv; ++ix) {
            memcpy(cdata[icmp][ix], trdata[icmp][ix], nt*FSIZE);
        }
    }
    for (i=0; i<nc; ++i) DoFFT(FFTW_FORWARD, planf, nrcv, ntfft, fftscale, cdata[i]);

    mig3d(nrcv, rxyz, xi, yi, zi, vf, ixmin, ixmax, iymin, iymax, izmin, izmax, ifs, ife, df,
       cdata, ttp, ttsv, ttsh, rap, rasv, rash, impsv, imgpsh);

    return 0;
}

int mig3d(int nrcv, float** rxyz, float* xi, float* yi, float* zi, float* vf,
   int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax, int ifs, int ife, float df,
   float*** cdata, float*** ttp, float*** ttsv, float*** ttsh, float*** rap, float*** rasv, float*** rash, float*** imgpsv, float*** imgpsh)
{
    int ix, iy, iz, iw, i, j;
    float e[3], phaseshift[3], rsig[3], isig[3];

    float xmin = xi[0];
    float   dx = xi[1] - xi[0];

    float**  rpsvh = ealloc2float(nrcv, 3);
    float**  ipsvh = ealloc2float(nrcv, 3);
    float*   dist  = ealloc1float(nrcv);
    int*     ixy   = ealloc1int(nrcv);

    // do the computation
    for (iy=iymin; iy<=iymax; ++iy) {
        for (ix=ixmin; ix<=ixmax; ++ix) {
            for (iz=izmin; iz<=izmax; ++iz) {
                float psv = 0.0, psh = 0.0;
                for (iw=ifs; iw<=ife; ++iw) {
                    float dmin = FLT_MAX;
                    e[0] = e[1] = e[2] = 0.0; // reset energies
                    for (i=0; i<nrcv; ++i) {
                        float x = xi[ix] - rxyz[i][0];
                        float y = yi[iy] - rxyz[i][1];
                        float z = zi[iz] - rxyz[i][2];
                        float hdist = sqrtf(x*x + y*y);
                        float hcos = (hdist < FLT_MIN)? 1.0 : x/hdist;
                        float hsin = (hdist < FLT_MIN)? 0.0 : y/hdist;

                        // project 3C receivers to radial/transversal resp P/SV/SH
                        // first rotate two horizontal components
                        rpsvh[1][i] =  hcos*cdata[1][i][2*iw]     + hsin*cdata[2][i][2*iw]; //real, radial horizontal
                        ipsvh[1][i] =  hcos*cdata[1][i][2*iw + 1] + hsin*cdata[2][i][2*iw + 1]; //imaginary, radial horizontal
                        rpsvh[2][i] = -hsin*cdata[1][i][2*iw]     + hcos*cdata[2][i][2*iw]; //real, trans, SH
                        ipsvh[2][i] = -hsin*cdata[1][i][2*iw + 1] + hcos*cdata[2][i][2*iw + 1];; //imaginary, trans, SH

                        dist[i]  = sqrtf(hdist*hdist + z*z);
                        if (dist[i] < dmin) dmin = dist[i];
                        ixy[i] = NINT((hdist - xmin)/dx);
                        float vcosp = ABS(cosf(rap[i][ixy[i]][iz]));
                        float vsinp = sinf(rap[i][ixy[i]][iz]);
                        float vcosv = ABS(cosf(rasv[i][ixy[i]][iz]));
                        float vsinv = sinf(rasv[i][ixy[i]][iz]);

                        // then project vertical to radial P and transversal SV
                        rpsvh[0][i] =  vcosv*rpsvh[1][i]           + vsinp*cdata[0][i][2*iw]; //real, P
                        ipsvh[0][i] =  vcosv*ipsvh[1][i]           + vsinp*cdata[0][i][2*iw + 1]; //imaginary, P
                        rpsvh[1][i] = -vsinv*rpsvh[1][i]           + vcosp*cdata[0][i][2*iw]; //real, SV
                        ipsvh[1][i] = -vsinv*ipsvh[1][i]           + vcosp*cdata[0][i][2*iw + 1]; //imaginary, SV
                    }

                    // phase shift in f --- time shift in t
                    float t0;
                    for (i=0, t0=FLT_MAX; i<nrcv; ++i) if (t0 > ttp[i][ixy[i]][iz]) t0 = ttp[i][ixy[i]][iz];
                    for (j=0; j<3; ++j) {
                        for (i=0; i<nrcv; ++i) {
                            phaseshift[0] = 2.0*PI*df*iw*( ttp[i][ixy[i]][iz]  - t0)/vf[0];
                            phaseshift[1] = 2.0*PI*df*iw*(ttsv[i][ixy[i]][iz]  - t0)/vf[1];
                            phaseshift[2] = 2.0*PI*df*iw*(ttsh[i][ixy[i]][iz]  - t0)/vf[2];
                            float psin = sinf(phaseshift[j]);
                            float pcos = cosf(phaseshift[j]);
                            // phase shift plus spherical divergence correction
                            float rv = (pcos*rpsvh[j][i] - psin*ipsvh[j][i])*dist[i];
                            float iv = (psin*rpsvh[j][i] + pcos*ipsvh[j][i])*dist[i];
                            e[j] += rv*rv + iv*iv;
                            rpsvh[j][i] = rv;
                            ipsvh[j][i] = iv;
                        }
                        e[j] /= (float) nrcv;
                        rsig[j] = quick_select(rpsvh[j], nrcv);  // median real for signature real
                        isig[j] = quick_select(ipsvh[j], nrcv);
                    }
                    float rpsvxcor = (rsig[1]*rsig[0] + isig[1]*isig[0])/(e[0] + e[1] + e[2]);
                    //float ipsvxcor = (isig[1]*rsig[0] - rsig[1]*isig[0])/(e[0] + e[1] + e[2]);
                    //psv += rpsvxcor*rpsvxcor + ipsvxcor*ipsvxcor;
                    psv += rpsvxcor;
                    float rpshxcor = (rsig[2]*rsig[0] + isig[2]*isig[0])/(e[0] + e[1] + e[2]);
                    //float ipshxcor = (isig[2]*rsig[0] - rsig[2]*isig[0])/(e[0] + e[1] + e[2]);
                    //psh += rpshxcor*rpshxcor + ipshxcor*ipshxcor;
                    psh += rpshxcor;
                }
                psv = psv*psv;
                imgpsv[iy][ix][iz] = psv;
                psh = psh*psh;
                imgpsh[iy][ix][iz] = psh;
            }
        }
    }

    return 0;
}
