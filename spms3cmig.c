/* Copyright (c) RWS, Read Well Services, Oslo, 2012.*/
/* All rights reserved.                       */

/* SPMS3CMIG: $Revision: 1.0 $ ; $Date: 2012/09/19 18:32:28 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include "segyhdr.h"
#include <fftw3.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPMS3CMIG - fast 3D/3C MIGration/location of Micro-Seismic events             ",
"                                                                               ",
" spms3cmig <stdin >stdout  [optional parameters]                               ",
"                                                                               ",
" Parameters:                                                                   ",
" key=ep                gather sorting key                                      ",
" nlevel=12             number of VSP receiver levels/stations                  ",
" nc=4                  number of components expected in a gather               ",
" opt=0                 =0 use mean to extract signature                        ",
"                       =1 use madian to extract signature                      ",
"                                                                               ",
" dz=10           [m]   depth z sampling interval of image                      ",
" dr=10           [m]   radial (lateral distance) sampling interval of image    ",
" da=           [deg]   azimuth sampling interval of image                      ",
"                       if not given, calculated internally                     ",
"                                                                               ",
" zwin=300        [m]   half length in z of search window centered at source point",
" rwin=300        [m]   half length in x of search window                       ",
" awin=30       [deg]   half length in azimuth of search window                 ",
"                       standard deviation of azimuth by SPLOCEVT (key=ungpow)  ",
"                       is used if bigger than the given one                    ",
"                                                                               ",
" zmax=sdepth+300 [m]   max. z of traveltime table                              ",
" zmin=sdepth-300 [m]   min. z of traveltime table                              ",
" xmax=1500       [m]   max. r of traveltime table                              ",
" xmin=-xmax      [m]   min. r of traveltime table                              ",
"                                                                               ",
" fmin=20.0      [Hz]   lower limit of frequency band to be filtered            ",
" fmax=240.0     [Hz]   upper limit of frequency band to be filtered            ",
" dfmax=10.0     [Hz]   upper limit of sampling interval in frequency domain    ",
"                                                                               ",
" vp=           [m/s]   velocity of P,  default fetched from header wevel       ",
" vsv=          [m/s]   velocity of SV, default fetched from header swevel      ",
" vsh=vsv       [m/s]   velocity of SH                                          ",
"                       S velicities set to square root of 3 of Vp if not given ",
"                                                                               ",
//" nxy=                  horizontal/radial dimension of traveltime/ray-angle table",
" ttp=                  file of travel time table for P                         ",
" ttsv=                 file of travel time table for SV                        ",
" ttsh=                 file of travel time table for SH                        ",
" rap=                  file of ray angle table for P                           ",
" rasv=                 file of ray angle table for SV                          ",
" vpf=1.0               velocity factor of P                                    ",
" vsvf=1.0              velocity factor of SV                                   ",
" vshf=1.0              velocity factor of SH                                   ",
"                                                                               ",
" output=0              =-2 output event location with data for next iteration  ",
"                       =-1 output input data rotated toward new source position",
"                       =0 output only slices of image subcube through event location",
"                       =1 output SH and SV image subcube                       ",
"                                                                               ",
" fsemb=keyNo-semb.su   filen where slices on max. semblence are saved (for QC) ",
"                       applied only if output=-1                               ",
"                                                                               ",
" verbose=0             >0 output info                                          ",
"                                                                               ",
" Notes:                                                                        ",
"                                                                               ",
" This is a general implementation of 3D/3C migration/location of microseismic  ",
" events with semblence weighted Wiener decon filter as proposed by Jakob Haldorsen",
" (Haldorsen et at., : Locating microseismic sources using migration-based ",
" deconvolution",
" The program computes object function within a cube/sector (z, r, azimuth) centered",
" at source or injection point and determines the location as point with max.   ",
" semblance/cross correlation between P and SH/SV.",
" Input data must be sorted by event, receiver and component in order of ZEN.   ",
" Depths of receivers and source are fetched from header keyword gelev and sdepth  ",
" respectively. The computational time is proportional to the size of search grid",
" of image cube. To narrow the search grid, the source position estmated by program",
" SPLOCEVT is used. Here keywords fx/fy/fz give the lateral source position while",
" unscale/ungpow hold the values for central azimuth and standard deviation (half",
" search window). If the values of these extended header are zero, injection point",
" (sx/sy/sdepth) is used as search center point.                                ",
" Result of max. semblance between P and SH is used because it gives more correct",
" azimuth estimation. Keyword sdel holds lateral distance r, unscale azimuth,   ",
" ungpow the max. semblance, fx/fy/fz the estimated source location, dx error in",
" lateral distance r, dz error in z, dy error in distance along azimuthal direction",
" (error in azimuth x radius). The errors are estimated as half width of the peak",
" around the max. semblence at its half amplitude in respective directions. The ",
" error eistimated this way along azimuthal direction is large constantly (> 30 ",
" degrees). So the error estimated by SPLOCEVT is used if it is available and smaller.",
" If external traveltime and ray angle tables are input, the must have dimensions",
" exactly as specified by xmin/xmax/dx. zmin/zmax/dz, sorted by receiver, x and z.",
" They are binary floating point files and can be easily QC with ZTerra zeos if",
" a corresponding simple .H file is created.",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
" 1. A usual run in office without result from SPLOCEVT",
"    use cascade to speed up for fine resolution",
"                                                                               ",
" susort ep duse < ZEN.su |\\",
" segyclean |\\",
" sushw key=fx,fy,fz a=0 |\\",
" spms3cmig key=ep nc=3 output=-2 vp=4800 vsv=2700 rwin=600 awin=90 da=3 |\\",
" spms3cmig key=ep nc=3 output=-1 rwin=100 zwin=100 awin=6 da=0.5 |\\",
" tee ZEN+PSVSH.su |\\"
" spsort ep duse |\\",
" suxwigb perc=99 windowtitle=\"Input + projected  P, SV, SH after location\" &",
"                                                                               ",
" 2. Typical field sequence for real time monitoring",
"                                                                               ",
" mirf2su length=8.1 overlap=0.1 file=$recno.rcd hfile=myfile.hdr |\\",
" spsort ep |\\",
" subfilt fstoplo=10 fpasslo=20 fpasshi=240 fstophi=360 |\\",
" sushw key=styp a=$StageNo |\\",
" sushw match=styp key=sx,sy,sdepth infile=injpoint.xyz |\\",
" sushw match=cdp key=gx,gy,gelev infile=receiver.xyz |\\",
" sushw match=cdp key=otrav infile=receiver-azimuth.tbl |\\",
" sp2crot key=cdp a=otrav mode=0 vsp=-1 nc=3 pol=-1,1,-1 |\\",
" spevent key=ep level=2 nc=3 verbose=1 |\\",
" spsort ep |\\",
" spevent key=ep level=3 nc=4 verbose=1 |\\",
" tee $recno.su |\\",
" suwind key=corr min=1 |\\",
" sushw key=wevel,swevel a=4800,2700 |\\",
" splocevt key=ep nc=4 level=3 ps=1 |\\",
" suwind key=corr min=1 |\\",
" spms3cmig key=ep nc=7 output=-1 |\\",
" tee ZEN+PSVSH.su |\\"
" spsort ep duse |\\",
" suxwigb perc=99 windowtitle=\"Input + projected P, SV, SH after double location\"   ",
"                                                                               ",
"                                                                               ",
" Version 1.0.0 last modified Dec. 2012 by Sanyu Ye                             ",
"                                                                               ",
 NULL};

/* Credits:
 *      RWS: Sanyu Ye, Nov. 2012
 *
 *
 * Notes:
 *
 */
/**************** end self doc *******************************************/

#include "sprinthlp.c"

// forward declare prototype
int  DoFFT(int fft, fftwf_plan plan, int nwfft, int ntfft, float fftscale, float **cdata);

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
    int ifs, ife, iw;      /* index of frequency window in samples	*/
    int ntfft;      /* length of time window in samples	*/
    int nf;      /* length of time window in samples	*/

    int verbose;
    int output, opt;
    int nrcv;
    int nxy, nzz; /* dimensions of traveltime table (two-sided)*/
    int nx, ny, nz; /* max dimensions of output trace array */
    int ixpsv, iypsv, izpsv, ixpsh, iypsh, izpsh, ixshv, iyshv, izshv; /* index of microseismic event/max. cross-correlation */

    float dt;       /* time sample interval (sec)	*/
    float df, dfmax, fmin, fmax, fNyq;
    float dx, da, dz, xmin, xmax, amin, amax, zmin, zmax;
    float xwin, awin, zwin;
    float a, azimc, astdd, azpsv, azpsh, azshv, rinj, zinj, xinj, yinj;
    float v[3], vf[3], e[3], rsig[3], isig[3], sxyz[3];
    float psin[3], pcos[3];

    char *file;                 /* name of time table file	*/
    FILE *filefp = NULL;        /* file pointer for file	*/

    char *sfile;                /* name of file	for semblance slices */
    FILE *sfp = NULL;           /* file pointer for semblance slices	*/

    segy tr, tro;

    int nc = 4;   // number of components

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* Get parameters  */
    if(!getparint("verbose", &verbose))    verbose = 0;
    if(!getparint("nlevel", &nrcv)) nrcv = 12;
    if(!getparint("nc",     &nc))   nc   = 4;
    if(!getparint("output", &output)) output = 0;
    if(!getparint("opt",    &opt))       opt = 0;


    /* read first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES ) err("Cannot get first trace");

    float scale = (float) tr.scalco;
    if (scale == 0.0) scale = 1.0;
    else if (scale <  0.0) scale = -1.0/scale;


    nt = tr.ns;
    dt = ((double) tr.dt) / 1000000.0;
    fNyq = 0.5/dt; // Nyquist frequency

    if(!getparfloat("dr",   &dx))     dx = 10.0;
    if(!getparfloat("dz",   &dz))     dz = 10.0;

    if(!getparfloat("rwin", &xwin)) xwin = 300.0;
    if(!getparfloat("zwin", &zwin)) zwin = 300.0;

    // get injection point coordinates
    xinj = scale*tr.sy;
    yinj = scale*tr.sx;
    zinj = NINT(scale*tr.sdepth/dz) * dz;
    azimc = atan2(yinj, xinj)*180.0/PI;
    rinj = sqrtf(xinj*xinj + yinj*yinj);
    
    if(!getparfloat("zmin", &zmin)) {
        zmin = zinj - 300;
    }
    if(!getparfloat("zmax", &zmax)) {
        zmax = zinj + 300;
    }
    if(!getparfloat("xmax", &xmax)) xmax = 1500.0;
    if(!getparfloat("xmin", &xmin)) xmin = -xmax;

    nxy = NINT((xmax - xmin)/dx) + 1;
    nzz = NINT((zmax - zmin)/dz) + 1;
    
    nx = 2*NINT(xwin/dx) + 1;
    nz = 2*NINT(zwin/dz) + 1;
    ny = nz;
    
    if(!getparfloat("da",     &da))    da = 0.0;
    if(!getparfloat("awin",   &awin))  awin = 30.0;
    if ( awin != 0.0 && da != 0.0) {    
        ny = 2*NINT(awin/da) + 1;
    }

    if(!getparfloat("vpf",  &vf[0]))  vf[0] = 1.0;
    if(!getparfloat("vsvf", &vf[1]))  vf[1] = 1.0;
    if(!getparfloat("vshf", &vf[2]))  vf[2] = 1.0;
    if(!getparfloat("vp",  &v[0])) {
        v[0] = (float) tr.wevel;
    }
    if (v[0] < 1500.0) err("P velocity must be properly given (Vp=%1.0f [m/s])", v[0]);
    if(!getparfloat("vsv", &v[1]))  {
        v[1] = (float) tr.swevel;
    }
    if ( v[1] < 1000.0 ) {
        warn("S velocity set to default from the given one Vs=%1.0f to new =%1.0f [m/s])", v[1], 0.58*v[0]);
        v[1] = 0.58*v[0];
    }
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


    // figure out samling interval and number in frequncy domain
    ntfft = NINT(2.0 * ( fNyq / dfmax + 1));
    if (nt > ntfft) {
        ntfft = 2*(nt/2 + 1);
    }
    nf = ntfft/2;
    df = fNyq/(ntfft/2 - 1);
    ifs = NINT(fmin/df); // starting sample of filter window in f
    ife = NINT(fmax/df); // end sample of filter window in f
    int nff = ife - ifs + 1;

    // allocate memory input data traces
    segyhdr** hdrs2d = (segyhdr**) ealloc2(nrcv, nc, HDRBYTES);
    float*** trdata = ealloc3float(nt, nrcv, nc);
    float**  rpsvh = ealloc2float(nrcv, 3);
    float**  ipsvh = ealloc2float(nrcv, 3);
    float*** cdata = ealloc3float(ntfft, nrcv, 3);
    float*** imgpsv= ealloc3float(nz, nx, ny);
    float*** imgpsh= ealloc3float(nz, nx, ny);
    float*** imgshv= ealloc3float(nz, nx, ny);
    float**  rxyz  = ealloc2float(3, nrcv);
    float*   dist  = ealloc1float(nrcv);
    float*   xi    = ealloc1float(nxy);
    float*   zi    = ealloc1float(nzz);
    float** prc  = ealloc2float(nff, nrcv);
    float** prs  = ealloc2float(nff, nrcv);
    float** pra  = ealloc2float(nff, nrcv);
    float** pic  = ealloc2float(nff, nrcv);
    float** pis  = ealloc2float(nff, nrcv);
    float** pia  = ealloc2float(nff, nrcv);
    float** svrc = ealloc2float(nff, nrcv);
    float** svrs = ealloc2float(nff, nrcv);
    float** svra = ealloc2float(nff, nrcv);
    float** svic = ealloc2float(nff, nrcv);
    float** svis = ealloc2float(nff, nrcv);
    float** svia = ealloc2float(nff, nrcv);
    float** shrc = ealloc2float(nff, nrcv);
    float** shrs = ealloc2float(nff, nrcv);
    float** shic = ealloc2float(nff, nrcv);
    float** shis = ealloc2float(nff, nrcv);

    /* zero out data memory */
    memset(**cdata,  0, 3*nrcv*ntfft*FSIZE);
    memset(**trdata, 0, nc*nrcv*nt*FSIZE);
    memset(*hdrs2d,  0, nc*nrcv*HDRBYTES);
    memset(**imgpsv, 0, ny*nx*nz*FSIZE);
    memset(**imgpsh, 0, ny*nx*nz*FSIZE);
    memset(**imgshv, 0, ny*nx*nz*FSIZE);
    memset(*prc,  0, nrcv*nff*FSIZE);
    memset(*prs,  0, nrcv*nff*FSIZE);
    memset(*pra,  0, nrcv*nff*FSIZE);
    memset(*pic,  0, nrcv*nff*FSIZE);
    memset(*pis,  0, nrcv*nff*FSIZE);
    memset(*pia,  0, nrcv*nff*FSIZE);
    memset(*svrc, 0, nrcv*nff*FSIZE);
    memset(*svrs, 0, nrcv*nff*FSIZE);
    memset(*svra, 0, nrcv*nff*FSIZE);
    memset(*svic, 0, nrcv*nff*FSIZE);
    memset(*svis, 0, nrcv*nff*FSIZE);
    memset(*svia, 0, nrcv*nff*FSIZE);
    memset(*shrc, 0, nrcv*nff*FSIZE);
    memset(*shrs, 0, nrcv*nff*FSIZE);
    memset(*shic, 0, nrcv*nff*FSIZE);
    memset(*shis, 0, nrcv*nff*FSIZE);

    float**** tt   = ealloc4float(nzz, nxy, nrcv, 3);
    float**** ra   = ealloc4float(nzz, nxy, nrcv, 3);

    memset(***tt,    0, 3*nrcv*nxy*nzz*FSIZE);
    memset(***ra,    0, 3*nrcv*nxy*nzz*FSIZE);

    for (ix=0; ix<nxy; ++ix) xi[ix] = xmin + ix*dx;
    for (iz=0; iz<nzz; ++iz) zi[iz] = zmin + iz*dz;

    int filegiven = 0;
    if (output == -1) {
        if (getparstring("fsemb", &sfile)) {
            sfp = efopen(sfile, "a");
            filegiven = 1;
        }
    }

    // read in traveltime table files
    size_t npdata = nrcv*nxy*nzz;
    size_t npread = 0;
    if (getparstring("ttp", &file)) {
        float* fv;
        filefp = efopen(file, "r");
        npread = fread((char *) &tt[0][0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
        fclose(filefp);
        if (npread != npdata) err("TT file %s has wrong size read in (%d != %d)", npdata, npread);
        if (vf[0] != 1.0) for(i=0, fv=&tt[0][0][0][0]; i<npdata; ++i, ++fv) *fv *= vf[0];
        
        if (getparstring("ttsv", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &tt[1][0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("TT file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (vf[1] != 1.0) for(i=0, fv=&tt[1][0][0][0]; i<npdata; ++i, ++fv) *fv *= vf[1];
        
        if (getparstring("ttsh", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &tt[2][0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("TT file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (vf[2] != 1.0) for(i=0, fv=&tt[2][0][0][0]; i<npdata; ++i, ++fv) *fv *= vf[2];
        if (getparstring("rap", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &ra[0][0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("Ray angle file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (getparstring("rasv", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &ra[1][0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("Ray angle file %s has wrong size read in (%d != %d)", npdata, npread);
        }
        if (getparstring("rash", &file)) {
           filefp = efopen(file, "r");
           npread = fread((char *) &ra[2][0][0][0], FSIZE, npdata + 1, filefp);  // read over the end to detect error
           fclose(filefp);
           if (npread != npdata) err("Ray angle file %s has wrong size read in (%d != %d)", npdata, npread);
        }
    }

    /* creat fftw plan*/
    fftwf_plan planf = NULL;
    //fftwf_plan planb = NULL;
    fftwf_complex* ctr = (fftwf_complex*) &cdata[0][0][0];
    int ntt   = nt; //ntfft - 1;
    int nhalf = ntfft/2;
    //int nreal = ntfft*nrcv;
    //int ncmpl = nreal/2;
    planf = fftwf_plan_many_dft_r2c(1, &ntt, nrcv, &cdata[0][0][0], &ntfft, 1, ntfft, ctr, &nhalf, 1, nhalf, FFTW_MEASURE);
    //planb = fftwf_plan_many_dft_c2r(1, &ntt, nrcv, ctr, &ncmpl, 1, nhalf, &cdata[0][0][0], &nreal, 1, ntfft, FFTW_MEASURE);
    float fftscale = 1.0/((float) (ntt));  // fft scaling factor

    /* Read headers and data while getting a count */
    int eof = 0;
    int icmp  = 0;
    ngather = ntr = ntotal = 0;
    float xc = 0.0, yc = 0.0, zc = 0.0;  // center coordinates of receiver array
    do { /* Main loop over traces/gather */
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data */
            icmp = ntr%nc;
            ix = ntr/nc;  // receiver level index
            if (ntr > nc*nrcv - 1) err("  array dimension too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nc*nrcv, ntr+1, ngather+1, key, vtoi(type, val));
            if (ix > nx - 1) err("  array dimension nx too small (%d < %d) traces input for %d-th gather (%s=%d)",
                nx, ix+1, ngather+1, key, vtoi(type, val));
            memcpy(&hdrs2d[icmp][ix], &tr, HDRBYTES);
            memcpy(trdata[icmp][ix], tr.data, FSIZE*tr.ns);

            if (!icmp) { // get receiver coordinates
                rxyz[ix][0] = scale*tr.gy;
                rxyz[ix][1] = scale*tr.gx;
                rxyz[ix][2] = scale*tr.gelev;
            }

            if ( ntr == 0 ) {
                // get source coordinates as estimated by splocevt
                sxyz[0] = tr.fy;
                sxyz[1] = tr.fx;
                sxyz[2] = tr.fz;
                // azimuth and its standard deviation
                if ( tr.unscale != 0.0 ) azimc = tr.unscale;
                astdd = tr.ungpow;
            }
            ++ntr;
        } else { // new gather or END_OF_FILE
            if (!tt[0][0][0][0] && !tt[0][nrcv -1][nxy - 1][nzz -1]) {
                // populate traveltime table and ray angle table
                for (i=0; i<nrcv; ++i) {
                    for (iz=0; iz<nzz; ++iz) {
                        for (ix=0; ix<nxy; ++ix) {
                            //float x = xi[ix] - rxyz[ircv][0];
                            //float y = yi[iy] - rxyz[i][1];
                            float z = zi[iz] - rxyz[i][2];
                            float x = xmin + ix*dx;
                            float d = sqrtf(z*z + x*x);
                            tt[0][i][ix][iz] = d / v[0];
                            tt[1][i][ix][iz] = d / v[1];
                            tt[2][i][ix][iz] = d / v[2];
                            float angle = atan2f(z, ABS(x));  // -x in accordance with Jakob's angle
                            ra[0][i][ix][iz] = angle;
                            ra[1][i][ix][iz] = angle;
                            ra[2][i][ix][iz] = angle;
                        }
                    }
                }
            }

            // calc center coordinates of receiver array
            if (!xc && !yc && !zc) {  // if not set
                for (i=0; i<nrcv; ++i) {
                    xc += rxyz[i][0];
                    yc += rxyz[i][1];
                    zc += rxyz[i][2];
                }
                xc /= nrcv;
                yc /= nrcv;
                zc /= nrcv;
            }
            
            ++ngather;
            ntotal += ntr;

            int keyno = vtoi(type, val);
            if (verbose) warn("  Processing %d traces of %d-th gather (%s=%d)...", ntr, ngather, key, keyno);

            // figure out indices for image subcube
            float z0 = dz*NINT(sxyz[2]/dz);
            if (z0 < dz) z0 = zinj;
            int izmin = MAX(0,  NINT((z0 - zwin - zmin)/dz));
            int izmax = MIN(nzz-1, NINT((z0 + zwin - zmin)/dz));

            float r0 = sqrtf(sxyz[0]*sxyz[0] + sxyz[1]*sxyz[1]);
            if ( r0 < dx ) r0 = rinj;
            int ixmin = MAX(0,  NINT((MAX(0.0, r0 - xwin) - xmin)/dx));
            int ixmax = MIN(nxy-1, NINT((r0 + xwin - xmin)/dx));

            if ( awin < astdd ) {
                awin = astdd;
                da = 2.0 * awin / (ny - 1);
            }
            if (da == 0.0) {
                da = 2.0 * awin / (ny - 1);
            }
            amin = azimc - awin;
            amax = azimc + awin;
            if (verbose > 1) {
                warn("Central Azimuth=%3.1f  Half Windows=%3.1f  Step=%3.1f [deg]", azimc, awin, da);
            }


            itr = 0;
            for (icmp=0; icmp<3; ++icmp) {
                for (ix=0; ix<nrcv; ++ix) {
                    memcpy(cdata[icmp][ix], trdata[icmp][ix], nt*FSIZE);
                }
            }
            for (i=0; i<3; ++i) DoFFT(FFTW_FORWARD, planf, nrcv, ntfft, fftscale, cdata[i]);

            // do the computation
            //ixpsv = iypsv = izpsv = ixpsh = iypsh = izpsh = (ixmax - ixmin)/2;
            float maxpsv = 0.0, maxpsh = 0.0, maxshv = 0.0;
            for (ix=ixmin; ix<=ixmax; ++ix) {
                for (iz=izmin; iz<=izmax; ++iz) {
                    
                    float t0 = FLT_MAX;  // minimum tracel time of P
                    for (i=0; i<nrcv; ++i) if (t0 > tt[0][i][ix][iz]) t0 = tt[0][i][ix][iz];

                    for (i=0; i<nrcv; ++i) {
                        float x = xi[ix];
                        float z = zi[iz] - rxyz[i][2];
                        //float hdist = ABS(x);

                        dist[i]  = sqrtf(x*x + z*z);

                        //ixy[i] = NINT((hdist - xmin)/dx);
                        float vcosp = ABS(cosf(ra[0][i][ix][iz]));
                        float vsinp = sinf(ra[0][i][ix][iz]);
                        float vcosv = ABS(cosf(ra[1][i][ix][iz]));
                        float vsinv = sinf(ra[1][i][ix][iz]);

                        for (iw=ifs; iw<=ife; ++iw) {
                            // phase shift in f --- time shift in t
                            for (j=0; j<3; ++j) {
                                float phaseshift = 2.0*PI*df*iw*(tt[j][i][ix][iz] - t0);
                                psin[j] = sinf(phaseshift);
                                pcos[j] = cosf(phaseshift);
                            }

                            int ifrq = iw - ifs;
                            prc[i][ifrq] =   vcosp*(cdata[2][i][2*iw]*pcos[0] - cdata[2][i][2*iw+1]*psin[0]);
                            prs[i][ifrq] =   vcosp*(cdata[1][i][2*iw]*pcos[0] - cdata[1][i][2*iw+1]*psin[0]);
                            pra[i][ifrq] =   vsinp*(cdata[0][i][2*iw]*pcos[0] - cdata[0][i][2*iw+1]*psin[0]);
                            svrc[i][ifrq] = -vsinv*(cdata[2][i][2*iw]*pcos[1] - cdata[2][i][2*iw+1]*psin[1]);
                            svrs[i][ifrq] = -vsinv*(cdata[1][i][2*iw]*pcos[1] - cdata[1][i][2*iw+1]*psin[1]);
                            svra[i][ifrq] =  vcosv*(cdata[0][i][2*iw]*pcos[1] - cdata[0][i][2*iw+1]*psin[1]);
                            shrc[i][ifrq] =         cdata[1][i][2*iw]*pcos[2] - cdata[1][i][2*iw+1]*psin[2];
                            shrs[i][ifrq] =       -(cdata[2][i][2*iw]*pcos[2] - cdata[2][i][2*iw+1]*psin[2]);
                            pic[i][ifrq] =   vcosp*(cdata[2][i][2*iw]*psin[0] + cdata[2][i][2*iw+1]*pcos[0]);
                            pis[i][ifrq] =   vcosp*(cdata[1][i][2*iw]*psin[0] + cdata[1][i][2*iw+1]*pcos[0]);
                            pia[i][ifrq] =   vsinp*(cdata[0][i][2*iw]*psin[0] + cdata[0][i][2*iw+1]*pcos[0]);
                            svic[i][ifrq] = -vsinv*(cdata[2][i][2*iw]*psin[1] + cdata[2][i][2*iw+1]*pcos[1]);
                            svis[i][ifrq] = -vsinv*(cdata[1][i][2*iw]*psin[1] + cdata[1][i][2*iw+1]*pcos[1]);
                            svia[i][ifrq] =  vcosv*(cdata[0][i][2*iw]*psin[1] + cdata[0][i][2*iw+1]*pcos[1]);
                            shic[i][ifrq] =         cdata[1][i][2*iw]*psin[2] + cdata[1][i][2*iw+1]*pcos[2];
                            shis[i][ifrq] =       -(cdata[2][i][2*iw]*psin[2] + cdata[2][i][2*iw+1]*pcos[2]);
                        }
                    }
                                    
                    for (iy=0, a = amin; iy < ny; ++iy, a += da) {
                        float hcos = cosf(a*PI/180.0);
                        float hsin = sinf(a*PI/180.0);
                        float psv = 0.0, psh = 0.0, shv = 0.0;
                        for (iw=0; iw < nff; ++iw) {
                            e[0] = e[1] = e[2] = 0.0; // reset energies
                            for (i=0; i<nrcv; ++i) {
                                rpsvh[0][i] = dist[i]*(hcos*prc[i][iw] + hsin*prs[i][iw] + pra[i][iw]); //real, P
                                ipsvh[0][i] = dist[i]*(hcos*pic[i][iw] + hsin*pis[i][iw] + pia[i][iw]); //imaginary, P
                                rpsvh[1][i] = dist[i]*(hcos*svrc[i][iw] + hsin*svrs[i][iw] + svra[i][iw]);//real, SV
                                ipsvh[1][i] = dist[i]*(hcos*svic[i][iw] + hsin*svis[i][iw] + svia[i][iw]); //imaginary, SV
                                rpsvh[2][i] = dist[i]*(hcos*shrc[i][iw] + hsin*shrs[i][iw]);//real, SV
                                ipsvh[2][i] = dist[i]*(hcos*shic[i][iw] + hsin*shis[i][iw]); //imaginary, SV

                                e[0] += rpsvh[0][i]*rpsvh[0][i] + ipsvh[0][i]*ipsvh[0][i];
                                e[1] += rpsvh[1][i]*rpsvh[1][i] + ipsvh[1][i]*ipsvh[1][i];
                                e[2] += rpsvh[2][i]*rpsvh[2][i] + ipsvh[2][i]*ipsvh[2][i];
                            }

                            for (j=0; j < 3; ++j) {
                                rsig[j] = isig[j] = 0.0;
                                if (opt == 1) {
                                    rsig[j] = quick_select(rpsvh[j], nrcv);  // median real for signature real
                                    isig[j] = quick_select(ipsvh[j], nrcv);
                                } else {
                                    for (i=0; i<nrcv; ++i) {
                                        rsig[j] += rpsvh[j][i];
                                        isig[j] += ipsvh[j][i];
                                    }
                                    rsig[j] /= (float) nrcv;
                                    isig[j] /= (float) nrcv;
                                }
                            }

                            float etotal = (e[0] + e[1] + e[2])/(float) nrcv;
                            float rpsvxcor = (rsig[1]*rsig[0] + isig[1]*isig[0])/etotal;
                            psv += rpsvxcor;
                            float rpshxcor = (rsig[2]*rsig[0] + isig[2]*isig[0])/etotal;
                            psh += rpshxcor;
                        }
                        shv = ABS(psv)*ABS(psh);
                        if (shv > maxshv) {
                            ixshv = ix;
                            iyshv = iy;
                            izshv = iz;
                            azshv = a;
                            maxshv= shv;
                        }
                        imgshv[iy][ix - ixmin][iz - izmin] = shv;

                        psv = psv*psv;
                        if (psv > maxpsv) {
                            ixpsv = ix;
                            iypsv = iy;
                            izpsv = iz;
                            azpsv = a;
                            maxpsv= psv;
                        }
                        imgpsv[iy][ix - ixmin][iz - izmin] = psv;

                        psh = psh*psh;
                        if (psh > maxpsh) {
                            ixpsh = ix;
                            iypsh = iy;
                            izpsh = iz;
                            azpsh = a;
                            maxpsh= psh;
                        }
                        imgpsh[iy][ix - ixmin][iz - izmin] = psh;

                    }  // end of loop over azimuth
                }
            }

            // compute source coordinates x/y
            if (azpsh <   0.0) azpsh += 360.0;
            if (azpsh > 360.0) azpsh -= 360.0;
            if (azpsv <   0.0) azpsv += 360.0;
            if (azpsv > 360.0) azpsv -= 360.0;
            if (azshv <   0.0) azshv += 360.0;
            if (azshv > 360.0) azshv -= 360.0;

            float sy = yc + xi[ixpsh]*sinf(azpsh*PI/180.0);  // event location from SH
            float sx = xc + xi[ixpsh]*cosf(azpsh*PI/180.0);
            float sz = zi[izpsh];
            //float sy1 = yc + xi[ixpsv]*sinf(azpsv*PI/180.0);  // event location from SV
            //float sx1 = xc + xi[ixpsv]*cosf(azpsv*PI/180.0);
            //float sz1 = zi[izpsv];
            //float sy2 = yc + xi[ixshv]*sinf(azshv*PI/180.0);  // event location from SV
            //float sx2 = xc + xi[ixshv]*cosf(azshv*PI/180.0);
            //float sz2 = zi[izshv];

            if(verbose) {
                warn("Respective position by PSH/PSV/PSHV: z=%1.0f/%1.0f/%1.0f  r=%1.0f/%1.0f/%1.0f  azim=%3.1f/%3.1f/%3.1f",
                      zi[izpsh], zi[izpsv], zi[izshv], xi[ixpsh], xi[ixpsv], xi[ixshv], azpsh, azpsv, azshv);
            }

            // try to figure out error (distance to half peak value)
            float threshold = 0.5 * imgpsh[iypsh][ixpsh - ixmin][izpsh - izmin];
            int ndx=0, ndz=0, nda=0;
            // azimuthal direction
            for (iy = iypsh - 1; iy >= 0; --iy, ++nda) {
                if (imgpsh[iy][ixpsh - ixmin][izpsh - izmin] < threshold ) break;
            }
            for (iy = iypsh + 1; iy < ny; ++iy, ++nda) {
                if (imgpsh[iy][ixpsh - ixmin][izpsh - izmin] < threshold ) break;
            }
            for (ix = ixpsh - 1; ix >= ixmin; --ix, ++ndx) {
                if (imgpsh[iypsh][ix - ixmin][izpsh - izmin] < threshold ) break;
            }
            for (ix = ixpsh + 1; ix < ixmax; ++ix, ++ndx) {
                if (imgpsh[iypsh][ix - ixmin][izpsh - izmin] < threshold ) break;
            }
            for (iz = izpsh - 1; iz >= izmin; --iz, ++ndz) {
                if (imgpsh[iypsh][ixpsh - ixmin][iz - izmin] < threshold ) break;
            }
            for (iz = izpsh + 1; iz < izmax; ++iz, ++ndz) {
                if (imgpsh[iypsh][ixpsh - ixmin][iz - izmin] < threshold ) break;
            }
            float rErr = 0.5*dx*ndx;
            float zErr = 0.5*dz*ndz;
            float aErr = 0.5*da*nda;
            if (astdd > 1.0 && astdd < aErr) aErr = astdd;
            float yErr = xi[ixpsh]*aErr*PI/180.0;

            //  output
            if (output == -2) {  // output original input data with estimated source for next run
                for (i=0; i<nrcv; ++i) {
                    for (icmp=0; icmp<nc; ++icmp) {
                        memcpy(&tro, &hdrs2d[icmp][i], HDRBYTES);
                        memcpy(tro.data, trdata[icmp][i], tro.ns*FSIZE);
                        tro.fx = sy - yc;  // event location from SH
                        tro.fy = sx - xc;
                        tro.fz = sz;
                        tro.dx = rErr;
                        tro.dy = yErr;
                        tro.dz = zErr;
                        tro.unscale = azpsh;
                        tro.ungpow  = 6.0;  // default half width nagle for recomputation
                        //tro.sx     = NINT(sx/scale);  // event location X from SH
                        //tro.sy     = NINT(sy/scale);  // event location Y from SH
                        //tro.sdepth = NINT(sz/scale);  // event location Z from SV
                        puttr(&tro);
                    }
                }
            }
            if (output == -1) {  // output original input data projected as P/SV/SH towards estimated source for QC
                for (i=0; i<nrcv; ++i) {
                    for (icmp=0; icmp<nc; ++icmp) {
                        memcpy(&tro, &hdrs2d[icmp][i], HDRBYTES);
                        memcpy(tro.data, trdata[icmp][i], tro.ns*FSIZE);
                        tro.fx = sx;  // event location from SH
                        tro.fy = sy;
                        tro.fz = sz;
                        tro.dx = rErr;
                        tro.dy = yErr;
                        tro.dz = zErr;
                        tro.gdel = NINT(xi[ixpsv]/scale);
                        tro.sdel = NINT(xi[ixpsh]/scale);
                        tro.otrav = NINT(10.0*azpsh);
                        tro.gaps  = NINT(10.0*aErr);
                        tro.unscale = azpsh;
                        tro.ungpow  = imgpsh[iypsv][ixpsv - ixmin][izpsv - izmin];
                        tro.duse = icmp + 1;
                        puttr(&tro);
                    }
                    float x = sx - rxyz[i][0];
                    float y = sy - rxyz[i][1];
                    float z = sz - rxyz[i][2];
                    float ah = (180.0/PI)*atan2(y, x);
                    float av = (180.0/PI)*atan2(z, sqrtf(x*x + y*y));
                    Rotate2C(0, nt-1, 0, ah, trdata[2][i], trdata[1][i], trdata[1][i], trdata[2][i]);
                    Rotate2C(0, nt-1, 0, av, trdata[1][i], trdata[0][i], trdata[0][i], trdata[1][i]);
                    for (icmp=0; icmp<3; ++icmp) {
                        tro.duse = nc + icmp + 1;
                        memcpy(tro.data, trdata[icmp][i], tro.ns*FSIZE);
                        puttr(&tro);
                    }
                }
            } 
            int ns = MAX(nz, ny);
            memcpy(&tro, &hdrs2d[0][0], HDRBYTES); // copy back header to output trace
            tro.ns = (output == 1)? nz : ns;
            tro.fx = sy;  // event location from SH
            tro.fy = sx;
            tro.fz = sz;
            tro.dx = rErr;
            tro.dy = yErr;
            tro.dz = zErr;
            tro.gdel = NINT(xi[ixpsv]/scale);
            tro.sdel = NINT(xi[ixpsh]/scale);
            tro.otrav = NINT(10.0*azpsh);
            tro.gaps  = NINT(10.0*aErr);
            tro.unscale = azpsh;
            tro.ungpow  = imgpsh[iypsv][ixpsv - ixmin][izpsv - izmin];
            tro.d1 = dz;
            tro.d2 = dx;
            tro.f1 = zi[izmin];
            tro.f2 = xi[ixmin];
            tro.dt = NINT(1000.0*dz);
            if (output == -1 || output == 0) { // output slices
                if (output == 0) sfp = stdout;
                else if (!filegiven) {
                    char filename[32];
                    sprintf(filename, "%d-semb.su", keyno);
                    sfp = efopen(filename, "w");
                }
                tro.duse = 1;
                tro.cdpt = 1;
                tro.cdp = 1;
                tro.nvs = iypsh;
                memset(tro.data, 0, ns*FSIZE);
                for (i=0; i < nx; ++i) {  // slice along X/SH
                    tro.nhs = i + 1;
                    tro.f2 = xi[i + ixmin];
                    memcpy(tro.data, &imgpsh[iypsh][i][0], nz*FSIZE);
                    fputtr(sfp, &tro);
                }
                tro.cdp = 2;
                tro.nvs = iypsv;
                for (i=0; i < nx; ++i) {  // slice along X/SV
                    tro.nhs = i + 1;
                    tro.f2 = xi[i + ixmin];
                    memcpy(tro.data, &imgpsv[iypsv][i][0], nz*FSIZE);
                    fputtr(sfp, &tro);
                }
                tro.cdp = 3;
                tro.nvs = iyshv;
                for (i=0; i < nx; ++i) {  // slice along X/SV
                    tro.nhs = i + 1;
                    tro.f2 = xi[i + ixmin];
                    memcpy(tro.data, &imgshv[iyshv][i][0], nz*FSIZE);
                    fputtr(sfp, &tro);
                }
                tro.duse = 2;
                tro.cdpt = 2;
                tro.dt = NINT(1000.0*da);
                tro.f1 = amin;
                tro.d1 = da;
                tro.cdp = 1;
                tro.nvs = izpsh;
                memset(tro.data, 0, ns*FSIZE);
                for (i=0; i < nx; ++i) {  // time slice along X/SH
                    tro.nhs = i;
                    tro.f2 = xi[i + ixmin];
                    for(j=0; j<ny; ++j) tro.data[j] = imgpsh[j][i][izpsh-izmin];
                    fputtr(sfp, &tro);
                }
                tro.cdp = 2;
                tro.nvs = izpsv;
                for (i=0; i < nx; ++i) {  // time slice along X/SV
                    tro.nhs = i;
                    tro.f2 = xi[i + ixmin];
                    for(j=0; j<ny; ++j) tro.data[j] = imgpsv[j][i][izpsv-izmin];
                    fputtr(sfp, &tro);
                }
                tro.cdp = 3;
                tro.nvs = izshv;
                for (i=0; i < nx; ++i) {  // time slice along X/SHV
                    tro.nhs = i;
                    tro.f2 = xi[i + ixmin];
                    for(j=0; j<ny; ++j) tro.data[j] = imgshv[j][i][izshv-izmin];
                    fputtr(sfp, &tro);
                }
                if (output == -1 && !filegiven) fclose(sfp);
            }
            if (output == 1) {
                tro.cdp = 1;
                for (iy = 0; iy < ny; ++iy) {
                    tro.nvs = iy;
                    for (ix = 0; ix < nx; ++ix) {
                        tro.f2 = xi[ix + ixmin];
                        tro.nhs = ix;
                        memcpy(tro.data, &imgpsh[iy][ix][iz], nz*FSIZE);
                        puttr(&tro);
                    }
                }
                tro.cdp = 2;
                for (iy = 0; iy < ny; ++iy) {
                    tro.nvs = iy;
                    for (ix = 0; ix < nx; ++ix) {
                        tro.f2 = xi[ix + ixmin];
                        tro.nhs = ix;
                        memcpy(tro.data, &imgpsv[iy][ix][iz], nz*FSIZE);
                        puttr(&tro);
                    }
                }
            }
            
            /* zero out data memory */
            memset(**imgpsv, 0, nx*ny*nz*FSIZE);
            memset(**imgpsh, 0, nx*ny*nz*FSIZE);
            memset(**imgshv, 0, nx*ny*nz*FSIZE);
            memset(**cdata,  0, 3*nrcv*ntfft*FSIZE);
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

