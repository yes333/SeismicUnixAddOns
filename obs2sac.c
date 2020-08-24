/* Copyright (c) Sanyu Ye <sanyu_ye@yahoo.com> 2012 */
/* All rights reserved.                       */

#include "su.h"
#include "sac.h"
#include "header.h"
#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "datetime.c"
#include "hlputil.c"


/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" OBS2SAC - read serial files of raw data of CAS OBS, convert/split to multiple ",
"              SAC files for each component                                     ",
"                                                                               ",
"   obs2sac files= rcvno=                                                      ",
"                                                                               ",
" Required parameter:                                                           ",
"                                                                               ",
" files=           filename pattern of raw data of CAS OBS                      ",
" rcvno=           OBS receiver station number                                  ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
" chan=3,4,2,1     sample order of components/channels of raw multiplexed data  ",
"                  1=P (hydrophone), 2=Z (vertical), 3=X (inline), 4=Y (crossline)",
"                  for broadband geophone 5=Z, 6=X, 7=Y                         ",
"                  channel number is appended to output file extension          ",
"                                                                               ",
" lat=             latitude of OBS, in degrees                                  ",
" long=            longitude of OBS, in degrees                                 ",
" gz=         [m]  water depth at OBS receiver station                          ",
"                                                                               ",
" fsac=obs         prefix of output SAC filename                                ",
"                  default filename prefixNN-yyyymmddhhmmss.sac#                ",
"                  where NN is receiver number (rcvno=)                         ",
"                  and # channel number (1=P/hydrophone, 2=V, 3=X, 4=Y geophone)",
" fhdr=obsNN.hdr   output ascii filename containing info of conversion          ",
"                                                                               ",
" gx=         [m]  local/UTM X coordinate (EASTING)                             ",
" gy=         [m]  local/UTM Y coordinate (NORTHING)                            ",
" scale=1          scaling factor for coodinates and elevations                 ",
"                  gx/gy will be stored as integer in header nzsize/nysize      ",
"                  multplied by scaling factor (scale=)                         ",
"                                                                               ",
" sps=125          sampling rate (samples per second)                           ",
" nbpc=3           numbers of bytes of data sample of single channel            ",
"                                                                               ",
" length=48 [hous] max. length of individual sac files                          ",
" t0=              start time of single SAC file, for active seismic only.      ",
"                  if given, only one SAC file with this starting time and time ",
"                  window 'length' is written                                   ",
"                  format =yyyy,mm,dd,hh,mm                                     ",
"                                                                               ",
" endian=0         =1 [no] byte swap to write SAC files in big endian format    ",
"                                                                               ",
" verbose=0        >0 echo info                                                 ",
"                                                                               ",
" Notes:                                                                        ",
"  SAC (Seismic Analysis Code) file format is primarily made for analysis of    ",
"  earthquake seismic data. Since OBS is widely used for passive seismic study, ",
"  SAC is the natural choice as file format.                                    ",
"  The real sampling interval is determined by time difference of first sample  ",
"  between second and first data files. Usually it is bit larger then the norminal",
"  sampling interval as input by 1/sps. In case of error messages for missing or",
"  overlapping sample, please check header file and remove first few files and  ",
"  start conversion at file which has normal/average duration/sampling interval.",
"  It is recommended to convert recorded raw data file by file without merge, as",
"  there may be a tiny time gap (usually smaller than one sample interval)      ",
"  between two consecutive raw data files.                                      ",
"                                                                               ",
" Examples:                                                                     ",
" 1. convert 4 component OBS data                                               ",
" obs2sac files=OBSRAW/2*.??? rcvno=27 t0=2009,4,17,17,40 length=24 verbose=1   ",
"                                                                               ",
" 2. convert 3 component broadband OBS data                                     ",
" obs2sac files=OBSRAW/A*.??? rcvno=27 chan=6,7,5 fsac=broadband fhdr=bb27.hdr  ",
"                                                                               ",
" Version 1.0.0   Last updated Jan, 2013 by Sanyu Ye                            ",
"                                                                               ",
NULL};

/*
 *
 * Credits:
 *	Sanyu Ye sanyu_ye@yahoo.com
 *      based on the codes originally written by Xuelin Qiu
 *      SoftSeis, Oslo
 */
/**************** end self doc ***********************************/


// forward declaration of function
void fname2time(const int nfiles, char** filenames, double* ftime);

int main(int argc, char **argv)
{
    char *fhdr;                 /* name of header info file	*/
    char *fsac;                 /* name of SAC file	*/
    char *fpattern;             /* file name pattern of OBS raw data files	*/

    FILE *fpsac[4] = {NULL, NULL, NULL, NULL};
    FILE *fpimg = NULL;
    FILE *fphdr = NULL;

    int i, j, k;        /* counter				*/
    int nc, nchan[4];   /* number of data components	*/
    int sps, nbpc;      /* sampling rate, number byte per channel sample*/
    int verbose;        /* echo info ...			*/

    int endian;
    int rcvno = 0;
    float length, scale;
    double gx, gy, gz, glat, glong, dt;

    SACHdr sac = sac_null;
    const size_t SAC_SIZE = sizeof(SACHdr);

     /* Initialize */
    initargs(argc, argv);
    requestdoc(2); /* stdin not used */

    /* Set filenames */
    MUSTGETPARSTRING("files", &fpattern);
    MUSTGETPARINT("rcvno", &rcvno);

    if (!getparstring("fhdr", &fhdr)) {
        fhdr = malloc(64);
        sprintf(fhdr, "obs%02d.hdr", rcvno);
    }
    if (!getparstring("fsac", &fsac)) fsac = "obs";

    /* Set parameters */
    if ( (nc = countparval("chan")) >= 3 )  {
        getparint("chan", nchan);
    } else {
        nc = 4;   // default 4 component
        nchan[0] = 3;
        nchan[1] = 4;
        nchan[2] = 2;
        nchan[3] = 1;
    }
    if (!getparint("nbpc", &nbpc)) nbpc = 3;
    if (!getparint("sps", &sps)) sps = 125;
    if (!getparint("endian", &endian)) endian = 0;
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparfloat("length", &length)) length = 48.0;
    if (!getparfloat("scale", &scale)) scale = 1.0;
    if (scale < 1.0) {
        scale = 1.0;
        warn("Scaling factor for distance (scale=) must be >= 1, reset to default scale=1");
    }
    const double DIST_SCALE = scale;  // scaling factor to distances

    if (!getpardouble("gx", &gx)) gx = 0.0;
    if (!getpardouble("gy", &gy)) gy = 0.0;
    if (!getpardouble("gz", &gz)) gz = 0.0;
    if (!getpardouble("lat",  &glat))   glat  = 0.0;
    if (!getpardouble("long", &glong))  glong = 0.0;

    cwp_Bool IsActiveSeismic = FALSE;
    int nfread = 0;
    int* arri = ealloc1int(8);
    memset(arri, 0, 8*ISIZE);
    double jdayStart = 0.0;   // start time of conversion
    if ((nfread = countparval("t0")) >= 3) {
        getparint("t0", arri);
        jdayStart = julian_day(arri[0], arri[1], arri[2], arri[3], arri[4], 0, 0.0);
        IsActiveSeismic = TRUE;
        warn("Starting Time: %d-%02d-%2d %02d:%02d:00   Duration=%4.1f [hour]",
                arri[0], arri[1], arri[2], arri[3], arri[4], length);
    }

    // list files with file pattern
    char cmdbuf[64];
    sprintf(cmdbuf, "ls %s | wc -w", fpattern);
    FILE* pipefp = epopen(cmdbuf, "r");
    int nfiles = 0;
    fscanf(pipefp, "%d", &nfiles);
    epclose(pipefp);
    if ( nfiles <= 0 ) {
        err(" No file found for pattern %s", fpattern);
    }
    if(verbose) warn("Total %d raw data files found for file pattern '%s'", nfiles, fpattern);

    sprintf(cmdbuf, "ls %s", fpattern);
    pipefp = epopen(cmdbuf, "r");
    int* nsf = ealloc1int(nfiles);      // array holding number of sample for each file
    char** filenames = (char**) ealloc2(64, nfiles, sizeof(char));
    for (i=0; i<nfiles; ++i) {
        fscanf(pipefp, "%s", filenames[i]);
    }
    epclose(pipefp);

    // figure out start time of every raw data file from filename
    double* ftime = ealloc1double(nfiles);
    fname2time(nfiles, filenames, ftime);

    /* Open output header file */
    fphdr = efopen(fhdr, "w+");
    if (jdayStart) fprintf(fphdr,   "Starting Time: %d-%02d-%2d %02d:%02d:00   Duration=%4.1f [hour]\n\n",
                        arri[0], arri[1], arri[2], arri[3], arri[4], length);
    fprintf(fphdr,   "Total %d raw data files found for file pattern '%s'\n", nfiles, fpattern);

    // print out info for raw data files
    int nsdiff;  // number of samples in a file
    double msec;
    struct tm tmf;
    struct stat statf;
    stat(filenames[0], &statf);
    int ns = statf.st_size / nc / 3;
    dt = (nfiles > 1)? (ftime[1] - ftime[0])*86400.0/(double)ns : 1.0/(double) sps;
    if (NINT(1.0/dt) != sps ) {
        fprintf(fphdr, "\nSampling rate = %1.0f  != %d as input\n", 1.0/dt, sps);
        warn("Sampling rate = %1.0f  != %d as input", 1.0/dt, sps);
    }
    for (i=0; i<nfiles; ++i) {
        jday2tm(ftime[i], &msec, &tmf);
        stat(filenames[i], &statf);
        nsf[i] = statf.st_size / nc / 3;
        double tdiff = (i < nfiles - 1)? (ftime[i+1] - ftime[i])*86400.0 : 0.0;
        fprintf(fphdr, "%3d '%s' NS=%d Start time %d-%02d-%02d %02d:%02d:%02d.%03.0f Tdiff=%6.3f\n",
            i + 1, filenames[i], nsf[i], tmf.tm_year + 1900, tmf.tm_mon + 1, tmf.tm_mday,
            tmf.tm_hour, tmf.tm_min, tmf.tm_sec, msec, tdiff);
        if (verbose) warn("%3d '%s' NS=%d Start time %d-%02d-%02d %02d:%02d:%02d.%03.0f Tdiff=%6.3f",
            i + 1, filenames[i], nsf[i], tmf.tm_year + 1900, tmf.tm_mon + 1, tmf.tm_mday,
            tmf.tm_hour, tmf.tm_min, tmf.tm_sec, msec, tdiff);
        nsdiff = (i<nfiles-1)? NINT( (tdiff - nsf[i]*dt)/dt ) : 0 ;
        if (ABS(nsdiff) != 0) {
            warn("    WARNING: Missing %d samples to next file!!!", nsdiff);
            fprintf(fphdr, "    WARNING: Missing %d samples to next file!!!\n", nsdiff);
        }
        if (statf.st_size % (nc * nbpc) != 0) {
            fprintf(fphdr, "    WARNING: File size not multiple of sample size!!!\n");
            warn("    WARNING: File size not multiple of sample size!!!");
        }
    }
    fflush(fphdr);

    char name[32];

    // set common header words
    sprintf(name, "%d", rcvno);
    sprintf(sac.kstnm, "OBS%d", rcvno);

    sac.norid = rcvno;
    sac.delta = (float) dt;
    sac.b     = 0.0;   // starting time always zero
    if ( glat  != 0.0) sac.stla = glat;
    if ( glong != 0.0) sac.stlo = glong;
    if ( gz != 0.0)    sac.stdp = gz;
    if ( gx != 0.0) sac.nxsize = (int) (DIST_SCALE * gx);
    if ( gy != 0.0) sac.nysize = (int) (DIST_SCALE * gy);
    sac.scale = DIST_SCALE;
    sac.leven = TRUE;
    sac.iftype = ITIME;
    sac.iztype = IB;
    sac.lpspol = FALSE;
    sac.lcalda = TRUE;
    sac.unused27 = FALSE;
    sac.idep = IVOLTS;

    
    unsigned char *p;
    int  idata;
    size_t nDataSampleSize = nc * nbpc;
    char* buf = ealloc1(nDataSampleSize, sizeof(char));
    float data = 0.0;  // data buffer
    double dtnew = 0.0;

    int fclosed = 0;
    int NS = 0;
    for ( k = 0; k < nfiles; ++k ) {  // loop over files
        if (ftime[k] + nsf[k]*dt/86400.0 < jdayStart) continue;  // skip
        if (IsActiveSeismic && (ftime[k] - jdayStart)*24.0 > length) break;  // end conversion

        if (verbose) warn("Converting %d-th data file %s", k+1, filenames[k]);

        /* open input image file containing raw OBS data */
        fpimg = efopen(filenames[k], "rb");

        if (!IsActiveSeismic) { // passive, convert file by file
            jday2tm(ftime[k], &msec, &tmf);

            sac.nzyear = tmf.tm_year + 1900;
            sac.nzjday = dayofyear(tmf.tm_year + 1900, tmf.tm_mon + 1, tmf.tm_mday);
            sac.nzhour = tmf.tm_hour;
            sac.nzmin  = tmf.tm_min;
            sac.nzsec  = tmf.tm_sec;
            sac.nzmsec = NINT(msec);
            sac.npts = nsf[k];  // update header word for number of samples

            // open sac files
            for (j = 0; j < nc; ++j)
            {
                sprintf(name, "%s%d-%02d%02d%02d%02d%02d%02d%03.0f.sac%d",
                    fsac, rcvno, tmf.tm_year - 100, tmf.tm_mon + 1, tmf.tm_mday,
                    tmf.tm_hour, tmf.tm_min, tmf.tm_sec, msec, nchan[j]);

                sac.iinst  = nchan[j];  // component
                setCompName(nchan[j], sac.kcmpnm);

                fpsac[j] = efopen(name, "wb");
                efwrite(&sac, 1, SAC_SIZE, fpsac[j]);
            }

            // read sample by sample, channel by channel
            for (j = 0; j < nsf[k]; ++j) {
                fread(buf, 1, nDataSampleSize, fpimg);
                for (i = 0, p = (unsigned char*) buf; i < nc; ++i) {
                    idata = 0;
                    memcpy(&idata, p + i*nbpc, nbpc); // assume big endian
                    swap_int_4(&idata);
                    if (nbpc == 3) idata = idata >> 8;
                    data = (float) idata;
                    //data[i] = (SF * data[i] - DC) * SSF;

                    if (endian) swap_float_4(&data);
                    efwrite(&data, 1, FSIZE, fpsac[i]);
                }
            }
            efclose(fpimg);
            for (i = 0; i < nc; ++i) efclose(fpsac[i]);
            fclosed = 1;
        } else { // convert multiple files into one single SAC file of given length
            double tdiff = (k < nfiles - 1)? (ftime[k+1] - ftime[k])*86400.0 : 0.0;
            nsdiff = (k < nfiles-1)? NINT( (tdiff - nsf[k]*dt)/dt ) : 0 ;
            if (!dtnew) {
                dtnew = tdiff/nsf[k];
                if (ABS(nsdiff) != 0) {
                    warn("Adjust sampling interval from old %10.8f to new %10.8f [ms]", dt*1000.0, dtnew*1000.0);
                    fprintf(fphdr, "Adjust sampling interval from old %10.8f to new %10.8f [ms]", dt*1000.0, dtnew*1000.0);
                }
                dt = dtnew;
                sac.delta = dt;
            } else if (ABS(nsdiff) != 0) {
                warn("WARNING: Missing %d samples from %d-th file %s to next file %s!!!",
                        nsdiff, k+1, filenames[k], filenames[k+1]);
                fprintf(fphdr, "WARNING: Missing %d samples from %d-th file %s to next file %s!!!",
                        nsdiff, k+1, filenames[k], filenames[k+1]);
            }
            
            int offset = NINT(86400.0*(jdayStart - ftime[k])/dt);

            if (NS == 0) {  // sac file not open yet
                if (offset > 0) {
                    fseek(fpimg, offset*nDataSampleSize, SEEK_SET);
                    jday2tm(jdayStart, &msec, &tmf);
                } else {
                    jday2tm(ftime[k], &msec, &tmf);
                }
                sac.nzyear = tmf.tm_year + 1900;
                sac.nzjday = dayofyear(tmf.tm_year + 1900, tmf.tm_mon + 1, tmf.tm_mday);
                sac.nzhour = tmf.tm_hour;
                sac.nzmin  = tmf.tm_min;
                sac.nzsec  = tmf.tm_sec;
                sac.nzmsec = NINT(msec);
                sac.npts = nsf[k];  // update header word for number of samples

                // open sac files
                for (j = 0; j < nc; ++j)
                {
                    sprintf(name, "%s%d+%02d%02d%02d%02d%02d%02d%03.0f.sac%d",
                        fsac, rcvno, tmf.tm_year - 100, tmf.tm_mon + 1, tmf.tm_mday,
                        tmf.tm_hour, tmf.tm_min, tmf.tm_sec, msec, nchan[j]);

                    sac.iinst  = nchan[j];  // component
                    setCompName(nchan[j], sac.kcmpnm);

                    fpsac[j] = efopen(name, "wb");
                    efwrite(&sac, 1, SAC_SIZE, fpsac[j]);
                }
            }

            while (fread(buf, 1, nDataSampleSize, fpimg) > 0) {
                for (i = 0, p = (unsigned char*) buf; i < nc; ++i) {
                    idata = 0;
                    memcpy(&idata, p + i*nbpc, nbpc); // assume big endian
                    swap_int_4(&idata);
                    if (nbpc == 3) idata = idata >> 8;
                    data = (float) idata;
                    //data[i] = (SF * data[i] - DC) * SSF;

                    if (endian) swap_float_4(&data);
                    efwrite(&data, 1, FSIZE, fpsac[i]);
                }

                ++NS;

                if (NS*dt > 3600.0 * length) { // close sac files
                    // update SAC header
                    for (j = 0; j < nc; ++j)
                    {
                        sac.npts = NS;
                        sac.iinst  = nchan[j];  // component
                        setCompName(nchan[j], sac.kcmpnm);
                        
                        rewind(fpsac[j]);
                        efwrite(&sac, 1, SAC_SIZE, fpsac[j]);
                        efclose(fpsac[j]);
                    }
                    fclosed = 1;
                    break;
                }
            }
            efclose(fpimg);
        }
    }
    
    if (!fclosed) {
        for (j = 0; j < nc; ++j)
        {
            sac.npts = NS;
            sac.iinst  = nchan[j];  // component
            setCompName(nchan[j], sac.kcmpnm);

            rewind(fpsac[j]);
            efwrite(&sac, 1, SAC_SIZE, fpsac[j]);
            efclose(fpsac[j]);
        }
    }

    efclose(fphdr);

    //if (verbose) warn("Totally %d minutes (%7.3f days) of data have been written",
    //       nmin, (double) nmin/86400.0);

    return (CWP_Exit());
}

void fname2time(const int nfiles, char **filenames, double *ftime)
{
    int i, j, fnlength;
    int year, mon, mday, hour, min, sec;
    double msec;
    unsigned int date[8], ms[3];
    unsigned char c;

    for (j=0; j<nfiles; ++j) {

        fnlength = strlen(filenames[j]);

        /* get the date and time information for the file name */

        for (i=0;i<8;i++) {
            c = filenames[j][fnlength-12+i];
            if ( !isxdigit(c) ) {
                err("Invalid (non-xdigit) character found in last 12 of filename %s", filenames[j]);
            }
            date[i] = (c < 58)? c - 48 : toupper(c) - 55;
        }

        year = (date[0] << 2 | date[1] >> 2)  & 0x0f;  /* 0x0f=00001111 */
        mon  = (date[1] << 2 | date[2] >> 2)  & 0x0f;  /* 0x0f=00001111 */
        mday = (date[2] << 3 | date[3] >> 1)  & 0x1f;  /* 0x1f=00011111 */
        hour = (date[3] << 4 | date[4])       & 0x1f;  /* 0x1f=00011111 */
        min  = (date[5] << 2 | date[6] >> 2)  & 0x3f;  /* 0x3f=00111111 */
        sec  = (date[6] << 4 | date[7])       & 0x3f;  /* 0x3f=00111111 */

        year += 2000;

        for (i=9;i<12;i++) {
            c = filenames[j][fnlength-12+i];
            if ( !isxdigit(c) ) {
                err("Invalid (non-xdigit) character found in last 11 of filename %s", filenames[j]);
            }
            ms[i-9] = (c < 58)? c - 48 : toupper(c) - 55;
        }

        msec = (double) (ms[0] << 8 | ms[1] << 4 | ms[2]);
        msec /= 2.048;

        ftime[j] = julian_day(year, mon, mday, hour, min, sec, msec);
    }
    return;
}
