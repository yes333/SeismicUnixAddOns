/* Copyright (c) Sanyu Ye <sanyu_ye@yahoo.com> 2012 */
/* All rights reserved.                       */

#include "su.h"
#include "sac.h"
#include "sedis4.h"
#include "header.h"
#include <time.h>
#include <stdlib.h>

#include "datetime.c"
#include "hlputil.c"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SEDIS2SAC -  read a SEDIS image file of raw OBS data, convert/split to multiple",
"              SAC files for each component                                     ",
"                                                                               ",
"   sedis2sac fimg= rcvno=                                                      ",
"                                                                               ",
" Required parameter:                                                           ",
"                                                                               ",
" fimg=            input SEDIS image filename                                   ",
" rcvno=           OBS receiver station number                                  ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
" t0=              start time of single SAC file, for active seismic only.      ",
"                  if given, only one SAC file with this starting time and time ",
"                  window 'length' is written                                   ",
"                  format =yyyy,mm,dd,hh,mm                                     ",
" length=48 [hrs]  max. record length of individual sac files (if t0 is given)  ",
"                                                                               ",
" lat=             latitude of OBS, in degrees                                  ",
" long=            longitude of OBS, in degrees                                 ",
" gz=         [m]  water depth at OBS receiver station                          ",
"                                                                               ",
" fsac=obs         prefix of output SAC filename                                ",
"                  default filename prefixNN-yyyymmddhhmmss.sac#                ",
"                  where NN is receiver number (rcvno=)                         ",
"                  and # component number (1=P/hydrophone, 2=Z, 3=X, 4=Y geophone)",
" fhdr=obsNN.hdr   output ascii filename containing SEDIS header info           ",
"                                                                               ",
" gx=         [m]  local/UTM X coordinate (EASTING)                             ",
" gy=         [m]  local/UTM Y coordinate (NORTHING)                            ",
" scale=1          scaling factor for coodinates and elevations                 ",
"                  gx/gy will be stored as integer in header nxsize/nysize      ",
"                  multplied by scaling factor (scale=)                         ",
"                                                                               ",
" chan=2,3,4,1     sample order of components/channels of raw multiplexed data  ",
"                  1=P (hydrophone), 2=Z (vertical), 3=X, 4=Y geophone          ",
"                                                                               ",
" endian=0         =1 [no] byte swap to write SAC files in big endian format    ",
"                                                                               ",
" verbose=0        >0 echo info                                                 ",
"                                                                               ",
" Notes:                                                                        ",
"  SAC (Seismic Analysis Code) file format is primarily made for analysis of    ",
"  earthquake seismic data. Since OBS is widely used for passive seismic study, ",
"  SAC is the natural choice as intermediate file format.                       ",
"  SAC uses only geographycal coordinates (lat/long) while SU the UTM. Also     ",
"  station ID and component are signated by letters in SAC, only numbers are    ",
"  accepted for SU header. Therefore some SAC headerwords are used to store     ",
"  station ID (NORID), component number (NMAGTYP) and UTM coodinates (NX/YSi)   ",
"                                                                               ",
" Version 1.1.0   Last updated Feb, 2013 by Sanyu Ye                            ",
"                                                                               ",
NULL};

/*
 *
 * Credits:
 *	Sanyu Ye <sanyu_ye@yahoo.com>
 *      based on the codes originally written by Xuelin Qiu
 *      SoftSeis, Oslo
 */
/**************** end self doc ***********************************/


// forward declaration of function

int TimeDate2tm(const TimeDate* td, struct tm *tm)
{
    tm->tm_year = td->YEAR;
    tm->tm_mon  = td->MONTH;
    tm->tm_mday = td->MDAY;
    tm->tm_hour = td->HOUR;
    tm->tm_min  = td->MIN;
    tm->tm_sec  = td->SEC;
    return tm->tm_yday;
}

int main(int argc, char **argv)
{
    char *fhdr;                 /* name of header info file	*/
    char *fsac;                 /* name of SAC file	*/
    char *fimg;                 /* name of ascii navigation file	*/

    FILE *fpsac[4] = {NULL, NULL, NULL, NULL};
    FILE *fpimg = NULL;
    FILE *fphdr = NULL;

    int i, j;      /* counter				*/
    //int ns;     /* number of data samples		*/
    int verbose;    /* echo info ...			*/

    int endian;
    int rcvno = 0;
    int nc, nchan[4];   /* number of data components	*/
    float dt, length, scale;
    double gx, gy, gz, glat, glong;

    Sedis_Header sh;
    SACHdr sac = sac_null;

    const size_t SHDR_SIZE = sizeof(Sedis_Header);
    const size_t SAC_SIZE  = sizeof(SACHdr);

     /* Initialize */
    initargs(argc, argv);
    requestdoc(2); /* stdin not used */

    /* Set filenames */
    MUSTGETPARSTRING("fimg", &fimg);
    MUSTGETPARINT("rcvno", &rcvno);

    if (!getparstring("fhdr", &fhdr)) {
        fhdr = malloc(64);
        sprintf(fhdr, "obs%02d.hdr", rcvno);
    }
    /* Open output header file */
    fphdr = efopen(fhdr, "wa");

    if (!getparstring("fsac", &fsac)) fsac = "obs";

    /* Set parameters */
    if ( (nc = countparval("chan")) >= 3 )  {
        getparint("chan", nchan);
    } else {
        nc = 4;   // default 4 component
        nchan[0] = 2;
        nchan[1] = 3;
        nchan[2] = 4;
        nchan[3] = 1;
    }
    if (!getparint("endian", &endian)) endian = 0;
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparfloat("dt", &dt)) dt = 0.004;
    if (!getparfloat("length", &length)) length = 48.0;
    if (!getparfloat("scale", &scale)) scale = 1.0;
    if (scale < 1.0) {
        scale = 1.0;
        fprintf(fphdr, "Scaling factor for distance (scale=) must be >= 1, reset to default scale=1");
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
    double jdayStart = 0.0;   // start time of conversion
    if ((nfread = countparval("t0")) >= 3) {
        int* arri = ealloc1int(8);
        memset(arri, 0, 8*ISIZE);
        getparint("t0", arri);
        jdayStart = julian_day(arri[0], arri[1], arri[2], arri[3], arri[4], 0, 0.0);
        IsActiveSeismic = TRUE;
        fprintf(fphdr, "Starting time for active seismic: %d-%02d-%02d %02d:%02d  Duration: %4.2f hours ~= %d minutes\n",
                arri[0], arri[1], arri[2], arri[3], arri[4], length, (int)(60*length));
        if(verbose) warn("Starting time for active seismic: %d-%02d-%02d %02d:%02d  Duration: %4.2f hours ~= %d minutes",
                arri[0], arri[1], arri[2], arri[3], arri[4], length, (int)(60*length));
    }


    /* open input image file containing raw OBS data */
    fpimg = efopen(fimg, "rb");
    if (verbose) warn("SEDIS image file (%s) opened successfully", fimg);
    fprintf(fphdr,    "SEDIS image file (%s) opened successfully\n", fimg);

    fseek(fpimg, 512, SEEK_SET);  // skip first 512 byte
    size_t nread = fread((char *) &sh, 1, SHDR_SIZE, fpimg);  // read SAC header block

    // check header info
    if (sh.HeaderSize != 80) {
        fprintf(fphdr, "Invalid header size (%d != 80)\n", sh.HeaderSize);
        err("Invalid header size (%d != 80)", sh.HeaderSize);
    }

    int nspb = (int) (sh.BlockSamples1 + 256 * sh.BlockSamples2);       // n of samples in a minute block
    unsigned short bps = sh.SampleBytes;    // byte per sample for all channel together
    double dtb = 60.0 / (double) nspb;      // sampling interval calculated from block samples

    fprintf(fphdr, "Data Header '%s'\n", sh.DataHeader);
    fprintf(fphdr, "Software Revision %d ; Board Serial Number %d\n", sh.Rev, sh.Board);
    fprintf(fphdr, "Available GPS Satellite %d\n", sh.NumberSV);
    fprintf(fphdr, "Number of samples in a minute block %d with sampling interval %5.3f [ms]\n", nspb, 1000.0*dtb);

    unsigned short confw = sh.ConfigWord;
    unsigned short nbpc = ( (confw >> 7) & 0x01 ) ? 4 : 3; // n byte per sample per channel

    unsigned short wbite =  confw & 0x0f;
    unsigned short dtms = 0;  // sampling interval in miliseconds
    switch (wbite) {
        case 0:
            dtms = 16;
            break;
        case 1:
            dtms = 8;
            break;
        case 2:
            dtms = 4;
            break;
        case 3:
            dtms = 2;
            break;
        case 4:
            dtms = 1;
            break;
        default:
            warn("Invalid sampling interval bit map coding %d", wbite);
            break;
    }
    if ( (confw >> 4) & 0x01 ) dtms *= 2;

    if (dtms !=  (unsigned short)(dtb*1000.0 + 0.5)) {
        warn("Sample interval from block %5.3f differ from bit map coding %d [ms]", 1000*dtb, dtms);
        fprintf(fphdr, "Sample interval from block %5.3f differ from bit map coding %d [ms]\n", 1000*dtb, dtms);
    }

    unsigned short chanbm = sh.ChannelBitMap;
    unsigned short cbm1 =  chanbm       & 0x01;
    unsigned short cbm2 = (chanbm >> 1) & 0x01;
    unsigned short cbm3 = (chanbm >> 2) & 0x01;
    unsigned short cbm4 = (chanbm >> 3) & 0x01;
    unsigned short cbm5 = (chanbm >> 4) & 0x01;
    unsigned short cbm6 = (chanbm >> 5) & 0x01;

    unsigned short nchn = cbm1 + cbm2 + cbm3 + cbm4 + cbm5 + cbm6;

    if ( nchn != nc) {
        warn("number channels input differ (nc=%d != %d=nchn) from and overridden by bit map coding", nc, nchn);
        fprintf(fphdr, "number channels input differ (nc=%d != %d=nchn) from and overridden by bit map coding\n", nc, nchn);
        nc = nchn;
    }
    if ( nchn*nbpc != bps) {
        warn("number channels differ (bps/nbpc=%d != %d=nchn) from bit map coding", bps/nbpc, nchn);
        fprintf(fphdr, "number channels differ (bps/nbpc=%d != %d=nchn) from bit map coding\n", bps/nbpc, nchn);
    }

    unsigned short year=sh.SampleTime.YEAR;
    unsigned short mon =sh.SampleTime.MONTH;
    unsigned short mday=sh.SampleTime.MDAY;
    unsigned short hour=sh.SampleTime.HOUR;
    unsigned short min =sh.SampleTime.MIN;
    unsigned short sec =sh.SampleTime.SEC;

    double jday = julian_day(year, mon, mday, hour, min, sec, 0.0);

    fprintf(fphdr,   "        Start block time %d-%02d-%02d %02d:%02d:%02d\n", year, mon, mday, hour, min, sec);
    if(verbose) warn("        Start block time %d-%02d-%02d %02d:%02d:%02d", year, mon, mday, hour, min, sec);
    struct tm tmGPS, tmSDS;
    TimeDate2tm(&sh.GPSTime, &tmGPS);
    TimeDate2tm(&sh.SedisTime, &tmSDS);
    fprintf(fphdr,   "         Drift %d [ms]  GPS Time %d-%02d-%02d %02d:%02d:%02d SEDIS Time %d-%02d-%02d %02d:%02d:%02d\n",
       sh.Drift, tmGPS.tm_year, tmGPS.tm_mon, tmGPS.tm_mday, tmGPS.tm_hour, tmGPS.tm_min, tmGPS.tm_sec,
                 tmSDS.tm_year, tmSDS.tm_mon, tmSDS.tm_mday, tmSDS.tm_hour, tmSDS.tm_min, tmSDS.tm_sec);

    size_t nDataBlockSize = nc*nbpc*nspb;
    char* buf = ealloc1(nDataBlockSize, sizeof(char));
    float *data = ealloc1float(nspb);  // data buffer
    //fread(buf, 1, nDataBlockSize, fpimg);

    double jday_current, jday_last = jday;

    int nmin = 0;
    fseek(fpimg, 512, SEEK_SET);  // rewind back to first header again by skipping first 512 byte
    // check run just read and print all parameters
    while( (nread = fread((char *) &sh, 1, SHDR_SIZE, fpimg)) == SHDR_SIZE
            && (sh.HeaderSize == 80) && (sh.SampleBytes == bps) ) {
        fread(buf, 1, nDataBlockSize, fpimg);
        year=sh.SampleTime.YEAR;
        mon =sh.SampleTime.MONTH;
        mday=sh.SampleTime.MDAY;
        hour=sh.SampleTime.HOUR;
        min =sh.SampleTime.MIN;
        sec =sh.SampleTime.SEC;

        if(nmin%60 == 0) {
            fprintf(fphdr,   "%5d-th block time %d-%02d-%02d %02d:%02d:%02d\n",
                nmin + 1, year, mon, mday, hour, min, sec);
            fprintf(fphdr,   "         Battery Voltage %4.1f ; Temperature %3.1f °C\n",
               (float) sh.Bat*50.0/1024.0, ( (float) sh.Temp*5000.0/1024.0 - 600.0)/10.0);
            fprintf(fphdr,   "         Compass flag =%d : %s\n", sh.Compass, sh.PositionStr);
        }

        jday_current = julian_day(year, mon, mday, hour, min, sec, 0.0);
        float tdiffinsec = (float) ((jday_current - jday_last)*86400.0 - 60.0);
        if ( nmin > 0 && ABS(tdiffinsec) > dtb ) {
            warn("WARNING: %5d-th block time %d-%02d-%02d %02d:%02d:%02d",   nmin + 1, year, mon, mday, hour, min, sec);
            warn("WARNING:   Time discrepancy to previous block = %8.6f [s]",   tdiffinsec);
            fprintf(fphdr, "%5d-th block time %d-%02d-%02d %02d:%02d:%02d\n", nmin + 1, year, mon, mday, hour, min, sec);
            fprintf(fphdr, "WARNING:   Time discrepancy to previous block = %8.6f [s]\n",   tdiffinsec);
        }
        //if(verbose > 10 && nmin%60 == 0) warn("%5d-th block time %d-%02d-%02d %02d:%02d:%02d",
        //        nmin + 1, year, mon, mday, hour, min, sec);
        ++nmin;
        jday_last = jday_current;
    }

    if(verbose) warn("%5d-th/last block time %d-%02d-%02d %02d:%02d:%02d",   nmin, year, mon, mday, hour, min, sec);
    fprintf(fphdr,   "%5d-th/last block time %d-%02d-%02d %02d:%02d:%02d\n", nmin, year, mon, mday, hour, min, sec);
    fprintf(fphdr,   "         Battery Voltage %4.1f ; Temperature %3.1f °C\n",
       (float) sh.Bat*50.0/1024.0, ( (float) sh.Temp*5000.0/1024.0 - 600.0)/10.0);
    fprintf(fphdr,   "         Compass flag =%d : %s\n", sh.Compass, sh.PositionStr);
    TimeDate2tm(&sh.GPSTime, &tmGPS);
    TimeDate2tm(&sh.SedisTime, &tmSDS);
    fprintf(fphdr,   "         Drift %d [ms]  GPS Time %d-%02d-%02d %02d:%02d:%02d SEDIS Time %d-%02d-%02d %02d:%02d:%02d\n",
       sh.Drift, tmGPS.tm_year, tmGPS.tm_mon, tmGPS.tm_mday, tmGPS.tm_hour, tmGPS.tm_min, tmGPS.tm_sec,
                 tmSDS.tm_year, tmSDS.tm_mon, tmSDS.tm_mday, tmSDS.tm_hour, tmSDS.tm_min, tmSDS.tm_sec);
    double jday_end = julian_day(year, mon, mday, hour, min + 1, sec, 0.0);
    double dsec = (jday_end - jday)*86400.0;
    if ( ABS(dsec - 60.0*nmin) > dtb ) {
        fprintf(fphdr,   "Discrepancy in end and start time %5.3f [s] = %d samples",
                dsec - 60.0*nmin, (int) ((dsec - 60.0*nmin)/dtb));
        warn("Discrepancy in end and start time %5.3f [s] = %d samples",
                dsec - 60.0*nmin, (int) ((dsec - 60.0*nmin)/dtb));
    }
    if(verbose) warn("Total record length %8.6f days = %3.2f hours",   jday_end - jday, dsec/3600.0);
    fprintf(fphdr,   "Total record length %8.6f days = %3.2f hours\n", jday_end - jday, dsec/3600.0);

    fflush(fphdr);  // flush write buffer to header file

    char name[32];
    
    // set common header words
    sprintf(name, "%d", rcvno);
    sprintf(sac.kstnm, "OBS%d", rcvno);
    sac.norid = rcvno;
    sac.delta = dtb;
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

    int NS = 0, nmin_read = 0, nfiles = 0;
    nmin = 0;
    cwp_Bool eof = cwp_false;
    // read from begining and convert to SAC file
    fseek(fpimg, 512, SEEK_SET);  // skip first 512 byte
    
    do {
        nread = fread((char *) &sh, 1, SHDR_SIZE, fpimg);
        eof = !( (nread == SHDR_SIZE) && (sh.HeaderSize == 80) && (sh.SampleBytes == bps) );
        if (!eof) {
            fread(buf, 1, nDataBlockSize, fpimg);
            ++nmin_read;

            year=sh.SampleTime.YEAR;
            mon =sh.SampleTime.MONTH;
            mday=sh.SampleTime.MDAY;
            hour=sh.SampleTime.HOUR;
            min =sh.SampleTime.MIN;
            sec =sh.SampleTime.SEC;
            jday = julian_day(year, mon, mday, hour, min, sec, 0.0);

            if (jday < jdayStart) continue;  // skip
            if ( IsActiveSeismic && (jday - jdayStart)*1440.0 >= floor((double) length * 60.0) ) eof = cwp_true;  // end conversion
        }

        if ( nmin % ((int) (60.0*length)) == 0 || eof ) { // close previous file and start a new sac file
            if (nmin > 0) {  // update SAC header and rewrite at beginning of file
                for (i = 0; i < nc; ++i) {
                    sac.imagtyp  = nchan[i];  // component
                    sac.npts = NS;  // update header word for number of samples
                    sac.e = (NS - 1) * dtb;
                    setCompName(nchan[i], sac.kcmpnm);

                    rewind(fpsac[i]);
                    if (endian) {
                        float* toSwap  = (float*) &sac;
                        for(j=0; j < 110; ++j, ++toSwap) swap_float_4(toSwap);
                    }
                    efwrite(&sac, 1, SAC_SIZE, fpsac[i]);
                    efclose(fpsac[i]);
                    fpsac[i] = NULL;
                }

                ++nfiles;
                NS = 0;
                if (eof) break;
            }

            if (!eof) {
                sac.nzyear = year;
                sac.nzjday = dayofyear(year, mon, mday);
                sac.nzhour = hour;
                sac.nzmin  = min;
                sac.nzsec  = sec;
                sac.nzmsec = 0.0;
                for (i = 0; i < nc; ++i) {
                    sprintf(name, "%s%02d-%d%02d%02d%02d%02d%02d.sac%d",
                          fsac, rcvno, year, mon, mday, hour, min, sec, nchan[i]);
                    fpsac[i] = efopen(name, "wb");

                    efwrite(&sac, 1, SAC_SIZE, fpsac[i]);
                }
            }
        }

        if (!eof) {
            ++nmin;
            NS += nspb;

            for (j = 0; j < nc; ++j) {
                for (i = 0, p = (unsigned char*) buf; i < nspb; ++i) {
                    idata = 0;
                    //memcpy(&idata, p + nbpc*(nc*i + j) - (4 - nbpc), 4); // assume little endian
                    //if (nbpc == 3) idata = idata >> 8;
                    memcpy(&idata, p + nbpc*(nc*i + j), nbpc); // assume big endian
                    swap_int_4(&idata);
                    if (nbpc == 3) idata = idata >> 8;
                    data[i] = (float) idata;
                    //data[i] = (SF * data[i] - DC) * SSF;

                    if (endian) swap_float_4(&data[i]);
                }
                efwrite(data, 1, nspb*FSIZE, fpsac[j]);
            }
        }
    } while (!eof);


    fprintf(fphdr,    "Totally %d x %d files of %d minutes (or %5.3f hours or %5.3f days) of data have been written",
            nc, nfiles, nmin, (double) nmin/60.0, (double) nmin/1440.0);
    if (verbose) warn("Totally %d x %d files of %d minutes (or %5.3f hours or %5.3f days) of data have been written",
            nc, nfiles, nmin, (double) nmin/60.0, (double) nmin/1440.0);

    efclose(fpimg);
    efclose(fphdr);
    for (i=0; fpsac[i] && i < nc; ++i) efclose(fpsac[i]);

    return (CWP_Exit());
}
