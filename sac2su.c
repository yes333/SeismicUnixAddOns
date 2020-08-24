/* Copyright (c) Sanyu Ye <sanyu_ye@yahoo.com> 2012 */
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"
#include "sac.h"
#include "header.h"
#include <time.h>
#include <stdlib.h>

#include "datetime.c"
#include "hlputil.c"


/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SAC2SU -  read an SAC file of single component OBS data, split into multiple  ",
"           su traces according to navigation (source & receiver) parameter file",
"                                                                               ",
"   sac2su fsac= fnav=  > stdout                                                ",
"                                                                               ",
" Required parameter:                                                           ",
"                                                                               ",
" fsac=            input SAC record filename                                    ",
" fnav=            input navigation filename                                    ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
" length=10   [s]  length of trace                                              ",
" t0=-1.0     [s]  starting time in reference to reduced zero time              ",
" rv=8.0   [km/s]  reduction velocity                                           ",
" maxd=200   [km]  maximum distance to be output                                ",
"                                                                               ",
" rcvno=           OBS station number, default exteacted from SAC header KSTNM  ",
" cmpno=           OBS component number, default read from SAC header IINST     ",
"                  1=P (hydrophone), 2=Z (vertical), 3=X, 4=Y geophone          ",
"                                                                               ",
" dt=         [s]  sampling interval if not given in SAC file                   ",
" gx=         [m]  local/UTM X coordinate (EASTING)                             ",
" gy=         [m]  local/UTM Y coordinate (NORTHING)                            ",
" gz=         [m]  water depth at OBS receiver station                          ",
" scale=1          scaling factor for coodinates and elevations                 ",
"                                                                               ",
" drift=     [ms]  clock drift of OBS station, positive if faster than GPS clock",
" td1=             date & time when OBS clock is synchronized before deployment ",
" td2=             date & time when OBS clock drift is checked after recovery   ",
"                  format =yyyy,mm,dd,hh,mm                                     ",
"                                                                               ",
" endian=0         =1 [no] byte swap due to different endian of SAC file        ",
"                                                                               ",
" verbose=0        >0 echo info                                                 ",
"                                                                               ",
" Notes:                                                                        ",
"  SAC (Seismic Analysis Code) file format is primarily made for analysis of    ",
"  earthquake seismic data. Since OBS is widely used for passive seismic study, ",
"  SAC is the natural choice as file format. For active OBS seismic studies, SU ",
"  or SEGY format is standard, so SAC data need to be split/converted to SU     ",
"  traces that correspond to individual shots.                                  ",
"  SAC uses only geographycal coordinates (lat/long) while SU the UTM. Also     ",
"  station ID and component are signated by letters in SAC, only numbers are    ",
"  accepted for SU header. Therefore some SAC headerwords are used to store     ",
"  station ID (KSTNM), component number (IINST) and UTM coodinates (INTERNAL7/8)",
"  Navigation file must contain source (shot) parameters like UTM coordinates,  ",
"  accurate shot time with miliseconds, and optionally water depth beneath shot ",
"  point. OBS receiver information, like UTM coodinates, water depth and clock  ",
"  drift can directly given through console parameters. However it is strongly  ",
"  recommended to put them in the navigation file for better tracability.       ",
"  Receiver info in navigation file overide those given by command line or      ",
"  stored in SAC file.                                                          ",
"  Navigation file is white-space-separated ascii text, not fixed format, easily",
"  to be created using unix utilities like awk. It contains three types of      ",
"  records: clock drift of receiver stations, their local/UTM coordinates, and  ",
"  coordinates and exact times of all shotpoints. The records are designated with",
"  capital letter 'D', 'R' and S' respectively at first column, followed by     ",
"  station or shot number. Comments can be added in the file, best with starting",
"  letter 'H' at begining, similar to ukooa/P190 file. the detailed format of   ",
"  navigation file looks as below:                                              ",
"                                                                               ",
" H info regarding survey parameters can be presented here                      ",
" H it is a good practice to put such info in one piece for later references    ",
" H for yourself, your cooperation partners and third paty users                ",
" H ......                                                                      ",
" D NNN time-drift-in-ms yyyy mm dd hh mm yyyy mm dd hh mm                      ",
" D ......                                                                      ",
" R NNN EASTING NORTHING WATER-DEPTH-IN-METER                                   ",
" R ......                                                                      ",
" S NNNN EASTING NORTHING WATER-DEPTH yyyy mm dd hh mm ss.ms                    ",
" S .......                                                                     ",
"                                                                               ",
" Version 1.0.0   Last updated Jan, 2013 by Sanyu Ye                            ",
"                                                                               ",
NULL};

/*
 *
 * Credits:
 *	Sanyu Ye <sanyu_ye@yahoo.com>
 *      SoftSeis, Oslo
 */
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
    cwp_String key;		/* header key word from segy.h		*/
    cwp_String type;            /* type of key				*/
    Value val;                  /* value of key			*/
    int index;
    char *fsac;                 /* name of SAC file	*/
    char *fnav;                 /* name of ascii navigation file	*/

    FILE *fpsac = NULL;
    FILE *fpnav = NULL;

    int i;      /* counter				*/
    int ns;     /* number of data samples		*/
    int verbose;    /* echo info ...			*/

    int endian;
    int rcvno = 0, srcno = 0;
    float dt, length, t0, rv, scale, maxd;
    
    segy tr;
    SACHdr rcdhdr;

     /* Initialize */
    initargs(argc, argv);
    requestdoc(2); /* stdin not used */

    /* Make sure stdout is a file or pipe */
    switch (filestat(STDOUT)) {
        case TTY:
            err("stdout can't be tty");
            break;
        case DIRECTORY:
            err("stdout must be a file, not a directory");
            break;
        case BADFILETYPE:
            err("stdout is illegal filetype");
            break;
        default: /* Others OK */
            break;
    }

    /* Set filenames */
    MUSTGETPARSTRING("fsac", &fsac);
    MUSTGETPARSTRING("fnav", &fnav);
    if (!getparstring("key", &key))	  key = "corr";
    type = hdtype(key);
    index = getindex(key);

    /* Set parameters */
    if (!getparint("rcvno", &rcvno)) rcvno = 0;
    if (!getparint("ns", &ns)) ns = 0;
    if (!getparint("endian", &endian)) endian = 0;
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparfloat("length", &length)) length = 10.0;
    if (!getparfloat("maxd", &maxd)) maxd = 200.0;
    if (!getparfloat("t0", &t0)) t0 = -1.0;
    if (!getparfloat("rv", &rv)) rv = 8.0;
    if (!getparfloat("scale", &scale)) scale = 1.0;
    if ( scale < 1.0) {
        scale = 1.0;
        warn("Scaling factor for distance (scale=) must be >= 0, reset to default scale=1");
    }
    const double DIST_SCALE = scale;  // scaling factor to distances


    /* Open input SAC file */
    fpsac = efopen(fsac, "r");
    int nread = (int) fread((char *) &rcdhdr, 1, sizeof(SACHdr), fpsac);  // read SAC header block
    /* open navigation file containing source info */
    fpnav = efopen(fnav, "ra");
    if (verbose) warn("navigation file (%s) opened successfully", fnav);

    // swap byte if SAC created on big endian machine
    if (endian) {
        float* toSwap  = (float*) &rcdhdr;
        for(i=0; i < 110; ++i, ++toSwap) swap_float_4(toSwap);
    }

    // decipher receiver station number and component
    int RCVNO = (rcdhdr.norid > 0)? rcdhdr.norid : a2i(rcdhdr.kstnm, "0123456789", 8);
    if (RCVNO > 0) {
        if(verbose) warn("Receiver station number %d given in SAC file", RCVNO);
    } else if (RCVNO == 0 && rcvno == 0) {
        err("Receiver station number not given in SAC file, must be specified by rcvno=");
    } else if (RCVNO > 0 && rcvno > 0) {
        if (RCVNO != rcvno) {
            warn("Receiver station number (=%d) in SAC header overridden to the specified (rcvno=%d)", RCVNO, rcvno);
            RCVNO = rcvno;
        }
    } else if (RCVNO == 0) RCVNO = rcvno;

    // get number of sampls and sampling rate
    int NPTS = rcdhdr.npts;
    float DT = rcdhdr.delta;
    if (DT <= 0.0) {
        if (!getparfloat("dt", &dt)) {
            err("Sampling interval dt= must be specified");
        }
        DT = dt;
    }
    int NS = NINT(length/DT + 0.5);  // calc number of sample of trace

    // read nav file to find clock drift info of corresponding station
    int DRIFT = 0, drift;
    int year1, mon1, day1, hour1, min1;
    int year2, mon2, day2, hour2, min2;
    char* lname = malloc(16);
    char* line = (char*) ealloc1(1,128);
    memset(line, 0, 128);
    int ntoread = 85, nbread = 0, nfread = 0;
    while ( (nbread = getline(&line, &ntoread, fpnav)) > 0 ) {
        if (line[0] == 'D') {  // clock drift line
            nfread = sscanf(line, "%s %d %d %d %d %d %d %d %d %d %d %d %d",
               lname, &rcvno, &drift, &year1, &mon1, &day1, &hour1, &min1,
                           &year2, &mon2, &day2, &hour2, &min2);
            if ( rcvno == RCVNO ) { // drift info found
                if (nfread != 13) {
                    warn("Missing or wrong number of parameters for clock drift info");
                    err("%s", line);
                } else {
                    DRIFT = drift;
                    break;
                }
            }
        } else if (line[0] == 'S') {  // already get to source line, stop
            break;
        }
    }
    if (DRIFT == 0) { // try to read from named parameters
        if (!getparint("drift", &drift)) {
            warn("Clock drift for station %d not given, assumed zero (=0)", RCVNO);
        } else {
            int* arri = ealloc1int(8);
            memset(arri, 0, 8*ISIZE);
            if ((nfread = countparval("td1")) >= 4) {
                getparint("td1", arri);
                year1 = arri[0]; mon1 = arri[1]; day1 = arri[2];
                hour1 = arri[3]; min1 = arri[4]; ;
            } else {
                err("Date/Time for clock sync must be given (td1=)");
            }
            if ((nfread = countparval("td2")) >= 4) {
                getparint("td2", arri);
                year2 = arri[0]; mon2 = arri[1]; day2 = arri[2];
                hour2 = arri[3]; min2 = arri[4]; ;
            } else {
                err("Date/Time for clock drift check must be given (td2=)");
            }
        }
        DRIFT = drift;
    }
    if (DRIFT == 0) {
        warn("Clock drift not given anywhere, assumed zero");
    } else {
        warn("Clock drift %d miliseconds for station %d between %d-%02d-%02d %02d:%02d ~ %d-%02d-%02d %02d:%02d",
            DRIFT, RCVNO, year1, mon1, day1, hour1, min1, year2, mon2, day2, hour2, min2);
    }

    // read nav file to find coordinate info of corresponding station
    double GX = 0.0, GY = 0.0, GZ = 0.0;
    double gx = 0.0, gy = 0.0, gz = 0.0;
    rewind(fpnav);
    memset(line, 0, 128);
    while ( (nbread = getline(&line, &ntoread, fpnav)) > 0 ) {
        if (line[0] == 'R') {  // receiver station line
            nfread = sscanf(line, "%s %d %lf %lf %lf", lname, &rcvno, &gx, &gy, &gz);
            if ( rcvno == RCVNO ) { // drift info found
                if (nfread < 4) {
                    warn("Missing parameters for receiver station info");
                    warn("%s", line);
                } else {
                    GX = gx;  GY = gy;
                    if (nfread > 4) GZ = gz;
                    break;
                }
            }
        } else if (line[0] == 'S') {  // already get to source line, stop
            break;
        }
    }
    if (GX == 0.0 && GY == 0.0) {
        if (verbose) warn("Receiver UTM coordinates not found in navigation file, trying to read from command line...");
        if (getpardouble("gx", &gx) && getpardouble("gy", &gy)) {
            GX = gx; GY = gy;
        } else {
            if (verbose) warn("Receiver UTM coordinates not specified, trying to retrieve from SAC file...");
            if (rcdhdr.nxsize == -12345 || rcdhdr.nysize == -12345) {
                err("Receiver UTM coordinates not specified anywhere, should be given in navigation file");
            } else {
                GX = rcdhdr.nxsize/DIST_SCALE;
                GY = rcdhdr.nysize/DIST_SCALE;
            }
        }
        if (GZ == 0) {
            if (verbose) warn("Water depth at receiver not found in navigation file, trying to read from command line...");
            if (getpardouble("gz", &gz)) {
                GZ = gz;
            } else {
                if (verbose) warn("Water depth at receiver not specified, trying to retrieve from SAC file...");
                if (rcdhdr.stdp == -12345.0) {
                    if (verbose) warn("Water depth at receiver not specified anywhere, should be given in navigation file");
                } else {
                    GZ = rcdhdr.stdp;
                }
            }
        }
    }

    if (verbose) {
        warn("Receiver station %d cordinates: X=%5.1f  Y=%5.1f  Z=%3.1f", RCVNO, GX, GY, GZ);
    }

    double trcv, tsrc, tdrift1, tdrift2;

    // construct starting date and time in julian day
    int year, mon, day, dayy, hour, min, sec, msec;
    year = rcdhdr.nzyear;
    dayy = rcdhdr.nzjday;
    hour = rcdhdr.nzhour;
    min  = rcdhdr.nzmin;
    sec  = rcdhdr.nzsec;
    msec = rcdhdr.nzmsec;
    trcv = julian_day(year, 1, 1, hour, min, sec, msec) + dayy - 1;

    if ( DRIFT) {
        // construct starting and ending date and time for clock drift
        tdrift1 = julian_day(year1, mon1, day1, hour1, min1, 0, 0);
        tdrift2 = julian_day(year2, mon2, day2, hour2, min2, 0, 0);
    }

    // set common trace header values
    memset(&tr, 0, HDRBYTES);  // reset output trace header

    tr.cdp = RCVNO;
    tr.trid = 1;
    tr.dt = (short) NINT(1000000*DT);
    tr.ns = NS;
    tr.scalco = tr.scalel = (short) ((DIST_SCALE <= 1.0)? 1/DIST_SCALE : -DIST_SCALE);
    tr.gx = NINT(DIST_SCALE*GX);
    tr.gy = NINT(DIST_SCALE*GY);
    tr.gwdep =  NINT(DIST_SCALE*GZ);
    tr.duse = rcdhdr.imagtyp;
    tr.cdpt = rcdhdr.imagtyp;
    tr.f1   = t0;


    int ntr = 0;  // global trace number counter

    // loop through nav file to find shot info
    //double SX = 0.0, SY = 0.0, SZ = 0.0;
    double x = 0.0, y = 0.0, z = 0.0;
    rewind(fpnav);
    memset(line, 0, 128);
    while ( (nbread = getline(&line, &ntoread, fpnav)) > 0 && !feof(fpnav)) {
        if (line[0] == 'S') {  // source point line
            nfread = sscanf(line, "%s %d %lf %lf %lf %d %d %d %d %d %d.%d",
                lname, &srcno, &x, &y, &z, &year, &mon, &day, &hour, &min, &sec, &msec);

            if (nfread != 12) {
                warn("Missing or wrong number of parameters forsshot point info");
                err("%s", line);
            }

            // calc distance and check if it is within range specified
            double offset = sqrt((GX - x)*(GX - x) + (GY - y)*(GY - y));
            if (ABS(offset) > 1000.0*maxd) {
                if (verbose > 9) warn("Shot %d out of range (%.3f > %.0f)", srcno, ABS(0.001*offset), maxd);
                continue;
            }
            // construct starting and ending date and time for clock drift
            tsrc = julian_day(year, mon, day, hour, min, sec, msec);

            // calc clock drift at source time
            double tdsrc = 0.0;
            if (DRIFT) {
                tdsrc = 0.001*((double) DRIFT)*(tsrc - tdrift1)/(tdrift2 - tdrift1);
            }
            // calc starting time and sample
            double tstart = 86400.0*(tsrc - trcv) - tdsrc + t0; //  in sec
            if (rv > 0.0) tstart += 0.001*offset/rv;


            int ISTART = NINT(tstart/DT);
            if ( (ISTART + NS) < 0 ) { // preceding of recording window
                if(verbose) warn("SHOT %d of %d-%02d-%02d %02d:%02d:%02d.%03d is preceding recording window",
                        srcno, year, mon, day, hour, min, sec, msec);
                continue;
            }
            if ( ISTART > NPTS - 1) { // out of recording window
                if(verbose) warn("SHOT %d of %d-%02d-%02d %02d:%02d:%02d.%03d is out of recording window. Terminating",
                        srcno, year, mon, day, hour, min, sec, msec);
                break;
            }

            ++ntr;

            if (verbose > 3) {
                if(verbose) warn("%4d-th: SHOT %d of distance %7.3f km of time %d-%02d-%02d %02d:%02d:%02d.%03d",
                        ntr, srcno, 0.001*offset, year, mon, day, hour, min, sec, msec);
            }

            int is1 = MAX(0, ISTART);
            int is2 = MIN(NPTS - 1, ISTART + NS);

            memset(tr.data, 0, NS*FSIZE);
            fseek(fpsac, is1*FSIZE  + sizeof(SACHdr), SEEK_SET);
            fread(&tr.data[ISTART >= 0 ? 0 : -ISTART], 1, (is2 - is1)*FSIZE, fpsac);

            if (endian) {
                for(i=0; i<NS; ++i) swap_float_4(&tr.data[i]);
            }

            // set trace header
            tr.tracf  = ntr;
            tr.fldr   = srcno;
            tr.sx     = NINT(DIST_SCALE*x);
            tr.sy     = NINT(DIST_SCALE*y);
            tr.swdep  = NINT(DIST_SCALE*z);
            tr.offset = NINT(DIST_SCALE*offset);

            tr.tracl = ntr;
            tr.tracr = ntr;  

            // save source time in header
            tr.year     = (short) year;
            tr.day      = (short) dayofyear(year, mon, day);
            tr.hour     = (short) hour;
            tr.minute   = (short) min;
            tr.sec      = (short) sec;
            tr.timbas   = (short) msec;  // milisecond

            /* Write the trace to disk */
            puttr(&tr);
        }
    }

    efclose(fpnav);
    efclose(fpsac);

    if (verbose) warn("Totally %d traces have been written", ntr);

    return (CWP_Exit());
}
