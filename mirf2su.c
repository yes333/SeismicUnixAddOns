/* Copyright (c) Read Well Services, 2012.*/
/* All rights reserved.                       */

/* MIRF2SU: $Revision: 1.0.0 $ ; $Date: 2012/10/27 03:03:03 $     */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>
#include <stdlib.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" MIRF2SU - read an Avalon MirF4 file, split into multiple traces of a given    ",
"           length with possible overlap, while map header info into trace      ",
"           header keys                                                         ",
"                                                                               ",
"   mirf2su file=   > stdout                                                    ",
"                                                                               ",
" Required parameter:                                                           ",
" file=            input MIRF record filename                                   ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
" length=8.1  [s]  length of trace.                                             ",
"                  =0 output as single trace provided it has up to 65535 samples",
"                                                                               ",
" overlap=0.1 [s]  time overlap between two consecutive traces                  ",
"                                                                               ",
" key=corr         key holding the ordering no. of splitt traces (1,2,3,...)    ",
"                                                                               ",
" hfile=header     file holding extra info about unmapped record/channel headers",
"                                                                               ",
" format=(header)  forcedly set format code (in case of header error)           ",
"                                                                               ",
" verbose=0        >0 echo info                                                 ",
"                                                                               ",
" Notes:                                                                        ",
"  every individual MIRF trace is cut to a length given by length=8.1 [s] as    ",
"  default and output consecutively with an ordering number (PanelNo) stored in ",
"  header word specified by key=. Time overlap between two consecutive traces   ",
"  can be given by overlap= (default 0.1 s).                                    ",
"  A new unique gather key is made up by multplying field record number by 10   ",
"  plus trace ordering number and store in keyword ep (ep=10*fldr + PanelNo).   ",
"  the output data must be sorted by cutted panel ep before further input for   ",
"  event detection.                                                             ",
"                                                                               ",
"  Unmapped header info in both record and channel headers are written to hfile ",
"  as ascii text.                                                               ",
"  Currently only format 2/5 (16bit Geochain IFP/24bit integer) are implemented.",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
"  mirf2su file=f_fldr-no.rcv verbose=1 hfile=f_fldr-no.hdr |\\",
"  spgsort key=fldr sort=ep |\\",
"  ......",
"                                                                               ",
"                                                                               ",
" Version 1.1.3   Last updated Nov, 2012 by Sanyu Ye                            ",
"                                                                               ",
NULL};

/*
 *
 * Credits:
 *	Sanyu Ye sanyu_ye@yahoo.com.com
 */
/**************** end self doc ***********************************/

typedef struct GeneralHeader // define general header structure
{
    int MIRF_version;
    int File_type;
    int Format_code;
    int Correlation_flag;
    int Controller_type;
    int Tool_system;
    int Channels_defined;
    int Test_mode;
    int Dataset_id;
    int Year;
    int Month;
    int Day;
    int Hour;
    int Minute;
    int Second;
    int Timezone_bias_seconds;
    int Source_id;
    int Number_of_receivers;
    int Receiver_number;
    int Number_in_stack;
    int Measurement_units;
    int Receiver_polarity;
    int Source_refenrence_channel;
    int reserved;
    int Record_number;
    int Stack_number;
    int Fix_number;
    int Tool_MD;
    int SIus;
    int Line_number;
    int Gun_pressure;
    int Software_version_x100;
    int SCX;
    int SCY;
    int TCX;
    int YCY;
    int WRE;
    int SRE;
    int SD;
    int S2M;
    int Source_elevation_error;
    int External_reference_delay_us;
    int Supplied_Ts_us;
    int TB_advance_us;
    int Tool_skew_us;
    int Raw_SCX_x10;
    int Raw_SCY_x10;
    int Controller_second;
    int Microsecond;
    int Timestamp_mode;
    int SDE;
    int dummy[76];
}  GenHdr;

typedef struct ChannelHeader
{
    int Legacy_Owner;
    int Descriptor;
    int Format_code;
    int NS;
    int Pointer_us;
    int Owner;
    int RCX;
    int RCY;
    int TVD;
    int MDO;
    int Reserved1;
    int Reserved2;
    float SSF;
    float DC;
    float SF;
    float Max_magnitude;
} ChnHdr;

#include "sprinthlp.c"

float convFormat(const int ns, const int format, const float SSF, const float DC,
        const float SF, char* buf, float* data, const int verbose);

float convFormat2(const int ns, const int format, const float SSF, const float DC,
        const float SF, char* buf, float* data, const int verbose);

int main(int argc, char **argv)
{
    cwp_String key;		/* header key word from segy.h		*/
    cwp_String type;            /* type of key				*/
    Value val;                  /* value of key			*/
    int index;
    char *file;                 /* name of seg2 file	*/
    char *hfile;                /* name of ascii header file	*/

    FILE *filefp = NULL;        /* file pointer for tape	*/
    FILE *headerfp;             /* file pointer for hfile	*/

    int i;      /* counter				*/
    int itr;    /* current trace number			*/
    int nread;
    int ns, NS = 0;     /* number of data samples		*/
    int dt;         /* sampling interval     		*/
    int verbose;    /* echo info ...			*/

    int format;     /* flag for to specify override format	*/
    float length, overlap;
    float* data = NULL;  // buffer for trace data

    char* buf = NULL;  // buffer for channel data
    int func;
    
    segy tr;
    GenHdr rcdhdr;
    ChnHdr chnhdr;

    /* Initialize */
    initargs(argc, argv);
    requestdoc(1); /* stdin not used */

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
    MUSTGETPARSTRING("file", &file);
    if (!getparstring("hfile", &hfile)) hfile = "header";
    if (!getparstring("key", &key))	  key = "corr";
    type = hdtype(key);
    index = getindex(key);

    /* Set parameters */
    if (!getparint("func", &func)) func = 0;
    if (!getparint("ns", &ns)) ns = 0;
    if (!getparint("dt", &dt)) dt = 0;
    if (!getparint("format", &format)) format = 0;
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparfloat("length", &length)) length = 8.1;
    if (!getparfloat("overlap", &overlap)) overlap = 0.1;


    /* Open input MIRF files */
    filefp = efopen(file, "r");
    nread = (int) fread((char *) &rcdhdr, 1, 512, filefp);  // read all pointer to trace blocks to buffer
    /* the header file containing header infor */
    headerfp = efopen(hfile, "w");
    if (verbose) warn("header file (%s) opened successfully", hfile);

    // decipher file header block and write to info header file
    int NCHN = rcdhdr.Channels_defined;
    int NRCV = rcdhdr.Number_of_receivers;

    fprintf(headerfp, "Record no: %d  Filename: %s  MIRF Version=%d  ",
        rcdhdr.Record_number, file, rcdhdr.MIRF_version);
    fprintf(headerfp, "Software_version: %4.2f\n", rcdhdr.Software_version_x100/100.0);
    fprintf(headerfp, "File Type: %s data  ", (rcdhdr.File_type == 0)? "raw" : "stacked");
    fprintf(headerfp, "%s\n", (rcdhdr.Correlation_flag == 0)? "uncorrelated" : "correlated");
    fprintf(headerfp, "Controller type: %d  Tool system: %d  Test mode: %d  Data ID: %d\n",
        rcdhdr.Controller_type, rcdhdr.Tool_system, rcdhdr.Test_mode, rcdhdr.Dataset_id);
    fprintf(headerfp, "(maximum) number of receivers =%d   ", NRCV);
    fprintf(headerfp, "Number of channels defined =%d\n", NCHN);
    fprintf(headerfp, "Timestamp mode=%d %s time  ", rcdhdr.Timestamp_mode,
       rcdhdr.Timestamp_mode == 0 ? "PC local" :  rcdhdr.Timestamp_mode == 1 ? "PC UTC" : "GPS UTC");
    fprintf(headerfp, "Timezone bias %d [s]\n", rcdhdr.Timezone_bias_seconds);
    fprintf(headerfp, "Starting Time: %4d-%02d-%02d %02d:%02d:%02d:%07.3f\n", rcdhdr.Year, rcdhdr.Month, rcdhdr.Day,
                       rcdhdr.Hour, rcdhdr.Minute, rcdhdr.Second, 0.001*rcdhdr.Microsecond);
    fprintf(headerfp, "Mesurement units: %s   ", rcdhdr.Measurement_units == 1 ? "millimeters" : "millifeet");
    fprintf(headerfp, "Receiver polarity: %d\n", rcdhdr.Receiver_polarity);
    fprintf(headerfp, "Source ID: %d  Source reference channel %d  Gun Pressure %d\n",
        rcdhdr.Source_id, rcdhdr.Source_refenrence_channel, rcdhdr.Gun_pressure);
    fprintf(headerfp, "Stack number: %d  Fix number %d\n", rcdhdr.Stack_number, rcdhdr.Fix_number);
    fprintf(headerfp, "Tool MD: %d  S2M: %d\n", rcdhdr.Tool_MD, rcdhdr.S2M);
    fprintf(headerfp, "External reference delay: %d  Supplied time correlation to seismic datum: %d\n",
        rcdhdr.External_reference_delay_us, rcdhdr.Supplied_Ts_us);
    fprintf(headerfp, "Delay TB to first sample: %d  Timing skew correction %d\n", rcdhdr.TB_advance_us, rcdhdr.Tool_skew_us);
    fprintf(headerfp, "Controller second: %d   Microsecond: %d\n", rcdhdr.Controller_second, rcdhdr.Microsecond);
    fprintf(headerfp, "Seismic data elevation: %d\n", rcdhdr.SDE);

    if ( rcdhdr.Format_code == -1) {
        fprintf(headerfp, "Variable data format specified by channel header\n");
        if (verbose) warn("Variable data format specified by channel header");
    } else {
        if (!format) format = rcdhdr.Format_code;
        else if (format != rcdhdr.Format_code) {
            if (verbose) warn("format (=%d) in record header is overriden by user input (=%d)",
                    rcdhdr.Format_code, format);
        }
    }
    /***
    // extract file number from record file name
    int digitFound = 0;
    int sufixFound = 0;
    int fileno = 0;
    for (i=0, file += strlen(file) - 2; i < 10; ++i, --file) {
        if ( *file == '.' ) { 
            *file = '\0';
            sufixFound = 1;
            continue;
        }
        if (sufixFound && isdigit(*file)) digitFound = 1;
        if (digitFound && !isdigit(*file)) break;
    }
    if(digitFound) fileno = atoi(++file);
    ****/

    // construct starting date and time
    struct tm tacq;
    tacq.tm_year = rcdhdr.Year - 1900;
    tacq.tm_mon  = rcdhdr.Month - 1;
    tacq.tm_mday = rcdhdr.Day;
    tacq.tm_hour = rcdhdr.Hour;
    tacq.tm_min  = rcdhdr.Minute;
    tacq.tm_sec  = rcdhdr.Second;
    //tacq.tm_zone = "UTC";
    //tacq.tm_gmtoff = 0;
    //tzset();
    setenv("TZ", "UTC", 1);
    time_t ttacq = mktime(&tacq);  // time in second
    double us = rcdhdr.Microsecond;

    // set common trace header values
    memset(&tr, 0, HDRBYTES);  // reset output trace header

    tr.fldr = rcdhdr.Record_number;
    tr.trid = 1;
    tr.scalco = tr.scalel = -1000;   // unit set to millimeter/feet
    tr.dt = (short) rcdhdr.SIus;
    tr.sx = rcdhdr.SCX;
    tr.sy = rcdhdr.SCY;
    tr.sdel = rcdhdr.SRE;
    tr.sdepth = rcdhdr.SD;
    tr.gdel = rcdhdr.WRE;

    int ctr = 0;  // global trace number counter
    int n_chn_missing = 0;
    int n_chn_live = 0;
    int nSize = 2;  // default number of bytes per sample
    int nOffsetData = 0; // offset of data block to the end of channel header

    /* main loop to read channel header and traces */
    for (itr = 0; itr < NCHN; ++itr) {
        fflush(stdout);
        fflush(headerfp);
        //read channel header
        fseek(filefp, 512 + itr * 64, SEEK_SET);
        nread = (int) fread((char *) &chnhdr, 1, 64, filefp);
        unsigned int ns_chn = chnhdr.NS;
        if (!NS && ns_chn > 0) NS = ns_chn;  // assign for first time
        if (ns_chn == 0) {
            if (verbose > 1) warn("no data sample for channel %d", itr + 1);
            fprintf(headerfp, "no data sample for channel %d\n", itr + 1);
            ++n_chn_missing;
            continue;
        }

        if ( ns_chn != NS) {
            err("inconsistent number of samples ns (%d != %d)\n", ns, NS);
        }

        // print info to header file
        fprintf(headerfp, "\nHeader info for channel %d:\n", itr + 1);
        format = chnhdr.Format_code;
        if (chnhdr.Format_code == 2) {
            fprintf(headerfp, "Format code 2, Geochain IFP 16 bit\n");
            nSize = 2;
        } else if (chnhdr.Format_code == 5) {
            fprintf(headerfp, "Format code 5, 24 bit integer\n");
            nSize = 3;
        } else {
            fprintf(headerfp, "Data format (=%d) is not implemented, abort", chnhdr.Format_code);
            efclose(filefp);
            efclose(headerfp);
            err("Data format (=%d) dis not implemented, abort", chnhdr.Format_code);
        }
        // print info to header file
        fprintf(headerfp, "Owner/Receiver=%d   Descriptor/Component=%d   Cursor position=%d [us]\n",
            chnhdr.Owner, chnhdr.Descriptor, chnhdr.Pointer_us);
        fprintf(headerfp, "SSF=%10.6g   SF=%10.6g  DC=%10.6g   maxAmpl=%10.6g\n",
            chnhdr.SF, chnhdr.SSF, chnhdr.DC, chnhdr.Max_magnitude);

        // set trace header
        tr.tracf = itr + 1;
        tr.cdp = chnhdr.Owner;
        tr.cdpt = chnhdr.Descriptor;
        tr.duse = n_chn_live % 3 + 1;
        tr.gx = chnhdr.RCX;
        tr.gy = chnhdr.RCY;
        tr.gelev = chnhdr.TVD;
        tr.offset = rcdhdr.Tool_MD + chnhdr.MDO;
        tr.unscale = chnhdr.SF;
        tr.ungpow = chnhdr.SSF;


        // read datata trace
        ++n_chn_live;

        fseek(filefp, 512 + NCHN * 64 + nOffsetData, SEEK_SET);
        if (!buf) buf = (char*) ealloc1(NS * 4, 1);  // big enough for all formats
        if (!data) data = ealloc1float(NS);
        nread = (int) fread((char *) buf, 1, ns_chn*nSize, filefp);
        nOffsetData += ns_chn * nSize;  // add offset for next file
        float medianV = convFormat(ns_chn, format, chnhdr.SSF, chnhdr.DC, chnhdr.SF, buf, data, verbose);
        if (verbose > 1) {
            fprintf(headerfp, "Median value (=%10.6g)\n", medianV);
            warn("Channel %d: median value (=%10.6g)", itr + 1, medianV);
        }
        int nt = NINT(1000000.0*length / tr.dt);  // sample number specified by user
        int no = NINT(1000000.0*overlap / tr.dt);  // number of overlapping samples of next trace
        
        if (NS > SU_NFLTS && nt == 0) {
            err("Number of samples (NS=%d) too big to be accommodated in single segy trace. Set length= to split trace", NS);
        }
        if (nt > SU_NFLTS) {
            err("Number of samples (tr.ns=%d) too big to be accommodated in single segy trace. Set length= to smaller value", nt);
        }

        if (nt == 0) { // output as single trace
            nt = NS;
            no = 0;
        }

        int ntr = 1;  // number of SU traces the original SEG2 trace are splitted into
        if (nt > 0) {
            ntr = NS / (nt - no);
            if ((int) NS - (ntr - 1)*(nt - no) - nt > 0) ++ntr;
        } else {
            nt = (int) NS;
        }
        tr.ns = (unsigned short) nt;

        for(i=0; i<ntr; ++i) { // set headers and copy data
            int it = i*(nt - no);
            int nsize = (NS > it + nt)? nt : NS - it + 1;
            memset(tr.data, 0, SU_NFLTS*FSIZE);
            memcpy(tr.data, &data[it], nsize*FSIZE);
            tr.tracl = ++ctr;
            tr.tracr = ctr;  // no. of splitted trace
            tr.ep = 10*tr.fldr + i +  1;
            setval(type, &val, (double) (i + 1));  // no. of splitted trace
            puthval(&tr, index, &val);

            // calc time of trace start
            int t_usec = (it*tr.dt + (int) us) % 1000000;
            int t_sec  = (it*tr.dt + (int) us) / 1000000;
            time_t ttstart = ttacq + t_sec;
            //struct tm* tstart = localtime(&ttstart);
            struct tm* tstart = gmtime(&ttstart);
            tr.year     = (short) tstart->tm_year + 1900;
            tr.day      = (short) tstart->tm_yday;
            tr.hour     = (short) tstart->tm_hour;
            tr.minute   = (short) tstart->tm_min;
            tr.sec      = (short) tstart->tm_sec;
            tr.timbas   = (short) ((int) t_usec / 1000);  // milisecond
            tr.trwf     = (short) ((int) t_usec % 1000);  // microsecond
            /* Write the trace to disk */
            puttr(&tr);
        }
    }

    if (n_chn_live != 3 * NRCV) {
        warn("No of live channels not in accord with no of receivers (%d != 3 x %d)",n_chn_live, NRCV);
    }
    efclose(filefp);
    efclose(headerfp);
    if (verbose) warn("Totally %d traces have been written", ctr);

    return (CWP_Exit());
}

float convFormat(const int ns, const int format, const float SSF, const float DC, 
        const float SF, char* buf, float* data, const int verbose)
{
    int i, j;
    short *sbuf;
    short exponent;
    int idata;
    float factor;
    char* p;
    
    memset(data, 0, ns * FSIZE);

    float minV = FLT_MAX;
    int lastValidSample = 0, InvalidSampleFound = 0;
    if (format == 2) {
        for (i = 0, sbuf = (short*) buf; i < ns; ++i, ++sbuf) {
            idata = (int) *sbuf;
            if (!InvalidSampleFound) lastValidSample = i - 1;
            if (idata == 0 || idata == 1 || idata == 2) {
                InvalidSampleFound = 1;
                if (verbose > 10) warn("sampling error occurs at %d sample with value=%d", i+1, idata);
                continue;
            }
            exponent = idata & 0x0003;
            idata = idata >> 2;
            switch (exponent) {
                case 3:
                    factor = 64.0;
                    break;
                case 2:
                    factor = 16.0;
                    break;
                case 1:
                    factor = 4.0;
                    break;
                default:
                    factor = 1.0;
                    break;
            }
                    
            data[i] = (float) idata / factor;
            if (data[i] < minV) minV = data[i];
            data[i] = (SF * data[i] - DC) * SSF;
            if (InvalidSampleFound) {
                // linear interpolate
                for (j =lastValidSample + 1; j < i; ++j) {
                    data[j] = (float) (j - lastValidSample) / (float) (i - lastValidSample)
                            * (data[i] - data[lastValidSample]) + data[lastValidSample];
                }
                InvalidSampleFound = 0;
            }
        }
    } else if (format == 5) {        
        for (i = 0, p = (char*) buf; i < ns; ++i) {
            idata = 0;
            //memcpy(&idata, p + 3*i - 1, 4); // assume little endian
            ///pi = &idata;
            ///memcpy(pi+1, p + 3*i, 3); // assume little endian
            //idata = idata >> 8;
            //pidata = &idata;
            memcpy(&idata, p + 3*i, 3); // assume big endian
            swap_int_4(&idata);
            //idata = idata << 8;
            //pidata = &uidata;
            data[i] = (float) (idata / 256);
            data[i] = (SF * data[i] - DC) * SSF;
        }
    } else {
        warn("data format (=%d) not implemented yet", format);
        return 0.0;
    }

    // remove DC
    //double dsum = 0;
    //for (i = 0; i < ns; ++i) dsum += data[i];
    //float avg = dsum / ns;
    memcpy(buf, data, ns * FSIZE);
    float median = quick_select((float*) buf, ns);
    for (i = 0; i < ns; ++i) {
        data[i] -=  median;
    }

    return median;
}
