/* Copyright (c) Read Well Services, 2011.*/
/* All rights reserved.                       */

/* SEG22SU: $Revision: 1.0.0 $ ; $Date: 2011/04/27 03:03:03 $     */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <time.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SEG22SU - read an SEG-2 file while map descriptor info into trace header keys  ",
"                                                                               ",
"   seg22su file=   > stdout                                                     ",
"                                                                               ",
" Required parameter:                                                           ",
" file=            input seg-2 filename                                         ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
" map=             =key,DESCRIPTOR_NAME,postion,scale                           ",
"                  map scaled nummerical value at position after the descriptior",
"                  into su trace header key. Multiple entries permitted.        ",
"                  if position and scale are omitted, default set to 1 and 1.0  ",
"                                                                               ",
" hfile=header     ascii file holding printable text of file/trace descriptor blocks",
"                                                                               ",
" verbose=0        >0 echo info                                                 ",
"                                                                               ",
" format=(header)  convert to native floating point according to format value   ",
"                    if not set (=0), assume =4 (IEEE floating point)           ",
" ns=(header)      number of samples (use if trace header ns is not set)        ",
" dt=        [us]  sampling interval (microsecond)                              ",
"                    use if failed to extract from descriptor                   ",
"                                                                               ",
" Notes:                                                                        ",
"  Descriptions in both file and trace blocks are dumped to hfile as ascii text.",
"  User should check hfile to remap/override and add missing parameters.        ",
"  ACQUISITION_DATE dd/mm/yyyy and ACQUISITION_TIME hh:mm:ss are mapped as usual.",
"  Following descriptor names are used to figure out the fraction of second:    ",
"  ACQUISITION_TIME_US, ACQUISITION_TIME_MICROSECONDS, ACQUISITION_TIME_NANOSECONDS.",
"  Kerwords timbas & trwf are used to store mili- and microseconds.             ",
"                                                                               ",
" Examples:                                                                     ",
"                                                                               ",
" seg22su file=my.sg2 map=ungpow,DESCALING_FACTOR map=gaps,SKEW,1,1000.0 >my.su ",
"                                                                               ",
" More built-in/default mapping:                                                ",
"                                                                               ",
"   map=dt,SAMPLE_INTERVAL,1,1000000.0",
"   map=sx,SOURCE_LOCATION,1,100.0",
"   map=sy,SOURCE_LOCATION,2,100.0",
"   map=sdepth,SOURCE_LOCATION,3,100.0",
"   map=gx,RECEIVER_LOCATION,1,100.0",
"   map=gy,RECEIVER_LOCATION,2,100.0",
"   map=gelev,RECEIVER_LOCATION,3,100.0",
"   map=gdel,DATUM,1,100.0",
"   map=tracf,CHANNEL_NUMBER,1,1.0",
"   map=cdp,RECEIVER_STATION_NUMBER,1,1.0",
"   map=fldr,SHOT_SEQUENCE_NUMBER,1,1.0",
"   map=ep,SOURCE_STATION_NUMBER,1,1.0",
"   map=unscale,DESCALING_FACTOR,1,1.0",
"                                                                               ",
" Format Code:                                                                  ",
" 1: 16 bit fix point                                                           ",
" 2: 32 bit fix point                                                           ",
" 3: 20 bit SEGD floating point, not supported                                  ",
" 4: 32 bit IEEE floating point                                                 ",
" 5: 64 bit IEEE floating point, not supported                                  ",
"                                                                               ",
" Version 1.1.0   Last updated April, 2011 by Sanyu Ye                          ",
"                                                                               ",
NULL};

/*
 *
 * Credits:
 *	Sanyu Ye sanyu.ye@readgroup.com
 *      Read Well Services, Oslo
 */
/**************** end self doc ***********************************/

// forward declaration

int setHdrKey(const char* key, const char* name, const int pos, const float scale, const char* strbuf, segy* ptr);
int setAddKeys(int nkeys, char** keys, char** names, int* pos, float* scales, const char* strbuf, segy* ptr);
char* getSubstr(const char* strbuf, const char* delim, const char* name, const int pos, int* nc);
double getValue(const char* strbuf, const char* name, const int pos, const float scale);
static void setval(cwp_String type, Value *valp, double dval);

const double Invalid_Value = -333.333;


int setDefKeys(const char* strbuf, segy* ptr)
{
    int nFailed = 0;

    nFailed += setHdrKey("dt",     "SAMPLE_INTERVAL",          1,  1000000.0, strbuf, ptr);
    nFailed += setHdrKey("sx",     "SOURCE_LOCATION",          1,      100.0, strbuf, ptr);
    nFailed += setHdrKey("sy",     "SOURCE_LOCATION",          2,      100.0, strbuf, ptr);
    nFailed += setHdrKey("sdepth", "SOURCE_LOCATION",          3,      100.0, strbuf, ptr);
    nFailed += setHdrKey("gx",     "RECEIVER_LOCATION",        1,      100.0, strbuf, ptr);
    nFailed += setHdrKey("gy",     "RECEIVER_LOCATION",        2,      100.0, strbuf, ptr);
    nFailed += setHdrKey("gelev",  "RECEIVER_LOCATION",        3,      100.0, strbuf, ptr);
    nFailed += setHdrKey("gdel",   "DATUM",                    1,      100.0, strbuf, ptr);
    nFailed += setHdrKey("tracf",  "CHANNEL_NUMBER",           1,        1.0, strbuf, ptr);
    nFailed += setHdrKey("cdp",    "RECEIVER_STATION_NUMBER",  1,        1.0, strbuf, ptr);
    nFailed += setHdrKey("fldr",   "SHOT_SEQUENCE_NUMBER",     1,        1.0, strbuf, ptr);
    nFailed += setHdrKey("ep",     "SOURCE_STATION_NUMBER",    1,        1.0, strbuf, ptr);
    nFailed += setHdrKey("unscale","DESCALING_FACTOR",         1,        1.0, strbuf, ptr);

    return nFailed;
}

int main(int argc, char **argv)
{
    char *file;     /* name of seg2 file	*/
    char *hfile;    /* name of ascii header file	*/

    FILE *filefp = NULL; /* file pointer for tape	*/
    FILE *headerfp;      /* file pointer for hfile	*/

    int i, j;   /* counter				*/
    int itr;    /* current trace number			*/
    //int trmin;  /* first trace to read			*/
    //int trmax;  /* last trace to read			*/
    int ns;     /* number of data samples		*/
    int dt;     /* sampling interval     		*/
    int verbose; /* echo info ...			*/

    int format; /* flag for to specify override format	*/
    int nkeys, pos[32];
    float scales[32];
    cwp_String pars[4];

    unsigned short usbuf;
    short           sbuf;
    unsigned int   uibuf;
    int             ibuf;
    float           fbuf;

    char buf[SU_NFLTS];     // general buffer
    char strbuf[SU_NFLTS];  // string buffer holding description text for parsing

    segy tr;

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

    /* Set parameters */
    if (!getparint("ns", &ns)) ns = 0;
    if (!getparint("dt", &dt)) dt = 0;
    if (!getparint("format", &format)) format = 0;
    //if (!getparint("trmin", &trmin)) trmin = 1;
    //if (!getparint("trmax", &trmax)) trmax = INT_MAX;
    if (!getparint("verbose", &verbose)) verbose = 0;

    nkeys = countparname("map");
    cwp_String keys[nkeys];
    cwp_String names[nkeys];
    for (i=0; i<nkeys; ++i) {
        int npar = countnparval(i, "map");
        if (npar < 2 || npar > 4) err("%-dth map=%s must have at least 2 but at most 4 parameters KEY,NAME", i+1, pars[0]);
        getnparstringarray(i, "map", pars);
        keys[i] = pars[0];
        names[i] = pars[1];
        pos[i]    = (npar > 2)? eatoi(pars[2]) : 1;
        scales[i] = (npar > 3)? eatof(pars[3]) : 1.0;
        if (verbose>1) {
            warn("%d-th mapping: key='%s'  name='%s'  pos=%d  scale=%10.4e", i+1, keys[i], names[i], pos[i], scales[i]);
        }
        if (pos[i] < 1 || pos[i] > 3 ) {
            warn(" Possibly wrong position number (=%d, usually 1<=pos<=3) for %d-th mapping", pos[i], i+1);
        }
    }

    /* Open files */
    filefp = efopen(file, "r");

    // read file block
    int nread = (int) fread((char *) buf, 1, 32, filefp);  // read first 32 bytes
    memcpy(&usbuf, &buf[0], 2);
    if (verbose>1) warn("File Descriptor Block ID = 0x%2x%2x", buf[0], buf[1]);
    if (usbuf != 0x3a55 && usbuf != 0x553a) err("Invalid File Descriptor Block ID");
    unsigned short ID = usbuf;
    int native = ( ID == 0x3a55 ) ? 1 : 0;   // if file is native, no byte swapping needed
    if (verbose) warn("SEG-2 File %s %s in native format", file, (native)? "is" : "is not");
    memcpy(&usbuf, &buf[2], 2);
    if (!native) swap_u_short_2(&usbuf);
    unsigned short RevNo = usbuf;
    memcpy(&usbuf, &buf[4], 2);
    if (!native) swap_u_short_2(&usbuf);
    unsigned short M = usbuf;   // size of Trace Pointer Sub-block
    memcpy(&usbuf, &buf[6], 2);
    if (!native) swap_u_short_2(&usbuf);
    unsigned short N = usbuf;   // number of Traces
    if (verbose)   warn("RevisionNo=%d  Number of traces=%d   Size of Trace Pointer Sub-block=%d", RevNo, N, M);
    if ( M < 4 * N) err("Inconsistent Number of traces (%d != %d) Size of Trace Pointer Sub-block=%d", N, M/4);

    // read all pointer to trace blocks
    nread = (int) fread((char *) buf, 1, M, filefp);  // read all pointer to trace blocks to buffer
    int* ptrb = ealloc1int(N);
    memset(ptrb, 0, N*ISIZE);
    for (i=0; i<N; ++i) {
        memcpy(&uibuf, &buf[i * ISIZE], 4);
        if (!native) swap_u_int_4(&uibuf);
        ptrb[i] = (int) uibuf;
    }

    /* the header file in ascii containing all printable text */
    headerfp = efopen(hfile, "w");
    if (verbose) warn("header file (%s) opened successfully", hfile);

    // read string block
    nread = ptrb[0] - M - 32;
    nread = (int) fread((char *) buf, 1, nread, filefp);
    int nc = 0, linefeed = 0;
    for (j=0, i=0; i<nread; ++i) {  // dump printable of block to header file
        if ( (isprint(buf[i]) || isspace(buf[i])) && !(iscntrl(buf[i])) ) {
            strbuf[j++] = buf[i];
            fputc(buf[i], headerfp);
            linefeed = 1;
            ++nc;
        } else if (linefeed) {  // insert space or line break at first non-printable character
            strbuf[j++] = ' ';  // space for simple inline parsing
            fputc('\n', headerfp);  // linefeed for easy grep / awk by unix utilities
            linefeed = 0;
        }
    }
    strbuf[j] = '\0';  // end string buffer

    // parse header descriptor block to find starting date and time
    int ndate, ntime;
    char* cs = getSubstr(strbuf, " ", "ACQUISITION_DATE", 1, &ndate);
    for(i=0; i<ndate; ++i) buf[i] = cs[i];
    buf[ndate] = ' ';
    cs = getSubstr(strbuf, " ", "ACQUISITION_TIME", 1, &ntime);
    for(i=0; i<ntime; ++i) buf[ndate + 1 + i] = cs[i];
    buf[ndate + ntime + 2] = '\0';
    struct tm tacq;
    strptime(buf, "%d/%m/%Y%t%T", &tacq);  // convert date/time to unix tm
    double us = getValue(strbuf, "ACQUISITION_TIME_US", 1, 1.0);
    if (us == Invalid_Value) us = getValue(strbuf, "ACQUISITION_TIME_MICROSECONDS", 1, 1.0);
    if (us == Invalid_Value) us = getValue(strbuf, "ACQUISITION_TIME_NANOSECONDS",  1, 0.001);

    memset(&tr, 0, HDRBYTES);  // reset output trace header
    tr.year     = (short) tacq.tm_year + 1900;
    tr.day      = (short) tacq.tm_yday;
    tr.hour     = (short) tacq.tm_hour;
    tr.minute   = (short) tacq.tm_min;
    tr.sec      = (short) tacq.tm_sec;
    tr.timbas   = (short) ((int) us / 1000);  // milisecond
    tr.trwf     = (short) ((int) us % 1000);  // microsecond

    tr.trid = 1;
    tr.scalco = tr.scalel = -100;   // unit set to cm
    // try to figure out header parameters if present from File Descriptor Block
    setDefKeys(strbuf, &tr);
    setAddKeys(nkeys, keys, names, pos, scales, strbuf, &tr);

    /* main loop to read traces */
    int nSizeTrace = ptrb[1] - ptrb[0];
    for (itr = 0; itr < N; ++itr) {
        fflush(stdout);
        fflush(headerfp);
        fseek(filefp, ptrb[itr], SEEK_SET);
        nread = (int) fread((char *) buf, 1, 32, filefp);
        memcpy(&usbuf, &buf[0], 2);
        if (!native) swap_u_short_2(&usbuf);
        ID = usbuf;
        if (ID != 0x4422) {
            warn("Invalid trace descriptor block ID=0x%4x (!=0x4422). Trace skipped", ID);
            continue;  // skip
        }
        memcpy(&usbuf, &buf[2], 2);
        if (!native) swap_u_short_2(&usbuf);
        unsigned short X = usbuf;
        memcpy(&uibuf, &buf[4], 4);
        if (!native) swap_u_int_4(&uibuf);
        unsigned int Y = uibuf;
        memcpy(&uibuf, &buf[8], 4);
        if (!native) swap_u_int_4(&uibuf);
        unsigned int NS = uibuf;
        if (NS == 0) {
            if (ns == 0) err("Number of samples not found in data, must be set by ns=");
            else NS = ns;
        }
        if (buf[12]) format = (int) buf[12];
        if (format == 3 || format ==5)
            err("format=%d, %d bit floating point not supported", format, (format == 3)? 20 : 64);
        int nSampleSize = 4;
        if (format == 1) nSampleSize = 2;  // 16 bit fix point
        if (verbose && !itr)
            warn("%d Data samples of format=%d --- %d bit %s point", NS, format, 8*nSampleSize, (format==4)? "floating" : "fix");
        if (itr < N - 1) nSizeTrace =  ptrb[itr+1] - ptrb[itr];
        if (nSizeTrace != X + Y)
            warn("Inconsistent trace block size ( %d != %d + %d ) for %d-th trace", nSizeTrace, X, Y, itr + 1);

        // read and dump string text to header file
        nread = (int) fread((char *) buf, 1, X - 32, filefp);
        fprintf(headerfp, "===%d===Description for Trace %d:\n", itr + 1, itr + 1);
        linefeed = 0;
        for (j=0, i=0; i<(X - 32); ++i) {
            if ( (isprint(buf[i]) || isspace(buf[i])) && !(iscntrl(buf[i])) ) {
                strbuf[j++] = buf[i];
                fputc(buf[i], headerfp);
                linefeed = 1;
                ++nc;
            } else if (linefeed) {
                strbuf[j++] = ' ';
                fputc('\n', headerfp);
                linefeed = 0;
            }
        }
        strbuf[j] = '\0';

        // find out NAME-VALUE pair in text buffer and set related values to header keys
        setDefKeys(strbuf, &tr);
        setAddKeys(nkeys, keys, names, pos, scales, strbuf, &tr);

        ns = MIN(NS, SU_NFLTS);
        tr.ns = ns;
        tr.tracl = tr.tracr = itr + 1;
        if (tr.dt <= 0  &&  dt > 0) tr.dt = dt;
        // read trace data block
        for (i=0; i<ns; ++i) {
            if (format == 1) {
                fread(&sbuf, nSampleSize, 1, filefp);
                if(!native) swap_short_2(&sbuf);
                tr.data[i] = (float) sbuf;
            }
            if (format == 2) {
                fread(&ibuf, nSampleSize, 1, filefp);
                if(!native) swap_int_4(&ibuf);
                tr.data[i] = (float) ibuf;
            }
            else { // 32 bit IEEE floating
                fread(&fbuf, nSampleSize, 1, filefp);
                if(!native) swap_float_4(&fbuf);
                tr.data[i] = fbuf;
            }
        }

        /* Write the trace to disk */
        puttr(&tr);
    }
    efclose(headerfp);
    if (verbose>1) warn("header file written with %d characters and closed successfully", nc);
    if (verbose) warn("Totally %d traces have been written", itr);


    return (CWP_Exit());
}

int setAddKeys(int nkeys, char** keys, char** names, int* pos, float* scales, const char* strbuf, segy* ptr)
{
    int i, nFailed = 0;

    for(i=0; i<nkeys; ++i) {
        nFailed += setHdrKey(keys[i], names[i], pos[i], scales[i], strbuf, ptr);
    }

    return nFailed;
}

int setHdrKey(const char* key, const char* name, const int pos, const float scale, const char* strbuf, segy* ptr)
{
    int index = getindex(key);
    cwp_String type = hdtype(key);
    double d = getValue(strbuf, name, pos, scale);
    if (d == Invalid_Value) return -1;
    Value val;
    setval(type, &val, d);
    puthval(ptr, index, &val);
    return 0;
}

char* getSubstr(const char* strbuf, const char* delim, const char* name, const int pos, int* nc)
{
    int i;
    const char* numset = "-.0123456789";

    char* cs = strstr(strbuf, name); // find name
    if (cs == NULL) {
        *nc = 0;
        return NULL;
    }

    for(i=0; i<pos; ++i) {
        cs = strpbrk(cs, delim);  // skip to end of name
        cs = strpbrk(cs, numset); // start of numbers
        *nc = strcspn(cs, delim);  // length of number string
    }
    
    return cs;
}

double getValue(const char* strbuf, const char* name, const int pos, const float scale)
{
    int nc;
    const char* delim = " ,/|#\t\n";

    char* cs = getSubstr(strbuf, delim, name, pos, &nc);
    
    if (!cs) return Invalid_Value;

    double d = eatod(cs);
    return d*scale;
}

static void setval(cwp_String type, Value *valp, double dval) 
{
    switch (*type) {
        case 's':
            err("can't set char header word");
            break;
        case 'h':
            valp->h = (short) dval;
            break;
        case 'u':
            valp->u = (unsigned short) dval;
            break;
        case 'l':
            valp->l = (long) dval;
            break;
        case 'v':
            valp->v = (unsigned long) dval;
            break;
        case 'i':
            valp->i = (int) dval;
            break;
        case 'p':
            valp->p = (unsigned int) dval;
            break;
        case 'f':
            valp->f = (float) dval;
            break;
        case 'd':
            valp->d = dval;
        default:
            err("unknown type %s", type);
            break;
    }
    return;
}
