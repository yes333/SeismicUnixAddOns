/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                      */

/* SUSETGATHER: $Revision: 1.0 $ ; $Date: 2008/02/07 $        */

#include <su.h>
#include <segy.h>
#include <header.h>
#include <segyhdr.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SUSETGATHER- set gather parameters upon key value change                      ",
"                                                                               ",
"   susetgather <infile >outfile  [optional parameters]                         ",
"                                                                               ",
" Optional Parameters:                                                          ",
"   key=cdpt        gather sorting key                                          ",
"                                                                               ",
"   xkey=           inline sorting key                                          ",
"                   if not specified, gather is treated as 2D and               ",
"                   trace is numbered sequencially                              ",
"   dir=0           direction of xkey value                                     ",
"                   =0 determined automatically (=1 if xkey[1] => xkey[0])      ",
"                   =1  xkey value increasing                                   ",
"                   =-1 xkey value decreasing                                   ",
"                                                                               ",
"   nmax=1024       max. number of traces of gathers                            ",
"                                                                               ",
"   transpose=0     =0 output in orignal inline then xline order                ",
"                   =1 output transposed first xline then inline                ",
"                                                                               ",
"   verbose=0       >0 echo information                                         ",
"                                                                               ",
" Notes:                                                                        ",
"   it counts inline (along X) and crossline (along Y) numbers of 2/3D gather   ",
"   with common 'key' value. The number of inline is determined by abrupt opposite",
"   change of 'xkey' value (e.g. from monotonously increase to sudden decrease or",
"   vice versa). The X and Y dimensions are stored in keywords tr.sfs and tr.sfe",
"   respectively, while inline and xline numbers (starting from 1) in tr.nhs and",
"   tr.nvs.                                                                     ",
"                                                                               ",
" Version 1.2 last modified Dec. 2011 by Sanyu Ye                               ",
NULL};

/*
* Credits: READ Well Service, Sanyu Ye, Feb. 2008.
*
*/
/**************** end self doc ********************************/

/** main **/
int main(int argc, char **argv)
{
    cwp_String key, xkey;    /* header key word from segy.h */
    cwp_String type, xtype;   /* type of key                 */
    int index, xindex;         /* index of key                */
    Value val, valnew, xval; /* value of key                */
    int nsegy;          /* Bytes read in for a trace   */
    int ns;             /* number of time samples of trace		*/
    int verbose;        /* flag for echoing information */
    int ntr, ntotal;    /* trace counters */
    int ngather;        /* gather counters */
    int nmax, nxmax;          /* max. number of traces expected in gather */
    int nx, ny;         /* actual dimensions of gather */
    int transpose;      /* flag	to output traces (grid) transposed */
    int eof = 0;        /* END OF File flag			*/
    int i, ix, dir;
    int dim = 2;            /* default 2D gather without xkey set */
    double dval, dvalnew;   /* double values for xkey */
    segy tr, outtr, *ptr;
    segyhdr *hdrs = NULL; /* placeholder for incoming trace headers	*/
    float **trdata = NULL; /* temporary array for trace data		*/

    /* hook up getpar to handle the parameters */
    initargs(argc,argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose= 0;
    if (!getparint("transpose", &transpose)) transpose = 0;
    if (!getparint("nmax", &nmax)) nmax = 1024;

    /* Get info from first trace */
    if (!(nsegy = gettr(&tr))) err("can't read first trace");
    ns = (int) tr.ns;

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdpt";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);
    if (getparstring("xkey", &xkey)) {
        dim = 3;
        xtype = hdtype(xkey);
        xindex = getindex(xkey);
    }
    if (!getparint("dir", &dir)) dir = 0;

    // allocate and reset memory for input and output traces
    hdrs = (segyhdr*) ealloc1(nmax, HDRBYTES);
    trdata = ealloc2float(ns, nmax);
    memset(hdrs, 0, nmax * HDRBYTES);
    memset(*trdata, 0, nmax * ns * FSIZE);

    /* Read headers and data while getting a count */
    ngather = 0;
    ntr = 0;
    ntotal = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > nmax - 1) err("\nNumber of traces exceeding nmax=%d\n", nmax);
            memcpy(&hdrs[ntr], &tr, HDRBYTES);
            memcpy(trdata[ntr], tr.data, FSIZE * ns);
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            ++ngather;
            if (verbose > 1) warn("  processing %d traces of %d-th gather with %s=%d",
                                     ntr, ngather, key, vtoi(type, val));

            ntotal += ntr;

            nxmax = nx = 0;
            ny = 1;
            for (i=0; i<ntr; ++i) {  // loop over gather to check xkey change
                if (ntr > 1 && dim == 3) {
                    ptr = (segy*) &hdrs[i];
                    gethval(ptr, xindex, &xval);
                    dvalnew = vtod(xtype, xval);
                    if (i == 0) {
                        dval = dvalnew;
                    } else if (i == 1) { // compare xkey value and set direction at second input trace
                        if (dir == 0) dir = (dvalnew >= dval)? 1 : -1;
                        else if (dir > 0 && dvalnew < dval) {
                            err("xkey direction conflict: xkey value specified increasing, but dereasing in data");
                        }
                        else if (dir < 0 && dvalnew > dval) {
                            err("xkey direction conflict: xkey value specified decreasing, but inreasing in data");
                        }
                    } 
                    if (i + 1 == ntr) { // last trace of gather
                        if (nxmax == 0) nxmax = nx + 1;
                        if (nx + 1 > nxmax) {
                            if (verbose) warn(" New max. inline number (%d > %d) on %d xline of %d gather with %s=%d",
                                nx, nxmax, ny, ngather, key, vtoi(type, val));
                            nxmax = nx + 1;
                        }
                    }
                    // compare xkey value and check change of direction or last trace in gather
                    else if ( (dir < 0 && dvalnew > dval) || (dir > 0 && dvalnew < dval) ) {
                        if (nxmax == 0) nxmax = nx;
                        if (nx > nxmax) {
                            if (verbose) warn(" New max. inline number (%d > %d) on %d xline of %d gather with %s=%d",
                                nx, nxmax, ny, ngather, key, vtoi(type, val));
                            nxmax = nx;
                        }
                        ++ny;
                        nx = 0;  // reset nx
                    }
                    ++nx;
                    hdrs[i].nhs = nx;
                    hdrs[i].nvs = ny;
                    dval = dvalnew;
                } else {  // 2D
                    hdrs[i].nhs = i + 1;
                    hdrs[i].nvs = 1;
                    if (nxmax == 0) nxmax = ntr;
                }
            }

            if (ntr == 1 && verbose > 1) { // single trace gather, warn
                 warn("  single trace gather for %d gather with %s=%d", ngather, key, vtoi(type, val));
            }
            
            if ( ntr != nxmax*ny) {
                 warn("  non-multiple number of traces (%d != %d x %d) on %d gather with %s=%d",
                                    ntr, nxmax, ny, ngather, key, vtoi(type, val));
                 if (transpose) err(" cannot output traces transposed");
            }
            for (i=0; i<ntr; ++i) {
                ix = (transpose)? (i%ny)*nxmax + i/ny : i;
                memcpy(&outtr, &hdrs[ix], HDRBYTES);
                memcpy(outtr.data, trdata[ix], FSIZE * ns);
                //outtr.ntr = ntr;  // disable because it cause suximage to fail
                outtr.shortpad = (transpose)? ny : nxmax;
                outtr.sfs = (transpose)? ny : nxmax;
                outtr.sfe = (transpose)? nxmax : ny;
                if (transpose) {
                    outtr.nhs = i%ny + 1;
                    outtr.nvs = i/nxmax + 1;
                }                
                puttr(&outtr);
            }

            val = valnew;
            // reset output data
            memset(hdrs, 0, nmax * HDRBYTES);
            memset(*trdata, 0, nmax * ns * FSIZE);
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose > 0) warn("  Totally %d gathers with %d traces processed", ngather, ntotal);

    return(CWP_Exit());
}


