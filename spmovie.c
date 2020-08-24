/* Copyright (c) READ Well Services, 2009.*/
/* All rights reserved.                      */

/* SPMOVIE: $Revision: 1.0 $ ; $Date: 2009/02/20 $        */

#include <su.h>
#include <segy.h>
#include <header.h>
#include <unistd.h>
#include <sys/types.h>
#include <string.h>
#include <dirent.h>
#include <unistd.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                               ",
" SPMOVIE - display data as movie frame by frame                ",
"                                                               ",
"   spmovie <infile >outfile  [optional parameters]             ",
"                                                               ",
" Optional Parameters:                                          ",
"   key=cdp         frame key                                   ",
"                                                               ",
"   names=          keys thier values append to file name       ",
"                   for better identification of individual su file",
"                                                               ",
"   desc=           descriptive name of movie directory         ",
"                   should not contain blanks and special characters",
"                                                               ",
"   view=0          =1; start spviewer automatically            ",
"                                                               ",
"   rhost=          start movie (spviewer) at remote host       ",
"                                                               ",
"   verbose=0       >0 echo information                         ",
"                                                               ",
" spmovie is a thin wrapper around spviewer. It breaks the incoming segy trace ",
" stream upon the change of the key value into separate files stored under a ",
" directory, while starting the spviewer in a background process ",
"                                                             ",
" Notes:                                                      ",
"  Incoming su stream must be sorted in right order           ",
"  A temparary directory named movieNNNNN-desc (NNNNN process id)  ",
"  created under the working directory. Don't forget to delete",
"  it after finishing the viewing"
"                                                             ",
" Exsamples:                                                  ",
"   spdbread dbpath=mydata.db select=\"cdpt+(1:240:10)|ep+(1000:2000:10)|sx+\" |\\ ",
"   sugain tpow=2 |\\                                         ",
"   spmovie key=ep names=cdpt desc=every10rec-every10shot       ",
"                                                             ",
"   susort cdp offset < mydata.su |\\ ",
"   sugain tpow=2 |\\                                         ",
"   spmovie rhost=my_pc verbose=10 ",
"                                                             ",
"   segyread tape=mydata.segy |\\ ",
"   susort2 cdp offset |\\                                    ",
"   sugain tpow=2 |\\                                         ",
"   spmovie rhost=my_pc verbose=1 ",
NULL};

/*
* Credits: READ Well Service, Sanyu Ye, Feb. 2009.
*
*/
/**************** end self doc ********************************/

/** main **/
int main(int argc, char **argv)
{
    cwp_String key;    /* header key word from segy.h */
    cwp_String type;   /* type of key                 */
    int index;         /* index of key                */
    Value val, valnew; /* value of key                */
    cwp_String mkey[8];  /* array of keywords to append		*/
    cwp_String mtype[8]; /* array of types name		*/
    int mindex[8];	/* index array of keywords to append	*/
    Value valm;         /* value array of name keywords from trace */
    int mnkeys, ikey;		/* number of header fields set for appending keys	*/
    int nsegy;         /* Bytes read in for a trace   */
    int view;          /* flag for starting spviewer  */
    int verbose;       /* flag for echoing information */
    int ntr, ntotal=0; /* trace counters */
    int ngathers=0;    /* gather counters */
    char *rhost;       /* remote host name	*/
    char tmpstr[BUFSIZ];  /* string buffer to hold key value */
    char tmpdir[BUFSIZ];  /* string buffer to hold movie directory */
    char tmpfile[BUFSIZ]; /* string buffer to hold filename */
    char cmdviewer[BUFSIZ]; /* string buffer to hold command starting spviewer */
    FILE *fpdata;           /* fp for su data file		*/
    cwp_String desc;    /* descriptive name append to mivie directory name */

    segy tr;

    /* hook up getpar to handle the parameters */
    initargs(argc,argv);
    requestdoc(0);

    if (!getparint("view", &view))          view=0;
    if (!getparint("verbose", &verbose))    verbose=0;
    if (!getparstring("desc",&desc))        desc="";
    if (!getparstring("rhost",&rhost))      rhost="";

    /* read first trace */
    nsegy = gettr(&tr);

    /* get SU sorting key */
    if (!getparstring("key", &key)) key = "cdp";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    /* Get name "key" values */
    if ((mnkeys=countparval("names"))!=0) {
        getparstringarray("names", mkey);

        /* get types and indexes corresponding to the matching keys */
        for (ikey=0; ikey<mnkeys; ++ikey) {
            mtype[ikey]=hdtype(mkey[ikey]);
            mindex[ikey]=getindex(mkey[ikey]);
        }

        if(verbose > 1) warn("  %d keys for name input", mnkeys);
    }

    // create movie directory
    sprintf(tmpdir, "%s/movie%d", getcwd(NULL,0), (int)getpid());
    strcat(tmpdir, "-");
    strcat(tmpdir, desc);

    if (verbose > 1) warn(" Data stored under %s", tmpdir);
    if (mkdir(tmpdir, 0777)) {
        err("\n\tFailed to create directory %s\n", tmpdir);
    }
    sync(); // ensure that the dir is created on disk

    // start spviewer
    if (view) {
        if (STREQ(rhost, "")) {
            //sprintf(cmdviewer, " source=movie%d ", getpid());
            //execlp("spviewer", "spviewer", cmdviewer, "&", NULL);
            sprintf(cmdviewer, "spviewer source=%s &", tmpdir);
            system(cmdviewer);
        }
        else {
            //sprintf(cmdviewer, " \"/soft/spdb/bin/spviewer source=%s \" ", tmpdir);
            //execlp("ssh", "ssh", rhost, cmdviewer, NULL);
            sprintf(cmdviewer, "ssh %s \"/soft/spdb/bin/spviewer source=%s \" &", rhost, tmpdir);
            system(cmdviewer);
        }
        if(verbose > 1) {
            warn("cmdline=%s\n", cmdviewer);
        }
    }

    // create first data file 
    sprintf(tmpfile, "%s/d%05d", tmpdir, ++ngathers);
    /* loop over name key fields and get values and append to name*/
    for (ikey=0; ikey<mnkeys; ++ikey) {
        gethval(&tr, mindex[ikey], &valm);
        sprintf(tmpstr, "-%d",  vtoi(mtype[ikey], valm));
        strcat(tmpfile, tmpstr);
    }
    sprintf(tmpstr, "-%d.su", vtoi(type, val));
    strcat(tmpfile, tmpstr);
    if (verbose > 1) warn(" First file name %s", tmpfile);
    fpdata = efopen(tmpfile, "w+");
    sprintf(cmdviewer, "spviewer source=%s &", tmpdir);
    
    ntr = 0;
    do {
        gethval(&tr, index, &valnew);
        //if(verbose >= 100) {
        //    warn("key old=%d  new=%d\n", vtoi(type, val), vtoi(type, valnew));
        //}
        if (!valcmp(type, val, valnew) ) { /* same key */
            ++ntr;
        } else {  /* new gather */
            ntotal += ntr;
            // close previous file and open a new file
            efclose(fpdata);
            sync(); // ensure data is written to disk 
            if (verbose > 2) warn(" %d-th gather has %d traces.", ngathers, ntr);

            // create new data file
            sprintf(tmpfile, "%s/d%05d", tmpdir, ++ngathers);
            /* loop over name key fields and get values and append to name*/
            for (ikey=0; ikey<mnkeys; ++ikey) {
                gethval(&tr, mindex[ikey], &valm);
                sprintf(tmpstr, "-%d",  vtoi(mtype[ikey], valm));
                strcat(tmpfile, tmpstr);
            }
            sprintf(tmpstr, "-%d.su", vtoi(type, valnew));
            strcat(tmpfile, tmpstr);
            fpdata = efopen(tmpfile, "w+");
            
            ntr = 1; // reset
        }

        fputtr(fpdata, &tr); // output trace

        val = valnew; // cache current key value for later comparison

        /*read next trace */
        nsegy = gettr(&tr);
    } while (nsegy > HDRBYTES);

    ntotal += ntr;
    efclose(fpdata);  // close last file
    if (verbose >= 10) warn(" %d-th gather has %d traces.", ngathers, ntr);

    if (verbose) warn(" Totally %d gathers defined with total %d traces.", ngathers, ntotal);

    return(CWP_Exit());
}


