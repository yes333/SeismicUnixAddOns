/* Copyright (c) READ Well Services 2008.*/
/* All rights reserved.                      */

/* SUMERGE: $Revision: 1.0 $ ; $Date: 2008/02/05 17:38:30 $    */

#include "su.h"
#include "segy.h"
#include "header.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                       ",
" SUMERGE - merge up to 4 su data files, trace by trace from each file  ",
"                                                                       ",
" sumerge files=data1,data2,... [ scale=1.0,1.0,... trid= ] >stdout     ",
"                                                                       ",
" Required parameters:                                                  ",
"   files=data1,data2,..   up to 4 su data files                        ",
"                                                                       ",
" Optional parameter:                                                   ",
"   scale=1.0,1.0,...      scaler for data files before output          ",
"                                                                       ",
"   trid=                  array of trace id from individual files      ",
"                          =1,2,3,4 e.g. set trid                       ",
"                                                                       ",
"   key=cdpt               keyword used to match input scaling factors  ",
"   nrcv=240               max. number of receivers                     ",
"   scfile=                file containing scaling/macthing factors     ",
"                                                                	",
"   verbose=1               verbose = 0 not echo info                   ",
"                                                                       ",
NULL};

/* Credits:
*    Sanyu Ye, READ Well Service, Feb. 2008
*
*/
/**************** end self doc ***********************************/

segy tr;

int main(int argc, char **argv)
{
    cwp_String key;		/* header key word from segy.h		*/
    cwp_String type;            /* type of key				*/
    Value val;                  /* value of key			*/
    double dval, dvallast;      /* double value of key			*/
    int index;                  /* index of key				*/
    FILE *fp[4]={NULL,NULL,NULL,NULL};   /* file pointer for up to 4 files    */
    cwp_String  filename[4];             /* filenames of input data */
    int nfiles=0;                        /* number of input files  */
    int nt=0;                            /* number of sample points on trace    */
    int nsegy=0;                         /* number of bytes on traces        */
    int ntr=0;                           /* number of trace being processed    */
    int i, j, n;                         /* loop counter  */
    int verbose;                         /* flag for echoing information */
    int trid[4];                         /* trace ids to be set  */
    float scale[4]={0.0,0.0,0.0,0.0};    /* weighting factors for files*/
    float scaler=1.0;                    /* readin global weighting factors for geophone comp */
    
    char *scfile="";	/* name of input file of reflection coefficients	*/
    FILE *scfp=NULL;	/* pointer to input file		*/
    cwp_Bool from_file=cwp_false;   /* is the data from infile?	*/
    int nrcv, nRow;                 /* max. and actual scaling factor values read in */
    double *mkeys = NULL;           /* pointer containing values of matching keys */
    double *vkeys = NULL;           /* pointer containing values of keys to be set */

    /* Initialize */
    initargs(argc, argv);
    requestdoc(1); /* 1 file args required */

    if (!getparint("verbose", &verbose)) verbose= 1;

    if ((nfiles=countparval("files"))!=0) {
        getparstringarray("files",filename);
    } else if (nfiles > 4) {
        err("Number of input data files cannot larger than 4");
    } else {
        err("Input data files must be given!");
    }

    /* get scaling factors */
    if ((n=countparval("scale"))!=0) {
        if (n>nfiles) {
            err("number of scale's cannot be large than file's!");
        } else {               
            getparfloat("scale",scale);         
            if (n<nfiles) {
                for (i=n; i<nfiles; ++i) {
                    scale[i]=scale[n-1];
                }
            }
        }
    } else { /* set default */
            for (i=0; i<nfiles; ++i)
                    scale[i]=1.;
    }
    if (verbose>1) {
        for (i=0; i<nfiles; ++i) fprintf(stderr, "scale[%d]=%.3f  ", i, scale[i]);
        fprintf(stderr, "\n");
    }

    /* matching key  */
    if (!getparstring("key", &key))	 key="cdpt";
    type = hdtype(key);
    index = getindex(key);

    /* get name of infile for reflection coefficients*/
    getparstring("scfile",&scfile);
    if (*scfile!='\0') {
        /* open infile */
        if((scfp=efopen(scfile,"r"))==NULL)
                err("cannot open scfile=%s\n",scfile);

        /* set from_file flag */
        from_file=cwp_true;
    }
    if (!getparint("nrcv", &nrcv)) nrcv=240;
    if (from_file) { 
        mkeys = ealloc1double(nrcv);
        vkeys = ealloc1double(nrcv);

        if (verbose) warn("  Match key=%s for reading scaling factor from file=%s", key, scfile);
        /* reading all data from ascii infile */
        /*  fscanf returns: 0   : characters there, but no conversion (error)
        *		  EOF : eof before conversion
        *		  else: number of conversions 
        */
        for (nRow=0; nRow < nrcv; ++nRow) {
            int ret = fscanf(scfp, "%lf", &mkeys[nRow]);
            if(ret == EOF || ret == 0) break;  /* else everything is okay: get out of the loop */
            ret = fscanf(scfp, "%lf", &vkeys[nRow] );
            if(ret == EOF || ret == 0) break;  /* else everything is okay: get out of the loop */
        }
        
        if (verbose) warn(" %d values are read in from scfile=%s", nRow, scfile);
    }

    /* get trace ids */
    if ((n=countparval("trid"))!=0) {
        if (n != nfiles) {
            err("number of trace IDs must be equal to file's!");
        } else {               
            getparint("trid",trid);
        }
    } else { /* set default */
            for (i=0; i<nfiles; ++i)
                    trid[i]=0;
    }
    
    /* Open files given for reading */
    for (i = 0; i < nfiles; ++i) {
        fp[i] = efopen(filename[i], "r");
    }

    dvallast = 0.0;
    do { /* loop over the traces */        
        for (i = 0; i < nfiles; ++i) { /* loop over the files*/
            nsegy =  fgettr(fp[i], &tr);
            if ( !nt && nsegy ) {
                nt = tr.ns;
                if (verbose) warn("trace info: ns=%d ", nt);
            }
            
            if (from_file) { // retrieve scaler
                gethval(&tr, index, &val);
                dval = vtod(type, val);
                if (dval != dvallast) {
                    for (j=0; j<nRow; ++j) { //loop over to find right scaler 
                        if ( dval == mkeys[j] ) { 
                            scaler = vkeys[j];
                            break;
                        }
                    }
                    dvallast = dval;
                }
            }

            if ( nsegy <= HDRBYTES) {
                if ( i > 0 ) err(" Inconsistent data:  %d-th trace for %-th file %s missing", ntr+1, i+1, filename[i]);
                else break; // EOF for first file, out of loop
            } 

            if ( nt != tr.ns ) err("Inconsistent data:  %d-th trace for %-th file %s --> nt=%d != %d", 
                    ntr, i+1, filename[i], tr.ns, nt);
            if ( trid[i] ) tr.trid = trid[i]; 
            if ( scale[i] != 1.0 ) {  /* apply scaling factors */
                for (j = 0; j < nt; ++j) tr.data[j] *= scale[i];               
            }
            if (from_file && i > 0) { /* apply scaling factor to geophone components */
                for (j = 0; j < nt; ++j) tr.data[j] *= scaler;               
            }
            puttr(&tr);
        }

        if (nsegy > HDRBYTES) ++ntr;

    } while (nsegy > HDRBYTES);

    for (i = 0; i < nfiles; ++i) efclose(fp[i]);
    
    if (verbose) warn("Totally %d traces for each of %d files merged", ntr, nfiles);

    free1double(mkeys);
    free1double(vkeys);
    
    return(CWP_Exit());
}
