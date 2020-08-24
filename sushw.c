/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* SUSHW: $Revision: 1.21 $ ; $Date: 2004/01/28 21:30:40 $		*/

#include <su.h>
#include <segy.h>
#include <header.h>

#ifndef VALUETYPE
#define VALUETYPE double
#define EALLOC1VALUETYPE(x) ealloc1double(x)
#define FREE1VALUETYPE(x) free1double(x)
#endif

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SUSHW - Set one or more Header Words using trace number, mod and              ",
"	 integer divide to compute the header word values or input              ",
"	 the header word values from a file                                     ",
"                                                                               ",
" Required Parameters for setting headers from binary infile:                   ",
" key=key1,key2 ... is the list of header fields as they appear in infile       ",
" infile= 	binary file of values for field specified by                    ",
" 		key1,key2,...                                                   ",
"                                                                               ",
" Required Parameters for setting headers from ascii infile/table:              ",
" match=key_1,key_2 ... is the list of matching header fields as they appear in infile",
" key=key1,key2     ... is the list of header fields as they appear in infile   ",
" infile=           ascii table of values for field specified by                ",
"                   key_1  key_2 ...  key1 key2 ...                             ",
"                   ...               ...                                       ",
"                                                                               ",
"                   values must be separated by white space                     ",
"                   the number of columns must be the same as                   ",
"                   that of keys to be matched and to be set                    ",
"                                                                               ",
" Required Parameters for setting headers from commandline value array:         ",
" match=key_1,key_2 ... is the list of matching header fields                   ",
" key=key1,key2     ... is the list of header fields                            ",
" values=           array of values for fields specified by                     ",
"                   key_1,key_2,...,key1 key2,...,...,...                       ",
"                                                                               ",
"                   the number of values must be the multiplum of               ",
"                   that of keys to be matched and to be set                    ",
"                                                                               ",
" Optional parameters ():                                                       ",
" key=cdp,...			header key word(s) to set                       ",
" a=0,...			value(s) on first trace                         ",
" b=0,...			increment(s) within group                       ",
" c=0,...			group increment(s)                              ",
" d=0,...			trace number shift(s)                           ",
" j=ULONG_MAX,ULONG_MAX,...	number of elements in group                     ",
" scale=1.0,...                 scaling factor(s) for key word(s) to set        ",
"                                                                               ",
" match=key_1,...	header key word(s) to match, no default                 ",
" maxlines=512          max. number of lines in ascii infile (if match=key(s) set)",
" infile=               ascii table of values for field specified by            ",
"                       key_1  key_2 ...  key1 key2 ...                         ",
        "                       or binary file containing values set to header words",
"                                                                               ",
" values=               array containing match keys and header words to be set  ",
" sync=0                =1 match key(s) in trace header are sorted as in ascii infile",
"                          this enable populating trace header keys in on run   ",
"                          thus greatly speed up for large dataset              ",
" interp=0		=1 linearly interpolate and constantly extrapolate      ",
" 			   value(s) according to matching keyword               ",
" sort=1,,		sort orders for match keys, used only if interp=1       ",
" 			e.g. =1,-1,1 asc for 1st & 3rd, desc for 2nd key        ",
"                                                                               ",
" verbose=0		=n print out info to stderr for every n-th trace matched",
"                                                                               ",
" There are 4 possibilities to set head words:                                  ",
"                                                                               ",
" ... compute header fields                                                     ",
"   sushw <stdin >stdout key=cdp,.. a=0,..  b=0,.. c=0,.. d=0,.. j=..,..        ",
"                                                                               ",
" ... or set headers from a binary file                                         ",
"   sushw <stdin > stdout  key=key1,..    infile=binary_file		",
" 									",
" ... or set headers from a ascii file				        ",
"   sushw match=key_1,..  scale=1.0,.. maxlines=1024 \\                 ",
"         key=key1,..    infile=ascii_file	<stdin > stdout         ",
" 									",
" ... or set headers from an array of values as an command line argument",
"   sushw match=key_1,.. key=key1,.. interp=1 <stdin >stdout \\         ",
"     values=k1,k2,..v1,v2,..,...                                       ",
" 									",
" Notes:								",
" Fields that are getparred must have the same number of entries as key	",
" words being set. Any field that is not getparred is set to the default",
" value(s) above. Explicitly setting j=0 will set j to ULONG_MAX.	",
" 									",
" The value of each header word key is computed using the formula:	",
" 	i = itr + d							",
" 	val(key) = ( a + b * (i % j) + c * (i / j) ) * scale		",
" where itr is the trace number (first trace has itr=0, NOT 1)		",
" 									",
" Examples:								",
" 1. set every dt field to 4ms						",
" 	sushw <indata key=dt a=4000 |...				",
" 2. set the sx field of the first 32 traces to 6400, the second 32 traces",
"    to 6300, decrementing by -100 for each 32 trace groups		",
"   ...| sushw key=sx a=6400 c=-100 j=32 |...				",
" 3. set the offset fields of each group of 32 traces to 200,400,...,6400",
"   ...| sushw key=offset a=200 b=200 j=32 |...				",
" 4. perform operations 1., 2., and 3. in one call			",
"  ..| sushw key=dt,sx,offset a=4000,6400,200 b=0,0,200 c=0,-100,0 j=0,32,32 |",
" 									",
" In this example, we set every dt field to 4ms.  Then we set the first	",
" 32 shotpoint fields to 6400, the second 32 shotpoint fields to 6300 and",
" so forth.  Next we set each group of 32 offset fields to 200, 400, ...,",
" 6400.									",
" 									",
" Example of a typical processing sequence using suchw:			",
"  sushw <indata key=dt a=4000 |					",
"  sushw key=sx a=6400 c=-100 j=32 |					",
"  sushw key=offset a=200 b=200 j=32 |			     		",
"  suchw key1=gx key2=offset key3=sx b=1 c=1 |		     		",
"  suchw key1=cdp key2=gx key3=sx b=1 c=1 d=2 >outdata	     		",
" 									",
" Again, it is possible to eliminate the multiple calls to both sushw and",
" sushw, as in Example 4.						",
" 									",
" Reading header values from a binary file:				",
" If the parameter infile=binary_file is set, then the values that are to",
" be set for the fields specified by key=key1,key2,... are read from that",
" file. The values are read sequentially from the file and assigned trace",
" by trace to the input SU data. The infile consists of C (unformated)	",
" binary floats in the form of an array of size (nkeys)*(ntraces) where	",
" nkeys is the number of floats in the first (fast) dimension and ntraces",
" is the number of traces.						",
" 									",
" Reading header values from a ascii table:				",
" If the parameter infile=ascii_file is set, then the entire table is read",
" in and cached in the memory. When a trace is read the values of trace ",
" headers specified by key=key_1,key_2,... are compared to the cached table.",
" If a match is found, the trace fields specified by key=key1,key2,... are",
" set to the values of the corresponding record. If no match is found, the",
" header values remain unchanged by default setting interp=0. However, if",
" interp=1, header values are linearly interpolated and constantly extrapolated",
" by matching headers. The interp=1 option is meant to load parameters",
" like muting table, scaling factors etc. Interpolation/extraplolation ",
" works only for up to three keys, e.g. inline, crossline, and offset.",
" The value is first interpolated between two adjacent offset, then crossline,",
" and last inline. The parameter to be interpolated should sampled densely",
" enough along every inline that are presented in the table.            ",
" 									",
" Comment: 								",
" Users wishing to edit one or more header fields (as in geometry setting)",
" may do this via the following sequence:				",
"     sugethw < sudata output=geom  key=key1,key2 ... > hdrfile 	",
" Now edit the ASCII file hdrfile with any editor, setting the fields	",
" appropriately. Convert hdrfile to a binary format via:		",
"     a2b < hdrfile n1=nfields > binary_file				",
" Then set the header fields via:					",
"     sushw < sudata infile=binary_file key1 key2 ... > sudata.edited	",
" Alternatively header fields can be output with key values, edited and ",
" read in again:                                                        ",
"     sugethw < sudata output=geom  key=key_1,...,key1,key2 ... > hdrfile",
"     sushw match=key_1,... key=key1,key2,... infile=hdrfile.edited     ",
"           < sudata > sudata.edited					",
" 									",
" 									",
" Caveat: 								",
" If the (number of traces)*(number of key words) exceeds the number of	",
" values in the infile then the user may still set a single header field",
" on the remaining traces via the parameters key=keyword a,b,c,d,j.	",
"  									",
" Example:								",
"    sushw < sudata key=key1,key2 ... infile=binary_file [Optional Parameters]",
"                                                                               ",
" Version 2.0.0   Last updated June, 2012 by Sanyu Ye                           ",
"                                                                               ",
NULL};

/* Credits:
 *	SEP: Einar Kajartansson
 *	CWP: Jack K. Cohen
 *      CWP: John Stockwell, added multiple fields and infile= options
 *      RWS: Sanyu Ye, added ascii file/array of multiple fields according to matching keys
 *
 * Caveat:
 *	All constants are cast to doubles.
 */
/**************** end self doc ****************************************/


segy tr;

int verbose = 0;

/* Prototypes */
double mod(double x, double y);
void setval(cwp_String type, Value *valp, double a, double b,
        double c, double i, double j, double scale);

int FindInterval(int ikey, VALUETYPE dv, VALUETYPE *dprev, VALUETYPE **keys, int *order,
        int nRows, int *nbegp,  int* n1p, int* n2p, VALUETYPE* d1p, VALUETYPE* d2p);

void interp1(int ikey,  int* sort, int n1, int n2, VALUETYPE d, VALUETYPE **keys, VALUETYPE **var, int nv, VALUETYPE* v);

void interp2(int sort, VALUETYPE d, VALUETYPE d1, VALUETYPE d2, VALUETYPE *v1, VALUETYPE *v2, int nv, VALUETYPE* v);

int interp(int nkeys, VALUETYPE *dv, VALUETYPE **keys, int *sort, int nRows, VALUETYPE **var, int nv, VALUETYPE* v);

int main(int argc, char **argv) {
    cwp_String mkey[SU_NKEYS];  /* array of keywords to match		*/
    cwp_String mtype[SU_NKEYS]; /* array of types name		*/
    int mindex[SU_NKEYS];	/* index array of keywords to match	*/
    Value valm[SU_NKEYS];	/* value array of matching keywords from trace */

    int mnkeys;		/* number of header fields set for matching keys	*/
    int mcount=0;		/* number of traces matched 	*/

    double *values = NULL;  /* pointer containing value array */
    double **mkeys = NULL;  /* pointer containing values of matching keys */
    double **vkeys = NULL;  /* pointer containing values of keys to be set */

    int maxRows = 0;        /* max. number of lines of ascii infile */
    int nRow = 0;           /* actual number of lines of ascii infile read in*/

    cwp_String key[SU_NKEYS];  /* array of keywords			*/
    cwp_String type[SU_NKEYS]; /* array of keywords			*/
    int index[SU_NKEYS];	/* name of type	of getparred key	*/
    double dintp[SU_NKEYS]; /* key values obtained through interpolation */

    int ikey;		/* key counter 				*/
    int nkeys;		/* number of header fields set		*/
    int count=0;		/* number of header fields from file	*/
    int nvalues=0;		/* number of elements in value array	*/

    double i;		/* parameters for computing fields	*/
    int itr = 0;		/* trace counter 			*/
    Value val;		/* value of key field to set		*/

    char *infile="";	/* name of input file of header values	*/
    FILE *infp=NULL;	/* pointer to input file		*/
    cwp_Bool from_file=cwp_false; /* is the data from infile?	*/

    float *afile=NULL;	/* array of "a" values from file	*/
    double *a=NULL;		/* array of "a" values			*/
    double *b=NULL;		/* array of "b" values			*/
    double *c=NULL;		/* array of "c" values			*/
    double *d=NULL;		/* array of "d" values			*/
    double *j=NULL;		/* array of "j" values			*/
    double *scale=NULL;	/* scaling factor for keys to be set	*/
    int *sort=NULL;	        /* sorting order for keys to be matched	*/
    int n;			/* number of a,b,c,d,j values		*/

    int interpolate = 0;    /* if set !=0, interpolate/extrapolate */
    int sync = 0;           /* if set !=0, assume data and input ascii have exactly same sort order in match keys */
    int ipos = 0;           /* global position in infile */

    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* Get "key" values */
    if ((nkeys=countparval("key"))!=0) {
        getparstringarray("key", key);

    } else {
        nkeys = 1;
        key[0]="cdp";
    }

    /* get types and indexes corresponding to the keys */
    for (ikey=0; ikey<nkeys; ++ikey) {
        type[ikey]=hdtype(key[ikey]);
        index[ikey]=getindex(key[ikey]);
    }

    /* Get matching "key" values */
    if ((mnkeys=countparval("match"))!=0) {
        getparstringarray("match", mkey);

        /* get types and indexes corresponding to the matching keys */
        for (ikey=0; ikey<mnkeys; ++ikey) {
            mtype[ikey]=hdtype(mkey[ikey]);
            mindex[ikey]=getindex(mkey[ikey]);
        }

        if (!getparint("maxlines", &maxRows))	maxRows = 512;
        if (!getparint("verbose", &verbose))	verbose = 0;
        if (!getparint("interp", &interpolate))     interpolate = 0;
    }

    /* get name of infile */
    getparstring("infile", &infile);

    if(!getparint("sync", &sync))  sync = 0;

    /* if infile is specified get specified keys from file */
    if (*infile!='\0') {

        /* open infile */
        if((infp=efopen(infile, "r"))==NULL)
            err("cannot open infile=%s\n", infile);

        /* set from_file flag */
        from_file=cwp_true;
    }

    // get value array
    if (mnkeys > 0 && !from_file) {
        nvalues = countparval("values");
        if ( nvalues % (mnkeys + nkeys) != 0) {
            err("Number of values must be multiplum of that of keys to be matched and set!");
        }
        values=ealloc1double(nvalues);
        getpardouble("values", values);
        from_file=cwp_true;
    }

    /* get sorting order of match keys if interpolation applies */
    if ( interpolate ) {
        if ((n=countparval("sort"))!=0) {
            if (n>mnkeys) {
                err("number of sort's cannot be large than key's!");
            } else {
                sort=ealloc1int(mnkeys);
                getparint("sort", sort);
                if (n<mnkeys) {
                    for (ikey=n; ikey<mnkeys; ++ikey) {
                        sort[ikey]=sort[n-1];
                    }
                }
            }
        } else { /* set default */
            sort=ealloc1int(mnkeys);
            for (ikey=0; ikey<mnkeys; ++ikey)
                sort[ikey]=1;
        }
    }

    /* get scaling factors */
    if ((n=countparval("scale"))!=0) {
        if (n>nkeys) {
            err("number of scale's cannot be large than key's!");
        } else {
            scale=ealloc1double(nkeys);
            getpardouble("scale", scale);
            if (n<nkeys) {
                for (ikey=n; ikey<nkeys; ++ikey) {
                    scale[ikey]=scale[n-1];
                }
            }
        }
    } else { /* set default */
        scale=ealloc1double(nkeys);
        for (ikey=0; ikey<nkeys; ++ikey)
            scale[ikey]=1.;
    }

    /* If not from file or value array, getpar a,b,c,d,j */
    if (!from_file && nvalues == 0) {
        /* get "a" values */
        a=ealloc1double(nkeys);
        if ((n=countparval("a"))!=0) {
            getpardouble("a",a);
            if (n<nkeys) {
                for (ikey=n; ikey<nkeys; ++ikey) {
                    a[ikey]=a[n-1];
                }
            }
        } else {
            for (ikey=0; ikey<nkeys; ++ikey) a[ikey]=0.;
        }

        /* get "b" values */
        b=ealloc1double(nkeys);
        if ((n=countparval("b"))!=0) {
            getpardouble("b",b);
            if (n<nkeys) {
                for (ikey=n; ikey<nkeys; ++ikey) {
                    b[ikey]=b[n-1];
                }
            }
        } else {
            for (ikey=0; ikey<nkeys; ++ikey) b[ikey]=0.;
        }

        /* get "c" values */
        c=ealloc1double(nkeys);
        if ((n=countparval("c"))!=0) {
            getpardouble("c",c);
            if (n<nkeys) {
                for (ikey=n; ikey<nkeys; ++ikey) {
                    c[ikey]=c[n-1];
                }
            }
        } else {
            for (ikey=0; ikey<nkeys; ++ikey) c[ikey]=0.;
        }

        /* get "d" values */
        d=ealloc1double(nkeys);
        if ((n=countparval("d"))!=0) {
            getpardouble("d",d);
            if (n<nkeys) {
                for (ikey=n; ikey<nkeys; ++ikey) {
                    d[ikey]=d[n-1];
                }
            }
        } else {
            for (ikey=0; ikey<nkeys; ++ikey) d[ikey]=0.;
        }

        /* get "j" values */
        if ((n=countparval("j"))!=0) {
            if (n!=nkeys)
                err("number of j values not equal to number of keys");

            j=ealloc1double(n);
            getpardouble("j", j);

            /* make sure that j!=0 */
            for (ikey=0; ikey<nkeys; ++ikey)
                if(j[ikey]==0) j[ikey]=ULONG_MAX;
        } else {
            j=ealloc1double(nkeys);
            for (ikey=0; ikey<nkeys; ++ikey) j[ikey]=ULONG_MAX;
        }
    } else if (mnkeys == 0) { /* if reading from a file */
        /* allocate space for afile */
        afile=ealloc1float(nkeys);
    } else {  /* if reading from a acsii file/table or value array*/
        if ( nvalues == 0 ) { // read from ascii file
            int bEnd = 0;

            mkeys = ealloc2double(mnkeys, maxRows);
            vkeys = ealloc2double(nkeys, maxRows);

            if (verbose) warn("m=%d n=%d file=%s" , mnkeys, nkeys, infile);
            /* reading all data from ascii infile */
            /*  fscanf returns: 0   : characters there, but no conversion (error)
             *		  EOF : eof before conversion
             *		  else: number of conversions
             */
            for (nRow=0; nRow < maxRows; ++nRow) {
                for (ikey = 0; ikey < mnkeys; ++ikey ) {
                    int ret = fscanf(infp, "%lf", &mkeys[nRow][ikey]);
                    if(ret == EOF) {
                        bEnd = 1; break;  /* else everything is okay: get out of the loop */
                    }
                    if(ret == 0) {
                        bEnd = 1; break;  /* else everything is okay: get out of the loop */
                    }
                }
                for (ikey = 0; ikey < nkeys && !bEnd; ++ikey ) {
                    int ret = fscanf(infp, "%lf", &vkeys[nRow][ikey] );
                    if(ret == EOF) {
                        bEnd = 1; break;  /* else everything is okay: get out of the loop */
                    }
                    if(ret == 0) {
                        bEnd = 1; break;  /* else everything is okay: get out of the loop */
                    }
                }

                if (bEnd) break;
            }
        } else { // read from value array
            maxRows = nvalues / (mnkeys + nkeys);
            mkeys = ealloc2double(mnkeys, maxRows);
            vkeys = ealloc2double(nkeys, maxRows);
            for (nRow = 0; nRow < maxRows; ++nRow) {
                for (ikey = 0; ikey < mnkeys; ++ikey ) {
                    mkeys[nRow][ikey] = (double) values[nRow*(mnkeys + nkeys) + ikey];
                }
                for (ikey = 0; ikey < nkeys; ++ikey ) {
                    vkeys[nRow][ikey] = (double) values[nRow*(mnkeys + nkeys) + mnkeys + ikey];
                }
            }
            nRow = maxRows;
        }

        if (verbose > 1) warn("%d lines read", nRow);
        if (verbose >= 10) {
            for (n=0; n < nRow; ++n) {
                fprintf(stderr, "%4d -> ", n);
                for (ikey = 0; ikey < mnkeys; ++ikey ) {
                    fprintf(stderr, " %6d ", (int)mkeys[n][ikey]);
                }
                for (ikey = 0; ikey < nkeys; ++ikey ) {
                    fprintf(stderr, " %.2f ", vkeys[n][ikey]);
                }
                fprintf(stderr, " END\n");
            }
        }
    }

    /* loop over traces */
    while (gettr(&tr)) {

        if (from_file && mnkeys > 0) {  /* match keys and set value from file to trace by trace */

            cwp_Bool isMatch = cwp_false;

            /* loop over matching key fields and get values */
            for (ikey=0; ikey<mnkeys; ++ikey) {
                /* get header values */
                gethval(&tr, mindex[ikey], &valm[ikey]);
            }

            if (!interpolate) { /* no interpolation, exact match of key values required */
                for (n = ipos; n < nRow && !isMatch; ++n) {
                    for (ikey=0, isMatch = cwp_true; isMatch && ikey<mnkeys; ++ikey) {
                        setval(mtype[ikey], &val, mkeys[n][ikey], 0, 0, 0, ULONG_MAX, 1.0);
                        isMatch = isMatch && !valcmp(mtype[ikey], valm[ikey], val);
                    }
                }

                if ( !sync ) ipos = 0;
                else if ( isMatch ) ipos = n - 1;
                else ipos = 0;

                if (isMatch) {
                    /*warn("found match for trace %d key[%d][0]=%d" , itr,  n-1,(int) mkeys[n-1][0]);*/
                    for (ikey=0; ikey<nkeys; ++ikey) {
                        setval(type[ikey], &val, vkeys[n-1][ikey], 0, 0, 0, ULONG_MAX, scale[ikey]);
                        puthval(&tr, index[ikey], &val);
                    }
                    ++mcount;
                }
            } else { /* interpolation/extrapolation */
                double d[3];

                if ( mnkeys > 3 ) err("Interpolation over >3 keys is not implemented");
                for (ikey=0; ikey<mnkeys; ++ikey){
                    d[ikey] = vtod(mtype[ikey], valm[ikey]);
                }

                interp(mnkeys, d, mkeys, sort, nRow, vkeys, nkeys, dintp);
                for (ikey=0; ikey<mnkeys && verbose>30; ++ikey){
                    fprintf(stderr, "%s=%d ", mkey[ikey], (int) d[ikey]);
                }
                for (ikey=0; ikey<nkeys && verbose>30; ++ikey){
                    fprintf(stderr, "dintp[%d]=%f ", ikey, dintp[ikey]);
                }
                if (verbose>30) fprintf(stderr, "\n");
                for (ikey=0; ikey<nkeys; ++ikey) {
                    setval(type[ikey], &val, dintp[ikey], 0, 0, 0, ULONG_MAX, scale[ikey]);
                    puthval(&tr, index[ikey], &val);
                }
                ++mcount;
            }
        } else if (from_file) {/* use the "a" value from file to trace by trace */
            if (efread(afile, FSIZE, nkeys, infp)!=0) {
                for (ikey=0; ikey<nkeys; ++ikey) {
                    double a_in;
                    a_in=(double) afile[ikey];
                    setval(type[ikey], &val, a_in,
                            0, 0, 0, ULONG_MAX, scale[ikey]);
                    puthval(&tr, index[ikey], &val);
                    ++count;
                }
            }
        } else { /* use getparred values of a,b,c,d,j */
            for (ikey=0; ikey<nkeys; ++ikey) {
                i = (double) itr + d[ikey];

                setval(type[ikey], &val, a[ikey], b[ikey],
                        c[ikey], i, j[ikey], scale[ikey]);
                puthval(&tr, index[ikey], &val);
            }
        }

        ++itr;
        puttr(&tr);
    }

    if (from_file && verbose) {
        if (infp) efclose(infp);
        if (mnkeys > 0) {
            warn("  %d traces of total %d set with new header values", mcount, itr);
        } else if (count < (int)(itr*nkeys) ) {
            warn("itr=%d > count=%d", (int) itr*count, count);
            warn("n traces=%d > data count =%d", (itr*nkeys), count);
        }
    }

    /* free allocated memory */
    if (a) free1double(a);
    if (b) free1double(b);
    if (c) free1double(c);
    if (d) free1double(d);
    if (j) free1double(j);
    if (sort) free1int(sort);
    if (scale) free1double(scale);
    if (afile) free1float(afile);
    if (mkeys) free2double(mkeys);
    if (vkeys) free2double(vkeys);
    if (values) free1double(values);

    return(CWP_Exit());
}


void setval( cwp_String type, Value *valp, double a, double b,
        double c, double i, double j, double scale) {
    switch (*type) {
        case 's':
            err("can't set char header word");
            break;
        case 'h':
            valp->h = scale * (a + b * mod(i, j) + c * ((int) (i/j)));
            break;
        case 'u':
            valp->u = scale * (a + b * mod(i, j) + c * ((int) (i/j)));
            break;
        case 'l':
            valp->l = (long) (scale * ((a + b * mod(i, j) + c * ((int) (i/j)))));
            break;
        case 'v':
            valp->v = (unsigned long) (scale * ((a + b * mod(i, j) + c * ((int) (i/j)))));
            break;
        case 'i':
            valp->i = scale * (a + b * mod(i, j) + c * ((int) (i/j)));
            break;
        case 'p':
            valp->p = scale * (a + b * mod(i, j) + c * ((int) (i/j)));
            break;
        case 'f':
            valp->f = scale * (a + b * mod(i, j) + c * ((int) (i/j)));
            break;
        case 'd':
            valp->d = scale * (a + b * mod(i, j) + c * ((int) (i/j)));
        default:
            err("unknown type %s", type);
            break;
    }
    return;
}


double mod(double x, double y)	/* As defined in Knuth, vol. 1	*/ {
    return y == 0.0 ? x : x - y * floor(x/y);
}

#include "interp_incl.c"

