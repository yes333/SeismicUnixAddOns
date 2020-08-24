/* Copyright (c) Colorado School of Mines, 2008.*/
/* All rights reserved.                       */

/* SUMUTE: $Revision: 1.31 $ ; $Date: 2006/11/07 22:58:42 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
char *sdoc[] = {
" 	   								",
" SUMUTE - mute above/below a polygonal curve given by user input or a	",
"          header word. Surgical mute can be performed between two time ",
"          values defined by two key header words                       ",
" 	   								",
" sumute <stdin >stdout [optional parameters]                           ",
" 									",
" Muting mode:                                                          ",
"   mode=0              =0 top mute (above)                            	",
"                       =1 bottom mute (below)                         	",
"                       =2 to mute below AND above a straight line	",
"                           in that case xmute,tmute describes the total",
"                           time length of the muted zone		",
"                       =3 surgical mute between keywords muts and mute	",
" 									",
" Optional parameters for muting curve, default via header key word     ",
"   tkey=mute           header word containing tmute value in ms	",
"   bias=0.0            constant shift in ms added to header value      ",
" 									",
"  ... or input via arrays:						",
"   xmute=              array of position values as specified by	",
"                       the `key' parameter				",
"   tmute=              array of corresponding time values (sec)	",
"                       in case of air wave muting, correspond to 	",
"                       air blast duration				",
" 	   								",
"  ... or input via files:						",
"   nmute=              number of x,t values defining mute		",
"   xfile=              file containing position values as specified by	",
"                       the `key' parameter				",
"   tfile=              file containing corresponding time values (sec)	",
" 									",
" More optional parameters:						",
"   key=offset          Key header word specifying trace offset 	",
"                       e.g. =tracl  use trace number instead           ",
"   scale=tr.scalco     scaling factor, default tr.scalco for key=offset",
"        =1.0           default 1.0 for other than default key=offset   ",
" 									",
"   ntaper=0            number of points to taper before hard		",
"                       mute (sine squared taper)			",
"   type=4              taper type					",
"                       =1 linear (default)                             ",
"                       =2 cosine					",
"                       =3 sine                                         ",
"                       =4 sine^2					",
"                       =5 gaussian (+/-2.0)                            ",
"                       =6 gaussian (+/-3.8)                            ",
" 									",
"   linvel=330          linear velocity					",
"   tm0=0               time shift of linear mute at \'key\'=0		",
" 									",
"   verbose=0           >0  echoes information				",
    " 									",
" Notes: 								",
" 									",
" The tmute interpolant is extrapolated to the left by the smallest time",
" sample on the trace and to the right by the last value given in the	",
" tmute array.								",
"									",
" The files tfile and xfile are files of binary (C-style) floats.	",
"									",
" In the context of this program \"above\" means earlier time and	",
" \"below\" means later time (above and below as seen on a seismic section.",
"									",
" The mode=2 option is intended for removing air waves. The mute is	",
" is over a narrow window above and below the polygonal line specified	",
" by the values in tmute, xmute or tfile and xfile.			",
"									",
" If data are spatial, such as the (z-x) output of a migration, then    ",
" depth values are used in place of times in tmute and tfile. The value ",
" of the depth sampling interval is given by the d1 header field	",
"									",
" If data are complex frequency, such as output of spfk and spfk2taup,  ",
" then frequency values in Hz are used in place of times  The value of  ",
" the frequency sampling interval is given by the d1 header field	",
"									",
" Mute time retrieved from header must be in milisecond, meter or 1/10 Hz",
"									",
" Caveat: if data are seismic time sections, then tr.dt must be set. If ",
" data are seismic depth/complex freqeuncy sections, then tr.trid must  ",
" be set to 30 and 64 respectively, and tr.d1 header field must be set.	",
"									",
NULL};
    
/* Credits:
 *
 *	SEP: Shuki Ronen
 *	CWP: Jack K. Cohen, Dave Hale, John Stockwell
 *	DELPHI: Alexander Koek
 *	RWS: Sanyu Ye
 *
 * Trace header fields accessed: ns, dt, delrt, key, tkey, trid, d1
 * Trace header fields modified: muts or mute
 */
/**************** end self doc ***********************************/
    
/* forward declaration of prototype functions */
float* MakeTaper(int n, int type, int verbose);

segy tr;

int main(int argc, char **argv)
{
    char *key;      /* header key word from segy.h		*/
    char *type;     /* ... its type				*/
    int index;      /* ... its index			*/
    Value val;      /* ... its value			*/
    float fval;     /* ... its value cast to float		*/
    char *tkey;     /* header key word for tmute		*/
    char *ttype;    /* ... its type				*/
    int tindex;     /* ... its index			*/
    Value tval;     /* ... its value			*/
    float bias;     /* ... constant shift of its value 		*/
    float scale;    /* scaling factor for key	*/
    float linvel;   /* linear velocity			*/
    float tm0;      /* time shift of linear mute for 'key'=0*/
    float *xmute=NULL;   /* array of key mute curve values	*/
    float *tmute=NULL;   /* ...		mute curve time values 	*/
    float *taper=NULL;   /* ...		taper values	*/
    int nxmute;         /* number of key mute values		*/
    int ntmute;         /* ...		mute time values 	*/
    int ntaper, tptype;	/* taper window length and type		*/
    int below;          /* mute below curve			*/
    int mode;           /* kind of mute (top, bottom, linear)	*/
    int absolute;       /* Take absolute value of key for mode=2 */
    int verbose;        /* flag for printing information	*/

    int ntr=0;          /* trace counter 		*/
    int nxtmute;	/* number of mute values 		*/
    cwp_String xfile="";	/* file containing positions by key	*/
    FILE *xfilep;		/* ... its file pointer			*/
    cwp_String tfile="";	/* file containing times	 	*/
    FILE *tfilep;		/* ... its file pointer			*/

    cwp_Bool seismic;	/* cwp_true if seismic, cwp_false not seismic */
    cwp_Bool cmplxtr;	/* cwp_true if complex trace data like FK */
    cwp_Bool hdr = cwp_false;	/* cwp_true if header value is used for muting */

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);


    /* Get parameters */
    if ( (nxmute = countparval("xmute")) && (ntmute = countparval("tmute")) ) {
        if (nxmute != ntmute)
            err("lengths of xmute, tmute must be the same");
        xmute = ealloc1float(nxmute);	getparfloat("xmute", xmute);
        tmute = ealloc1float(nxmute);	getparfloat("tmute", tmute);
    } else if (getparint("nmute", &nxtmute)){
        nxmute = nxtmute;
        xmute = ealloc1float(nxtmute);
        tmute = ealloc1float(nxtmute);
        if (!getparstring("tfile", &tfile)) err("tfile must be given");
        if (!getparstring("xfile", &xfile)) err("xfile must be given");
        if((xfilep=fopen(xfile, "r"))==NULL)
            err("cannot open xfile=%s\n", xfile);
        if (fread(xmute, sizeof(float), nxtmute, xfilep)!=nxtmute)
            err("error reading xfile=%s\n", xfile);
        fclose(xfilep);

        if((tfilep=fopen(tfile, "r"))==NULL)
            err("cannot open tfile=%s\n", tfile);
        if (fread(tmute, sizeof(float), nxtmute, tfilep)!=nxtmute)
            err("error reading tfile=%s\n", tfile);
        fclose(tfilep);
    } else {
        if (!getparstring("tkey", &tkey)) tkey = "mute";
        ttype = hdtype(tkey);
        tindex = getindex(tkey);
        if (!getparfloat("bias", &bias))  bias = 0;
        hdr = cwp_true;
    }

    if (!getparint("ntaper", &ntaper))	ntaper = 0;
    if(!getparint("type", &tptype)) tptype = 4;
    if(tptype < 1 || tptype > 6) err( " Unkonwn taper type (=%d)", tptype);
    if (!getparint("mode", &mode))		mode = 0;
    if (getparint("below", &below))	{
        mode = below;
        warn("use of below parameter is obsolete. mode value set to %d \n", mode);
    }
    if (!getparint("verbose", &verbose))	verbose = 0;
    if (!getparint("absolute", &absolute))	absolute = 1;
    if (!getparfloat("linvel", &linvel))	linvel = 330;
    if (!getparfloat("tm0", &tm0))		tm0 = 0;
    if (linvel==0) err("linear velocity can't be 0");

    /* get key type and index */
    if (!getparstring("key", &key))     key = "offset";
    type = hdtype(key);
    index = getindex(key);

    /* Set up taper weights if tapering requested */
    if (ntaper > 0) taper = MakeTaper(ntaper, tptype, verbose);

    /* Get info from first trace */
    if (!gettr(&tr)) err("can't read first trace");

    if (!getparfloat("scale", &scale)) {
        if ( !strcmp(key, "offset") ) {
            scale = ( tr.scalco < 0 )? -1.0/tr.scalco : (tr.scalco > 0)? tr.scalco : 1.0;
            if ( verbose > 0) warn(" Offset scaling factor=%.4f", scale);
        } else {
            scale = 1.0;
        }
    }

    float mscale = 0.001; // defauft mute time is given in ms
    int ndata = 1; // number of float for one sample , 1 for seismic, 2 for complex trace
    cmplxtr = tr.trid == FUNPACKNYQ;
    seismic = ISSEISMIC(tr.trid);
    if (seismic) {
        if (!tr.dt) err("dt header field must be set");
    } else if (cmplxtr) { // complex traces, eg. frequency data
        if (!tr.d1) err("  d1 header field must be set for complex traces");
        ndata = 2;
        mscale = 0.1; // mute time given in 0.1 Hz
    } else if (tr.trid==30) {   /* depth section */
        if (!tr.d1) err("d1 header field must be set");
        mscale = 1.0; // mute time is given in meter
    } else {
        err("tr.trid = %d, unsupported trace id", tr.trid);
    }

    /* Loop over traces */
    do {
        int nt     = (int) tr.ns;
        float tmin = (seismic)? tr.delrt/1000.0 : tr.f1;
        float dt = (seismic)? ((double) tr.dt)/1000000.0 : tr.d1;
        float t, t1, t2;
        int nmute;
        int itaper;
        int topmute;
        int botmute;
        register int i;

        ++ntr;

        if ( mode == 3 ) {
            t1 = mscale*((float)tr.muts);
            t2 = mscale*((float)tr.mute);
            if (t1 > t2 && verbose > 1) warn("Mute start %d > %d mute end for %d-th trace cdp=%d fldr=%d offset=%d !",
                    tr.muts, tr.mute, ntr, tr.cdp, tr.fldr, tr.offset);
        } else if ( hdr ) {
            /* get t value of key */
            gethval(&tr, tindex, &tval);
            t = mscale*(vtof(ttype, tval) + bias);
        } else {
            /* get value of key and convert to float */
            gethval(&tr, index, &val);
            fval = vtof(type, val)*scale;

            /* linearly interpolate between (xmute,tmute) values */
            intlin(nxmute, xmute, tmute, tmin, tmute[nxmute-1], 1, &fval, &t);
            if (absolute) fval = abs(fval);
            if ( mode == 2 ) {
                t1=NINT(tm0 + fval/linvel - 0.5*t);
                t2=NINT(tm0 + fval/linvel + 0.5*t);
            }
        }

        /* do the mute */
        if (mode==0) {	/* mute above */
            nmute = NINT((t - tmin)/dt);
            if (nmute > 0) memset( (void *) tr.data, 0, nmute*FSIZE*ndata);
            for (i = 0; i < ntaper; ++i) {
                itaper = ndata*(nmute + i);
                if ( itaper >= 0 && itaper < nt) {
                    tr.data[itaper] *= taper[i];
                    if (cmplxtr) tr.data[itaper + 1] *= taper[i];
                }
            }
            tr.muts = NINT(t/mscale);
        } else if (mode==1){	/* mute below */
            nmute = NINT((tmin + (nt/ndata)*dt - t)/dt);
            if (nmute > 0) memset( (void *) (tr.data+nt-nmute*ndata), 0, nmute*FSIZE*ndata);
            for (i = 0; i < ntaper; ++i) {
                itaper = nt-(nmute-1)*ndata-i;
                if ( itaper >= 0 && itaper < nt) {
                    tr.data[itaper] *= taper[i];
                    if (cmplxtr) tr.data[itaper + 1] *= taper[i];
                }
            }
            tr.mute = NINT(t/mscale);
        } else if (t1 < t2 ) { // mode==2, 3, air wave mute or surgical mute
            topmute=NINT((t1 - tmin)/dt);
            botmute=NINT((t2 - tmin)/dt);
            topmute=MIN(MAX(0, topmute), nt);
            botmute=MIN(nt, botmute);
            nmute=botmute - topmute;
            if (nmute > 0) memset( (void *) (tr.data+topmute*ndata), 0, nmute*FSIZE*ndata);
            for (i = 0; i < ntaper; ++i){
                itaper=(topmute - i)*ndata;
                if (itaper > 0) {
                    tr.data[itaper] *=taper[i];
                    if (cmplxtr) tr.data[itaper + 1] *= taper[i];
                }
            }
            for (i = 0; i < ntaper; ++i){
                itaper=(botmute + i)*ndata;
                if (itaper<nt) {
                    tr.data[itaper] *=taper[i];
                    if (cmplxtr) tr.data[itaper + 1] *= taper[i];
                }
            }

        }
        puttr(&tr);
    } while (gettr(&tr));

    if (xmute) free1float(xmute);
    if (tmute) free1float(tmute);
    if (taper) free1float(taper);

    return(CWP_Exit());
}

float* MakeTaper(int n, int type, int verbose)
{
    int i;
    float env = 0.0, f, x;
    const float min = 0.0, max = 1.0;
    const float EPS = 3.8090232;    /* exp(-EPS*EPS) = 5e-7, "noise" level  */

    if (n < 1) return NULL;

    float* w = ealloc1float(n);

    for (i = 0; i < n; i++) {
        f = (float) (i+1) / n;
        switch (type) {
            case 1: env = min + (max - min) * f;
                break;
            case 2: env = 0.5 * (1.0 - cos(PI * f));
                break;
            case 3: env = sin(PI * f / 2.);
                break;
            case 4: env = sin(PI * f / 2.)*sin(PI * f / 2.);
                break;
            case 5: x = 2.0 * (1 - f);
                env = exp(-(x * x));
                break;
            case 6: x = EPS * (1 - f);
                env = exp(-(x * x));
                break;
            default:err(" unknown taper taper type (=%d)", type);
        }
        w[i] = env;
    }

    if (verbose > 1) {
        fprintf(stderr, "  %d weighting factors for taper: ", n);
        for (i=0; i<n; ++i) {
            if ( !(i%20) ) fprintf(stderr, "\n");
            fprintf(stderr, "%.3f ", w[i]);
        }
        fprintf(stderr, "\n");
    }

    return w;
}

