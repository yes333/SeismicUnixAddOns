/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* SUMEAN: $Revision: 1.3 $ ; $Date: 2003/08/20 18:32:49 $	*/

#include "su.h"
#include "segy.h"

/*************************** self documentation ******************************/
char *sdoc[] = {
"									",
" SUMEAN - calculate the mean of sample values for a time window        ",	
"									",
" sumean < stdin > stdout [optional parameters] 			",
"									",
" Optional parameters: 							",
"									",
"   mode=1              =1 calculate mean                               ",
"                       =2 find maximum sample value                    ",
"									",
"   tmin=0.0            start time of the window			",
"   tmax=(from header)  end time of the window				",
"   itmin=0             min time sample of the window			",
"   itmax=(from header) max time sample of the window			",
"   nt=itmax-itmin+1    number of time samples of the window		",
"									",
"   key=unscale         header key to store the calculated mean		",
"									",
"   power=2.0           mean to the power                               ",
"                         (e.g. =1.0 mean amplitude, =2.0 mean energy)  ",
"   sqrt=1              take square root of the calculated mean         ",
"                         (normally used in conjunction with power=2)   ",
"   reci=0              =1 take reciprocity of mean (1/mean)            ",
"                                                                       ",
"   abs=0               preserve sign if power=1.0                      ",
"                       =1 use absolute value if power=1.0              ",
"									",
"   verbose=0           >0 display value of each trace/section          ",
"									",
" Description:			 					",
" Each sample is raised to the requested power, and the sum of all those",
" values is averaged for the time window.                               ",
"									",
" Notes:			 					",
" Hard zeros (muted samples) are excluded from calculation.             ",
"			 						",
" Examples:                                                             ",
" 1. Normalize trace amplitude to its absolute maximum                  ",
" sumean tmin=0.9 tmax=2.9 mode=2 abs=1 key=unscale < mydata.su |\\",
" suhtmath op=div key=unscale |\\",
" suxwigb perc=99 windowtitle=\"trace normalized to its max. amplitude\" &",
"			 						",
NULL};

/* Credits:
 *  Bjoern E. Rommel, IKU, Petroleumsforskning / October 1997
 *		    bjorn.rommel@iku.sintef.no
 *  Sanyu Ye, RWS, Read Well Service, Oslo, Norway / Feb. 2009
 *	added time window and key to store the calculated mean
 */

/**************** end self doc ***********************************************/


/* Globals */
segy tr;

int main (int argc, char **argv)
{
    int verbose;	   /* flag for printing extra information	*/
    int itr;               /* trace number				*/
    int it;		   /* sample number				*/
    int nt;		   /* number of time samples of the time window */
    int itmin;             /* smallest time sample (zero-based)    */
    int itmax;             /* largest time sample (zero-based)     */
    float tmin;            /* minimum time to calculate        	*/
    float tmax;            /* maximum time to calculate		*/
    float dt;              /* sampling interval, from tr.dt for seismic data (sec) */
    float power;           /* mean to the power of			*/
    double tmean;	   /* average mean of trace			*/
    double gmean;	   /* average mean of section			*/
    int sqroot;            /* square root flag */
    int reci;              /* reciprocity flag */
    int ab;                /* absolute value flag */
    int mode;              /* absolute value flag */

    cwp_String key;	/* key storing calculated mean		*/
    cwp_String type;	/* type for key	*/
    int index;		/* indexes for key 	*/

    Value val;		/* value of key		*/

    /* Initialize */
    initargs (argc, argv);
    requestdoc (0);

    /* Get optional parameters */
    if (!getparint ("verbose", &verbose))   verbose = 0;
    if (!getparint ("mode", &mode))     mode = 1;
    if (!getparint ("sqrt", &sqroot))   sqroot = 1;
    if (!getparint ("reci", &reci))     reci = 0;
    if (!getparint ("abs", &ab))        ab = 0;
    if (!getparfloat ("power", &power))   power = 2.0;

    if (!getparstring ("key", &key))   key = "unscale" ;
    type  = hdtype(key);
    index = getindex(key);

    /* Get info from first trace */
    if (!gettr(&tr))  err ("can't get first trace");

    /* Time gating parameters */
    dt = ((double) tr.dt)/1000000.0;
    if (getparint("itmin", &itmin)) {
            tmin = itmin*dt;
    } else if (getparfloat("tmin", &tmin)) {
            itmin = NINT(tmin / dt);
    } else {
            itmin = 0;
            tmin = 0;
    }
    if (getparint("itmax", &itmax)) {
            tmax = itmax*dt;
            nt = itmax - itmin + 1;
    } else if (getparfloat("tmax", &tmax)) {
            itmax = NINT(tmax / dt);
            nt = itmax - itmin + 1;
    } else if (getparint("nt", &nt)) {
            itmax = itmin + nt - 1;
            tmax = itmax*dt;
    } else {
            itmax = tr.ns - 1;
            tmax = itmax*dt;
            nt = itmax - itmin + 1;
    }


    /* Check time gating values */
    if (itmin < 0)
            err("itmin=%d should be positive", itmin);
    if (nt > SU_NFLTS)
            err("nt=%d exceeds SU_NFLTS=%d", nt, SU_NFLTS);
    if (itmin > itmax)
            err("itmin=%d, itmax=%d conflict", itmin, itmax);
    if (tr.ns <= itmax)
            err("tr.ns=%d, itmax=%d window cannot extend over the trace length", tr.ns, itmax);

    /* Initialize mean value of section */
    gmean = 0.0;

    itr = 0;
    /* Loop through traces */
    do {
        tmean = 0.0;   /* reset mean value of trace */

        /* Loop through samples */
        for (it = itmin, nt = 0; it <= itmax; ++it) {
            if (tr.data[it] == 0.0) continue;  // skip hard zero

            if (mode == 2) {
                float v = tr.data[it]; 
                if( ab ) v = ABS(v);
                if (v > tmean) tmean = v;
            } else {
                /* Raise sample to the requested power, add to mean value of trace */
                if( ab ){
                   tmean += pow (fabs (tr.data[it]), power);
                }else{
                   tmean += pow ((tr.data[it]), power);
                }
                ++nt;
            }
        }

        if (mode != 2) {
            /* Average mean value of trace */
            tmean = (nt == 0)? 0 : tmean / nt;
            if (sqroot) tmean = sqrt(fabs(tmean));
        }    
        if (reci) tmean = (tmean == 0)? 0 : 1.0/tmean;

        val.f = (float) tmean;
        puthval(&tr, index, &val);
        puttr(&tr);
        
        if ( mode == 2 ) {
            if (tmean > gmean) gmean = tmean; 
        } else {
            /* Add to the mean value of section */
            gmean += tmean;
        }

        itr++ ;	/* count traces */
        /* Print mean value of trace */
        if (verbose > 0)   
            warn("trace: %i  %s value: %e", itr, (mode == 2)? "max" : "mean", tmean);

  } while (gettr(&tr)); 

  /* Average mean value of section */
  if ( mode != 2) gmean = gmean / itr;

  /* Print mean value of section */
  if (verbose > 0 ) warn("global %s: %e\n", (mode == 2)? "max" : "mean", gmean);
  
  return(CWP_Exit());
}
