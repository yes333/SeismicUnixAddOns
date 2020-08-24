/* Copyright (c) Read Well Services, 2008-2011.*/
/* All rights reserved.                       */

/* SPWFS */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"                                                                               ",
" SPWFS - Wave Field Seperation in Tau-P domain for OBC seismic data            ",
"                                                                               ",
"    spwfs <stdin >stdout [options]                                             ",
"                                                                               ",
" Required Parameters:                                                          ",
"    none                                                                       ",
"                                                                               ",
" Optional Parameters:                                                          ",
"    nc=3           =4 four component wave field separation (P Z X Y)           ",
"                   =3 three component (P Z X)                                  ",
"                   =2 two component (P Z)                                      ",
"                                                                               ",
"    modes=-2       codes designate output wave mode                            ",
"                   =-1 P   upgoing pressure on P (hyddrophone)                 ",
"                   =1  P downgoing pressure on P (hyddrophone)                 ",
"                   =-2 P   upgoing velocity on Z                               ",
"                   =2  P downgoing velocity on Z                               ",
"                   =-3 P   upgoing velocity on X                               ",
"                   =3  P downgoing velocity on X                               ",
"                   =-4 P   upgoing velocity on Y                               ",
"                   =4  P downgoing velocity on Y                               ",
"                   =-5 S   upgoing pressure on P (hyddrophone)                 ",
"                   =5  S downgoing pressure on P (hyddrophone)                 ",
"                   =-6 S   upgoing velocity on Z                               ",
"                   =6  S downgoing velocity on Z                               ",
"                   =-7 S   upgoing velocity on X                               ",
"                   =7  S downgoing velocity on X                               ",
"                   =-8 S   upgoing velocity on Y                               ",
"                   =8  S downgoing velocity on Y                               ",
"                                                                               ",
"    subtract=0     =1 compute all complimentary wave modes                     ",
"                                                                	",
"    key=cdpt       keyword used to match input reflection coefficient 	",
"                                                                	",
"    vp=1.6         P-velocity of seafloor [km/s]			",
"    vs=0.2         S-velocity of seafloor [km/s]			",
"    db=1.6         density of seafloor [g/cm**2]                       ",
"    vw=1.48        P-velocity in sea water [km/s]                      ",
"    dw=1.03        density of sea water [g/cm**2]                      ",
"                   the values above give a default Rpp = 0.255        	",
"                                                                	",
"    matched=1      geophone data is matched                          	",
"                   =0 not matched, but calibrated                      ",
"                     (geophone to be backscaled by impedance)	        ",
"                                                                	",
"    nrcv=240       max. number of receivers                            ",
"    rcfile=        file containing reflection coefficients         	",
"                                                                	",
"    verbose=0      >0 output info                                      ",
"                                                                	",
" Note:                                                                 ",
"  Input tau-p data must be sorted trace by trace in order P, Z, X, Y.  ",
"    Their polarities should be as recorded, e.g. >>>> for shots recorded",
"    in first quadrant (X --> E, Y --> N, Z --> Down). In other words the",
"    polarities tau-p transform of P, Z, X, Y should be the same like >>>>",
"    in first quadrant (px>0 & py>0).                                   ",
"  Output data are sorted trace by trace in the order as specified by   ",
"    the modes. keyword duse code is set to be mode number              ",
"    duse=mode                                                          ",
"                                                                	",
" Version 2.0.0  last modified May 2011 by Sanyu Ye                     ",
"                                                                	",
NULL};

/* Credits:
 *	Sanyu Ye: READ Well Service, Feb. 2008
 */
/**************** end self doc *******************************************/

/* define output wave mode */
#define MODE_P_TOTOL_UPGOING            1
#define MODE_P_TOTOL_DOWNGOING          2
#define MODE_S_TOTOL_UPGOING            3
#define MODE_S_TOTOL_DOWNGOING          4
#define MODE_P_POTENTIAL_UPGOING        5
#define MODE_P_POTENTIAL_DOWNGOING      6
#define MODE_S_POTENTIAL_UPGOING        7
#define MODE_S_POTENTIAL_DOWNGOING      8
#define MODE_P_V_Z_UPGOING              9
#define MODE_P_V_Z_DOWNGOING           10
#define MODE_S_V_H_UPGOING             11
#define MODE_S_V_H_DOWNGOING           12
#define MODE_S_V_V_UPGOING             13
#define MODE_S_S_V_UPGOING             14
#define MODE_P_UPGOING_2C              15
#define MODE_P_UPGOING_4C              16
 
/* prototypes */
float InferV(float Vw, float Dw, float Vp0, float D0, float Vs0, float Rpp, float* Vp, float* Vs);
int GetCMatrix(int subtract, int matched, float Vp, float Vs, float Imp, float px, float py, float**** c4d);

#define  MaxModes 8

int main(int argc, char **argv) {
    cwp_String key;		/* header key word from segy.h		*/
    cwp_String type;            /* type of key				*/
    Value val;                  /* value of key			*/
    double dval, dvallast;      /* double value of key			*/
    int index;                  /* index of key				*/
    int nt;                     /* numsamps as int			*/
    int nc;                     /* number of components			*/
    int nsegy;                  /* length bytes read for segy trace	*/
    int nm, modes[12];          /* number of output wave modes and its values */
    segyhdr *hdrs1d = NULL;     /* buffer array of input trace headers  */
    float *trbuf=NULL;          /* trace buffer	containing trace data, padded with zero to ntfft  */
    float **trdata=NULL;	/* tau-p-domain for all input components   */
    float **trout=NULL;         /* tau-p-domain traces of all output components   */
    float Vw, dw;               /* water velocity and density           */
    float Vp0, Vs0, d0;         /* default P-, S-velocity and density of sea bottom */
    float Vp, Vs, Rpp, Imp;     /* P-, S-velocity, RC and impedance at receiver site */
    float px, py    ;           /* input trace px py header values */

    int verbose;		/* =0 no info, >0 output info  */
    int subtract;               /* compute complimentary modes	*/
    int matched=1;              /* is geophone data mached?	*/
    int i, j, m, ntr;                 /* loop and trace counter */
    int nrcv, nRow;             /* max. and actual rc values read in */
    double *mkeys = NULL;           /* pointer containing values of matching keys */
    double *vkeys = NULL;           /* pointer containing values of keys to be set */

    char *rcfile="";                /* name of input file of reflection coefficients	*/
    FILE *rcfp=NULL;                /* pointer to input file		*/
    cwp_Bool from_file=cwp_false;   /* is the data from infile?	*/
    cwp_Bool RppFound=cwp_false;    /* is reflectivity value found for that receiver?	*/
    
    segy tr;
    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* get match key */
    if (!getparstring("key", &key))	 key="cdpt";
    type = hdtype(key);
    index = getindex(key);

    if (!getparint("nc", &nc)) nc=3;
 
    /* get output wave modes */
    if ((nm=countparval("modes")) > 0) {
        getparint("modes", modes);         
        for (i=0; i<nm; ++i) {
             if (modes[i] > MaxModes || modes[i] < -MaxModes) err("invalid wave mode modes[%d]=%d", i, modes[i]);
        }
    } else { /* set default */
        nm = 1;
        modes[0] = -3;
    }
    
    if (!getparint("verbose", &verbose))    verbose=0;
    if (!getparint("matched", &matched))    matched=1;
    if (!getparint("subtract", &subtract))  subtract=0;

    if (!getparfloat("vw", &Vw)) Vw=1.48;
    if (!getparfloat("dw", &dw)) dw=1.03;
    if (!getparfloat("vp", &Vp0)) Vp0=1.60;
    if (!getparfloat("vs", &Vs0)) Vs0=0.20;
    if (!getparfloat("db", &d0)) d0=1.60;

    /* get name of infile for reflection coefficients*/
    getparstring("rcfile",&rcfile);

    /* if infile is specified get specified keys from file */
    if (*rcfile!='\0') {
            /* open infile */
            if((rcfp=efopen(rcfile,"r"))==NULL) err("cannot open rcfile=%s\n",rcfile);
            from_file=cwp_true; /* set from_file flag */
    }
   
    
    /* Evaluate time bounds from getpars and first header */
    if ( (nsegy=gettr(&tr)) == 0) err("can't get first trace");

    nt = tr.ns;

    hdrs1d = (segyhdr *) ealloc1(nc, HDRBYTES);
    trdata = ealloc2float(nt, nc);
    trout = ealloc2float(nt, nm);
    // matrix coifficients c[l][k][j][i], c[modes][subtract][-/+][column]
    float**** c = ealloc4float(4, 2, 2, MaxModes);

    /* read rc values from file */
    if (from_file) { 
        if (!getparint("nrcv", &nrcv)) nrcv=240;
        mkeys = ealloc1double(nrcv);
        vkeys = ealloc1double(nrcv);

        /* reading all data from ascii infile */
        /*  fscanf returns: 0   : characters there, but no conversion (error)
        *		  EOF : eof before conversion
        *		  else: number of conversions 
        */
        for (nRow=0; nRow < nrcv; ++nRow) {
            int ret = fscanf(rcfp, "%lf", &mkeys[nRow]);
            if(ret == EOF || ret == 0) break;  /* else everything is okay: get out of the loop */
            ret = fscanf(rcfp, "%lf", &vkeys[nRow] );
            if(ret == EOF || ret == 0) break;  /* else everything is okay: get out of the loop */
        }
        
        if (verbose) warn(" %d values are read in from rcfile=%s", nRow, rcfile);
    }
    
    //decomposing wavefield
    ntr = 0;
    dvallast = 0.0;
    /* Main loop over traces */
    while (nsegy > HDRBYTES ) {
        i = ntr % nc;
        memcpy(&hdrs1d[i], &tr, HDRBYTES);  //cache header
        memcpy(trdata[i], tr.data, nt*FSIZE); // cache data

        if ( ++ntr%nc == 0 ) { // nc traces read in, start separation
            if (from_file) { // retrieve rc
                gethval(&tr, index, &val);
                dval = vtod(type, val);
                if (dvallast != dval) {
                    RppFound = cwp_false;
                    for (i=0; i<nRow; ++i) { //loop over to find right rc 
                        if ( dval == mkeys[i] ) { 
                            Rpp = vkeys[i];
                            RppFound = cwp_true;
                            dvallast = dval;
                            break;
                        }
                    }
                    if (RppFound) { 
                        Imp = InferV(Vw, dw, Vp0, d0, Vs0, Rpp, &Vp, &Vs);
                    } else { //assume default values
                        Vp = Vp0;
                        Vs = Vs0;
                        Imp = d0*Vp;
                    }
                }
            } else { //assume default values
                Vp = Vp0;
                Vs = Vs0;
                Imp = d0*Vp;
            }
            
            px = tr.fx;
            py = tr.fy;
            
            //Decompose(nt, matched, Vp, Vs, Imp, px, py, nc, trdata, nm, modes, trout);

            memset(***c, 0, MaxModes*2*2*4*FSIZE);
            GetCMatrix(subtract, matched, Vp, Vs, Imp, px, py, c);
            for (m = 0; m < nm; ++m) {
                int mode = ABS(modes[m]) - 1;
                int updown = (modes[m] < 0) ? 0 : 1;
                for (i = 0; i < nt; ++i) { // loop over samples
                    for (trout[m][i] = 0, j = 0; j < nc; ++j) { // loop over input components
                        trout[m][i] += c[mode][subtract][updown][j] * trdata[j][i];
                    }
                }
            }
            
            // back to t and output
            for (i=0; i<nm; ++i) {
                memcpy(&tr, &hdrs1d[(i<nc)? i : nc-1], HDRBYTES);
                tr.duse = modes[i];  // mark component with mode
                memcpy(tr.data, trout[i], nt*FSIZE);
                puttr(&tr);
            }
            if (verbose && (ntr/nc - 1)%(1000/verbose) == 0) { //output info
                gethval(&tr, index, &val);
                warn(" %d-th trace: %s=%d  px=%.6f py=%.6f of %d compomnents are processed", 
                        ntr/nc, key, vtoi(type, val), px, py, nc);
                // print out coefficients for inspection
                for (i=0; i<nm; ++i) {
                    int mode = ABS(modes[i]) - 1;
                    int updown = (modes[i] < 0) ? 0 : 1;
                    fprintf(stderr, "        Mode=%2d  P=%8.6f Z=%8.6f X=%8.6f Y=%8.6f \n",
                        modes[i], c[mode][subtract][updown][0],c[mode][subtract][updown][1],c[mode][subtract][updown][2],c[mode][subtract][updown][3]);
                }
            }
        }

        //read next trace
        nsegy = gettr(&tr);
    }

    if (verbose) warn(" Totally %d traces for each of %d compomnents are processed", ntr/nc, nc);
    
    if (hdrs1d) free(hdrs1d);
    if (trbuf) free1float(trbuf);
    if (trdata) free2float(trdata);
    if (trout) free2float(trout);
    
    return(CWP_Exit());
}
    
float InferV(float Vw, float Dw, float Vp0, float D0, float Vs0, float Rpp, float* Vp, float* Vs)
{
    // assuming percentage change of density is n1 times than that (a) of Vp , change (b) of Vs is n2 times than Vp
    // upon percentage change c in reflection coefficient, there is formula
    // a = 2 * c * R0/ ( ( 1 + R0) * ( 1 - R0 ) * ( 1 + n1 ) )  ~= 2cR0/(1+n1);    b = n2 * a;
    
    float n1 = 3.0, n2 = 5.0; // assumption
    float a, b, c, R0;
    
    R0 = (Vp0 * D0 - Vw * Dw) / (Vp0 * D0 + Vw * Dw);
    c = Rpp / R0 - 1.0;
    a = 2 * c * R0 / ( ( 1 + R0) * ( 1 - R0 ) * ( 1 + n1 ) );
    b = n2 * a;
    *Vp = (1.0 + a) * Vp0;
    *Vs = (1.0 + b) * Vs0;
    return Dw*Vw*(1.0 + Rpp)/(1.0 - Rpp);  // impedance of seafloor
}

int GetCMatrix(int subtract, int matched, float Vp, float Vs, float Imp, float px, float py, float**** c4d)
{
    int i, j;  // loop counter
    float a, b;
    float p, VpQp, divVpQp, VsQs, divVsQs;
    float COStheta2;
    double minCOStheta2 = 0.04; // about 78 degrees
    cwp_Bool overcritical = cwp_false;
    float v11,v12,v13,v14,v21,v22,v23,v24,v31,v32,v33,v34;  // voif for P without sign
    float v41,v42,v43,v44,v51,v52,v53,v54,v61,v62,v63,v64;  // for SV
    float t31,t32,t33,t34,t61,t62,t63,t64;

    p = sqrt(px*px + py*py);

    //Cap the minimum of V*Q = cos(theta) thus division not go to unlimited
    COStheta2 = 1.0 - p*p*Vp*Vp;
    if (COStheta2 < -minCOStheta2  ) {
        VpQp = sqrt(-COStheta2);
        divVpQp = 1.0/VpQp;
        overcritical = cwp_true;
    } else if(COStheta2 < 0.0) {
        VpQp = sqrt(-COStheta2);
        divVpQp = 1.0/sqrt(minCOStheta2);
        overcritical = cwp_true;
    } else if (COStheta2 <= minCOStheta2) {
        VpQp = sqrt(COStheta2);
        divVpQp = 1.0/sqrt(minCOStheta2);
    } else {
        VpQp = sqrt(COStheta2);
        divVpQp = 1.0/VpQp;
    }
    COStheta2 = 1.0 - p*p*Vs*Vs;
    if (COStheta2 < minCOStheta2  ) {
        VsQs = -sqrt(-COStheta2);
        divVsQs = 1.0/VsQs;
    } else if(COStheta2 < 0.0) {
        VsQs = -sqrt(-COStheta2);
        divVsQs = 1.0/sqrt(minCOStheta2);
    } else if (COStheta2 <= minCOStheta2) {
        VsQs = sqrt(COStheta2);
        divVsQs = 1.0/sqrt(minCOStheta2);
    } else {
        VsQs = sqrt(COStheta2);
        divVsQs = 1.0/VsQs;
    }

    a = p*p*Vs*Vs;
    b = 1.0 - 2.0*p*p*Vs*Vs;

    //c = (!p)?  0.0: px/p;
    //d = (!p)?  0.0: py/p;

    if (!subtract) {
        v11 = px*px*Vs*Vs;
        v12 = px*py*Vs*Vs;
        v13 = px*b*Vp*divVpQp/2;
        v14 = px*Vp/2;
        v21 = v12;
        v22 = py*py*Vs*Vs;
        v23 = py*b*Vp*divVpQp/2;
        v24 = py*Vp/2;
        v31 = px*Vs*Vs*VpQp/Vp;
        v32 = py*Vs*Vs*VpQp/Vp;
        v33 = b/2;
        v34 = VpQp/2;
        v41 = (1 - 2*px*px*Vs*Vs)/2;
        v42 = px*py*Vs*Vs;
        v43 = px*Vs*VsQs;
        v44 = px*Vp/2;
        v51 = v42;
        v52 = (1 - 2*py*py*Vs*Vs)/2;
        v53 = py*Vs*VsQs;
        v54 = py*Vp/2;
        v61 = px*b*Vs*divVsQs/2;
        v62 = py*b*Vs*divVsQs/2;
        v63 = a;
        v64 = p*p*Vp*Vs*divVsQs/2;

        // nagative sign for S3=-P, or P=-S3
        t31 = -px*Vs*Vs*b/Vp;
        t32 = -py*Vs*Vs*b/Vp;
        t33 = -b*b*divVpQp/2;
        t34 = -b/2;
        t61 = t31;
        t62 = t32;
        t63 = -2*a*Vs*VsQs/Vp;
        t64 = -a;

        // c4d[modes][subtract][-/+][column]
        // up / down P Vz
        c4d[1][0][0][0] = -v34;
        c4d[1][0][1][0] =  v34;
        c4d[1][0][0][1] =  v33;
        c4d[1][0][1][1] =  v33;
        c4d[1][0][0][2] = -v31;
        c4d[1][0][1][2] =  v31;
        c4d[1][0][0][3] = -v32;
        c4d[1][0][1][3] =  v32;
        // up / down P Vx
        c4d[2][0][0][0] =  v14;
        c4d[2][0][1][0] =  v14;
        c4d[2][0][0][1] = -v13;
        c4d[2][0][1][1] =  v13;
        c4d[2][0][0][2] =  v11;
        c4d[2][0][1][2] =  v11;
        c4d[2][0][0][3] =  v12;
        c4d[2][0][1][3] =  v12;
        // up / down P Vy
        c4d[3][0][0][0] =  v24;
        c4d[3][0][1][0] =  v24;
        c4d[3][0][0][1] = -v23;
        c4d[3][0][1][1] =  v23;
        c4d[3][0][0][2] =  v21;
        c4d[3][0][1][2] =  v21;
        c4d[3][0][0][3] =  v22;
        c4d[3][0][1][3] =  v22;

        // up / down SV Vz
        c4d[5][0][0][0] = -v64;
        c4d[5][0][1][0] =  v64;
        c4d[5][0][0][1] =  v63;
        c4d[5][0][1][1] =  v63;
        c4d[5][0][0][2] =  v61;
        c4d[5][0][1][2] = -v61;
        c4d[5][0][0][3] =  v62;
        c4d[5][0][1][3] = -v62;
        // up / down SV Vx
        c4d[6][0][0][0] = -v44;
        c4d[6][0][1][0] = -v44;
        c4d[6][0][0][1] =  v43;
        c4d[6][0][1][1] = -v43;
        c4d[6][0][0][2] =  v41;
        c4d[6][0][1][2] =  v41;
        c4d[6][0][0][3] = -v42;
        c4d[6][0][1][3] = -v42;
        // up / down SV Vy
        c4d[7][0][0][0] = -v54;
        c4d[7][0][1][0] = -v54;
        c4d[7][0][0][1] =  v53;
        c4d[7][0][1][1] = -v53;
        c4d[7][0][0][2] = -v51;
        c4d[7][0][1][2] = -v51;
        c4d[7][0][0][3] =  v52;
        c4d[7][0][1][3] =  v52;

        // up / down P pressure on Z
        c4d[0][0][0][0] = -t34;
        c4d[0][0][1][0] = -t34;
        c4d[0][0][0][1] =  t33;
        c4d[0][0][1][1] = -t33;
        c4d[0][0][0][2] = -t31;
        c4d[0][0][1][2] = -t31;
        c4d[0][0][0][3] = -t32;
        c4d[0][0][1][3] = -t32;

        // up / down S pressure on Z
        c4d[4][0][0][0] = -t64;
        c4d[4][0][1][0] = -t64;
        c4d[4][0][0][1] =  t63;
        c4d[4][0][1][1] = -t63;
        c4d[4][0][0][2] =  t61;
        c4d[4][0][1][2] =  t61;
        c4d[4][0][0][3] =  t62;
        c4d[4][0][1][3] =  t62;


    } else {
        v11 = (1 - px*px*Vs*Vs);
        v12 = px*py*Vs*Vs;
        v13 = px*b*Vp*divVpQp/2;
        v14 = px*Vp/2;
        v21 = v12;
        v22 = (1 - py*py*Vs*Vs);
        v23 = py*b*Vp*divVpQp/2;
        v24 = py*Vp/2;
        v31 = px*Vs*Vs*VpQp/Vp;
        v32 = py*Vs*Vs*VpQp/Vp;
        v33 = (1 + 2*a)/2;
        v34 = VpQp/2;
        v41 = (1 + 2*px*px*Vs*Vs)/2;
        v42 = px*py*Vs*Vs;
        v43 = px*Vs*VsQs;
        v44 = px*Vp/2;
        v51 = v42;
        v52 = (1 + 2*py*py*Vs*Vs)/2;
        v53 = py*Vs*VsQs;
        v54 = py*Vp/2;
        v61 = px*b*Vs*divVsQs/2;
        v62 = py*b*Vs*divVsQs/2;
        v63 = (1 - a);
        v64 = p*p*Vp*Vs*divVsQs/2;

        // nagative sign for S3=-P, or P=-S3
        t31 = -px*Vs*Vs*b/Vp;
        t32 = -py*Vs*Vs*b/Vp;
        t33 = -b*b*divVpQp/2;
        t34 = -(1 + 2*a)/2;
        t61 = t31;
        t62 = t32;
        t63 = -2*p*p*Vs*Vs*Vs*VsQs/Vp;
        t64 = -(1 - a);

        // c4d[modes][subtract][-/+][column]
        // up / down P Vz
        c4d[1][1][0][0] =  v34;
        c4d[1][1][1][0] = -v34;
        c4d[1][1][0][1] =  v33;
        c4d[1][1][1][1] =  v33;
        c4d[1][1][0][2] =  v31;
        c4d[1][1][1][2] = -v31;
        c4d[1][1][0][3] =  v32;
        c4d[1][1][1][3] = -v32;
        // up / down P Vx
        c4d[2][1][0][0] = -v14;
        c4d[2][1][1][0] = -v14;
        c4d[2][1][0][1] =  v13;
        c4d[2][1][1][1] = -v13;
        c4d[2][1][0][2] =  v11;
        c4d[2][1][1][2] =  v11;
        c4d[2][1][0][3] = -v12;
        c4d[2][1][1][3] = -v12;
        // up / down P Vy
        c4d[3][1][0][0] = -v24;
        c4d[3][1][1][0] = -v24;
        c4d[3][1][0][1] =  v23;
        c4d[3][1][1][1] = -v23;
        c4d[3][1][0][2] = -v21;
        c4d[3][1][1][2] = -v21;
        c4d[3][1][0][3] =  v22;
        c4d[3][1][1][3] =  v22;
        // up / down SV Vz
        c4d[5][1][0][0] =  v64;
        c4d[5][1][1][0] = -v64;
        c4d[5][1][0][1] =  v63;
        c4d[5][1][1][1] =  v63;
        c4d[5][1][0][2] = -v61;
        c4d[5][1][1][2] = +v61;
        c4d[5][1][0][3] = -v62;
        c4d[5][1][1][3] = +v62;
        // up / down SV Vx
        c4d[6][1][0][0] =  v44;
        c4d[6][1][1][0] =  v44;
        c4d[6][1][0][1] = -v43;
        c4d[6][1][1][1] =  v43;
        c4d[6][1][0][2] =  v41;
        c4d[6][1][1][2] =  v41;
        c4d[6][1][0][3] =  v42;
        c4d[6][1][1][3] =  v42;
        // up / down SV Vy
        c4d[7][1][0][0] =  v54;
        c4d[7][1][1][0] =  v54;
        c4d[7][1][0][1] = -v53;
        c4d[7][1][1][1] =  v53;
        c4d[7][1][0][2] =  v51;
        c4d[7][1][1][2] =  v51;
        c4d[7][1][0][3] =  v52;
        c4d[7][1][1][3] =  v52;

        // up / down P pressure on Z
        c4d[0][1][0][0] = -t34;
        c4d[0][1][1][0] = -t34;
        c4d[0][1][0][1] = -t33;
        c4d[0][1][1][1] =  t33;
        c4d[0][1][0][2] =  t31;
        c4d[0][1][1][2] =  t31;
        c4d[0][1][0][3] =  t32;
        c4d[0][1][1][3] =  t32;

        // up / down S pressure on Z
        c4d[4][1][0][0] = -t64;
        c4d[4][1][1][0] = -t64;
        c4d[4][1][0][1] = -t63;
        c4d[4][1][1][1] =  t63;
        c4d[4][1][0][2] = -t61;
        c4d[4][1][1][2] = -t61;
        c4d[4][1][0][3] = -t62;
        c4d[4][1][1][3] = -t62;
    }

    if (!matched) {
        for (i=0; i<MaxModes; ++i) {
            for (j=1; j<=3; ++j) {
                c4d[i][0][0][j] *= Imp;
                c4d[i][0][1][j] *= Imp;
            }
        }
    }

    return 0;
}

