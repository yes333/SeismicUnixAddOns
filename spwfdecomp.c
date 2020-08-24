/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPWFDECOMP v. 2.0 */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"									",
" SPWFDECOMP - Wave Field Decomposition in F-P domain for OBC data      ",
"									",
"    spwfdecomp <stdin >stdout [options]				",
"									",
" Required Parameters:                                                  ",
"    none								",
"									",
" Optional Parameters:							",
"    nc=3           =4 four component wave field separation (P Z X Y)   ",
"                   =3 three component (P Z X)                          ",
"                   =2 two component (P Z)                              ",
"                                                                	",
"    modes=1,...    codes designate output wave mode             	",
"                   1 = P total velocity field upgoing          	",
"                   2 = P total velocity field downgoing           	",
"                   3 = S total velocity field upgoing           	",
"                   4 = S total velocity field downgoing           	",
"                   5 = P potential field upgoing                       ",
"                   6 = P potential field downgoing                     ",
"                   7 = S potential field upgoing                       ",
"                   8 = S potential field downgoing                     ",
"                   9 = P vertical velocity field upgoing           	",
"                  10 = P vertical velocity field downgoing           	",
"                  11 = S horizontal velocity field upgoing           	",
"                  12 = S horizontal velocity field downgoing          	",
"                                                                	",
"    matched=1      geophone data is matched                          	",
"                   =0 not matched, but calibrated                      ",
"                     (geophone to be backscaled by impedance)	        ",
"                                                                	",
"    rckey=         keyword where reflection coefficient (RC) is stored	",
"    scale=0.001    scaling factor to be applied on RC                  ",
"                                                                	",
"    vp=1.6         P-velocity of seafloor [km/s]			",
"    vs=0.2         S-velocity of seafloor [km/s]			",
"    db=1.6         density of seafloor [g/cm**2]                       ",
"    vw=1.48        P-velocity in sea water [km/s]                      ",
"    dw=1.03        density of sea water [g/cm**2]                      ",
"                   the values above give a default Rpp = 0.255        	",
"                                                                	",
"    verbose=0      =n output info for every 1000/n-th trace            ",
"                                                                	",
" Note:                                                                 ",
"    because decomposition is performed in F-P domain, input data must  ",
"    be complex frequency (spfk2taup ftau=0 ...).                       ",
"    input data must be sorted trace by trace in order P, Z, X, Y    	",
"    only first trace (P-Hydrophone) need to have RC stored in header   ",
"    output data are sorted trace by trace in the order as specified   	",
"    by the modes. keyword duse is set to be mode number duse=mode	",
"                                                                       ",
" Version 1.0 last modified Oct. 2010 by Sanyu Ye                       ",
NULL
};

/* Credits:
 *	Sanyu Ye: READ Well Service, April. 2010
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

#define  MaxModes 12


/* prototypes */
float InferV(float Vw, float Dw, float Vp0, float D0, float Vs0, float Rpp, float* Vp, float* Vs);

void  wfs_(int* nf, float* df, float* fmin, float* fmax, int* IsMatched, int* IsV, float* Vp, float* Vs,
      float* Imp, float* px, float* py, int* nc, void** trdata, int* nm, int* modes, void** trout, void** coeff);

int verbose; /* =0 no info, >0 output info  */

int main(int argc, char **argv) {
    cwp_String key; /* header key word from segy.h		*/
    cwp_String type; /* type of key				*/
    Value val; /* value of key			*/
    int index; /* index of key				*/
    int nt; /* numsamps as int			*/
    int nc; /* number of components			*/
    int nsegy; /* length bytes read for segy trace	*/
    int nm, modes[12]; /* number of output wave modes and its values */
    int nf; /* number of frequencies */
    segyhdr *hdrs1d = NULL; /* buffer array of input trace headers  */
    float Vw, dw; /* water velocity and density           */
    float Vp0, Vs0, d0; /* default P-, S-velocity and density of sea bottom */
    float Vp, Vs, Rpp, Imp; /* P-, S-velocity, RC and impedance at receiver site */
    float px, py; /* input trace px py header values */
    float dt, fmin, fmax, df, scale; /*  */

    int IsMatched = 1; /* is geophone data mached?	*/
    int IsVSensor= 0; /* is geophone velocity sensor?	*/
    int i, ntr; /* loop and trace counter */
    int verbose; /* =0 no info, >0 output info  */

    cwp_Bool from_hdr = cwp_false; /* is RC data stored in header key?	*/

    segy tr;

    /* Initialize */
    initargs(argc, argv);
    requestdoc(1);

    /* get first trace */
    if ((nsegy = gettr(&tr)) == 0) err("can't get first trace");
    if(tr.trid != FUNPACKNYQ) // input always f-p
        err("input not complex freq data, trid=%d (!=%d)", tr.trid, FUNPACKNYQ);

    if (!getparint("nc", &nc)) nc = 3;

    /* get output wave modes */
    if ((nm = countparval("modes")) > 0) {
        getparint("modes", modes);
        for (i = 0; i < nm; ++i) {
            if (modes[i] > MaxModes || modes[i] < 1) err("invalid wave mode (modes[%d]=%d)", i, modes[i]);
        }
    } else { /* set default */
        nm = 1;
        modes[0] = 1;
    }

    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparint("matched", &IsMatched)) IsMatched = 1;
    if (!getparint("vsensor", &IsVSensor)) IsVSensor = 1;

    if (!getparfloat("vw", &Vw)) Vw = 1.48;
    if (!getparfloat("dw", &dw)) dw = 1.03;
    if (!getparfloat("vp", &Vp0)) Vp0 = 1.60;
    if (!getparfloat("vs", &Vs0)) Vs0 = 0.20;
    if (!getparfloat("db", &d0)) d0 = 1.60;

    /* get key where reflection coifficiet is stored*/
    if (getparstring("rckey", &key)) {
        from_hdr = cwp_true;
        type = hdtype(key);
        index = getindex(key);
    }
    if (!getparfloat("scale", &scale)) scale = 0.001;

    if (!getparfloat("fmin", &fmin)) fmin = 3.0;
    if (!getparfloat("fmax", &fmax)) fmax = 0.5/dt;

    df = tr.d1;
    nt = tr.ns;
    nf = nt/2;

    hdrs1d = (segyhdr *) ealloc1(nc, HDRBYTES);
    complex** coeff = ealloc2complex(4, MaxModes);
    memset(*coeff, 0, 2*4*MaxModes*FSIZE);
    complex** trdata = ealloc2complex(nf, 4);
    memset(*trdata, 0, 2*4*nf*FSIZE);
    complex** trout = ealloc2complex(nf, nm);
    memset(*trout, 0, 2*nm*nf*FSIZE);

    //decomposing wavefield
    ntr = 0;
    /* Main loop over traces */
    while (nsegy > HDRBYTES) {
        i = ntr % nc;
        //transform to f and cache
        memcpy(&hdrs1d[i], &tr, HDRBYTES); //cache header
        memcpy(trdata[i], tr.data, nt * FSIZE); // cache data

        if (++ntr % nc == 0) { // nc traces read in, start separation
            if (from_hdr) { // retrieve rc
                gethval((segy*) &hdrs1d[0], index, &val);  // fetch RC from first trace (H)
                Rpp = scale*vtod(type, val);
                Imp = InferV(Vw, dw, Vp0, d0, Vs0, Rpp, &Vp, &Vs);
            } else { //assume default values
                Vp = Vp0;
                Vs = Vs0;
                Imp = d0*Vp;
            }

            px = tr.fx;
            py = tr.fy;

            // do seperation
            //Decompose(nt, IsMatched, Vp, Vs, Imp, px, py, nc, trdata, nm, modes, trout);
            wfs_(&nf, &df, &fmin, &fmax, &IsMatched, &IsVSensor, &Vp, &Vs, &Imp, &px, &py, 
                    &nc, &trdata[0][0], &nm, &modes[0], &trout[0][0], &coeff[0][0]);

            // output
            for (i = 0; i < nm; ++i) {
                memcpy(&tr, &hdrs1d[(i < nc) ? i : nc - 1], HDRBYTES);
                tr.duse = modes[i]; // mark component with mode
                memcpy(tr.data, trout[i], nt * FSIZE);
                puttr(&tr);
            }
            if (verbose && (ntr/nc - 1)%(1000/verbose) == 0) { //output info
                gethval(&tr, index, &val);
                warn(" %d-th trace:  px=%.6f py=%.6f of %d compomnents are processed",
                        ntr/nc, px, py, nc);
                // print out coefficients for inspection
                for (i=0; i<nm; ++i) {
                    fprintf(stderr, "   Mode=%2d  P=(%8.5f , %8.5f) Z=(%8.5f , %8.5f) X=(%8.5f , %8.5f) Y=(%8.5f , %8.5f) \n",
                        modes[i], coeff[modes[i]-1][0].r, coeff[modes[i]-1][0].i, coeff[modes[i]-1][1].r, coeff[modes[i]-1][1].i,
                                  coeff[modes[i]-1][2].r, coeff[modes[i]-1][2].i, coeff[modes[i]-1][3].r, coeff[modes[i]-1][3].i);
                }
            }

            memset(*coeff, 0, 2*4*MaxModes*FSIZE);
            memset(*trdata, 0, 4*nt*FSIZE);
            memset(*trout, 0, 2*nm*nf*FSIZE);
        }

        //read next trace
        nsegy = gettr(&tr);
    }

    if (verbose) warn(" Totally %d traces for each of %d compomnents are processed", ntr / nc, nc);

    return (CWP_Exit());
}

float InferV(float Vw, float Dw, float Vp0, float D0, float Vs0, float Rpp, float* Vp, float* Vs) {
    // assuming percentage change of density is n1 times than that (a) of Vp , change (b) of Vs is n2 times than Vp
    // upon percentage change c in reflection coefficient, there is formula
    // a = 2 * c * R0/ ( ( 1 + R0) * ( 1 - R0 ) * ( 1 + n1 ) )  ~= 2cR0/(1+n1);    b = n2 * a;

    float n1 = 3.0, n2 = 5.0; // assumption
    float a, b, c, R0;

    R0 = (Vp0 * D0 - Vw * Dw) / (Vp0 * D0 + Vw * Dw);
    c = Rpp / R0 - 1.0;
    a = 2 * c * R0 / ((1 + R0) * (1 - R0) * (1 + n1));
    b = n2 * a;
    *Vp = (1.0 + a) * Vp0;
    *Vs = (1.0 + b) * Vs0;
    return Dw * Vw * (1.0 + Rpp) / (1.0 - Rpp); // impedance of seafloor
}

