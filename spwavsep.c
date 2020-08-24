/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPWAVSEP v. 1.0 written by Sanyu Ye */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"									",
" SPWAVSEP - Wave Seperation for VSP data in tau-p domain               ",
"									",
"    spwavsep <stdin >stdout [options]					",
"									",
" Optional Parameters:							",
"                                                                	",
"    vp=wevel       keyword containing P-velocity at receiver [m/s]	",
"    vs=swevel      keyword containing S-velocity at receiver [m/s]	",
"                                                                	",
"    verbose=0      =1 output info                                      ",
"                                                                	",
" Note:                                                                 ",
"    ray parameter p must be stored on tr.fx in [s/km] as output by SPFK2TAUP",
"    input data must be sorted trace by trace in order Z, X (H radial)	",
"    output data are sorted trace by trace in order P, S           	",
NULL};

/* Credits:
 *	Sanyu Ye: READ Well Services, April. 2010
 */
/**************** end self doc *******************************************/

// forward declaration
void Decompose(int nt, float Vp, float Vs,float p, float** trdata);


int main(int argc, char **argv) {
    cwp_String key_p, key_s;	/* header key word from segy.h		*/
    cwp_String type_p, type_s;  /* type of key				*/
    Value val;                  /* value of key			*/
    int index_p, index_s;       /* index of key				*/
    int nt;                     /* numsamps as int			*/
    int nc = 2;                 /* number of components			*/
    int nsegy;                  /* length bytes read for segy trace	*/
    segyhdr *hdrs1d = NULL;     /* buffer array of input trace headers  */
    float **trdata=NULL;	/* tau-p-domain for all input components   */
    float Vp, Vs;               /* P-, S-velocity at receiver site */
    float p;                    /* input trace px py header values */
    int i, ntr;                 /* loop and trace counter */
    int verbose;		/* =0 no info, >0 output info  */
    segy tr;

    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* get keys */
    if (!getparstring("vp", &key_p))	 key_p="wevel";
    type_p = hdtype(key_p);
    index_p = getindex(key_p);
    if (!getparstring("vs", &key_s))	 key_s="swevel";
    type_s = hdtype(key_s);
    index_s = getindex(key_s);
    
    if (!getparint("verbose", &verbose)) verbose=0;
    
    /* get first trace */
    if ( (nsegy=gettr(&tr)) < HDRBYTES) err("can't get first trace");

    nt = tr.ns;

    hdrs1d = (segyhdr *) ealloc1(nc, HDRBYTES);
    trdata = ealloc2float(nt, nc);
    
    //decomposing wavefield
    ntr = 0;
    /* Main loop over traces */
    while (nsegy > HDRBYTES ) {
        i = ntr % nc;
        memcpy(&hdrs1d[i], &tr, HDRBYTES);  //cache header
        memcpy(trdata[i], tr.data, nt*FSIZE); // cache data

        if ( ++ntr%nc == 0 ) { // nc traces read in, start separation
            gethval(&tr, index_p, &val);
            Vp = vtod(type_p, val)/1000.0;
            gethval(&tr, index_s, &val);
            Vs = vtod(type_s, val)/1000.0;
            
            p = tr.fx;
            
            Decompose(nt, Vp, Vs, p, trdata); 
                       
            // output
            for (i=0; i<nc; ++i) {
                memcpy(&tr, &hdrs1d[i], HDRBYTES);
                memcpy(tr.data, trdata[i], nt*FSIZE);
                puttr(&tr);
            }
        }

        //read next trace
        nsegy = gettr(&tr);
    }

    if (verbose) warn(" Totally %d traces for each of %d compomnents are processed", ntr/nc, nc);
    
    if (hdrs1d) free(hdrs1d);
    if (trdata) free2float(trdata);
    
    return(CWP_Exit());
}

// back project X and Z to P and S according their incident angles
void Decompose(int nt, float Vp, float Vs,float p, float** trdata)
{
    int i;  // loop counter
    float ax, az, sinp, cosp, sins, coss, f;
    //float const maxsintheta = 0.9798; // correspond max angle of 78 degree
    float const maxsintheta = 1.0;

    sinp = p*Vp;  // sin(THETA-P)
    if ( fabs(sinp) >= maxsintheta ) {
        sinp = (sinp > 0.0) ? maxsintheta : -maxsintheta;
        cosp = 0.0;
    } else {
        cosp = sqrt(1.0 - sinp*sinp);  // cos(THETA-P)  <= 0.2
    }
    sins = p*Vs;  // sin(THETA-S)
    if ( fabs(sins) >= maxsintheta ) {
        sins = (sins > 0.0) ? maxsintheta : -maxsintheta;
        coss = 0.0;
    } else {
        coss = sqrt(1.0 - sins*sins);  // cos(THETA-S) <= 0.2
    }
    f = sinp*sins + cosp*coss;

    for (i=0; i<nt; ++i) { // loop over samples
        ax = trdata[1][i];
        az = trdata[0][i];
        trdata[0][i] = (ax*sins + az*coss)/f;
        trdata[1][i] = (ax*cosp - az*sinp)/f;
    }
}
