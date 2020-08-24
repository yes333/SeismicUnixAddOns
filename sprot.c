/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPCOMPROT v. 1.0 written by Sanyu Ye */


#include <cwp.h>
#include <su.h>
#include <header.h>
#include <segy.h>
#include <segyhdr.h>

/*********************** self documentation ******************************/
char *sdoc[] = {
"									",
" SPCOMPROT - Compute and apply component rotation along given axis     ",
"             for either 4C OBC or 3C VSP data                          ",
"									",
"   spcomprot <stdin >stdout [options]					",
"									",
" Optional Parameters:							",
"                                                                	",
"   nc=4          number of OBC input components, in order H, C3, C2 and C1 ",
"                 usually Hydrophone(H), Vertical(Z), inline(X) and crossline(Y)"
"                 =3  for VSP in order Z, H2, H1                        ",
"                 C1, C2, C3 or H1, H2, Z should build a right hand coordinate system",
"                                                                       ",
"   axis=-1       rotation axis, count backward from last component     ",
"                                                                       ",                                                                       ",
"   a=wevel       keyword containing rotation angle                     ",
"   scale=0.01    scaling factor to convert angle in degrees (-180 ~ 180)",
"                                                                	",
"   verbose=0      =1 output info                                       ",
"                                                                	",
" Note:                                                                 ",
"   hydrophone(H) component will be output without any change           ",
"   output data are sorted trace by trace in order                      ",
"   H, Vertical(Z), inline(X), crossline(Y)                             ",
NULL};

/* Credits:
 *	Sanyu Ye: READ Well Services, April. 2010
 */
/**************** end self doc *******************************************/

// forward declaration
void Rotate(int nt, float a, float b, float** trdata);


int main(int argc, char **argv) {
    cwp_String  key_a;          /* header key word from segy.h		*/
    cwp_String type_a;          /* type of key				*/
    Value val;                  /* value of key			*/
    int index_a;                /* index of key				*/
    int nt;                     /* numsamps as int			*/
    int nc;                     /* number of components			*/
    int nsegy;                  /* length bytes read for segy trace	*/
    segyhdr *hdrs1d = NULL;     /* buffer array of input trace headers  */
    float **trdata=NULL;	/* in/out components   */
    float a;                    /* rotation angle            applied   */
    float scale;                /* input trace px py header values */
    int i, ntr;                 /* loop and trace counter */
    int verbose;		/* =0 no info, >0 output info  */
    segy tr;

    
    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    /* get keys */
    if (!getparstring("a", &key_p))	 key_p="wevel";
    type_p = hdtype(key_p);
    index_p = getindex(key_p);
    if (!getparstring("key", &key))	 key="cdpt";
    type = hdtype(key);
    index = getindex(key);
    
    if (!getparint("nc", &nc)) nc=4;
    if (!getparint("axis", &axis)) axis=-1;
    else {
        if (axis >=0 || axis < -3 || axis < -nc)
            err("invalid axis (=%d), should be between -1 and -3, count backwards", axis);
    }
    if (!getparfloat("scale", &scale)) scale=0.01;
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
            a = vtod(type_p, val)*scale;
            gethval(&tr, index_s, &val);
            b = vtod(type_s, val)*scale;
                        
            Rotate(nt, a, b, trdata);
                       
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
void Rotate(int nt, float a, float b, float** trdata)
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
