/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* SUCHW: $Revision: 1.16 $ ; $Date: 2007/04/04 19:31:34 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"

/*********************** self documentation **********************/
char *sdoc[] = {
"									",
" SUCHW - Change Header Word using one, two or three header word fields	",
"									",
"  suchw <stdin >stdout [optional parameters]				",
"									",
" Required parameters:							",
" none									",
"									",
" Optional parameters:							",
" key1=cdp,...	output key(s) 						",
" key2=cdp,...	input key(s) 						",
" key3=cdp,...	input key(s)  						",
" key4=....	optional input key(s), ignored if not given             ",
" a=0,...		overall shift(s)				",
" b=1,...		scale(s) on first input key(s) 			",
" c=0,...		scale on second input key(s) 			",
" d=1,...		overall scale(s)				",
" e=1,...		exponent on first input key(s)                  ",
" f=1,...		exponent on second input key(s)                 ",
" g=1,...		exponent on third input key(s)                  ",
" abs=0,...		=1,... use absolute value of key2(s)            ",
"									",
" Note: if the number of a constant array is less than key1�s,          ",
"       the last one is used for the rest of the array			",
"									",
" The value of header word key1 is computed from the values of		",
" key2, key3 and optional key4 by:					",
"									",
"  val(key1) = (a + b * ([abs]val(key2))^e + c * val(key3)^f)[*val(key4)^g] / d",
"									",
" Examples:								",
" Shift cdp numbers by -1:						",
"	suchw <data >outdata a=-1					",
"									",
" Add 1000 to tracr value:						",
" 	suchw key1=tracr key2=tracr a=1000 <infile >outfile		",
"									",
" We set the receiver point (gx) field by summing the offset		",
" and shot point (sx) fields and then we set the cdp field by		",
" averaging the sx and gx fields (we choose to use the actual		",
" locations for the cdp fields instead of the conventional		",
" 1, 2, 3, ... enumeration):						",
"									",
"   suchw <indata key1=gx key2=offset key3=sx b=1 c=1 |			",
"   suchw key1=cdp key2=gx key3=sx b=1 c=1 d=2 >outdata			",
"									",
" Do both operations in one call:					",
"									",
" suchw<indata key1=gx,cdp key2=offset,gx key3=sx,sx b=1,1 c=1,1 d=1,2 >outdata",
"									",
" Divide water depth at receiver by scaling factor and print out        ",
"									",
" suchw <indata key1=scalel key2=scalel abs=1 |				",
" suchw key1=unscale key2=gwdep key4=scalel g=-1 |			",
" sugethw output=geom key=cdp,cdpt,unscale				",
"									",
NULL};

/* Credits:
 *	SEP: Einar Kjartansson
 *	CWP: Jack K. Cohen
 *      CWP: John Stockwell, 7 July 1995, added array of keys feature
 *      Delphi: Alexander Koek, 6 November 1995, changed calculation so
 *              headers of different types can be expressed in each other
 *      RWS: Sanyu Ye, Oct. 2007, READ Well Service, Oslo,
 *              added absolute value for key2 and  
 *              optional key4 array to enable multiplication/division
 */
/**************** end self doc ***********************************/

/* prototype for function used internally */
void changeval(cwp_String type1, Value *valp1, cwp_String type2,
	       Value *valp2, cwp_String type3, Value *valp3,
		double a, double b, double c, double d, double e, double f, double g,
                int ab, int nkey4, cwp_String type4, Value *valp4);

segy tr;

int
main(int argc, char **argv)
{
	cwp_String key1[SU_NKEYS];	/* output key(s)		*/
	cwp_String key2[SU_NKEYS];	/* first input key(s)		*/
	cwp_String key3[SU_NKEYS];	/* second input key(s)		*/
	cwp_String key4[SU_NKEYS];	/* third input key(s)		*/
	cwp_String type1[SU_NKEYS];	/* array of types for key1	*/
	cwp_String type2[SU_NKEYS];	/* array of types for key2	*/
	cwp_String type3[SU_NKEYS];	/* array of types for key3	*/
	cwp_String type4[SU_NKEYS];	/* array of types for key4	*/
	int nkeys;			/* number of keys to be computed*/
	int nkey4;			/* number of key4s =0 ignored   */
	int n;				/* counter of keys getparred	*/
	int ikey;			/* loop counter of keys 	*/
	int index1[SU_NKEYS];		/* array of indexes for key1 	*/
	int index2[SU_NKEYS];		/*      ....        for key2	*/
	int index3[SU_NKEYS];		/*      ....        for key3	*/
	int index4[SU_NKEYS];		/*      ....        for key4	*/

	Value val1;			/* value of key1		*/
	Value val2;			/* value of key2		*/
	Value val3;			/* value of key3		*/
	Value val4;			/* value of key4		*/

	double *a=NULL;			/* array of "a" values		*/
	double *b=NULL;			/* array of "b" values		*/
	double *c=NULL;			/* array of "c" values		*/
	double *d=NULL;			/* array of "d" values		*/
	double *e=NULL;			/* array of "e" values		*/
	double *f=NULL;			/* array of "f" values		*/
	double *g=NULL;			/* array of "g" values		*/
        int *ab=NULL;                   /* array of "abs" values	*/
        
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);

	/* Get parameters */
	/* get key1's */
	if ((n=countparval("key1"))!=0){
		nkeys=n;
		getparstringarray("key1",key1);
	} else { /* set default */
		nkeys=1;
		key1[0]="cdp";	
	}

	/* get key2's */
	if ((n=countparval("key2"))!=0){
		if (n!=nkeys)
			err("number of key2's and key1's must be equal!");

		getparstringarray("key2",key2);
	} else { /* set default */
		if (nkeys!=1)
			err("number of key2's and key1's must be equal!");

		key2[0]="cdp";	
	}

	/* get key3's */
	if ((n=countparval("key3"))!=0){
		if (n!=nkeys)
			err("number of key3's and key1's must be equal!");

		getparstringarray("key3",key3);
	} else { /* set default */
            for (ikey=0; ikey < nkeys; ++ikey) key3[ikey]="cdp";
	}
        
	/* get key4's */
	if ((nkey4=countparval("key4"))!=0){
		if (nkey4!=nkeys)
			err("number of key4's and key1's must be equal!");

		getparstringarray("key4",key4);
	} 
        
	/* get a's */
	if ((n=countparval("a"))!=0) {
            if (n>nkeys) {
		err("number of a's cannot be large than key1's!");
            } else {               
        	a=ealloc1double(nkeys);
                getpardouble("a",a);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        a[ikey]=a[n-1];
                    }
                }
            }
	} else { /* set default */
        	a=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			a[ikey]=0.;
	}
        
	/* get b's */
	if ((n=countparval("b"))!=0) {
            if (n>nkeys) {
		err("number of b's cannot be large than key1's!");
            } else {               
        	b=ealloc1double(nkeys);
                getpardouble("b",b);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        b[ikey]=b[n-1];
                    }
                }
            }
	} else { /* set default */
        	b=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			b[ikey]=1.;
	}
        
	/* get c's */
	if ((n=countparval("c"))!=0) {
            if (n>nkeys) {
		err("number of c's cannot be large than key1's!");
            } else {               
        	c=ealloc1double(nkeys);
                getpardouble("c",c);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        c[ikey]=c[n-1];
                    }
                }
            }
	} else { /* set default */
        	c=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			c[ikey]=0.;
	}
        
	/* get d's */
	if ((n=countparval("d"))!=0) {
            if (n>nkeys) {
		err("number of d's cannot be large than key1's!");
            } else {               
        	d=ealloc1double(nkeys);
                getpardouble("d",d);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        d[ikey]=d[n-1];
                    }
                }
            }
	} else { /* set default */
        	d=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			d[ikey]=1.;
	}
        
	/* get e's */
	if ((n=countparval("e"))!=0) {
            if (n>nkeys) {
		err("number of e's cannot be large than key1's!");
            } else {               
        	e=ealloc1double(nkeys);
                getpardouble("e",e);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        e[ikey]=e[n-1];
                    }
                }
            }
	} else { /* set default */
        	e=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			e[ikey]=1.;
	}
        
	/* get f's */
	if ((n=countparval("f"))!=0) {
            if (n>nkeys) {
		err("number of f's cannot be large than key1's!");
            } else {               
        	f=ealloc1double(nkeys);
                getpardouble("f",f);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        f[ikey]=f[n-1];
                    }
                }
            }
	} else { /* set default */
        	f=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			f[ikey]=1.;
	}
        
	/* get g's */
	if ((n=countparval("g"))!=0) {
            if (n>nkeys) {
		err("number of g's cannot be large than key1's!");
            } else {               
        	g=ealloc1double(nkeys);
                getpardouble("g",g);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        g[ikey]=g[n-1];
                    }
                }
            }
	} else { /* set default */
        	g=ealloc1double(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			g[ikey]=1.;
	}
        
	/* get abs's */
	if ((n=countparval("abs"))!=0) {
            if (n>nkeys) {
		err("number of abs's cannot be large than key1's!");
            } else {               
        	ab=ealloc1int(nkeys);
                getparint("abs",ab);
                if (n<nkeys) {
                    for (ikey=n; ikey<nkeys; ++ikey) {
                        ab[ikey]=ab[n-1];
                    }
                }
            }
	} else { /* set default */
        	ab=ealloc1int(nkeys);
		for (ikey=0; ikey<nkeys; ++ikey)
			ab[ikey]=0;
	}
        
	for (ikey=0; ikey<nkeys; ++ikey) {

            /* get types and index values */
            type1[ikey]  = hdtype(key1[ikey]);
            type2[ikey]  = hdtype(key2[ikey]);
            type3[ikey]  = hdtype(key3[ikey]);
            index1[ikey] = getindex(key1[ikey]);
            index2[ikey] = getindex(key2[ikey]);
            index3[ikey] = getindex(key3[ikey]);
            if (nkey4>0) {
                type4[ikey]  = hdtype(key4[ikey]);
                index4[ikey] = getindex(key4[ikey]);                
            }
	}

	/* loop over traces */
	while (gettr(&tr)) {

		/* loop over key fields */
		for (ikey=0; ikey<nkeys; ++ikey) {
			
			/* get header values */
			gethval(&tr, index2[ikey], &val2);
			gethval(&tr, index3[ikey], &val3);
                        if (nkey4>0) {
                            gethval(&tr, index4[ikey], &val4);
                            type4[ikey]  = hdtype(key4[ikey]);
                            index4[ikey] = getindex(key4[ikey]);                
                        }

			changeval(type1[ikey], &val1, type2[ikey], &val2,
				type3[ikey], &val3, a[ikey], b[ikey], c[ikey],
				d[ikey], e[ikey], f[ikey], g[ikey], ab[ikey], 
                                nkey4, type4[ikey], &val4);
			puthval(&tr, index1[ikey], &val1);
		}
		puttr(&tr);
	}
        
        if (a) free1double(a);
        if (b) free1double(b);
        if (c) free1double(c);
        if (d) free1double(d);
        if (e) free1double(e);
        if (f) free1double(f);
        if (g) free1double(g);
        if (ab) free1int(ab);

	return(CWP_Exit());
}


void changeval(cwp_String type1, Value *valp1, cwp_String type2,
	       Value *valp2, cwp_String type3, Value *valp3,
		double a, double b, double c, double d, double e, double f, double g,
                int ab, int nkey4, cwp_String type4, Value *valp4)
{
        double dval1=0.0, dval2=0.0, dval3=0.0, dval4=0.0;
	dval2=vtod( type2, *valp2);
	if (ab) dval2=fabs(dval2);
	dval3=vtod( type3, *valp3);
	dval1=(a+b*pow(dval2,e)+c*pow(dval3,f))/d;
        if (nkey4>0) {
            dval4=vtod( type4, *valp4);
            if (dval4 == 0 && g < 0.0) {
                warn("val(key1) = 0 set due to division by zero : val(key4)=0");
                dval1 = 0.0;
            } else {
                dval1 *= pow(dval4,g);
            }
        }
        
	switch (*type1) {
	case 's':
		err("can't change char header word");
	break;
	case 'h':
		valp1->h = (short) dval1;
	break;
	case 'u':
		valp1->u = (unsigned short) dval1;
	break;
	case 'l':
		valp1->l = (long) dval1;
	break;
	case 'v':
		valp1->v = (unsigned long) dval1;
	break;
	case 'i':
		valp1->i = (int) dval1;
	break;
	case 'p':
		valp1->p = (unsigned int) dval1;
	break;
	case 'f':
		valp1->f = (float) dval1;
	break;
	case 'd':
		valp1->d = (double) dval1;
	break;
	default:
		err("unknown type %s", type1);
	break;
	}
}
