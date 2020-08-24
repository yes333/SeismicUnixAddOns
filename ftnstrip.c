/* Copyright (c) Colorado School of Mines, 2007.*/
/* All rights reserved.                       */

/* FTNSTRIP: $Revision: 1.5 $ ; $Date: 2003/08/19 21:24:44 $		*/

#include "par.h"

/*********************** self documentation ******************************/
char *sdoc[] = {
" 									",
" FTNSTRIP - convert a file of binary data plus record delimiters created",
"      via Fortran to a file containing only binary values (as created via C)",
" 									",
" ftnstrip <ftn_data >c_data [optional parameters]			",
"									",
" Optional Parameters:							",
" swap=0  =1 swap bytes of record delimiters                           ",
" 									",
" Caveat: this code assumes the conventional Fortran format of header	",
"         and trailer integer containing the number of byte in the	",
"         record.  This is overwhelmingly common, but not universal.	",
" 									",
NULL};

/* Credits:
 *	RWS: Sanyu Ye  (READ Well Service, Oslo, Norway)
 *           added swap flag as optional parameter
 *	CWP: Jack K. Cohen
 */
/**************** end self doc *******************************************/


int
main(int argc, char **argv)
{
        int swap;
	int n1bytes;
	char *buf;


	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);
        
        if (!getparint("swap",&swap)) swap = 0;

	while (efread(&n1bytes, ISIZE, 1, stdin)) {
                if ( swap ) swap_int_4(&n1bytes);
		buf = ealloc1(n1bytes, 1);
		efread(buf, n1bytes, 1, stdin);
		efwrite(buf, n1bytes, 1, stdout);
		free1(buf);
		efread(&n1bytes, ISIZE, 1, stdin);
	}

	return(CWP_Exit());
}
