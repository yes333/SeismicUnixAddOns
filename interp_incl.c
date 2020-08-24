
/* Copyright (c) READ Well Service, 2007.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/******************************************************************************
* Functions with prototypes for functions used internally to interpolate
* header values based on up to 3 keywrods
* Written by Sanyu Ye, sanyu.ye@readgroup.com
******************************************************************************
*
******************************************************************************/
/**************** end self doc ********************************/

int FindInterval(
    int ikey,        /* key level/index to search (which match key) */
    VALUETYPE dv,       /* value to be searched */
    VALUETYPE *dprev,   /* values of previous keys that must be matched during the search */
    VALUETYPE **keys,  /* array containing the input table of match keys */
    int *sort,      /* sorting order of the input match keys, >0 ascending, <0 descending */
    int nRows,       /* Total number of the rows/lines of input match keys */
    int *nbegp,      /* input: starting search line; output: new starting line for next search */
    int* n1p, int* n2p, /* interval of line number between which the searched value located */
    VALUETYPE* d1p, VALUETYPE* d2p /* interval of key values between which the search value located */
    /* return value =-2 beyond first row; =+-1 between two rows; =0 exact match; =2 beyond last row */
    ) 
{
    int startIntervalFound = 0, match = sort[ikey];
    VALUETYPE prevKeyValue = (sort[ikey] > 0 )? (VALUETYPE) INT_MIN : (VALUETYPE) INT_MAX;
    for (int n = *nbegp; n < nRows; ++n) {  // loop over rows of match keys
        int isMatchPrevKeys = 1;
        // in case of search for 2nd or third key, find first match in previous keys
        for (int i=0; i<ikey && isMatchPrevKeys; ++i) {
            isMatchPrevKeys = isMatchPrevKeys && dprev[i] == keys[n][i];
        }

        if (!isMatchPrevKeys ) continue; /* skip to first line that previous key values matched*/

        if ( keys[n][ikey] == prevKeyValue ) continue;  // skip line with same key value
        else { // key value change
            prevKeyValue = keys[n][ikey];

            if ( dv == keys[n][ikey] ) { // exact match found, break
                *nbegp = n;  /* new start point for next level search */
                *d1p = *d2p = keys[n][ikey];
                *n1p = *n2p = n;
                match = 0;
                break;
            } else if ( (sort[ikey] > 0 && dv > keys[n][ikey]) 
                     || (sort[ikey] < 0 && dv < keys[n][ikey]) ) { // first value of interval found
                *nbegp = n;  /* new start point for next level search */
                *d1p = *d2p = keys[n][ikey];
                *n1p = *n2p = n;
                startIntervalFound = 1;
            } else if ( (sort[ikey] > 0 && dv < keys[n][ikey]) /*check if second value of interval found */ 
                     || (sort[ikey] < 0 && dv > keys[n][ikey]) ) {
                *d2p = keys[n][ikey];
                *n2p = n;
                
                if (!startIntervalFound) { // even the start is not found
                    *nbegp = n;  /* new start point for next level search */
                    *d1p = *d2p;
                    *n1p = *n2p;
                    match = 2;
                }
                break;
            }
        } 
    }
    return match;
}

void interp1(
    int     ikey,   /* index of match key */
    int*   sort,    /* sorting array */
    int     n1,     /* first row */
    int     n2,     /* last row between them the values are interpolated */
    VALUETYPE  d,      /* key value of trace to which the interpolation is done */
    VALUETYPE  **keys, /* 2-D array of input table values of match keys */
    VALUETYPE  **var, /* 2-D array of input table values */
    int     nv,     /* array size */
    VALUETYPE* v   /* output values interpolated */
    ) 
{
    int i;
    for ( i=0; i<nv; ++i) {
        if ( (sort[ikey] > 0 && d <= keys[n1][ikey]) || (sort[ikey] < 0 && d >= keys[n1][ikey]) )  v[i] = var[n1][i];
        else if ( (sort[ikey] > 0 && d >= keys[n2][ikey]) || (sort[ikey] < 0 && d <= keys[n2][ikey]) ) v[i] = var[n2][i];
        else if( n1 == n2 || keys[n2][ikey] == keys[n1][ikey] ) v[i] = (var[n1][i] + var[n2][i])/2;
        else v[i] = var[n1][i] + (var[n2][i] - var[n1][i])*(d - keys[n1][ikey])/(keys[n2][ikey] - keys[n1][ikey]);
    }
    return;
}

void interp2(
    int         sort,  /* sort direction */
    VALUETYPE  d,      /* key value of trace to which the interpolation is done */
    VALUETYPE  d1,      /* first key value  */
    VALUETYPE  d2,      /* second key value  */
    VALUETYPE* v1,  /* first 1-D array to interpolation between */
    VALUETYPE* v2,  /* second 1-D array to interpolation between */
    int     nv,    /* array size */
    VALUETYPE* v    /* output values interpolated */
    ) 
{
    int i;
    for ( i=0; i<nv; ++i) {
        if ( (sort > 0 && d <= d1) || (sort < 0 && d >= d1) ) v[i] = v1[i];
        else if ( (sort > 0 && d >= d2) || (sort < 0 && d <= d2) ) v[i] = v2[i];
        else if( d1 == d2 ) v[i] = (v1[i] + v2[i])/2;
        else v[i] = v1[i] + (v2[i] - v1[i])*(d - d1)/(d2 - d1);
    }
    return;
}

int interp(
    int nkeys,      /* Number of the match keys */
    VALUETYPE *d,      /* 1D array of key values to be searched */
    VALUETYPE **keys,  /* 2D array containing the input table of match keys */
    int *sort,      /* sorting order of the input match keys, >0 ascending, <0 descending */
    int nRows,      /* Total number of the rows/lines of input match keys */
    VALUETYPE **var,   /* 2D array containing the input table of values */
    int nv,         /* Total number of the keys to be set */
    VALUETYPE* v       /* interval of key values between which the search value located */
    /* return value =-2 beyond first row; =+-1 between two rows; =0 exact match; =2 beyond last row */
    ) 
{
    int n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14;
    int nbeg1=0, nbeg3, nbeg5, nbeg7, nbeg9, nbeg11, nbeg13;
    VALUETYPE d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14;
    VALUETYPE *v3=NULL, *v5=NULL, *v7=NULL, *v9=NULL, *v11=NULL, *v13=NULL;
    VALUETYPE dprev[2] = {0, 0};
    int match;

    if ( nkeys > 3 ) err("interpolation not implemented for more than 1D");

    match = FindInterval(0, d[0], dprev, keys, sort, nRows, &nbeg1, &n1, &n2, &d1, &d2);
    if (verbose > 100) fprintf(stderr, "%d %f %f %f %f %d %d %d %f %f\n", 0, d[0], dprev[0], keys[n1][0], keys[n2][0], nbeg1, n1, n2, d1, d2);
    if ( 1 == nkeys ) { /* just 1-D interpolation */
        interp1(0, sort, n1, n2, d[0], keys, var, nv, v );
        return 1;
    } else { /* one level more */
        dprev[0] = d1; nbeg3 = n1;
        match = FindInterval(1, d[1], dprev, keys, sort, nRows, &nbeg3, &n3, &n4, &d3, &d4);
        if (verbose>100) fprintf(stderr, "%d %f %f %f %f %d %d %d %f %f\n", 1, d[1], dprev[0], keys[n3][1], keys[n4][1], nbeg3, n3, n4, d3, d4);
        if ( d1 != d2 ) {
            dprev[0] = d2; nbeg5 = n2;
            match = FindInterval(1, d[1], dprev, keys, sort, nRows, &nbeg5, &n5, &n6, &d5, &d6);
            if (verbose>100) fprintf(stderr, "%d %f %f %f %f %d %d %d %f %f\n", 1, d[1], dprev[0], keys[n5][1], keys[n6][1], nbeg5, n5, n6, d5, d6);
        } else {
            n5 = n3; n6 = n4; d5 = d3; d6 = d4;
        }
        if ( 2 == nkeys ) { /* 2D */
            v3 = EALLOC1VALUETYPE(nv);
            v5 = EALLOC1VALUETYPE(nv);
            interp1(1, sort, n3, n4, d[1], keys, var, nv, v3 );
            interp1(1, sort, n5, n6, d[1], keys, var, nv, v5 );
            interp2(sort[0], d[0], d1, d2, v3, v5, nv, v );
            FREE1VALUETYPE(v3);
            FREE1VALUETYPE(v5);
            return 2;
        } else { /* 3D, one level more */
            dprev[0] = d1; dprev[1] = d3; nbeg7 = n3;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg7, &n7, &n8, &d7, &d8);
            dprev[0] = d1; dprev[1] = d4; nbeg9 = n4;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg7, &n9, &n10, &d9, &d10);
            dprev[0] = d2; dprev[1] = d5; nbeg11 = n5;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg11, &n11, &n12, &d11, &d12);
            dprev[0] = d2; dprev[1] = d6; nbeg13 = n6;
            match = FindInterval(2, d[2], dprev, keys, sort, nRows, &nbeg13, &n13, &n14, &d13, &d14);
            v3 = EALLOC1VALUETYPE(nv);
            v5 = EALLOC1VALUETYPE(nv);
            v7 = EALLOC1VALUETYPE(nv);
            v9 = EALLOC1VALUETYPE(nv);
            v11 = EALLOC1VALUETYPE(nv);
            v13 = EALLOC1VALUETYPE(nv);
            interp1(2, sort, n7, n8, d[2], keys, var, nv, v7 );
            interp1(2, sort, n9, n10, d[2], keys, var, nv, v9 );
            interp1(2, sort, n11, n12, d[2], keys, var, nv, v11 );
            interp1(2, sort, n13, n14, d[2], keys, var, nv, v13 );
            interp2(sort[1], d[1], d3, d4, v7, v9, nv, v3 );
            interp2(sort[1], d[1], d5, d6, v11, v13, nv, v5 );
            interp2(sort[0], d[0], d1, d2, v3, v5, nv, v );
            FREE1VALUETYPE(v3);
            FREE1VALUETYPE(v5);
            FREE1VALUETYPE(v7);
            FREE1VALUETYPE(v9);
            FREE1VALUETYPE(v11);
            FREE1VALUETYPE(v13);
            return 3;
        }
    }
}
