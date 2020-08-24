// collection of frequently used help functions

const double Invalid_Value = -333.333;

int a2i(const char* strbuf, const char* numset, const int npos)
{
    int i, j, k = 0;
    char *buf = malloc(npos+1);
    memset(buf, 0, npos + 1);
    int l = strlen(numset);

    for(i=0; i < npos; ++i) {
        for (j=0; j < l; ++j) {
            if (strbuf[i] == numset[j]) {
                buf[k++] = strbuf[i];
                break;
            }
        }
    }

    return k > 0 ? atoi(buf) : 0;
}

void setCompName(const int nchn, char* cmpnm)
{
    switch (nchn) {
        case 1:
            strcpy(cmpnm, "P Hydro");
            break;
        case 2:
            strcpy(cmpnm, "Z UD");
            break;
        case 3:
            strcpy(cmpnm, "X EW");
            break;
        case 4:
            strcpy(cmpnm, "Y NS");
            break;
        default:
            strcpy(cmpnm, "UNKNOWN");
            break;
    }
    return;
}

char* getSubstr(const char* strbuf, const char* delim, const char* name, const int pos, int* nc)
{
    int i;
    const char* numset = "-.0123456789";

    char* cs = strstr(strbuf, name); // find name
    if (cs == NULL) {
        *nc = 0;
        return NULL;
    }

    for(i=0; i<pos; ++i) {
        cs = strpbrk(cs, delim);  // skip to end of name
        cs = strpbrk(cs, numset); // start of numbers
        *nc = strcspn(cs, delim);  // length of number string
    }

    return cs;
}

double getValue(const char* strbuf, const char* name, const int pos, const float scale)
{
    int nc;
    const char* delim = " ,/|#\t\n";

    char* cs = getSubstr(strbuf, delim, name, pos, &nc);
    
    if (!cs) return Invalid_Value;

    double d = eatod(cs);
    return d*scale;
}

static void setval(cwp_String type, Value *valp, double dval)
{
    switch (*type) {
        case 's':
            err("can't set char header word");
            break;
        case 'h':
            valp->h = (short) dval;
            break;
        case 'u':
            valp->u = (unsigned short) dval;
            break;
        case 'l':
            valp->l = (long) dval;
            break;
        case 'v':
            valp->v = (unsigned long) dval;
            break;
        case 'i':
            valp->i = (int) dval;
            break;
        case 'p':
            valp->p = (unsigned int) dval;
            break;
        case 'f':
            valp->f = (float) dval;
            break;
        case 'd':
            valp->d = dval;
        default:
            err("unknown type %s", type);
            break;
    }
    return;
}

void Rotate2C(int itmin, int itmax, int skip, float angle, float* x, float* y, float* x1, float* y1)
// x1/y1  rotate  anticlockwise to x/y  | or x/y project to x1/y1 
{
    int it;  // loop counter
    float sina, cosa;

    sina = sin(angle*PI/180.0);
    cosa = cos(angle*PI/180.0);

    for (it=itmin; it<=itmax; ++it) { // loop over samples
        if (skip && (x[it] == 0.0 && y[it] == 0.0)) {
            x1[it - itmin] = 0.0;
            y1[it - itmin] = 0.0;
            continue;
        }
        float a =  cosa*x[it] + sina*y[it];
        float b = -sina*x[it] + cosa*y[it];
        x1[it - itmin] = a;
        y1[it - itmin] = b;
    }
    return;
}


/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

static float quick_select(float *arr, int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    if (n == 1) return arr[0];
    if (n == 2) return 0.5*(arr[0] + arr[1]);

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return (n%2)? arr[median] : 0.5*(arr[median] + arr[median + 1]);

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return (n%2)? arr[median] : 0.5*(arr[median] + arr[median + 1]);
        }

		/* Find median of low, middle and high items; swap into position low */
		middle = (low + high) / 2;
		if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
		if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
		if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

		/* Swap low item (now in position middle) into position (low+1) */
		ELEM_SWAP(arr[middle], arr[low+1]) ;

		/* Nibble from each end towards middle, swapping items when stuck */
		ll = low + 1;
		hh = high;
		for (;;) {
		    do ll++; while (arr[low] > arr[ll]) ;
		    do hh--; while (arr[hh]  > arr[low]) ;

		    if (hh < ll)
		    break;

		    ELEM_SWAP(arr[ll], arr[hh]) ;
		}

		/* Swap middle item (in position low) back into correct position */
		ELEM_SWAP(arr[low], arr[hh]) ;

		/* Re-set active partition */
		if (hh <= median) low = ll;
		if (hh >= median) high = hh - 1;
	}

    return (n%2)? arr[median] : 0.5*(arr[median] + arr[median + 1]);
}

#undef ELEM_SWAP


