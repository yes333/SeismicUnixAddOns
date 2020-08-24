#include <math.h>
const int MAXIT=100; //Maximum allowed number of iterations.
const float INVALID_VALUE = -999.0;

void CalcAzim(float x, float* a, float* f, float* df)
{
    float rad = x * PI / 180.0;
    float cosa = cosf(rad);
    float sina = sinf(rad);

    *f = a[0]*cosa*cosa + a[1]*sina*sina + a[2]*cosa*sina + a[3]*cosa + a[4]*sina + a[5];

    *df = a[2]*(cosa*cosa - sina*sina) + 2.0*(a[1] - a[0])*cosa*sina - a[3]*sina + a[4]*cosa;

    return;
}

float rtsafe(void (*funcp)(float, float*, float *, float *), float* a, float x1, float x2,
        float xacc)
//Using a combination of Newton-Raphson and bisection, find the root of a function bracketed
//between x1 and x2. The root, returned as the function value rtsafe, will be refined until
//its accuracy is known within ±xacc. funcp is a user-supplied routine that returns both the
//function value and the first derivative of the function.
{
    int j;
    float df, dx, dxold, f, fh, fl;
    float temp, xh, xl, rts;

    (*funcp)(x1, a, &fl, &df);
    (*funcp)(x2, a, &fh, &df);

    if ( fl*fh > 0.0 ) {
        warn("Opposite values over zero crossing expected: fl(%4.2f)=%7.5f fh(%4.2f)=%7.5f", x1, fl, x2, fh);
        return INVALID_VALUE;
    }

    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    if (fl < 0.0) { //Orient the search so that f(xl) < 0.
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
    }
    rts = 0.5 * (x1 + x2); //Initialize the guess for root,
    dxold = fabs(x2 - x1); //the “stepsize before last,”
    dx = dxold; //and the last step.
    (*funcp)(rts, a, &f, &df);
    for (j = 1; j <= MAXIT; j++) { //Loop over allowed iterations.
        if ((((rts - xh) * df - f)*((rts - xl) * df - f) > 0.0) //Bisect if Newton out of range,
                || (fabs(2.0 * f) > fabs(dxold * df))) { //or not decreasing fast enough.
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rts = xl + dx;
            if (xl == rts) return rts; //Change in root is negligible.
        } else { //Newton step acceptable. Take it.
            dxold = dx;
            dx = f / df;
            temp = rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (fabs(dx) < xacc) return rts; //Convergence criterion.
        (*funcp)(rts, a, &f, &df);
        //The one new function evaluation per iteration.
        if (f < 0.0) //Maintain the bracket on the root.
            xl = rts;
        else
            xh = rts;
    }
    err("Maximum number of iterations exceeded in rtsafe");
    return 0.0; //Never get here.
}
