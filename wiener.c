/**-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
C
C
C               WIEN: A SUBROUTINE TO SOLVE WIENER'S EQUATIONS USING THE
C               LEVINSON RECURSIVE ALGORITHM. TAKEN FROM 'DIGITAL
C               FILTERING' BY A.MESKO.
C
C               IT PRODUCES A FILTER (C) WHICH WHEN CONVOLVED WITH
C               THE INPUT SIGNAL X WILL PRODUCE THE DESIRED
C               OUTPUT Z.
C
C               PARAMETERS:
C
C               R(NXC) REAL VECTOR CONTAINING THE AUTOCORRELATION OF
C               THE SIGNAL X .IT SHOULD BE PADDED WITH ZEROS BEYOND
C               THE LENGTH OF X.
C               THE VALUES SHOULD FIRST BE DIVIDED BY THE AUTO-
C               CORRELATION OF THE Z SIGNAL EVALUATED AT TIME ZERO.
C
C               G(NXC) : THE CROSS CORRELATION OF THE INPUT SIGNAL X
C               ,WITH THE OUTPUT Z,DIVIDED BY THE AUTOCORRELATION OF
C               Z EVALUATED AT ZERO TIME.
C
C               NXC: INTEGER, THE LENGTH OF THE CROSS CORRELATION
C               VECTOR.
C
C               C(NXC) :REAL,THE OUTPUT FILTER COEFFICIENTS.
C
 ************************************************************************/
int WIEN(int NXC, float* R, float* G, float* C)
{
    int I, M, MP1;
    float TEMP[NXC], A[NXC];
    //float K, L;


    float RMAX = 0.72E33;

    // Set recursion start values.C
    A[0] = 1.0;
    C[0] = G[0] / R[0];
    float ALPHA = R[0];
    float BETA = R[1];
    float GAMMA = C[0] * R[1];
    if (C[0] > RMAX) return -1;


    //   Main loop
    for (MP1 = 1; MP1 < NXC - 1; ++MP1) {
        float K = -BETA / ALPHA;
        if (K > RMAX) return -1;
        M = MP1 - 1;

        // Part 2 (from mesko)
        for (I = 1; I < M; ++I) TEMP[I] = A[I] + K * A[M - I + 2];
        for (I = 1; I < M; ++I) A[I] = TEMP[I];
        A[MP1] = K * A[0];

        //  Part 3C
        ALPHA = ALPHA + K * BETA;
        BETA = 0.0;
        for (I = 1; I < MP1; ++I) BETA = BETA + A[I] * R[M + 3 - I];

        //  Part 4
        float L = (G[MP1] - GAMMA) / ALPHA;

        //  Part 5
        for (I = 1; I < M; ++I) C[I] = C[I] + L * A[M - I + 2];
        C[MP1] = L * A[1];

        // Part 6
        GAMMA = 0.0;
        for (I = 1; I < MP1; ++I) GAMMA = GAMMA + C[I] * R[M - I + 3];
    }

    return 0;
}

