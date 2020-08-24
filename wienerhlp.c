
int CalcMix(int mix, int nwx, int nwy, float** weight)
{
    int i, j;
    // only mix=1 implemented yet

    for (i=0; i<nwy; ++i) {
        for (j=0; j<nwx; ++j) {
            weight[i][j] = 1.0;
        }
    }
    return nwx*nwy;
}

int CalcAvgAC(int nxmax, int nymax, int** trid, float*** acorr,
        int nwx, int nwy, float** weight, int ix, int iy, int ncorr, float* AvgAC)
{
    int i, j, ixx, iyy, it, nw=0;

    for (it=0; it<ncorr; ++it) {
        float w = 0.0; // total weighting factor
        float z = 0.0; // weighted acorr at sample it
        for (nw=0, i=0; i<nwy; ++i) {
            for (j=0; j<nwx; ++j) {
                iyy = iy - nwy/2 + i;
                ixx = ix - nwx/2 + j;
                if (ixx < 0 || ixx >= nxmax || iyy < 0 || iyy >= nymax) continue;
                if (!trid[iyy][ixx] || acorr[iyy][ixx][it] == 0.0) continue; // skip sample with zero acorr
                ++nw;
                w += weight[i][j];
                z += weight[i][j]*acorr[iyy][ixx][it];
            }
        }
        AvgAC[it] = (w > 1.0 )? z / w : z;
    }
    return nw;
}
