/*
 *
 **********************************************************************************************
 *
 * Compute Inverse global forward slant stack (tau-p transform) via
 * FK  transform
 *
 * Input is:
 *
 * ---------------------------------------------------------------------------------------------*
 * float *traces:   Input tau-p trace data (npx*npy*ntau)                                       *
 * int *ntt:        length of data trace in time domain                                         *
 * int *ntmax:      length of time vector (padded with zeroes from trace length)                *
 * int *nxxmax:     Max number of in-line traces                                                *
 * int *nyymax:     Max number of cross-line traces                                             *
 * int *nxx:        number of in-line traces                                                    *
 * int *nyy:        number of x-line traces                                                     *
 * float *dtt:      time Sampling in seconds                                                    *
 * float *dxx:      in-line spatial trace sampling (meters)                                     *
 * float *dyy:      x-line spatial trace sampling (meters)                                      *
 * float *xxmin:    in-line offset of first trace                                               *
 * float *yymin:    x-line offset of first trace                                                *
 * int *npxx:       number of output px-traces                                                  *
 * int *npyy:       number of output py-traces                                                  *
 * float *pxxmin:   minimum slowness px                                                         *
 * float *pyymin:   minimum slowness py                                                         *
 * float *maxfreqq: maximum temporal frequency                                                  *
 * float *work:     Work array of size -> 4*(*ntmax) + 4*nxfft + 6*nxfft*nyfft (spatial ffts)   *
 *
 * ---------------------------------------------------------------------------------------------*
 *
 * Output is:                                                                                   *
 *
 * float *traces:   inverse tau-p transformed data i.e (x,y,t)- doamain                         *
 *
 * ---------------------------------------------------------------------------------------------*
 *
 *
 */





/******************************************************************************
 *
 * Compute inverse global slant stack in FK domain
 *
 ******************************************************************************/


void inv_taup_3d_(float *traces, int *ntt, int *nttmax, int *nxxmax, int *nyymax, int *nxx, int *nyy, float *dtt, float *dxx,
        float *dyy, float *xxmin, float *yymin, int *npxx, int *npyy, float *pxxmin, float *pyymin,
        // float *maxfreqq, int *nxxfft, int *nyyfft, float *work)
        float *pxxmax, float *pyymax, float *maxfreqq, int *nxxfft, int *nyyfft, float *slice, float *work)
        
{
    int i, ix, ip, itau, iw, ik, it, iy; /* loop counters */
    float maxfreq, fw;       /* first and maximum frequency */
    float dw, dk;            /* frequency and wavenumber sampling interval */
    int nw, nkx, nky;         /* number of frequencies and wavenumbers */
    int ntaufft;            /* number of padded time samples */
    int nkyfft;             /* number of padded wavenumber samples along Y-axis*/
    int nkxfft;             /* number of padded wavenumber samples along X-axis*/
    int nkmax;
    int ntau = *ntt;
    float pxmax, pymax;
    float dtau= *dtt;	/* time intercept sampling ingerval */
    float fpy=0.;           /* first slope in px */
    float fpx=0.;           /* first slope in py */
    float fk;               /* first wavenumber */
    float w, k;		/* angular frequency, slope and wavenumber */
    float fftscl;  		/* scale factor for FFT */
    float *pk;		/* k-interpolated w trace */
    float *tr_fft;          /* padded trace for FFT */
    complex *tr_p;
    float *ctr;          /* */
    complex *ctr_x;        /* */
    complex *hk;            /* slant stacked single trace */
    complex *ole;
    complex czero;		/* complex zero number */
    float fky, fkx, dky, dkx;
    float xshift, c, s, temp, phase;
    int nkymax, nkxmax, ipx, ipy, ikx, iky, nkxx;
    complex *cwork, *help;
    int ntmax, nxmax, nymax, ntw, nnw, nwork;
    float wscal, dteta, dteta2, teta, vpxmax, vpymax, ddpx;
    int it0, it1, it2, itx, ipx2, ipy2, iptx, ipty;
    
    
    float dt, xmin, ymin, dx, dy, pxmin, pymin, dpx, dpy, pxx, pyy;
    int nt, nx, ny, npx, npy;
    
    /* Convert from input pointers */
    
    nxmax = *nxxmax;
    nymax = *nyymax;
    dt = *dtt;
    nt = *ntt;
    ntmax = *nttmax;
    ny = *npyy;
    nx = *npxx;
    ymin = *yymin*0.001;
    xmin = *xxmin*0.001;
    dy = *dyy*0.001;
    dx = *dxx*0.001;
    
    
    npy = *nyy;
    npx = *nxx;
    pymin = *pyymin; // ignored *1000.;
    pxmin = *pxxmin; // ignored *1000.;
    dpy = 2.*fabs(pymin)/(npy-1);
    dpx = 2.*fabs(pxmin)/(npx-1);
    maxfreq = *maxfreqq;
    
    vpxmax = *pxxmax; // ignored*1000;
    vpymax = *pyymax; // ignored*1000;
    if(vpxmax < pxmin)  vpxmax = pxmin;
    if(vpymax < pymin)  vpymax = pymin;
    
    
    /* compute slope sampling interval */
    czero = cmplx(0.0, 0.0);
    pymax=pymin+(npy-1)*dpy;
    pxmax=pxmin+(npx-1)*dpx;
    
    /* determine lengths and scale factors for FFT's */
    fpy=pymin+0.5*(pymax-pymin-(npy-1)*dpy);
    fpx=pxmin+0.5*(pxmax-pxmin-(npx-1)*dpx);
    pymax=fpy+(npy-1)*dpy;
    pxmax=fpx+(npx-1)*dpx;
    ntaufft = npfar(1.5*ntau);
    fftscl = 1.0/ntaufft;
    nkyfft =1;
    if(ny > 1) nkyfft = npfa(npy);
    nkxfft = npfa(npx);
    
    nkxfft = *nxxfft;
    nkyfft = *nyyfft;
    
    ntw = ntmax/2;
    
    /* determine frequency and wavenumber sampling */
    nw = ntaufft/2+1;
    dw = 2.0*PI/(ntaufft*dtau);
    fw = 0.;
    fky = -PI/dy;
    fkx = -PI/dx;
    
    dky = 2.0*PI/(nkyfft*dy);
    dkx = 2.0*PI/(nkxfft*dx);
    nkymax = nkyfft;
    nkxmax = nkxfft;
    nky = nkyfft;
    nkx = nkxfft;
    
    /* allocate working space for FFT's */
    
    // nwork = 3*ntaufft + nkxmax*nkymax + 2*(npx*npy + nkxmax*nkymax + nkxfft*nkyfft);
    
    tr_fft = (float*)work;
    nwork = ntaufft;
    pk = (float*)(work+nwork);
    nwork = nwork + nkxmax;
    
    ctr = (float*)traces;
    ctr_x= (complex*)traces;
    
    cwork = (complex*)(work+nwork);
    
    tr_p = (complex*)(cwork);
    nwork = npx + npy;
    hk = (complex*)(cwork+nwork);
    nwork = nwork + nkxmax + nkymax;
    help = (complex*)(slice);
    
    
    nnw = 2.*PI*maxfreq/dw + 1;
    
    
    // Do tau-p muting defined from input parameters pxmax & pymac
    
    
    ddpx = vpxmax - pxmin;
    itx = nt*0.04;
    iptx = npx*0.1;
    ipty = npy*0.1;
    dteta = 0.5* PI/itx;
    dteta2 = PI/iptx;
    
    if(fabs(ddpx) < dpx) ddpx = dpx;
    for (ipy=0; ipy<npy; ipy++) {
        for (ipx=0; ipx<npx; ipx++) {
            pxx = pxmin + ipx*dpx;
            if(fabs(pxx) < fabs(vpxmax)) continue;
            it0 = nt*((-fabs(pxx) - pxmin)/ddpx);
            it1 = it0-itx;
            if(it1 < 0) it1 = 0;
            it2 = it0+itx;
            if(it2 > nt ) it2 = nt;
            
            for (it=it1; it<it2; it++) {
                teta = (it-it1)*dteta;
                wscal = 0.5 + 0.5*cos(teta);
                traces[ipy*nxmax*ntmax+ipx*ntmax+it] = traces[ipy*nxmax*ntmax+ipx*ntmax+it]*wscal;
            }
            for (it=it2; it<nt; it++) {
                traces[ipy*nxmax*ntmax+ipx*ntmax+it] = 0.;
            }
        }
        
        for (it=0;it<nt;it++) {
            for (ipx=0; ipx<npx; ipx++) {
                if(fabs(traces[ipy*nxmax*ntmax+ipx*ntmax+it]) > 0.) {
                    ipx2 = ipx + iptx;
                    if(ipx2 > npx) ipx2 = npx;
                    for (ix=ipx; ix<ipx2; ix++) {
                        teta = (ix-ipx)*dteta2;
                        wscal = 0.5 - 0.5*cos(teta);
                        traces[ipy*nxmax*ntmax+ix*ntmax+it] = traces[ipy*nxmax*ntmax+ix*ntmax+it]*wscal;
                    }
                    break;
                }
            }
        }
        
        for (it=0;it<nt;it++) {
            for (ipx=npx-1; ipx>=0; ipx--) {
                if(fabs(traces[ipy*nxmax*ntmax+ipx*ntmax+it]) > 0.) {
                    ipx2 = ipx - iptx;
                    if(ipx2 < 0) ipx2 = 0;
                    for (ix=ipx; ix>ipx2; ix--) {
                        teta = (ipx-ix)*dteta2;
                        wscal = 0.5 - 0.5*cos(teta);
                        traces[ipy*nxmax*ntmax+ix*ntmax+it] = traces[ipy*nxmax*ntmax+ix*ntmax+it]*wscal;
                    }
                    break;
                }
            }
        }
    }
    
    
    
    dteta2 = PI/ipty;
    ddpx = vpymax - pymin;
    if(fabs(ddpx) < dpy) ddpx = dpy;
    for (ipx=0; ipx<npx; ipx++) {
        for (ipy=0; ipy<npy; ipy++) {
            pyy = pymin + ipy*dpy;
            if(fabs(pyy) < fabs(vpymax)) continue;
            it0 = nt*((-fabs(pyy) - pymin)/ddpx);
            it1 = it0-itx;
            if(it1 < 0) it1 = 0;
            it2 = it0+itx;
            if(it2 > nt ) it2 = nt;
            
            for (it=it1; it<it2; it++) {
                teta = (it-it1)*dteta;
                wscal = 0.5 + 0.5*cos(teta);
                traces[ipy*nxmax*ntmax+ipx*ntmax+it] = traces[ipy*nxmax*ntmax+ipx*ntmax+it]*wscal;	     }
            for (it=it2; it<nt; it++) {
                traces[ipy*nxmax*ntmax+ipx*ntmax+it] = 0.;
            }
        }
        
        for (it=0;it<nt;it++) {
            for (ipy=0; ipy<npy; ipy++) {
                if(fabs(traces[ipy*nxmax*ntmax+ipx*ntmax+it]) > 0.) {
                    ipy2 = ipy + ipty;
                    if(ipy2 > npy) ipy2 = npy;
                    for (iy=ipy; iy<ipy2; iy++) {
                        teta = (iy-ipy)*dteta2;
                        wscal = 0.5 - 0.5*cos(teta);
                        traces[iy*nxmax*ntmax+ipx*ntmax+it] = traces[iy*nxmax*ntmax+ipx*ntmax+it]*wscal;
                    }
                    break;
                }
            }
        }
        
        for (it=0;it<nt;it++) {
            for (ipy=npy-1; ipy>=0; ipy--) {
                if(fabs(traces[ipy*nxmax*ntmax+ipx*ntmax+it]) > 0.) {
                    ipy2 = ipy - ipty;
                    if(ipy2 < 0) ipy2 = 0;
                    for (iy=ipy; iy>ipy2; iy--) {
                        teta = (ipy-iy)*dteta2;
                        wscal = 0.5 - 0.5*cos(teta);
                        traces[iy*nxmax*ntmax+ipx*ntmax+it] = traces[iy*nxmax*ntmax+ipx*ntmax+it]*wscal;
                    }
                    break;
                }
            }
        }
    }
    
    
    // *********************************** End Muting .....
    
    
    
    /* loop over traces in tau-p domain */
    
    for (ipy=0; ipy<npy; ipy++) {
        
        for (ipx=0; ipx<npx; ipx++) {
            
            /* pad tau axis with zeros */
            for (itau=0; itau<ntau; itau++)
                tr_fft[itau]=traces[ipy*ntmax*nxmax+ipx*ntmax+itau];
            for (itau=ntau; itau<ntaufft; itau++)
                tr_fft[itau]=0.0;
            
            /* Fourier transform tau to frequency */
            
            pfarc(1, ntaufft, tr_fft, help);
            
            for (it=0; it<ntw; it++) {
                w = it*dw;
                ctr[ipy*nxmax*ntmax+ipx*ntmax+2*it]   = help[it].r;
                ctr[ipy*nxmax*ntmax+ipx*ntmax+2*it+1] = help[it].i;
            }
            
        }
    }
    
    
    // Second in Crossline Direction
    
    for(ipx=0;ipx<npx;ipx++) {
        
        /* loop over w */
        for (iw=0, w=fw; iw<ntw; ++iw, w+=dw) {
            
            /* scale ctr by p sampling interval */
            for (ipy=0; ipy<npy; ipy++) {
                tr_p[ipy].r = ctr[ipy*nxmax*ntmax+ipx*ntmax+2*iw];
                tr_p[ipy].i = ctr[ipy*nxmax*ntmax+ipx*ntmax+2*iw+1];
            }
            
            /* compute number of k's for interpolation */
            //nkx = (pmax-pmin)*w/dk + 1;
            nky = nkyfft;
            /* pk = alloc1float(nk); */
            
            /* compute p values at which to interpolate tr_p */
            for (iky=0, k=fky; iky<nky; iky++, k+=dky) {
                if (w==0.0) pk[iky]= k/dw;
                else pk[iky] = k/w;
            }
            
            /* interpolate tr_pa to obtain tr_k */
            ints8c(npy, dpy, fpy, tr_p, czero, czero, nky, pk, hk);
            if(w==0.0) {
                for (iky=0;iky<nky; iky++)
                    hk[iky] = czero;
            }
            
            xshift = ymin;
            for (iky=0, k=fky; iky<nky; iky++, k+=dky) {
                phase = k*xshift;
                c = cos(phase);
                s = sin(phase);
                temp = hk[iky].r*c-hk[iky].i*s;
                hk[iky].i = hk[iky].r*s+hk[iky].i*c;
                hk[iky].r=temp;
            }
            
            
            /* inverse Fourier transform from k to x */
            pfacc(1, nkyfft, hk);
            
            for (iky=1; iky<nkyfft; iky+=2) {
                hk[iky].r = -hk[iky].r;
                hk[iky].i = -hk[iky].i;
            }
            
            
            /* copy 1-D interpolated array to 2-D array */
            for (ipy=0; ipy<ny; ipy++) {
                ctr[ipy*nxmax*ntmax+ipx*ntmax+2*iw] =   hk[ipy].r;
                ctr[ipy*nxmax*ntmax+ipx*ntmax+2*iw+1] = hk[ipy].i;
            }
        }
    }
    
    
    
    // First in Inline Direction
    
    for(ipy=0;ipy<ny;ipy++) {
        
        /* loop over w */
        for (iw=0, w=fw; iw<ntw; ++iw, w+=dw) {
            
            /* scale ctr by p sampling interval */
            for (ipx=0; ipx<npx; ipx++)
                //	tr_p[ipx] = crmul(ctr[ipx*nw+iw],dpx*fftscl);
            {
                tr_p[ipx].r = ctr[ipy*nxmax*ntmax+ipx*ntmax+2*iw];
                tr_p[ipx].i = ctr[ipy*nxmax*ntmax+ipx*ntmax+2*iw+1];
            }
            
            /* compute number of k's for interpolation */
            //nkx = (pmax-pmin)*w/dk + 1;
            nkx = nkxfft;
            /* pk = alloc1float(nk); */
            
            /* compute p values at which to interpolate tr_p */
            for (ikx=0, k=fkx; ikx<nkx; ikx++, k+=dkx) {
                if (w==0.0) pk[ikx]= k/dw;
                else pk[ikx] = k/w;
            }
            
            /* interpolate tr_pa to obtain tr_k */
            ints8c(npx, dpx, fpx, tr_p, czero, czero, nkx, pk, hk);
            if(w==0.0) {
                for (ikx=0;ikx<nkx; ikx++)
                    hk[ikx] = czero;
            }
            
            xshift = xmin;
            for (ikx=0, k=fkx; ikx<nkx; ikx++, k+=dkx) {
                phase = k*xshift;
                c = cos(phase);
                s = sin(phase);
                temp = hk[ikx].r*c-hk[ikx].i*s;
                hk[ikx].i = hk[ikx].r*s+hk[ikx].i*c;
                hk[ikx].r=temp;
            }
            
            
            /* inverse Fourier transform from k to x */
            pfacc(1, nkxfft, hk);
            
            for (ikx=1; ikx<nkxfft; ikx+=2) {
                hk[ikx].r = -hk[ikx].r;
                hk[ikx].i = -hk[ikx].i;
            }
            
            
            /* copy 1-D interpolated array to 2-D array */
            for (ix=0; ix<nx; ix++) {  //33 changed npx to nx 
                ctr[ipy*nxmax*ntmax+ix*ntmax+2*iw]=hk[ix].r;
                ctr[ipy*nxmax*ntmax+ix*ntmax+2*iw+1]=hk[ix].i;
            }
        }
    }
    
    
    
    /* inverse Fourier transform from frequency to time */
    
    for (iy=0; iy<ny; iy++) {
        for (ix=0; ix<nx; ix++) {
            for (it=0; it<ntw; it++) {
                help[it].r = ctr[iy*nxmax*ntmax+ix*ntmax+2*it];
                help[it].i = ctr[iy*nxmax*ntmax+ix*ntmax+2*it+1];
            }
            
            for (it=ntw; it<nw; it++)
                help[it] = czero;
            
            pfacr(-1, ntaufft, help, tr_fft);
            for (it=0; it<nt; it++)
                traces[iy*nxmax*ntmax+ix*ntmax+it]=tr_fft[it]/ntaufft;
        }
    }
    
    return;
}
