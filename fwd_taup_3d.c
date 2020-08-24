/*
 *
 **********************************************************************************************
 *
 * Compute global forward slant stack (tau-p transform) via
 * FK  transform
 *
 * Input is:
 *
 * ---------------------------------------------------------------------------------------------*
 * float *traces:    Input trace data (nx*ny*nt)                                                *
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
 * float *trace:   tau-p transformed data                                                       *
 *
 * ---------------------------------------------------------------------------------------------*
 *
 *
 */

void fwd_taup_3d_(float *traces, int *ntt, int *nttmax, int *nxxmax, int *nyymax, int *nxx, int *nyy, float *dtt, float *dxx,
        float *dyy, float *xxmin, float *yymin, int *npxx, int *npyy, float *pxxmin, float *pyymin,
        //           float *maxfreqq, int *nxxfft, int *nyyfft, float *work)
        float *pxxmax, float *pyymax, float *maxfreqq, int *nxxfft, int *nyyfft, float *slice, float *work)
        
{
    int i, it, itau, iy, ix, iw, ik, ip;  /* loop counters */
    float fx, fw;		/* first frequency */
    float fpx;            	/* first slope */
    float dw, dkx;            /* frequency and wavenumber sampling interval */
    int nnw, nw, nk;           /* number of frequencies and wavenumbers */
    int ntpad;              /* number of time samples to pad */
    int ntfft;              /* number of padded time samples */
    int ntau;		/*  */
    /* int lk; */		/* last wavenumber */
    int nxfft;		/* number of x samples to pad */
    float fftscl;		/* scale factor for FFT */
    float w, p, k;            /* frequency, slope and wavenumber */
    float phase;            /* phase angle */
    float xshift;           /* shift to account for non-centered x-axis */
    float c, s;              /* auxiliary variables */
    float xmax;             /* maximum horizontal value */
    float fkx;          	/* first wavenumber */
    float pxmax;             /* maximum slope */
    int lwrap;              /* samples required to make k periodic */
    float fkax;              /* first wavenumber with k periodic */
    float temp;		/* auxiliary variable */
    int nka;                /* number of wavenumbers with k perioic */
    float *tr_fft;          /* padded trace for FFT */
    complex *tr_x;          /* W-transformed input tracescaled  by dx */
    float *kp;		/* K-transformed slant stacked single trace */
    complex *hp;            /* slant stacked single trace */
    float *ctr;          /* F-K transformed trace */
    float *ctr_p;        /*  */
    complex *tr_k;          /* K-transformed trace  */
    complex *tr_ka;         /*  */
    complex czero;		/* complex number zero */
    complex *tr_xy;
    int ipx, ipy, nkx, nyfft;
    int nky, nkax, nkay, ntxpad, ntypad;
    float ymax, pymax, fky, fkay, dky, fpy, fy, maxfreq;
    complex *cwork, *help, *ole;         /* complex working area */
    int nxmax, nymax, ntmax, ntw, nwork;
    
    float dt, xmin, ymin, dx, dy, pxmin, pymin, dpx, dpy, pxx, pyy;
    int nt, nx, ny, npx, npy, nnw2, it0, it1, it2, itx, ipx2, ipy2, iptx, ipty;
    float wscal, maxfreq2, dteta, dteta2, teta, vpxmax, vpymax, ddpx;
    
    /* Convert from input pointers */
    
    nxmax = *nxxmax;
    nymax = *nyymax;
    dt = *dtt;
    nt = *ntt;
    ntmax = *nttmax;
    ny = *nyy;
    nx = *nxx;
    ymin = *yymin*0.001;
    xmin = *xxmin*0.001;
    dy = *dyy*0.001;
    dx = *dxx*0.001;
    npy = *npyy;
    npx = *npxx;
    pymin = *pyymin; // *1000.;
    pxmin = *pxxmin; //*1000.;
    dpy = 2.*fabs(pymin)/(npy-1);
    dpx = 2.*fabs(pxmin)/(npx-1);
    maxfreq = *maxfreqq;
    ntw = ntmax/2;
    vpxmax = *pxxmax; //*1000;
    vpymax = *pyymax; //*1000;
    if(vpxmax < pxmin)  vpxmax = pxmin;
    if(vpymax < pymin)  vpymax = pymin;
    
    czero = cmplx(0.0, 0.0);
    
    
    /* compute slope sampling interval */
    pymax = pymin+(npy-1)*dpy;
    pxmax = pxmin+(npx-1)*dpx;
    czero = cmplx(0.0, 0.0);
    
    
    /* determine lengths and scale factors for FFT's */
    fx=0.0;
    fy=0.0;
    fpy=pymin+0.5*(pymax-pymin-(npy-1)*dpy);
    fpx=pxmin+0.5*(pxmax-pxmin-(npx-1)*dpx);
    pymax = (dpy<0.0)?fpy:fpy+(npy-1)*dpy;
    pxmax = (dpx<0.0)?fpx:fpx+(npx-1)*dpx;
    ymax = (dy<0.0)?fy:fy+(ny-1)*dy;
    xmax = (dx<0.0)?fx:fx+(nx-1)*dx;
    ntau = nt;
    ntxpad = ABS(pxmax*xmax)/dt;
    ntypad = ABS(pymax*ymax)/dt;
    ntpad = ntxpad;
    if(ntpad < ntypad) ntpad = ntypad;
    ntfft = npfar(MAX(nt+ntpad, ntau));
    fftscl = 1.0/ntfft;
    lwrap = 16;      /* samples required to make k periodic */
    nxfft = npfa((nx+lwrap));
    nyfft = ny;
    if(ny > 1) nyfft = npfa((ny+lwrap));
    // ntfft = npfar(MAX(nt+ntpad,ntau));
    ntfft = npfar(1.5*nt);
    
    nxfft = *nxxfft;
    nyfft = *nyyfft;
    
    
    /* determine frequency and wavenumber sampling */
    nw = ntfft/2+1;
    dw = 2.0*PI/(ntfft*dt);
    fw = 0.0;
    nkx = nxfft;
    nky = nyfft;
    nk = nxfft;
    if(nky > nkx) nk = nyfft;
    nky = nyfft;
    dky = 2.0*PI/(nyfft*dy);
    dkx = 2.0*PI/(nxfft*dx);
    fky = -PI/dy;
    fkx = -PI/dx;
    /* lk = PI/dx; */
    fkay = fky-lwrap*dky;
    fkax = fkx-lwrap*dkx;
    nkay = lwrap+nky+lwrap;
    nkax = lwrap+nkx+lwrap;
    
    
    if(ntw > nw) ntw = nw;
    
    nwork = 4*ntfft + 4*nkax + 6*nxfft*nyfft;
    
    tr_fft = (float*)work;
    nwork = ntfft;
    kp = (float*)work + nwork;
    nwork = ntfft + npx;
    
    cwork = (complex*)(work+nwork);
    ctr = (float*)traces;
    ctr_p = (float *)traces;
    tr_ka = (complex*)cwork;
    tr_x = tr_k=tr_ka+lwrap;        /* pointers to tr_ka[lwrap] */
    nwork = nwork + nkax;
    hp = (complex*)(cwork+nwork);
    nwork = nwork + nx + ny;
    help = (complex*)(cwork+nwork);
    nwork = nwork + 2*ntfft;
    
    nwork = 2*nxfft*nyfft;
    tr_xy = (complex*)(slice);
    ole = (complex*)(slice+nwork);
    
    
    /* loop over traces */
    
    nnw = 2*PI*maxfreq/dw + 1;
    maxfreq2 = maxfreq*1.2;
    if(0.2*maxfreq < 10.) maxfreq2 = maxfreq + 10.;
    nnw2 = 2.*PI*maxfreq2/dw + 1;
    dteta = PI/(nnw2-nnw-1);
    
    if(nnw > ntw) nnw = ntw;
    if(nnw2 > ntw) nnw2 = ntw;
    
    
    for (iy=0; iy<ny; iy++) {
        for (ix=0; ix<nx; ix++) {
            
            /* pad time axis with zeros */
            for (it=0; it<nt; it++)
                tr_fft[it]=traces[iy*nxmax*ntmax+ix*ntmax+it];
            for (it=nt; it<ntfft; it++)
                tr_fft[it]=0.0;
            
            /* Fourier transform time to frequency */
            pfarc(1, ntfft, tr_fft, help);
            
            for (it=0; it<nnw; it++) {
                ctr[iy*nxmax*ntmax+ix*ntmax+2*it]   = help[it].r;
                ctr[iy*nxmax*ntmax+ix*ntmax+2*it+1] = help[it].i;
            }
            
            for (it=nnw; it<nnw2; it++) {
                teta = (it-nnw)*dteta;
                wscal = 0.5 + 0.5*cos(teta);
                ctr[iy*nxmax*ntmax+ix*ntmax+2*it]   = wscal*help[it].r;
                ctr[iy*nxmax*ntmax+ix*ntmax+2*it+1] = wscal*help[it].i;
            }
            
            for (it=nnw2; it<ntw; it++) {
                ctr[iy*nxmax*ntmax+ix*ntmax+2*it]   = 0.;
                ctr[iy*nxmax*ntmax+ix*ntmax+2*it+1] = 0.;
            }
            
        }
    }
    
    // loop over w
    
    for (iw=0, w=fw; iw<nnw2; ++iw, w+=dw) {
        
        //  scale tr_x by x sampling interval
        
        for(iy=0; iy<ny; iy++) {
            for (ix=0; ix<nx; ix++) {
                tr_xy[iy*nxfft+ix].r = ctr[iy*nxmax*ntmax+ix*ntmax+2*iw];
                tr_xy[iy*nxfft+ix].i = ctr[iy*nxmax*ntmax+ix*ntmax+2*iw+1];
            }
            
            if(npx > 1) {
                for (ix=nx; ix<nxfft; ix++) {
                    tr_xy[iy*nxfft+ix].r = 0.;
                    tr_xy[iy*nxfft+ix].i = 0.;
                }
            }
        }
        
        if(npy > 1) {
            for(iy=ny; iy<nyfft; iy++) {
                for (ix=0; ix<nxfft; ix++) {
                    tr_xy[iy*nxfft+ix].r = 0.;
                    tr_xy[iy*nxfft+ix].i = 0.;
                }
            }
        }
        
        
        fk_trans(ctr_p, tr_xy, tr_x, tr_k, tr_ka, hp, ole, kp, ny, nx, nyfft, nxfft, dky, dkx, fkay, fkax, nky, nkx,
                dy, dx, npy, npx, ntmax, fpy, fpx, dpy, dpx, fy, fx, iw, w, nkay, nkax, ntfft, ymin, xmin, nxmax);
        
    }
    
    // loop over p
    
    
    for (ipy=0; ipy<npy; ipy++) {
        
        for (ipx=0; ipx<npx; ipx++) {
            
            // Fourier transform frequency to time
            
            for (it=0; it<ntw; it++) {
                help[it].r = ctr_p[ipy*nxmax*ntmax+ipx*ntmax+2*it];
                help[it].i = ctr_p[ipy*nxmax*ntmax+ipx*ntmax+2*it+1];
            }
            
            for (it=ntw; it<nw; it++)
                help[it] = czero;
            
            pfacr(-1, ntfft, help, tr_fft);
            
            
            // copy to output array
            for (itau=0; itau<ntau; itau++)
                traces[ipy*nxmax*ntmax+ipx*ntmax+itau] = tr_fft[itau]/ntfft;
        }
    }
    
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
    
    return;
}



void fk_trans(float *ctr_p, complex *tr_xy, complex *tr_x, complex *tr_k, complex *tr_ka, complex *hp, complex *ole,
        float *kp, int ny, int nx, int nyfft, int nxfft, float dky, float dkx, float fkay, float fkax,
        int nky, int nkx, float dy, float dx, int npy, int npx, int ntmax, float fpy, float fpx,
        float dpy, float dpx, float fy, float fx, int iw, float w, int nkay, int nkax, int ntfft,
        float ymin, float xmin, int nxmax) {
    int iy, ip, ik, ix;
    int lwrap = 16;
    float xshift, phase, temp, k, fftscl;
    float c, s, p;
    complex czero;
    czero = cmplx(0., 0.);
    fftscl = 1.0/ntfft;
    
    // Fourier transform tr_x to tr_k
    
    for(iy=0;iy<nyfft;iy++) {
        for(ix=0;ix<nxfft;ix++) {
            tr_x[ix] = tr_xy[iy*nxfft+ix];
            if(ix%2 == 1) {
                tr_x[ix].r = -tr_x[ix].r;
                tr_x[ix].i = -tr_x[ix].i;
            }
            
        }
        
        pfacc(-1, nxfft, tr_x);
        
        for (ix=0;ix<nxfft;ix++) {
            ole[iy*nxfft+ix] = tr_x[ix];
        }
    }
    
    
    // Running over inline line  direction first .....
    
    for(iy=0;iy<nyfft;iy++) {
        for(ix=0;ix<nxfft;ix++) {
            tr_x[ix] = ole[iy*nxfft+ix];
        }
        
        
        // wrap-around tr_k to avoid interpol end effects
        for (ik=0; ik<lwrap; ik++)
            tr_ka[ik] = tr_k[ik+nkx-lwrap];
        for (ik=lwrap+nkx; ik<lwrap+nkx+lwrap; ik++)
            tr_ka[ik] = tr_k[ik-lwrap-nkx];
        
        
        xshift = -xmin;
        for (ik=0, k=fkax; ik<nkax; ik++, k+=dkx) {
            phase = k*xshift;
            c = cos(phase);
            s = sin(phase);
            temp = tr_ka[ik].r*c-tr_ka[ik].i*s;
            tr_ka[ik].i = tr_ka[ik].r*s+tr_ka[ik].i*c;
            tr_ka[ik].r=temp;
        }
        
        
        //  compute k values at which to interpolate tr_k
        for (ip=0, p=fpx; ip<npx; ip++, p+=dpx) {
            kp[ip] = w*p;
        }
        
        
        // 8-point sinc interpolation of tr_k to obtain h(p)
        ints8c(nkax, dkx, fkax, tr_ka, czero, czero, npx, kp, hp);
        
        // phase shift to account for non-centered x-axis
        for (ip=0; ip<npx; ip++) {
            ole[iy*nxfft+ip].r = hp[ip].r;
            ole[iy*nxfft+ip].i = hp[ip].i;
        }
        
    }
    
    
    // Running over cross line  direction second .....
    
    for(ix=0;ix<npx;ix++) {
        for(iy=0;iy<nyfft;iy++) {
            tr_x[iy].r = ole[iy*nxfft+ix].r;
            tr_x[iy].i = ole[iy*nxfft+ix].i;
            if(iy%2 == 1) {
                tr_x[iy].r = -tr_x[iy].r;
                tr_x[iy].i = -tr_x[iy].i;
            }		      }
        
        pfacc(-1, nyfft, tr_x);
        
        for (iy=0;iy<nyfft;iy++) {
            ole[iy*nxfft+ix] = tr_x[iy];
        }
    }
    
    
    
    // Running over cross line  direction second .....
    
    for(ix=0;ix<npx;ix++) {
        for(iy=0;iy<nyfft;iy++) {
            tr_x[iy] = ole[iy*nxfft+ix];
        }
        
        
        // wrap-around tr_k to avoid interpol end effects
        for (ik=0; ik<lwrap; ik++)
            tr_ka[ik] = tr_k[ik+nky-lwrap];
        for (ik=lwrap+nky; ik<lwrap+nky+lwrap; ik++)
            tr_ka[ik] = tr_k[ik-lwrap-nky];
        
        
        xshift = -ymin;
        for (ik=0, k=fkay; ik<nkay; ik++, k+=dky) {
            phase = k*xshift;
            c = cos(phase);
            s = sin(phase);
            temp = tr_ka[ik].r*c-tr_ka[ik].i*s;
            tr_ka[ik].i = tr_ka[ik].r*s+tr_ka[ik].i*c;
            tr_ka[ik].r=temp;
        }
        
        
        // compute k values at which to interpolate tr_k
        for (ip=0, p=fpy; ip<npy; ip++, p+=dpy) {
            kp[ip] = w*p;
        }
        
        // 8-point sinc interpolation of tr_k to obtain h(p)
        ints8c(nkay, dky, fkay, tr_ka, czero, czero, npy, kp, hp);
        
        for (ip=0; ip<npy; ip++) {
            ctr_p[ip*nxmax*ntmax+ix*ntmax+2*iw] = hp[ip].r;
            ctr_p[ip*nxmax*ntmax+ix*ntmax+2*iw+1] = hp[ip].i;
        }
        
    }
    
    return;
}



