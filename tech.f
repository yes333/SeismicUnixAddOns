c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  ccorr(a,b,length,cmax)
      real  a(length), b(length)
c
c     compute crosscorrelation 
c

      cmax = 0
      ccsum=0
      if(length.lt.1)return

      do i=1,length
         ccsum=ccsum+a(i)*b(i)
      enddo
      cmax=ccsum/length
      return
      END



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  acorr(a,b,length,cmax)
      real  a(length), b(length)
c
c     compute crosscorrelation dvided by autocorrelation
c
      cmax = 0
      ccsum=0
      aasum=0
      do i=1,length
         ccsum=ccsum+a(i)*b(i)
         aasum=aasum+b(i)*b(i)
      enddo
      if(aasum.gt.0.0)cmax=ccsum/aasum
      return
      END

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  sdev(a,n,aver,stdev)
      real  a(n)
c
c     Compute standard deviation
c
      if(n.lt.2)return
      stdev = 0

      aasum=0
      do i=1,n
         aasum=aasum + (a(i)-aver)**2
      enddo
      if(aasum.gt.0.0)stdev=sqrt(aasum/(n-1))
      return
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c 
       subroutine hamm_weights(coswg , nktaper)
       real coswg(nktaper)
c
c      compute Hamming weights
c
c
       pi   = 4*atan(1.0)
      

       dv = pi/nktaper
       v = pi
       do itap = 1, nktaper
          v    = v - dv
          coswg(itap) = (0.5 + 0.5*cos(v))
       enddo
       R E T U R N
       END





c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c 
      subroutine fwdfft( dtabfr, nsamp , fftdata , ksamp )
c
c
c
c     Compute Fourier transforms of input data 
c
c
c     dtabfr -  Contains the original X,t domain data
c     fftdata  -  Will contain z-comp fourier transform for different windows
 

c     ** input time array
c
      real  dtabfr(nsamp)
      real fftdata(ksamp)
 

      do isamp = 1,nsamp
         fftdata(isamp) = dtabfr(isamp)
      enddo
      
      
      do isamp = nsamp+1,ksamp
         fftdata(isamp) = 0.0
      enddo
      
      call rwsrfft(fftdata,ksamp,+1)

      return
      END
c     end of fwdfft



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
c
      subroutine set_sinclen(method,nfilt)
      character*(*) method


c     ******************  SINC BLOCK  ******************
      parameter (maxfunc=32 , nsinc=32 )
      real  tab(-nsinc+1 : nsinc , 1:maxfunc)
      real  sp (-nsinc+1:nsinc)
      character*8 cmethod
      common /csinc/ cmethod , nsinc2 , nfunc , tab , sp
c*******************************************************


      nfunc  = min(nfilt,maxfunc)
      nsinc2 = min(nfilt,nsinc)/2

      cmethod = method

      if(cmethod.eq.'TABLE') then
        
        call sinctable( tab(-nsinc2+1,1),nsinc2,nfunc) 

      endif

      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine sinctable( tab , nsinc2,nfunc) 
      real tab(-nsinc2+1:nsinc2 , 1:nfunc)

      real*8 aint,adif ,value1,value2,pi,fi,w

      data pi/3.141592654d00/

      nfilt = nsinc2*2
      aint  = dble(1) / dble(nfilt-1)

      do j=1,nfunc

         adif = (j-1)*aint

         do i = -nsinc2+1 , nsinc2                 
            fi  = ( adif - float(i) ) * pi
            if(fi.eq.0.0) then
               tab(i,j) = 1
            else
               w  = 0.5+0.5*dcos(fi/dble(nsinc2)) ! hamming weight
               tab(i,j) =  w*dsin(fi)/fi         
            endif
         end do 

      enddo

      return
      end







c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  invsinc(value,arr1,lsamp,it0,tdif)
      real arr1(lsamp) 

c     it0= The closest sample on the left side
c     tdif - The time difference {ms}/sr from the value to it0
c
c     Do the inverse sinc operation 
c     
c     Given VAlUE at time it0+tdif 
c     compute the response at it0-3 , it0-2 ... it0+4
c
c     (Currently hard coded to 16 points)
c
      real*8 fi,pi

      data pi/3.141592654d00/

      data small/0.001/

      data nsinc2/8/


      if (tdif .lt.small) then
         arr1(it0) = arr1(it0) + value  
         return
      elseif(1.0 - tdif.lt.small) then
         arr1(it0+1) = arr1(it0+1) + value 
         return 
      endif



c     ** reduce sinc length
c
      il1=nsinc2-1
 20   if (it0.lt.il1) then
         il1=il1-1
         goto 20
      endif

      il2=nsinc2
 21   if (it0+il2.gt.lsamp) then
         il2=il2-1
         goto 21
      endif
 

      do 10 i=-il1,il2
         fi = ( tdif - float(i) ) * pi
         w  = 0.5+0.5*dcos(fi/dble(nsinc2)) ! hamming weight

         arr1(it0+i) = arr1(it0+i) + value *w*dsin(fi)/fi
   10 continue

      return
      end





c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  rsam(aux , lsamp , sr , arr1 , maxlen , si)
      real aux(lsamp)    ! Input array 
      real arr1(maxlen)  ! output array
      real sr            ! sample rate in
      real si            ! sample rate out 


c     ******************  SINC BLOCK  ******************
      parameter (maxfunc=32 , nsinc=32 )
      real  tab(-nsinc+1 : nsinc , 1:maxfunc)
      real  sp (-nsinc+1:nsinc)
      character*8 cmethod
      common /csinc/ cmethod , nsinc2 , nfunc , tab , sp
c*******************************************************


      double precision ti1,tk , dsi , dsr

      logical m

      if( (cmethod.ne.'TABLE'.and.cmethod.ne.'EXACT').or.
     $    (nsinc2.lt.4.or.nsinc2.gt.16) ) then

c        ** default        
         call set_sinclen('TABLE',16)
      endif




      m = cmethod.eq.'TABLE'


      dsi = dble(si)
      dsr = dble(sr)


c     ** Resample trace
c     ** Time domain process
c
 
      k=0 

   50 continue
      k = k+1
      if(k.gt.maxlen) goto 100
      tk = (k-1)  * dsi

      i1=  tk/sr + 1.

      if (i1.gt.lsamp) goto 100

      ti1 = (i1-1)*dsr

      tdif = (tk-ti1) / dsr
      if(m) then
         arr1(k) = tabsinc(aux,lsamp,i1,tdif)
      else  
         arr1(k) = sinc(aux,lsamp,i1,tdif)
      endif
      goto 50
  100 continue

      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  tabrsam(inparr,lsamp,sr,outarr,maxlen,si,
     $                    tab , nsinc2 , nfunc)
      real inparr(lsamp) , outarr(maxlen)
      real tab(-nsinc2+1:nsinc2,1:nfunc) 

      double precision ti1,tk , dsi , dsr,value1,value2,w , tdif

      data small/0.0001/


      dsi = dble(si)
      dsr = dble(sr)

      do k = 1 , maxlen 
         tk = (k-1)  * dsi

         i1=  INT (tk/dsr + dble(1) )

         if (i1.le.lsamp) then

            ti1 = (i1-1)*dsr

            tdif = (tk-ti1) / dsr
            if(tdif.gt.dble(1)) tdif=dble(1)
            if(tdif.lt.dble(0)) tdif=dble(0)

            it0 = i1

            if (tdif .lt.small) then
               outarr(k) = inparr(i1)
            elseif(1.0 - tdif.lt.small) then
               outarr(k) = inparr(i1+1)
            else 

               il1=nsinc2-1
               il2=nsinc2

               do while (it0.le.il1)
                  il1 = il1 - 1
               enddo

               do while (it0+il2.gt.lsamp)
                  il2 = il2 - 1
               enddo

               aint = dble(1)/dble(nsinc2*2-1)

               KK = tdif/aint


c              ** weight at right  sample
               w = (tdif - KK*aint ) /aint

               jval1 = KK +1
               jval2 = jval1 + 1
               value1 = 0.0
               value2 = 0.0

c              ** interpolate
c     
               il1=-1*il1
               do i=il1,il2
                  value1 = value1 + inparr( it0 + i) * tab(i,jval1)
                  value2 = value2 + inparr( it0 + i) * tab(i,jval2)
               end do

               outarr(k) = value1*(1.-w) + value2*w
            endif
         endif

      enddo



  100 continue

      return
      end


c
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      real function  tabsinc(arr1,lsamp,it0,tdif)
      real arr1(lsamp) 




cc
c     tabulated interpolators
c

c     it0= The closest sample on the left side
c     tdif - The time difference {ms}/sr from the value to it0
c     Precomputation of interpolators

c     ******************  SINC BLOCK  ******************
      parameter (maxfunc=32 , nsinc=32 )
      real  tab(-nsinc+1 : nsinc , 1:maxfunc)
      real  sp (-nsinc+1:nsinc)
      character*8 cmethod
      common /csinc/ cmethod , nsinc2 , nfunc , tab , sp
c*******************************************************


      real*8 aint, value1,value2


      data small/0.0001/
      data icall/0/
      data nerr/0/

      tabsinc = -999.25
      if(nsinc2.lt.4.or.nsinc2.gt.16) then
         write(6,*) 'Error in tabsinc: Not initialised'
         return
      endif


      if (tdif.gt.1.0.or.tdif.lt.0.or.it0.lt.1.or.it0.gt.lsamp) then
         nerr=nerr+1
         write (6,*)  ' Sinc error:',tdif,it0,lsamp
         if (nerr.gt.10) then
           stop ' 10 consecutive errors'
         endif
      elseif (tdif .lt.small) then
         tabsinc=arr1(it0)
         return
      elseif(1.0 - tdif.lt.small) then
         tabsinc=arr1(it0+1)
         return 
      endif

      nerr  =   0

      il1=nsinc2-1
      il2=nsinc2

      do while (it0.le.il1)
        il1 = il1 - 1
      enddo


      do while (it0+il2.gt.lsamp)
        il2 = il2 - 1
      enddo


      aint = dble(1)/dble(nsinc2*2-1)

      KK = tdif/aint


c     ** weight at right  sample
      w = (tdif - KK*aint ) /aint

      jval1 = KK +1
      jval2 = jval1 + 1
      value1 = 0.0
      value2 = 0.0

c     ** interpolate
c
      il1=-1*il1
      do i=il1,il2
        value1 = value1 + arr1( it0 + i) * tab(i,jval1)
        value2 = value2 + arr1( it0 + i) * tab(i,jval2)
      end do

      tabsinc = value1*(1.-w) + value2*w
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      real function  sinc(arr1,lsamp,it0,tdif)
      real arr1(lsamp) 
c
c
c     "EXACT" interpolator
c
c
c

c
c     it0= The closest sample on the left side
c     tdif - The time difference {ms}/sr from the value to it0
c

c     ******************  SINC BLOCK  ******************
      parameter (maxfunc=32 , nsinc=32 )
      real  tab(-nsinc+1 : nsinc , 1:maxfunc)
      real  sp (-nsinc+1:nsinc)
      character*8 cmethod
      common /csinc/ cmethod , nsinc2 , nfunc , tab , sp
c*******************************************************




 
      real*8 pi , fi , w , value
      data pi/3.141592654d00/

      data small/0.001/
      data otdif/-99.0/
      
      data nerr/0/

c     ** default
      if(nsinc2.lt.4.or.nsinc2.gt.16)call set_sinclen('EXACT',16)


      if (tdif .lt.small) then
         sinc=arr1(it0)
         return
      elseif(1.0 - tdif.lt.small) then
         sinc=arr1(it0+1)
         return 
      elseif ( ABS(tdif).gt.1.0 .or. it0.lt.1.or.it0.gt.lsamp) then
         nerr=nerr+1
         if (nerr.gt.10) then
           write (6,*)  ' Sinc error:',tdif,it0,lsamp
           stop ' sinc 10 consecutive errors'
         endif
      endif

      nerr  =   0
      value = 0.0

      il1 = nsinc2-1
      il2 = nsinc2


c     ** adjust start and stop 
c

      do while (it0.le.il1)
        il1 = il1 - 1
      enddo


      do while (it0+il2.gt.lsamp)
        il2 = il2 - 1
      enddo


c     ** compute sinc interpolator
c        
      if(tdif .ne. otdif) then        
        otdif = tdif                  
        do i = -nsinc2 + 1 , nsinc2                  
         w  = 0.5+0.5*dcos(fi/dble(nsinc2))   ! hamming weight
         fi  = ( tdif - dble(i) ) * pi
         sp(i) =  w * dsin(fi)/fi         
        end do 
      endif

c     ** interpolate
c
      il1=-1*il1
      do i=il1,il2
        value = value + arr1( it0 + i) * sp(i)
      end do

      sinc = value
      return
      end




c******************************************************************
c******************************************************************
c     title :           prfft
c
c     keywords :        apsimlib
c
c     language :        fortran 77
c
c     description :     real to complex forward fft
c                       complex to real inverse fft
c                       both in-place
c
c     calling sequence :
c
c       call rwsrfft(sig, nreal, flag)
c
c     arguments :
c
c       sig   = complex : input / output array
c
c       nreal = integer : no. of real samples in array
c                         (must be a power of 2)
c
c       flag  = integer : direction flag :    1 = forward
c                                            -1 = inverse
c
c     entry points : prfft only
c
c     common blocks : none
c
c     method : see hockney & jesshope, parallel computers, 302-304
c
c     external references : pcfft
c
c     system utility references : none
c
c     notes : the standard packed format is used for the cfft
c             ie. since the imaginary parts of the zero and
c             nyquist frequency components are zero, real part of
c             the nyquist component is placed in the imaginary part
c             of the zero  frequency component. this allows the
c             transform to occupy the same amount of space as the
c             original data.
c
c     revision history : pc01 - based on vprfft
c                        jh02   add scaling                  11/08/87
c
c
c--------------------------------------------------------------------
c
      subroutine rwsrfft(sig, nreal, flag)

c
c     **** data structure ****
c
      real  pi                   
      parameter  (pi = 3.14159265)
c
c     **** input arguments ****
c
      integer nreal               
      integer flag                
      complex sig(0:nreal/2-1)    
      complex ccc
c
c     **** internal variables ****
c
      integer  ncomp              
      complex  omega              
c                                 
      complex  cc                 
      complex  fplus              
      complex  fminus             
c
c     **** executable code ****
c
      ncomp = nreal / 2           
c
      if(flag.eq.1) then          
c
         call rwscfft(sig(0), ncomp, 1)  
c
ccc      ccc =  exp( -2.0*pi*(0.0,1.0)/real(nreal) )
         do 100 k=1,ncomp/2-1
             omega = (0.0,1.0)*exp(-2.0*pi*(0.0,1.0)*real(k)/
     +               real(nreal))
ccc         omega  =  (0.0 , 1.0)  * ccc
            cc = conjg( sig(ncomp-k) )
            fplus = sig(k) + cc
            fminus = omega * (sig(k) - cc)
            sig(k) = fplus - fminus
            sig(ncomp-k) = conjg(fplus + fminus)
ccc         ccc = ccc*ccc
100      continue
c
         sig(ncomp/2) = 2.0 * conjg( sig(ncomp/2) )
         sig(0) = 2.0 * cmplx( real(sig(0))+aimag(sig(0)),
     +                         real(sig(0))-aimag(sig(0)) )
c
      elseif(flag.eq.-1) then                    
c
ccc      ccc =  exp( +2.0*pi*(0.0,1.0)/real(nreal) )
         do 200 k=1,ncomp/2-1
             omega = (0.0,1.0)*exp(+2.0*pi*(0.0,1.0)*real(k)/
     +               real(nreal))

ccc         omega  =  (0.0 , 1.0)  * ccc
            cc = conjg( sig(ncomp-k) )
            fplus = sig(k) + cc
            fminus = omega * (sig(k) - cc)
            sig(k) = fplus + fminus
            sig(ncomp-k) = conjg(fplus - fminus)
ccc         ccc = ccc*ccc
  200    continue
c
         sig(ncomp/2) = 2.0 * conjg( sig(ncomp/2) )
         sig(0) = cmplx( real(sig(0))+aimag(sig(0)),
     +                   real(sig(0))-aimag(sig(0)) )
c
         call rwscfft(sig(0), ncomp, -1) 
         do 300 i=0,ncomp-1
            sig(i) = sig(i) / 4.0
  300    continue
c
      endif
c
      return
      end


c******************************************************************
c******************************************************************
      subroutine rwscfft(cx,lx,isign)
c----------------------------------------------------------------------
c
c     a subrotine to perform a complex fast fourier transform
c                         lx
c     cx(x) = sqrt(1/lx) sum (cx(j)*exp(2*pi*isign*i*(j-1)*(k-1)/lx) 
c                        j=1            for k=1,2,.....,(lx=2**integer)
c
c     ref. fundamentals of geophysical data processing (claerbout)
c
c----------------------------------------------------------------------
c
c     cx = complex array containing wavelet
c     lx = length of wavelet in samples
c     isign = sign of fft       1=forward fft
c                              -1=inverse fft
c
c----------------------------------------------------------------------
c
c     version id01 written bi i.davis as vpcfft
c             pc02 name changed to pcfft to avoid duplicate vpcfft subr.
c             cj03 vectorise by taking constants from loops
c             pc04 add extra comments
c             jh05 add minus sign in catmp calculation         11/08/87
c
c----------------------------------------------------------------------
      complex cx(*),carg,cexp,cw,ctemp,catmp
*
      integer lx,isign
      integer i,j,l,m,istep
*
      real sc
c
c             sc = scale factor to correct amplitudes for two-way fft
c
      if(isign.eq.-1) then
        sc=1./lx
c
c            scale input array
c
        do 5 i=1,lx
          cx(i)=cx(i)*sc
   5    continue
      endif

c     **  bit reverse order of data
c
      j=1

      do 30 i=1,lx
        if(i.gt.j) go to 10
          ctemp=cx(j)
          cx(j)=cx(i)
          cx(i)=ctemp
  10    m=lx/2
  20    if(j.le.m) goto 25
          j=j-m
          m=m/2
          if(m.ge.1) goto 20
  25        j=j+m
  30  continue

c     ** perform complex fft
c
      l=1
      catmp=-(0.0,1.0)*3.14159265*isign

  40  istep=2*l

      carg=catmp/l
      do 60 m=1,l
        cw=cexp(carg*(m-1))
        do 50 i=m,lx,istep
           ctemp=cw*cx(i+l)
           cx(i+l)=cx(i)-ctemp
           cx(i)=cx(i)+ctemp
  50    continue
  60  continue

      l=istep
      if(l.lt.lx) goto 40
      return
      end
 
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE LJVPWIEN (WAVLT,  LWAVLT, DOUT,LDOUT, 
     +                     WSTAB,  WGAIN , LFIL,  inpanti,
     +                     FILT ,  IOPT  , IER, filmax,
     +                     XC,AC,APWORK,ERFILT,TEMP,nmax)

C-------------------------------------------------------------------
C  THIS SUBROUTINE DESIGNS AND APPLIES A GENERAL PURPOSE
C  WIENER LEAST SQUARES INVERSE FILTER TO SHAPE A GIVEN INPUT
C  (WAVLT) INTO SOME DESIRED OUTPUT (DOUT) OR SIMPLY TO INVERT
C  THE INPUT.
C-------------------------------------------------------------------
C
C  CALLING SEQUENCE :
C
C           CALL VPWIEN( WAVLT, LWAVLT,  IZW,  DOUT,LDOUT, IZD,
C    1                   WSTAB,  WGAIN, LFIL,   IZF,
C    2                    FILT, ACTOUT,  IZA,        IOPT, IER)
C-------------------------------------------------------------------
C  INPUT ARGUMENTS :
C
C  WAVLT    ARRAY CONTAINING THE INPUT WAVELET DIMENSIONED
C           WAVLT(1,..,LWAVLT).
C
C  LWAVLT   LENGTH OF ARRAY WAVLT.
C
C  DOUT     ARRAY CONTAINING THE DESIRED OUTPUT, DIMENSIONED
C           DOUT(1,..,LDOUT).
C             NOTE THAT IF DOUT IS ZERO PADDED AT THE BEGINNING,
C           AN UNNECESSARILY LONG FILTER MAY BE REQUIRED AS THE
C           FILTER MUST PERFORM TRANSLATION AS WELL AS SHAPING.
C           NOT ACCESSED IF IOPT = 1 OR 3. ASSUMES A DELTA FUNCTION.
C
C  LDOUT    LENGTH OF ARRAY DOUT.
C           = 1 FOR IOPT = 1 AND 3.
C
C
C  WSTAB    WHITE LIGHT ADDED TO STABILIZE MATRIX INVERSION IF IT IS
C           NECESSARY AS A PERCENTAGE OF THE ZERO LAG AUTOCORRELATION
C           VALUE.  .OO2 IS SUGGESTED.
C
C  WGAIN    PER CENT WHITE LIGHT ADDED TO ZERO LAG AUTOCORRELATION
C           VALUE TO CONTROL GAIN AT FREQUENCIES WITH LOW SIGNAL TO
C           NOISE RATIO.  SUGGESTED VALUES :
C                          WGAIN = O.5    IF IOPT = 1
C                          WGAIN = 0.     IF IOPT = 2
C
C  LFIL     TOTAL LENGTH OF THE WIENER INVERSE FILTER TO COMPUTE IN
C           SAMPLES.
C
C  INPANTI  THIS IS EQUIVALENT TO THE ANTICIPATION COMPONENT.
C           A GUIDELINE IS: TIME OF PEAK VALUE OF WAVLT - TIME
C           OF PEAK VALUE OF DOUT, IN SAMPLES.  IF NEGATIVE, USE 1.
C           MUST BE LESS THAN LFIL.
C           IANTI=1 CORRESPONDS TO AN ALL MEMORY FILTER.
C
C  IOPT     TYPE OF INVERSE FILTER TO COMPUTE.
C
C           IOPT = 1 FOR AN INVERSE FILTER IN WHICH CASE DOUT SHOULD
C           BE ZERO SAVE FOR DOUT(1) = 1.
C           IOPT = 2 FOR A SHAPING FILTER TO THE DESIRED OUTPUT IN
C           ARRAY DOUT.
C           IOPT = 3 .  AS IOPT = 1 BUT DC IS NOT REMOVED FROM INPUT.
C           IOPT = 4 .  AS IOPT = 2 BUT DC IS NOT REMOVED FROM INPUT.
C-------------------------------------------------------------------
C  OUTPUT ARGUMENTS :
C
C  FILT     ARRAY CONTAINING THE COMPUTED WIENER FILTER, DIMENSIONED
C           FILT(1,..,LFIL).  TIME ZERO CORRESPONDS TO INDEX IANTI
C
C
C  IER      ERROR RETURN FLAG.
C           IER = 1 : ARRAY PROCESSOR SOFT OR HARD FAILURE.
C           IER = 2 : IOPT NOT SET.
C           IER = 5 : IANTI SHOULD BE SUCH THAT 1.LE.IANTI.LE.LFIL
C
      LOGICAL SEARCH

      DIMENSION WAVLT(lwavlt), DOUT(ldout), FILT(lfil)

      REAL XC(NMAX),AC(NMAX),APWORK(NMAX),ERFILT(NMAX),TEMP(NMAX)

 
      DATA ONE /1./

      IF (LFIL.GT.NMAX.OR.LWAVLT+LDOUT.GT.NMAX) THEN
          WRITE(6,*)' THE FILTER IS TOO LONG , OR '
          WRITE(6,*)' THE WAVELETS ARE TOO LONG.'
          ier=-1 
          return 
      END IF




      do i=1,lfil
         filt(i) = 0.0
      enddo

      IER = 0

      ianti = inpanti 


      IF (ianti.eq.0) ianti = 1
      IF (ianti .lt. 0 .OR. ianti .gt. LFIL  ) then    
        IER = 5
        RETURN
      endif


c     ** check for zero input window
c
      do i=1,lwavlt
        if(wavlt(i).ne.0.0) goto 55
      enddo

      ier=66
      return

   55 continue


C       ** DC REMOVAL REQUIRED ?
C
      IF ( IOPT .eq. 2 ) then        
        RLWVLT = LWAVLT
        CALL VPSVE(WAVLT,1,1,SUM,LWAVLT)
        CALL VPDIV(RLWVLT,1,1,SUM,1,1,SUM,1,1, 1)
        CALL VPSUB(SUM,1,0,WAVLT,1,1,WAVLT,1,1,LWAVLT)
      endif



C       ** CLEAR OUT CROSS-CORRELATION VECTOR.
C
      CALL VPCLR(XC,1,1,NMAX)



C       **  SHAPING OPTION. IOPT = 2 OR 4
C
C       PUT CROSS-CORRELATION OF INPUT (WAVLT) WITH DESIRED OUTPUT (DOUT)
C       INTO XC.
C
        IF ( IOPT .EQ. 2 .OR. IOPT .EQ. 4 ) then
           CALL VPCLR(APWORK,1,1,LWAVLT+LFIL)
           CALL VPMOV(DOUT,1,1,APWORK,ianti,1,LDOUT)
           CALL VPCNV(APWORK,1,1,WAVLT,1,1,XC,1,1, LFIL,LWAVLT)

C        ** SPIKING OPTION, IOPT = 1 OR 3
C
C        FILL THE CROSSCORRELATION VECTOR WITH THE FIRST 'IANTI'
C        SAMPLES OF THE INPUT IN REVERSE ORDER, CORRESPONDING TO A
C        TWO-SIDED WIENER INVERSE FILTER DESIGN.
C
        ELSE IF ( IOPT .EQ. 1 .OR. IOPT .EQ. 3 ) then 
           CALL VPMOV(WAVLT,ianti,-1,XC,1,1,ianti)
        else
           IER = 2                               
           RETURN
        endif



C       ** CALCULATE AUTOCORRELATION OF INPUT FOR LAG=0 TO LFIL-1.IN CASE
C          LFIL>LWAVLT WE MUST POST-PAD WAVLT WITH ZEROS TO A TOTAL LENGTH
C          OF LWAVLT+LFIL-1 = LCON SAMPLES.
C
      CALL VPCLR(APWORK,1,1,LFIL+LWAVLT)
      CALL VPMOV(WAVLT,1,1,APWORK,1,1,LWAVLT)
      CALL VPCNV(APWORK,1,1,WAVLT,1,1,AC,1,1,LFIL,LWAVLT)



C       **  ADD WHITE LIGHT TO DIAGONAL FOR GAIN CONTROL.
C
      DIAG = 1. + WGAIN*.01
c     ac(1) = ac(1)*(1. + wgain*0.01

      CALL VPSMUL(AC,1,1,DIAG,AC,1,1,1)



C       **  NOW CALCULATE WIENER INVERSE FILTER. DO NOT SCALE ITS MAXIMUM
C           ABSOLUTE VALUE TO 1.
C
      IFLAG = 1
      CALL VPCLR(TEMP,1,1,NMAX)
      KTRY = 10
      CALL WIEN(AC,TEMP,XC,LFIL,ERFILT,FILT,KTRY,WSTAB,LFIL,1,1)
              
      
      CALL VPMXMG(FILT, 1, 1, ERFILT, LFIL)
      CALL VPDIV(ERFILT,1,1,ONE,1,1,FMAX,1,1,1)
      CALL VPABS(FMAX,1,1, FMAX,1,1, 1)
      CALL VPSMUL(FILT,1,1,FMAX,FILT,1,1,LFIL)

      filmax = fmax
      return
      end









c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine shiftrc(trcin,ipntr,nsmp,itshift,trcout)
c
c     shift in quarter of a sample rate
c
c
c     shift in 1/4 sr
c     ---------------
c
c     itshift = 4*altime*1000/isr - ftbfact*iftbmms*4/isr
c
      dimension trcin(1),trcout(1)
c
c     dimension isshift(20),ipshift(20),twork(20)
      dimension rsinc1(8),rsinc2(8),rsinc3(8)
c
      data rsinc1/ -0.00583580, 0.04035960,-0.14005053, 0.89166689,
     +              0.27481657,-0.07685405, 0.01818759,-0.00057665/
      data rsinc2/ -0.00346141, 0.03929961,-0.14670730, 0.61239052,
     +              0.61239052,-0.14670730, 0.03929961,-0.00346141/
      data rsinc3/ -0.00057665, 0.01818759,-0.07685405, 0.27481735,
     +              0.89166594,-0.14005071, 0.04025969,-0.00583582/
      data nsinc /8/
c
      nsinc2 = nsinc/2
c
      do 2100 i = 1,nsmp
         trcout(i) = 0.0
 2100 continue
c
c
c
c
         if(itshift.ge.0) then
            ismpshf = (itshift+3)/4
            nsmp1 = 1 + ismpshf
            nsmp2 = nsmp
            ij    = ipntr
         else
            ismpshf = itshift/4
            nsmp1 = 1
            nsmp2 = nsmp + ismpshf
            ij = ipntr - ismpshf
         endif
c
         isubsamp = ismpshf*4  - itshift
c
         if(isubsamp.eq.0) then
            do 2200 j = nsmp1,nsmp2
                trcout(j) = trcout(j) + trcin(ij)
                ij = ij + 1
 2200       continue
c
         elseif(isubsamp.eq.1) then
            do 2300 j = nsmp1,nsmp2
               iptr1 = ij - nsinc2
               sincsum = 0.0
               do 2250 k = 1,nsinc
                  sincsum = sincsum + trcin(iptr1+k) * rsinc1(k)
 2250          continue
               trcout(j) = trcout(j) + sincsum
               ij = ij + 1
 2300       continue
         elseif(isubsamp.eq.2) then
            do 2400 j = nsmp1,nsmp2
               iptr1 = ij - nsinc2
               sincsum = 0.0
               do 2350 k = 1,nsinc
                  sincsum = sincsum + trcin(iptr1+k) * rsinc2(k)
 2350          continue
               trcout(j) = trcout(j) + sincsum
               ij = ij + 1
 2400       continue
c
         elseif(isubsamp.eq.3) then
            do 2500 j = nsmp1,nsmp2
               iptr1 = ij - nsinc2
               sincsum = 0.0
               do 2450 k = 1,nsinc
                  sincsum = sincsum + trcin(iptr1+k) * rsinc3(k)
 2450          continue
               trcout(j) = trcout(j) + sincsum
               ij = ij + 1
 2500       continue
         endif
      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine ljshiftrc(trcin,lsamp,sr,shift,trcout)

c
c     Shift trcin(lsamp)  shift [ms] --> trcout(lsamp)
c
c     shift  -  Amount shift in  [ms]  (a positive shift means shift down)
c     sr     -  Sample rate in ms
c
      real trcin(lsamp),trcout(lsamp)

      double precision dshift,dsr

      call set_sinclen('EXACT',16) 


      do 2100 i = 1,lsamp
         trcout(i) = 0.0
 2100 continue
 
      
      dsr    = sr
      dshift = -1.0 * shift



      kmod = INT (dshift/dsr)
      tdif = (dshift - kmod*dsr)/dsr


      if( dshift.le.0.0 ) then
          kmod = kmod - 1
          tdif = 1.0 + tdif
      endif


      do 750 ksmp = 1,lsamp
          ismp = ksmp + kmod
          if(ismp.gt.0.and.ismp.le.lsamp-1) then
            trcout(ksmp) = sinc(trcin,lsamp,ismp,tdif)
          endif
  750 continue
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine single_median(arr,aux,lsamp,nmed,alimit,
     $     vnull_replace,vnull,medbuf)
      real arr(lsamp) , medbuf(nmed) , aux(lsamp)
      character vnull_replace*(*)
c
c     arr(1:lsamp)  ---> Input buffer
c     aux(1:lsamp)  ---> Output buffer  
c     medbuf(1:*)   ---> Work buffer
c
c    

      logical vnl
      logical compute 

      vnull = -999.25

      vnl = vnull_replace.eq.'Y'

      nhits = 0

    
      do i=1,lsamp
        aux(i) = arr(i)    
      enddo
  
      nm2 = nmed/2
      do i=nm2+1,lsamp-nm2
         k=0 

         avsum = 0.0
         aamin = 1e30
         aamax = -1e30
         do j=i-nm2,i+nm2
            if(arr(j).ne.vnull)then
              k=k+1
              medbuf(k) = arr(j)
              if(j.ne.0)avsum = avsum + arr(j)
              aamin = min(aamin,arr(j))
              aamax = max(aamax,arr(j))
            endif
         enddo      

c        ** the average of the other values in the window 
         average = avsum/(k-1)


         compute = k.ge.3


c        ** replace-logic with strong constraints !!!
c        ** should be extremal value
c         if(alimit.ne.1.0) then
c            if(arr(i).eq.aamin.or.arr(i).eq.aamax) then
c               testval = abs( (average-arr(j))/average)
c               compute = testval.gt.alimit
c            endif
c         endif

 
         amed = arr(i)
         if(compute)call tmed(medbuf,1,k,amed)

c         ** simple relative difference
c
         
         testval = abs((amed-arr(i))/average)

         if(testval.gt.alimit) then
            aux(i) = amed
            if(vnl) aux(i) = vnull
            nhits = nhits+1
         endif
         
      enddo

      write(6,*) 'MeDiAn hits:' , nhits
      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine medianfilt(arr,aux,lsamp,nmed,ired,medpct,
     $     vnull_replace , vnull,mode)
      real arr(lsamp) , aux(lsamp)
      character*(*) vnull_replace
      real medpct
c
c     mode          -- 0 = simple median filtering, 1=blocking 
c
c     arr(1:lsamp)  -- Input array     
c     aux(1:lsamp)  -- Output array
c
c
      parameter (maxnmed=121)
      real medbuf(maxnmed)


c     ** Median filtering
      
      logical replace 

      if(nmed.lt.3)return



c    ** just do an interpolation of rindt at the moment 
c
      if(mode.ne.0)call interpolat(arr,lsamp,vnull,1)



      do i=1,lsamp
        aux(i)=arr(i)
      enddo


      if(nmed.gt.maxnmed) then
c       call qmese('Median array too small')
        return
      endif
      

      if(mode.eq.0) then
         apct = medpct/100.
         call single_median(arr,aux,lsamp,nmed,apct,
     $        vnull_replace,vnull,medbuf)
         return
      endif

      llsamp = lsamp



      if(ired.gt.1) then
        k=0
        kksamp = (lsamp+ired-1)/ired
        do i = 1 , kksamp
          sum=0.0
          nsq = 0
          iad = (i-1)*ired
          do j=1,ired           
            if(iad+j.le.lsamp)then 
              sum = sum + arr(iad+j)
              nsq = nsq + 1
            endif
          enddo
         
          if(nsq.gt.0) then 
              k=k+1
              aux(k) = sum/float(nsq)
          endif
        enddo 

        llsamp = k 
        do i=1,llsamp
           arr(i)=aux(i)
        enddo
      endif

      
c     ** do the blocking
c
      call  medblock(arr,aux,medbuf,llsamp,nmed)
      

      if(ired.gt.1) then
        k=0
        do i=1,llsamp
          do j=0,ired-1           
            k=k+1
            if(k.le.lsamp)arr(k) = aux(i)
          enddo
        enddo 
        do i=1,lsamp
          aux(i) = arr(i)
        enddo 
      endif


      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine medblock(arr,aux,medbuf,llsamp,nmed)
      real arr(llsamp) , aux(llsamp)
      real medbuf(llsamp)
c
c
c     Blocking using median filter techniques
c
c
c     arr      -   Input array
c     aux      -   Output array
c     medbuf   -   Work array
c     llsamp    -   Length of input array
c     nmed     -   Length of blocking filter in points
c     ired     -   Reduce the computations involved [ired=1 ==> No reduction] 

      logical replace 

      

      m1 = 3
      m2 = nmed
    
      do imed = m1 , m2 , 2 
        niter = 0   
        do i=1,llsamp
          aux(i)=arr(i)
        enddo

 100    niter=niter+1
        replace = .false.
        nrep = 0
        
        nm2 = imed/2
    
        do 500 i=nm2+1,llsamp-nm2
         k=0 
         do j=i-nm2,i+nm2
            k=k+1
            medbuf(k) = arr(j)
         enddo      


c        ** optimise by comparison with old value
         if(niter.gt.3) then
              nequal=0
              nless =0
              ngt   =0
              do m = 1 , imed
                 if( medbuf(m).eq.aux(i)) nequal=nequal+1
                 if( medbuf(m).lt.aux(i)) nless = nless+1
                 if( medbuf(m).gt.aux(i)) ngt = ngt +1
                 if(nequal.gt.nm2 ) goto 500
              enddo   
            

           if(nless+nequal.gt.nm2.and.nequal.ge.nless) goto 500
           if(ngt  +nequal.gt.nm2.and.nequal.ge.ngt)   goto 500
           if(nless.eq.ngt) goto 500 

         endif

         call tmed(medbuf,1,imed,value)
               
           
         if( value.ne. arr(i) ) then
             nrep=nrep+1
             replace=.true.
         endif

         aux(i) = value

 500   continue
c      ** next i

c       write(6,*) 'mlen:',imed,' Iter:',niter,'  repl:',nrep

       if(replace) then 
           do i=1,llsamp
            arr(i)=aux(i)  
           enddo 
           goto 100
       endif  
       write(6,*) 'median length:',imed,' # Iterations:',niter
      enddo
     
 
      return
      end




c-c-c-c-c-c-c-c-c-c-c--cc--c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
cc    ** simplified median function
c
      real function rmed(array,n)
      real array(n)
      call tmed(array,1,n,amed)
      rmed = amed
      return
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c    Taken from /soft/RWS/source/misc/tmed.f
c
c
c    7>                                                                <72
      Subroutine tmed(x,n,m,y)
c======================================================================c
c
c  Function   : Find median values in m direction of a
c               nxm matrix.
c
c  Operation  : y <--- MEDIAN(x,n,m)
c
c  Parameters :
c          In :
c               x    - input matrix                      (real)
c               n    - number of rows in x               (integer)
c               m    - number of colomns in x            (integer)
c          Out:
c               y    - n dimensional result
c
c  Note       :
c
c  Calls      :  none
c
c  History    :  Date of coding : 26/11/1987
c                Company        : PKR/EIF
c                Modifications  : <none>
c
c
c======================================================================c
c
      real*4       x(1)
c-lj  real*4       y(1)

c
c     Determine if m is odd or even
c
      modd = (m/2)*2-m+1
c
c     For odd number m the median value is found by choosing the
c     first element that has equal number of elements greater/equal
c     and less/equal to itself. If one or more elements are equal
c     the best choice is the element with the least difference between
c     greater/equal and less/equal elements compared.
c     For even number m the procedure is the same but ther will be
c     at least a difference of 1  - this is compensated for in modd.
c
      do 300 i=1,n
         nbest = 999999
         do 200 j=1,m
            ng = 0
            nl = 0
            do 100 l=1,m
               if (x((l-1)*n+i).ge.x((j-1)*n+i)) ng=ng+1
               if (x((l-1)*n+i).le.x((j-1)*n+i)) nl=nl+1
  100       continue
            if ((ng.eq.(nl+modd)).or.(ng.eq.(nl-modd))) then
               nbest = abs(ng-nl)
               nj = j
               go to 300
            end if
            if (abs(ng-nl).lt.nbest) then
               nbest = abs(ng-nl)
               nj = j
            end if
  200    continue
c 300 y(i) = x((nj-1)*n+i)
  300 y    = x((nj-1)*n+i)
      Return
      End



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
c
      subroutine maxcorr(a1 , la1 , a2 , la2 , idelay , subsamp ,cmax)
c
c      correlate sequences a1 and a2 , Return sample for best correlation
c      The function argument is max correlation value
c      The sequence a1 is shorter than a2 .only  positive lags considered
c


       integer la1,la2 
       integer idelay 
       real cmax

       real a1(la1),a2(la2)


c      Local variables
c
       double precision  aapos(-1:1) , sump ,sum1,sum2,aver,psum
       double precision  rollon , rollof , avern
       double precision  qa , qb , qc , xpeak


       ll1 = la1  
       ll2 = la2

       aapos(-1)  =  0.0
       aapos(0)   = -1.0
       aapos(1)   =  0.0 
       avern      =  0.0


       idelpos = 0


       nlags = la2 - la1     

       if(nlags.lt.0) then
         cmax = 0.0
         idelay = - 1
         return
       endif 



c      ** Find energy in channel 1
c      ** and channel 2
       sum1=0.0
       sum2=0.0
       do 100 i = 1,la1  
           sum1 = sum1 + a1(i) * a1(i)
           sum2 = sum2 + a2(i) * a2(i)
  100  continue                                  



       do ilag = 0 , nlags  

         sump=0.0
         do i=1,la1  
           sump = sump + a1(i) * a2 (i+ilag)
         end do   

         if(ilag.gt.0) then
            rollof = a2(ilag)       * a2(ilag)
            rollon = a2(ilag + la1) * a2(ilag + la1) 
            sum2 = sum2 - rollof + rollon 
         endif

         psum = sum1 * sum2 
         if(psum.gt.0.0) then
            psum = dsqrt(psum)
            aver = sump / psum 
            if(aver.gt.aapos(0) ) then
               aapos(0)   = aver
               idelpos = ilag
               aapos(-1) = avern
            endif  
            avern = aver
            if(ilag.eq.idelpos+1) aapos(1) = aver
         endif
       end do   


c      ** find subsample (second order polynomial)
c

       if(idelpos.gt.0.and.idelpos.lt.nlags) then

         qa = 0.5*(aapos(-1)+ aapos(1) - 2*aapos(0))
         qb = 0.5*(aapos(1) - aapos(-1)) 
         qc = aapos(0)

c        **  x= -b/(2a)
c
         xpeak  = 0.5*( aapos(1)  - aapos (-1)  )/
     +          (2*aapos(0) - aapos(1) - aapos(-1))

         cmax = qa*xpeak*xpeak + qb*xpeak + qc
         
       else 
         cmax    = aapos(0)
         xpeak   = 0.0
       endif
       subsamp = xpeak
       idelay = idelpos
       return
       end 
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
c
      subroutine ccorr_stack( rtrace , lsamp , ntrace  ,
     +                        stack,medbuf,stktyp,iftrace,iltrace,
     +                        cdelay , ccorr , flag  ,
     +                        maxdel , maxiter , rlimit , niter)

      SAVE
      character*(*)   stktyp           ! Straight / median stack
      real     rtrace(1)               ! Data buffer 
      integer  lsamp                   ! Nr of samples
      integer  ntrace                  ! Nr. of traces
      real     stack (lsamp)           ! stack buffer        
      real     cdelay (ntrace)         ! The individual delays
      real     ccorr (ntrace)          ! Cross corr coefficients
      integer  flag  (ntrace)          ! Live trace flags 
      real     medbuf(ntrace)          ! Median buffer
      integer  maxdel                  ! Max delay to try
      integer  maxiter                 ! Max nr. of iterations
      real     rlimit                  ! Stop iterating
      integer  niter                   ! Return the number of iterations

c     Local variables 
c

      integer iter      ! Iteration count
      integer lwin      ! Length of trace window
      integer iad_data  ! Address in data buffer 
      integer iret      ! Return value from check routine 
      integer idelay    ! Last delay in samples




      real    cmax      ! Last corr. coef
      real    sumcorr
      real    oldave
      real    avecorr

c
c     The purpose of this routine is twofold 
c
c     1) Produce a crosscorrelation stack
c     2) Compute the delays of the individual traces
c

      nsamp = lsamp




c     ** The correlating windows of the traces will have this length
c

      lwin=lsamp-2*maxdel

      iter = 0
      avecorr = 0.0



c     ** iteration start point 
  100 continue
      iter = iter + 1
      oldave = avecorr 

      if (iter.gt.maxiter) goto 9999


c     ** compute straight or median stack   
c     ** 

      do 102 i = 1 , lsamp

          stack(i) = 0.0
          kfold    = 0

          do 101 itrace = iftrace , iltrace 

            if(flag(itrace) .eq.1) then

              iad_data = (itrace-1)*lsamp
              icdel =  INT(cdelay(itrace))
              j = i -  icdel
              tdif = abs ( cdelay(itrace) - icdel) 

              if( cdelay(itrace).le.0.0 ) then
                j=j-1
                tdif = 1.0 - tdif
              endif

              if(j.eq.0) then
                 j=1
                 tdif=0.0
              elseif(j.eq.lsamp) then
                 j=lsamp-1
                 tdif=1.0
              endif


              if(j.ge.1.and.j.le.lsamp) then
c               val =  tabsinc(rtrace(iad_data+1),lsamp,j,tdif)
                val =     sinc(rtrace(iad_data+1),lsamp,j,tdif)
                kfold = kfold + 1
                if(stktyp.eq.'STRAIGHT') then 
                  stack(i)=stack(i) + val
                else
                  medbuf(kfold) = val
                endif
             endif

            endif

  101     continue          

c - lj    3/11.99
          if(kfold.ge.1)then
          if(stktyp.eq.'STRAIGHT') then 
            stack(i) = stack(i)/float(kfold)  
          else
            call tmed (medbuf,1,kfold, stack(i) )
          endif
          endif
  102 continue          






c     ** correlate all traces with the stack 
c
      sumcorr = 0.0

      kfold = 0
      do 911 itrace = 1 , ntrace 

          if(flag(itrace) .eq. 1) then


          iad_data = (itrace-1)*lsamp + maxdel +1

c         ** Find delay idelay for maximum correlation
c

          call maxcorr(rtrace(iad_data),lwin,stack,lsamp,idelay,ss,cmax)


c           ** delay in samples
          if(idelay.eq.-1) then
              write(6,*) ' ccorr_stack : LOGICAL ERROR'
              idelay = maxdel
          elseif(idelay.eq.0.or.idelay.eq.2*maxdel) then
c             write(6,*) ' ccorr_stack : No peak in correlation'
              idelay = maxdel
          endif

          kfold = kfold + 1
          cdelay(itrace) = idelay   -  maxdel + ss
          ccorr(itrace) = cmax  

          sumcorr = sumcorr + ccorr(itrace)

          else
            cdelay(itrace) = 0.0
            ccorr(itrace)  = 0.0
          endif

  911 continue                           

c      ** the average correlation for this trace
c
       avecorr = 0.0
       if(kfold.ge.1) avecorr = sumcorr / float(kfold)

       if(iter.gt.1 ) then
         if(avecorr.lt.oldave) then
ccc        write(6,*) ' ccorr_stack : correlation decreasing'
           goto 9999
         elseif(avecorr-oldave.lt.avecorr*rlimit) then
           goto 9999
         endif
       endif
       goto 100

 9999 continue
      niter = min(iter,maxiter)

      return
      end





c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine frconv(c1,c2,trace,ctrace,ksamp,ishift,sr)

c     ** Convolution by multiplication in frequency domain
c

      complex c1(ksamp/2)  , c2(ksamp/2)
      complex ctrace(ksamp/2)              ! NB ,same adress as trace 

      real trace(ksamp)

      pi = atan (1.0) * 4.0

      nw = ksamp/2
      df      = 500.0 / (sr*nw)
      freq    = 0.0
      ts = -1.0*ishift*sr/1000.0

      do iw = 1, nw
           w = 2.0 * pi * freq
           ctrace (iw) = c1(iw)*c2(iw) * cexp( (0.0,-1.0) * w * ts )
           freq = freq + df
      enddo

      call rwsrfft(trace,ksamp,-1)
 
      return
      end




c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  ctrconv(aux,lsamp,filt,nfilt,arr1,narr1,ishift)
      complex aux(lsamp) , filt(nfilt) , arr1(narr1)

c
c     Convolution of a trace and a filter 
c     Input and output has the same length
c     Aux is the input trace
c
c     Replacement installed 26 OCT 95 by lj

      complex*16 tem

      do i=1,narr1
         arr1(1) = cmplx( 0.0 , 0.0)
      enddo

     
      do k =  1 ,   max(lsamp,narr1)

            ipoint = k + ishift
            iback = max (1,ipoint-nfilt+1)
            ito   = min (lsamp,ipoint)

            tem = cmplx( 0. , 0.)
            do i = iback,ito
               if(ipoint-i+1.gt.nfilt) stop 'ctrconv'
               tem = tem + aux(i) * filt(ipoint-i+1)
            enddo        

            arr1(k) = tem 
      enddo

      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  trconv(aux,lsamp,filt,nfilt,arr1,ishift)
      real aux(lsamp) , filt(nfilt) , arr1(lsamp)
c
c     Convolution of a trace and a filter 
c     Input and output has the same length
c     Aux is the input trace
c

c     Replacement installed 26 OCT 95 by lj

     
      do k =  1 ,   lsamp 

            ipoint = k + ishift
            iback = max (1,ipoint-nfilt+1)
            ito   = min (lsamp,ipoint)

            tem=0.
            do i = iback,ito
               tem = tem + aux(i) * filt(ipoint-i+1)
            enddo        

            arr1(k) = tem 
      enddo

      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPFILR(IPHASE,SAMP,B,FILT,NFILT,maxfilt,ctemp,mctemp)
C
C----------------------------------------------------------------------
C        ROUTINE TO DESIGN A RICKER WAVELET
C----------------------------------------------------------------------
C
      real  FILT(maxfilt)

      parameter (np=2048)

      complex ctemp(mctemp)

      LOGICAL    GTMAX

      DATA   PI / 3.1415927 /
      DATA   GTMAX / .FALSE. /
      data icall/0/
      icall=icall+1
c     if(icall.eq.1) ip = malloc(np*8)


C
c     TERMR=5
c     TERMW=6
      NLEN = NINT(B/SAMP)
      IF(NLEN.GT.256) THEN
         NSAMPS = 256
         GTMAX = .TRUE.
      ELSE
         NSAMPS = NLEN
      ENDIF

c-lj  IF(IPHASE .EQ. 2)THEN
c-lj    NSAMPS=1024
c-lj  ENDIF
      NFILT = 2*NSAMPS
C         CHANGED FROM 2*NFILT+1 28/10/86
C
      TIME = 0.0
      DO 10 I=NSAMPS+1, 2*NSAMPS+1
         FT = 6.0 * (TIME**2.0) / (B**2.0)
         FILT(I) = SQRT( PI ) / 2.0 * ( FT-0.5 ) * EXP( -1.0*FT )
         FILT(I) = -4.0 / SQRT( PI ) * FILT(I)
         J = 2*NSAMPS + 2 - I
         FILT(J) = FILT(I)
         TIME = TIME + SAMP
   10 CONTINUE
C
C          MINIMUM PHASE OPTION
C
      IF (IPHASE.EQ.2) THEN
         IF (IREMP2(NFILT).NE.0)THEN
           NFILT=NFILT+IREMP2(NFILT)
         ENDIF

c-lj
           nfilt = 2048
c-lj
           DO 7711 II=1,NFILT          
 7711        CTEMP(II)=CMPLX(FILT(II),0.0)
           CALL MINPHS(CTEMP,NFILT)
           NFILT=NLEN * 4
           DO 7713 II=1,NFILT
 7713         FILT(II)=REAL(CTEMP(II))
           CALL VPTAPR(FILT,1,NFILT,0,3,IERR)
      END IF
      RETURN
      END



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine vpfilp(itype,iphase,samp,sl,hppl,sh,hpph,
     $     filt,maxfilt,ctemp,mctemp)

      real filt(maxfilt) 
c
c
c----------------------------------------------------------------------
c        routine to design a butterworth low-pass, high-pass or
c        band-pass filter
c
c        itype=1  ==> Bandpass
c        itype=2  ==> Highpass
c        itype=3  ==> Lowpass
c
c        iphase=1  ==> Zero phase
c        iphase=2  ==> Minimum phase
c
c        samp     - Sample rate in second
c
      parameter (np=4096)
      complex ctemp(mctemp)

c-lj  integer    f, fnyq
       real      f, fnyq

c-lj
      double precision dp1,dp2,one,a,pi

      real       n
      integer    m

      data icall/0/
 

c-lj  data  pi / 3.1415927 /
      one=1.0
      pi = datan(one)*4.0  
c-lj end 
 
      nsamps=256
      fnyq = 1.0 / (2.0 * samp)


c-lj  
      finc = 1.0

c-lj       
c     ** fnyq should be an integer according to algo. design
      if(fnyq.lt.124) finc = fnyq/125.



      n = sh / (20.0 * log10(2.0))
      m = sl / (20.0 * log10(2.0))




      if (iphase.eq.2) then
        icall=icall+1
c       if(icall.eq.1) ip = malloc(np*8)
        call vpclr(filt, 1, 1, 512)
      else
        call vpclr(filt, 1, 1, 257)
      end if



      if(itype.eq.1) then
c-lj     do 20 f=1,fnyq
c-lj           a = 1.0 / ( sqrt(1.0 + (f/hpph)**(2.0*n)) *
c-lj +                     sqrt(1.0 + (hppl/f)**(2.0*m)) )

         do 20 f=finc,fnyq+0.1,finc

               dp1 = one*f/hpph
               dp2 = one*hppl/f

               a = 1.0 / ( dsqrt(one + dp1**(2.0*n)) *
     +                     dsqrt(one + dp2**(2.0*m)) )

            do 10 i=129,257
               filt(i) = filt(i) + a*cos(2.0 * pi * f * (i-129) * samp)
   10       continue

   20    continue
 
      elseif(itype.eq.2) then

         do 40 f=finc,fnyq+0.1,finc
            do 30 i=129,257
               a = 1.0 /  sqrt(1.0 + (hppl/f)**(2.0*m))
               filt(i) = filt(i) + a*cos(2.0 * pi * f * (i-129) * samp)
   30       continue
   40    continue
 
      elseif(itype.eq.3) then
         do 60 f=finc,fnyq+0.1,finc
            do 50 i=129,257
               a = 1.0 / sqrt(1.0 + (f/hpph)**(2.0*n))
               filt(i) = filt(i) + a*cos(2.0 * pi * f * (i-129) * samp)
   50       continue
   60    continue
      endif
 
      ilast=257
 
c     ** scale maximum value to 1
c
      call vpmax(filt, 257, ampmax, idum)
      do 100 i=129,257
         filt(i) = filt(i) / ampmax
         filt(258-i) = filt(i)
  100 continue
      nfilt=257
 
c      **   minimum phase option
c
      if (iphase.eq.2) then
         if (iremp2(nfilt).ne.0)then
           nfilt=nfilt+iremp2(nfilt)

           do 8811 ii=1,nfilt          
 8811        ctemp(ii)=cmplx(filt(ii),0.0)

           call minphs(ctemp,nfilt)

           do 8813 ii=1,nfilt
 8813         filt(ii)=real(ctemp(ii))
         endif
      endif
 
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine minphs(cx,lx)
c-----------------------------------------------------------------------------
c
c           a subroutine to convert a complex wavelet into its minimum
c           phase equivalent by use of a hilbert transform technique
c
c           ref. fundamentals of geophysical data processing (claerbout)
c
c------------------------------------------------------------------------------
c
c           cx = complex array containing input wavelet and returning 
c                min. phs. wavelet
c           lx = no. of elements in complex array (must be a power  of 2)
c
c------------------------------------------------------------------------------
      complex cx(*)
      integer lx
c
      k=lx/2
      l=k/2
c
c     forward fft
c
      call rwscfft(cx,lx,1)
c
c     check for complex zero
c
      do  5 i=1,lx
        if(cx(i) .eq. cmplx(0.0,0.0))then
          cx(i) = cmplx(.00001,0.0)
        endif
    5 continue
c
c     form log of spectrum
c
      do 10 i=1,lx
  10  cx(i)=.5*clog(cx(i)*conjg(cx(i)))
c
c     begin hilbert transform
c
      call rwscfft(cx,lx,-1)
c
c     leave t=0, double values at positive time, zero neg. time values
c
      do 20 j=2,k
      cx(j)=cx(j)+cx(j)
  20  cx(lx+2-j)=cmplx(0.)
c
c     end hilbert transform
c
      call rwscfft(cx,lx,1)
c
c     exponentiate
c
      do 30 i=1,lx
  30  cx(i)=cexp(cx(i))
c
c     inverse fft
c
      call rwscfft(cx,lx,-1)
      return
      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      function iremp2(iarg)
c-------------------------------------------------------------
c                  return the remainder of the argument
c                  when the largest power of 2 that is
c                  greater than or equal to the argument
c                  is subtracted
c
c-------------------------------------------------------------
      ipow=1
 500  if(ipow.ge.iarg.or.ipow.gt.10000)goto 1000 
	 ipow=ipow*2
      goto 500
 1000 iremp2= ipow-iarg
      return

      end
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine av3smooth(arr,lsamp)
      real arr(lsamp)

      nsm2 = 1
      do i = 2 ,lsamp
         j1=max(i-nsm2,1)
         j2=min(i+nsm2,lsamp)
         if (arr(i).gt.arr(j1).and.arr(i).gt.arr(j2) ) then
           arr(i) = (arr(j2) + arr(j1))/2.0
         endif
      enddo
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c       3d velocity smoothing
c
c
	subroutine  vsmooth(velbuf,lsamp,ntrace,nsx,nsz,work,
     &                      type,vnull)
        character*(*) type

        real work(*)

	real velbuf(lsamp,ntrace) 


c      ** horizontal direction       
c
        if(nsx.gt.2)then
           do isamp = 1,lsamp
              do itrace = 1  , ntrace
                 work(itrace) = velbuf(isamp,itrace)
              enddo
              if(type.eq.'AVERAGE') then
                 call avsmooth(work, work(ntrace+1),ntrace,nsx,vnull)
              else
                 call medsmooth(work,work(ntrace+1),
     &           work(2*ntrace+1),ntrace,nsx,vnull)
              endif

              do itrace = 1  , ntrace
                 velbuf(isamp,itrace) = work(ntrace+itrace)
              enddo
           enddo
        endif


c      ** vertical direction       
c
        if(nsz.gt.2)then
           do itrace = 1,ntrace
              do isamp = 1  , lsamp
                 work(isamp) = velbuf(isamp,itrace)
              enddo
              if(type.eq.'AVERAGE') then
                 call avsmooth(work, work(lsamp+1),lsamp,nsz,vnull)
              else
                 call medsmooth(work,work(lsamp+1),
     &           work(2*lsamp+1),lsamp,nsz,vnull)
              endif
                 

              do isamp = 1  , lsamp
                 velbuf(isamp,itrace) = work(lsamp+isamp)
              enddo
           enddo
        endif
        return
        end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
c      subroutine avsmooth(arr,oarr,lsamp,nsmooth)
c      real arr(lsamp) ,oarr(lsamp)
c      nsm2 = nsmooth/2
c      do i = 1 ,lsamp
c         j1=max(i-nsm2,1)
c         j2=min(i+nsm2,lsamp)
c         sum=0
c         nv=0
c         do j=j1,j2
c           nv = nv+1
c           sum=sum+arr(j)
c         enddo
c         oarr(i) = sum/nv
c      enddo
c      return
c      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c  from calib
c
      subroutine avsmooth(arr,oarr,lsamp,nsmooth,vnull)
      real arr(lsamp) ,oarr(lsamp)

      nsm2 = nsmooth/2

      do i = 1 ,lsamp
         oarr(i) = vnull
         if(arr(i).ne.vnull) then 
            j1=max(i-nsm2,1)
            j2=min(i+nsm2,lsamp)
            sum=0
            nv=0
            do j=j1,j2
               if( arr(j).ne. vnull) then
                  nv = nv+1
                  sum = sum + arr(j)
               endif
            enddo
            oarr(i) = sum/nv
         endif
      enddo
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine medsmooth(arr,oarr,work,lsamp,nsmooth,vnull)
      real arr(lsamp) ,oarr(lsamp)
      real work(nsmooth)


      if(nsmooth.gt.lsamp) then
         ic = 0
         do k = 1,lsamp
            if( arr(k).ne.vnull) then
               ic = ic + 1
               work(ic) = arr(k)
            endif
         enddo
         call tmed(work,1,ic,amed)
c         print *,' amed= ',amed
         do i=1,lsamp
            oarr(i) = amed
         enddo
         R E T U R N
      endif



      nsm2 = nsmooth/2
      do i = 1 ,lsamp

         j1=max(i-nsm2,1)
         j2=min(i+nsm2,lsamp)

         np = j2 - j1 + 1

         ic = 0
         do k = j1,j2
            if( arr(k).ne.vnull) then
               ic = ic + 1
               work(ic) = arr(k)
            endif
         enddo
         if(ic.gt.0) then
            call tmed(work,1,ic,oarr(i))
         else
            oarr(i) = vnull
         endif
      enddo
      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine get_extrm(arr,lsamp,rindt,ifsamp,ilsamp,nrindt)
      real arr(lsamp)
c
c     Get first and last samplke not equal to rindt 
c

      ifsamp=0
      ilsamp=-1

      do i=1,lsamp
        if (arr(i).ne.rindt) then
          ifsamp=i
          goto 10
        endif
      enddo
      return

   10 continue

c     ** Find last  sample
c

      do i=lsamp,1,-1
        if (arr(i).ne.rindt) then
          ilsamp=i
          goto 11
        endif
      enddo

   11 continue


      nrindt=0
      do i=ifsamp,ilsamp
        if (arr(i).eq.rindt) then
          nrindt=nrindt+1
        endif
      enddo


      return
      end 

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine interpolat(arr,lsamp,rindt,emode)
      real arr(lsamp)
      integer emode  
c
c     emode = 0  No extrapolation
c     emode = 1  Linear extrapolation
c     emode = 2  Constant extrapolation
c
      double precision val,delta,c,x1,x2,y1,y2

      logical copymode

c     ** Interpolate values not equal to rindt
c     ** Find first sample
c     ** extrapolate later


      call get_extrm(arr,lsamp,rindt,ifsamp,ilsamp,nrindt)
      if(ifsamp.le.0.or.ilsamp.le.0)return
 
      y1 = arr(ifsamp)
      x1 = ifsamp
      x2 = 0.0

      do i = ifsamp , ilsamp

        if(arr(i) .eq. rindt) then
           if(copymode) then
             ix1 = i-1
              y1 = arr(i-1)

              do 20 j= i+1,lsamp
               if (arr(j).eq.rindt) goto 20
               ix2 = j
               y2 = arr(j)
                c = (y2 - y1)/float(ix2 - ix1)
               goto 21
   20         continue 
   21         copymode=.false.
           endif
           arr(i) = y1 + c*(i-ix1)
        else
           copymode=.true.
        endif

      enddo

c     **  extrapolation
c     
      if(ifsamp.gt.1.and.emode.ne.0) then

         if(ifsamp.lt.lsamp) then
            delta = arr(ifsamp) - arr(ifsamp+1)
            if(arr(ifsamp+1).eq.rindt) delta = 0.0
         else
            delta = 0.0
         endif

         if(emode.eq.2) delta=0.0

         val = arr(ifsamp) 
         do i= ifsamp-1,1,-1
            val = val + delta
            arr(i) = val
         enddo
      endif

      if(ilsamp.lt.lsamp.and.emode.ne.0) then
         if(ilsamp.gt.1) then
            delta = arr(ilsamp) - arr(ilsamp-1)
            if(arr(ilsamp-1).eq.rindt) delta = 0.0
         else
            delta = 0.0
         endif
         if(emode.eq.2)delta=0.0 
         val = arr(ilsamp)
         do i= ilsamp+1,lsamp
            val = val + delta
            arr(i) = val
         enddo
      endif

      return
      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine interpolate_old(arr,lsamp,rindt)
      real arr(lsamp)
      
      double precision val,delta

      logical copymode

c     ** Interpolate values not equal to rindt
c     ** Find first sample
c     ** extrapolate later

      ifsamp=0
      do i=1,lsamp
        if (arr(i).ne.rindt) then
          ifsamp=i
          goto 10
        endif
      enddo
      return

   10 continue

c     ** Find last  sample
c
      ilsamp=-1
      do i=lsamp,1,-1
        if (arr(i).ne.rindt) then
          ilsamp=i
          goto 11
        endif
      enddo
      return
   11 continue


      y1 = arr(ifsamp)
      x1 = ifsamp
      x2 = 0.0

      do i = ifsamp , ilsamp

        if(arr(i) .eq. rindt) then
           if(copymode) then
             ix1 = i-1
              y1 = arr(i-1)

              do 20 j= i+1,lsamp
               if (arr(j).eq.rindt) goto 20
               ix2 = j
               y2 = arr(j)
                c = (y2 - y1)/float(ix2 - ix1)
               goto 21
   20         continue 
   21         copymode=.false.
           endif
           arr(i) = y1 + c*(i-ix1)
        else
           copymode=.true.
        endif

      enddo

c  **  extrapolation
c
      if(ifsamp.gt.1) then
       if(ifsamp.lt.lsamp) then
         delta = arr(ifsamp) - arr(ifsamp+1)
         if(arr(ifsamp+1).eq.rindt) delta = 0.0
       else
         delta = 0.0
       endif

        val = arr(ifsamp) 
        do i= ifsamp-1,1,-1
         val = val + delta
         arr(i) = val
        enddo
      endif

      if(ilsamp.lt.lsamp) then
        if(ilsamp.gt.1) then
          delta = arr(ilsamp) - arr(ilsamp-1)
          if(arr(ilsamp-1).eq.rindt) delta = 0.0
        else
          delta = 0.0
        endif
        val = arr(ilsamp)
        do i= ilsamp+1,lsamp
         val = val + delta
         arr(i) = val
        enddo
      endif

      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
      subroutine intpol(values,nl,val,r)
c
c
c     Return fractional sample  
c
c     the values array must be monotonous [rising]
c
      real values(nl)

      r = 0.0
      do i=1,nl-1
         if(values(i).lt.val.and.values(i+1).gt.val) then
            diff = values(i+1) - values(i)
            r = i + (val  - values(i)) / diff
            return
         endif
      enddo
      return
      end





c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c     ** given a monotonous function varr(n) , return index of value v
c
      real function arrfind(v,varr,n)
      real varr(n)

      arrfind = 0.0

      if(varr(1).lt.varr(n))then
       if(v.lt.varr(1)) return
       
       do i=1,n-1
         if(v.ge.varr(i).and.v.lt.varr(i+1))then
            arrfind = i + (v-varr(i))/(varr(i+1)-varr(i))
            return
         endif
       enddo
      
       else  

       if(v.lt.varr(n))then
          arrfind = float(n)
          return
       endif
       
       do i=1,n-1
         if(v.le.varr(i).and.v.gt.varr(i+1))then
            arrfind = i + (v-varr(i+1))/(varr(i)-varr(i+1))
            return
         endif
       enddo
      endif

      return
      END

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine i2pol(cnt,arr,n,vcnt, varr , iret )
      real cnt(n) , arr(n) , vcnt , varr
c
c     ** Interpolate values in arr              
c     ** Use the value vcnt (element in cnt) as control
c
c     vcnt (I)        -   sample value , x-coordinate corresponding to 
c                                  output value varr 
c     cnt (1:n)  (I)   -   sample array (x coordinates)
c     arr (1:n)  (I)   -   Value array  (y coordinates)
c     varr             -   Output value at x-coordinate vcnt
c
      iret=-1 
      do 99 i=1,n-1
         if(vcnt .ge. cnt(i) .and.
     +      vcnt .le. cnt(i+1) ) then
            it1=i
            it2=i+1
            diff = cnt(it2) - cnt(it1)
            gscale = (vcnt  - cnt(it1)) / diff
            varr = arr(it1)*(1.-gscale) + arr(it2)*gscale
            goto 199
         endif
   99 continue
c     call qerrmsg('i2pol outside limits')
      return
  199 iret=0
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine rws_taper(arr1,np,itaper,aux,iopt) 
      real arr1(np),aux(np)
c
c     itaper - taper length in samples 
c
c     iopt = 0  taper on
c     iopt = 1  taper off
c     iopt = 2  on and off 

      pi=atan(1.0)*4.0

      do 778 i=1,np
          aux(i)=arr1(i)
  778 continue

      if(itaper.le.0) return


c     ** Taper on
c

      if (iopt.eq.0.or.iopt.eq.2) then
        ifsamp=1
        ilsamp=itaper
        npp= ilsamp-ifsamp+1
        dv= pi/npp
        v=pi
        do 777 i=ifsamp , ilsamp 
           v = v-dv
           weight = 0.5 + 0.5*cos(v)
           aux (i) = arr1(i) * weight           
  777   continue 
      endif


c     ** Taper off
c
      if(iopt.eq.1.or.iopt.eq.2) then
        ifsamp = np-itaper+1
        ilsamp = np
        npp= ilsamp-ifsamp+1
        dv= pi/npp
        v=0.0
        do 877 i = ifsamp , ilsamp
           v = v+dv
           weight = 0.5 + 0.5*cos(v)
           aux (i) = arr1(i) *  weight
  877   continue 
      endif
      return
      end











      SUBROUTINE VPTAPR(ARRAY,ISTART,IEND,IUP,ITYPE,IERR)
C-----------------------------------------------------------------------
C
C     A SUBROUTINE TO APPLY A TAPER TO A DATA ARRAY OVER A SPECIFIED
C     RANGE OF SAMPLES.
C
C-----------------------------------------------------------------------
C
C             ARRAY  = ARRAY OF REAL VALUES TO BE TAPERED
C             ISTART = START SAMPLE (ARRAY INDEX) TO BE TAPERED
C             IEND   =   END SAMPLE (ARRAY INDEX) TO BE TAPERED
C             IUP    = FLAG : 0 = TAPER DOWN, 1 = TAPER UP
C                              2 = DOUBLE SIDED TAPER (UP THEN DOWN)
C             ITYPE  = TAPER TYPE : 0 = LINEAR
C                                   1 = COSINE
C                                   2 = COSINE SQUARED (HANNING)
C                                   3 = HAMMING
C                                   4 = COSINE ** 0.5
C             IERR    = ERROR FLAG: 0 = O.K.
C                                   1 = END BEFORE START
C                                   2 = IUP OUTSIDE RANGE PERMITTED
C                                   3 = ITYPE NOT RECOGNISED
C
C
C-----------------------------------------------------------------------
C
C             VERSION: PC01  DATE 24/01/85
C             UPDATE :
C
C-----------------------------------------------------------------------
C
      REAL ARRAY(*)
C
      INTEGER ISTART,IEND,IUP,ITYPE,IERR
      INTEGER TERMW
C
      PARAMETER (PI=3.1415926)
C
      DATA  TERMW  / 6 /
C---------------------------------------------------------------------
C        FORMATS
C---------------------------------------------------------------------
 100     FORMAT(/11X,A,$)
C-----------------------------------------------------------------------
C             SET DEFAULTS
C-----------------------------------------------------------------------
      IERR = 0
      HPI = PI/2.0
C-----------------------------------------------------------------------
C             ERROR CHECKS
C-----------------------------------------------------------------------
      IF(ISTART .GT. IEND)THEN
        WRITE(TERMW, 100)'*** TAPER WINDOW ERROR ***'
        IERR = 1
        GOTO 999
      ELSEIF(IUP .LT. 0 .OR. IUP .GT. 2)THEN
        WRITE(TERMW, 100)'*** TAPER DEFINITION ERROR ***'
        IERR = 2
        GOTO 999
      ELSEIF(ITYPE .LT. 0 .OR. ITYPE .GT. 4)THEN
        WRITE(TERMW, 100)'*** TAPER TYPE UNKNOWN ***'
        IERR = 3
        GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C             SET UP SCALE FACTOR INCREMENT
C-----------------------------------------------------------------------
      NINC = IEND - ISTART
      NSAMP = NINC + 1
C-----------------------------------------------------------------------
C             TYPE = 0, LINEAR
C-----------------------------------------------------------------------
      IF(ITYPE .EQ. 0)THEN
        SCRINC = 1.0/NINC
        IF(IUP .EQ. 0)THEN
          DO 2000 ISAMP=ISTART,IEND
 2000       ARRAY(ISAMP) = ARRAY(ISAMP) * (1.0 - SCRINC*(ISAMP-ISTART))
        ELSEIF(IUP .EQ. 1)THEN
          DO 2010 ISAMP=ISTART,IEND
 2010       ARRAY(ISAMP) = ARRAY(ISAMP) * (0.0 + SCRINC*(ISAMP-ISTART))
        ELSEIF(IUP .EQ. 2)THEN
          NINC = (NSAMP+1)/2 - 1
          SCRINC = 1.0/NINC
          IMID   = (NSAMP+1)/2 + ISTART - 1
          DO 2020 ISAMP=ISTART,IMID
 2020       ARRAY(ISAMP) = ARRAY(ISAMP) * (0.0 + SCRINC*(ISAMP-ISTART))
          IF(MOD(IMID,2) .EQ. 0)THEN
            IMID = IMID + 1
          ENDIF
          DO 2021 ISAMP=IMID,IEND
 2021       ARRAY(ISAMP) = ARRAY(ISAMP) * (1.0 - SCRINC*(ISAMP-IMID))
        ENDIF
C-----------------------------------------------------------------------
C             TYPE = 1, COSINE
C----------------------------------------------------------------------
      ELSEIF(ITYPE .EQ. 1)THEN
        ANGINC = PI/(2.0*NINC)
        IF(IUP .EQ. 0)THEN
          DO 2100 ISAMP=ISTART,IEND
 2100       ARRAY(ISAMP) = ARRAY(ISAMP) * COS(0.0+ANGINC*(ISAMP-ISTART))
        ELSEIF(IUP .EQ. 1)THEN
          DO 2110 ISAMP=ISTART,IEND
 2110       ARRAY(ISAMP) = ARRAY(ISAMP) * COS(HPI-ANGINC*(ISAMP-ISTART))
        ELSEIF(IUP .EQ. 2)THEN
          ANGINC =   PI / NINC
          DO 2120 ISAMP=ISTART,IEND
 2120       ARRAY(ISAMP) = ARRAY(ISAMP)*COS(HPI-ANGINC*(ISAMP-ISTART))
        ENDIF
C-----------------------------------------------------------------------
C             TYPE = 2, COS ** 2.0
C-----------------------------------------------------------------------
      ELSEIF(ITYPE .EQ. 2)THEN
        ANGINC = PI/(2.0*NINC)
        IF(IUP .EQ. 0)THEN
          DO 2200 ISAMP=ISTART,IEND
 2200       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &                    (COS(0.0 + ANGINC*(ISAMP-ISTART))) ** 2.0
        ELSEIF(IUP .EQ. 1)THEN
          DO 2210 ISAMP=ISTART,IEND
 2210       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &                    (COS(HPI - ANGINC*(ISAMP-ISTART))) ** 2.0
        ELSEIF(IUP .EQ. 2)THEN
          ANGINC =  PI  / NINC
          DO 2220 ISAMP=ISTART,IEND
 2220       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &                    (COS(HPI - ANGINC*(ISAMP-ISTART))) ** 2.0
        ENDIF
C-----------------------------------------------------------------------
C             TYPE = 3, HAMMING
C-----------------------------------------------------------------------
      ELSEIF(ITYPE .EQ. 3)THEN
        ANGINC = PI/(  NINC  )
        IF(IUP .EQ. 0)THEN
          DO 2300 ISAMP=ISTART,IEND
 2300       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &            (0.46 * (COS(0.0 + ANGINC*(ISAMP-ISTART))) + 0.54)
        ELSEIF(IUP .EQ. 1)THEN
          DO 2310 ISAMP=ISTART,IEND
 2310       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &            (0.46 * (COS( PI - ANGINC*(ISAMP-ISTART))) + 0.54)
        ELSEIF(IUP .EQ. 2)THEN
          ANGINC = 2.0*PI/NINC
          DO 2320 ISAMP=ISTART,IEND
 2320       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &            (0.46 * (COS( PI - ANGINC*(ISAMP-ISTART))) + 0.54)
        ENDIF
C-----------------------------------------------------------------------
C             TYPE = 4, COS ** 0.5
C-----------------------------------------------------------------------
      ELSEIF(ITYPE .EQ. 4)THEN
        ANGINC = PI/(2.0*NINC)
        IF(IUP .EQ. 0)THEN
          DO 2500 ISAMP=ISTART,IEND
 2500       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &                    (COS(0.0 + ANGINC*(ISAMP-ISTART))) ** 0.5
        ELSEIF(IUP .EQ. 1)THEN
          DO 2510 ISAMP=ISTART,IEND
 2510       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &                    (COS(HPI - ANGINC*(ISAMP-ISTART))) ** 0.5
        ELSEIF(IUP .EQ. 2)THEN
          ANGINC =  PI /  NINC
          DO 2520 ISAMP=ISTART,IEND
 2520       ARRAY(ISAMP) = ARRAY(ISAMP) *
     &                    (COS(HPI - ANGINC*(ISAMP-ISTART))) ** 0.5
        ENDIF
      ENDIF
 999  RETURN
      END




c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c

      subroutine phaserot(arr1,nsamps,degrees,work,maxwork)
      real arr1 (nsamps)
      real work(maxwork)

c
c     Phase rotation
c
c
c     arr1 - array 
c     nsamps - length of array
c     degrees - Degres to rotate
c
      parameter(nc=8192)
      pointer(ip_power , power)
      pointer(ip_phase , phase)
      pointer(ip_trace , trace)
      pointer(ip_ctrace , ctrace)
      real  power(nc)
      real  phase(nc)
      real trace(nc*2)
      complex ctrace(nc)

      ksamps=2
      do while(ksamps.lt.nsamps)
         ksamps=ksamps*2
      enddo
 
      nw=ksamps/2

      if(maxwork.lt.ksamps) then
         write(6,*) 'phaserot: Work buffer too small'
         return
      endif

      ip_trace = loc(work(1))
      ip_ctrace = ip_trace
 
      do i = 1,ksamps
        trace(i) = 0.0
      enddo

      do i = 1,nsamps
        trace(i) = arr1(i)
      enddo

c     ** forward transform
c

      call rwsrfft(trace,ksamps,+1)
 
      pi = atan (1.0) * 4.0

      fs = 2*pi*degrees/360.
      do iw = 1, nw
           ctrace (iw) = ctrace(iw) * cexp((0.0,-1.0) * fs)
      enddo


c     ** inverse transform
c

      call rwsrfft(trace,ksamps,-1)

      do i=1,nsamps
        arr1(i)=trace(i)
      enddo

 
      return
      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
c
      subroutine fildes(filt,maxfilt,work,maxwork,nfilt,sr,destype,
     $                  phase, ifrick,ilfreq,ihfreq,ildb,ihdb)

      character*(*) destype , phase
c
c     ** butterworth wavelet filterdesign
c     ** and ricker

      real filt(maxfilt) , work(maxwork)

      real flow,lowdb,fhigh,highdb

        logical bandpass
        logical lowpass
        logical highpass
        logical ricker
        logical zero
        logical minimum


      character elem*12,flag*4

      iret = -1


      bandpass = destype.eq.'BANDPASS'
      highpass = destype.eq.'HIGHPASS'
       lowpass = destype.eq.'LOWPASS'
      ricker   = destype.eq.'RICKER'

      iphase=1
      zero=.true.
      minimum=phase.eq.'MINIMUM'
      if (minimum) iphase=2


      if (.not. ricker) then

         if (bandpass) then
             flow   = ilfreq  
             lowdb  = ildb    
             fhigh  = ihfreq      
             highdb = ihdb         
             ifilt  = 1
         else if (lowpass) then
             fhigh  = ihfreq   
             highdb = ihdb     
             ifilt  = 3

         else if (highpass) then
             flow   = ilfreq   
             lowdb  = ildb    
             ifilt  = 2
         else
           return
         endif

         SAMP = sr / 1000.0
         CALL VPFILP(IFILT,IPHASE,SAMP,lowdb,flow,highdb,
     +               fhigh,filt,maxfilt,work,maxwork)


c        ** Determine the nr. of filter points necessary
c
         rlim=0.0001

         nfilt=257
         IF (IPHASE.EQ.2) nfilt=256
         call truncfilt(filt,nfilt,iphase,rlim,newnfilt)
         nfilt = newnfilt


         IF (IPHASE.EQ.2) THEN
           IZERO=1
         ELSE
           IZERO=nfilt/2+1
         ENDIF

      else
         rwlfrq=ifrick
         ifilt    = 4
         RWL = (SQRT(6.0)*1000.0) / (3.1415927*RWLFRQ)
         SAMP = sr     / 1000.0
         RWL  = RWL  / 1000.0
         CALL VPFILR(IPHASE,SAMP,RWL,filt,NUM,maxfilt,work,maxwork)
         nfilt=num
      endif

c     ** Normalise so that the rms value is unchanged
c
      call normfilt( filt,nfilt) 

      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine  normfilt( filt,nfilt) 
      real filt(nfilt)
      double precision energy , fac

      energy=0.0
      do i=1,nfilt
        energy = energy+filt(i)*filt(i)
      enddo

      fac = 1
      if(energy.gt.0.0) fac = dsqrt( dble(1)/energy)

      do i=1,nfilt
        filt(i) = filt(i)*fac
      enddo

      return
      END



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
c
         subroutine  truncfilt(filt,nfilt,iphase,rlim,newnfilt)
         real filt(nfilt)
c
c        Truncate filter
c
c        iphase=0 ++>   Mixed phase 
c        iphase=1 ==>   Zero phase
c        iphase=2 ===>  Minimum phase filter 
c        zero phase implies an odd nr. of filterpoints

         double precision energy,sum

         if(iphase.eq.0) then
           energy = 0.0
           do  i = 1 , nfilt
            energy = energy +  filt(i) * filt(i)
           enddo
           sum    = 0.0
           do i  = 1 , nfilt
             sum  = sum + filt(i)*filt(i)
             if (energy-sum.le.rlim*energy) then 
               newnfilt=i
               return
             endif
           enddo
         endif




         ifac=2
         icenter = nfilt/2+1
         IF (IPHASE.NE.1) then
           icenter=1
           ifac=1
         endif

c        ** Total energy
c
         energy = filt(icenter) *filt(icenter)
         do 333 i = icenter+1   , nfilt
           energy = energy + ifac * filt(i) * filt(i)
  333    continue

         kmax = 0   
         sum    = filt(icenter) * filt(icenter)
         do 444 i = icenter+1 , nfilt
           sum  = sum + ifac * filt(i)*filt(i)
           if (energy-sum.le.rlim*energy) goto 445
           kmax = kmax+1
  444    continue
  445    continue

         
         if (iphase.eq.1) then
           k=0
           do 555 i = icenter-kmax ,icenter+ kmax 
             k=k+1
             filt(k)=filt(i)
  555      continue
         else
           k = kmax
         endif

         newnfilt=k
         return
         end





c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      Subroutine ftbpick(trace,nstart,nsln,nrtyp,lv,nctyp,
     +              itrshp,iftbo,subsamp,ierr)
c======================================================================c
c
c Function   : Determine FTB detecting first amplitude value above
c              a certain percentage of rms value
c
c Parameters :
c         In :
c              trace  - Trace buffer                    (real(*))
c              nstart - Start sample for search         (integer)
c              nsln   - Search window length in samples (integer)
c              nrtyp  - Type of rough determination     (integer)
c                       1 - Peak
c                       2 - Trough
c                       3 - Threshold (absolute value)
c              lv     - Amplitude level to trig         (integer)
c              nrtyp  - Type of refinement              (integer)
c                       1 - No refinement
c                       2 - Closest extremum
c                       3 - Backtrack to opposite extremum
c                       4 - Backtrack to onset
c                       5 - Backtrack to zero
c                       6 - Backtrack through opposite extremum
c                           to onset
c              itrchp - Threshold level in per cent     (integer)
c         Out:
c              iftbo  - Estimated ftb in samples        (integer)
c              subsamp- Estimated ftb subsample         (real   )
c              ierr   - Error return                    (integer)
c                       0  - Ok
c                       >0 - see below
c
c Notes      :
c
c Subr-calls : none
c Func-calls : none
c W-AP-calls : none
c
c History    : Made from ftblv routine.
c              Date of coding: 16-11-1989
c              Company       : RWS/PKR
c              Modification  : None
c
c======================================================================c
c
      real*4 trace(1)
c
      ierr  = 0
      pct   = lv/100.0
      fsamp = nsln + 1
c
c     Check previous ftb
c
c     if ((iftb.lt.0).or.(iftb.gt.nsamp)) then
c        ierr = 1
c        go to 9999
c     end if
c
c     Compute rms value for search window
c
      rms   = 0.0
      isamp = 0
c
      do 100 i=nstart,nstart+nsln
         rms = rms + trace(i)**2
  100 continue
c     ** no live samples
      if(rms.eq.0.0) then
         ierr = 9
         iftbo = 1
         subsamp = 0.0
         go to 9999
      endif

      rms = sqrt(rms/fsamp)
c
c
c     Search for first amplitude value above pct*rms
c     ----------------------------------------------
c
      if(nrtyp.eq.1) then
         thrsh = pct*rms
         do 210 i=nstart,nstart+nsln
            if (trace(i).gt.thrsh) then
               isamp = i
               go to 300
            endif
  210    continue
c
      elseif(nrtyp.eq.2) then
         thrsh = -pct*rms
         do 220 i=nstart,nstart+nsln
            if (trace(i).lt.thrsh) then
               isamp = i
               go to 300
            endif
  220    continue
c
      elseif(nrtyp.eq.3) then
         thrsh = pct*rms
         do 230 i=nstart,nstart+nsln
            if (abs(trace(i)).gt.thrsh) then
               isamp = i
               go to 300
            endif
  230    continue
c
      endif
c
  300 continue
c
c     Check if isamp is reasonable
c
      if (isamp.eq.0) then
         ierr = 2
         iftbo = nstart
         subsamp = 0.0
c        call qerrmsg('No sample above triglevel')
         go to 9999
      else if (isamp.eq.nstart+nsln) then
         ierr = 3
         iftbo = nstart
         subsamp = 0.0
c        call qerrmsg('No samples above triglevel')
         go to 9999
      end if

c
c
c     FTB is now roughly detected
c     do refined determination according to nctyp
c     -------------------------------------------
c
      if (nctyp.eq.1) then
         iftbo = isamp
         subsamp = 0.0
      elseif (nctyp.eq.2.or.nctyp.eq.4) then
c
c        Search for first local extremum
c        -------------------------------
c
         a = abs(trace(isamp))
         iftbo = isamp
         do 1200 i=isamp+1,nstart+nsln
            if (abs(trace(i)).gt.a) then
               iftbo = i
               a = abs(trace(i))
            else
               go to 1210
            end if
 1200    continue
 1210    continue
         isamp = iftbo
c
      else if (nctyp.eq.3.or.nctyp.eq.6) then
c
c        Search for first previous extremum
c        ----------------------------------
c
c
c        Backtrack to zero-crossing
c
         asgn = sign(1.0,trace(isamp))
c
         do 1300 i=isamp,nstart,-1
            if (asgn*trace(i).lt.0.0) go to 1310
 1300    continue
c
 1310    isamp = i + 1
c
         a = asgn*trace(isamp)
         iftbo = isamp
         do 1330 i=isamp-1,nstart,-1
            if (asgn*trace(i).lt.a) then
               iftbo = i
               a = asgn*trace(i)
            else
               go to 1340
            end if
 1330    continue
 1340    continue
         isamp = iftbo
c
c
      else if (nctyp.eq.5) then
c
c        Backtrack to zero
c        -----------------
c
         eps = 0.1*abs(trace(isamp))
c
         do 1500 i=isamp,1,-1
 1500       if (abs(trace(i)).lt.eps) go to 1510
c
 1510    iftbz1 = i
         if(eps*trace(iftbz1).lt.0.0) iftbz1 = iftbz1 + 1
         iftbz2 = 0
c
         signb = trace(isamp)
         do 1530 i=isamp,1,-1
            if (signb*trace(i).lt.0.0) then
               iftbz2 = i
               go to 1540
            endif
 1530    continue
c
 1540    continue
         if(iftbz1.le.0) then
            if(iftbz2.le.0) then
               iftbo = isamp
            else
               iftbo = iftbz2
            endif
         else
            if(iftbz2.le.0) then
               iftbo = iftbz1 - 1
            else
               if(abs(iftbz1-iftbz2).lt.4) then
                  iftbo = iftbz2
               else
                  iftbo = iftbz1 - 1
                  if(iftbz1.gt.iftbz2) iftbo = iftbz1
               endif
            endif
         endif
c
      end if
c
      if (nctyp.eq.6.or.nctyp.eq.4) then
c
c        Backtrack to onset
c        ------------------
c
         eps = itrshp*abs(trace(isamp))/100.0
         iftbo = isamp
         asgn = sign(1.0,trace(isamp))
c
         do 1600 i=isamp,nstart,-1
            if (asgn*trace(i).lt.eps) go to 1650
            iftbo = i
 1600    continue
c
 1650    continue
         iftbo = iftbo - 1
      endif
c
c     Calculate subsample
c     -------------------
c
      if(nctyp.eq.4.or.nctyp.eq.6) then
         subsamp = (asgn*eps-trace(iftbo))/(trace(iftbo+1)-trace(iftbo))
      elseif(nctyp.eq.5) then
         subsamp = trace(iftbo)/(trace(iftbo)-trace(iftbo+1))
      elseif(nctyp.eq.2.or.nctyp.eq.3) then
         subsamp = 0.5*(trace(iftbo+1) - trace(iftbo-1))/
     +          (2*trace(iftbo) - trace(iftbo+1) - trace(iftbo-1))
      endif


 9999 continue

      Return
      End

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
c
      subroutine balance(input,output,scal,lsamp,ifsamp,winlen,
     $     eff,nwin,fb)

      real  input(lsamp)
      real  output(lsamp)
      real  scal(lsamp)
      integer winlen(nwin)
      real eff(nwin)
      real fb
c
c     input  (I) -  Input array
c     output (O) - output array
c     scal   (O) - Computed continous scaling function 
c     lsamp  (I) - The Number of samples in the arrays
c     ifsamp (I) - First sample to use in the input array
c     winlen (I) 
c     eff    (I)
c     nwin   (I)
c     fb     (I)

c
c     50 % overlapping window scaling of input data
c     The output will be a superposition of nwin windows
c


c     ** find feedback level 
c

      call  scalfact(input,lsamp,ifsamp,lsamp,yfb)


      do i=1,lsamp
         scal(i) = 0.0
      enddo


c     ** Do the computation
c
      do 305 iw = 1,nwin

              if (winlen(iw).eq.0.or.eff(iw).eq.0.0) then
                 write(6,*) 'balance : wrong wnd params'
                 return
              endif
                
c                       ** nr of samples pr window  / overlap length
c
              nsw = winlen(iw)
              nswolap = nsw/2

c                       ** nr of windows
c
              nswin  = (lsamp-ifsamp)/nswolap

              im2=0 
              y2 =0
              l1 = ifsamp
              l2 = l1 + nsw - 1
              l2 = min(l2,lsamp)
              if(l2.eq.lsamp) im2=lsamp+1

              c  = 0
c-lj 01.12.98 y1 = 1

              do 310 i = ifsamp , lsamp

c               ** compute new scaling factors
c
                if (i.ge.im2) then
                  call  scalfact(input,lsamp,l1,l2,y)
                  y1=y2
                  y2=y
                  if (y1.eq.0) y1=y2

                  c =(y2-y1)/nswolap

                  if (im2.eq.0) then
                    im1=ifsamp
                    im2=l1 + (l2-l1)/2.
                  else
                    im1 = im2
                    im2 = im1 + nswolap
                  endif

                  l1 = l1+nswolap
                  l1 = min(l1,lsamp)

                  l2 = l1 + nsw - 1
                  l2 = min(l2,lsamp)
                endif


                v = y1 + c * (i - im1)
                if(y.eq.0) v=0.0
                scal(i) = scal(i) + v * eff(iw)
310           continue
305   continue


c     ** do the scaling
c
      do  i=1,lsamp
        if(i.lt.ifsamp) scal(i) = scal(ifsamp)
        output(i)=input(i)*scal(i) + fb * yfb * input(i) 
      enddo

      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
       subroutine  scalfact(input,lsamp,l1,l2,y)
       real input(lsamp)

c     ** Average level of output trace
c
      parameter (avelevel=1000.)

      double precision ave,sum,small

c     data small/0.00001/
c     data big/9999999/

      sum = 0.

      y=0
      if(l2.le.l1) return 

      pavel = avelevel
      do 430 ii=l1,l2
           sum = sum +  abs(input(ii))
430   continue

      ave = sum/(l2-l1+1)

      if (ave.le.0.0) then
         y=0  
      else
        y    = pavel/ave
      endif

      return
      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine accuamp(inp,ninp,amp,nyq,sr,work)
      real inp(ninp) , amp(nyq) , work(1)
c
c
c     compute accumulated amplitude spectrum
c     integer frequencies only 

      maxl = 2048

1     if(ninp/2.lt.maxl) then
        maxl=maxl/2
        if(maxl.gt.256) goto 1
      endif


      ksamp = maxl*2

      k=0
      do 100 i=1,ninp 
         work(i) = inp(i)
  100 continue

      do 110 i = ninp+1, ksamp  
           work(i) = 0.0
  110 continue


c     ** Perform forward real to complex fft transform
c

      call rwsrfft(work , ksamp  , 1)


      do i=1,nyq
        amp(i) = 0.0
      enddo


c       ** Find the total energy 
c
      df      = 500./(sr*maxl)
      freq   = 0
      energy=0.0 
      knext = 1
      do i = 1,ksamp,2
              a = work(i)
              b = work(i+1)
          energy = energy + sqrt(a*a+b*b)
          freq = freq + df
          if (freq.ge.knext) then
            amp(knext) = energy
            knext = knext + 1
            if(knext.gt.nyq) goto 1000
          endif
      enddo         

1000  continue

      do i=1,nyq
        amp(i) = amp(i)/energy
      enddo
      amp(nyq) = 1
      return
      end
c


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine curve_fit( inp , ninp , curves , ncurves , icur)
c
      real inp(ninp) , curves(ninp*ncurves)
cc
c
c     Find best curve (least squares criteria)
c
      best = 1e30
      do i=1,ncurves
        iad = (i-1)*ninp 
        sum = 0.0
        do ip = 1,ninp 
          sum = sum + (curves(ip+iad) - inp(ip) )**2
        enddo
        ave = sqrt ( sum / float(ninp))
        if(ave.lt.best) then
            best = ave
            ibest = i
        endif
      enddo
      icur = ibest
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c     
c     
      subroutine projpnt(xa,ya,xb,yb, xc,yc, shotinc , 
     +                xnew,ynew,cerror,error,r,rnew,dx,dy,isp)

      double precision dist  , xsnew , ysnew

c       
c     Project the point (xc,yc) onto the line between
c     (xa,ya) and (xb,yb) The output point will be (xout,yout)
c       
c       
      double precision a1,a2,b1,b2,c1,c2,ba1,ba2,d1,d2,ab


c-lj 25.05.99  More North south problems . 
      logical northsouth , interior_point
c-lj  END


      real small
      data small/1e-08/  

c    |
c    |                            B
c    |          C
c    |
c    |
c    |                       (xout,yout)
c    |
c    |                    A
c    |
c    | ___________________________________________
c     

c    CD.AB=0   (Dot product)
c     
c    (b2-a2)/(b1-a1) = (d2-a2)/(d1-a1)  (Triangles)
c     
c         
      a1 = xa
      a2 = ya
 
      b1 = xb
      b2 = yb

      ba1 = b1-a1 
      ba2 = b2-a2 

      c1 = xc
      c2 = yc

      if(ba1.eq.0.0.and.ba2.eq.0.0) then
       ab = 0.0 
       d1=0.0
       d2=0.0
      else 
       ab = dsqrt(ba1**2.0 + ba2**2.0)
       d2  = (ba1*c1 + ba2*c2) * ba2 - (ba2*a1 -ba1*a2)*ba1
       d2  = d2/(ba2*ba2+ba1*ba1)

c-lj 20.05.99 North south produced NaN
       if( abs(ba1).lt.small) then
          d1 = (b1+a1)/2.
          northsouth=.true. 
       else
          d1  = (ba1*c1 + ba2*c2 -ba2*d2)/ba1
          northsouth=.false. 
       endif
c-lj
      endif 
      
     
      xsnew = d1
      ysnew = d2

c
c        xnew   -  x value of (sx,sy) point projected onto line
c        ynew   -  y value of (sx,sy) point projected onto line
c        cerror  -  Crossline error
c        error  -  Total error
c        r      - Radial distance from SOL to projected point
c         rnew   - Radial distance after regridding
c         dx     - Change in x-coordinate in line direction as a result of regridding
c         dy     - Change in y-coordinate in line direction as a result of regridding
c        isp     - Shotpoint
c




c          ** compute distance to first point
c
      dist = dsqrt ( (xsnew-xa)*(xsnew-xa) + (ysnew-ya)*(ysnew-ya) )
      r=dist

c-lj 25.05.99
cWhat about this ???
c     use_y = abs(ya-yb).gt.abs(xa-xb) 
c End whatabout
c     if(xsnew.lt.xa) r=-1.0*r
      if(northsouth) then
         aaamin = min(ya,yb)
         aaamax = max(ya,yb)
         interior_point = ysnew.ge.aaamin.and.ysnew.le.aaamax

         if(.not.interior_point)then
             if(abs(ysnew-ya).lt.abs(ysnew-yb)) r=-1.0*r 
         endif
      else
         aaamin = min(xa,xb)
         aaamax = max(xa,xb)
         interior_point = xsnew.ge.aaamin.and. xsnew.le.aaamax

c        ** LEAVE it as it always has been
c         if(.not.interior_point)then
c             if(abs(xsnew-xa).lt.abs(xsnew-xb)) r=-1.0*r 
c         endif

          if(xsnew.lt.xa) r=-1.0*r
          
      endif
c-lj End      


      isp = nint (r/shotinc + 1.0)

c     ** radial distance after regridding
c
      rnew =  (isp-1)*shotinc 




c     ** determine crossline sign
c     ** (use vector arithm!!)
c
c-LJ 13.02.01 USE Y-coordinates only , except when NORTH-SOUTH
      idirc = 1 
      if (northsouth) then
        if(xc.gt.xsnew) idirc = -1
      else
        if(yc.lt.ysnew) idirc = -1
      endif
c      abx = abs(xc-xsnew)
c      aby = abs(yc-ysnew)
c      if(abx.gt.aby) signc = (xc-xsnew) 
c      if(aby.gt.abx) signc = (yc-ysnew) 
c      if(signc.lt.0.0) idirc=-1
c-lj end      


      

c     ** distance moved when projected 
c
      cerror  =  (xsnew-xc)**2 + (ysnew-yc)**2 
      if (cerror.gt.small) then
          cerror = idirc * sqrt(cerror)
      else 
          cerror =0.0
      endif 

      rmov = rnew-r


c     ** total distance moved
c
      idiri = 1 
      if(r.gt.rnew) idiri=-1
      error=cerror*cerror + rmov*rmov
      if(error.gt.small) then 
         error = idiri*sqrt(error)
      else
         error=0.0 
      endif 

c     ** x and y increments after regridding
c
      if(r.ne.0.0) then
        dx = (xsnew-xa)*(rnew/r-1.)
        dy = (ysnew-ya)*(rnew/r-1.)
      else
        dx = 0.0 
        dy = 0.0
      endif

      xnew = xsnew
      ynew = ysnew 

      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine ljpolfit(x,y,nr,xv,yv,icm,istat)
      real x(nr),y(nr)
c
c
c     Polynomial fitting by overlapping second degree polynomials
c
c
c     If icm>1 redo coefficient computation
c     if icm=2 Test on i1,i2 
c     if icm=3 Linear interpolation



c     ** X must be sorted ( The smallest value first)
c

      real*8 c1(3) , c2(3)
      real*8 a0,a1,a2
      equivalence (c1(1) , a2)
      equivalence (c1(2) , a1)
      equivalence (c1(3) , a0)

ccc   real warr(100),aerr(100)
      data a0,a1,a2/3*0.0/

      data i1old/0/


ccc   data icall,ndeg/0,3/
ccc   if(icall.eq.0) then
ccc      do i=1,100
ccc        warr(i)=1.0
ccc      enddo
ccc      icall=1
ccc   endif



      istat=-1

      ix1=-1

      if(nr.le.1) return

      if(icm.ge.1) then


c     ** Linear
c
        if(nr.eq.2) then

          a1 = ( y(1) - y(2) )/ ( x(1) - x(2) )
          a0 =   y(1) - a1*x(1)

        elseif (nr.ge.3) then

c       ** Find the two three point polynomials to merge
c
          i1=-1
          i2=-1
          i=0
          do while (i.le.nr-2.and.i1.eq.-1)
            i=i+1
            if(x(i).gt.x(i+1)) then
cc            write(6,*)   'ljpolfit x-array not sorted'
cc            write(6,*) i,x(i),x(i+1)
              istat=-2              
              return 
            endif   

            if(xv.ge.x(i).and.xv.le.x(i+1) ) then
              i1=i-1
              i2=i
            endif

          enddo

          if(i1.eq.-1) then
            write(6,*) 'ljpolfit x-value out of bounds',xv
            if(xv.lt.x(1))  i1 = 1
            if(xv.gt.x(nr)) i1 = nr-2
          endif

          i1 = max(1,i1)
          i2 = max(1,i2)

          i1 = min(i1,nr-2)
          i2 = min(i2,nr-2)


c         ** conditional recomputation
c
          if (icm.eq.1.or.i1.ne.i1old) then

ccc         call polfit(x(i1) , y(i1) ,3, ndeg, warr,aerr,rmserr,xv,yv1)
             call lj3d( x(i1) , y(i1) , c1)
             if(i2.ne.i1) then
ccc            call polfit(x(i2),y(i2),3,ndeg,warr,aerr,rmserr,xv,yv2)
               call lj3d( x(i2) , y(i2) , c2)
             else
               do i=1,3
                c2(i)=c1(i)
              enddo
ccc           yv2 = yv1
             endif
          endif
          i1old = i1
        endif

      endif

      if(nr.le.3) then
        yv = a2*xv*xv + a1*xv + a0
      else

        yv1 = c1(1)*xv*xv + c1(2)*xv + c1(3)
        yv2 = c2(1)*xv*xv + c2(2)*xv + c2(3)  

cx      i1    i2  
c       X     X  xv X    X
c                ^
c                |
c  P1   .............
c  P2         ............


        ccoef = (xv - x(i2))/(x(i2+1)-x(i2))

        yv = yv1*(1.-ccoef) + yv2*ccoef
      endif

 9999 istat=0
      return
      end


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine lj3d(xx,yy,c)
      real xx(3),yy(3)
      real*8 c(3),x(3),y(3)
c
c     Solve equation system
c
c     a2*x1**2 + a1*x1 + a0 = y1
c     a2*x2**2 + a1*x2 + a0 = y2
c     a2*x3**2 + a1*x3 + a0 = y3
c
      real*8 aa,bb1,bb2,a2,a1,a0

      do i=1,3
        x(i) = xx(i)
        y(i) = yy(i)
      enddo


c     aa =(y(1)-y(2))*(x(1)-x(3))-(y(1)-y(3))*(x(1)-x(2))
c     bb1 =(x(1)*x(1) - x(2)*x(2)) * (x(1) - x(3))
c     bb2 =(x(1)*x(1) - x(3)*x(3)) * (x(1) - x(2))

      aa =(y(1)-y(2))*(x(3)-x(2))-(y(3)-y(2))*(x(1)-x(2))
      bb1 =(x(1)*x(1) - x(2)*x(2)) * (x(3) - x(2))
      bb2 =(x(3)*x(3) - x(2)*x(2)) * (x(1) - x(2))



      if(bb1.eq.bb2.or.x(1).eq.x(2) ) then
        write(6,*) 'ljpolfit undetermined'
        c(1)=0.0
        c(2)=0.0
        c(3)=0.0
      else
        a2 = aa/(bb1 - bb2)
        a1 = ( y(1) - y(2) - a2*(x(1)*x(1)-x(2)*x(2))  ) /( x(1)-x(2) )
        a0 = y(1) - a2*x(1)*x(1) - a1*x(1)
        c(1) = a2 
        c(2) = a1 
        c(3) = a0 
      endif
      return
      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine linfit(x,y,nr,xv,yv,istat)
      real x(nr),y(nr)
c
c
c     Linear interpolation
c
c


c     ** X must be sorted 
c


      istat=-1

      ix1=-1
      if(nr.le.1) return

      idir=1
      if(x(nr).lt.x(1))idir=-1


c       ** Find the two points
c
      ifound=0
      do i=1,nr-1
         if(idir.eq.1) then
          if(xv.ge.x(i).and.xv.le.x(i+1))  ifound=i
         else
          if(xv.ge.x(i+1).and.xv.le.x(i))  ifound=i
         endif
         if(ifound.gt.0) goto 50
      enddo

   50 if(ifound.eq.0) then
         write(6,*) 'linfit x-value out of bounds',xv
         return
      endif

c         ** linear interpolation
c

       x1 = x(ifound)
       x2 = x(ifound+1)
       y1 = y(ifound)
       y2 = y(ifound+1)
       if (x2.eq.x1) then
           c = 0.0
       else
           c  = (y2-y1)/(x2-x1)
       endif
       yv = y1 + c*(xv-x1)
       istat=0
       return
       end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine hilbt_frq( cchb , ksamp , sr )
      complex cchb(ksamp)
c
c     ** Frequency domain hilbert transformation
c
c      
c     sr sample rate in seconds
c
      double precision dsr 
      data pi/0.0/
      if(pi.eq.0.0) pi= atan(1.0)*4

c         ** complex FFT 
c
      call rwscfft(cchb,ksamp,1)

      nw      = ksamp/2
      dsr     = sr


      do iw = 1, nw
         cchb(iw) = cmplx(0.0,1.0)*cchb(iw)
      end do
      
      do iw = nw+1, nw*2-2
         cchb(iw) = cmplx(0.0,-1.0)*cchb(iw)
      end do
 
c     ** Inverse complex FFT 
c     ** Time-domain hilbert transform
c
      call rwscfft(cchb,ksamp,-1)

      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine cmplx_ddt( ccdt , ksamp , sr )
      complex ccdt(ksamp)
c
c     ** complex derivation by multiplication with jw
c
c      
c     sr sample rate in seconds
c
      double precision dsr , dw  
      data pi/0.0/
      if(pi.eq.0.0) pi= atan(1.0)*4



c         ** complex forward FFT 
c
      call rwscfft(ccdt,ksamp,1)

      nw      = ksamp/2
      dsr     = sr
      dw =     pi/(dsr*nw)

c         ** Derivation. Positive frequencies first
c
      do iw = 1, nw  
           w = -1.0*(iw-1)*dw
           ccdt(iw)  =  ccdt(iw) * cmplx(0.0 , w)
      enddo

c         ** Derivation.Negative frequencies
c
      do iw = 2, nw  
           w = (iw-1)*dw
           ccdt(ksamp+2-iw)    =  ccdt(ksamp+2-iw)* cmplx(0.0,w)
      enddo

      ccdt(1)    = cmplx(0.0,0.0) 
      ccdt(nw+1) = cmplx(0.0,0.0) 


c     ** amplitude correction
c
      do i=1,ksamp
          ccdt(i) = ccdt(i)/float(nw)
      enddo


c         ** Inverse complex FFT 
c         ** Time derivative of the complex trace
c
      call rwscfft(ccdt,ksamp,-1)

      return
      end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine cmplx_intg( ccdt , ksamp , sr )
      complex ccdt(ksamp)
c
c     ** integration by multiplication with 1.0/jw
c
      
c     sr sample rate in seconds
      double precision dsr , dw  
      data pi/0.0/

      if(pi.eq.0.0) pi= atan(1.0)*4


c         ** complex FFT 
c
      call rwscfft(ccdt,ksamp,1)


      nw      = ksamp/2
      anw     = float(nw)
      dsr     = sr
      dw =     pi/(dsr*nw)


c         ** positive frequencies first
c
      do iw = 2, nw  
           w = -1.0*(iw-1)*dw
           ccdt(iw)    = anw* ccdt(iw) / cmplx(0.0 , w)
      enddo

c         ** then negative frequencies

      do iw = 2, nw  
          w = (iw-1)*dw
          ccdt(ksamp+2-iw)    =  anw*ccdt(ksamp+2-iw)/ cmplx(0.0,w)
      enddo

      anyq = 0.5/dsr
      ihz1 = NINT(1.0*anw/anyq)
      ihz4 = NINT(4.0*anw/anyq)

      do iw=1,ihz1
        ccdt(iw)        = cmplx(0.0 , 0.0)
        ccdt(ksamp-iw+1)  = cmplx(0.0 , 0.0)
      enddo

      npp = ihz4-ihz1
      dv= pi/npp
       v=pi
      do iw = ihz1+1, ihz4
          v = v-dv
          weight = 0.5 + 0.5*cos(v)
          ccdt (iw)        = ccdt(iw)       * weight           
          ccdt (ksamp-iw+1)  = ccdt(ksamp-iw) * weight           
      enddo


c         ** Inverse complex FFT 
c
      call rwscfft(ccdt,ksamp,-1)

      return
      end

c
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-cc-c-c-c-c-c-c-c-c-
c     
c     
      subroutine hilbt_construct(hdif,nhilb)
c     
c     Construct hilbert differentiator
c     

      real hdif(nhilb)

      pi     = atan(1.0)*4.0
      nhilb2 = nhilb/2

      do  i=1,nhilb
         hdif(i)=0.0
      enddo  

      do i = 1 , nhilb2-1 , 2
         hdif(i+nhilb2)=-2.0/(pi*i)
      enddo

c     ** anti-symmetrical around t=0
c     
      do i=1,nhilb2-1,2
         hdif(nhilb2-i)=-1.0*hdif(i+nhilb2)
      enddo
      return
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
c
      subroutine regression(x,y,w,n,x1,y1,x2,y2,aa,bb)
c
c
c     Fit a straight line to x , y values
c     w is the weight function
c
      real x(n) , y(n) , w(n)

      double precision sumx,sumy,sumxx,sumxy,xx,yy,a,b
      double precision sumw , ww

      sumw=0.0
      sumx=0.0
      sumy=0.0
      sumxx=0.0
      sumxy=0.0

      xmn = 1e30
      xmx =-1e30

      do i=1,n
        xx=x(i)
        yy=y(i)
        ww=w(i)
        sumw = sumw + ww
        sumx = sumx + xx*ww
        sumy = sumy + yy*ww
        sumxx= sumxx + xx * xx * ww
        sumxy= sumxy + xx * yy * ww

        if(xx.gt.xmx) xmx=xx
        if(xx.lt.xmn) xmn=xx
      enddo

c     y = ax+b  nb !!!!!!!!

      b =  (sumxy*sumx - sumy*sumxx)/(sumx*sumx-sumw*sumxx)
      a =  (sumy - sumw*b)/sumx


      x1 = xmn
      y1 = a*xmn+b 
      x2 = xmx
      y2 = a*xmx+b 

      aa=a
      bb=b
      return
      end




c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
       subroutine autocorr(a,na,c,nc)
       real a(na),c(nc)
       real*8 sum
       do k = 1 , min(na,nc)
          sum = 0.0
          do j = 1 , na
            sum = sum + a(j)*a(k-1+j)
          end do
          c(k) = sum
       end do
       return
       end

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine ccrosscorr(a,nl,b,nb,c,ialen)
      complex a(nl) , b(nb) , c(ialen)  
c
c     Complex same domain crosscorrelation
c
c
      complex*16 sum

      do 502 i=1,ialen
         sum = cmplx(0.0 , 0.0)
         k  =0
         do 602 j= 1 , nb
            if (k+i.gt.nl) goto 603
            sum = sum + b(j) * CONJG(a(k+i))
            k = k + 1
 602     continue
 603     c(i) = sum 
 502  CONTINUE     
      return
      end
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine crosscorr(a,nl,b,nb,c,ialen)
      real a(nl) , b(nb) , c(ialen)  
c
c     Real same domain crosscorrelation
c
c
      real*8 sum
      do 502 i=1,ialen
         sum=0.
         k  =0
         do 602 j= 1 , nb
            if (k+i.gt.nl) goto 603
            sum = sum + b(j) * a(k+i)
            k = k + 1
 602     continue
 603     c(i) = sum 
 502  CONTINUE     
      return
      end



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-cc--c-c-c-cc-c-c-c-c-c-c-c-c-c-c-c-cc-c-c-
c
c
          subroutine tcrosscorr(trace1,nsamp,pulse,npulse,cc_trpu,ncc,
     $                      ifound , subsamp , cmax)
          real trace1(nsamp) , pulse(npulse) , cc_trpu(ncc)
          real*8 qa,qb,qc,xpeak
c
c
c         Crosscorrelate two time series in the time domain.
c         Return crosscorrelation function and sample,subsample and max correlation
c
c
          call crosscorr(trace1,nsamp,pulse,npulse,cc_trpu,ncc)


          call findpeak(cc_trpu , nsamp , cmax , ifound,subsamp)

c-------------------------------------------------------------------------
c        alternative way of finding subsamp , but cmax is normalized
c
         if(ifound.gt.maxdel.and.ifound.lt.nsamp-npulse-maxdel) then
         maxdel = 10
         call  maxcorr (pulse,npulse,trace1(ifound-maxdel),
     $                  npulse+2*maxdel,idelay,subsamp,cmax_dummy)
         endif
c-------------------------------------------------------------------------



         return
         END



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
          subroutine findpeak(arr,nsamp,ampmax,ifound,subsamp)
          real arr(nsamp)
c
c         find maximum value of a time series
c
c       
c
          small = 1e-18  
          ifound = 0
          ampmax = 0
          subsamp = 0
          amsmax = 0.0

          if(nsamp.lt.3)return

          do  jk = 2 , nsamp-1
             amsb = abs(arr(jk-1))
             ams  = abs(arr(jk))
             amsa = abs(arr(jk+1))
             if (ams.gt.amsmax.and.ams.ge.amsb.and.ams.ge.amsa) then
                amsmax = ams
                ifound = jk
             endif
          enddo
          if(ifound.eq.0) return


c-----------------------------------------------------------------------------
c        ** find subsample (second order polynomial)
c     

          qa = 0.5*
     $        (arr(ifound-1)+arr(ifound+1)-2*arr(ifound))
          qb = 0.5*(arr(ifound+1)-arr(ifound-1))
          qc = arr(ifound)

c        **  x= -b/(2a)
c     
          denom= 2*arr(ifound)-arr(ifound+1)-arr(ifound-1)

          if(abs(denom).lt.small)then
             ifound=0
             return
          ENDIF

          xpeak=dble(0.5)*( arr(ifound+1)  - arr (ifound-1) )/
     $        (2*arr(ifound) - arr(ifound+1)-arr(ifound-1))

          subsamp = xpeak
          ampmax = qa*xpeak*xpeak + qb*xpeak + qc

          return
          END




c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
c
      subroutine fcrosscorr(a,na,b,nb,c,nc,work,nwork)
      real a(na), b(nb) , c(nc) 
      real work(nwork)
      real*8 sum
c
c     frequency domain crosscorrelation
c
c
      pointer (ip1,cp1)
      pointer (ip2,cp2)
      complex cp1(100000) , cp2(100000)      
c
c
c     ** work should be big enogh to hold 
c     ** two fourier transforms  
c     ** on output work will contain the fourier transform of the crosscorr 
c
      nn = nwork/2+1

      ip1 = loc (work(1))
      ip2 = loc (work(nn))

      do i=1,nwork
       work(i)=0.0
      enddo

      do i=1,na
       work(i)=a(i)
      enddo

      do i=1,nb
       work(i+nn-1) = b(i)
      enddo

      nfft = max(na,nb)
      ksamp=2
      do while(ksamp.lt.nfft)
       ksamp=ksamp*2
      enddo

      call rwsrfft(cp1,ksamp,1)
      call rwsrfft(cp2,ksamp,1)

      do i=1,ksamp/2
        cp2(i) = cp1(i) * CONJG( cp2(i))
      enddo

      call rwsrfft(cp2,ksamp,-1)

      do i=1,nc   
        c(i) = work(i+nn-1)
      enddo

      return
      end




c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      subroutine corrppsf(a,b,c,nn,work,nwork,srsec,
     $                    fl,fh,tstart,tend,ppsfk,cmd)
c
c     Pure phase shift correction (to remove harmonics)
c     before crosscorrelation   
c
c
c     a      - array containing geophone
c     b      - Sweep signal
c     c      - Output array
c     fl     - Sweep start frequency
c     fh     - Sweep end frequency
c     tstart - Sweep start time
c     tend   - Sweep end time 
c     ppsfk  - kfactor (Li et. al Geophysics march/april 1995)
c     cmd    - What kind of data should be output

      character *(*) cmd   

      real a(nn), b(nn), c(nn)
      real work(nwork)


      pointer (ip1,cp1)
      pointer (ip2,cp2)
      pointer (ip3,cp3)
      complex cp1(100000) , cp2(100000) ,cp3(100000)

 
      complex ci
      real rpart
      real k

      k = ppsfk 

      ci = cmplx (0.0 , 1.0)

      ksamp=2 
      do while (ksamp.lt.nn)
         ksamp = ksamp*2
      enddo

      ksamp = ksamp*2

      ip1 = loc (work(1))
      ip2 = loc (work(ksamp*2+2))
      ip3 = loc (work(ksamp*4+4))


      do i=1,ksamp
         cp1(i) = cmplx(0.0,0.0)
         cp2(i) = cmplx(0.0,0.0)
         cp3(i) = cmplx(0.0,0.0)
      enddo


      pi = 4*atan(1.0)


      ifs = tstart/srsec + 1
      ils = tend/srsec   + 1 
      nps = ils - ifs    + 1 

      TT = (nps+1)*srsec   ! sweep length


c     ** t=0 at ksamp/2

      iad = ksamp/2



c     ** construct the sk time function
c
      do i = ifs , ils
         t = (i-ifs)*srsec
         
c        ** instantaneous phase
c
         fi = (k+1)*fl*t +(k+1)*(fh-fl)*t*t/(2*TT)


c        ** time representation of phase filter
c
         rpart = exp(2*pi*ci*(k+1)*fi)
         cp1(iad+i-ifs+1) = cmplx ( rpart , 0.0)

      enddo


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     ** Work on sweep only ??  
c
      do i=1,nps
       cp2 (i+iad) = cmplx ( a(i+ifs-1) ,0.0)
       cp3 (i+iad) = cmplx ( b(i+ifs-1) ,0.0)
      enddo


 

      isign = 1
      call rwscfft(cp1,ksamp,isign)
      call rwscfft(cp2,ksamp,isign)
      call rwscfft(cp3,ksamp,isign)

 
c     ** do the phase shift operation
c
      do i=1,ksamp
        cp2(i)= cp2(i) * CONJG( cp1(i))
        cp3(i)= cp3(i) * CONJG( cp1(i))
      enddo


c     ** reverse complex fft
c
      isign=-1
      call rwscfft(cp2,ksamp,isign)
      call rwscfft(cp3,ksamp,isign) 



c     ** zero the "negative time noise"
c     
      if(cmd.eq.'NOISE') then

         do i= iad+1 , ksamp
            cp2(i+iad)     = cmplx(0.0,0.0)
            cp3(i+iad)     = cmplx(0.0,0.0)
         enddo

      else

c     ** zero the "negative time noise"
c     
         do i= 1 , iad
            cp2(i+iad)     = cmplx(0.0,0.0)
            cp3(i+iad)     = cmplx(0.0,0.0)
         enddo
      endif



c     ** forward fourier transform again
c
      isign=1
      call rwscfft(cp2,ksamp,isign)
      call rwscfft(cp3,ksamp,isign) 


c     ** inverse ppsf
c
      do i=1,ksamp
        cp2(i)= cp2(i) *  cp1(i)
        cp3(i)= cp3(i) *  cp1(i)
      enddo


c     ** crosscorrelate trace and sweep
c
      do i=1,ksamp
        cp2(i) = cp2(i)*CONJG(cp3(i))
      enddo 


c     ** reverse complex fft
c
      isign=-1
      call rwscfft(cp2,ksamp,isign)


      do i=1,nn
         c(i) = real ( cp2(i))
      enddo

      return
      end
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-cc-c-c-c-c-c-c-c-c-c-c-
c     
c     
      subroutine moveout_c(time_in,dts,velocity,
     $     sx,rx,sy,ry,rz,time_out)

      dxsr   = sx-rx
      dysr   = sy-ry
      ofst   = sqrt(dxsr*dxsr+dysr*dysr)
      vel    = velocity*dts

      z      = time_in*vel*0.5

      z0     = amax1(z,rz+10.0)
      t0     = (z0/vel)

      if (rz.le.10.0) then 
         t = sqrt(t0*t0+(ofst*ofst*0.25)/(vel*vel))*2.0
      else 
         t0r    = ((z0-rz)/vel)
         tshift =(z0-z)/vel
         dzr    = z0-rz
         xs     = 0.0
         xr     = ofst
         x0     = ((dzr*xs+z0*xr)/(dzr+z0))
         xdiff  = xr-x0
         t      = sqrt(t0*t0+((x0*x0)/(vel*vel)))
     &        + sqrt(t0r*t0r+((xdiff*xdiff)/(vel*vel)))
         t      = amax1(t - 2.0*tshift,1.00001)

      end if

      time_out = t

      return
      END
c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
c
      subroutine monk(mt, mdt, ht, hdt, dt,  nt,
     $     itstart,itend,w)


c
c
c
c     Compute the four MONK weights 
c
c

c     ** input arguments [all arrays in time domain]
c
      real  mt(nt)            ! model trace [mode2]
      real  mdt(nt)           ! time deriv of model trace 
      real  ht(nt)            ! hilbert of model trace 
      real  hdt(nt)           ! derivative of hilbert
      real  dt(nt)            ! original trace [mode1] 

      real  w(4)              ! monk coefs (Output)
      integer nt              ! number of time samples
      integer itstart         ! time window start 
      integer itend           ! time window end 


c     ** local parameters
c
      double precision a(0:4), e(0:4)
      double precision zeroval , r

      zeroval = 1.0e-35

      do i=0,4
         a(i)=0.0
         e(i)=0.0
      enddo
      

C     ** Zero lag correlation of model trace
c
      do it = itstart, itend
         a(0) = a(0) + mt(it)*mt(it)
      end do 

C     ** Zero lag correlation of model trace and Hilbert transform derivative
      do it = itstart, itend
         a(3) = a(3) + mt(it)*hdt(it)
      end do 

C     ** Zero lag autocorrelation of derivative of model trace
c
      do it = itstart, itend
         a(4) = a(4) + mdt(it)*mdt(it)
      end do 


c     ** crosscorrelations between data trace and the four model
c     ** derived traces
c
      do it = itstart, itend
         e(1) = e(1) + mt(it)*dt(it)
      end do 

      do it = itstart, itend
         e(2) = e(2) + mdt(it)*dt(it)
      end do 

      do it = itstart, itend
         e(3) = e(3) + ht(it)*dt(it)
      end do 

      do it = itstart, itend
         e(4) = e(4) + hdt(it)*dt(it)
      end do 

      r  = a(0)*a(4) - a(3)*a(3)

      if ((a(0).lt.zeroval).or.(abs(r).lt.zeroval)) then
         w(1) = 1.0
         w(2) = 0.0
         w(3) = 0.0
         w(4) = 0.0
      else
         w(1) = (e(1)*a(4)-e(4)*a(3))/r
         w(2) = (e(2)*a(0)+e(3)*a(3))/r
         w(3) = (a(4)*e(3)+a(3)*e(2))/r
         w(4) = (a(0)*e(4)-a(3)*e(1))/r
      end if
 
      return
      end 



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c     
      subroutine orig_corrppsf(a,b,c,nn,work,nwork,srsec,
     $                    tstart,tend,fl,fh,ppsfk)

 

c    ** a - contains original trace
c    ** b - Will be pass signal
c    ** c - Will be reject signal

      
         
      real a(nn) , b(nn) , c(nn)
      real work(nwork)
      
      pointer (ipf,cpf)
      pointer (ipa,cpa)
      pointer (ipb,cpb)
      pointer (ipc,cpc)
        
      parameter (iasize=100000)
      complex cpf(iasize) , cpa(iasize), cpb(iasize), cpc(iasize)
      

      complex ci
      
      real k
      
      k = ppsfk
      
      ci = cmplx (0.0 , 1.0)
         
      ksamp=512
      do while (ksamp.lt.nn)
         ksamp = ksamp*2
      enddo
      
c     ksamp=ksamp*2
      
      if(nwork.lt.ksamp*8) then
         write(6,*) ' Increase work buffers'
         return
      endif

        
      ipf = loc (work(1))
      ipa = loc (work(ksamp*2+2))
      ipb = loc (work(ksamp*4+2))
      ipc = loc (work(ksamp*6+2))
      

      do i=1,ksamp
         cpf(i) = cmplx(0.0,0.0)
         cpa(i) = cmplx(0.0,0.0)
         cpb(i) = cmplx(0.0,0.0)
         cpc(i) = cmplx(0.0,0.0)
      enddo 

      pi     = 4*atan(1.0)

      ifs = tstart/srsec + 1
      ils = tend/srsec   + 1
      nps = ils - ifs    + 1 
      TT = (nps+1)*srsec



      
c     ** t=0 at ksamp/2
      
      iad = ksamp/2       
         
      

c     ** construct the sk time function
c    
      do i=1,nps
         t = (i-1)*srsec
      
c        ** instantaneous phase
c     
         fi = (k+1)*fl*t +(k+1)*(fh-fl)*t*t/(2*TT)
      
c        ** time representation of phase filter
c     
        rdel =  exp(2*pi*ci*(k+1)*fi)
c       cpf(i+iad) = cmplx ( exp(2*pi*ci*(k+1)*fi) , 0.0)
        cpf(i+iad) = cmplx ( rdel , 0.0) 
      enddo
         
      
      do i=1,nps
       cpa (i+iad) = cmplx ( a(i+ifs-1) ,0.0)
      enddo
      
      do i=1,nn
         a(i) = 0.0
         b(i) = 0.0
         c(i) = 0.0
      enddo
         
      
      isign = 1
      call rwscfft(cpf,ksamp,isign)
      call rwscfft(cpa,ksamp,isign)
      

c     ** do the phase shift operation
c        
      do i=1,ksamp
        cpa(i) = cpa(i) * CONJG( cpf(i))
      enddo 
      

c     ** reverse complex fft
c     
      isign=-1
      call rwscfft(cpa,ksamp,isign)


c     ** zero negative time
c     
      do i= 1,iad
      
         cpc(i)      = cmplx( 0.0 , 0.0)
c        cpc(i+iad)  = cmplx( real(cpa(i+iad)) , 0.0)
         cpc(i+iad)  = cpa(i+iad)
      
         cpb(i)      = cpa(i)
c        cpb(i)      = cmplx( real(cpa(i)) , 0.0)
         cpb(i+iad)  = cmplx( 0.0,0.0)
       
      enddo
      
      
c     ** forward fourier transform again
c     
     
      isign=1  
      call rwscfft(cpa,ksamp,isign)
      call rwscfft(cpb,ksamp,isign)
      call rwscfft(cpc,ksamp,isign)
      
      
c     ** inverse ppsf
c      
      do i=1,ksamp
        cpa(i)= cpa(i) *  cpf(i)
        cpb(i)= cpb(i) *  cpf(i)
        cpc(i)= cpc(i) *  cpf(i)
      enddo


      
c     ** reverse complex fft
c     
      isign=-1
      call rwscfft(cpa,ksamp,isign)
      call rwscfft(cpb,ksamp,isign)
      call rwscfft(cpc,ksamp,isign)
      
       
      do i=1,nps
        a(i+ifs-1) = real ( cpa(i+iad))
        b(i+ifs-1) = real ( cpb(i+iad))
        c(i+ifs-1) = real ( cpc(i+iad))
      enddo
     
      
      return
      end  
