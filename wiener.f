c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPWIEN (WAVLT,  LWAVLT, DOUT,LDOUT,
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
          !WRITE(6,*)' THE FILTER IS TOO LONG , OR '
          !WRITE(6,*)' THE WAVELETS ARE TOO LONG.'
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

      return

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



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPMAX(ARRAY, NSAMPS, MAX, NMAX)
C
C        ROUTINE TO RETURN THE MAX. VALUE IN AN ARRAY AND ITS INDEX
C        (IF NOT UNIQUE, RETURNS FIRST INDEX FOR MAX. VALUE)
C
      DIMENSION  ARRAY(*)
      REAL       MAX
C
      MAX = -1.0E+33
C
      DO 10 I=1,NSAMPS
         IF(ARRAY(I).GT.MAX) THEN
            MAX = ARRAY(I)
            NMAX = I
         ENDIF
   10 CONTINUE
      RETURN
      END


c-c-c-c-c-c-c-c-c-c--cc-c--cc-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c

      SUBROUTINE VPCLR(A, IAINDX, IAINC, NSAMPS)
C
C     ROUTINE TO CLEAR NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT AN
C     INCREMENT IAINC
C     REPLACES AP ROUTINE VCLR
C
C     VERSION: JB01 26/02/85  INITIAL RELEASE                J. BARNES
C----------------------------------------------------------------------
C
      DIMENSION A(*)
C
      DO 10 I=0,NSAMPS-1
         A( IAINDX+(I*IAINC) ) = 0.0
   10 CONTINUE
      RETURN
      END



c-c-c-c-c-c-c-c-c-c--cc-c--cc-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPSMUL(A, IAINDX, IAINC, SCALAR, R, IRINDX, IRINC,
     +                  NSAMPS)
C
C     ROUTINE TO MULTIPLY NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT AN
C     INCREMENT IAINC BY SCALAR.  RESULT RETURNED IN ARRAY R (STARTING
C     AT ELEMENT IRINDX AT AN INCREMENT IRINC).
C     REPLACES AP ROUTINE VSMUL
C
C     VERSION: JB01 26/02/85  INITIAL RELEASE                J. BARNES
C----------------------------------------------------------------------
C
      DIMENSION A(*)
      DIMENSION R(*)
C
      DO 10 I=0,NSAMPS-1
         R( IRINDX+(I*IRINC) ) = SCALAR * A( IAINDX+(I*IAINC) )
   10 CONTINUE
      RETURN
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPCNV(A, IAINDX, IAINC, B, IBINDX, IBINC,
     +                 R, IRINDX, IRINC, NSAMPS, NOPTS)
C
C     ROUTINE TO CONVOLVE NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT AN
C     INCREMENT IAINC WITH ARRAY B STARTING AT ELEMENT IBINDX AT AN
C     INCREMENT IBINC. RESULT RETURNED IN ARRAY R (STARTING AT ELEMENT
C     IRINDX AT AN INCREMENT IRINC).
C     REPLACES AP ROUTINE VADD
C
C     VERSION: JB01 26/02/85  INITIAL RELEASE                J. BARNES
C     VERSION: ID02 09/10/86   put in block if to speed it up.
C----------------------------------------------------------------------
C
      DIMENSION A(*)
      DIMENSION B(*)
      DIMENSION R(*)
C
      IF(IBINC.EQ.1.AND.IAINC.EQ.1)THEN
        DO 20 I=0,NSAMPS-1
          RSUM = 0.0
          DO 10 J=0,NOPTS-1
   10       RSUM = RSUM + A( IAINDX+I+J) * B( IBINDX+J )
   20     R( IRINDX+(I*IRINC) ) = RSUM
      ELSE IF(IBINC.EQ.-1.AND.IAINC.EQ.1)THEN
        DO 21 I=0,NSAMPS-1
          RSUM = 0.0
          DO 11 J=0,NOPTS-1
   11       RSUM = RSUM + A( IAINDX+I+J ) * B( IBINDX-J )
   21     R( IRINDX+(I*IRINC) ) = RSUM
      ELSE
        DO 22 I=0,NSAMPS-1
          RSUM = 0.0
          DO 12 J=0,NOPTS-1
   12       RSUM = RSUM + A(IAINDX+(I+J)*IAINC) * B( IBINDX+(J*IBINC) )
   22     R( IRINDX+(I*IRINC) ) = RSUM
      ENDIF
      RETURN
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPSVE(A, IAINDX, IAINC, R, NSAMPS)
C
C     ROUTINE TO SUM NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT
C     AN INCREMENT IAINC. THE RESULT IS RETURNED AS R.
C     REPLACES AP ROUTINE SVE
C
C     VERSION: YS01 31/12/85  INITIAL RELEASE                Y. SMITH
C----------------------------------------------------------------------
C
      DIMENSION A(*)
C
      R = 0.0
C
      DO 10 I=0,NSAMPS-1
        R = R + A(IAINDX+(I*IAINC))
   10 CONTINUE
      RETURN
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPMXMG(ARRAY,ISTART,IINC, MAX,NSAMPS)
C
C        ROUTINE TO RETURN THE MAX. VALUE IN AN ARRAY AND ITS INDEX
C        (IF NOT UNIQUE, RETURNS FIRST INDEX FOR MAX. VALUE)
C
      DIMENSION  ARRAY(*)
      REAL       MAX(*)
C
      MAX(1) = -1.0E+33
C
      DO 10 I=ISTART,NSAMPS,IINC
	 TEMP=ARRAY(I)
	 IF (TEMP.LT.0.0) THEN
	    TEMP=-TEMP
         ENDIF
         IF(TEMP.GT.MAX(1)) THEN
            MAX(1) = TEMP
            MAX(2) = I
         ENDIF
   10 CONTINUE
      RETURN
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c
c
c
      SUBROUTINE VPDIV(A, IAINDX, IAINC, B, IBINDX, IBINC,
     +                 R, IRINDX, IRINC, NSAMPS)
C
C     ROUTINE TO DIVIDE NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT
C     AN INCREMENT IAINC INTO ARRAY B STARTING AT ELEMENT IBINDX AT AN
C     INCREMENT IBINC. RESULT RETURNED IN ARRAY R (STARTING AT ELEMENT
C     IRINDX AT AN INCREMENT IRINC).
C     REPLACES AP ROUTINE VDIV
C
C     VERSION: YS01 24/12/85  INITIAL RELEASE                Y. SMITH
C----------------------------------------------------------------------
C
      DIMENSION A(*)
      DIMENSION B(*)
      DIMENSION R(*)
C
      DO 10 I=0,NSAMPS-1
        R(IRINDX+(I*IRINC)) = B(IBINDX+(I*IBINC)) / A(IAINDX+(I*IAINC))
   10 CONTINUE
      RETURN
      END






C-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
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
C               TEMP(NXC) REAL WORKING STORAGE.
C
C               G(NXC) : THE CROSS CORRELATION OF THE INPUT SIGNAL X
C               ,WITH THE OUTPUT Z,DIVIDED BY THE AUTOCORRELATION OF
C               Z EVALUATED AT ZERO TIME.
C
C               NXC: INTEGER, THE LENGTH OF THE CROSS CORRELATION
C               VECTOR.
C
C               A(NXC) : REAL, WORKING STORAGE.
C
C               C(NFIL) :REAL,THE OUTPUT FILTER COEFFICIENTS.
C
C               NTRY:INTEGER, THE NUMBER OF TIMES TO TRY
C               ADDING WHITE NOISE IF THE AUTOCORRELATION MATRIX IS
C               SINGULAR.
C
C               WPC: REAL, THE PERCENTAGE OF WHITE NOISE TO ADD ON
C               EACH TRY (TYPICALLY 0.5).
C
C               IXC1: INTEGER. THE CROSS CORRELATION VECTOR MUST
C               CONTAIN ALL THE VALUES, EVEN IF THEIR INDICES ARE
C               NEGATIVE. G(1) SHOULD CONTAIN THE FIRST NON ZERO
C               VALUE . IXC1 PASSES THE INDEX OF THE G VECTOR
C               THAT CONTAINS THE NON SHIFTED CROSS CORRELATION
C               VALUE ,I.E., AT T=0.
C
C               IFIL1: INTEGER, THE INDEX OF THE C VECTOR WHERE
C               THE TIME=0 VALUE OF THE FILTER IS TO BE DISPLACED TO.
C
C
C
      SUBROUTINE WIEN(R,TEMP,G,NXC,A,C,NTRY,WPC,NFIL,IXC1,IFIL1)
      REAL R(*),TEMP(*),G(*),A(*),C(*)
      REAL K,L


      ITRY=00
      RMAX=.72E33
      IOFF=IFIL1-IXC1

  10  ITRY=ITRY+1
      IF (ITRY.GT.NTRY) GOTO 1000
      IF (ITRY.GT.1) THEN
          IERR = ITRY
      END IF


C     ** Set recursion start values.
C
      A(1)=1.0
      C(IOFF+1)=G(1)/R(1)
      ALPHA=R(1)*(1.0+ITRY*WPC/100.)
      BETA=R(2)
      GAMMA=C(IOFF+1)*R(2)
      IF (C(IOFF+1).GT.RMAX) GOTO 10


C     **   Main loop
C
      DO 600 MP1=2,NXC-1
          K=-BETA/ALPHA
          IF (K.GT.RMAX) GOTO 10
          M=MP1-1

C         ** Part 2 (from mesko)
C
          DO 20 I=2,M
             TEMP(I)=A(I)+K*A(M-I+2)
   20     CONTINUE
          DO 25 I=2,M
             A(I)=TEMP(I)
   25     CONTINUE
          A(MP1)=K*A(1)


C         **   Part 3
C
          ALPHA=ALPHA+K*BETA
          BETA=0.0
          DO 30 I=1,MP1
             BETA=BETA+A(I)*R(M+3-I)
   30     CONTINUE


C         **  Part 4
C
C
          L=(G(MP1)-GAMMA)/ALPHA

C         **  Part 5
C
          DO 40 I=1,M
             C(IOFF+I)=C(IOFF+I)+L*A(M-I+2)
   40     CONTINUE
          C(MP1+IOFF)=L*A(1)


C         **  Part 6
C
          GAMMA=0.0
          DO 50 I=1,MP1
             GAMMA=GAMMA+C(I+IOFF)*R(M-I+3)
   50     CONTINUE

  600 CONTINUE




 1000 RETURN
      END


c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPMOV(A, IAINDX, IAINC, R, IRINDX, IRINC, NSAMPS)
C
C     ROUTINE TO MOVE NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT AN
C     INCREMENT IAINC TO ARRAY R STARTING AT ELEMENT IRINDX AT AN
C     INCREMENT IRINC.
C     REPLACES AP ROUTINE VMOV
C
C     VERSION: JB01 26/02/85  INITIAL RELEASE                J. BARNES
C----------------------------------------------------------------------
C
      DIMENSION A(*)
      DIMENSION R(*)
C
      DO 10 I=0,NSAMPS-1
         R( IRINDX+(I*IRINC) ) = A( IAINDX+(I*IAINC) )
   10 CONTINUE
      RETURN
      END



c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPABS(A, IAINDX, IAINC,
     +                  R, IRINDX, IRINC, NSAMPS)
C
C     ROUTINE TO TAKE THE ABSOLUTE VALUE OF NSAMPS OF ARRAY A FROM
C     ELEMENT IAINDX AT INCREMENT IAINC. RESULT RETURNED IN ARRAY R
C     (STARTING AT ELEMENT IRINDX AT AN INCREMENT IRINC).
C     REPLACES AP ROUTINE VABS
C
C     VERSION: YS01 31/12/85  INITIAL RELEASE                Y. SMITH
C----------------------------------------------------------------------
C
      DIMENSION A(*)
      DIMENSION R(*)
C
      DO 10 I=0,NSAMPS-1
        R(IRINDX+(I*IRINC)) = ABS ( A(IAINDX+(I*IAINC)) )
   10 CONTINUE
      RETURN
      END

c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-c-
c
c
      SUBROUTINE VPSUB(A, IAINDX, IAINC, B, IBINDX, IBINC,
     +                 R, IRINDX, IRINC, NSAMPS)
C
C     ROUTINE TO SUBTRACT NSAMPS OF ARRAY A FROM ELEMENT IAINDX AT
C     AN INCREMENT IAINC FROM ARRAY B STARTING AT ELEMENT IBINDX AT AN
C     INCREMENT IBINC. RESULT RETURNED IN ARRAY R (STARTING AT ELEMENT
C     IRINDX AT AN INCREMENT IRINC).
C     REPLACES AP ROUTINE VSUB
C
C     VERSION: YS01 31/12/85  INITIAL RELEASE                Y. SMITH
C----------------------------------------------------------------------
C
      DIMENSION A(*)
      DIMENSION B(*)
      DIMENSION R(*)
C
      DO 10 I=0,NSAMPS-1
        R(IRINDX+(I*IRINC)) = B(IBINDX+(I*IBINC)) - A(IAINDX+(I*IAINC))
   10 CONTINUE
      RETURN
      END
