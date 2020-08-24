!
! File:   wfdecomp.F95
! Author: yes
!
! Created on November 30, 2009, 1:38 PM
!
subroutine wfs(nf, df, fmin, fmax, IsMatched, IsV, Vp, Vs, Imp, px, py, nc, trdata, nm, modes, trout, HZXY)
   
    implicit none
    integer, parameter :: SP = 4, MaxModes = 12      ! use only single precision
    integer(SP),  intent(in)    :: nf, nc, nm        ! number of sample per trace, components input and output
    integer(SP),  intent(in)    :: IsMatched         ! if traces are matched
    integer(SP),  intent(in)    :: IsV               ! geophones are velocity sensor instead accelerometer
    integer(SP),  intent(in)    :: modes(nm)         ! what kinds of separated waves output
    real(SP),     intent(in)    :: df, fmin, fmax    ! sampling interval, min and max freq
    real(SP),     intent(in)    :: Vp, Vs, Imp
    real(SP),     intent(in)    :: px, py
    complex(SP),  intent(inout) :: trdata(nf, nc)    ! input traces of nc components
    complex(SP),  intent(out)   :: trout(nf, nm)     ! output traces of nm components
    complex(SP),  intent(out)   :: HZXY(4, MaxModes) ! coefficient matrix

    integer(SP) :: m
    real(SP)    :: f
    real(SP)    :: a, b, c, d, p
    real(SP)    :: COStheta2
    real(SP)    :: minCOStheta2 = 0.01
    complex(SP) :: Fa
    complex(SP) :: VpQp, divVpQp, VsQs, divVsQs
    complex(SP) :: scale          ! 1/iw ( 1/(j*2*PI*F) ) applied to accelerometer (convert to velocity)

    real(SP), parameter  :: PI = 3.14159

    p = sqrt(px*px + py*py)

!   Cap the minimum of V*Q = cos(theta) thus division not go to unlimited
    COStheta2 = 1.0 - p*p*Vp*Vp
    if (COStheta2 <= -minCOStheta2  ) then
        VpQp = -cmplx(0, sqrt(-COStheta2))
        divVpQp = 1.0/VpQp
    else if(COStheta2 <= 0.0) then
        VpQp = -cmplx(0, sqrt(-COStheta2))
        divVpQp = 1.0/cmplx(0, sqrt(minCOStheta2))
    else if (COStheta2 <= minCOStheta2) then
        VpQp = cmplx(sqrt(COStheta2), 0)
        divVpQp = 1.0/cmplx(sqrt(minCOStheta2), 0)
    else
        VpQp = cmplx(sqrt(COStheta2), 0)
        divVpQp = 1.0/VpQp
    endif
    COStheta2 = 1.0 - p*p*Vs*Vs
    if (COStheta2 <= -minCOStheta2  ) then
        VsQs = -cmplx(0, sqrt(-COStheta2))
        divVsQs = 1.0/VsQs
    else if(COStheta2 <= 0.0) then
        VsQs = -cmplx(0, sqrt(-COStheta2))
        divVsQs = 1.0/cmplx(0, sqrt(minCOStheta2))
    else if (COStheta2 <= minCOStheta2) then
        VsQs = cmplx(sqrt(COStheta2), 0)
        divVsQs = 1.0/cmplx(sqrt(minCOStheta2), 0)
    else
        VsQs = cmplx(sqrt(COStheta2), 0)
        divVsQs = 1.0/VsQs
    endif

    a = 1.0 - 2.0*p*p*Vp*Vp
    b = 1.0 - 2.0*p*p*Vs*Vs

    if ( p == 0.0 ) then
        c = 0.0
        d = 0.0
    else
        c = px/p
        d = py/p
    endif

    Fa = 4.0*p*p*Vs*Vs*Vs*VsQs*VpQp/Vp + b*b

    HZXY(3,1) = -px*a*Vs*Vs/Vp/Imp
    HZXY(4,1) = -py*a*Vs*Vs/Vp/Imp
    HZXY(2,1) =  0.5*a*b*divVpQp/Imp
    HZXY(1,1) = -0.5*a/Imp
    HZXY(3,2) = px*Vs*Vs/Vp/Imp
    HZXY(4,2) = py*Vs*Vs/Vp/Imp
    HZXY(2,2) = 0.5*b*divVpQp/Imp
    HZXY(1,2) = 0.5/Imp
    HZXY(3,3) = -0.5*c*b*b*Vp*divVpQp/Vs/Imp
    HZXY(4,3) = -0.5*d*b*b*Vp*divVpQp/Vs/Imp
    HZXY(2,3) = -0.5*p*Vs*b/Imp
    HZXY(1,3) = -0.5*p*Vp*b*divVsQs/Imp
    HZXY(3,4) = -0.5*c*b*Vp*divVpQp/Vs/Imp
    HZXY(4,4) = -0.5*d*b*Vp*divVpQp/Vs/Imp
    HZXY(2,4) = 0.5*p*Vs/Imp
    HZXY(1,4) = 0.5*p*Vp*divVsQs/Imp
    HZXY(3,5) = px*Vs*Vs*Vs*VpQp*divVsQs/Vp/Vp
    HZXY(4,5) = py*Vs*Vs*Vs*VpQp*divVsQs/Vp/Vp
    HZXY(2,5) = 0.5*b*divVpQp
    HZXY(1,5) = (2.0*p*p*Vs*Vs*Vs*Vs*VpQp*VpQp/Vp/Vp - 0.5*b*(p*p*Vs*Vs*VpQp/Vp + VsQs/Vs)*Vs*divVsQs)/Fa
    HZXY(3,6) =  HZXY(3,5)
    HZXY(4,6) =  HZXY(4,5)
    HZXY(2,6) = -HZXY(2,5)
    HZXY(1,6) =  HZXY(1,5)
    HZXY(3,7) = -c*b*Vs*divVsQs/Vp
    HZXY(4,7) = -d*b*Vs*divVsQs/Vp
    HZXY(2,7) = p*Vs*Vs/Vp
    HZXY(1,7) = -p*(2.0*b*Vs*Vs*VpQp/Vp - 0.5*b*b*Vs*divVsQs)/Fa
    HZXY(3,8) = -HZXY(3,7)
    HZXY(4,8) = -HZXY(4,7)
    HZXY(2,8) =  HZXY(2,7)
    HZXY(1,8) = -HZXY(1,7)
    HZXY(3,9) = -px*Vs*Vs*VpQp/Vp/Imp
    HZXY(4,9) = -py*Vs*Vs*VpQp/Vp/Imp
    HZXY(2,9) = 0.5*b/Imp
    HZXY(1,9) = -0.5*VpQp/Imp
    HZXY(3,10) = -HZXY(3,9)
    HZXY(4,10) = -HZXY(4,9)
    HZXY(2,10) =  HZXY(2,9)
    HZXY(1,10) = -HZXY(1,9)
    HZXY(3,11) = 0.5*c*b/Imp
    HZXY(4,11) = 0.5*d*b/Imp
    HZXY(2,11) = p*Vs*VsQs/Imp
    HZXY(1,11) = -0.5*p*Vp/Imp
    HZXY(3,12) =  HZXY(3,11)
    HZXY(4,12) =  HZXY(4,11)
    HZXY(2,12) = -HZXY(2,11)
    HZXY(1,12) =  HZXY(1,11)

    if ( IsMatched == 0 )  HZXY(2:4, :) = HZXY(2:4, :)*Imp
    
    if (IsV == 0 ) then ! accelerometer geophone, V=A/(img*Omiga) = A/(img*2*PI*f)
        do m=0, nf
            f=m*df
            if (f <= fmin .or. f > fmax ) then
                scale = cmplx(0.0, 0.0)
            else
                scale = 1.0 / cmplx(0, 1.0/(2.0*PI*f))
            endif
            trdata(m, 1:nc) = trdata(m, 1:nc) * scale
        enddo
    endif

    do m=1, nm
        trout(:,m) = MATMUL(trdata(:, 1:nc), HZXY(1:nc, modes(m)))
    enddo

end subroutine wfs

