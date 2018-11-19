
!----------------------------------

! changed the obsolecent f77 features in the two routines below
! now still awful Fortran, but at least conforms to f90 standard

       real(kind=8) function rsple(I1,I2,X,Y,Q,S)
       implicit none

! rsple returns the value of the function y(x) evaluated at point S
! using the cubic spline coefficients computed by rspln and saved in Q.
! If S is outside the interval (x(i1),x(i2)) rsple extrapolates
! using the first or last interpolation polynomial. The arrays must
! be dimensioned at least - x(i2), y(i2), and q(3,i2).

       integer i1,i2
       real(kind=8) X(*),Y(*),Q(3,*),s
       integer i,ii
       real(kind=8) h

       i = 1
       II=I2-1

! GUARANTEE I WITHIN BOUNDS.
       I=MAX0(I,I1)
       I=MIN0(I,II)

! SEE IF X IS INCREASING OR DECREASING.
       IF(X(I2)-X(I1) < 0) goto 1
       IF(X(I2)-X(I1) >= 0) goto 2

! X IS DECREASING. CHANGE I AS NECESSARY.
 1         IF(S-X(I) <= 0) goto 3
       IF(S-X(I) > 0) goto 4

 4     I=I-1

       IF(I-I1 < 0) goto 11
       IF(I-I1 == 0) goto 6
       IF(I-I1 > 0) goto 1

 3     IF(S-X(I+1) < 0) goto 5
       IF(S-X(I+1) >= 0) goto 6

 5     I=I+1

      IF(I-II < 0) goto 3
      IF(I-II == 0) goto 6
      IF(I-II > 0) goto 7

! X IS INCREASING. CHANGE I AS NECESSARY.
 2     IF(S-X(I+1) <= 0) goto 8
       IF(S-X(I+1) > 0) goto 9

 9     I=I+1

       IF(I-II < 0) goto 2
       IF(I-II == 0) goto 6
       IF(I-II > 0) goto 7

 8     IF(S-X(I) < 0) goto 10
       IF(S-X(I) >= 0) goto 6

 10    I=I-1
       IF(I-I1 < 0) goto 11
       IF(I-I1 == 0) goto 6
       IF(I-I1 > 0) goto 8

 7     I=II
       GOTO 6
 11    I=I1

! CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
 6     H=S-X(I)
       RSPLE=Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))

       end function rsple

!----------------------------------

       subroutine rspln(I1,I2,X,Y,Q,F)
       implicit none

! Subroutine rspln computes cubic spline interpolation coefficients
! for y(x) between grid points i1 and i2 saving them in q.The
! interpolation is continuous with continuous first and second
! derivatives. It agrees exactly with y at grid points and with the
! three point first derivatives at both end points (i1 and i2).
! X must be monotonic but if two successive values of x are equal
! a discontinuity is assumed and separate interpolation is done on
! each strictly monotonic segment. The arrays must be dimensioned at
! least - x(i2), y(i2), q(3,i2), and f(3,i2).
! F is working storage for rspln.

       integer i1,i2
       real(kind=8) X(*),Y(*),Q(3,*),F(3,*)
       integer i,j,k,j1,j2
       real(kind=8) y0,a0,b0,b1,h,h2,ha,h2a,h3a,h2b
       real(kind=8) YY(3),small
       equivalence (YY(1),Y0)
       data SMALL/1.0d-08/,YY/0.0d0,0.0d0,0.0d0/

       J1=I1+1
       Y0=0.0d0

! BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL
       IF(I2-I1 < 0) return
       IF(I2-I1 == 0) goto 17
       IF(I2-I1 > 0) goto 8

 8     A0=X(J1-1)
! SEARCH FOR DISCONTINUITIES.
       DO 3 I=J1,I2
       B0=A0
       A0=X(I)
       IF(DABS((A0-B0)/DMAX1(A0,B0)).LT.SMALL) GOTO 4
 3     CONTINUE
 17    J1=J1-1
       J2=I2-2
       GOTO 5
 4     J1=J1-1
       J2=I-3
! SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
 5     IF(J2+1-J1 < 0) goto 9
       IF(J2+1-J1 == 0) goto 10
       IF(J2+1-J1 > 0) goto 11

! ONLY TWO POINTS. USE LINEAR INTERPOLATION.
 10    J2=J2+2
       Y0=(Y(J2)-Y(J1))/(X(J2)-X(J1))
       DO J=1,3
       Q(J,J1)=YY(J)
       Q(J,J2)=YY(J)
       enddo
       GOTO 12

! MORE THAN TWO POINTS. DO SPLINE INTERPOLATION.
 11    A0=0.
       H=X(J1+1)-X(J1)
       H2=X(J1+2)-X(J1)
       Y0=H*H2*(H2-H)
       H=H*H
       H2=H2*H2
! CALCULATE DERIVITIVE AT NEAR END.
       B0=(Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
       B1=B0

! EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
       DO I=J1,J2
       H=X(I+1)-X(I)
       Y0=Y(I+1)-Y(I)
       H2=H*H
       HA=H-A0
       H2A=H-2.0d0*A0
       H3A=2.0d0*H-3.0d0*A0
       H2B=H2*B0
       Q(1,I)=H2/HA
       Q(2,I)=-HA/(H2A*H2)
       Q(3,I)=-H*H2A/H3A
       F(1,I)=(Y0-H*B0)/(H*HA)
       F(2,I)=(H2B-Y0*(2.0d0*H-A0))/(H*H2*H2A)
       F(3,I)=-(H2B-3.0d0*Y0*HA)/(H*H3A)
       A0=Q(3,I)
       B0=F(3,I)
       enddo

! TAKE CARE OF LAST TWO ROWS.
       I=J2+1
       H=X(I+1)-X(I)
       Y0=Y(I+1)-Y(I)
       H2=H*H
       HA=H-A0
       H2A=H*HA
       H2B=H2*B0-Y0*(2.0d0*H-A0)
       Q(1,I)=H2/HA
       F(1,I)=(Y0-H*B0)/H2A
       HA=X(J2)-X(I+1)
       Y0=-H*HA*(HA+H)
       HA=HA*HA

! CALCULATE DERIVATIVE AT FAR END.
       Y0=(Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
       Q(3,I)=(Y0*H2A+H2B)/(H*H2*(H-2.0d0*A0))
       Q(2,I)=F(1,I)-Q(1,I)*Q(3,I)

! SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
       DO J=J1,J2
       K=I-1
       Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
       Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
       Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
       I=K
       enddo
       Q(1,I)=B1
! FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
 9     J2=J2+2
       DO J=1,3
       Q(J,J2)=YY(J)
       enddo

! SEE IF THIS DISCONTINUITY IS THE LAST.
 12    IF(J2-I2 < 0) then
       goto 6
       else
       return
       endif

! NO. GO BACK FOR MORE.
 6     J1=J2+2
       IF(J1-I2 <= 0) goto 8
       IF(J1-I2 > 0) goto 7

! THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
 7     DO J=1,3
       Q(J,I2)=YY(J)
       enddo

       end subroutine rspln

