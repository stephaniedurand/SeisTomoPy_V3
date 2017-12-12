ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
       subroutine sph2geo22(NS,p,np,radial_basis,
     +            filename,filename2,filename3,NN,model,rad,r)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Transforme les coefficients d'harmoniques spheriques en coordonnees geographiques


       implicit NONE
       integer model,NN,np,ind,mm,ll,NS,k,kk,radius,npp,ligne(13),nb,l
       integer ttheta,pphi
       real(kind=8) radial_basis(model,0:NN-1),yraw(3721)
       real(kind=8) phi,theta,cost,sint,splitf,splitfd
       real(kind=8) x(2*NS+1),dx(2*NS+1)
       real(kind=8) p(np),dg2rad,tpi,dp(np)
       real(kind=8) deg,rad,r(model),rsple,Ralpha,Rrho
       real*8 q(3,model),f(3,model),func(0:NN-1),func3(model)
       character*(*) filename,filename2,filename3
       external rsple

!      Definition de la loi de sclaing Ralpha
       Ralpha=55.d-2
       Rrho=20.d-2

       tpi=8.d0*datan(1.d0)
       dg2rad=tpi/360.d0

       npp=0
       do k=1,NS
       npp=(2*k+1)+npp
       enddo 

       do k=0,NN-1
        do l=1,model
        func3(l)=radial_basis(l,k)
        enddo
        call rspln(1,model,r,func3,q,f)
        func(k)=rsple(1,model,r,func3,q,rad)
       enddo

       do mm=1,npp
       yraw(mm)=0.d0
       do k=0,NN-1
       yraw(mm)=yraw(mm)
     +            +func(k)*p((mm-1)*NN+k+1)
       enddo
       enddo
       
       open(13,file='../output_files_map/output_map_SEISGLOB2.out')
       open(10,file=trim(adjustl(filename)))
       open(11,file=trim(adjustl(filename2)))
       open(12,file=trim(adjustl(filename3)))
       do pphi=0,359
       do ttheta=-89,89 !89
       theta=(90.d0-dfloat(ttheta))*dg2rad
       phi=dfloat(pphi)*dg2rad


       ind=1
       splitf=0.d0
       do ll=1,NS

       sint=dsin(theta)
       cost=dcos(theta)
       call lgndr(ll,cost,sint,x,dx)

       splitf=splitf+1.d0*x(1)*yraw(ind)
       ind=ind+1

       if (ll.ne.0) then
       do mm=1,ll
       splitf=splitf
     +            +(2.d0)*yraw(ind)*(dcos(dble(mm)*phi))*x(mm+1)
       ind=ind+1
       splitf=splitf
     +           +(2.d0)*yraw(ind)*(dsin(dble(mm)*phi))*x(mm+1)
       ind=ind+1
       enddo
       endif
       enddo
       write(10,*) ttheta,pphi,Ralpha*splitf*100.d0
       write(11,*) ttheta,pphi,splitf*100.d0
       write(12,*) ttheta,pphi,Rrho*splitf*100.d0
       write(13,*) ttheta,pphi,Ralpha*splitf*100.d0,
     + splitf*100.d0,Rrho*splitf*100.d0
       enddo
       enddo
       close(10)
       close(11)
       close(12)
       close(13)        

       end subroutine sph2geo22

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
       subroutine sph2geo21(NS,p,np,radial_basis,
     +            filename,filename2,filename3,NN,model,rad,r)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Transforme les coefficients d'harmoniques spheriques en coordonnees geographiques


       implicit NONE
       integer model,NN,np,ind,mm,ll,NS,k,kk,radius,npp,ligne(13),nb,l
       integer ttheta,pphi
       real(kind=8) radial_basis(model,0:NN-1),yraw(3721)
       real(kind=8) phi,theta,cost,sint,splitf,splitfd
       real(kind=8) x(2*NS+1),dx(2*NS+1)
       real(kind=8) p(np),dg2rad,tpi,dp(np)
       real(kind=8) deg,rad,r(model),rsple,Ralpha,Rrho
       real*8 q(3,model),f(3,model),func(0:NN-1),func3(model)
       character*(*) filename,filename2,filename3
       external rsple

!      Definition de la loi de sclaing Ralpha
       Ralpha=55.d-2
       Rrho=20.d-2

       tpi=8.d0*datan(1.d0)
       dg2rad=tpi/360.d0

       npp=0
       do k=1,NS
       npp=(2*k+1)+npp
       enddo 

       do k=0,NN-1
        do l=1,model
        func3(l)=radial_basis(l,k)
        enddo
        call rspln(1,model,r,func3,q,f)
        func(k)=rsple(1,model,r,func3,q,rad)
       enddo

       do mm=1,npp
       yraw(mm)=0.d0
       do k=0,NN-1
       yraw(mm)=yraw(mm)
     +            +func(k)*p((mm-1)*NN+k+1)
       enddo
       enddo
       
       open(13,file='../output_files_map/output_map_SEISGLOB1.out')
       open(10,file=trim(adjustl(filename)))
       open(11,file=trim(adjustl(filename2)))
       open(12,file=trim(adjustl(filename3)))
       do pphi=0,359
       do ttheta=-89,89 !89
       theta=(90.d0-dfloat(ttheta))*dg2rad
       phi=dfloat(pphi)*dg2rad


       ind=1
       splitf=0.d0
       do ll=1,NS

       sint=dsin(theta)
       cost=dcos(theta)
       call lgndr(ll,cost,sint,x,dx)

       splitf=splitf+1.d0*x(1)*yraw(ind)
       ind=ind+1

       if (ll.ne.0) then
       do mm=1,ll
       splitf=splitf
     +            +(2.d0)*yraw(ind)*(dcos(dble(mm)*phi))*x(mm+1)
       ind=ind+1
       splitf=splitf
     +           +(2.d0)*yraw(ind)*(dsin(dble(mm)*phi))*x(mm+1)
       ind=ind+1
       enddo
       endif
       enddo
       write(10,*) ttheta,pphi,Ralpha*splitf*100.d0
       write(11,*) ttheta,pphi,splitf*100.d0
       write(12,*) ttheta,pphi,Rrho*splitf*100.d0
       write(13,*) ttheta,pphi,Ralpha*splitf*100.d0,
     + splitf*100.d0,Rrho*splitf*100.d0
       enddo
       enddo
       close(10)
       close(11)
       close(12)
       close(13)        

       end subroutine sph2geo21

!--------------------------------------------------------

!--------------------------------------------------------
        subroutine lgndr(l,c,s,x,dx)

        ! computes Legendre function x(l,m,theta)
        ! theta=colatitude,c=cos(theta),s=sin(theta),l=angular order,
        ! sin(theta) restricted so that sin(theta) > 1.e-7
        ! x(1) contains m=0, x(2) contains m=1, x(k+1) contains m=k
        ! m=azimuthal (longitudinal) order 0 <= m <= l
        ! dx=dx/dtheta
        !
        ! subroutine originally came from Physics Dept. Princeton through
        ! Peter Davis, modified by Jeffrey Park

        implicit none

        ! argument variables
        integer l
        double precision x(2*l+1),dx(2*l+1)
        double precision c,s

        ! local variables
        integer i,lp1,lpsafe,lsave
        integer m,maxsin,mmm,mp1

        double precision sqroot2over2,c1,c2,cot
        double precision ct,d,f1,f2
        double precision f3,fac,g1,g2
        double precision g3,rfpi,sqroot3,sos
        double precision ss,stom,t,tol
        double precision v,y

        tol = 1.d-05
        rfpi = 0.282094791773880d0
        sqroot3 = 1.73205080756890d0
        sqroot2over2 = 0.707106781186550d0

        if(s >= 1.0d0-tol) s=1.0d0-tol
        lsave=l
        if(l<0) l=-1-l
        if(l>0) goto 1
        x(1)=rfpi
        dx(1)=0.0d0
        l=lsave
        return
 1      if (l /= 1) goto 2
        c1=sqroot3*rfpi
        c2=sqroot2over2*c1
        x(1)=c1*c
        x(2)=-c2*s
        dx(1)=-c1*s
        dx(2)=-c2*c
        l=lsave
        return
 2          sos=s
        if (s<tol) s=tol
        cot=c/s
        ct=2.0d0*c
        ss=s*s
        lp1=l+1
        g3=0.0d0
        g2=1.0d0
        f3=0.0d0

! evaluate m=l value, sans (sin(theta))**l
        do i=1,l
          g2=g2*(1.0d0-1.0d0/(2.0d0*i))
        enddo
        g2=rfpi*dsqrt((2*l+1)*g2)
        f2=l*cot*g2
        x(lp1)=g2
        dx(lp1)=f2
        v=1.0d0
        y=2.0d0*l
        d=dsqrt(v*y)
        t=0.0d0
        mp1=l
        m=l-1

! these recursions are similar to ordinary m-recursions, but since we
! have taken the s**m factor out of the xlm's, the recursion has the powers
! of sin(theta) instead
 3      g1=-(ct*mp1*g2+ss*t*g3)/d
        f1=(mp1*(2.0d0*s*g2-ct*f2)-t*ss*(f3+cot*g3))/d-cot*g1
        x(mp1)=g1
        dx(mp1)=f1
        if (m == 0) goto 4
        mp1=m
        m=m-1
        v=v+1.0d0
        y=y-1.0d0
        t=d
        d=dsqrt(v*y)
        g3=g2
        g2=g1
        f3=f2
        f2=f1
        goto 3
! explicit conversion to integer added
 4      maxsin=int(-72.0d0/log10(s))

! maxsin is the max exponent of sin(theta) without underflow
        lpsafe=min0(lp1,maxsin)
        stom=1.0d0
        fac=sign(1.0d0,dble((l/2)*2-l) + 0.50d0)

! multiply xlm by sin**m
        do m=1,lpsafe
        x(m)=fac*x(m)*stom
        dx(m)=fac*dx(m)*stom
        stom=stom*s
        enddo

! set any remaining xlm to zero
        if (maxsin <= l) then
        mmm=maxsin+1
        do m=mmm,lp1
        x(m)=0.0d0
        dx(m)=0.0d0
        enddo
        endif

        s=sos
        l=lsave

        end subroutine lgndr         


