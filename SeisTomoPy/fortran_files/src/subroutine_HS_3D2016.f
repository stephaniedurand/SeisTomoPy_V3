ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
       subroutine sph2geo23D2016(NS,p,np,radial_basis,
     +            filename,filename2,filename3,pnts,NN,model,
     +           dvvst,cd,dep,width)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Transforme les coefficients d'harmoniques spheriques en coordonnees geographiques


       implicit NONE
       integer model,NN,np,ind,mm,ll,NS,k,kk,radius,npp,ligne(13),nb
       real(kind=8) radial_basis(model,0:NN-1),yraw(1681)
       real(kind=8) phi,theta,cost,sint,splitf,splitfd
       real(kind=8) x(2*NS+1),dx(2*NS+1)
       real(kind=8) p(np),dg2rad,tpi,dp(np)
       real(kind=8) Ralpha,Rrho 
       real(kind=8) ttheta,pphi,deg,pnts(3,10000)
       character*(*) filename,filename2,filename3
       integer depth,s,cd,dep,width
       real(kind=8) dvvst(3721,579)
       character(len=255) homedir,tampon

       call getenv("HOME", tampon)
       homedir=trim(adjustl(tampon))//'/SeisTomoPy_files/'
       
       tpi=8.d0*datan(1.d0)
       dg2rad=tpi/360.d0

       npp=0
       do k=1,NS
       npp=(2*k+1)+npp
       enddo 

c       open(1,file='../models/S40RTS_5km.sph')
c       cd=1
c       do depth=2890,0,-5
c
c!      Comme dans S40RTS
c       Ralpha(cd)=((2.d0-3.d0)/(6371.d0-3480.d0))
c     + *(6371.d0-dfloat(depth))+
c     + 2.d0-((2.d0-3.d0)/(6371.d0-3480.d0))*6371.d0
c       Ralpha(cd)=1.d0/Ralpha(cd)
c       Rrho=50.d-2
c
c       read(1,'(3721f15.8)') (dvvst(s,cd),s=1,3721) 
c!       write(*,'(3721f15.8)') (dvvst(s,cd),s=1,3721) 
c       cd=cd+1
c       enddo
c       close(1)

       open(10,file=trim(adjustl(filename)))
       open(11,file=trim(adjustl(filename2)))
       open(12,file=trim(adjustl(filename3)))
       open(13,file=trim(adjustl(homedir))//
     + 'output_files_cross/output_cross_3D2016.out')
       open(14,file=trim(adjustl(homedir))//
     + 'output_files_cross/3D2016_input_AxiSEM.sph')
       open(15,file=trim(adjustl(homedir))//
     + 'output_files_cross/3D2016_PathPy.sph')
       write(14,*) (width+1)*cd
       do kk=1,width+1
       deg=(180.d0-(90.d0-pnts(3,kk)))
       ttheta=pnts(2,kk)
       pphi=pnts(1,kk)
       theta=(90.d0-ttheta)*dg2rad
       phi=pphi*dg2rad

       radius=1
       do depth=2890,0,-5

       ind=1
       splitf=0.d0
       do ll=0,NS
c       write(*,*) ll,ind
       sint=dsin(theta)
       cost=dcos(theta)
       call lgndr(ll,cost,sint,x,dx)

       splitf=splitf+1.d0*x(1)*dvvst(ind,radius)
       ind=ind+1

       if (ll.ne.0) then
       do mm=1,ll
       splitf=splitf
     +            +(2.d0)*dvvst(ind,radius)*(dcos(dble(mm)*phi))*x(mm+1)
       ind=ind+1
       splitf=splitf
     +           +(2.d0)*dvvst(ind,radius)*(dsin(dble(mm)*phi))*x(mm+1)
       ind=ind+1
c       write(*,*) ll,ind-1
       enddo
       endif
       enddo
       Ralpha = 0.d0
       Rrho = 0.d0
       write(10,*) deg,6371-depth,Ralpha*splitf*100.d0
       write(11,*) deg,6371-depth,splitf*100.d0
       write(12,*) deg,6371-depth,Rrho*splitf*100.d0
       write(13,*) ttheta,pphi,6371-depth,Ralpha*splitf*100.d0,
     + splitf*100.d0,Rrho*splitf*100.d0
       write(14,*) 6371-depth,deg,Ralpha*splitf*100.d0,
     + splitf*100.d0,Rrho*splitf*100.d0
       write(15,*) 6371-depth,kk-1,Ralpha*splitf*100.d0,
     + splitf*100.d0,Rrho*splitf*100.d0

       radius=radius+1
       enddo
       enddo
       close(14)
       close(13)
       close(12)
       close(11)
       close(10)
       close(15)

       end subroutine sph2geo23D2016

