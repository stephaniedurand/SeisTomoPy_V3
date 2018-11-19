ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
       subroutine sph2geo2s40rts(NS,p,np,radial_basis,
     +            filename,filename2,filename3,NN,model,rad,r)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      Transforme les coefficients d'harmoniques spheriques en coordonnees geographiques


       implicit NONE
       integer model,NN,np,ind,mm,ll,NS,k,kk,radius,npp,ligne(13),nb,l
       integer ttheta,pphi,rad,depth,s
       real(kind=8) radial_basis(model,0:NN-1),yraw(1680)
       real(kind=8) phi,theta,cost,sint,splitf,splitfd
       real(kind=8) x(2*NS+1),dx(2*NS+1)
       real(kind=8) p(np),dg2rad,tpi,dp(np)
       real(kind=8) deg,r(model),rsple,Ralpha,Rrho
       real*8 q(3,model),f(3,model),func(0:NN-1),func3(model)
       character*(*) filename,filename2,filename3
       external rsple
       real(kind=8) dvvst(3721)
       character(len=255) homedir,tampon

       call getenv("HOME", tampon)
       homedir=trim(adjustl(tampon))//'/SeisTomoPy_files/'

       tpi=8.d0*datan(1.d0)
       dg2rad=tpi/360.d0

!      Comme dans S40RTS
       Ralpha=((2.d0-3.d0)/(6371.d0-3480.d0))
     + *(6371.d0-dfloat(rad))+
     + 2.d0-((2.d0-3.d0)/(6371.d0-3480.d0))*6371.d0
       Ralpha=1.d0/Ralpha
       Rrho=30.d-2

       npp=0
       do k=1,NS
       npp=(2*k+1)+npp
       enddo

       open(1,file='../models/S40RTS_1km.sph')
       do depth=2890,0,-1
       if (depth.eq.rad) then
       read(1,'(3721f15.8)') (dvvst(s),s=1,3721)
!       write(*,'(1681f15.8)') (dvvst(s,cd),s=1,1681)
       else
       read(1,'(3721f15.8)')
       endif
       enddo
       close(1)

       open(13,file=trim(adjustl(homedir))//
     + 'output_files_map/output_map_S40RTS.out')
       open(10,file=trim(adjustl(filename)))
       open(11,file=trim(adjustl(filename2)))
       open(12,file=trim(adjustl(filename3)))
       do pphi=0,359
       do ttheta=-89,89 !89
       theta=(90.d0-dfloat(ttheta))*dg2rad
       phi=dfloat(pphi)*dg2rad

       ind=1
       splitf=0.d0
       do ll=0,NS

       sint=dsin(theta)
       cost=dcos(theta)
       call lgndr(ll,cost,sint,x,dx)

       splitf=splitf+1.d0*x(1)*dvvst(ind)
       ind=ind+1

       if (ll.ne.0) then
       do mm=1,ll
       splitf=splitf
     +            +(2.d0)*dvvst(ind)*(dcos(dble(mm)*phi))*x(mm+1)
       ind=ind+1
       splitf=splitf
     +           +(2.d0)*dvvst(ind)*(dsin(dble(mm)*phi))*x(mm+1)
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

       end subroutine sph2geo2s40rts

