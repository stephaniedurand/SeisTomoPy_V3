!=====================================================================
!                          PROFILE-MAKER PROGRAM
!       S. Durand, May 2016, V0.0 
!=====================================================================
       implicit none
!    Nombres de degres inverse pour chaque donnees
       integer c,NSmax,np,model,NN,spline,k,l,ddeg,kk,cclat,cclon,aaz,
     + np1,NSmax1,lllmax
       integer s,cd,cY,depth,cd1,dep,width,degstart,tomo
       parameter (np=35280,model=579,NN=21,spline=1,np1=9240,
     + NSmax1=20)
!   earth model
       real(kind=8) r(model),rho(model),kappa(model),mu(model),
     +        rhos2(model),qkap(model),qmu(model),pg(model)
       real(kind=8) radial_basis(model,0:NN-1)
       real(kind=8) p(np1),xx,yy,dp(np)
       real(kind=8) p1(np),dp1(np1)
       real(kind=8) clat,clon,az,az2,elat,elon,slat,slon,dist
       real(kind=8) pnts(3,10000),az12,az21,Gc,blat,blon,plat(55),
     + plon(55),step,eqlon(4995),eqlat(4995),eqdep(4995)
       integer nrp,igrc
!   Nom fichiers
       character*(80) filename2
       character str1*12,str2*12
!   Normalization
       real(kind=8) pig,a,g,pi,rhobar,wnorm,vnorm,wes
       real(kind=8) dg2rad 
!   Read model files
       real(kind=8) Ralpha(model),Rrho 
c       real(kind=8) dvvst1(1861,573),dvvst2(1861,573),dvvst3(1861,573)
       real(kind=8) dvvst1(3721,579),dvvst2(3721,579),dvvst3(3721,579),
     + dvvst4(3721,579),dvvst5(3721,579),dvvpt5(3721,579),
     + dvvst6(3721,579),dvvst7(3721,579),dvvst8(3721,579)


       common pig,a,g,pi,rhobar,wnorm,vnorm,wes
       common dg2rad

ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

!       Normalisation as minos: 
       pig=1.d0
       a=6371000.d0                    ! mean earth radius (m)
       g=6.6723d-11                    ! Newtonian gravity constant
       pi=dacos(-1.d0)
       rhobar=5515.d0                  ! mean earth density
       wnorm=dsqrt(pi*g*rhobar)       ! frequency normalization =environ 1.07 mHz
       vnorm=a*wnorm                   ! velocity normalization
       wes=1.d3*wnorm/(2.d0*pi)        ! frequence en mHz
       dg2rad=pi/360.d0
       step=4.d0


!   read center lat, center lon and az
        read(*,*) cclat,cclon,aaz,dep,width,tomo,NSmax
        clat=dfloat(cclat)
        clon=dfloat(cclon)
        az=dfloat(aaz)
        if (mod(width,2).eq.0) then
          degstart=width/2
        endif
        if (mod(width,2).ne.0) then 
          degstart=(width+1)/2
        endif

!----------------------------------------------------------------------
!        read spline functions 
!        read model coefficients
!----------------------------------------------------------------------

       open(3,file='../data/spline_basis.dat')
       do l=0,NN-1
       do k=1,model
       read(3,*) xx,yy,radial_basis(model+1-k,l)
       enddo
       enddo
       close(3)

       write(111,*) 'Spline base computed' 

       open(1,file='../models/SEISGLOB2.sph')
       do k=1,np 
       read(1,*) p(k),dp(k)       
       enddo
       close(1)

       open(1,file='../models/SEISGLOB1.sph')
       do k=1,np1 
       read(1,*) p1(k),dp1(k)       
       enddo
       close(1)

       write(111,*) 'Model coefficients read' 

!----------------------------------------------------------------------
!        Compute the Great-Circle-Path 
!----------------------------------------------------------------------

       az2=az+180.d0
       if (az2.ge.360.d0) az2=az2-360.d0
       dist=dfloat(degstart)*111.195 !90.d0*111.195
       ddeg=1       
       pnts=0.d0

c       write(*,*) clat,clon,az
       call MYGDS(clat,clon,az,dist,slat,slon)       
c       write(*,*) clat,clon,az2
       call MYGDS(clat,clon,az2,dist,elat,elon)   
c       write(*,*)    
c       write(*,*) elat,elon
c       write(*,*) slat,slon
       call MYGRT(elat,elon,slat,slon,dist,az12,az21,Gc)
c      write(*,*) dist,az12
c       call cgrc(elat,elon,slat,slon,ddeg,pnts,nrp,igrc)

       open(1,file='../coupes_NEW.xyz')
       open(2,file='../coupes_NEW.xy')
       open(3,file='../points_ligne.xy')
       open(6,file='../points_chauds_coupe.xy')
       open(4,file='../points_chauds_map.xy')
       open(5,file='../data/points_chauds.xy')
       open(7,file='../data/catalogue.xy')  
       open(8,file='../eq_coupe.xy')       
       write(6,*) "NaN ", "NaN"
       write(8,*) "NaN ", "NaN"
       pnts(1,1)=elon
       pnts(2,1)=elat
!        pnts(3,1)=dfloat(degstart)-0.d0
! Modified May 2017 : Problem plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pnts(3,1)=-dfloat(degstart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(1,*) pnts(1,1),pnts(2,1),90.d0-pnts(3,1)
       write(2,*) pnts(1,1),pnts(2,1)
       do k=1,55
       read(5,*) plon(k),plat(k)
       if ((elon.le.plon(k)+step).and.(elon.ge.plon(k)-step)
     +  .and.(elat.le.plat(k)+step).and.(elat.ge.plat(k)-step)) then
c       write(6,'(a)') '>'
       write(6,*) (90.d0-pnts(3,1)), 6281.
c       write(4,'(a)') '>'
       write(4,*) plon(k),plat(k)
       plon(k)=1e16
       plat(k)=1e16
       endif
       enddo
       close(5)
       do k=1,4995
       read(7,*) eqlon(k),eqlat(k),eqdep(k)
       if ((elon.le.eqlon(k)+step).and.(elon.ge.eqlon(k)-step)
     +  .and.(elat.le.eqlat(k)+step).and.(elat.ge.eqlat(k)-step)
     + .and.(eqdep(k).ge.50)) then
c       write(8,'(a)') '>'
       write(8,*) (90.d0-pnts(3,1)), 6371-eqdep(k)
       eqlon(k)=1e16
       eqlat(k)=1e16
       endif
       enddo
       close(7)
 
       c=2
       do k=degstart-1,1,-1
       dist=dfloat(k)*111.195
       call MYGDS(clat,clon,az2,dist,blat,blon)
       call MYGRT(elat,elon,blat,blon,dist,az12,az21,Gc)
       pnts(1,c)=blon
       pnts(2,c)=blat
!        pnts(3,c)=dfloat(degstart)-dist
!
! Modified May 2017 : Problem plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pnts(3,c)=dfloat(degstart-degstart-k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(1,*) pnts(1,c),pnts(2,c),90.d0-pnts(3,c)
       write(2,*) pnts(1,c),pnts(2,c)
       if (k.eq.degstart) then
c       write(3,'(a)') '> -Ggreen'
       write(3,*) pnts(1,c),pnts(2,c) 
       endif
c       if (k.eq.30) then
c       write(3,'(a)') '> -Gred'
c       write(3,*) pnts(1,c),pnts(2,c) 
c       endif
       do kk=1,55
       if ((pnts(1,c).le.plon(kk)+step).and.
     + (pnts(1,c).ge.plon(kk)-step)
     +  .and.(pnts(2,c).le.plat(kk)+step).and.
     + (pnts(2,c).ge.plat(kk)-step)) then
c       write(*,*) plon(kk)-step,pnts(1,c),plon(kk)+step
c       write(*,*) plat(kk)-step,pnts(2,c),plat(kk)+step
c       write(4,'(a)') '>'
       write(4,*) plon(kk),plat(kk)
c       write(6,'(a)') '>'
       write(6,*) (90.d0-pnts(3,c)), 6281.
       plon(kk)=1e16
       plat(kk)=1e16
       endif
       enddo
       do kk=1,4995
       if ((pnts(1,c).le.eqlon(kk)+step).and.
     + (pnts(1,c).ge.eqlon(kk)-step)
     +  .and.(pnts(2,c).le.eqlat(kk)+step).and.
     + (pnts(2,c).ge.eqlat(kk)-step)
     + .and.(eqdep(kk).ge.50)) then
c       write(8,'(a)') '>'
       write(8,*) (90.d0-pnts(3,c)), 6371-eqdep(kk)
       eqlon(kk)=1e16
       eqlat(kk)=1e16
       endif
       enddo
       c=c+1
       enddo
       pnts(1,c)=clon
       pnts(2,c)=clat
!        pnts(3,c)=0.d0
! Modified May 2017 : Problem plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pnts(3,c)=dfloat(degstart-degstart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(1,*) pnts(1,c),pnts(2,c),90.d0
       write(2,*) pnts(1,c),pnts(2,c)
c       write(3,'(a)') '> -Gyellow'
       write(3,*) pnts(1,c),pnts(2,c) 
       do kk=1,55
       if ((pnts(1,c).le.plon(kk)+step).and.
     + (pnts(1,c).ge.plon(kk)-step)
     +  .and.(pnts(2,c).le.plat(kk)+step).and.
     + (pnts(2,c).ge.plat(kk)-step)) then
c       write(4,'(a)') '>'
       write(4,*) plon(kk),plat(kk)
c       write(6,'(a)') '>'
       write(6,*) (90.d0), 6281.
       plon(kk)=1e16
       plat(kk)=1e16
       endif
       enddo
       do kk=1,4995
       if ((pnts(1,c).le.eqlon(kk)+step).and.
     + (pnts(1,c).ge.eqlon(kk)-step)
     +  .and.(pnts(2,c).le.eqlat(kk)+step).and.
     + (pnts(2,c).ge.eqlat(kk)-step)
     + .and.(eqdep(kk).ge.50)) then
c       write(8,'(a)') '>'
       write(8,*) (90.d0-pnts(3,c)), 6371-eqdep(kk)
       eqlon(kk)=1e16
       eqlat(kk)=1e16
       endif
       enddo
       c=c+1
       do k=1,degstart-1
       dist=dfloat(k)*111.195
       call MYGDS(clat,clon,az,dist,blat,blon)
       call MYGRT(elat,elon,blat,blon,dist,az12,az21,Gc)
       pnts(1,c)=blon
       pnts(2,c)=blat
!        pnts(3,c)=-(-degstart+dfloat(degstart+k))
! Modified May 2017 : Problem plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pnts(3,c)=dfloat(k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(1,*) pnts(1,c),pnts(2,c),90.d0-pnts(3,c)
       write(2,*) pnts(1,c),pnts(2,c)
c       if ((k.eq.30).or.(k.eq.60)) then
c       write(3,'(a)') '> -Gred'
c       write(3,*) pnts(1,c),pnts(2,c) 
c       endif
       do kk=1,55
       if ((pnts(1,c).le.plon(kk)+step).and.
     + (pnts(1,c).ge.plon(kk)-step)
     +  .and.(pnts(2,c).le.plat(kk)+step).and.
     + (pnts(2,c).ge.plat(kk)-step)) then
c       write(4,'(a)') '>'
       write(4,*) plon(kk),plat(kk)
c       write(6,'(a)') '>'
       write(6,*) (90.d0-pnts(3,c)), 6281.
       plon(kk)=1e16
       plat(kk)=1e16
       endif
       enddo
       do kk=1,4995
       if ((pnts(1,c).le.eqlon(kk)+step).and.
     + (pnts(1,c).ge.eqlon(kk)-step)
     +  .and.(pnts(2,c).le.eqlat(kk)+step).and.
     + (pnts(2,c).ge.eqlat(kk)-step)
     + .and.(eqdep(kk).ge.50)) then
c       write(8,'(a)') '>'
       write(8,*) (90.d0-pnts(3,c)), 6371-eqdep(kk)
       eqlon(kk)=1e16
       eqlat(kk)=1e16
       endif
       enddo
       c=c+1
       enddo
       pnts(1,width+1)=slon
       pnts(2,width+1)=slat
       dist=dfloat(degstart)*111.195
       call MYGRT(elat,elon,slat,slon,dist,az12,az21,Gc)
!        pnts(3,width+1)=-(-degstart+dfloat(degstart+k))
! Modified May 2017 : Problem plot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pnts(3,width+1)=dfloat(degstart)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(1,*) pnts(1,width+1),pnts(2,width+1),90.d0-pnts(3,width+1)
       write(2,*) pnts(1,width+1),pnts(2,width+1)
c       write(3,'(a)') '> -Gred'
       write(3,*) pnts(1,width+1),pnts(2,width+1) 
       do kk=1,55
       if ((pnts(1,width+1).le.plon(kk)+step).and.
     + (pnts(1,width+1).ge.plon(kk)-step)
     +  .and.(pnts(2,width+1).le.plat(kk)+step).and.
     + (pnts(2,width+1).ge.plat(kk)-step)) then
c       write(4,'(a)') '>'
       write(4,*) plon(kk),plat(kk)
c       write(6,'(a)') '>'
       write(6,*) (90.d0-pnts(3,width+1)), 6281.
       plon(kk)=1e16
       plat(kk)=1e16
       endif
       enddo
       do kk=1,4995
       if ((pnts(1,width+1).le.eqlon(kk)+step).and.
     + (pnts(1,width+1).ge.eqlon(kk)-step)
     +  .and.(pnts(2,width+1).le.eqlat(kk)+step).and.
     + (pnts(2,width+1).ge.eqlat(kk)-step)
     + .and.(eqdep(kk).ge.50)) then
c       write(8,'(a)') '>'
       write(8,*) (90.d0-pnts(3,width+1)), 6371-eqdep(kk)
       eqlon(kk)=1e16
       eqlat(kk)=1e16
       endif
       enddo
       close(1)
       close(2)
       close(3)
       close(4)
       close(6)
       close(8)

!----------------------------------------------------------------------
!        Make profile .xyz file 
!----------------------------------------------------------------------
   
       open(1,file='../models/S40RTS_5km.sph')
       open(2,file='../models/S362WMANIM_5km.sph')
       open(3,file='../models/SEMUCB_WM1_5km.sph')   
       open(5,file='../models/S12RTS_5km.sph')   
       open(6,file='../models/P12RTS_5km.sph')  
       open(7,file='../models/SGLOBE_5km.sph')   
       open(8,file='../models/3D2016_5km.sph') 
       cd=1
       do depth=2890,0,-5

!      Comme dans S40RTS
       Ralpha(cd)=((2.d0-3.d0)/(6371.d0-3480.d0))
     + *(6371.d0-dfloat(depth))+
     + 2.d0-((2.d0-3.d0)/(6371.d0-3480.d0))*6371.d0
       Ralpha(cd)=1.d0/Ralpha(cd)
       Rrho=50.d-2

       read(1,'(3721f15.8)') (dvvst1(s,cd),s=1,3721) 
       read(2,'(3721f15.8)') (dvvst2(s,cd),s=1,3721)        
       read(3,'(3721f15.8)') (dvvst3(s,cd),s=1,3721) 
       read(5,'(3721f15.8)') (dvvst5(s,cd),s=1,3721) 
       read(6,'(3721f15.8)') (dvvpt5(s,cd),s=1,3721) 
       read(7,'(3721f15.8)') (dvvst6(s,cd),s=1,3721) 
       read(8,'(3721f15.8)') (dvvst7(s,cd),s=1,3721) 


       cd=cd+1
       enddo

       close(1)
       close(2)
       close(3)
       close(5)
       close(6)
       close(7)
       close(8)

c       cd=1
c       do depth=2890,30,-50
c       write(4,'(1681f15.8)') (dvvst1(s,cd),s=1,1681) 
c       write(5,'(1681f15.8)') (dvvst2(s,cd),s=1,1681) 
c       if (depth.ge.50) then
c        write(6,'(1681f15.8)') (dvvst3(s,cd),s=1,1681) 
c       endif
c
c       cd=cd+10
c       enddo
c       close(4)
c       close(5)
c       close(6)
c
c       STOP

       if (tomo.eq.1) then
       call sph2geo22(40,p,np,radial_basis,
     +       '../coupes_NEW_SEISGLOB2_VP.res',
     +   '../coupes_NEW_SEISGLOB2_VS.res',
     + '../coupes_NEW_SEISGLOB2_RHO.res',pnts,NN,model,dep,width)
       endif

       if (tomo.eq.2) then
       call sph2geo2S40RTS(NSmax,p,np,radial_basis,
     + '../coupes_NEW_S40RTS_VP.res',
     + '../coupes_NEW_S40RTS_VS.res',
     + '../coupes_NEW_S40RTS_RHO.res',pnts,NN,model,
     +   dvvst1,Ralpha,Rrho,579,dep,width)
       endif

       if (tomo.eq.4) then
       call sph2geo2S362WMANIM(NSmax,p,np,radial_basis,
     + '../coupes_NEW_S362WMANIM_VP.res',
     + '../coupes_NEW_S362WMANIM_VS.res',
     + '../coupes_NEW_S362WMANIM_RHO.res',pnts,NN,model,
     + dvvst2,579,dep,width)
       endif

       if (tomo.eq.3) then 
       call sph2geo2SEMUCB(NSmax,p,np,radial_basis,
     + '../coupes_NEW_SEMUCBWM1_VP.res',
     + '../coupes_NEW_SEMUCBWM1_VS.res',
     + '../coupes_NEW_SEMUCBWM1_RHO.res',pnts,NN,model,
     + dvvst3,579,dep,width)
       endif

       if (tomo.eq.5) then
       call sph2geo21(NSmax1,p1,np1,radial_basis,
     +       '../coupes_NEW_SEISGLOB1_VP.res',
     +   '../coupes_NEW_SEISGLOB1_VS.res',
     + '../coupes_NEW_SEISGLOB1_RHO.res',pnts,NN,model,dep,width)
       endif

       if (tomo.eq.6) then
       call sph2geo2SP12RTS(NSmax,p,np,radial_basis,
     + '../coupes_NEW_SP12RTS_VP.res',
     + '../coupes_NEW_SP12RTS_VS.res',
     + '../coupes_NEW_SP12RTS_RHO.res',pnts,NN,model,
     +   dvvst5,dvvpt5,579,dep,width)
       endif

       if (tomo.eq.7) then
       call sph2geo2SGLOBE(NSmax,p,np,radial_basis,
     + '../coupes_NEW_SGLOBE_VP.res',
     + '../coupes_NEW_SGLOBE_VS.res',
     + '../coupes_NEW_SGLOBE_RHO.res',pnts,NN,model,
     +   dvvst6,579,dep,width)
       endif

       if (tomo.eq.8) then
       call sph2geo23D2016(NSmax,p,np,radial_basis,
     + '../coupes_NEW_3D2016_VP.res',
     + '../coupes_NEW_3D2016_VS.res',
     + '../coupes_NEW_3D2016_RHO.res',pnts,NN,model,
     +   dvvst7,579,dep,width)
       endif

       if (tomo.eq.9) then
       open(9,file='../models/YOURMODEL_5km.sph')   
       cd=1
       do depth=2890,0,-5
       read(9,'(3721f15.8)') (dvvst8(s,cd),s=1,3721) 
       cd=cd+1
       enddo
       close(9)

       call sph2geo2YOURMODEL(NSmax,p,np,radial_basis,
     + '../coupes_NEW_MYMODEL_VP.res',
     + '../coupes_NEW_MYMODEL_VS.res',
     + '../coupes_NEW_MYMODEL_RHO.res',pnts,NN,model,
     +   dvvst8,579,dep,width)
       endif

!----------------------------------------------------------------------
!        Files for GMT plots 
!----------------------------------------------------------------------
      
       open(1,file='../lignes_depth.xy')
c       write(1,'(a)') '> -W0.5p,black,--'
       do k=-900,2700
       write(1,'(f8.3,1x,i4)') dfloat(k)/10.d0,5961
       enddo
c       write(1,'(a)') '> -W0.5p,black,--'
       do k=-900,2700
       write(1,'(f8.3,1x,i4)') dfloat(k)/10.d0,5701
       enddo
c       write(1,'(a)') '> -W0.5p,black,--'
       do k=-900,2700
       write(1,'(f8.3,1x,i4)') dfloat(k)/10.d0,5371
       enddo
       close(1)

       open(1,file='../points_ligne_coupe.xy')
c       write(1,'(a)') '> -Ggreen'
       write(1,*) 90-width/2, '   ',6371
c       write(1,'(a)') '> -Gyellow'
       write(1,*) 90, '   ',6371
c       write(1,'(a)') '> -Gred'
       write(1,*) 90+width/2, '   ',6371
       close(1)



       end 

!---------------------------------------------------------------------
!      subfonction Great-circle-path 
!---------------------------------------------------------------------

      SUBROUTINE MYGDS(ALAT,ALON,AZ,DIST,BLAT,BLON)

c Geodesic program. Given epicenter (alat,alon), azimuth of travel and distance
c traveled, finds arrival point (blat,blon);  dist in km

      real(kind=8) ALAT,ALON,AZ,DIST,BLAT,BLON
      real(kind=8) ALT,ALN,A0,B0,F,OPI,RAD,
     + E2,EPS,SIN1,V,SINA,COSA,COS1,TC2,C2,EPS0,TB1,
     + B1,TAN1,S1,G0,G2,G4,G6,SIGM,U1P,SIN2P,SIN4P,TSS,
     +S12,SIGP,T1,T2,T3,U2P,SINU1,U2,SINP1,A1,Q1,Q2,
     + X1,AMU,DLAMB
 
      ALT=ALAT
      ALN=ALON
      IF (ALAT.EQ.90.d0) ALAT=ALAT-0.000001
      IF (AZ.EQ.90.d0.OR.AZ.EQ.270.d0) AZ=AZ-0.000001
      A0=6378.388
      B0=6356.912
      F=(A0-B0)/A0
      PI=4.*DATAN(1.d0)
      RAD=PI/180.d0
      ALAT=ALAT*RAD
      ALON=ALON*RAD
      AZ=AZ*RAD
      E2=1.d0-(B0**2.d0)/(A0**2.d0)
      EPS=E2/(1.d0-E2)
      SIN1=DSIN(ALAT)
      V=A0/(DSQRT(1.d0-E2*(SIN1**2.d0)))
      SINA=DSIN(AZ)
      COSA=DCOS(AZ)
      COS1=DCOS(ALAT)
      TC2=(COSA*COSA)*(COS1*COS1)
      C2=TC2+(SIN1*SIN1)
      EPS0=C2*EPS
      TB1=DSQRT(1.d0+EPS*TC2)
      B1=(V*TB1)/(1.d0+EPS0)
      TAN1=DTAN(ALAT)
      IF (COSA.EQ.0.d0) COSA=0.00001
      S1=DSQRT(1.d0+EPS0)
      G0=1.d0-(EPS0/4.d0)+((7.d0*EPS0*EPS0)/64.d0)-
     + ((15.d0*(EPS0**3.d0))/256.d0)
      G2=(EPS0/8.d0)-(.0625*EPS0*EPS0)+((145.d0*(EPS0**3))/2048.d0)
      G4=((5.d0*EPS0*EPS0)/256.d0)-((5.d0*(EPS0**3.d0))/256.d0)
      G6=(29.d0*(EPS0**3.d0))/6144.d0
      SIGM=(DIST*G0)/B1
      U1P=DATAN2(TAN1,(COSA*S1))
      SIN2P=DSIN(2.d0*U1P)
      SIN4P=DSIN(4.d0*U1P)
      TSS=((EPS0/4.d0)-((EPS0*EPS0)/8.d0))
      S12=2.d0*U1P-TSS*SIN2P-((EPS0*EPS0)/128.d0)*SIN4P
      SIGP=S12+SIGM
      T1=SIGM+(2.d0*G2*DSIN(SIGM))*DCOS(SIGP)
      T2=(2.d0*G4*DSIN(2.d0*SIGM))*DCOS(2.d0*SIGP)
      T3=(2.d0*G6*DSIN(3.d0*SIGM))*DCOS(3.d0*SIGP)
      U2P=U1P+T1+T2+T3
      SINU1=TAN1/(DSQRT(1.d0+EPS+(TAN1*TAN1)))
      C=DSQRT(C2)
      SINU2=(((B1*C)/B0)*DSIN(U2P))-((EPS-EPS0)/(1.d0+EPS0))*SINU1
      U2=ASIN(SINU2)
      SINP1=SINU2/(DSQRT(1.d0-E2*(DCOS(U2)*DCOS(U2))))
      BLAT=ASIN(SINP1)
      A1=B1*(DSQRT(1.d0+EPS0))
      IF (DCOS(U2).EQ.0.d0) U2=U2-0.00001
      Q1=(A1*DCOS(U2P))/(A0*DCOS(U2))
      if(q1.gt.-1.d0.and.q1.lt.1.d0) Q2=ACOS(Q1)
      if(q1.eq.1.d0) q2=0.d0
      if(q1.eq.-1.d0) q2=pi
      if(q1.gt.1.d0) go to 100
      if(q1.lt.-1.d0) go to 200
300   X1=SIN1*SINA
      AMU=DATAN2(X1,COSA)
      AZ=AZ/RAD
      U2P=U2P/RAD
      IF (AZ.GT.180.d0) Q2=-Q2
      IF (U2P.GT.180.d0.OR.U2P.LT.0.d0) Q2=-Q2
      DLAMB=Q2-AMU
      BLON=DLAMB+ALON
      BLAT=BLAT/RAD
      BLON=BLON/RAD
      IF (DABS(BLON).GT.180.d0) BLON=BLON-SIGN(360.d0,BLON)
      ALAT=ALT
      ALON=ALN
      RETURN
c100   write(6,*)'Flag in "Q2=ACOS(Q1)",Q1 = ',q1,';Q2 taken as 0'
100   q2=0.d0
      go to 300
c200   write(6,*)'     Flag in "Q2=ACOS(Q1)", Q1= ',q1,';Q2 taken as PI'
200   q2=pi
      go to 300
      END SUBROUTINE MYGDS

      SUBROUTINE MYGRT(ALAT1,ALON1,ALAT2,ALON2,DISD,AZ12,AZ21,GC)

c Great circle program. Given epicenter and station, finds distance,
c take-off and back- azimuths, and length of great circle

      real(kind=8) ALAT1,ALON1,ALAT2,ALON2,DISD,AZ12,AZ21,GC
      real(kind=8) PI,ATH,BTH,RAD,H,P,GR,TR,SINTR,COSTR,R1,Z1,G,T,
     + SINT,COST,R2,DG,COSDG,SINDG,DGR,DT,X,Y,Z,Q,COS12,P0,B0,E0,
     + C0,C2,C4,U0,U,DIST

      PI = 4.d0*DATAN(1.d0)
      ATH=6378.388
      BTH=6356.912
      RAD = PI/180.d0
      H = 1.d0 - BTH*BTH/(ATH*ATH)
      P = H/(1.d0 - H)
      GR = ALON1*RAD
      TR = ALAT1*RAD
      SINTR = DSIN(TR)
      COSTR = DCOS(TR)
      IF (SINTR .EQ. 0.d0) SINTR = .000001
      IF (COSTR .EQ. 0.d0) COSTR = .000001
      R1 = ATH/DSQRT(1.d0 - H*SINTR*SINTR)
      Z1 = R1*(1.d0 - H)*SINTR
      G = ALON2*RAD
      T = ALAT2*RAD
      IF (T .EQ. 0.d0) T = .00001
      SINT = DSIN(T)
      COST = DCOS(T)
      R2 = ATH/DSQRT(1.d0 - H*SINT*SINT)
      DG = G - GR
      COSDG = DCOS(DG)
      COSDG = DCOS(DG)
      SINDG = DSIN(DG)
      DGR = GR - G
      DT = T - TR
      Q = SINT*COSTR/((1.d0 + P)*COST*SINTR) + H*R1*COSTR/(R2*COST)
      X = R2*COST*COSDG
      Y = R2*COST*SINDG
      Z = R2*(1.d0 - H)*SINT
      AZ12 = DATAN2(SINDG,(Q - COSDG)*SINTR)
      Q = SINTR*COST/(COSTR*SINT*(1.d0 + P)) + H*R2*COST/(R1*COSTR)
      AZ21 = DATAN2(DSIN(DGR),SINT*(Q - DCOS(DGR)))
      COS12 = DCOS(AZ12)
      CTA2 = COSTR*COSTR*COS12*COS12
      P0 = P*(CTA2 + SINTR*SINTR)
      B0 = (R1/(1.d0 + P0))*DSQRT(1.d0 + P*CTA2)
      E0 = P0/(1.d0 + P0)
      GC = 2.d0*PI*B0*DSQRT(1.d0 + P0)*(1.d0 - E0*(.25 + E0*(3.d0/64.d0
     *                                          + 5.d0*E0/256.d0)))
      C0 = 1.d0 + P0*(.25 - P0*(3.d0/64.d0 - 5.d0*P0/256.d0))
      C2 = P0*(-.125 + P0*(1.d0/32.d0 - 15.d0*P0/1024.d0))
      C4 = (-1.d0/256.d0 + 3.d0*P0/1024.d0)*P0*P0
      U0 = DATAN2(SINTR,COSTR*COS12*DSQRT(1.d0 + P0))
      U = DATAN2(R1*SINTR + (1.d0 + P0)*(Z - Z1),(X*COS12 - Y*SINTR*
     *    DSIN(AZ12))*DSQRT(1.d0 + P0))
      DISD = U - U0
      IF (U .LT. U0) DISD = PI + PI + DISD
      DIST = B0*(C0*( DISD ) + C2*(DSIN(U + U) - DSIN(U0 + U0))
     *                       + C4*(DSIN(4.d0*U) - DSIN(4.d0*U0)))
      DISD = DISD/RAD
      AZ12 = AZ12/RAD
      AZ21 = AZ21/RAD
      IF (AZ12 .LT. 0.d0) AZ12 = 360.d0 + AZ12
      IF (AZ21 .LT. 0.d0) AZ21 = 360.d0 + AZ21
        if(disd.gt.355.d0) disd=360.d0-disd
      RETURN
      END SUBROUTINE MYGRT



