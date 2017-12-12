!=====================================================================
!                          PROFILE-MAKER PROGRAM
!       S. Durand, May 2016, V0.0 
!=====================================================================
       implicit none
!    Nombres de degres inverse pour chaque donnees
       integer c,NSmax,np,model,NN,spline,k,l,ddeg,kk,ddepth,lay,
     + tomo,para,np1,NSmax1,lllmax
       parameter (np=35280,model=2891,NN=21,spline=1,
     + np1=9240,NSmax1=20)
!   earth model
       real(kind=8) r(model),rho(model),kappa(model),mu(model),
     +        rhos2(model),qkap(model),qmu(model),pg(model)
       real(kind=8) radial_basis(model,0:NN-1)
       real(kind=8) p(np),xx,yy,dp(np)
       real(kind=8) p1(np1),dp1(np1)
       real(kind=8) depth,rad(model)
       integer nrp,igrc
!   Nom fichiers
       character*(80) filename2
       character str1*12,str2*12
!   Normalization
       real(kind=8) pig,a,g,pi,rhobar,wnorm,vnorm,wes
       real(kind=8) dg2rad 

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


!   read center lat, center lon and az
        read(*,*) ddepth,tomo,NSmax
        depth=dfloat(ddepth)

!----------------------------------------------------------------------
!        read spline functions 
!        read model coefficients
!----------------------------------------------------------------------

       open(3,file='../data/spline_basis_new.dat')
       do l=0,NN-1
       do k=1,model
       read(3,*) rad(k),yy,radial_basis(k,l)
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
!        Make map .xyz file 
!----------------------------------------------------------------------
       if (tomo.eq.1) then
       call sph2geo22(40,p,np,radial_basis,
     +       '../output_files_map/map_NEW_SEISGLOB2_VP.xyz',
     + '../output_files_map/map_NEW_SEISGLOB2_VS.xyz',
     + '../output_files_map/map_NEW_SEISGLOB2_RHO.xyz',NN,model,
     + depth,rad)
       endif

       if (tomo.eq.2) then
       call sph2geo2s40rts(NSmax,p,np,radial_basis,
     +       '../output_files_map/map_NEW_S40RTS_VP.xyz',
     +  '../output_files_map/map_NEW_S40RTS_VS.xyz',
     +  '../output_files_map/map_NEW_S40RTS_RHO.xyz',NN,model,
     + ddepth,rad)
       endif

       if (tomo.eq.3) then
       call sph2geo2semucb(NSmax,p,np,radial_basis,
     + '../output_files_map/map_NEW_SEMUCBWM1_VP.xyz',
     +       '../output_files_map/map_NEW_SEMUCBWM1_VS.xyz',
     + '../output_files_map/map_NEW_SEMUCBWM1_RHO.xyz',NN,model,
     + ddepth,rad)
       endif

       if (tomo.eq.4) then
       call sph2geo2s362wmanim(NSmax,p,np,radial_basis,
     + '../output_files_map/map_NEW_S362WMANIM_VP.xyz',
     +       '../output_files_map/map_NEW_S362WMANIM_VS.xyz',
     + '../output_files_map/map_NEW_S362WMANIM_RHO.xyz',NN,model,
     + ddepth,rad)
       endif

       if (tomo.eq.5) then
       call sph2geo21(NSmax1,p1,np1,radial_basis,
     +       '../output_files_map/map_NEW_SEISGLOB1_VP.xyz',
     + '../output_files_map/map_NEW_SEISGLOB1_VS.xyz',
     + '../output_files_map/map_NEW_SEISGLOB1_RHO.xyz',NN,model,
     + depth,rad)
       endif

       if (tomo.eq.6) then
       call sph2geo2sp12rts(NSmax,p,np,radial_basis,
     +       '../output_files_map/map_NEW_SP12RTS_VP.xyz',
     +  '../output_files_map/map_NEW_SP12RTS_VS.xyz',
     +  '../output_files_map/map_NEW_SP12RTS_RHO.xyz',NN,model,
     + ddepth,rad)
       endif

       if (tomo.eq.7) then
       call sph2geo2sglobe(NSmax,p,np,radial_basis,
     +       '../output_files_map/map_NEW_SGLOBE_VP.xyz',
     +  '../output_files_map/map_NEW_SGLOBE_VS.xyz',
     +  '../output_files_map/map_NEW_SGLOBE_RHO.xyz',NN,model,
     + ddepth,rad)
       endif

       if (tomo.eq.8) then
       call sph2geo23d2016(NSmax,p,np,radial_basis,
     +       '../output_files_map/map_NEW_3D2016_VP.xyz',
     +  '../output_files_map/map_NEW_3D2016_VS.xyz',
     +  '../output_files_map/map_NEW_3D2016_RHO.xyz',NN,model,
     + ddepth,rad)
       endif

       if (tomo.eq.9) then
       call sph2geo2yourmodel(NSmax,p,np,radial_basis,
     +       '../output_files_map/map_NEW_MYMODEL_VP.xyz',
     +  '../output_files_map/map_NEW_MYMODEL_VS.xyz',
     +  '../output_files_map/map_NEW_MYODEL_RHO.xyz',NN,model,
     + ddepth,rad)
       endif
       end 



