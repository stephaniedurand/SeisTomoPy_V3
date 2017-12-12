cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                      c
c      Compute the model file needed for SeisTomoPy    c
c             from S. Durand March 2017                c
c                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                      c
c    In the directory where all the *xyz files are     c
c    - Run make clean, make                            c
c    - This code will create two model files:          c
c      YOURMODEL_1km.sph  and YOURMODEL_5km.sph        c
c    - Then you just have to copy this file in         c
c      ../fortran_files/models                         c
c                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       implicit none
       include 'in_HarmSpher'
c    
       integer nx,ny,l,m,s,t,x,y,ind,j,i,cY,cd,depth,NS,depth2
       integer model,modetype,modetype2,cc,kk,depth1,ndepth,ndepth2
       parameter (model=788,ndepth=579,NS=3721,ndepth2=2891)
       integer indx,form,lm,xxx,compt,k,np,lllmax
       parameter (nx=179,ny=360)
       real*8   sp2(0:60),sp(0:60),sptotal,sptotal2
       real*8 dataf(nx,ny),dataf2(nx,ny)
       real*8 dvvst(NS,ndepth),pi
       real*8 dvvst2(NS,ndepth2)
       real*8 racin(nlat),poids(nlat),xx,yyy
       real*8 ft(nlat,nlong),ft2(nlat,nlong)
       complex*16 yy(lmmax),yy2(lmmax)
       character ind2*4
c
        pi=dacos(-1.d0)                

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      decomposition en HS de chaque carte
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       open(10,file='../YOURMODEL_5km.sph')
       write(*,*) "Compute first file YOURMODEL_5km.sph"
       cd=1
       do depth=2890,0,-5
       write(ind2,'(i4)') depth
       open(1,file='../model/vs.'//trim(adjustl(ind2))//'.xyz')

       do y=0,180
       do x=-89,89
       read(1,*) xx,yyy,dataf(((90-x)),(y+180))
c       write(*,*) dataf(((90-x)),(y+180))
       enddo
       enddo
       do y=181,359
       do x=-89,89
       read(1,*) xx,yyy,dataf(((90-x)),(y-180))
c       write(*,*) dataf(((90-x)),(y-180))
       enddo
       enddo
       close(1)
       dataf=dataf/1.d2

       call gauss(racin,poids)
       call poly(racin)
       call norme
       call interpol(dataf,nx,ny,ft,nlat,nlong,racin,1)
       call spat_spec(ft,poids,yy,1) 

       lllmax=60
       call spec(yy,lllmax,sp,sptotal)

       cY=1
       do s=0,lllmax
       dvvst(cY,cd)=dsqrt(4.d0*pi)*real(yy(indx(s,0)))
       cY=cY+1
       do k=1,s
       dvvst(cY,cd)=dsqrt(2.d0*pi)
     + *real(yy(indx(s,k)))
       cY=cY+1
       dvvst(cY,cd)=dsqrt(2.d0*pi)
     +  *aimag(yy(indx(s,k)))
       cY=cY+1
       enddo
       enddo

       write(10,'(3721f15.8)') (dvvst(s,cd),s=1,3721)

       cd=cd+1
       enddo
       close(10)

       write(*,*) "Finish"
       
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      decomposition en HS de chaque carte
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       open(10,file='../YOURMODEL_1km.sph')
       cd=1
       write(*,*) "Compute second file YOURMODEL_1km.sph"
       do depth=2890,0,-1
       write(ind2,'(i4)') depth
       open(1,file='../model/vs.'//trim(adjustl(ind2))//'.xyz')

       do y=0,180
       do x=-89,89
       read(1,*) xx,yyy,dataf2(((90-x)),(y+180))
c       write(*,*) dataf(((90-x)),(y+180))
       enddo
       enddo
       do y=181,359
       do x=-89,89
       read(1,*) xx,yyy,dataf2(((90-x)),(y-180))
c       write(*,*) dataf(((90-x)),(y-180))
       enddo
       enddo
       dataf2=dataf2/1.d2
       close(1)

       call gauss(racin,poids)
       call poly(racin)
       call norme
       call interpol(dataf2,nx,ny,ft2,nlat,nlong,racin,1)
       call spat_spec(ft2,poids,yy2,1) 

       lllmax=60
       call spec(yy2,lllmax,sp2,sptotal2)

       cY=1
       do s=0,lllmax
       dvvst2(cY,cd)=dsqrt(4.d0*pi)*real(yy2(indx(s,0)))
       cY=cY+1
       do k=1,s
       dvvst2(cY,cd)=dsqrt(2.d0*pi)
     + *real(yy2(indx(s,k)))
       cY=cY+1
       dvvst2(cY,cd)=dsqrt(2.d0*pi)
     +  *aimag(yy2(indx(s,k)))
       cY=cY+1
       enddo
       enddo

       write(10,'(3721f15.8)') (dvvst2(s,cd),s=1,3721)

       cd=cd+1
       enddo
       close(10)
       write(*,*) "Finish"

       end
