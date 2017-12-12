
c
c SPECTRAL SUBROUTINES (Yanick Ricard) 
c passees en double precision le 20 fevrier 95(sylvain)

        
c*************************************************************************
      subroutine gauss(racin,poids)
c
c     Computes the Gauss points (racin) and their weight (poids)
c
c*************************************************************************
       implicit none
      include 'in_HarmSpher'
      real*8 pi,preci,pz,pza,pzaa,zz,dpz,zint
      real*8 poids(nlat),racin(nlat),co
      integer i,j,mr
c
      parameter(preci=1.d-12)
c
      pi=4.d0*datan2(1.d0,1.d0)
c
      mr=nlat/2
      co=pi/(nlat+.5d0)
c
      do 115 i=1,mr  
            zz=dcos(co*(dble(i)-0.25d0))
100       continue 
            pza=1
            pz=zz
            do 110 j=2,nlat
              pzaa=pza
              pza=pz
              pz=((2.*dble(j)-1.)/dble(j))*zz*pza
     &        -((dble(j)-1.)/dble(j))
     &        *pzaa
110        continue
              dpz=nlat*(zz*pz-pza)/(zz*zz-1.d0)
              zint=zz
              zz=zint-pz/dpz
      if (dabs(zz-zint).gt.preci) go to 100
        racin(i)=zz
        racin(nlat+1-i)=-zz
        poids(i)=2.d0/((1.d0-zz*zz)*dpz*dpz)
        poids(nlat+1-i)=poids(i)
115   continue
        return
      end
c*************************************************************************
       integer function indx(l,m)
c
c      Cumulative index of degrees and orders
c
c*************************************************************************
       implicit none
       integer l,m
       indx=l*(l+1)/2+m+1
       return
       end

c*************************************************************************
       integer function par(l,m)             
c
c      Parity of Legendre functions
c
c*************************************************************************
       implicit none
       integer l,m
       par=1-2*mod(l-m,2)
       return
       end
c*************************************************************************
       subroutine poly(racin)      
c
c      Values of Legendre functions on Gauss points
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher' 
c                 
       real*8 racin(nlat)
       real*8 plm(lmmax)
       real*8 xx
       integer i,lm,mr
c
c
c       write(*,*) 'Poly has changed the 1st march 2005'
c       write(*,*) 'Poly is now normalized'
c
       mr=nlat/2
c
       do 1 i=1,mr  
               xx=racin(i)
c           call legendre(plm,xx)
c
c  The normalization is included in harm (not in legendre)
c  Higher degrees can now be reached
            call harm(plm,xx)
c
       do lm=1,lmmax
           ylm(lm,i)=plm(lm)
       enddo
c
c
 1     continue
       return
       end
c**********************************************************************
       subroutine dpoly(racin)
c
c     First derivatives of Legendre functions
c     (d Plm/d teta) and (m Plm/sin teta)
c
c**********************************************************************
      implicit none
      include 'in_HarmSpher'
c                 
      real*8 racin(nlat)
      real*8 ct,st,cot,co
      integer l,m,mr,lm,i,indx
c     
c
      mr=nlat/2
       
      do i=1,mr 
         ct=racin(i)
         st=dsqrt(1.d0-ct*ct)
         cot=ct/st           
            do l=0,lmax 
               do m=0,l 
               lm=indx(l,m)
               co=0.d0
               if(m.ne.l) co=ylm(lm+1,i) 
               dtylm(lm,i)=m*cot*ylm(lm,i)+co
               dpylm(lm,i)=m*ylm(lm,i)/st 
               enddo
            enddo
      enddo  
      return
      end
c**********************************************************************
c      subroutine ddpoly(racin)
c
c      Second derivatives of Legendre functions
c      (d2 Plm/ d teta 2), d (m Plm/sin teta)/d teta
c
c**********************************************************************
c      
c       include 'in_HarmSpher'
c                 
c      real*8 racin(nlat)
c      real*8 ct,st,cot,co,st2
c
c      mr=nlat/2
c      
c      do i=1,mr 
c        ct=racin(i)
c        st2=1.d0-ct*ct
c        st=dsqrt(st2)
c        cot=ct/st
c           do l=0,lmax 
c              ll=l*(l+1)
c              do m=0,l 
c              lm=indx(l,m)
c              dttylm(lm,i)=-ll*ylm(lm,i)-cot*dtylm(lm,i)+
c     &                      m*m/st2*ylm(lm,i)
c              dtpylm(lm,i)=m*dtylm(lm,i)/st-cot*dpylm(lm,i) 
c              enddo
c           enddo
c      enddo  
c      return
c      end
c*************************************************************************
       subroutine norme                      
c
c      Normalization
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher'
c
       integer l,m,lm,indx
       znorm(1)=1.d0
         do 1 l=0,lmax
          znorm(indx(l,0))=dsqrt(2.d0*l+1.d0)
            do 1 m=1,l
             lm=indx(l,m)
        znorm(lm)=-znorm(lm-1)/dsqrt(dble((l+m)*(l+1-m)))
 1     continue
c
       return
       end
c*************************************************************************
c*************************************************************************
       subroutine legendre(plm,x) 
c
c      Values de plm for all degrees and orders lower than lmax at x
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher'
       real*8 plm(lmmax),x 
       real*8 xx,somx2,ff,x0,xx0,x1,x2
       integer l,m,lm,mp1mp1,indx
c
               xx=x
               xx0=1.d0
               plm(1)=xx0
               somx2=dsqrt(1.d0-xx*xx)
               ff=1.d0
          do 2 m=0,lmax-1
                  x0=xx0
                  xx0=-x0*ff*somx2
                  mp1mp1=indx(m+1,m+1)
                  plm(mp1mp1)=xx0
                  x1=xx*ff*x0
                  plm(mp1mp1-1)=x1
c          write(*,*) m+1,m+1
c          write(*,*) m+1,m
                  ff=ff+2.d0
c             do 3 l=min(m+2,lmax),lmax
c          write (*,*) 'do',m+2,lmax
             do 3 l=m+2,lmax
                     x2=(xx*(2.d0*dble(l)-1.d0)*x1-
     &               dble((l+m-1))*x0)/dble(l-m)
c         write(*,*) l,m
                     lm=indx(l,m)
                     plm(lm)=x2
                     x0=x1
                     x1=x2
 3           continue
 2         continue
       return
       end
c  subroutine tire de spectral et qui peuvent etre compile 
c   sans l'option -Kc
c SPECTRAL SUBROUTINES (Yanick Ricard) 
ci passee en double precision par Sylvain
c
c*************************************************************************
       subroutine interpol(donnee,nd,nnd,data,nr,nnr,racin,isigne)
c
c     Interpolation from a regular grid to a gauss grid (example)
c
c*************************************************************************
      implicit none
      integer i,j,it,nr,nd,nnd,nnr,isigne,k
      real*8 donnee(nd,nnd),data(nr,nnr),teta(0:257)
      real*8 provi
      real*8 racin(nr),ip,ipp1,itp1,xt
      real*8 t,pi,p,umdt1,umdt2,dt2,dt1,c1,c2,c3,x1,x2
c
      common/x/provi(0:300,550)
c
      pi=4.d0*datan2(1.d0,1.d0)     
c
      do 2 i=1,nr
  2   teta(i)=dacos(racin(i))
      teta(0)=0.d0
      teta(nr+1)=pi
      c1=real(nnd)/real(nnr)
      c2=nd/pi
      c3=.5d0
c
      if(isigne.eq.1) then
       x1=0.d0
       x2=0.d0
       do 3 j=1,nnd
       x1=x1+donnee(1,j)
       x2=x2+donnee(nd,j)
       do 3 i=1,nd
3      provi(i,j)=donnee(i,j)
       do 4 j=1,nnd
       provi(0,j)=x1/nnd
4      provi(nd+1,j)=x2/nnd
c
       do 1 j=1,nnr
          p=(j-1.)*c1+1.
          ip=p
          ipp1=ip+1
          dt2=p-ip
          umdt2=1.-dt2
             if(ipp1.eq.nnd+1) ipp1=1
       do 1 i=1,nr
          t=teta(i)*c2+c3
          it=t
          itp1=it+1
          dt1=t-it
          umdt1=1.-dt1
c
                 data(i,j)=(provi(it,ip)*umdt1+
     &                      provi(itp1,ip)*dt1)*umdt2+
     &                     (provi(it,ipp1)*umdt1+
     &                      provi(itp1,ipp1)*dt1)*dt2
c
 1     continue
                   endif
c
       if(isigne.eq.-1) then
       x1=0.
       x2=0.
       do 6 j=1,nnr
       x1=x1+data(1,j)
       x2=x2+data(nr,j)
       do 6 i=1,nr
6      provi(i,j)=data(i,j)
       do 7 j=1,nnr
       provi(0,j)=x1/nnr
7      provi(nr+1,j)=x2/nnr
c
       it=1
       do 10 i=1,nd
          t=(i-c3)/c2
          do 8 k=it,nr+1
          xt=teta(k)
          if(xt.gt.t) go to 5
8         continue
5         itp1=k
          it=k-1
             dt1=(t-teta(it))/(teta(itp1)-teta(it))
             umdt1=1-dt1
       do 10 j=1,nnd
          p=(j-1.)/c1+1.
          ip=p
          ipp1=ip+1
          dt2=p-ip
          umdt2=1.-dt2
             if(ipp1.eq.nnr+1) ipp1=1
c
               donnee(i,j)=(provi(it,ip)*umdt1+
     &                      provi(itp1,ip)*dt1)*umdt2+
     &                     (provi(it,ipp1)*umdt1+
     &                      provi(itp1,ipp1)*dt1)*dt2
 10    continue
                   endif
       return
       end
c*************************************************************************
       subroutine spat_spec(data,poids,yy,isigne)            
c
c      expansion (isigne=1) of "data(lat,long)" into "yy(lm)"
c      construction of data from yy (isigne=-1)
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher'
       character*8 func
C     
       real*8 poids(nlat)
       real*8 data(nlat,nlong)
       complex*16 yy(lmmax)
       integer ipar,isigne  
c                  
      ipar=1
c
      if(isigne.eq.1)  then 
c........................................................................
c......Direct Fourier transform
c........................................................................
      call fourier(data,isigne)
c........................................................................
c......Gauss quadrature
c........................................................................
      func='ylm'
      call gaussquad(data,poids,yy,isigne,ipar,func) 
                      endif
       if(isigne.eq.-1) then 
c........................................................................
c......products by the Plm
c........................................................................
      func='ylm'
      call gaussquad(data,poids,yy,isigne,ipar,func)
c........................................................................
c......Inverse Fourier Transform
c........................................................................
      call fourier(data,isigne)
                       endif
       return
       end
c*************************************************************************
       subroutine fourier(data,isigne)
c
c      replace data by its FFT
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher'
c
       real*8 data(nlat,nlong),pippi
       real*8 data1(1024),data2(1024),z
       complex*16 fft1(1024),fft2(1024)
       integer k,i,j,ni,isigne
       parameter(pippi=6.28318530717959d0)
c
       if(isigne.eq.1) then
c                 
       z=pippi/dble(nlong)
c
       do 1 i=1,nlat,2
       ni=i+1
       do 2 j=1,nlong
       data1(j)=data(i,j)
       data2(j)=data(ni,j)
 2     continue
       call twofft(data1,data2,fft1,fft2,nlong)
       do 3 j=1,nlong,2
       k=(j+1)/2
       data(i,j)=dreal(fft1(k))*z
       data(i,j+1)=dimag(fft1(k))*z
       data(ni,j)=dreal(fft2(k))*z
       data(ni,j+1)=dimag(fft2(k))*z
 3     continue
 1     continue
c
                     endif
c
       if(isigne.eq.-1) then
c
       do 10 i=1,nlat,2
       ni=i+1
       do 20 j=1,nlong,2
       k=(j+1)/2
       fft1(k)=dcmplx(data(i,j)-data(ni,j+1),data(i,j+1)+data(ni,j))
       fft1(nlong+2-k)=dcmplx(data(i,j)+data(ni,j+1),
     &                     -data(i,j+1)+data(ni,j))
 20    continue
       fft1(nlong/2+1)=dcmplx(0.d0,0.d0)
       call four1(fft1,nlong,isigne)
       do 40 j=1,nlong
       data(i,j)=dreal(fft1(j))
       data(ni,j)=dimag(fft1(j))
 40    continue
 10    continue
c
                       endif               
c
       return
       end
c*************************************************************************
      subroutine twofft(data1,data2,fft1,fft2,n)
c
c      Direct or inverse FFT of a real function, two lines
c      are processed at the same time.
c
c*************************************************************************
      implicit none
      integer n,j,n2
      real*8 data1(n),data2(n)  
      complex*16 fft1(n),fft2(n),h1,h2,c1,c2
      c1=dcmplx(0.5d0,0.0d0)
      c2=dcmplx(0.0d0,-0.5d0)
      do 11 j=1,n
      fft1(j)=dcmplx(data1(j),data2(j))
11    continue
      call four1(fft1,n,1)
c
      fft2(1)=dcmplx(dimag(fft1(1)),0.0d0)
      fft1(1)=dcmplx(dreal(fft1(1)),0.0d0)
c
      n2=n+2
      do 12 j=2,n/2+1
        h1=c1*(fft1(j)+dconjg(fft1(n2-j)))
        h2=c2*(fft1(j)-dconjg(fft1(n2-j)))
      fft1(j)=h1
      fft2(j)=h2
      fft1(n2-j)=dconjg(h1)
      fft2(n2-j)=dconjg(h2)
12    continue
c
      return
      end
c*************************************************************************
      subroutine four1(data,nn,isigne)
c*************************************************************************
      implicit none
      real*8 wr,wi,wpr,wpi,wtemp,theta
      real*8 data(*),tempr,tempi
      integer i,j,n,nn,m,mmax,istep,isigne
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        go to 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isigne*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=dble(wr)*data(j)-dble(wi)*data(j+1)
            tempi=dble(wr)*data(j+1)+dble(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      go to 2
      endif
      return
      end
c*************************************************************************
       subroutine gaussquad(data,poids,yy,isigne,ipar,fct)
c
c     Gauss quadrature
c
c*************************************************************************
      implicit none
      include 'in_HarmSpher'
c
      character*8 fct
c
      real*8 data(nlat,nlong)
      real*8 poids(nlat)
      complex*16 yy(lmmax)
      integer isigne,ipar,i,j
c
      if (fct.eq.   'ylm') 
     & call gaussquad1(data,poids,yy,isigne,ipar)
      if (fct.eq. 'dtylm')
     & call gaussquad2(data,poids,yy,isigne,ipar)
      if (fct.eq. 'dpylm') 
     & call gaussquad3(data,poids,yy,isigne,ipar)
c
       return
       end
c*************************************************************************
      subroutine vit_spec(datat,datap,poids,yys,yyt,isigne)
c
c      expansion (isigne=1) of a vector field with components
c      datat and datap into poloidal and toroidal parts yys and yyt
c      The opposite transform is performed if (isigne=-1)
c
c*************************************************************************
      implicit none
      include 'in_HarmSpher'
       character *8 fct
c
      real*8 poids(nlat)
      real*8 datat(nlat,nlong),datap(nlat,nlong),aux_spat(nlat,nlong)
      complex*16 yys(lmmax),yyt(lmmax),aux_spec(lmmax)
      complex*16 im
      integer ipar,isigne,l,m,ia,indx,j,i
                            
      im=(0.d0,1.d0)   

      if(isigne.eq.1) then
c........................................................................
c......FFT...............................................................
c........................................................................
c
      call fourier(datat,isigne)
      call fourier(datap,isigne)
c

c........................................................................
c......poloidal expansion isigne=1........................................
c........................................................................                     
c
      ipar=-1
      fct='dtylm'                    
      call gaussquad  (datat,poids,aux_spec,isigne,ipar,fct)               
c
      ipar=1
      fct='dpylm'
      call gaussquad  (datap,poids,   yys,  isigne,ipar,fct)               
c
      do 1 l=1,lmax
            do 1 m=0,l
               ia=indx(l,m)
               yys(ia)=(aux_spec(ia)+im*yys(ia))/dble(l)/dble(l+1.)
1     continue 
c........................................................................                     
c......toroidal expansion isigne=1........................................
c........................................................................                     
      ipar=-1         
      fct='dtylm'                    
      call gaussquad  (datap,poids,aux_spec,isigne,ipar,fct)               
c
      ipar=1
      fct='dpylm'                    
      call gaussquad  (datat,poids,  yyt   ,isigne,ipar,fct)               
c
      do 3 l=1,lmax
            do 3 m=0,l
               ia=indx(l,m)
               yyt(ia)=-(aux_spec(ia)-im*yyt(ia))/dble(l)/dble(l+1.)
3     continue             
                      endif
c
       if(isigne.eq.-1) then
c........................................................................
c......poloidal reconstruction isigne=-1..................................
c........................................................................                     
      ipar=-1                    
      fct='dtylm'                    
      call gaussquad  (datat,poids,  yys  ,isigne,ipar,fct)               
c
      ipar=1
      do 4 l=0,lmax
            do 4 m=0,l
               ia=indx(l,m)
               aux_spec(ia)=im*yys(ia)
4     continue             
      fct='dpylm'                    
      call gaussquad  (datap,poids,aux_spec,isigne,ipar,fct) 
c
c........................................................................
c......toroidal reconstruction isigne=-1..................................
c........................................................................                     
      ipar=-1
      fct='dtylm'                    
      call gaussquad  (aux_spat,poids,  yyt  ,isigne,ipar,fct)               
c
      do j=1,nlong
           do i=1,nlat
                datap(i,j)=-datap(i,j)-aux_spat(i,j)
           enddo
      enddo   
      ipar=1
      do 5 l=0,lmax
            do 5 m=0,l
               ia=indx(l,m)
               aux_spec(ia)=im*yyt(ia)
5     continue             
      fct='dpylm'                    
      call gaussquad  (aux_spat,poids,aux_spec,isigne,ipar,fct)               
c
      do j=1,nlong
           do i=1,nlat
               datat(i,j)=datat(i,j)-aux_spat(i,j)
           enddo
      enddo
c
c........................................................................
c......FFT...............................................................
c........................................................................
      call fourier(datat,isigne) 
      call fourier(datap,isigne) 
                      endif

      return 
      end
c*************************************************************************
       subroutine gaussquad1(data,poids,yy,isigne,ipar)
c
c     Gauss quadrature with Plm
c
c*************************************************************************
      implicit none
      include 'in_HarmSpher'
c
      real*8 xx,pi,x1,x2,yy1,yy2
      real*8 data(nlat,nlong)
      real*8 poids(nlat),d1,d2,d3,d4
      complex*16 yy(lmmax),xfct
      integer co,par,isigne,ipar,mr,m,l,indx,lm,i,ik,ikp1,kp1,ni,k,km
c
      pi=4.d0*datan2(1.d0,1.d0)     
c
       mr=nlat/2 
c
      if (isigne.eq.1) then
c
       do 1 m=0,lmax
          ik=2*m+1
          ikp1=ik+1
          do 2 l=m,lmax
          lm=indx(l,m)
          yy(lm)=dcmplx(0.d0,0.d0)
          yy1=0.d0
          yy2=0.d0
c
                      co=par(l,m)*ipar
                      do 31 i=1,mr 
                      ni=nlat+1-i
c
       xfct=ylm(lm,i)/4.d0/pi
c
                        xx=xfct*poids(i)
                        x1=data(i,ik)+dble(co)*data(ni,ik)
                        x2=data(i,ikp1)+dble(co)*data(ni,ikp1)
                      yy1=yy1+xx*x1
                      yy2=yy2+xx*x2
31                    continue
        yy(lm)=dcmplx(dble(yy1),dble(yy2))
        if(m.ne.0) yy(lm)=yy(lm)*dsqrt(2.d0)
 2    continue
 1    continue
c
                     endif
c
       if(isigne.eq.-1) then
c
       do 10 i=1,mr
         ni=nlat+1-i
         do 20 k=1,nlong,2
          km=(k+1)/2
           kp1=k+1
           m=km-1
            d1=0.d0
            d2=0.d0
            d3=0.d0
            d4=0.d0
            if(m.gt.lmax) go to 40
        do 30 l=max0(lmin,m),lmax
             lm=indx(l,m) 
             co=par(l,m)*ipar 
c
        xfct=    ylm(lm,i)
c                                  
             if(m.ne.0) xfct=xfct/dsqrt(2.d0)
c
              x1= dreal(yy(lm))*xfct
              x2=dimag(yy(lm))*xfct
c 
             d1=d1+x1
             d2=d2+x2
             d3=d3+x1*dble(co)
             d4=d4+x2*dble(co)
 30     continue
 40     continue
             data(i,k)=d1
             data(i,kp1)=d2
             data(ni,k)=d3
             data(ni,kp1)=d4
 20     continue
 10     continue
c
                      endif
c
       return
       end
c*************************************************************************
       subroutine gaussquad2(data,poids,yy,isigne,ipar)
c
c     Gauss quadrature with d Plm/ d teta
c
c*************************************************************************
      implicit none
      include 'in_HarmSpher'
c
      real*8 xx,pi,x1,x2,yy1,yy2
      real*8 data(nlat,nlong)
      real*8 poids(nlat),d1,d2,d3,d4
      complex*16 yy(lmmax),xfct
      integer co,par,isigne,ipar,mr,m,indx,ik,ikp1,l,lm,i,ni,k,km,kp1
c
      pi=4.d0*datan2(1.d0,1.d0)           
c
       mr=nlat/2 
c
      if (isigne.eq.1) then
c
       do 1 m=0,lmax
          ik=2*m+1
          ikp1=ik+1
          do 2 l=m,lmax
          lm=indx(l,m)
          yy(lm)=dcmplx(0.d0,0.d0)
          yy1=0.d0
          yy2=0.d0
c
                      co=par(l,m)*ipar
                      do 31 i=1,mr 
                      ni=nlat+1-i
c
          xfct=  dtylm(lm,i)/4.d0/pi
c
                        xx=xfct*poids(i)
                        x1=data(i,ik)+dble(co)*data(ni,ik)
                        x2=data(i,ikp1)+dble(co)*data(ni,ikp1)
                      yy1=yy1+xx*x1
                      yy2=yy2+xx*x2
31                    continue
        yy(lm)=dcmplx(dble(yy1),dble(yy2))
        if(m.ne.0) yy(lm)=yy(lm)*dsqrt(2.d0)
 2    continue
 1    continue
c
                     endif
c
       if(isigne.eq.-1) then
c
       do 10 i=1,mr
         ni=nlat+1-i
         do 20 k=1,nlong,2
          km=(k+1)/2
           kp1=k+1
           m=km-1
            d1=0.d0
            d2=0.d0
            d3=0.d0
            d4=0.d0
            if(m.gt.lmax) go to 40
        do 30 l=max0(lmin,m),lmax
             lm=indx(l,m) 
             co=par(l,m)*ipar 
c
        xfct=  dtylm(lm,i)
c                                  
             if(m.ne.0) xfct=xfct/dsqrt(2.d0)
c
              x1= dreal(yy(lm))*xfct
              x2=dimag(yy(lm))*xfct
c 
             d1=d1+x1
             d2=d2+x2
             d3=d3+x1*dble(co)
             d4=d4+x2*dble(co)
 30     continue
 40     continue
             data(i,k)=d1
             data(i,kp1)=d2
             data(ni,k)=d3
             data(ni,kp1)=d4
 20     continue
 10     continue
c
                      endif
c
       return
       end
c******************************************************************
       subroutine gaussquad3(data,poids,yy,isigne,ipar)
c
c     Gauss quadrature with 1/sin teta d Plm/ d phi
c
c*************************************************************************
      implicit none
      include 'in_HarmSpher'
c
      real*8 xx,pi,x1,x2,yy1,yy2
      real*8 data(nlat,nlong)
      real*8 poids(nlat),d1,d2,d3,d4
      complex*16 yy(lmmax),xfct
      integer co,par,isigne,ipar,m,ikp1,ik,mr,l,i,ni,lm,indx,k,km,kp1
c
      pi=4.d0*datan2(1.d0,1.d0)
c
       mr=nlat/2 
c
      if (isigne.eq.1) then
c
       do 1 m=0,lmax
          ik=2*m+1
          ikp1=ik+1
          do 2 l=m,lmax
          lm=indx(l,m)
          yy(lm)=dcmplx(0.d0,0.d0)
          yy1=0.d0
          yy2=0.d0
c
                      co=par(l,m)*ipar
                      do 31 i=1,mr 
                      ni=nlat+1-i
c
          xfct=  dpylm(lm,i)/4.d0/pi
c
                        xx=xfct*poids(i)
                        x1=data(i,ik)+dble(co)*data(ni,ik)
                        x2=data(i,ikp1)+dble(co)*data(ni,ikp1)
                      yy1=yy1+xx*x1
                      yy2=yy2+xx*x2
31                    continue
        yy(lm)=dcmplx(dble(yy1),dble(yy2))
        if(m.ne.0) yy(lm)=yy(lm)*dsqrt(2.d0)
 2    continue
 1    continue
c
                     endif
c
       if(isigne.eq.-1) then
c
       do 10 i=1,mr
         ni=nlat+1-i
         do 20 k=1,nlong,2
          km=(k+1)/2
           kp1=k+1
           m=km-1
            d1=0.d0
            d2=0.d0
            d3=0.d0
            d4=0.d0
            if(m.gt.lmax) go to 40
        do 30 l=max0(lmin,m),lmax
             lm=indx(l,m) 
             co=par(l,m)*ipar 
c
        xfct=  dpylm(lm,i)
c                                  
             if(m.ne.0) xfct=xfct/dsqrt(2.d0)
c
              x1= dreal(yy(lm))*xfct
              x2=dimag(yy(lm))*xfct
c 
             d1=d1+x1
             d2=d2+x2
             d3=d3+x1*dble(co)
             d4=d4+x2*dble(co)
 30     continue
 40     continue
             data(i,k)=d1
             data(i,kp1)=d2
             data(ni,k)=d3
             data(ni,kp1)=d4
 20     continue
 10     continue
c
                      endif
c
       return
       end
c*************************************************************************
        subroutine expdelta(colat,long,func,yy) 
c
c     Express in spherical harmonics the delta function
c     delta(colat,long)*func.
c     colat and long are the colatitude and longitude
c     in degrees
c
c*************************************************************************

      include 'in_HarmSpher'
c
       real*8 colat,long,func,x,plm(lmmax),teta,phi
       complex*16 yy(lmmax)
       real*8 zn,pi,cp,ss,cc,sp,c,s
       integer lm,l,indx,m
c
        pi=4.d0*datan2(1.d0,1.d0)                   
c
        teta=colat*pi/180.0d0
        phi=long*pi/180.0d0
        x=dcos(teta)
c
        call norme
        call legendre(plm,x)
c
        do lm=1,lmmax
          zn=znorm(lm)
          plm(lm)=plm(lm)*zn
        enddo
c         
        do l=0,lmax
        lm=indx(l,0)
        yy(lm)=func*plm(lm)*dcmplx(1.d0,0.d0)
        enddo
c
        cp=dcos(phi)
        sp=dsin(phi)
        c=1.
        s=0.
c
        do m=1,lmax
        cc=c*cp-s*sp
        ss=s*cp+c*sp
        do l=m,lmax
        lm=indx(l,m)
        yy(lm)=func*plm(lm)*cmplx(cc,ss)*sqrt(2.)
        enddo
        c=cc
        s=ss
        enddo
c
        return
        end

c*************************************************************************
       subroutine spec(yy,ll,sp,sptotal)
c
c      Computes spectra
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher' 
c                   
       complex*16 yy(lmmax)
       real*8 sp(0:lmax),sptotal
       integer ll,l,m,lm,indx
c                
       sptotal=0.d0
       do l=0,ll            
       sp(l)=0.d0
        do m=0,l 
        lm=indx(l,m)
        sp(l)=sp(l)+cdabs(yy(lm))**2
        enddo
        sptotal=sptotal+sp(l) 
       enddo                  
c
       sptotal=dsqrt(sptotal)
       do l=0,ll
       sp(l)=dsqrt(sp(l))
       enddo
c                   
       return
       end
c*************************************************************************
       subroutine correl(yy1,yy2,ll,co,cototal)
c
c      computes correlations
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher'                    
c
       complex*16 yy1(lmmax),yy2(lmmax)
       real*8 co(0:lmax),cototal,xxt,yyt,xyt,xx,xy,yy
       integer lm,l,m,ll,indx
c
       cototal=0.d0
       xxt=0.d0
       yyt=0.d0
       xyt=0.d0
       do l=0,ll            
       co(l)=0.d0
       xx=0.d0
       yy=0.d0
       xy=0.d0
        do m=0,l 
        lm=indx(l,m) 
        xx=xx+cdabs(yy1(lm))**2
        yy=yy+cdabs(yy2(lm))**2
        xy=xy+dreal(yy1(lm))*dreal(yy2(lm))+dimag(yy1(lm))*
     .   dimag(yy2(lm))
        enddo
       xxt=xxt+xx
       yyt=yyt+yy
       xyt=xyt+xy
       if(xx*yy.eq.0.) then
                       co(l)=0.d0
                       else
                       co(l)=xy/dsqrt(xx*yy)
                       endif
       enddo                  
c
       if(xxt*yyt.eq.0.) then
                       cototal=0.d0
                       else
                       cototal=xyt/dsqrt(xxt*yyt)
                       endif
c
       return
       end
c*************************************************************************
c*************************************************************************
       subroutine harm(pplm,x) 
c
c      Values de plm for all degrees and orders lower than lmax at x
c
c*************************************************************************
       implicit none
       include 'in_HarmSpher'
       real*8 pplm(lmmax),x 
       real*8 xx,somx2,ff,x0,xx0,x1,x2
       integer l,m,lm,mp1mp1,indx
c
               xx=x
               xx0=1.d0
               pplm(1)=xx0
               somx2=dsqrt(1.d0-xx*xx)
          do 2 m=0,lmax-1
                  x0=xx0
                  xx0=x0*somx2*dsqrt(dble(2*m+3)/dble(2*m+2))
                  mp1mp1=indx(m+1,m+1)
                  pplm(mp1mp1)=xx0
                  x1=xx*dsqrt(dble(2*m+3))*x0
                  pplm(mp1mp1-1)=x1
             do 3 l=m+2,lmax
                     x2=xx*dsqrt(dble(2*l-1)/dble(l-m-1))*x1-
     &               dsqrt(dble((l+m-1)/dble(2*l-3)))*x0
                    x2=x2*dsqrt(dble((2*l+1)*(l-m-1))/dble((l-m)*(l+m)))
                     lm=indx(l,m)
                     pplm(lm)=x2
                     x0=x1
                     x1=x2
 3           continue
 2         continue
       return
       end


