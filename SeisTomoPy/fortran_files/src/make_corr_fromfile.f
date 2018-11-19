c      This program is based on the spectral routines from Yanick Ricard
c      It :
c        - reads two 'xy' type spatial gmt files (Debayle and Ritsema models)
c        - interpolate the xy file on a Gauss grid 
c        - expand the spatial file in spherical harmonic coefficients.
c          The harmonics are real (cosine and sine) and fully normalized.
c        - computes the spectal of the two files 
c        - computes correlations by spherical harmonic degree 
c
       include 'in_HarmSpher'
c    
c       implicit real*8 (a-h,o-z)
       parameter(nx=179,ny=360)
       integer llmax,xxx,xxxx,i,j,ll,l,ind,mm,
     + tomo1,para1,depth1,tomo2,para2,depth2
       real*8 racin(nlat),poids(nlat)
       real*8 fichier1(nx,ny),fichier2(nx,ny),fichier3(nx,ny),
     + fichier4(nx,ny),cor(0:lmax)
       real*8 ft1(nlat,nlong),ftt1(nlat,nlong)
       real*8 ft2(nlat,nlong),ftt2(nlat,nlong)
       real*8 ft3(nlat,nlong),ftt3(nlat,nlong)
       real*8 ft4(nlat,nlong),ftt4(nlat,nlong)
       complex*16 yy1(lmmax),yyd1(lmmax)
       complex*16 yy2(lmmax),yyd2(lmmax)
       complex*16 yy3(lmmax),yyd3(lmmax)
       complex*16 yy4(lmmax),yyd4(lmmax)
       real*8 zn,rien,sptotal,pi,coefR,coefI
       character*15 namemodel1,namepara1
       character*15 namemodel2,namepara2
       character*4 toto,namedepth1,namedepth2
       character*300 filemoi1,filemoi2,filemoi3,filemoi4
       real*8 phi,theta,cost,sint,splitf1,splitfd
       real*8 splitf2
       real*8 splitf3
       real*8 splitf4,sp(0:lmax)
       real*8 x(2*40+1),dx(2*40+1)
       real*8 rsple
       external rsple
       character(len=255) homedir,tampon

       call getenv("HOME", tampon)
       homedir=trim(adjustl(tampon))//'/SeisTomoPy_files/'
c
       pi=dacos(-1.d0)
c
c   read filename
       read(*,'(1000a)') filemoi1
       read(*,'(1000a)') filemoi3
       read(*,*) llmax
       filemoi2=trim(adjustl(homedir))//
     + "output_files_corr/corr_fromfile.xy"

c   Work on Debayle's model___________________________________
c   1.read
      if (filemoi1.eq.filemoi3) then
      open(1,file=filemoi1)  
      else 
      open(1,file=filemoi1)  
      open(2,file=filemoi3)   
      endif
      do i=0,359
      do j=-89,89
      if (filemoi1.eq.filemoi3) then
         read(1,*) rien,rien,fichier1(j+90,i+1)   
         fichier2(j+90,i+1)=fichier1(j+90,i+1)   
      else
         read(1,*) rien,rien,fichier1(j+90,i+1)   
         read(2,*) rien,rien,fichier2(j+90,i+1)  
      endif  
      enddo
      enddo
      if (filemoi1.eq.filemoi3) then
      close(1)
      else
      close(1)
      close(2)
      endif

c   2.interpol on Gauss grid
      call gauss(racin,poids)
      call poly(racin)
      call norme
      call interpol(fichier1,nx,ny,ft1,nlat,nlong,racin,1)
      call interpol(fichier2,nx,ny,ft2,nlat,nlong,racin,1)
c   4. Expand in spher. harmonics
       call spat_spec(ft1,poids,yy1,1)
       call spat_spec(ft2,poids,yy2,1)
       do l=0,llmax
       do m=0,l
       lm=indx(l,m)
       yyd1(lm)=yy1(lm)
       yyd2(lm)=yy2(lm)
       enddo
       enddo

c   5. Compute correlation
       call correl(yy1,yy2,llmax,cor,cort)
       open(3,file=filemoi2)
       do i=1,llmax
       write(3,*) i,cor(i)
       enddo
       close(3)

       end

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
       

