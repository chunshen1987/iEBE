c $Id: cascinit.f,v 1.5 1998/06/15 13:35:08 weber Exp $
      subroutine cascinit(ZZ,AA,nucleus)
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'

      integer Z,A,i,getspin,fchg,ZZ,AA,antia, j,k, nrjc,izcnt,nucleus
      real*8 R,P,R2,ranf,sint,phi,nucrad,cost,drr
      real*8 xcm,ycm,zcm,pxcm,pycm,pzcm,densnorm,add,avd
c
c...law: mass correction according to binding energy
      real*8 meff

CJK averaged density of one Gaussian in units of central value
      avd = 0.5**1.5 

      pxcm=0.d0
      pycm=0.d0
      pzcm=0.d0

      A=abs(AA)
      Z=abs(ZZ)
      antia=A/AA

      PT_AA(nucleus)=A

      if(A.gt.AAmax) then
         write(6,*)'***(E): Mass ',A,' exceeds initialization limit!'
         write(6,*)'        -> increase parameter AAmax=',AAmax
         write(6,*)'           in include-file inputs.f '
         write(6,*)'        -> uqmd aborting ... '
         stop
      endif

      if (A.eq.1) then  ! proton or neutron
	i = 1
        PT_iso3(i,nucleus) = antia*(-1+2*Z)
	PT_ityp(i,nucleus) = antia*1
	PT_spin(i,nucleus) = 1
        PT_dectime(i,nucleus)=1.d31
        PT_charge(i,nucleus)=fchg(PT_iso3(i,nucleus),
     &                            PT_ityp(i,nucleus))
	PT_fmass(i,nucleus) = EMNUC
        PT_p0(i,nucleus)=sqrt(PT_fmass(i,nucleus)**2)
        return
      end if

C
C Initialisation of A nucleons in configuration space
C CAUTION: any changes have to be reported to JK, since this part has
C          severe implications for other parts of the code as well
C
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0
      R2 = nucrad(A)
      nrjc = 0
      do 1 j=1,A
        PT_dectime(j,nucleus)=1.d31
111     nrjc = nrjc+1
        R = R2*ranf(0)**(1./3.)
        cost = 1.-2.*ranf(0)
        sint = sqrt(1.-cost**2)
        phi = 2.*Pi*ranf(0)
        PT_r0(j,nucleus) = 0.d0
        PT_rx(j,nucleus) = R*sint*cos(phi)
        PT_ry(j,nucleus) = R*sint*sin(phi)
        PT_rz(j,nucleus) = R*cost
        do 11 k=1,j-1
          drr=(PT_rx(j,nucleus)-PT_rx(k,nucleus))**2
     &       +(PT_ry(j,nucleus)-PT_ry(k,nucleus))**2
     &       +(PT_rz(j,nucleus)-PT_rz(k,nucleus))**2
          if (drr.lt.2.6.and.nrjc.lt.CTParam(46)) goto 111
11      continue         
        xcm=xcm+PT_rx(j,nucleus)
        ycm=ycm+PT_ry(j,nucleus)
        zcm=zcm+PT_rz(j,nucleus)
1     continue
      if (nrjc.ge.CTParam(46)) then
         write(6,*)'*** warning: initialisation corrupt '
      end if

      xcm = xcm/dble(A)
      ycm = ycm/dble(A)
      zcm = zcm/dble(A)
      do 13 j=1,A
        PT_rx(j,nucleus) = PT_rx(j,nucleus)-xcm 
        PT_ry(j,nucleus) = PT_ry(j,nucleus)-ycm 
        PT_rz(j,nucleus) = PT_rz(j,nucleus)-zcm 
        PT_rho(j,nucleus) = avd
13    continue

C local proton density in nucleus A,Z
      do 14 j=1,Z
        do 15 k=j+1,Z
          drr=(PT_rx(j,nucleus)-PT_rx(k,nucleus))**2
     &       +(PT_ry(j,nucleus)-PT_ry(k,nucleus))**2
     &       +(PT_rz(j,nucleus)-PT_rz(k,nucleus))**2
          add=exp(-(2.0*gw*drr))
          PT_rho(j,nucleus) = PT_rho(j,nucleus)+add
          PT_rho(k,nucleus) = PT_rho(k,nucleus)+add
15      continue
14    continue

C local neutron density in nucleus A,Z
      do 16 j=Z+1,A
        do 17 k=j+1,A
          drr=(PT_rx(j,nucleus)-PT_rx(k,nucleus))**2
     &       +(PT_ry(j,nucleus)-PT_ry(k,nucleus))**2
     &       +(PT_rz(j,nucleus)-PT_rz(k,nucleus))**2
          add=exp(-(2.0*gw*drr))
          PT_rho(j,nucleus) = PT_rho(j,nucleus)+add
          PT_rho(k,nucleus) = PT_rho(k,nucleus)+add
17      continue
16    continue

      densnorm = (2.0*gw/pi)**1.5
      do 18 j=1,A
        PT_rho(j,nucleus) = PT_rho(j,nucleus)*densnorm
        PT_pmax(j,nucleus) = hqc*(3.0*pi*pi*PT_rho(j,nucleus))**(1./3.)
18    continue

      izcnt=0
      do 12 j=1,A
         P = PT_pmax(j,nucleus)*ranf(0)**(1./3.)
cdebug,sab use old init
c         P = 0.27*ranf(0)**(1./3.)
         cost = 1.-2.*ranf(0)
         sint = sqrt(1.-cost**2)
         phi = 2.*Pi*ranf(0)
         PT_px(j,nucleus) = P*sint*cos(phi)
         PT_py(j,nucleus) = P*sint*sin(phi)
         PT_pz(j,nucleus) = P*cost
         pxcm=pxcm+PT_px(j,nucleus)
         pycm=pycm+PT_py(j,nucleus)
         pzcm=pzcm+PT_pz(j,nucleus)
         if (j.le.Z) then
            PT_iso3(j,nucleus)= 1*antia
            PT_charge(j,nucleus)=1*antia
         else
            PT_iso3(j,nucleus)= -(1*antia)
            PT_charge(j,nucleus)=0
         endif

         PT_spin(j,nucleus) = getspin(1,-1)
         PT_ityp(j,nucleus) = 1*antia
12    continue

c perform CM-correction
      pxcm=pxcm/A
      pycm=pycm/A
      pzcm=pzcm/A
c      write(6,*)'pcm/A:',pxcm,pycm,pzcm
      do 2 i=1,A
c         bx = pxcm/ecm
c         by = pycm/PT_p0(i,nucleus)
c         bz = pzcm/PT_p0(i,nucleus)
c         call rotbos(0d0,0d0,bx,by,bz,
c     @                 PT_px(i,nucleus),PT_py(i,nucleus),
c     @                 PT_pz(i,nucleus),PT_p0(i,nucleus))
         PT_px(i,nucleus)=PT_px(i,nucleus)-pxcm
         PT_py(i,nucleus)=PT_py(i,nucleus)-pycm
         PT_pz(i,nucleus)=PT_pz(i,nucleus)-pzcm
c...law: effective masses for initial energy corr. (CTOption(11).eq.0)
         r=sqrt(PT_rx(i,nucleus)**2+PT_ry(i,nucleus)**2
     &         +PT_rz(i,nucleus)**2)
         p=sqrt(PT_px(i,nucleus)**2+PT_py(i,nucleus)**2
     &         +PT_pz(i,nucleus)**2)
         PT_fmass(i,nucleus) = meff(z,a,r,p)
         PT_p0(i,nucleus)=sqrt(PT_px(i,nucleus)**2+PT_py(i,nucleus)**2
     &                    +PT_pz(i,nucleus)**2+PT_fmass(i,nucleus)**2)
 2    continue
c end of CM-correction

cdebug,sab
c      write(6,*)'***(M) Nucleus ',nucleus,' initialized in CASCINIT ***'

      return
      end


      function nucrad(AA)
      implicit none
      real*8 nucrad, r_0, r_e
      integer A,AA
      include 'coms.f'

      A=abs(AA)
c  root mean square radius of nucleus of mass A 
c r_0 corresponding to rho0
      r_0 = (0.75/pi/rho0)**(1./3.) 
c  (Mayer-Kuckuck: "Kernphysik", Teubner-Verlag)
      r_e = r_0*A**(1./3.) !-0.89*A**(-1./3.)

c substract gaussian tails, for distributing centroids correctly
c      nucrad = sqrt(r_e**2-1.25/gw)
      nucrad = r_0*(0.5*(a + (a**(1./3.)-1.)**3.))**(1./3.)
      return
      end

      subroutine boostnuc(i1,i2,pin,b,dst)
      implicit none
      include 'coms.f'
      include 'options.f'
      integer i1,i2,i
      real*8 b,dst,ei,ti
      real*8 pin,beta,gamma

      do 1 i=i1,i2

      beta = pin/sqrt(pin**2+fmass(i)**2)
      gamma = 1.d0/sqrt(1.d0-beta**2)

c  Gallilei-Trafo in x-direction (impact parameter)
c  projectile hits at POSITIVE x
         rx(i) = rx(i) + b
c  distance between nuclei: projectile at NEGATIVE z for dst < 0
         if(CTOption(23).eq.0)then
           ti = r0(i)
c           r0(i) = gamma*(r0(i) - beta*rz(i))
           rz(i) = rz(i)/gamma+dst/gamma
         else
          rz(i) = (rz(i) + dst)
         end if


         Ei = p0(i)
         p0(i) = gamma*(p0(i) - beta*pz(i))
         pz(i) = gamma*(pz(i) - beta*Ei) 

 1    continue
      return
      end

      real*8 function meff(z,a,r,p)
c mean binding energy of a nucleon in a nucleus according to weizaecker
      implicit none
      include 'options.f'
      real*8 av,as,ac,aa,ap,mdef,r,p,e,EMNUC
      integer z,a
      parameter (av=0.01587,as=0.01834,ac=0.00071)
      parameter (aa=0.09286,ap=11.46,EMNUC=0.938)
      if(CTOption(11).ne.0.or.a.eq.1)then
        meff=EMNUC
        return
      end if
c...mass defect
      mdef=-(av*A)+as*A**0.66667+ac*z*z/a**0.33333+aa*(z-a/2.)**2/a
c...energy per nucleon = binding energy + nucleon mass 
      e=min(0d0,mdef/a)+EMNUC
      meff=sqrt(e**2-p**2)      
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getnucleus(nucleus,offset)
c     Unit    : Initialization
c     Author  : Steffen A. Bass, Jens Konopka
c     Date    : 10/19/96
c     Revision: 1.0 
c
cinput nucleus : 1=projectile, 2=target
cinput offset  : offset for location of nucleus in particle vectors
c
c output : via common blocks
c 
c This subroutine read in a nucleus which has been initialized
c by {\tt cascinit} and stored in the {\tt PT\_ *(i,nucleus)} arrays.
c The respective nucleus is then rotated randomly in configuration
c and momentum space to yield a new initial state.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'

c local variables
      real*8 eul1,eul2,eul3,ceul1,seul1,ceul2,seul2,ceul3,seul3
      real*8 vecx0,vecy0,vecz0,vecx1,vecy1,vecz1,vecx2
      real*8 vecy2,vecz2
      integer k,nucleus,offset

c functions
      real*8 ranf

csab the following section has been provided by Jens Konopka
csab only the variable names have been changed
c ***
c *** rotation around Euler-angles
c ***
      eul1 = ranf(0)*2.0*pi
      eul2 = ranf(0)*2.0*pi
      eul3 = ranf(0)*2.0*pi
cdebug
c      write(6,*)'Euler angels are ',eul1,eul2,eul3
 
      ceul1 = cos(eul1)
      seul1 = sin(eul1)
      ceul2 = cos(eul2)
      seul2 = sin(eul2)
      ceul3 = cos(eul3)
      seul3 = sin(eul3)

      do 178 k=1,PT_AA(nucleus)

c rotate in configuration space

         vecx0 = PT_rx(k,nucleus)
         vecy0 = PT_ry(k,nucleus)
         vecz0 = PT_rz(k,nucleus)


         vecx1 =   ceul1*vecx0 + seul1*vecy0
         vecy1 = -(seul1*vecx0)+ ceul1*vecy0
         vecz1 =                                     vecz0

         vecx2 =   ceul2*vecx1               + seul2*vecz1
         vecy2 =                       vecy1
         vecz2 = -(seul2*vecx1)              + ceul2*vecz1

         vecx0 =   ceul3*vecx2 + seul3*vecy2
         vecy0 = -(seul3*vecx2)+ ceul3*vecy2
         vecz0 =                                     vecz2

         rx(k+offset) = vecx0 
         ry(k+offset) = vecy0
         rz(k+offset) = vecz0 

c rotate in momentum space

         vecx0 = PT_px(k,nucleus)
         vecy0 = PT_py(k,nucleus)
         vecz0 = PT_pz(k,nucleus)


         vecx1 =   ceul1*vecx0 + seul1*vecy0
         vecy1 = -(seul1*vecx0)+ ceul1*vecy0
         vecz1 =                                     vecz0

         vecx2 =   ceul2*vecx1               + seul2*vecz1
         vecy2 =                       vecy1
         vecz2 = -(seul2*vecx1)              + ceul2*vecz1

         vecx0 =   ceul3*vecx2 + seul3*vecy2
         vecy0 = -(seul3*vecx2) + ceul3*vecy2
         vecz0 =                                     vecz2

         px(k+offset) = vecx0 
         py(k+offset) = vecy0
         pz(k+offset) = vecz0 

csab end of JK part

c initialize the other quantum numbers

         iso3(k+offset)=PT_iso3(k,nucleus)
         ityp(k+offset)=PT_ityp(k,nucleus)
         spin(k+offset)=PT_spin(k,nucleus)
         dectime(k+offset)=PT_dectime(k,nucleus)
         charge(k+offset)=PT_charge(k,nucleus)
         fmass(k+offset)=PT_fmass(k,nucleus)
         r0(k+offset)=PT_r0(k,nucleus)
         p0(k+offset)=PT_p0(k,nucleus)

 178  continue
      
cdebug,sab
c      write(6,*)'***(M) Nucleus ',nucleus,' rotated in GETNUCLEUS ***'

      return
      end

C####C##1#########2#########3#########4#########5#########6#########7##
      real*8 function rnfWSX(AA,zmin,zmax)
c  yields a $x^n$ distributet value for $x$ between mmin and mmax
cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      real*8 zmin,zmax,ranf,a,rr,yf,z,nucrad
      integer AA
      parameter (a=0.54d0)
c      nucrad = 1.128*A**(1./3.)-0.89*A**(-1./3.)
c      rr=1.15d0*AA**(0.3333d0)-0.161d0*AA**(-0.3333d0)
        rr=nucrad(aa)
c      rr=zmax-1d0
      if(zmax.lt.rr)then
        write(6,*)'rnfwsx: maximum radius seems too low'
        stop
      end if

 108  continue 
      z=zmin+ranf(0)*(zmax-zmin)
      yf=z*z/((zmax-zmin)**3)*0.5d0/(1d0+exp(z-rr)/a)
      rnfWSX=z 
      if(yf.gt.1d0)stop'rnfWSX: wrong normalisaton:' 
c write(6,*)'rnfWSX: wrong normalisaton:',yf
      if(yf.gt.ranf(0))    return
      goto 108                        
       
c      return
      end


      subroutine cascinit_old(ZZ,AA)
      implicit none
      include 'coms.f'
      include 'options.f'
      integer Z,A,i,getspin,fchg,ZZ,AA,antia, j,k, nrjc,izcnt
      real*8 R,P,R2,ranf,sint,phi,nucrad,cost,drr
      real*8 xcm,ycm,zcm,pxcm,pycm,pzcm,densnorm,add,avd
      real*8 rho(nmax), pmax(nmax)
c...law: mass correction according to binding energy
      real*8 meff

CJK averaged density of one Gaussian in units of central value
      avd = 0.5**1.5 

      pxcm=0.d0
      pycm=0.d0
      pzcm=0.d0

      A=abs(AA)
      Z=abs(ZZ)
      antia=A/AA

      if (A.eq.1) then  ! proton or neutron
	i = npart+1
        iso3(i) = antia*(-1+2*Z)
	ityp(i) = antia*1
	spin(i) = 1
        dectime(i)=1.d31
        charge(i)=fchg(iso3(i),ityp(i))
	fmass(i) = EMNUC
        call setonshell(i)
        return
      end if

C
C Initialisation of A nucleons in configuration space
C CAUTION: any changes have to be reported to JK, since this part has
C          severe implications for other parts of the code as well
C
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0
      R2 = nucrad(A)
      nrjc = 0
      do 1 j=npart+1,npart+A
        dectime(j)=1.d31
111     nrjc = nrjc+1
        R = R2*ranf(0)**(1./3.)
        cost = 1.-2.*ranf(0)
        sint = sqrt(1.-cost**2)
        phi = 2.*Pi*ranf(0)
        r0(j) = 0.d0
        rx(j) = R*sint*cos(phi)
        ry(j) = R*sint*sin(phi)
        rz(j) = R*cost
        do 11 k=npart+1,j-1
          drr=(rx(j)-rx(k))**2+(ry(j)-ry(k))**2+(rz(j)-rz(k))**2
          if (drr.lt.2.6.and.nrjc.lt.CTParam(46)) goto 111
11      continue         
        xcm=xcm+rx(j)
        ycm=ycm+ry(j)
        zcm=zcm+rz(j)
1     continue
      if (nrjc.ge.CTParam(46)) then
         write(6,*)'*** warning: initialisation corrupt '
      end if

      xcm = xcm/dble(A)
      ycm = ycm/dble(A)
      zcm = zcm/dble(A)
      do 13 j=npart+1,npart+A
        rx(j) = rx(j)-xcm 
        ry(j) = ry(j)-ycm 
        rz(j) = rz(j)-zcm 
        rho(j) = avd
13    continue

C local proton density in nucleus A,Z
      do 14 j=npart+1,npart+Z
        do 15 k=j+1,npart+Z
          drr=(rx(j)-rx(k))**2+(ry(j)-ry(k))**2+(rz(j)-rz(k))**2
          add=exp(-(2.0*gw*drr))
          rho(j) = rho(j)+add
          rho(k) = rho(k)+add
15      continue
14    continue

C local neutron density in nucleus A,Z
      do 16 j=npart+Z+1,npart+A
        do 17 k=j+1,npart+A
          drr=(rx(j)-rx(k))**2+(ry(j)-ry(k))**2+(rz(j)-rz(k))**2
          add=exp(-(2.0*gw*drr))
          rho(j) = rho(j)+add
          rho(k) = rho(k)+add
17      continue
16    continue

      densnorm = (2.0*gw/pi)**1.5
      do 18 j=npart+1,npart+A
        rho(j) = rho(j)*densnorm
        pmax(j) = hqc*(3.0*pi*pi*rho(j))**(1./3.)
18    continue

      izcnt=0
      do 12 j=npart+1,npart+A
         P = pmax(j)*ranf(0)**(1./3.)
cdebug,sab use old init
c         P = 0.27*ranf(0)**(1./3.)
         cost = 1.-2.*ranf(0)
         sint = sqrt(1.-cost**2)
         phi = 2.*Pi*ranf(0)
         px(j) = P*sint*cos(phi)
         py(j) = P*sint*sin(phi)
         pz(j) = P*cost
         pxcm=pxcm+px(j)
         pycm=pycm+py(j)
         pzcm=pzcm+pz(j)
         if (j.le.npart+Z) then
            iso3(j)= 1*antia
            charge(j)=1*antia
         else
            iso3(j)= -(1*antia)
            charge(j)=0
         endif

csab
c alternative sequence to initialize p,n randomly in the tables
c         if(izcnt.lt.Z.and.ranf(0).lt.(dble(Z)/dble(A))) then
c            iso3(j)= 1*antia
c            charge(j)=1*antia
c            izcnt=izcnt+1
c         else
c            if(j-npart-izcnt.le.A-Z) then
c               iso3(j)= -1*antia
c               charge(j)=0
c            else
c               iso3(j)= 1*antia
c               charge(j)=1*antia
c               izcnt=izcnt+1
c            endif
c         endif


c         fmass(j) = EMNUC
         spin(j) = getspin(1,-1)
         ityp(j) = 1*antia
c         call setonshell(j)
12    continue

c perform CM-correction
      pxcm=pxcm/A
      pycm=pycm/A
      pzcm=pzcm/A
c      write(6,*)'pcm/A:',pxcm,pycm,pzcm
c...law: that would be a correction of the first nuclei by the momentum
c of the second after a boost to proj./targ.-system 
c      do 2 i=1,A
      do 2 i=npart+1,npart+A
c         bx = pxcm/ecm
c         by = pycm/p0(i)
c         bz = pzcm/p0(i)
c         call rotbos(0d0,0d0,bx,by,bz,
c     @                 px(i),py(i),pz(i),p0(i))
         px(i)=px(i)-pxcm
         py(i)=py(i)-pycm
         pz(i)=pz(i)-pzcm
c...law: effective masses for initial energy corr. (CTOption(11).eq.0)
         r=sqrt(rx(i)*rx(i)+ry(i)*ry(i)+rz(i)*rz(i))
         p=sqrt(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))
         fmass(i) = meff(z,a,r,p)
c         write(6,*)'cascinit:',i,fmass(i),sqrt(p**2+fmass(i)**2),p
c  ,z,a,p,px(i),py(i),pz(i)
         call setonshell(i)
c         write(6,*) i,p0(i),px(i),py(i),pz(i),sqrt(p0(i)**2-p**2)
 2    continue
c end of CM-correction

      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine setonshell(i)
c     Unit     : Initialization
c     Author   : Markus Hofmann
c     Date     : 03/19/95
c     Revision : 1.0
c     This subroutine set particle i on-shell
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      include 'coms.f'
      integer i
      p0(i) = sqrt(px(i)**2+py(i)**2+pz(i)**2+fmass(i)**2)
      return
      end

