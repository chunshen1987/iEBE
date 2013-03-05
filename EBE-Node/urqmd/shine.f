c $Id: shine.f,v 1.4 1998/06/15 13:35:33 weber Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      integer function channel(itp,izp,mass,gam)
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 07/02/97
c  Revision : 1.0
c
cinput  itp    : ityp of particle to check
cinput  izp    : iso3-spin component
cinput  mass   : mass of actual particle
coutput channel: Integer to classify the dilepton channel
coutput gam    : Width of particle
c
c This function checks if a particle has a dilepton channel. In case
c it returns a number to characterize the channel and the width of 
c this particle. If {\rm mass} is .le. 0 the width is not calculated
c to save runtime if the width is not needed.
c 
c Coding for dilepton channels:
c \begin{displaymath}
c \begin{array}{rl}
c channel &   mode \\
c  1& \omega\to ll \\
c  2& \rho^0\to ll \\
c  3& \phi\to ll \\
c 11& \pi^0\to \gamma ll \\
c 12& \eta\to \gamma ll \\
c 13& \omega\to \pi^0 ll \\
c 14& \eta'\to \gamma ll \\
c 15& \Delta\to N ll \\
c 21& \pi^+\pi^-\to ll \\
c 22& KK\to ll \\
c \end{array}
c \end{displaymath}
c 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
	include 'comres.f'
	integer itp,izp
	real*8 gam,mass,widit,fwidth,piwidth,etwidth
	logical pitrue,deltrue,et1true
c some widths are set to zero in urqmd, but here we need them
	parameter(piwidth=8.02d-9, etwidth=1.18d-6)

	parameter(pitrue=.true.)
c	parameter(pitrue=.false.)

	parameter(deltrue=.true.)
c	parameter(deltrue=.false.)

	parameter(et1true=.true.)
c	parameter(et1true=.false.)

ccccc
	channel=0
	gam=1d0

c speed up the lookup
	if(itp.lt.mindel.or.itp.gt.itphi.or.iabs(izp).gt.1)return
	if(.not.deltrue.and.(itp.lt.pimeson.or.izp.ne.0))return

	if(itp.eq.itome.and.izp.eq.0)then
	  channel=1
c	  if(mass.gt.0d0)gam=widit(itp)
	  if(mass.gt.0d0)gam=fwidth(itp,izp,mass)
	elseif(itp.eq.itrho.and.izp.eq.0)then
	  channel=2
c	  if(mass.gt.0d0)gam=widit(itp)
	  if(mass.gt.0d0)gam=fwidth(itp,izp,mass)
	elseif(itp.eq.itphi.and.izp.eq.0)then
	  channel=3
	  if(mass.gt.0d0)gam=widit(itp)
c	  if(mass.gt.0d0)gam=fwidth(itp,izp,mass)
	elseif(itp.eq.pimeson.and.izp.eq.0.and.pitrue)then
	  channel=11
	  if(mass.gt.0d0)gam=piwidth
	elseif(itp.eq.iteta.and.izp.eq.0)then
	  channel=12
	  if(mass.gt.0d0)gam=etwidth
	elseif(itp.eq.itetapr.and.izp.eq.0.and.et1true)then
	  channel=14 
	  if(mass.gt.0d0)gam=widit(itp)
c	  gam=fwidth(itp,izp,mass)
	elseif(itp.eq.mindel.and.iabs(izp).eq.1.and.deltrue)then
	  channel=15 
cc	  if(mass.gt.0d0)gam=widit(itp)
	  if(mass.gt.0d0)gam=fwidth(itp,izp,mass)
	endif
	if(mass.lt.0d0)gam=1d0	
	
	return
	end
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine preplep
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 07/02/97
c  Revision : 1.0
c
cinput none
coutput via common block lepcom
c
c This routine is called before each propagation step.
c When first called it does some initialization of the dilepton routines.
c It also writes a header to the dielpton output file at each new event.
c It scans for vector mesons that may shine a dilepton, counts them
c and stores their densities.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'coms.f'
c      include 'comres.f'
      include 'options.f'
      integer i,id,nee,ich,channel,index(nmax),thisevent
      real*8 bardens,dens1(nmax),gtot,stot,ptsigtot
      character*77 file18
      logical firstcall,muflag
      save firstcall,thisevent
      common /lepcom/ dens1,nee,index,muflag
      data firstcall /.true./
      data thisevent /0/
51    format(A1,i6,4e12.4,f7.3)
52    format(A1,i6,5e12.4,f7.3)

      if(bf18)return

c initialize
      if(firstcall)then
        write(*,*)'Dileptons are written to unit 18'
        firstcall=.false.
        muflag=CTOption(43).ne.0
c open file 18  
        call getenv('ftn18',file18)
        if(file18(1:4).ne.'    ') then
          open(UNIT=18,FILE=file18,STATUS='unknown',FORM='FORMATTED')
        endif
      endif                                  

c new event->write header
        if(thisevent.ne.event)then
	    do 11 i=1,nmax
	      dens1(i)=0d0
	      index(i)=0
11	    continue
          thisevent=event
          stot=ptsigtot()
          write(18,51)'!',thisevent,ebeam,
     6	dble(CTOption(27)),stot
        endif

c scan particles for sources of radiated leptons
      nee=0      
c this loop can start with nbar+1 to save runtime
      do 1 id=1,npart
         ich=channel(ityp(id),iso3(id),-1d0,gtot)

c store what has been found (only vector mesons shine continously)
         if(ich.gt.0.and.ich.lt.4)then
           nee=nee+1
           index(nee)=id
           dens1(nee)=bardens(id)
         endif           
1     Continue

      Return
      END
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine shinelep(detet,tltime)
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 07/02/97
c  Revision : 1.0
c
cinput  detet  : timestep
cinput  tltime : computational frame time at the beginning of the timestep
c
c {\rm shinelep} calculates the direct and dalitz decays and writes 
c them to output. Of course between a call of preplep and shinelep the order
c of the particles may not be changed, because shinelep only treats particles
c stored by preplep.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'coms.f'
      include 'options.f'
      integer i,j,nee,ich,id,channel,index(nmax),Maxint,itp,idt
      real*8 dens1(nmax),dens,Me,Mmu,M,dete,detet,tltime
      real*8 dx,qm,em2,bardens,ranf,ddlzd2,srtmax,ltime
      real*8 eg,pg(4),betax,betay,betaz,gg
      real*8 pp0,ppx,ppy,ppz,ppt,pgt,gtot,geeconst
	real*8 gmuconst
      logical decflg,muflag
      common /lepcom/ dens1,nee,index,muflag
      parameter(Me=511d-6)
      parameter(Mmu=106d-3)
c MonteCarlo Intervall
	parameter(srtmax=1d0)
51    format(A1,i6,4e12.4,f7.3)
52    format(A1,i6,5e12.4,f7.3)

      if(bf18)return

	decflg=.false.
	ltime=tltime
	
	goto 123

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	entry declep(idt)
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 07/02/97
c  Revision : 0.1 beta - untested and uncompleted
c
cinput  idt  : itype of decaying particle
c
c This is an entry to shinelep, that is called if a particle
c decays. It checks wether it is a dilepton source and in case flows into the
c rest of {\rm shinelep}. For a decay the settings are gtot=gamma and
c dete=1, while a shinestep assumes gtot=1d0 and dete=dt.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if(bf18)return

c check if particle has a dilepton channel and determine its width
	ich=channel(ityp(idt),iso3(idt),fmass(idt),gtot)
	if(ich.eq.0)return

c found something
	nee=1
	decflg=.true.
	dens1(1)=0d0
	index(1)=idt
	ltime=acttime
c warning
c	if(ich.lt.4)write(*,*)'(W) lep: decay forced of ',ityp(idt),
c     6  'at time', ltime
	

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c here begins the loop over all known sources
123	continue

	if(nee.eq.0)return

c loop over all shining sources
      Do 11 i=1,nee
c assign needed quantities
        id=index(i)
        itp=ityp(id)
        M=fmass(id)
        pp0=p0(id)
        ppx=px(id)+ffermpx(id)
        ppy=py(id)+ffermpy(id)
        ppz=pz(id)+ffermpz(id)
        ppt=sqrt(ppx*ppx+ppy*ppy)
c? take the density as the average of before/after cascstep
c? (this is however no very good approximation unless dete=small, gamma->1)
c? only the densities of vector mesons need to be calculated
c? still quite unclear..
	  if(dens1(i).gt.1d-3)then           !dens1 is known for vec-mesons
          dens=.5d0*(dens1(i)+bardens(id))
        elseif(ich.lt.4)then               !dens1 of vec-meson was small/unknown
          dens=bardens(id)
        else                               !all pure dalitz decays
          dens=0d0
        endif
c assign ich, dt and width (if decflg=true dete,ich and gtot are already known)
	  if(.not.decflg)then
c channel returns gtot=1d0 when called with neg. masses
          ich=channel(itp,0,-1d0,gtot)
	    gg=pp0/M
	    dete=detet/hqc/gg
	  else
	    dete=1d0
        endif

        
c\\\\\\\\\|||||||||///////////    
c    shine vector mesons

c direct decays
        if(ich.lt.10)then
c integration step of the decay width \Gamma_{ll}	
	     if(muflag)then
		 dx=dete*gmuconst(itp,M,Mmu)/gtot
	     else
           	 dx=dete*geeconst(itp,M,Me)/gtot
           endif
c medium dependent output..
	     if(dens.gt.1d-3)then
             write(18,51)' ',ich,M,dx,ppt,ppz,dens
c5x             write(18,52)' ',ich,M,dx,ppt,ppz,ltime,dens
           else
             write(18,51)' ',ich,M,dx,ppt,ppz
c5x             write(18,52)' ',ich,M,dx,ppt,ppz,ltime
           endif
        endif !direct decay


c dalitz decays

c set maxint: the number of Monte-Carlo integrations
      if(ich.eq.1)then
c the omega meson may both produce shine dileptons and dalitz pairs
c omega-dalitz decays: how treat offshell-omegas?
c        if(dabs(M-massit(itp)).lt.widit(itp))then
cc	    write(*,*)ich,itp,M,massit(itp),widit(itp)
          ich=13
          maxint=2        !can be small because the omega is sampled 
c        endif
	elseif(ich.eq.11.or.ich.eq.15)then 
	  maxint=5          !usually there are many pi nd Deltas
	else
	  maxint=20
      endif 

      if((ich.gt.10.and.decflg).or.ich.eq.13)then
          betax=ppx/pp0
          betay=ppy/pp0
          betaz=ppz/pp0
c each meson yields up to Maxint dileptons
c first the decay is treated in the local rest frame
c then a boost according to beta is performed
          do 200 j=1,Maxint
c qm=mass of the virtual gamma* -- selected randomly between 0 and srtmax
             qm=ranf(0)*srtmax          
c compute the dalitz rate and get mass of the stable partner em2
	       if(muflag)then
              dx=dete*ddlzd2(itp,qm,M,em2,Mmu)/dble(Maxint)/gtot*srtmax
             else
              dx=dete*ddlzd2(itp,qm,M,em2,Me)/dble(Maxint)/gtot*srtmax
	       endif
             if(dx.gt.0d0)then
c energy and 4-momemtum of the \gamma^*=e+e- choosen randomly via rboost
c eg=energy of the gamma 
               eg=(M*M+qm*qm-em2*em2)/(2d0*M)
c now the mass qm and the enrgy eg of the g* are known->choose momentum
               call rboost(pg,eg,qm,betax,betay,betaz)
		   pgt=sqrt(pg(1)*pg(1)+pg(2)*pg(2))
		   write(18,51)' ',ich,qm,dx,pgt,pg(3)
c5x		   write(18,52)' ',ich,qm,dx,pgt,pg(3),ltime
             endif !dx.gt.0

200       continue      !Maxint
        endif           !dalitz
11    Continue          !shining loop
      Return
      End

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine collep(id1,id2,sigs,sqs,tltime)
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 19/06/97
c  Revision : 0.1 beta - untested and uncompleted
c
cinput   id1   :  index of first particle
cinput   id2   :  index of second scattering particle
cinput   sqs   :  sqrt(s) of the collison
cinput  sigs   :  cross section
cinput tltime  :  collision time
c
c Writes the QED part of the pi+pi-->e+e- and the k+k-->ee
c cross section to the dilepton
c file. The Formfactor has to be multiplied (and maybe modified)
c in the analysis routines. This will also be the place, where bremsstrahlung,
c Drell-Yan and other processes could be included.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
	include 'coms.f'
	include 'comres.f'
	include 'options.f'
      integer id1,id2,it1,it2,iz1,iz2,ich,nee,index(nmax)
      logical muflag
      real*8 tltime,sqs,sigs
      real*8 cc,alpha,Me,Mmu,emlep,manni
	real*8 M,pp0,ppx,ppy,ppz,ppt,dens,dx,dens1(nmax)
      common /lepcom/ dens1,nee,index,muflag      
      PARAMETER(cc=0.38937966d0,alpha=1d0/137d0)           
      parameter(Me=511d-6)
      parameter(Mmu=106d-3)
      real*8 bardens,massit
51    format(A1,i6,4e12.4,f7.3)
52    format(A1,i6,5e12.4,f7.3)

	if(bf18)return
	
c set itypes and iso3's
	it1=ityp(id1)
	it2=ityp(id2)
	iz1=iso3(id1)
	iz2=iso3(id2)

c check if dilepton production is possible
c this is the analogous to integer function channel
c in later revisions it should be excluded
	if(it1.eq.pimeson.and.it2.eq.pimeson.and.
     n      min(iz1,iz2).eq.-2.and.max(iz1,iz2).eq.2)then
	   ich=21
	   manni=massit(pimeson)
	elseif(min(it1,it2).eq.-itkaon.and.max(it1,it2).eq.itkaon.and.
     n      min(iz1,iz2).eq.-1.and.max(iz1,iz2).eq.1)then
	   ich=22
	   manni=massit(itkaon)
	else
	   return
	endif
	

c set kinematics
	M=sqs
	pp0=p0(id1)+p0(id2)
	ppx=px(id1)+px(id2)+ffermpx(id1)+ffermpx(id2)
	ppy=py(id1)+py(id2)+ffermpy(id1)+ffermpy(id2)
	ppz=pz(id1)+pz(id2)+ffermpz(id1)+ffermpz(id2)
	ppt=sqrt(ppx*ppx+ppy*ppy)


c take density as average
	dens=0.5d0*(bardens(id1)+bardens(id2))

	if(muflag)then
	  emlep=Mmu
	else
	  emlep=Me
	endif
c calculate the i1+i2->e+e- cross section divided by the Formfactor
	dx=4d0*pi/3d0*(alpha/M)**2*sqrt(1d0-(2d0*emlep/M)**2)
     6    *(1d0+2d0*(emlep/M)**2)*sqrt(1d0-(2d0*manni/M)**2)*cc/sigs

c write
	if(dens.gt.1d-3)then
        write(18,51)' ',ich,M,dx,ppt,ppz,dens
c5x     write(18,52)' ',ich,M,dx,ppt,ppz,tltime,dens
      else
        write(18,51)' ',ich,M,dx,ppx,ppy,ppz
c5x     write(18,52)' ',ich,M,dx,ppx,ppy,ppz,tltime
      endif
	
	return
	end	

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function DDLZD2(ITP,EMDL,em1,em2,emlep)        
c  Unit     : Propagation, Output
c  Author   : Luke Winckelmann, Christoph Ernst
c  Date     : 07/02/97
c  Revision : 0.1 beta - untested and uncompleted
c
cinput  itp   : itype of decaying resonsnce
cinput  emdl  : rest mass M of this resonance
cinput  em1   : mass of the virtual photon
cinput  emlep : lepton mass ($m_e$ or $m_\mu$ resp.)
coutput em2   : mass of the stable particle
c
c This function returns the differential probablity for Dalitz decays.
c
c Important local quantities:
c  G  : The width of the decay into a virtual photon (PR D54)
c  T  : Kinematic expression
c  F  : Formfactor acc.  to Landsberg (PL 128 (1985) 301) et.al.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C... output: em2 mass of second decay product for this dalitz-dec 
      implicit none
      include 'comres.f'
      real*8 emdl,em1,em2,T,emlep,emn,empi,massit
	real*8 g,F,g0,alpha,pi,elambd
	real*8 em,dglldm,sqrtla,elam,x,y,z,cut
      integer itp
      PARAMETER(PI=3.1415926535d0,ALPHA=1d0/137d0)
      parameter (cut=1d-6)
      ELAM(X,Y,Z)=X**2+Y**2+Z**2-2.*(X*Y+X*Z+Y*Z)                               
      SQRTLA(X,Y,Z)=SQRT(ELAM(X**2,Y**2,Z**2))                                  

      DDLZD2=0D0            

c snap off undefined input
	if(emdl.lt.1d-20)return
	if(itp.gt.107)return

c Dalitz decays of mesons
	if(itp.gt.minmes)then
C...ETA(549)
          IF(Itp.EQ.iteta)THEN
            EM2=0D0
            IF(EMDL.GT.em1-cut)RETURN
            ELAMBD=0.72d0
            G=0.463d-6
            T=2d0*(1d0-(emdl/em1)**2)**3
            F=max(0d0,elambd*elambd/(elambd*elambd-emdl*emdl))
C...OMEGA(783)       
          ELSE IF(itp.EQ.itome)THEN       
            EM2=massit(pimeson) 
            IF(EMDL.GT.EM1-EM2-cut)RETURN          
ce            ELAMBD=0.66d0          
            ELAMBD=0.65d0          
ce offshell omegas may produce high mass pairs that are not rejected
ce by the em1-em2 condition->throw such pairs away here
	      if(emdl.ge.elambd-cut)return
            G=0.717d-3          
            T=((1d0+emdl**2/(em1**2-em2**2))**2-
     &        (2d0*em1*emdl/(em1**2-em2**2))**2)**1.5
            F=max(0d0,elambd*elambd/(elambd*elambd-emdl*emdl))
c cassing-like
c		F=max(0d0,elambd**4/((elambd**2-emdl**2)**2
c     6       +(75d-3*elambd)**2))

c...eta'(958)
          ELSE IF(itp.EQ.itetapr)THEN       
            EM2=0d0       
            IF(EMDL.GT.EM1-cut)RETURN
            G=4.26d-6
            T=2d0*(1d0-(emdl/em1)**2)**3
c            elambd=0.77d0
c            F=max(0d0,elambd*elambd/(elambd*elambd-emdl*emdl))
c cassing-like
c            elambd=0.75d0
c		F=max(0d0,elambd**4/( (elambd**2-emdl**2)**2
c     6       +(0.14d0*elambd)**2) )
c Set Formfactor=1 somewhow crude but works better than other two
	      F=1d0  
c...pi0
          ELSE IF(itp.EQ.pimeson)THEN       
            EM2=0d0       
            IF(EMDL.GT.EM1-cut)RETURN
            G=7.7d-9
            T=2d0*(1d0-(emdl/em1)**2)**3          
            F=1d0+5.5d0*emdl*emdl       
          else
            return
          END IF            
c Dalitz decay rate d\Gamma/dM acc. to Landsberg et.al.
         ddlzd2=2d0/3d0*ALPHA/PI/EMDL*G*F*F*T
     *    *sqrt((1d0-(2d0*EMLEP/EMDL)**2))*(1d0+2d0*(EMLEP/EMDL)**2)
	   return
	endif !itp.ge.minmes


ccccccccccccccccccccccccccccccccccccccccc
c...Delta(1232)
c .. (according to Luke..)
        IF(itp.eq.mindel)THEN                                                
	em=em1
	emn=massit(nucleon)
	empi=massit(pimeson)
          IF(EM.LE.EMN+EMPI.OR.EM.LE.EMN+EMDL)RETURN 

          F=(-1.5)*(EM+EMN)/(EMN*((EMN+EM)**2-EMDL**2))                           
CKO F.-FCTG=3./(1+(0.725/1.232)**2)                                             
          EM2=EMN                                                               
          G=2.5d0                                                                
C...VDM GVM CHENG O'NEIL-FIT FOR M**2 DEPENDENCE P. 111                         
C     ONLY THE DOMINANT MAGNETIC ISO-VECTOR PART:                               
C         G=G*( 4./(1.-EMDL**2/0.77**2) - 1. )/3.                               
C         G=G*( FPIMOD(EMDL)-FPIMOD(0) + 1. )                                   
          T=-(1./6.*EMDL**6-1./2.*EMDL**4*EM**2-1./3.*EMDL**4                   
     . *EM*EMN-1./2.*EMDL**4*EMN**2+1./2.*EMDL**2*EM**4+2./                     
     . 3.*EMDL**2*EM**3*EMN+1./3.*EMDL**2*EM**2*EMN**2+2./3.                    
     . *EMDL**2*EM*EMN**3+1./2.*EMDL**2*EMN**4-1./6.*EM**6-                     
     . 1./3.*EM**5*EMN+1./6.*EM**4*EMN**2+2./3.*EM**3*EMN**3                    
     . +1./6.*EM**2*EMN**4-1./3.*EM*EMN**5-1./6.*EMN**6)
c          P=SQRTLA(EM,EMN,EMPI)/(2.*EM)
c          GTOT=0.47*P**3/((EMPI**2+0.6*P**2))
        G0=4.*PI*ALPHA*(G*F)**2/(16.*PI*EM1**3)*
     *    SQRTLA(EM1,EM2,EMDL)*T                                                
          DGLLDM=0.666667*ALPHA/PI/EMDL*sqrt(1.-(2.*EMLEP/EMDL)**2)*
     *      (1+2.*(EMLEP/EMDL)**2)*G0
        DDLZD2=MAX(0D0,DGLLDM) !/gtot
      RETURN         
        END IF
      END            


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine rboost(p,e,em,bex,bey,bez)   
c  Unit     : Propagation, Output
c  Author   : L. Winckelmann, Christoph Ernst
c  Date     : 07/02/97
c  Revision : 0.9 beta 
c
c
cinput   e           : 4. component of p before boost
cinput   em          : mass
cinput   bex,bey,bez : transformation betas
coutput  p(4)        : boosted 4 impuls
c
c Boost of particle with energy e, mass em and randomly
c selected $p_1, p_2, p_3$ according to $\vec \beta$. The particle is 
c first rotated in its local restframe and then boosted.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      real*8 pi,e,em,bex,bey,bez,the,phi
      real*8 pa
      parameter(pi=3.1415926535d0) 
      real*8 p(4)
      real*8 dga,dbep,dgabep,dbe2,ranf    
      
c.. randomly selected angles
      the=ranf(0)*pi       
      phi=ranf(0)*pi*2d0
  
c..determine the absolute of the space components
      pa=sqrt(e**2-em**2)  

      p(4)=e 
c...rotate (typically from z axis to direction theta,phi)      
      p(1)=pa*sin(the)*cos(phi)    
      p(2)=pa*sin(the)*sin(phi)    
      p(3)=pa*cos(the)       

c... now the boost-part:
c...lorentz boost (typically from rest to momentum/energy=beta)
      dbe2=bex*bex+bey*bey+bez*bez
      if(dbe2.gt.1d-20)then   
         dga=1d0/dsqrt(1d0-dbe2)   
         dbep=bex*p(1)+bey*p(2)+bez*p(3) 
         dgabep=dga*(dga*dbep/(1d0+dga)+p(4)) 

         p(1)=p(1)+dgabep*bex      
         p(2)=p(2)+dgabep*bey      
         p(3)=p(3)+dgabep*bez      
         p(4)=dga*(p(4)+dbep) 
      endif  
      return 
      end    

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function geeconst(itp,emdl,emlep)
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 07/02/97
c  Revision : 0.1 beta - untested and uncompleted
c
cinput  itp    :  itype of resonance
cinput emdl    :  restmass of resonance
cinput emlep   :  lepton mass
c
c Returns the width for the decay {\tt itp} $\to e^+e^-$ (ref. PR D54).
c The width is set constant. Any mass dependance has to be multiplied in
c the analysis routines.
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'
      integer itp
      real*8 emdl,emlep
      real*8 gee,massit

c cut on double lepton restmass
	if(emdl.lt.2d0*emlep)then
	  geeconst=0d0
	  return
	endif

      IF(ITP.EQ.itrho)THEN           
	  gee=6.77d-6
c cut on goldstone mass
        if(emdl.lt.2d0*massit(pimeson))then
          geeconst=0d0
          return
        endif
      ELSE IF(ITP.EQ.itome)THEN       
	  gee=0.6d-6
      ELSE IF(ITP.EQ.itphi)THEN           
	  gee=1.37d-6 
      ELSE           
        geeconst=0d0
        return
      END IF         

      geeconst=gee
c*SQRT(1d0-(2d0*EMLEP/EMDL)**2)*(1d0+2d0*(EMLEP/EMDL)**2)

      RETURN         
      END            

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      real*8 function gmuconst(ITP,EMDL,emlep)
c  Unit     : Propagation, Output
c  Author   : Christoph Ernst
c  Date     : 07/02/97
c  Revision : 0.1 beta - untested and uncompleted
c
cinput  itp    :  itype of resonance
cinput emdl    :  restmass of resonance
cinput emlep   :  lepton mass 
c
c Returns the width for the decay {\tt itp} $\to \mu^+\mu^-$ (ref. PR D54).
c The width is set constant. Any mass dependance has to be multiplied in
c the analysis routines.
c 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none
      include 'comres.f'
      integer itp
      real*8 emdl,emlep
      real*8 gmu,massit

	if(emdl.lt.2d0*emlep)then
	  gmuconst=0d0
	  return
	endif
	
      IF(ITP.EQ.itrho)THEN           
	  gmu=6.9d-6
        if(emdl.lt.2d0*massit(pimeson))gmu=0d0
c!!!unknown mumu-width of omega! choose ee value
      ELSE IF(ITP.EQ.itome)THEN       
	  gmu=0.6d-6
      ELSE IF(ITP.EQ.itphi)THEN           
	  gmu=1.0984d-6 
      ELSE           
        gmuconst=0d0
        return
      END IF         

      gmuconst=gmu*SQRT(1d0-(2d0*EMLEP/EMDL)**2)
     6   *(1d0+2d0*(EMLEP/EMDL)**2)

      RETURN         
      END            


