c $Id: boxprg.f,v 1.7 1997/08/25 08:17:13 weber Exp $
ccccccccccccccccccc
c
c	Last Update: 04/24/97
c
ccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine bptinit(ibox)
c
c 
c     Unit     : Initis all the particles setted by the bpt command
c     Author   : Mathias Brandstetter
c     Date     : 01/26/96
c     Update   : 04/24/97
c     Version  : 1.06
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit none

c Includs
	include 'coms.f'
	include 'comres.f'
	include 'boxinc.f'
	include 'options.f'

c Var
c	Zaehler, BoxNr., ,Spin, 
	integer i,ibox,fchg,getspin


c	RandomNumbergenerator, Particlemass, deacy times
	real*8 ranf,massit,dectim

c	p=Impuls, angle distribution
	real*8 P,cost,sint,phi	
c	dummy var
	real*8 dummy


	ecm=0.d0
	ebeam=0.d0

c  main program


c loop over all particles
	do 42 i=npart+1,npart+bptpart(ibox)
c		configuration space 
		r0(i)=0.d0	
            	rx(i)=lboxhalbe*(1-2*ranf(0)) 
            	ry(i)=lboxhalbe*(1-2*ranf(0)) 
            	rz(i)=lboxhalbe*(1-2*ranf(0))            	

c		Isospin and Ityp
		iso3(i)=bptiso3(ibox)
		ityp(i)=bptityp(ibox)

c set baryon and   meson numbers	        
                if(abs(ityp(i)).le.maxbar) then
                  nbar=nbar+1
                else
                  nmes=nmes+1
                endif	

c		charge
		charge(i)=fchg(iso3(i),ityp(i))
c		massarray
		fmass(i)=massit(ityp(i))
c		Spin
		spin(i)=getspin(ityp(i),-1)
c		decaytime
		dectime(i)=dectim(i,1)
	
42	continue
c	End of loop



	if (edensflag.le.0) then
c	homogenious momentum distribution, randomly distributed
c	max momentum is a parameter
		do 45 i=npart+1,npart+bptpart(ibox)
                        P=bptpmax(ibox)*ranf(0)**(1./3.)
			cost = 1.-2.*ranf(0)
	         	sint = sqrt(1.-cost**2)
	         	phi = 2.*Pi*ranf(0)
	         	px(i) = P*sint*cos(phi)
	         	py(i) = P*sint*sin(phi)
	         	pz(i) = P*cost
			call setonshell(i)
45		continue
	elseif (edensflag.ge.1) then

c	energiedensity

		dummy=0.d0

c loop over all particles
		do 60 i=npart+1,npart+bptpart(ibox)
c			write(*,*) 'dummy= ',dummy
			P=bptpmax(ibox)/bptpart(ibox)*ranf(0)**(1./3.)

			cost = 1.-2.*ranf(0)
        	 	sint = sqrt(1.-cost**2)
         		phi = 2.*Pi*ranf(0)
 
 
c different initialisations

			if (para.eq.0) then
         		
c Boxmode         		         		
				if (i.eq.1) write(*,*) 'Boxmode'
         		px(i) = P*sint*cos(phi)
         		py(i) = P*sint*sin(phi)
         		pz(i) = P*cost
			
			elseif(para.eq.1) then
c stream over stream          
				if (i.eq.1) write(*,*) 'streammode'	
                        px(i) = 0.d0
                        py(i) = 0.d0
                        pz(i) = P*(-1.d0)**i

			elseif(para.eq.2) then

c slab on slab
				if (i.eq.1) write(*,*) 'slabmode'
			
			px(i)=0.d0
			py(i)=0.d0
			if (rx(i).gt.0) then
				pz(i)=(-1.0d0)*P
			else
				pz(i)=P	

			endif
         		
         		endif

c                       Write(*,*)'0 ', rx(i), px(i), py(i),pz(i)

60		continue
	endif

c	Addition der Teilchenzahl zur Gesamtzahl
	npart=npart+bptpart(ibox)
c debug
	Write(*,*) 'Particles = ',npart
cc


        return

        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc MB'96



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      	function swapi(x,dx)
c Version: 1.0
c Date:    08/26/96
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 swapi, x, dx  

      	swapi = x
1 	if (swapi.lt.-dx) then
          swapi = swapi + 2.0d0*dx
c          write(6,*) 'particle swapped', swap/dx
          goto 1
      	end if  
2 	if (swapi.gt.dx) then
          swapi = swapi - 2.0d0*dx
c           write(6,*) 'particle swapped', swap/dx
          goto 2
      	end if  
      	end     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccJK

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	Function Energie(alpha,max)
c
c     Unit     : calculat the energy
c     Author   : Mathias Brandstetter
c     Date     : 08/27/96
c     Update   : 08/27/96
c     Version  : 0.99 b 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



	implicit none

	include 'coms.f'
	integer i

c arrays mit variabler laenge deklarieren.
	real*8 alpha,max
	real*8 Energie, E

	E=0
	Do 42 i=1,npart
           E=E+sqrt((alpha**2)*(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))+
     &      fmass(i)**2)
42	continue
	Energie=E-max	
	Return
	End
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc MB'96

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	Function Regula(me)
c
c     Unit     : Searches for the zero of the function
c     Author   : Mathias Brandstetter
c     Date     : 08/27/96
c     Update   : 08/27/96
c     Version  : 0.99 b 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	implicit none
	
	include 'coms.f'

c Variablendklaration
	real*8 Regula,xn,xu,x0,me,Energie
	integer i
	real*8 E1,E2,E3

c Initialisierung
       xu=0.d0
       x0=3.d0     
c Hauptprogramm
	Write(*,*) 'Regula is running!'
	i=0
10     Continue
	i=i+1
	E1=Energie(x0,me)        
	E2=Energie(xu,me)	

        xn=x0-(E1*(x0-xu))/(E1-E2)
        E3=Energie(xn,me)
        IF ((E2*E3).LE.0) then
            x0=xn
        else
       	    xu=xn
        EndIF
c debug        
c	Write(*,*) 'Energie1 = ',E1, ' Energie 2 = ',E2
c	Write(*,*) 'Energie3 = ',E3
c	Write(*,*) 'x = ',x0,' y = ',xu	,' xn= ',xn
c	Write(*,*) 'Energie = ', Energie(x0,0.) 
c c
        IF ((ABS(x0-xu).GE.1.D-12).and.(i.le.1000).and.(
     &      ((E3.ge.1.D-12).or.(-E3.ge.1.D-12)))) goto 10 
	Regula=xn	
	End
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc MB'96




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function wallcoll(ipart,wall)
c
c 
c     Unit     : Collisions with an imaginary wall
c     Author   : Mathias Brandstetter
c     Date     : 09/25/96
c     Update   : 11/29/96
c     Version  : 0.5
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	implicit none

c includes
	include 'coms.f'
	include 'boxinc.f'
	include 'options.f'
ccccccccc

c var
	real*8  betax,betay,betaz
	real*8  ty,tz
	real*8	tn
cc

	integer wall,ipart
	integer wally,wallz
	
	
c Mainprogram
c init the variables
	wall=0
	tn=0

c velocity 
	betax=px(ipart)/p0(ipart)
	betay=py(ipart)/p0(ipart)
	betaz=pz(ipart)/p0(ipart)

c check which wall is reached next and sort by impact time
c wall presents the wall and tn the time

	if (betax.lt.0) then
                wall=-4
                tn=(-lboxhalbe-rx(ipart))/(-max(-betax,1.d-13))
        else
                wall=-1
                tn=((lboxhalbe-rx(ipart))/max(betax,1.d-13))
        endif

		
	if (betay.lt.0) then
                wally=-5
                ty=(-lboxhalbe-ry(ipart))/(-max(-betay,1.d-13))
        else
                wally=-2
                ty=((lboxhalbe-ry(ipart))/max(betay,1.d-13))
        endif

	if(ty.lt.tn) then 
	   tn=ty
	   wall=wally
        endif
	 
	if (betaz.lt.0) then
                wallz=-6
                tz=(-lboxhalbe-rz(ipart))/(-max(-betaz,1.d-13))
        else
                wallz=-3
                tz=((lboxhalbe-rz(ipart))/max(betaz,1.d-13))
        endif

	if(tz.lt.tn) then
	   tn=tz
	   wall=wallz
	endif


c sets the time for the earliest wall collision
	wallcoll=tn
	return

	End

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc MB'96



ccccccccccccccccccccccccccccccccccccccc
cc Output 

ccccccccccccccc open Output cccccccccccccccccc

	subroutine MBopen

c
c     Unit     : fort.20 is createt                   
c     Author   : Mathias Brandstetter
c     Date     : 06/19/97
c     Update   : 06/19/97
c     Version  : 1.00   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        character*77 file20


cmb open a file for output               
		 call getenv('ftn20',file20)
                if (file20(1:4).ne.'    ') then
                    OPEN(UNIT=20,FILE=file20)
                endif
               
        Write(20,*) '! Number of Particles versus Time from MB96'

	end	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	subroutine MBHeader

c     Unit     : Header for fort.20 is created
c     Author   : Mathias Brandstetter
c     Date     : 06/19/97
c     Update   : 06/19/97
c     Version  : 1.00   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      include 'coms.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'


c       Header fuer ftn20.dat
cmb for each new event print the event number surrounded by
c two empty lines  
	write(20,*)
	write(20,*) '!Event: ',event,' NrEv: ',nevents,' randomseed: ',
     &	ranseed
	write(20,*) '!tsteps: ',nsteps
        Write(20,*) '!t     nucs   delta  	pi 	K+  	K-',
     &  ' Sigma/Lambda 	eta	omega	rho'
	write(20,*)
c       End Header	 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	end	

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine MBoutput

c     Unit     : Write the multiplicity for a few particles
c		 to fort.20 vs. the time.
c     Author   : Mathias Brandstetter
c     Date     : 06/19/97
c     Update   : 06/19/97
c     Version  : 1.00   
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      include 'coms.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'

c Var
	Integer nucs, delta, pionen, kaonen, kaonenm
	Integer sigma, lambda, sumsl, eta, omega, rhos
	Integer counter,ihad
	real*8 etot,rap
c        real*8  mbmass, ratio1, ratio2, wurzs, mwurzs

76     format(2X,11f12.6)


c
c	here a new output file is createt!
c	Integer nucs, delta, pionen, kaonen, kaonenm
c	Integer sigma, lambda, sumsl, eta, omega, rhos

c	inital values of the counter variables
	nucs   	= 0
	delta 	= 0
	pionen 	= 0
	kaonen 	= 0
	kaonenm	= 0
	sigma	= 0
	lambda	= 0
	sumsl	= 0
	eta	= 0
	omega	= 0
	rhos	= 0
ccccccccccccccccccccccc
c       loop over all particles, different species are counted
c	at the moment we count
c	nucleons ---> ityp = 1
c	deltas	 ---> ityp = 17 (Delta 1232)
c	pions    ---> ityp = 101 
c	kaons	 ---> ityp = 106
c       kaons+   ---> ityp =-106
c	sigma	 ---> ityp = 40
c	lambda	 ---> ityp = 27
c	eta	 ---> ityp = 102
c	omega	 ---> ityp = 103
c	rhos     ---> ityp = 104 (rho 770 )

c number of formed hadrons
	ihad=0

	do 99 counter=1,npart

	   i=counter
	   etot=sqrt(fmass(i)*fmass(i)+px(i)*px(i)
     &	  +py(i)*py(i)+pz(i)*pz(i))
	   rap=0.5*log( (etot+pz(i)) / (etot-pz(i)) )


csab count only formed hadrons around mid-rapidity
	   if(xtotfac(counter).gt..99d0.and.abs(rap).le.1d0) then
	      ihad=ihad+1
	      if (ityp(counter).eq.1)    nucs   = nucs   +1 
	      if (ityp(counter).eq.17)   delta  = delta  +1
	      if (ityp(counter).eq.101)  pionen = pionen +1
	      if (ityp(counter).eq.106)  kaonen = kaonen +1		
	      if (ityp(counter).eq.-106) kaonenm = kaonenm +1		
	      if (ityp(counter).eq.40)   sigma = sigma +1
	      if (ityp(counter).eq.27)   lambda = lambda +1
	      if (ityp(counter).eq.102)   eta = eta +1
	      if (ityp(counter).eq.103)   omega = omega +1
	      if (ityp(counter).eq.104)   rhos = rhos + 1
	   endif
 99	continue
	sumsl = sigma + lambda
	
c now write the values in to a file
	write(20,76) time, dble(ihad),
     &    dble(nucs), dble(delta), dble(pionen)
     &  , dble(kaonen), dble(kaonenm), dble(sumsl), dble(eta)
     &  , dble(omega), dble(rhos) 	


	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
