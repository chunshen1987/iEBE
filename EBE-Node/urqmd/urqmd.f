c $Id: urqmd.f,v 1.17 1998/06/15 13:35:36 weber Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      program UrQMD
c	USE DFLIB
c	USE DFPORT, TIME2 => TIME
c
c
c      Authors : The UrQMD collaboration 
c
c     Date    : 03/19/95
c     Revision: 1.1
c
cc    contact address:
cc
cc                     uqmd@th.physik.uni-frankfurt.de
cc
c               
c This is the main module of {\tt uqmd}. It servers as a connection between
c the initialization, the propagation (including the real part of the 
c optical potential) and the collision term (imaginary part of the optical
c potential).
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      include 'coms.f'
	include 'comres.f'
      include 'options.f'
      include 'colltab.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'boxinc.f'
	include 'complot.f'
c
      integer i,j,k,steps,ii,ocharge,ncharge, mc, mp, it1,it2
      integer dummy, timing, tottime
      real*8 sqrts,otime,xdummy,st
      
      real*8 Ekinbar, Ekinmes, ESky2, ESky3,EYuk, ECb, EPau
      common /energies/ Ekinbar, Ekinmes,ESky2,ESky3,
     $                  EYuk,ECb,EPau
      real*8 etot, etotjk

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

 
c
c     numerical/technical initialisation
c
	
      call uinit(0)

      call osc_header

c     profiling
      dummy = timing(0)
      tottime = 0
	call plotinit


cmb ccccccccccccccccccccccccccccccccc
c If you wish a file including the parical 
c multiplicities versus time just uncomment this.
c creats  fort.20.
c	call mbopen
cccccccccccccccccccccccccccccccccc


c
c  Main program
c

c     "dirty" counters for collisions and particles
      mc=0
      mp=0
c
c loop over all events
c
      do 10 event=1,nevents

c     start event here
c    
ccsws
cc  reset parameters for new event 
cc  start flag 
cc  reset time of last plotting to some negative number
      tplot=-30.0
ccsws

c     time is the system time at the BEGINNING of every timestep
      time = 0.0

c     initialize random number generator
c     call auto-seed generator only for first event and if no seed was fixed
      if(.not.firstseed.and.(.not.fixedseed)) then
         ranseed=-(1*abs(ranseed))
         call sseed(ranseed)
      else
         firstseed=.false.
      endif

cdebug
         write(6,*)'event# ',event,ranseed

c
c     initialisation of physics quantities 
c
      call init

cdebug
c      betax=0d0
c      betay=0d0
c      betaz=-13.27d0/sqrt(13.27d0**2+0.938d0**2)
c      do i=1,2
c         call rotbos(0d0,0d0,betax,betay,betaz,
c     &        px(i),py(i),pz(i),p0(i))
c      enddo


c old time if an old fort.14 is used 
      if(CTOption(40).eq.1)time=acttime

c output preparation

c write header to screen
c      call output(6)

c write headers to file
      call output(13)
      call output(14)
      call output(15)
      call output(16)
      if(event.eq.1)call output(17)

cmb cccccccccccccccccccccccccccccccccccc
c Writes the Haeder to file fort.20
c	call mbheader
ccccccccccccccccccccccccccccccccccccc


c for CTOption(4)=1 : output of initialization configuration
      if(CTOption(4).eq.1)call file14out(0)

c     participant/spectator model:
      if(CTOption(28).ne.0) call rmspec(0.5d0*bimp,-(0.5d0*bimp))

c     compute time of output
      otime = outsteps*dtimestep

cdebug
c      call readitin
c      call denscheck
c      call angmom
c      write(6,*) 'initial Etot ', Etot(), EtotJK()
c      write(6,*) 'Ekinbar', Ekinbar
c      write(6,*) 'Ekinmes', Ekinmes
c      write(6,*) 'ESky2  ', ESky2
c      write(6,*) 'ESky3  ', ESky3
c      write(6,*) 'EYuk   ', EYuk
c      write(6,*) 'ECb    ', ECb
c      write(6,*) 'EPau   ', EPau

c reset time step counter
      steps = 0

c profiling
      dummy = timing(1)
      tottime = tottime+dummy
c      write(6,'(i5,2f10.2)') steps,0.01*dble(dummy),0.01*dble(tottime)
c      write(77,'(i5,2f10.2)')steps,0.01*dble(dummy),0.01*dble(tottime)

c      pause
c      call pcheck('initial')

ccsws
cc  set minimum time interval between plotting
       dtplot=(nsteps*dtimestep)/dble(ISTEPS)
ccsws

c  loop over all timesteps

      do 20  steps=1,nsteps
c store coordinates in arrays with *_t
c this is needed for MD type propagation
         if (eos.ne.0) then
            do 23 j=1,npart
               r0_t(j) = r0(j)
               rx_t(j) = rx(j)
               ry_t(j) = ry(j)
               rz_t(j) = rz(j)
 23         continue
         end if

c we are at the beginning of the timestep, set current time (acttime) 
         acttime = time
	   call plotlive
c  option for MD without collision term
         if(CTOption(16).gt.0) goto 103

c  Load collision table with next collisions in current timestep
         call colload
c     check for collisions in time-step, nct = # of collisions in table
         if (nct.gt.0) then

c     entry-point for collision loop in case of full colload after every coll. 
 101        continue
            k = 0
c     normal entry-point for collision loop 
 100        continue
		  call plotlive
c     get next collision
            call getnext(k)

c     exit collision loop if no collisions are left
            if (k.eq.0) goto 102

c  propagate all particles to next collision time
c  store actual time in acttime, propagation time st=cttime(k)-acttime
		call preplep
		st=cttime(k)-acttime
            call cascstep(acttime,st)
c  new actual time (for upcoming collision)
            acttime = cttime(k)
            call shinelep(st,acttime)

c  perform collision 
cdebug
c           write(6,*)'scatter',ctag+1,etot(),cttime(k),npart
c           write(6,*)'scatter',ctag+1,npart,cttime(k),cti1(k),cti2(k)
c           write(6,*)'scatter',ctag+1,npart,cttime(k),ctsqrts(k)
c           call printtab

            if(cti2(k).gt.0.and.
     .           abs(sqrts(cti1(k),cti2(k))-ctsqrts(k)).gt.1d-3)then
               write(6,*)' ***(E) wrong collision update (col) ***'
               write(6,*)cti1(k),cti2(k),
     .              ctsqrts(k),sqrts(cti1(k),cti2(k))
c              stop
            else if(cti2(k).eq.0.and.
     .              abs(fmass(cti1(k))-ctsqrts(k)).gt.1d-3) then
c              ctsqrts(k)=fmass(cti1(k))
               write(6,*)' *** main(W) wrong collision update (decay)'
               write(6,*)ctag,cti1(k),ityp(cti1(k)),dectime(cti1(k)),
     .              fmass(cti1(k)),ctsqrts(k)
c              stop
            endif

cdebug,sab
c            etoto=etot()
            ocharge=charge(cti1(k))
            if(cti2(k).gt.0) ocharge=ocharge+charge(cti2(k))

c     store quantities in local variables for charge conservation check
            it1= ityp(cti1(k))
            if(cti2(k).gt.0)it2= ityp(cti2(k))

c increment "dirty" collision counter
            if(cti2(k).gt.0)then	!scatter
              mc=mc+1
c direct dilepton-production
         call collep(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),cttime(k))
		else				!decay
c the decayed particle might radiate a dilepton
c if it is a shining particle it has already done so
		  if(ityp(cti1(k)).ne.itome.and.ityp(cti1(k)).ne.itrho
     &	 	.and.ityp(cti1(k)).ne.itphi)
     &          call declep(cti1(k))
		endif
c     perform scattering/decay
            call scatter(cti1(k),cti2(k),ctsigtot(k),ctsqrts(k),
     &                   ctcolfluc(k))



cdebug,sab
c            if(dabs(etoto-etot()).gt.1d-5) then
c               write(6,*)'***(E) ',
c     &              'energy conservation violated in scatter coll/dec '
c     &              ,ctag,' delta(E)= ',etoto-etot(),ctsqrts(k)
c               write(6,*)'  please check channel io=',iline,'in make22'
cc               stop
c            endif
          
c
c  update collision table 
c
c     normal update mode
            if(CTOption(17).eq.0) then
               if(nexit.eq.0) then
c     new collision partners for pauli-blocked states (nexit=0)

                  call collupd(cti1(k),1)
c mbbox
c ne <-> gt 
                  if(cti2(k).gt.0) call collupd(cti2(k),1)
               else
                  ncharge=0
c     new collision partners for scattered/produced particles (nexit><0)
                  do 30 i=1,nexit
c     ncharge is used for charge conservation check
                     ncharge=ncharge+charge(inew(i))

                     call collupd(inew(i),1)
 30               continue

cdebug,sab 
c     charge conservation check
                  if(ocharge.ne.ncharge) then
                     write(6,*)'ch-conservation error coll/dec ',ctag
                     write(6,*)'   it1:',it1,'   it2:',it2
                     write(6,*)'   ch:',ocharge,ncharge
                     write(6,*)'cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)'
                     write(6,*)cti1(k),cti2(k),ctsigtot(k),ctsqrts(k)
c                     stop
                  endif
               endif

c     update collisions for partners of annihilated particles
               do 55 ii=1,nsav
                  call collupd(ctsav(ii),1)
 55            continue
               nsav=0

            else ! (CTOption(17).ne.0)
c     full collision load
               call colload
            endif

            if (CTOption(17).eq.0) goto 100
            goto 101

c     this is the point to jump to after all collisions in the timestep
c     have been taken care of
 102        continue

         endif ! (nct.gt.0)

c  After all collisions in the timestep are done, propagate to end of 
c  the timestep.

c     point to jump to in case of MD without collision term
 103     continue

c     increment timestep
         time = time+dtimestep

c  After all collisions in the timestep are done, propagate to end of 
c  the timestep.
	   call preplep
         call cascstep(acttime,time-acttime)
         call shinelep(time-acttime,time)

c     in case of potential interaction, do MD propagation step
         if (eos.ne.0) then

c set initial conditions for MD propagation-step
            do 24 j=1,npart
               r0(j) = r0_t(j)
               rx(j) = rx_t(j)
               ry(j) = ry_t(j)
               rz(j) = rz_t(j)
 24         continue

c now molecular dynamics trajectories
            call proprk(time,dtimestep)

         end if ! (eos.ne.0)

c     perform output if desired
         if(mod(steps,outsteps).eq.0.and.steps.lt.nsteps)then 
            if(CTOption(28).eq.2)call spectrans(otime)
cmb modifiy
c to prevent big output files comment this line
            call file14out(steps)
            call file13out(steps)
cmb ccccccccccccccccccccccccccc
c writes the particle multipicities from a few particles
c in the file: fort.20
c	call mboutput
cccccccccccccccccccccccccccc

cdebug
c           call denscheck
c         write(6,*)'inter Etot   ',Etot(),EtotJK()
c      write(6,*) 'Ekinbar', Ekinbar
c      write(6,*) 'Ekinmes', Ekinmes
c      write(6,*) 'ESky2  ', ESky2  
c      write(6,*) 'ESky3  ', ESky3  
c      write(6,*) 'EYuk   ', EYuk   
c      write(6,*) 'ECb    ', ECb    
c      write(6,*) 'EPau   ', EPau   
c           call angmom

c     profiling
c            dummy = timing(1)
c            tottime = tottime+dummy
c      write(6,'(i5,2f10.2)') steps,0.01*dble(dummy),0.01*dble(tottime)
c      write(77,'(i5,2f10.2)')steps,0.01*dble(dummy),0.01*dble(tottime)

         endif ! output handling

 20   continue ! time step loop


ce
	acttime=time
c optional decay of all unstable particles before final output
c DANGER: pauli-blocked decays are not performed !!!
         if(CTOption(18).eq.0) then
c no do-loop is used because npart changes in loop-structure
            i=0
            nct=0
            actcol=0
c disable Pauli-Blocker for final decays
            CTOption(10)=1
c decay loop structure starts here
 40         continue
            i=i+1

c is particle unstable
            if(dectime(i).lt.1.d30) then
 41            continue
c the decayed particle might radiate a dilepton
		   call declep(i)
c     perform decay
               call scatter(i,0,0.d0,fmass(i),xdummy)
c     backtracing if decay-product is unstable itself
               if(dectime(i).lt.1.d30) goto 41
            endif
c     check next particle
            if(i.lt.npart) goto 40
         endif ! final decay
c final output
c let the remaining particles (pions and etas) decay to dileptons
	do 44 i=1,npart
	  call declep(i)
44	continue

	   if(CTOption(28).eq.2)call spectrans(otime)

cmb ccccccccccccccccccc
c write the particle multiplicities from the last timestep 
c to fort.20
c	call mboutput
ccccccccccccccccccccc	


         call file13out(nsteps)
         call file14out(nsteps)
         call file16out

         call osc_event
cdebug
c         call denscheck
c     profiling
      dummy = timing(1)
      tottime = tottime+dummy
c      write(6,'(i5,2f10.2)') nsteps,0.01*dble(dummy),0.01*dble(tottime)
c      write(77,'(i5,2f10.2)')nsteps,0.01*dble(dummy),0.01*dble(tottime)
c         write(6,*)'final Etot   ',Etot(),EtotJK()
c      write(6,*) 'Ekinbar', Ekinbar
c      write(6,*) 'Ekinmes', Ekinmes
c      write(6,*) 'ESky2  ', ESky2
c      write(6,*) 'ESky3  ', ESky3
c      write(6,*) 'EYuk   ', EYuk
c      write(6,*) 'ECb    ', ECb
c      write(6,*) 'EPau   ', EPau

c      call pcheck('final  ')

c     increment "dirty" particle counter
         mp=mp+npart

c     end of event loop
 10   continue

      write(6,*)'no. of collisions = ',mc/dble(nevents), ' (per event)'
      write(6,*)'final particles   = ',mp/dble(nevents), ' (per event)'

      end

C####C##1#########2#########3#########4#########5#########6#########7##
      subroutine pcheck(s)
c
c  Author: L.A.Winckelmann
c
cinput  s : info-string  
c
c This routine prints out 4-momentum and charge of all particles
c and is used for debug purposes.

cccccCcc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      implicit none
      integer i,nch
      real*8 pcx,pcy,pcz,pc0
      character  s*7
      include 'coms.f'

      pcx=0d0
      pcy=0d0
      pcz=0d0
      pc0=0d0 
      nch=0
      do 108 i=1,npart
        pcx=pcx+px(i)+ffermpx(i)
        pcy=pcy+py(i)+ffermpy(i)
        pcz=pcz+pz(i)+ffermpz(i)
        pc0=pc0+p0(i)
        nch=nch+charge(i)
 108  continue 
      write(6,1008)s,' p_mu, charge',pc0,pcx,pcy,pcz,nch
      return
 1008 format(1X,A7,A12,4e14.6,i3)

      end

 

