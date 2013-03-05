c $Id: output.f,v 1.11 1997/08/26 12:31:17 weber Exp $
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      subroutine output(iunit)
c
c     Unit     : Collision Term
c     Author   : Steffen A. Bass, Luke A. Winckelmann
c     Date     : 01/09/95
c     Revision : 0.9 beta - uncompleted
c
c     This subroutine writes the event-header to file(iunit)
C
c                    
c     input : 
c             iunit  : output-unit
c
c     output: 
c
c
c important global variables:
c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      implicit none

      include 'comres.f'
      include 'coms.f'
      include 'options.f'
      include 'inputs.f'
      include 'newpart.f'
      include 'freezeout.f'
      include 'boxinc.f'

c
      integer iunit,i,ttime,iu,app,att,zpp,ztt
      integer iiunit,isunit, id, pdgid
      integer timestep,itotcoll,iinelcoll
      real*8 sigmatot,ptsigtot,stot,otime
      common /outco2/sigmatot


      character*4 reffram
      character*20 aa,ah,ai,ak
      character*36 ae,abt
      character*31 aee
      character*15 ab,aj,al,am
      character*13 ac,ag,pds,tds
      character*12 ad
      character*7 af
      character*9 ag2
      character*1 add
      character*171 apa14,apa15,apav
      character*2 apa,aop

c file15out
      integer ind,ind1,ind2,nin
      integer istr,ich,ii

      real*8 sqrts, sigpart, colldens, cdens,cdens_
      logical bdum,paulibl

c temporary storage for in-channel
      include 'outcom.f'
c      real*8 tsqrts,tstot,tsigpart
c      real*8 tr0(3),trx(3),try(3),trz(3),ttform(3),txtotfac(3),
c     @     tp0(3),tpx(3),tpy(3),tpz(3),tm(3),torigin(3)
c      integer tind(3),tityp(3),tiso3(3),tcoll(3),tstrange(3)
c      integer tlcoll(3),tcharge(3),torigin(3),tstrid(3)   


     
c functions
      integer fchg,strit
      character*1 echar
      
      integer iou(13:17)

      save ! important, do not erase 

      data iou/13,14,15,16,17/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              output formats
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c fileheader
 101  format(a20,3i7,a15,i2)
 301  format(a13,a13,i4,i4,a12,a13,i4,i4,a1)
c 305  format(a36,3f10.7)
 304  format(a36,3f6.2,a31,1f9.2)
 302  format (a7,i9,a13,i12,a9,a20,i4,a20,f5.1)
c 303  format(a20,i3,a15,e10.4,a15,e10.4,a15,e10.4)
 102  format(a2,15(i3,a2))
c 103  format(a2,12(e10.4,a2))
 306  format(a171)

 305  format(a36,3f11.7)
 303  format(a20,i3,a15,e11.4,a15,e11.4,a15,e11.4)
 103  format(a2,12(e11.4,a2))

corig
c standard particle information vector
c 201  format(9e16.8,i5,2i3,i6,i5,i4)
c special output for cto40 (restart of old event)
c 210  format(9e16.8,i5,2i3,i6,i5,i4,3e16.8)
c special output for mmaker
c 203  format(9e16.8,i5,2i3,i6,i5,2i4,2e16.8)
c same with index for file15
c 501  format(i5,9e16.8,i5,2i3,i6,i5,i3,i15)
c enhanced file16
c 503  format(9e15.7,i5,2i3,i6,i5,i4,2i4)
c same including freeze-out coordinates
c 213  format(9e16.8,i5,2i3,i6,i5,i4,8e16.8)

c standard particle information vector
 201  format(9D24.16,i5,2i3,i6,i5,i4)
c special output for cto40 (restart of old event)
 210  format(9D24.16,i5,2i3,i6,i5,i4,3D24.16)
c special output for mmaker
 203  format(9D24.16,i5,2i3,i6,i5,2i4,2D24.16)
c same with index for file15
 501  format(i5,9D24.16,i5,2i3,i6,i5,i3,i15)
c enhanced file16
 503  format(9e15.7,i5,2i3,i6,i5,i4,2i4)
c same including freeze-out coordinates
 213  format(9D24.16,i5,2i3,i6,i5,i4,8D24.16)

c collsision stats for file14
 202  format(8i8)
c same with EndOfEvent tag for file16
 602  format(a1,8i8)

c header-line for each collision in file15
 502  format(i1,i8,i4,i7,f8.3,4e12.4)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(iunit.eq.17)return
      if(bf13.and.(iunit.eq.13)) return
      if(bf14.and.(iunit.eq.14)) return
      if(bf15.and.(iunit.eq.15)) return
      if(bf16.and.(iunit.eq.16)) return

c     copy projectile/target info to local vars
      app=ap
      zpp=zp
      att=at
      ztt=zt

      if(iunit.eq.19) return
c
      aa='UQMD   version:     '
      ab='  output_file '
      abt='transformation betas (NN,lab,pro) '
      ac='projectile:  '
      ad='   target: '
      add=' '
      ae='impact_parameter_real/min/max(fm):  '
      aee='  total_cross_section(mbarn):  '
      af='event# '
      ag=' random seed:' 
      ah='equation_of_state: '
      ai=' total_time(fm/c): '
      aj='  E_lab(GeV/u):'
      ak=' Delta(t)_O(fm/c): '
      al='  sqrt(s)(GeV):'
      am='  p_lab(GeV/u):'
      apa='pa'
      aop='op'

      apa14='pvec: '//
     & 'r0              rx              ry              rz          '//
     & '    p0              px              py              pz      '//
     & '        m          ityp 2i3 chg lcl#  ncl or'    
      apa15='pvec:ind   '//
     & 'r0              rx              ry              rz          '//
     & '    p0              px              py              pz      '//
     & '        m          ityp 2i3 chg lcl#  ncl st'    
      if(iunit.eq.15) then
         apav=apa15
      else
         apav=apa14
      endif

      if(fixedseed) then
         ag2=' (fixed) '
      else
         ag2=' (auto)  '
      endif
      if(prspflg.eq.1) then
         pds='(ityp, char) '
         app=spityp(1)
         zpp=fchg(spiso3(1),app)
      else
         pds='(mass, char) '
      endif
      if(trspflg.eq.1) then
         tds='(ityp, char) '
         att=spityp(2)
         ztt=fchg(spiso3(2),att)
      else
         tds='(mass, char) '
      endif

c determine cross section of the projectile-target system
      sigmatot = ptsigtot()
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      otime=outsteps*dtimestep
      ttime=int(nsteps*dtimestep+0.01)

c... dont change the output of files you dont need without contacting
c	those who need the files !!!! 
      if(iunit.eq.15)then
       write(iou(15),502)0,event,Ap,At,bimp,ecm
     ,     ,sigmatot,ebeam,pbeam       
      else
      write(iou(iunit),101) aa,version, sigver, laires, ab,iunit
      write(iou(iunit),301) ac,pds, App, Zpp, ad,tds, Att, Ztt,add       
      write(iou(iunit),305) abt,betann,betatar,betapro
      write(iou(iunit),304) ae,bimp,bmin,bdist,aee,sigmatot
      write(iou(iunit),303) ah,eos,aj,ebeam,al,ecm,am,pbeam
      write(iou(iunit),302) af,event,ag,ranseed,ag2,ai,ttime,ak,otime
      write(iou(iunit),102) aop,(CTOption(i),CTOdc(i),i=1,15)
      write(iou(iunit),102) aop,(CTOption(i),CTOdc(i),i=16,30)
      write(iou(iunit),102) aop,(CTOption(i),CTOdc(i),i=31,45)
      write(iou(iunit),103) apa,(CTParam(i),CTPdc(i),i=1,12)
      write(iou(iunit),103) apa,(CTParam(i),CTPdc(i),i=13,24)
      write(iou(iunit),103) apa,(CTParam(i),CTPdc(i),i=25,36)
      write(iou(iunit),103) apa,(CTParam(i),CTPdc(i),i=37,48)
      write(iou(iunit),306) apav
      end if

c 
      return
c.....
      entry uounit(iiunit,isunit)
      iou(iiunit)=isunit
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file14out(timestep)
c
c     Unit     : Collision Term
c     Author   : Steffen A. Bass (new source)
c     Date     : 01/09/95
c     Revision : 0.1 beta - uncompleted
c
c     This subroutine writes the standard output-file (unit 14)
c                    
c     input : 
c             timestep  : timestep of output
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

c
      if(bf14)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      write(iou(14),*) npart,ttime
      write(iou(14),202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes

c now write particle-output

cernst write spectators
      if(CTOption(28).eq.2)then
         if(CTOption(41).eq.0) then
            do 141 i=1,nspec
               write(iou(14),201) r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @              pxs(i),pys(i),pzs(i),sfmass(i),
     @              sityp(i),siso3(i),scharge(i),
     @              -1,-1,0
 141        continue
         else
            do 142 i=1,nspec
               write(iou(14),210) r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @              pxs(i),pys(i),pzs(i),sfmass(i),
     @              sityp(i),siso3(i),scharge(i),
     @              -1,-1,0,1d34,0d0,1d0
 142        continue
         endif
      endif
chw special output for mmaker
      if (CTOption(25).eq.1) then
	 write (6,*) 'Timestep:',timestep
         do 10 i=1,npart
            write(iou(14),203) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),lstcoll(i),
     @           ncoll(i),mod(origin(i),100),strid(i),tform(i),
     &           xtotfac(i)
 10      continue
      else
         if(CTOption(41).eq.0) then
chw normal output
            do 13 i=1,npart
               write(iou(14),201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @              px(i)+ffermpx(i),py(i)+ffermpy(i),
     @              pz(i)+ffermpz(i),fmass(i),
     @              ityp(i),iso3(i),charge(i),
     @              lstcoll(i),ncoll(i),mod(origin(i),100)
 13         continue
         else
            do 31 i=1,npart
               write(iou(14),210) r0(i),rx(i),ry(i),rz(i),p0(i),
     @              px(i)+ffermpx(i),py(i)+ffermpy(i),
     @              pz(i)+ffermpz(i),fmass(i),
     @              ityp(i),iso3(i),charge(i),
     @              lstcoll(i),ncoll(i),mod(origin(i),100),
     @              dectime(i),tform(i),xtotfac(i)
 31         continue

         endif

      endif
c 
      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file13out(timestep)
c
c     Unit     : Collision Term
c     Author   : Steffen A. Bass (new source)
c     Date     : 01/09/95
c     Revision : 0.1 beta - uncompleted
c
c     This subroutine writes the standard output-file (unit 13),
c     including the freeze-out configuration of the particles
c                    
c     input : 
c             timestep  : timestep of output
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c

c
      if(bf13)return
      ttime=int(timestep*dtimestep+0.01)
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl
      write(iou(13),*) npart,ttime
      write(iou(13),202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes

c now write particle-output

cernst write spectators
        if(CTOption(28).eq.2)then
          do 191 i=1,nspec
            write(iou(13),213) r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @        pxs(i),pys(i),pzs(i),sfmass(i),
     @        sityp(i),siso3(i),scharge(i),
     @        -1,-1,0,r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @        pxs(i),pys(i),pzs(i)
 191      continue
        endif


      do 90 i=1,npart
         if(ncoll(i).eq.0) then
            write(iou(13),213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           r0(i),rx(i),ry(i),rz(i),p0(i),px(i)+ffermpx(i),
     @           py(i)+ffermpy(i),pz(i)+ffermpz(i)
         else
            write(iou(13),213) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           lstcoll(i),ncoll(i),mod(origin(i),100),
     @           frr0(i),frrx(i),frry(i),frrz(i),frp0(i),frpx(i),
     @           frpy(i),frpz(i)
         endif
 90   continue
c 
      return

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      entry file15out(ind1,ind2,sqrts,stot,sigpart)
c     Unit     : Collision Term
c     Author   : L. A. Winckelmann, S. A. Bass
c     Date     : 01/09/95
c     Revision : 0.1 beta - uncompleted
c     This subroutine writes information about the in-channel to file15
c     (the collision statistics file)
c     input : 
c             ind1    : index of particle 1
c             ind2    : index of particle 2 (=0 for decay of ind1)
c             sqrts   : sqrt(s) of collision
c             stot	: total cross section
c             sigpart	: partial cross section
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c we need some of these quantities even if bf15.eq.false
c      if(bf15.and.(CTOption(13).eq.0))return

c determine tag for scatter-input or decay-input
c and store entry channel in temporary observables
      bdum=paulibl(ind1,cdens)
      tsqrts=sqrts
      tstot=stot
      tsigpart=sigpart

      tind(1)=ind1
      tr0(1)=r0(ind1)
      trx(1)=rx(ind1)
      try(1)=ry(ind1)
      trz(1)=rz(ind1)
      tp0(1)=p0(ind1)
      tpx(1)=px(ind1)
      tpy(1)=py(ind1)
      tpz(1)=pz(ind1)
      tm(1)=fmass(ind1)
      tityp(1)=ityp(ind1)
      tiso3(1)=iso3(ind1)
      tstrange(1) = strit(tityp(1))    !  insert strangeness here
      tcoll(1) = ncoll(ind1)
      tlcoll(1)=lstcoll(ind1)
chw
      torigin(1)=origin(ind1)
c      ttform(1)=tform(ind1)
c      txtotfac(1)=xtotfac(ind1)
c      tstrid(1)=strid(ind1)
c mbbox
      if(ind2.le.0) then
         nin=1
      elseif(ind2.gt.0) then
         bdum=paulibl(ind1,cdens_)
         cdens=5d-1*(cdens+cdens_)
         nin=2
         tind(2)=ind2
         tr0(2)=r0(ind2)
         trx(2)=rx(ind2)
         try(2)=ry(ind2)
         trz(2)=rz(ind2)
         tp0(2)=p0(ind2)
         tpx(2)=px(ind2)
         tpy(2)=py(ind2)
         tpz(2)=pz(ind2)
         tm(2)=fmass(ind2)
         tityp(2)=ityp(ind2)
         tiso3(2)=iso3(ind2)
         tstrange(2)=strit(tityp(2))    ! insert strangeness here
         tcoll(2) = ncoll(ind2)
         tlcoll(2)=lstcoll(ind2)
chw
         torigin(2)=origin(ind2)
c         ttform(2)=tform(ind2)
c         txtotfac(2)=xtotfac(ind2)
c         tstrid(2)=strid(ind2)
      endif

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f15outch(colldens)

c      call dnscnt(nin,nexit,iline,acttime,colldens)

c      write(6,502) nin,nexit,iline,ctag,acttime,sqrts
c     ,     ,stot,sigpart,colldens

      if (bf15) return
c     This entry writes information about the collision to file15
c     one line for each particle:
c     format: x y z px py pz ityp iso3...

      write(iou(15),502) nin,nexit,iline,ctag,acttime,tsqrts
     ,     ,tstot,tsigpart,cdens
c      write(6,*)torigin(1)
      do 11 i=1,nin
         istr=strit(tityp(i))
         ich = fchg(tiso3(i),tityp(i))
         
         write(iou(15),501) tind(i),tr0(i),trx(i),try(i),trz(i),
     @                   tp0(i),tpx(i),tpy(i),tpz(i),tm(i),
     @                   tityp(i),tiso3(i),ich,tlcoll(i),
     @                   tcoll(i),istr,torigin(i)
 11   continue
      do 20 ii=1,nexit
         i=inew(ii)
         istr=strit(ityp(i))
c         ich = fchg(iso3(i),ityp(i))
         write(iou(15),501) i,r0(i),rx(i),ry(i),rz(i),
     @                   p0(i),px(i),py(i),pz(i),fmass(i),
     @                   ityp(i),iso3(i),charge(i),lstcoll(i),
     @                   ncoll(i),istr,origin(i) 
 20   continue

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry f16outch
c      write (6,*) 'bf16 cto',bf16,CTOption(13)
      if (bf16.or.(CTOption(13).eq.0)) return
c      write (6,*) 'f16outch'
      if (nin.eq.1) then
         tityp(2)=0
      endif
      
      do 22 ii=1,nexit
         i=inew(ii)
         write(iou(16),503) r0(i),rx(i),ry(i),rz(i),
     @        p0(i),px(i)+ffermpx(i),py(i)+ffermpy(i),
     @        pz(i)+ffermpz(i),fmass(i),
     @        ityp(i),iso3(i),charge(i),lstcoll(i),
     @        ncoll(i),mod(origin(i),100),tityp(1),tityp(2)
 22   continue

      return

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16out

      echar='E'
      itotcoll=ctag-dectag
      iinelcoll=itotcoll-NBlColl-NElColl

c   
      if(bf16) return

c now write particle-output
      if (CToption(13).eq.0) then
c ernst write spectators
         if (CTOption(28).eq.2)then
            do 151 i=1,nspec
               write(iou(16),201)r0s(i),rxs(i),rys(i),rzs(i),p0s(i),
     @              pxs(i),pys(i),pzs(i),sfmass(i),
     @              sityp(i),siso3(i),scharge(i),
     @              -1,-1,0
 151        continue
         endif

         do 12 i=1,npart
            write(iou(16),201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           dectag+lstcoll(i),ncoll(i),mod(origin(i),100)

 12      continue
      else
         do 14 i=1,npart
            write(iou(16),503) r0(i),rx(i),ry(i),rz(i),p0(i),
     @           px(i)+ffermpx(i),py(i)+ffermpy(i),
     @           pz(i)+ffermpz(i),fmass(i),
     @           ityp(i),iso3(i),charge(i),
     @           dectag+lstcoll(i),ncoll(i),mod(origin(i),100),-99,-99
 14      continue
      endif
c
c write collision counters etc.	
       write(iou(16),602) echar,itotcoll,NElColl,iinelcoll,NBlColl,
     @     dectag,NHardRes,NSoftRes,NDecRes
      return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16entry(ind)
c
c     This entry stores one decay for later output (must be done, in case
c     of pauli-blocked decay)
c
      tr0(3)=r0(ind)
      trx(3)=rx(ind)
      try(3)=ry(ind)
      trz(3)=rz(ind)
      tp0(3)=p0(ind)
      tpx(3)=px(ind)
      tpy(3)=py(ind)
      tpz(3)=pz(ind)
      tm(3)=fmass(ind)
      tityp(3)=ityp(ind)
      tind(3)=ind
      tiso3(3)=iso3(ind)
      tcharge(3)=charge(ind)

c     lstcoll is negative to identify decayed particles 

      tlcoll(3)=-(1*lstcoll(ind))
      tcoll(3)=ncoll(ind)
      torigin(3)=origin(ind)

      return
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry file16write

c     This entry writes the decay to file
      i=3

      if(bf16)return
      if (CTOption(13).eq.0) then
       
      write(iou(16),201) tr0(i),trx(i),try(i),trz(i),tp0(i),tpx(i),
     @        tpy(i),tpz(i),tm(i),tityp(i),tiso3(i),tcharge(i),
     @        tlcoll(i),tcoll(i),mod(torigin(i),100)
      else
      write(iou(16),503) tr0(i),trx(i),try(i),trz(i),tp0(i),tpx(i),
     @        tpy(i),tpz(i),tm(i),tityp(i),tiso3(i),tcharge(i),
     @        tlcoll(i),tcoll(i),mod(torigin(i),100),-98,-98
      endif

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_header

      if (bf19) return
      
      write (19,901) 'OSC1997A    '
      write (19,901) 'final_id_p_x'

 901  format (a12)

      if (CTOption(27).eq.0) then
         reffram='eqsp'
      elseif (CTOption(27).eq.1) then
         reffram='tar'
      elseif (CTOption(27).eq.2) then
         reffram='pro'
      else
         call error ('osc_header','Unknown Ref-Frame',
     .        dble(CTOption(27)),2)
         reffram='----'
      endif

      write (19,902) 'UrQMD', '1.0', app, zpp, att, ztt,
     .     reffram, ebeam, 1

 902  format (2(a8,2x),'(',i3,',',i6,')+(',i3,',',i6,')',2x,a4,2x,
     &     e10.4,2x,i8)

      return

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry osc_event

      if (bf19) return

      write (19,903) event, npart, bimp, 0D0

 903  format (i10,2x,i10,2x,f8.3,2x,f8.3)

c particles

      do 99 i=1,npart
         id = pdgid(ityp(i), iso3(i))
         write(19,904) i, id, 
     .        px(i), py(i), pz(i), p0(i), fmass(i),     
     .        rx(i), ry(i), rz(i), r0(i)
 99   continue

 904  format (i10,2x,i10,2x,9(d24.16,2x))

      return


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry getoldevent

c     read event header
      read(10,*,end=666) 
     @ aa

cklappt so nicht,version, sigver, laires, ab,iunit
      read(10,301) ac,pds, App, Zpp, ad,tds, Att, Ztt,add       
      read(10,305) abt,betann,betatar,betapro
      read(10,304) ae,bimp,bmin,bdist,aee,sigmatot
      read(10,303) ah,eos,aj,ebeam,al,ecm,am,pbeam
      read(10,302) af,event,ag,ranseed,ag2,ai,ttime,ak,otime
      read(10,102) aop,(CTOption(i),CTOdc(i),i=1,15)
      read(10,102) aop,(CTOption(i),CTOdc(i),i=16,30)
      read(10,102) aop,(CTOption(i),CTOdc(i),i=31,45)
      read(10,103) apa,(CTParam(i),CTPdc(i),i=1,12)
      read(10,103) apa,(CTParam(i),CTPdc(i),i=13,24)
      read(10,103) apa,(CTParam(i),CTPdc(i),i=25,36)
      read(10,103) apa,(CTParam(i),CTPdc(i),i=37,48)
      read(10,306) apav
c reset option 40
      CTOption(40)=1

c read event body
      read(10,*) npart,ttime
      read(10,202) itotcoll,NElColl,iinelcoll,NBlColl,dectag,
     @     NHardRes,NSoftRes,NDecRes
c      timestep=dble(ttime)/dtimestep
      ctag=itotcoll+dectag
c now read particle-output
      nbar=0
      do 39 i=1,npart
         read(10,210) r0(i),rx(i),ry(i),rz(i),p0(i),
     @        px(i),py(i),pz(i),fmass(i),
     @        ityp(i),iso3(i),charge(i),
     @        lstcoll(i),ncoll(i),origin(i),
     @        dectime(i),tform(i),xtotfac(i)
      if(abs(ityp(i)).le.maxbar)nbar=nbar+1
 39   continue
      nmes=npart-nbar
      acttime=r0(1) 
c     read options-file
      call getparams
      return
c stop in case of EoF
 666  stop


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      entry specout(ind,iu)
      i=ind
      if (CTOption(28).lt.0) return 
      if (iu.eq.16.and.bf16) return
      if (iu.eq.14.and.bf14) return
      write(iu,201) r0(i),rx(i),ry(i),rz(i),p0(i),
     @     px(i)+ffermpx(i),py(i)+ffermpy(i),
     @     pz(i)+ffermpz(i),fmass(i),
     @     ityp(i),iso3(i),charge(i),
     @     dectag+lstcoll(i),ncoll(i),mod(origin(i),100)
 
      return

      end



       
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine spectrans(tstep)
c  (when cto 28 is set to 2 this subroutine is called
c  to propagate the spectators along straight lines)
c  
c  Author: C. Ernst 1/96 (modified source)
c  Input tstep : timestep
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 dtime,energ,tstep
      integer j
      include 'coms.f'

      dtime=tstep

      do 1 j=1,nspec
         energ = p0s(j)
         r0s(j) = r0s(j) + dtime
         rxs(j) = rxs(j) + pxs(j)/energ*dtime
         rys(j) = rys(j) + pys(j)/energ*dtime
         rzs(j) = rzs(j) + pzs(j)/energ*dtime
1     continue

      return
      end
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
	real*8 function bardens (k)
c	Unit	  :	Output
c	Author  :	Henning Weber, Christoph Ernst
c
cinput k  : index of particle 
c
c Calculates the baryon density at location of particle k.
c ref. H. Weber,  diploma thesis
c ATTENTION: the loop extends only to nbar. You should put nmax instead
c   if you are not sure wether the particles are properly ordered.
c ATTENTION2: a cut-off for very distant baryons is included to save 
c   runtime, however this might underestimate the density.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	implicit none
	include 'coms.f'
	
	integer k,i,collk
	real*8 rhob, qij(4), pij(4),q2,p2,qp
	real*8 v1,v2,beta2,gamma2,dr,re,sr,cutoff
ccut	real*8 rap,rapk
      parameter (re = 8.66/4.0)
      parameter (cutoff=20d0)

cc	rap(k)=.5d0*log((p0(k)+pz(k))/(p0(k)-pz(k))
ccut	rapk=rap(k)
	collk=lstcoll(k)
	
	rhob = 0d0
      do 37 i=1,nbar
          if (abs(ityp(i)).lt.100.and.i.ne.collk
ccut     &     .and.dabs(rap(i)-rapk).lt.1.5D0
     &     )then    
c calculate space q_{ij}=q_i-q_j
            qij(4)=r0(i)-r0(k)
            qij(1)=rx(i)-rx(k)
            qij(2)=ry(i)-ry(k)
            qij(3)=rz(i)-rz(k)
            q2=qij(4)*qij(4)-qij(1)*qij(1)-qij(2)*qij(2)
     6           -qij(3)*qij(3)
     		sr=sqrt(-q2)
ce!! cutoff for very distant particles to save runtime
ce            if(sr.gt.cutoff)goto 37

c calculate momentum p_{ij}=p_i+p_j
            pij(4)=p0(i)+p0(k)
            pij(1)=px(i)+ffermpx(i)+px(k)+ffermpx(k)
            pij(2)=py(i)+ffermpy(i)+py(k)+ffermpy(k)
            pij(3)=pz(i)+ffermpz(i)+pz(k)+ffermpz(k)  
            
            p2=pij(4)*pij(4)-pij(1)*pij(1)-pij(2)*pij(2)
     6           -pij(3)*pij(3)
            qp=qij(4)*pij(4)-qij(1)*pij(1)-qij(2)*pij(2)
     .          -qij(3)*pij(3)
            dr=qp*qp/p2 - q2
ce!! cutoff for large dr 
ce!! for gamma=1 dr=20 corresponds to rhob \approx 1e-5 \approx 0
		if(dr.gt.cutoff)goto 37
c calculate beta^2 and gamma^2 in r-direction
	      if (i.eq.k.or.dr.lt.1E-10) then
	        beta2=0.0
	      else
	        v1=(qij(1)*(px(i)+ffermpx(i))+qij(2)*(py(i)+ffermpy(i))
     .	             +qij(3)*(pz(i)+ffermpz(i)))/p0(i)/sr
	        v2=(qij(1)*(px(k)+ffermpx(k))+qij(2)*(py(k)+ffermpy(k))
     .	             +qij(3)*(pz(k)+ffermpz(k)))/p0(k)/sr
              beta2=((v1-v2)/(1-v1*v2))**2
            endif
            gamma2=1.0/(1.0-beta2)
            if (beta2.gt.1.or.beta2.lt.0) then
               write (6,*) '(W) bardens: Wrong Beta!',beta2
            endif
            rhob=rhob+sign(1,ityp(k))*
     .         (pi*re)**(-1.5)*exp(-(dr/re*gamma2))/rho0
          endif
 37   continue
	bardens = rhob
	return
	end

ce	real*8 function rap(i)
ce	implicit none
ce	integer i
ce	include 'coms.f'
ce	rap=.5d0*dlog( (p0(i)+pz(i))/ (p0(i)-pz(i)))
ce	return
ce	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	real*8 function ptsigtot()
c
c     Unit     : Collision Term
c     Author   : Steffen A. Bass, Luke A. Winckelmann, Christoph Ernst
c     Date     : 07/07/97
c     Revision : 0.9 beta - uncompleted
c
c     This function caculates the total cross section of the reaction.
c	(Projectile - target total cross section)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      implicit none

      include 'comres.f'
      include 'coms.f'
      include 'options.f'

	integer indmn,indmx,itypmn,iso3mn,itypmx,iso3mx
	integer isigline,iline,collclass
	real*8 stot,sigel
	real*8 sigtot

c determine total cross section for reaction:
      if(Ap+At.gt.2) then
         stot=10.d0*pi*(bdist**2-bmin**2)
      elseif(Ap+At.eq.2) then
         stot=sigtot(1,2,ecm)
cccccccc for CTOption(7)=1 no elastic cross section:
         if(CTOption(7).eq.1) then
c first sort the two itypes for call to collclass and anndec
            if(abs(ityp(1)).lt.abs(ityp(2))) then
               indmn=1
               indmx=2
            else
               indmn=2
               indmx=1
            endif

            itypmn=ityp(indmn)
            iso3mn=iso3(indmn)
            itypmx=ityp(indmx)
            iso3mx=iso3(indmx)
            isigline=collclass(itypmx,iso3mx,itypmn,iso3mn)
c     the elastic cross section is always the first entry (#3)
            iline=SigmaLn(3,1,isigline)
c!!!!DANGER: does not work for unstable particles (-> detailed balance)
            call crossx(iline,ecm,ityp(1),iso3(1),
     &              fmass(1),ityp(2),iso3(2),fmass(2),sigel)
c
            if(stot-sigel.gt.0) then
               stot=stot-sigel
            else
               stot=sigel
            endif
         endif
      else
         stot=0.d0
      endif
c
	ptsigtot=stot
	return
	end

