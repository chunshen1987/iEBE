cc  plot the particles directly onto the screen while calculating
	subroutine plotlive
      implicit none
      include 'coms.f'
	include 'complot.f'
	integer j
        real cradx
        integer ilam,ikaon,idelta
c  check whether plotting flag is set
      if(iplot.eq.0) return	
      ilam=0
      ikaon=0
      idelta=0
cc  check whether at least dtplot time passed after last plotting
	if(acttime.lt.(tplot+dtplot)) return
cc  put all plotting commands in a buffer before executing
cc  remove (most) flickering
	CALL PGBBUF
cc  erase last plot 
	call pgeras
cc  convert time to string
      CALL PGNUMB(int(acttime*10000),-4,1,STTIME,NUMC)  
cc  plot time string
      CALL PGSCI(1)
      call PGTEXT(zmin,XMIN+1.0,'TIME: '//STTIME)
cc  set color from palette
            CALL PGSCI(4)
cc  plot projectile nucleons
	do 700 j=1,min(npart,abs(ap))
		rvalz = rz(j)
		rvalx = rx(j)
cc  use a filled circle for plotting
               cradx=amax1(.01e0,sngl(cradius*(1.+ry(j)/20.)))
		call PGCIRC(rvalz,rvalx,cradx)
700	continue
cc  same with target
	CALL PGSCI(3)
	do 702 j=abs(ap)+1,min(npart,abs(at)+abs(ap))
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(cradius*(1.+ry(j)/20.)))
		call PGCIRC(rvalz,rvalx,cradx)
702	continue
cc  same with produced particles
	CALL PGSCI(8)
	do 703 j=abs(ap)+abs(at)+1,npart
cc  select mesons only
		if(abs(ityp(j)).le.99) goto 703
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(crad2*(1.+ry(j)/20.)))
cc  use smaller radius crad2
		call PGCIRC(rvalz,rvalx,cradx)
703	continue
	CALL PGSCI(1)
	do 704 j=abs(ap)+abs(at)+1,npart
cc  plot produced baryons only 
		if((ityp(j).lt.0).or.(ityp(j).gt.99)) goto 704
		if(abs(ityp(j)).gt.99) goto 704
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(cradius*(1.+ry(j)/20.)))
		call PGCIRC(rvalz,rvalx,cradx);
704	continue
	CALL PGSCI(2)
	do 705 j=abs(ap)+abs(at)+1,npart
cc  plot produced antibaryons only
		if((ityp(j).ge.0).or.(ityp(j).lt.-99)) goto 705
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(cradius*(1.+ry(j)/20.)))
		call PGCIRC(rvalz,rvalx,cradx)
705	continue
cc
cc  now specific particles
cc
	CALL PGSCI(5)
	do 706 j=1,npart
cc lambda 
		if((ityp(j).lt.27).or.(ityp(j).gt.38)) goto 706 
                ilam=ilam+1
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(cradius*(1.+ry(j)/20.)))
		call PGCIRC(rvalz,rvalx,cradx)
706	continue
	CALL PGSCI(6)
	do 707 j=1,npart
cc kaon 
		if(abs(ityp(j)).ne.106) goto 707
                ikaon=ikaon+1
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(crad2*(1.+ry(j)/20.)))
		call PGCIRC(rvalz,rvalx,cradx)
707	continue
      call pgsci(1) 
	CALL PGSCI(7)
	do 708 j=1,npart
cc delta 
		if((ityp(j).lt.17).or.(ityp(j).gt.26)) goto 708
                idelta=idelta+1
		rvalz = rz(j)
		rvalx = rx(j)
                cradx=amax1(.01,sngl(cradius*(1.+ry(j)/20.)))*1.3
		call PGCIRC(rvalz,rvalx,cradx)
708	continue
cc  convert number of lambdas to string
      CALL PGNUMB(ilam,0,1,STTIME,NUMC)  
cc  plot  string
      CALL PGSCI(5)
      call PGTEXT(zmin+15.0,XMIN+1.0,'Lambdas: '//STTIME)
cc  convert number of kaons to string
      CALL PGNUMB(ikaon,0,1,STTIME,NUMC)  
cc  plot  string
      CALL PGSCI(6)
      call PGTEXT(zmin+30.,XMIN+1.0,'Kaons: '//STTIME)
cc  convert number of deltas to string
      CALL PGNUMB(idelta,0,1,STTIME,NUMC)  
cc  plot  string
      CALL PGSCI(7)
      call PGTEXT(zmin+40.,XMIN+1.0,'Deltas: '//STTIME)
cc  close buffer and putt everything on screen
      call PGEBUF
cc  remember time of last plotting
cc  close buffer and putt everything on screen
      call PGEBUF
cc  remember time of last plotting
	tplot=acttime
	return
	end


cc  initialize plotting window
cc
	subroutine plotinit
	implicit none
	include 'coms.f'
	include 'complot.f'
	include 'options.f'
	integer pgbeg
cc  check if plotting should be initialized
cc	if(iplot.eq.0) return
	ISTEPS=1000
      cradius=0.5
      crad2=0.4
      IPLOT=1
cc  size of output area
      ZMIN=-28.0
      ZMAX=28.0
      XMIN=-24.0
      XMAX=24.0
cc  output device (/ps, /xw (x windows), /wv (win95)
	device='/xw'
cc  open device
      if(PGBEG(0,device,1,1).NE.1) then
		iplot=0
		return
	end if 
cc  set output mode
      CALL PGVSTD
cc  set size of output area (inches, horiz. to vert. ratio)
      CALL PGPAP(8.0,1.0)
cc  change horizontal range depending on rest frame of calculation
      if(CTOption(27).eq.1) ZMIN=-ZMAX/4.0
      if(CTOption(27).eq.2) ZMAX=-ZMIN/4.0
cc  set x/y size of output area (in simulation coordinates)
      CALL PGENV(ZMIN,ZMAX,XMIN,XMAX,1,-2)
	return
	end

