ccsws
cc  some variable for plotting
      real*8 dtplot,tplot 
cc  size of plot area in horiz.direction zmin,zmax
cc  radius of plotting symbol baryons: cradius, mesons: crad2
      real zmin,zmax,xmin,xmax,rvalx,rvalz,cradius,crad2
cc  IPLOT flag for plotting 
      integer IPLOT,numc
cc  sttime contains time string
      character*32 sttime
cc  output device
	character*3 device
cc    isteps = max number of plotting events 
      integer ISTEPS
	common /plotit/dtplot,tplot,cradius,crad2,
     $               zmin,zmax,xmin,xmax,iplot,isteps,device
ccsws