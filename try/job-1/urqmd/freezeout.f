c $Id: freezeout.f,v 1.3 1998/06/15 13:35:20 weber Exp $
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c freezeout common block for uQMD
c     Unit     : output/collision term
c     Author   : Steffen A. Bass
c     Date     : 01/31/96
c     Revision : 1.0 
c
c

      real*8 frr0(nmax), frrx(nmax), frry(nmax), frrz(nmax),
     +     frp0(nmax), frpx(nmax), frpy(nmax), frpz(nmax)

c 8*nmax*nmax + 40*nmax real*8        
      
      common /frcoor/ frr0, frrx, frry, frrz, frp0, frpx, frpy, frpz 
