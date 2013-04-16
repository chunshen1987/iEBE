c $Id: colltab.f,v 1.3 1997/06/16 20:52:37 bleicher Exp $
c
c     Autor : Markus Hofmann, Steffen A. Bass
c     Date  : 09/16/94
c
cdes  This file contains the  uqmd collision tables
c
      integer ncollmax
      parameter (ncollmax = 800000) ! maximum number of entries in collision table
      integer nct,actcol,nsav,apt
      real*8 cttime(0:ncollmax),ctsqrts(ncollmax),ctsigtot(ncollmax)
      real*8 ctcolfluc(ncollmax)
      logical ctvalid(ncollmax)
      real*8 tmin
      integer cti1(ncollmax),cti2(ncollmax),ctsav(ncollmax)
c      integer updi1(ncollmax),updi2(ncollmax)
c
c     cttime  : collision time
c     ctsqrts : $sqrt{s}$ of collision
c     ctsigtot: total cross section in mbarn
c     tmin    : paramteter for {\tt collupd}
c     cti1    : index of particle 1
c     cti2    : index of particle 2
c     nct     : number of collisions in the table
c     actcol  : current collision
c     ctvalid : tag whether collision is {\em true} or {\em false}
c     ctsav   : list of particles which lost their collision partner
c     nsav    : number of entries in {\tt ctsav}
c     apt     : mass of first particle/composite in the part. arrays 
      common /colltab/cttime,ctsqrts,ctsigtot,tmin,cti1,cti2,nct,actcol,
     &     ctvalid,ctsav,nsav,apt,ctcolfluc
