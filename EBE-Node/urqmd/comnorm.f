c $Id: comnorm.f,v 1.4 1997/06/26 11:36:38 konopka Exp $
c JK the quadratically weighted x_i's were not implemented correctly and 
c    e.g. lead to false normalization factors of the Breit-Wigners. 
c    This is a preliminary and quick bug-fix to maintain a physicswise 
c    senseful UrQMD. I hope that I found all places, where this kind of 
c    bug appeared. A more detailed update will follow up.       
      integer n
      parameter (n = 400)
      real*8 x_norm(0:3,1:n),y_norm(0:3,1:n)
      real*8 y2a(0:3,1:n),y2b(0:3,1:n),dx
      common /normsplin/ x_norm,y_norm,y2a,y2b,dx

