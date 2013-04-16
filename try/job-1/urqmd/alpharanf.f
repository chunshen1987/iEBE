c $Id: alpharanf.f,v 1.1.1.1 1996/12/02 15:56:03 weber Exp $
C*****************   R A N F   *******************************
c sab 05.03.1992
      real*8 FUNCTION RANF(ix)
      integer ix
      real*8 rand
      external rand
       ranf=dble(rand())
c       write(7,*) ranf
      RETURN
      end
c
c
      subroutine SSEED(ranseed)
      integer ranseed,oldseed,time
      external time,srand
      save

      itt=0


      if(ranseed.gt.0) then
         call srand(abs(ranseed))
c         WRITE(6,*)'FIXED SEED = ',ranseed
         return
      endif
      ranseed=time()
c      write(7,*)'seed',ranseed
c      do 123 j=1,9
c         itt=itt*ts(10-j)+timerr(j) 
c 123  continue
c      itt=itt+timerr(6)*timerr(5)*timerr(1)**2
c      ranseed=abs(2*itt+1)
      if(abs(oldseed).eq.ranseed) then
         ranseed=-1*abs(ranseed)
         return
      else
         oldseed=ranseed
      endif
c      write(6,*)'AUTO-SEED = ',ranseed
      call srand(ranseed)
      RETURN
      END
