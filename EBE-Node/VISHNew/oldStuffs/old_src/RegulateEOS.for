!***************************************************************************!
!     This is a program for reading in the equation of state from           !
!     the original table and then interpolate and convert it into           !
!     a standard form for VISH2+1 can read in.                              !
!     The equation of state is got from fitting the lattice QCD             !
!     data and smooth cross over to the hadronic phase.                     !
!                                                                           !
!     Author: Chun Shen                                                     !
!     date: Oct. 29, 2010                                                   !
!***************************************************************************!

      Program Main

      Implicit none

!=======list of parameters===================================================
      Integer, Parameter :: EOSdatasize=501 !Original EOS table data size
      Integer, Parameter :: datasize=155500   !converted EOS table data size
!=======list of parameters end===============================================

!=======declear variables for reading equation of state======================
      Integer :: I, J, K, L, M, N   !loop variables
      Integer :: Itemp              !temperary variables
      Integer :: IEOS
      double precision :: Pe01, Pde1  !energy density data structure parametersfor data part 1
      double precision :: Pe02, Pde2  !energy density data structure parameters for data part 2
      double precision :: Pe03, Pde3  !energy density data structure parameters for data part 3
      double precision :: Pe04, Pde4  !energy density data structure parameters for data part 4
      double precision :: Se01, Sde1  !energy density data structure parametersfor data part 1
      double precision :: Se02, Sde2  !energy density data structure parameters for data part 2
      double precision :: Se03, Sde3  !energy density data structure parameters for data part 3
      double precision :: Se04, Sde4  !energy density data structure parameters for data part 4
      double precision :: Te01, Tde1  !energy density data structure parametersfor data part 1
      double precision :: Te02, Tde2  !energy density data structure parameters for data part 2
      double precision :: Te03, Tde3  !energy density data structure parameters for data part 3
      double precision :: Te04, Tde4  !energy density data structure parameters for data part 4
      double precision :: Peend, Seend, Teend   !the maximum energy density in the table
      double precision :: PEOSdata1(EOSdatasize),PEOSdata2(EOSdatasize),
     &                  PEOSdata3(EOSdatasize),PEOSdata4(EOSdatasize),
     &                  SEOSdata1(EOSdatasize),SEOSdata2(EOSdatasize),
     &                  SEOSdata3(EOSdatasize),SEOSdata4(EOSdatasize),
     &                  TEOSdata1(EOSdatasize),TEOSdata2(EOSdatasize),
     &                  TEOSdata3(EOSdatasize),TEOSdata4(EOSdatasize) !Pressure, Entropy, Temperature
      double precision :: ee, de, pp, ss, TT, ee0
      double precision :: PEOSL7, TEOSL7, SEOSL7
!=======declear variables for reading equation of state end==================

!=======common blocks========================================================
      common /EOSdataP/PEOSdata1,PEOSdata2,PEOSdata3,PEOSdata4
      common /EOSdataS/SEOSdata1,SEOSdata2,SEOSdata3,SEOSdata4
      common /EOSdataT/TEOSdata1,TEOSdata2,TEOSdata3,TEOSdata4

      common /EOSdatastructureP/ Pe01, Pde1, Pe02, Pde2, Pe03, Pde3,
     &                           Pe04, Pde4, Peend
      common /EOSdatastructureS/ Se01, Sde1, Se02, Sde2, Se03, Sde3,
     &                           Se04, Sde4, Seend
      common /EOSdatastructureT/ Te01, Tde1, Te02, Tde2, Te03, Tde3,
     &                           Te04, Tde4, Teend

      common /EOSSEL/ IEOS
!=======common blocks end====================================================

      call InputEOS

      open(2,FILE='EOS/s95p-PCE/EOS_Original.dat',FORM='Formatted', 
     &     Status='new')

      do I = 1, 501
            ee = Pe01 + Pde1*(I-1)
            write(2,1001)ee, PEOSdata1(I), SEOSdata1(I), TEOSdata1(I)
      enddo
      do I = 1, 501
            ee = Pe02 + Pde2*(I-1)
            write(2,1001)ee, PEOSdata2(I), SEOSdata2(I), TEOSdata2(I)
      enddo
      do I = 1, 501
            ee = Pe03 + Pde3*(I-1)
            write(2,1001)ee, PEOSdata3(I), SEOSdata3(I), TEOSdata3(I)
      enddo
      do I = 1, 501
            ee = Pe04 + Pde4*(I-1)
            write(2,1001)ee, PEOSdata4(I), SEOSdata4(I), TEOSdata4(I)
      enddo
      close(2)
1001  Format(4E15.6)

      open(1,FILE='EOS/s95p-PCE/EOS_converted.dat',FORM='Formatted',
     &     STATUS='new')

      ee0 = Pe01
      de  = Pde1
      do I = 1, datasize
            ee = ee0 + (I-1)*de
            pp = PEOSL7(ee)
            ss = SEOSL7(ee)
            TT = TEOSL7(ee)
            write(1,1000) ee, pp, ss, TT
      enddo
1000  Format(4E15.6)


      close(1)
      end

      Subroutine InputEOS

      Implicit none

!=======list of parameters===================================================
      Integer, Parameter :: EOSdatasize=501
!=======list of parameters end===============================================

!=======declear variables for reading equation of state======================
      Integer :: I, J, K, L, M, N   !loop variables
      double precision :: Itemp              !temperary variables
      Integer :: IEOS
      double precision :: nb0, e0   !lowest baryon density, lowest energy density
      double precision :: dnb
      Integer :: nnb          !grid spacing of baryon density, number of columns
      double precision :: de
      Integer :: ne, netot    !spacing in energy density, number of rows
      double precision :: Pe01, Pde1  !energy density data structure parametersfor data part 1
      double precision :: Pe02, Pde2  !energy density data structure parameters for data part 2
      double precision :: Pe03, Pde3  !energy density data structure parameters for data part 3
      double precision :: Pe04, Pde4  !energy density data structure parameters for data part 4
      double precision :: Se01, Sde1  !energy density data structure parametersfor data part 1
      double precision :: Se02, Sde2  !energy density data structure parameters for data part 2
      double precision :: Se03, Sde3  !energy density data structure parameters for data part 3
      double precision :: Se04, Sde4  !energy density data structure parameters for data part 4
      double precision :: Te01, Tde1  !energy density data structure parametersfor data part 1
      double precision :: Te02, Tde2  !energy density data structure parameters for data part 2
      double precision :: Te03, Tde3  !energy density data structure parameters for data part 3
      double precision :: Te04, Tde4  !energy density data structure parameters for data part 4
      double precision :: Peend, Seend, Teend   !the maximum energy density in the table
      double precision :: PEOSdata1(EOSdatasize),PEOSdata2(EOSdatasize),
     &                  PEOSdata3(EOSdatasize),PEOSdata4(EOSdatasize),
     &                  SEOSdata1(EOSdatasize),SEOSdata2(EOSdatasize),
     &                  SEOSdata3(EOSdatasize),SEOSdata4(EOSdatasize),
     &                  TEOSdata1(EOSdatasize),TEOSdata2(EOSdatasize),
     &                  TEOSdata3(EOSdatasize),TEOSdata4(EOSdatasize) !Pressure, Entropy, Temperature
!
!=======declear variables for reading equation of state end==================

!=======common blocks========================================================
      common /EOSdataP/PEOSdata1,PEOSdata2,PEOSdata3,PEOSdata4
      common /EOSdataS/SEOSdata1,SEOSdata2,SEOSdata3,SEOSdata4
      common /EOSdataT/TEOSdata1,TEOSdata2,TEOSdata3,TEOSdata4

      common /EOSdatastructureP/ Pe01, Pde1, Pe02, Pde2, Pe03, Pde3,
     &                           Pe04, Pde4, Peend
      common /EOSdatastructureS/ Se01, Sde1, Se02, Sde2, Se03, Sde3,
     &                           Se04, Sde4, Seend
      common /EOSdatastructureT/ Te01, Tde1, Te02, Tde2, Te03, Tde3,
     &                           Te04, Tde4, Teend

      common /EOSSEL/ IEOS
!=======common blocks end====================================================

!     Pressure(e)
C      PEOSdata1(1) = 0.0d0          !add a point P(0) = 0.0d0
      open(5,FILE='EOS/s95p-PCE/p1.dat',STATUS='OLD')
      read(5,*) nb0, Pe01
      read(5,*) dnb, nnb, Pde1, ne
      do I=1,ne
        read(5,*) PEOSdata1(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/p2.dat',STATUS='OLD')
      read(5,*) nb0, Pe02
      read(5,*) dnb, nnb, Pde2, ne
      do I=1,ne
        read(5,*) PEOSdata2(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/p3.dat',STATUS='OLD')
      read(5,*) nb0, Pe03
      read(5,*) dnb, nnb, Pde3, ne
      do I=1,ne
        read(5,*) PEOSdata3(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/p4.dat',STATUS='OLD')
      read(5,*) nb0, Pe04
      read(5,*) dnb, nnb, Pde4, ne
      Peend = Pe04+Pde4*(ne-1)            !the maximum energy density in the table
      do I=1,ne
        read(5,*) PEOSdata4(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

!     Entropy(e)
C      SEOSdata1(1) = 0.0d0          !add a point S(0) = 0.0d0
      open(5,FILE='EOS/s95p-PCE/s1.dat',STATUS='OLD')
      read(5,*) nb0, Se01
      read(5,*) dnb, nnb, Sde1, ne
      do I=1,ne
        read(5,*) SEOSdata1(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/s2.dat',STATUS='OLD')
      read(5,*) nb0, Se02
      read(5,*) dnb, nnb, Sde2, ne
      do I=1,ne
        read(5,*) SEOSdata2(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/s3.dat',STATUS='OLD')
      read(5,*) nb0, Se03
      read(5,*) dnb, nnb, Sde3, ne
      do I=1,ne
        read(5,*) SEOSdata3(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/s4.dat',STATUS='OLD')
      read(5,*) nb0, Se04
      read(5,*) dnb, nnb, Sde4, ne
      Seend = Se04+Sde4*(ne-1)            !the maximum energy density in the table
      do I=1,ne
        read(5,*) SEOSdata4(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

!     Temperature(e)
C      TEOSdata1(1) = 0.0d0          !add a point T(0) = 0.0d0
      open(5,FILE='EOS/s95p-PCE/t1.dat',STATUS='OLD')
      read(5,*) nb0, Te01
      read(5,*) dnb, nnb, Tde1, ne
      do I=1,ne
        read(5,*) TEOSdata1(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/t2.dat',STATUS='OLD')
      read(5,*) nb0, Te02
      read(5,*) dnb, nnb, Tde2, ne
      do I=1,ne
        read(5,*) TEOSdata2(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/t3.dat',STATUS='OLD')
      read(5,*) nb0, Te03
      read(5,*) dnb, nnb, Tde3, ne
      do I=1,ne
        read(5,*) TEOSdata3(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)

      open(5,FILE='EOS/s95p-PCE/t4.dat',STATUS='OLD')
      read(5,*) nb0, Te04
      read(5,*) dnb, nnb, Tde4, ne
      Teend = Te04+Tde4*(ne-1)            !the maximum energy density in the table
      do I=1,ne
        read(5,*) TEOSdata4(ne-I+1),Itemp  !read table and turn the order
      enddo 
      close(5)


      return
      end

CSHEN=================================================================
C====EOS from table===================================================
      Double Precision Function PEOSL7(ee)  ! for lattice P(e)
      Implicit none
       
!=======list of parameters===================================================
      Integer, Parameter :: EOSdatasize=501
!=======list of parameters end===============================================
      double precision :: ee
      double precision :: Pe01, Pde1  !energy density data structure parametersfor data part 1
      double precision :: Pe02, Pde2  !energy density data structure parameters for data part 2
      double precision :: Pe03, Pde3  !energy density data structure parameters for data part 3
      double precision :: Pe04, Pde4  !energy density data structure parameters for data part 4
      double precision :: Peend   !the maximum energy density in the table
      double precision :: PEOSdata1(EOSdatasize),PEOSdata2(EOSdatasize),
     &                    PEOSdata3(EOSdatasize),PEOSdata4(EOSdatasize) !Pressure

      common /EOSdataP/PEOSdata1,PEOSdata2,PEOSdata3,PEOSdata4

      common /EOSdatastructureP/ Pe01, Pde1, Pe02, Pde2, Pe03, Pde3,
     &                          Pe04, Pde4, Peend

      if(ee.lt.Pe01) then
C       PEOSL7 = ee*PEOSdata1(1)/Pe01
       PEOSL7 = (-3*PEOSdata1(1)+PEOSdata1(2))/(1.5d0*Pde1*Pde1)*ee*ee
     &          - (-9*PEOSdata1(1)+PEOSdata1(2))/(3*Pde1)*ee
       write(*,*)'Warning: out of the beginning of the table, PEOSL'
       write(*,*)'current e ', ee, ' smaller than Min e, ', Pe01
       write(*,*)'use quardratic exterpolation!'
      else if(ee.lt.Pe02) then
       call interpCubic(PEOSdata1, EOSdatasize, Pe01, Pde1, ee, PEOSL7)
      else if (ee.lt.Pe03) then
       call interpCubic(PEOSdata2, EOSdatasize, Pe02, Pde2, ee, PEOSL7)
      else if (ee.lt.Pe04) then
       call interpCubic(PEOSdata3, EOSdatasize, Pe03, Pde3, ee, PEOSL7)
      else if (ee.lt.Peend) then
       call interpCubic(PEOSdata4, EOSdatasize, Pe04, Pde4, ee, PEOSL7)
      else
       PEOSL7 = 0.0d0
       write(*,*) 'Warning: out of the end of table, PEOSL'
       write(*,*) 'current e ', ee, ' larger than  Max e, ', Peend
      endif

      return
      end
      
      Double Precision Function SEOSL7(ee)  ! for lattice P(e)
      Implicit none
       
!=======list of parameters===================================================
      Integer, Parameter :: EOSdatasize=501
!=======list of parameters end===============================================
      double precision :: ee
      double precision :: Se01, Sde1  !energy density data structure parametersfor data part 1
      double precision :: Se02, Sde2  !energy density data structure parameters for data part 2
      double precision :: Se03, Sde3  !energy density data structure parameters for data part 3
      double precision :: Se04, Sde4  !energy density data structure parameters for data part 4
      double precision :: Seend       !the maximum energy density in the table
      double precision :: SEOSdata1(EOSdatasize),SEOSdata2(EOSdatasize),
     &                    SEOSdata3(EOSdatasize),SEOSdata4(EOSdatasize) !entropy

      common /EOSdataS/SEOSdata1,SEOSdata2,SEOSdata3,SEOSdata4


      common /EOSdatastructureS/ Se01, Sde1, Se02, Sde2, Se03, Sde3,
     &                          Se04, Sde4, Seend

      if (ee.lt.Se01) then
       SEOSL7 = ee*SEOSdata1(1)/Se01
       SEOSL7 = (-3*SEOSdata1(1)+SEOSdata1(2))/(1.5d0*Sde1*Sde1)*ee*ee
     &          - (-9*SEOSdata1(1)+SEOSdata1(2))/(3*Sde1)*ee
       write(*,*)'Warning: out of the beginning of the table, SEOSL'
       write(*,*)'current e ', ee, ' smaller than Min e, ', Se01
       write(*,*)'use quardratic exterpolation!'
      else if (ee.lt.Se02) then
       call interpCubic(SEOSdata1, EOSdatasize, Se01, Sde1, ee, SEOSL7)
      else if (ee.lt.Se03) then
       call interpCubic(SEOSdata2, EOSdatasize, Se02, Sde2, ee, SEOSL7)
      else if (ee.lt.Se04) then
       call interpCubic(SEOSdata3, EOSdatasize, Se03, Sde3, ee, SEOSL7)
      else if (ee.lt.Seend) then
       call interpCubic(SEOSdata4, EOSdatasize, Se04, Sde4, ee, SEOSL7)
      else
       SEOSL7 = 0.0d0
       write(*,*) 'Warning: out of the end of the table, SEOSL'
       write(*,*) 'current e ', ee, ' larger than  Max e, ', Seend
      endif

      return
      end
      
      Double Precision Function TEOSL7(ee)  ! for lattice P(e)
      Implicit none
       
!=======list of parameters===================================================
      Integer, Parameter :: EOSdatasize=501
!=======list of parameters end===============================================
      double precision :: ee
      double precision :: Te01, Tde1  !energy density data structure parametersfor data part 1
      double precision :: Te02, Tde2  !energy density data structure parameters for data part 2
      double precision :: Te03, Tde3  !energy density data structure parameters for data part 3
      double precision :: Te04, Tde4  !energy density data structure parameters for data part 4
      double precision :: Teend   !the maximum energy density in the table
      double precision :: TEOSdata1(EOSdatasize),TEOSdata2(EOSdatasize),
     &                    TEOSdata3(EOSdatasize),TEOSdata4(EOSdatasize) !Temperature

      common /EOSdataT/ TEOSdata1,TEOSdata2,TEOSdata3,TEOSdata4

      common /EOSdatastructureT/ Te01, Tde1, Te02, Tde2, Te03, Tde3,
     &                          Te04, Tde4, Teend

      if (ee.lt.Te01) then
C       TEOSL7 = ee*TEOSdata1(1)/Te01
       TEOSL7 = (-3*TEOSdata1(1)+TEOSdata1(2))/(1.5d0*Tde1*Tde1)*ee*ee
     &          - (-9*TEOSdata1(1)+TEOSdata1(2))/(3*Tde1)*ee
       write(*,*)'Warning: out of the beginning of the table, TEOSL'
       write(*,*)'current e ', ee, ' smaller than Min e, ', Te01
       write(*,*)'use quardratic exterpolation!'
      else if (ee.lt.Te02) then
       call interpCubic(TEOSdata1, EOSdatasize, Te01, Tde1, ee, TEOSL7)
      else if (ee.lt.Te03) then
       call interpCubic(TEOSdata2, EOSdatasize, Te02, Tde2, ee, TEOSL7)
      else if (ee.lt.Te04) then
       call interpCubic(TEOSdata3, EOSdatasize, Te03, Tde3, ee, TEOSL7)
      else if (ee.lt.Teend) then
       call interpCubic(TEOSdata4, EOSdatasize, Te04, Tde4, ee, TEOSL7)
      else
       TEOSL7 = 0.0d0
       write(*,*) 'Warning: out of the end of the table, TEOSL'
       write(*,*) 'current e ', ee, ' larger than  Max e, ', Teend
      endif

      return
      end
  
CSHEN===end===========================================================
     
