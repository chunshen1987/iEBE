!==============================================================================!
!           Subroutines for OSCAR output format                                !
!     For the detail information about the OSCAR format, please                !
!     reference TECHQM website.                                                !
!                                                                              !
!     OSCAR output file for one hydro run is about 11 GB                       !
!                                                                              !
!     Author: Chun Shen                                                        !
!     Data: Oct. 30, 2010                                                      !
!                                                                              !
!     Bug fixed: add the missing Pi02 to the output and changed the            !
!                the output format                                             !
!                Date: Aug. 18, 2011                                           !
!                                                                              !
!                                                                              !
!==============================================================================!

      Subroutine OSCARheaderoutput(IWrite, NXPhy0, NXPhy, NYPhy0, NYPhy,
     &                             T0, DX, DY)

      Implicit none

CSHEN======parameters declearation===============================================
      double precision :: ViscousC, VisBeta
      double precision :: T0
      double precision :: DX, DY
      Integer :: IVisflag
      Integer :: NXPhy0, NXPhy, NYPhy0, NYPhy
      Integer :: IWrite
C==========OSCAR2008Header related parameters==============================
      character*12 OHVerN,OHHydType,OHhypsurf
      character*80 OHInit1,OHInit2,OHEOS,OHCharge,OHHyper,
     &             OHGEOM,OHGrid,OHViscosity1,OHViscosity2,
     &             OHComm,OHCOMM1,OHEND
      Integer :: IOHNt, IOHNx, IOHNy, IOHNz, IOHC, IOHD, IOHT
      double precision :: OT0, OT1, OX0, OX1, OY0, OY1
CSHEN======parameters declearation end===========================================


!====common blocks=========================================================
      Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Shear Viscosity

!====common blocks end=====================================================

C============Output OSCAR2008H HEADER====================================
       OHVerN = "OSCAR2008H"
       if (ViscousC .ne. 0.0d0) then
          OHHydType = "viscous"
       else
          OHHydType = "ideal"
       endif
       OHhypsurf = "history"
C       OHInit1   = "INIT: Glauber WN 70% Au+Au"
       OHInit1   = "INIT: CGC Initialization Au+Au"
       OHInit2   = "INIT: b=11.1, s0=109, t0=0.6, eta/s=0.20"
       OHEOS     = "EOS: EOSL-PCE"
       OHCharge  = "CHARGES: baryon"
       OHHyper   = "HYPER: full evolution"
       OHGEOM    = "GEOM: scaling2d"
       OHGrid    = "GRID: Euler"
       IOHNt = 801   !to be the same as ItimeOSCAR
       IOHNx = NXPhy-NXPhy0+1
       IOHNy = NYPhy-NYPhy0+1
       IOHNz = 0
       IOHC  = 0              !number of conserved charges
       IOHD  = 8              !number of dissipation parameters
       IOHT  = 3              !number of tranportation parameters
       OT0   = T0             !initial time
       OT1   = 16.6           !finial time(to be change in the end)
       OX0   = NXPhy0*DX      !lower bound of x grid (fm)
       OX1   = NXPhy*DX       !upper bound of x grid (fm)
       OY0   = NYPhy0*Dy      !lower bound of y grid (fm)
       OY1   = NYPhy*Dy       !upper bound of y grid (fm)
       OHViscosity1 = "VISCOSITY: shear viscousity only"
C       OHViscosity2 = "VISCOSITY: "
       OHCOMM    = "COMM: Diss: Pi00,Pi01,Pi11,Pi12,Pi22,Pi33,PPi"
       OHCOMM1   = "COMM: Transport: \eta/s, tau_pi, tau_Pi"
       OHEND     = "END_OF_HEADER"
       write(IWrite,"(3a12)")OHVerN, OHHydType, OHhypsurf  !Line 1
       write(IWrite,"(a80)")OHInit1                        !Line 2
       write(IWrite,"(a80)")OHInit2                        !Line 3
       write(IWrite,"(a80)")OHEOS                          !Line 4
       write(IWrite,"(a80)")OHcharge                       !Line 5
       write(IWrite,"(a80)")OHHyper                        !Line 6
       write(IWrite,"(a80)")OHGEOM                         !LINE 7
       write(IWrite,"(a80)")OHGrid                         !Line 8
       write(IWrite,"(7i5)")IOHNt,IOHNx,IOHNy,IOHNz,IOHC,IOHD,IOHT  !Line 9
       write(IWrite,"(6f8.3)")OT0,OT1,OX0,OX1,OY0,OY1               !Line 10
       write(IWrite,"(a80)")OHViscosity1                  !Line 11
C       write(IWrite,"(a80)")OHViscosity2
       write(IWrite,"(a80)")OHCOMM                         !Line 12
       write(IWrite,"(a80)")OHCOMM1                        !Line 13
       write(IWrite,"(a80)")OHEND                          !Line 14
C     stop
CSHEN========end=========================================================
      end

      Subroutine OSCARbodyoutput(IWrite, NX0, NX, NY0, NY, NZ0, NZ,
     &                           NXPhy0, NXPhy, NYPhy0, NYPhy, ITime,
     &                           Ed, PL, Temp, Vx, Vy,
     &                           Pi00, Pi01, Pi02,
     &                           Pi11, Pi12, Pi22, Pi33, PPi,
     &                           VRelaxT, VRelaxT0)

      Implicit none

CSHEN======parameters declearation===============================================
      Integer :: NXPhy0, NXPhy, NYPhy0, NYPhy
      Integer :: NX0, NX, NY0, NY, NZ0, NZ
      double precision, parameter :: HbarC=0.19733d0
      Integer :: I, J, K, L

      double precision :: ViscousC, VisBeta
      double precision :: ViscousCTemp
      Integer :: IVisflag
      Integer :: IWrite
      Integer :: ITime

C==========OSCAR2008H related parameters===================================
      double precision :: OEd, OPL, OTemp,
     &                    OPi00, OPi01, OPi02, OPi11, OPi12,
     &                    OPi22, OPi33, OPPi,
     &                    OEtas, OTaupi, OTaubPi
      Integer :: R_qgp
CSHEN=========end==========================================================


      double precision :: Ed(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: PL(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Vx(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Vy(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Temp(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi00(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi01(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi02(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi11(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi12(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi22(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: Pi33(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: PPi(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: VRelaxT(NX0:NX, NY0:NY, NZ0:NZ)
      double precision :: VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ)

CSHEN======parameters declearation end===========================================

      Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Shear Viscosity

      do 3007 K=NZ0,NZ
      do 3007 I=NXPhy0,NXPhy
      do 3007 J=NYPhy0,NYPhy
       OEd   = Ed(I,J,K)*HbarC
       OPL   = PL(I,J,K)*HbarC
       OTemp = Temp(I,J,K)*HbarC
       OPi00 = Pi00(I,J,K)*HbarC
       OPi01 = Pi01(I,J,K)*HbarC
       OPi02 = Pi02(I,J,K)*HbarC
       OPi11 = Pi11(I,J,K)*HbarC
       OPi12 = Pi12(I,J,K)*HbarC
       OPi22 = Pi22(I,J,K)*HbarC
       OPi33 = Pi33(I,J,K)*HbarC
       OPPi  = PPi(I,J,K)*HbarC
       if(IVisflag.eq.1) then
         OEtas = ViscousCTemp(Temp(I,J,K))
       else
         OEtas = ViscousC
       endif
       OTaupi  = VRelaxT(I,J,K)
       OTaubPi = VRelaxT0(I,J,K)
       if(OTemp .gt. 0.183) then
            R_qgp = 1
       else
            R_qgp = 0
       endif
       write(IWrite,3008)ITime-1, (I-NXPhy0), (J-NYPhy0),
     &                OEd, OPL, OTemp, R_qgp, Vx(I,J,K), Vy(I,J,K),
     &                OPi00, OPi01, OPi02, OPi11, OPi12,
     &                OPi22, OPi33,OPPi,OEtas,OTaupi,OTaubPi
3007  continue
3008  FORMAT(3i5,3f18.8,i5,13f18.8)

      end

      Subroutine OSCARredundantoutput(IWrite, NZ0, NZ, ITime,
     &                                NXPhy0, NXPhy, NYPhy0, NYPhy)
      Implicit none

      Integer :: NXPhy0, NXPhy, NYPhy0, NYPhy
      Integer :: NZ0, NZ
      Integer :: I, J, K
      Integer :: IWrite
      Integer :: ITime

C==========OSCAR2008H related parameters===================================
      double precision :: OEd, OPL, OTemp, OVx, OVy,
     &                    OPi00, OPi01, OPi02, OPi11, OPi12,
     &                    OPi22, OPi33, OPPi,
     &                    OEtas, OTaupi, OTaubPi
      Integer :: R_qgp
CSHEN=========end==========================================================


      do 3009 K=NZ0,NZ
      do 3009 I=NXPhy0,NXPhy
      do 3009 J=NYPhy0,NYPhy
       OEd   = 0.0d0
       OPL   = 0.0d0
       OTemp = 0.0d0
       OVx   = 0.0d0
       OVy   = 0.0d0
       OPi00 = 0.0d0
       OPi01 = 0.0d0
       OPi02 = 0.0d0
       OPi11 = 0.0d0
       OPi12 = 0.0d0
       OPi22 = 0.0d0
       OPi33 = 0.0d0
       OPPi  = 0.0d0
       OTaupi  = 0.0d0
       OTaubPi = 0.0d0
       OEtas = 0.0
       R_qgp = 0
       write(IWrite,3010)ITime-1, (I-NXPhy0), (J-NYPhy0),
     &                OEd, OPL, OTemp, R_qgp, OVx, OVy,
     &                OPi00, OPi01, OPi02, OPi11, OPi12,
     &                OPi22, OPi33,OPPi,OEtas,OTaupi,OTaubPi
3009  continue
3010  FORMAT(3i5,3f18.8,i5,13f18.8)
      end
