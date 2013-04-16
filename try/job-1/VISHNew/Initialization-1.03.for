! 02-17-2011: The code "Ed = Ed + 1D-10" added at the end of energy initialization.

!======== Pre-defined parameters =======================================

#define EOSDATALENGTH 155500 ! compile with -cpp option
#define CONSTPI 3.14159265

!=======================================================================

      Subroutine InitializeAll(NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0,NXPhy,NYPhy,T0,DX,DY,DZ,DT,MaxT,NDX,NDY,NDT,
     &  TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,PPI,PISc,XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Bd,Sd,Time,Temp0,Temp,CMu,T00,T01,T02,IAA,CofAA,PNEW,
     &  TEM0,ATEM0,Rj,EPS0,V10,V20,AEPS0,AV10,AV20,TFREEZ,TFLAG)

      Implicit Double Precision (A-H, O-Z)
      Dimension X(NX0:NX), Y(NY0:NY), Z(NZ0:NZ)

      Dimension TT00(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension ScT00(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension TT01(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension ScT01(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension TT02(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension ScT02(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Rj(NX0:NX, NY0:NY, NZ0:NZ)   ! density

      Dimension T00(NX0:NX, NY0:NY, NZ0:NZ)! ideal T00  energy momentum tensor
      Dimension T01(NX0:NX, NY0:NY, NZ0:NZ)! ideal T01
      Dimension T02(NX0:NX, NY0:NY, NZ0:NZ)! ideal T02


      Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT00(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT01(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT02(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT33(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT11(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT12(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT22(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)    !Bulk pressure
      Dimension PISc(NX0:NX, NY0:NY, NZ0:NZ)


      Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
      Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure
      Dimension Bd(NX0:NX, NY0:NY, NZ0:NZ) !net baryon density
      Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Dimension CMu(NX0:NX, NY0:NY, NZ0:NZ) !Local chemical potential
      Dimension IAA(NX0:NX, NY0:NY, NZ0:NZ)
      Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

      Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008
      Dimension etaTtp0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

      Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008
      Dimension XiTtP0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi in last time step


      Dimension BEd(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension CC(NX0:NX, NY0:NY, NZ0:NZ) !Check

CSHEN==========================================================================
C======output relaxation time for both shear and bulk viscosity================
      Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
      Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
CSHEN==========================================================================

CSHEN===EOS from tables========================================================
      Integer, Parameter :: RegEOSdatasize = EOSDATALENGTH  !converted EOS table size
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                 !total rows of energy density

CSHEN===EOS from tables end====================================================


C----------------------------------------------------------------------------------------------
      DIMENSION EPS0(NX0:NX,NY0:NY),EPS1(NX0:NX,NY0:NY) ! Energy density in previous and current step
      DIMENSION TEM0(NX0:NX,NY0:NY),TEM1(NX0:NX,NY0:NY) ! Temperature density in previous and current step

      DIMENSION V10(NX0:NX,NY0:NY),V20(NX0:NX,NY0:NY)   !velocity in X Y in Previous step
      DIMENSION V11(NX0:NX,NY0:NY),V21(NX0:NX,NY0:NY)   !velocity in X Y in current step

      DIMENSION AEPS0(NX0:NX, NY0:NY),AEPS1(NX0:NX, NY0:NY) ! Energy density in previous and current step
      DIMENSION ATEM0(NX0:NX, NY0:NY),ATEM1(NX0:NX, NY0:NY) ! Temperature density in previous and current step
      DIMENSION AV10(NX0:NX, NY0:NY),AV20(NX0:NX, NY0:NY)   !velocity in X Y in Previous step
      DIMENSION AV11(NX0:NX, NY0:NY),AV21(NX0:NX, NY0:NY)   !velocity in X Y in current step

      DIMENSION F0Pi00(NX0:NX,NY0:NY),FPi00(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi01(NX0:NX,NY0:NY),FPi01(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi02(NX0:NX,NY0:NY),FPi02(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi11(NX0:NX,NY0:NY),FPi11(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi12(NX0:NX,NY0:NY),FPi12(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi22(NX0:NX,NY0:NY),FPi22(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi33(NX0:NX,NY0:NY),FPi33(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step


      CHARACTER*60 EARTERM
      INTEGER TFLAG, EINS
      Parameter (EARTERM='EARLY')  ! 'EARLY' parameter for early ended the program after decoupling
      Parameter (pi=3.1415926d0)
      Parameter (gt=169.0d0/4.0d0)!total freedom of Quarks and Gluond  Nf=2.5  !change another in InitialES
      Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm


C--------------------------------------------------------------------------------------------------------
      parameter(NNEW=4)
      DIMENSION PNEW(NNEW)    !related to root finding
      COMMON /EOSSEL/ IEOS   !Type of EOS
      Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling redious
      common/Edec/Edec
      common/Edec1/Edec1

      Common /Nsm/ Nsm
      Common /Accu/Accu

      COMMON /Initialization/ IInit     !
      Common/IEOS2dec/ IEOS2dec
      Common/R0Aeps/ R0,Aeps

      Common /Timestep/ DT_1, DT_2
      Double Precision Time

      COMMON /IEin/ IEin     !  type of initialization  entropy/enrgy
! ---Zhi-Changes---
      Common /RxyBlock/ Rx2, Ry2
      Common /EK/ EK, HWN
      Double Precision sFactor
      Common /sFactor/ sFactor
! ---Zhi-End---

      Integer NN ! For initial anisotropy calcualtion
      Double Precision XX, YY, TotalE, XC, YC, angle, RR
      Double Precision XN(0:9), YN(0:9), Weight ! For initial anisotropy calcualtion
      Double Precision XNP(0:9), YNP(0:9), WeightP ! XN' and YN', the one using r^n in the weight

      Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Shear Viscosity

      Double Precision SEOSL7, PEOSL7, TEOSL7, SEOSL6
      Double Precision ss, ddt1, ddt2, ee1, ee2
      External SEOSL7


!   ---Zhi-End---

CSHEN======================================================================
C=============add chemical potential for EOS6
      Common /chempot/ XEKeos6, XDEeos6, ICOLNUMeos6,
     &                 INEeos6, XMUeos6(1:501,1:40),
     &                 XEeos6(1:501), XMe2eos6(1:501,1:40)
C=========XMUeos6 store the chemical potential data
C=========XEeos6 store the coresponding energy density
C=========XMe2eos6 store the second derivation of the chemical potential
C=========INEeos6 store the number of the row of chemical potential
CSHEN======================================================================

CSHEN======================================================================
C==========OSCAR2008H related parameters===================================
      Integer :: ItimeOSCAR=800     !output time steps for OSCAR2008H file
      Integer :: IOSCARWrite
      Logical :: IOSCAR            ! trigger for output OSCAR format hydro results, False: not output, True: output
      Common/OSCAR/ IOSCAR,IOSCARWrite
CSHEN=========end==========================================================
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata !CSHEN: for EOS from tables
      common /EOSdatastructure/ EOSe0, EOSde, EOSne


      Accu=3.0                  ! 3pt formula for derivative


!------------- Freezeout energy density and temperature -----------------

      IF (IEOS .EQ. 2.and.TFLAG.eq.0) THEN
         if(IEOS2dec.eq.0) then ! decouple to gluons
         TDEc   = ((EDEC-0.0)*120./(PI**2.)/169.
     &            *(HBARC**3.))**(.25) !changed by Evan
         TFREEZ = TDEc
         else ! decouple to pions
         TFREEZ=SQRT(HBARC*SQRT(10.0*HBARC*EDEC)/PI)
         end if
      Else if (IEOS.eq.4) Then
         ee     = EDEC/HbarC
         TFREEZ = TEIEOS4(ee)*HbarC
         eed    = EDec/Hbarc              !1/fm^4
      else if (IEOS.eq.5) then            !CSHEN for SM-EOS Q
         ee     = EDEC/HbarC              !1/fm^4
         TFREEZ = TEIEOS5pp(ee)*HbarC     !GeV
         eed    = EDEC/HbarC              !1/fm^4
      else if (IEOS.eq.6) then            !CSHEN New EOSL
         ee     = EDEC                    !GeV/fm^3
         TFREEZ = TEOSL6(ee)              !GeV
         eed    = EDec                    !GeV/fm^3
      else if (IEOS.eq.7) then            !CSHEN EOS from tables
         ee     = EDEC                    !GeV/fm^3
         TFREEZ = TEOSL7(ee)              !GeV
         eed    = Edec                    !GeV/fm^3
      end if

         eed    = EDec/Hbarc              !fm^(-4)

      Print *,'TFREEZ, Tde',TFREEZ,'GEV', TFREEZ/HbarC,'fm-1',
     &         Tde ,'fm-1'
      Print *,'EDec, ',EDec,'GEV', eed,'fm-1'

       Time = T0

CSHEN===================================================================
CSHEN=====Using a smaller time step for short initialization time \tau_0
       if (Time.lt.0.59) then
            DT = DT_2
            write(*,*) 'Using a smaller time step for tau<0.6, DT=',DT
       else
            DT = DT_1
       endif
CSHEN====END============================================================

!---------------- Four flow velocity initialization---------------------

        do 2560 K = NZ0,NZ
        do 2560 I = NXPhy0-3,NXPhy+3
        do 2560 J = NYPhy0-3,NYPhy+3

          U1(I,J,K)  = 0.0d0
          U2(I,J,K)  = 0.0d0
          U0(I,J,K)  = sqrt(1.0+U1(I,J,K)**2+U2(I,J,K)**2)

          PU1(I,J,K) = 0.0d0
          PU2(I,J,K) = 0.0d0
          PU0(I,J,K) = sqrt(1.0+PU1(I,J,K)**2+PU2(I,J,K)**2)
2560   continue


!------------------- Energy initialization -----------------------------

        If(IInit.eq.0) then ! Gaussian initial condition
            call initialES2Gaussian(Ed,Bd,Sd, DX,DY,DZ,DT,      !Gaussain initial condition
     &      NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
        elseif(IInit.eq.1) then
            call initialES2(Ed,Bd,Sd, DX,DY,DZ,DT,      !Origin Glauber initial condition
     &      NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
        elseif (IInit.eq.2) then
C====Input the initial condition from file====
!           Unit: fm^-4
          If (IEin==0) Then
            OPEN(2,file='Initial/InitialEd.dat',status='old')
            do 2562 I = NXPhy0,NXPhy
              read(2,*)  (Ed(I,J,NZ0), J=NYPhy0,NYPhy)
2562        continue
            Ed=Ed/HbarC
            close(2)
          Else If (IEin==1 .AND. IEOS==7) Then
            OPEN(2,file='Initial/InitialSd.dat',status='old') ! read from file first
            Do I = NXPhy0,NXPhy
              read(2,*)  (Sd(I,J,NZ0), J=NYPhy0,NYPhy)
            End Do
            Close(2)
            Do I = NXPhy0,NXPhy
            Do J = NYPhy0,NYPhy
              Call invertFunctionD(SEOSL7, 0D0, 315D0, 1D-3, 0D0,
     &                        Sd(I,J,NZ0), resultingEd)
              Ed(I,J,NZ0) = resultingEd/HbarC ! to fm^(-4)
            End Do
            End Do
          End If ! IEin==0
        else
            Print*,'Init=',IInit,' is not supported by this version.'! ---Zhi-Changes---
            Stop
      end if


      Ed = Ed + 1D-10


!---------- Then convert energy to pressure, entropy, temperature ------
      call EntropyTemp3 (Ed,PL, Temp,CMu,Sd,
     &         NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

!---------- Use sFactor ------------------------------------------------
! VER-1.03: add support for sFactor parameter
      If (IEOS==7) Then ! use sFactor
        Sd = Sd*sFactor
        Print*, "sFactor=", sFactor
        Do I = NXPhy0,NXPhy
        Do J = NYPhy0,NYPhy
          Call invertFunctionD(SEOSL7, 0D0, 315D0, 1D-3, 0D0,
     &                        Sd(I,J,NZ0), resultingEd)
          Ed(I,J,NZ0) = resultingEd/HbarC ! to fm^(-4)
        End Do
        End Do
        call EntropyTemp3 (Ed,PL, Temp,CMu,Sd,
     &         NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
      End If

      do 2571 K = NZ0,NZ
      do 2571 I = NX0,NX
      do 2571 J = NY0,NY
        Temp0(I,J,K)   = Temp(I,J,K)
        etaTtp0(I,J,K) = (Ed(I,J,K)+PL(I,J,K))*Temp(I,J,K)/(6.0*VisBeta) !extra (eta T)/tau_pi terms in I-S eqn 02/2008
 2571 continue


!--------------- Viscous terms initialization --------------------------

      If (ViscousC>1D-6) Then

        call TransportPi6( Pi00,Pi01,Pi02,Pi33, Pi11,Pi12,Pi22,
     &  PPI,Ed, Sd, PL, Temp, Temp0, U0,U1,U2,PU0,PU1,PU2, DX,DY,DZ,DT,
     &  NX0,NY0,NZ0, NX,NY,NZ, Time, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &  VRelaxT,VRelaxT0)


      End If

!CHANGES
!   ---Zhi-Changes---
!-------Regulate Pi(mu,nu) before adding it to T tensor
!      If (ViscousC>1D-6) Then
!        call regulatePi(Time,NX0,NY0,NZ0,NX,NY,NZ,
!     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
!     &  Ed,PL,PPI,
!     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33,Vx,Vy)
!      End If
!-------End of regulation---------
!       ---Zhi-End---

!------------- T(mu,nu) initialization ---------------------------------

      do 2570 K = NZ0,NZ
      do 2570 I = NXPhy0,NXPhy
      do 2570 J = NYPhy0,NYPhy

         ee = Ed(i,j,k)*HbarC ! [ee] = GeV/fm^3
         Cn = 0.0
         pp = PEPS(Cn,ee)
         ee = Ed(i,j,k)
         pp = pp/HbarC ! [pp] -> fm^(-4)
         BulkPr = PPI(I,J,K)
         epU0   = (ee+pp+BulkPr)*U0(i,j,k)

         TT00(i,j,k) = (epU0*U0(i,j,k)+Pi00(i,j,k)-pp-BulkPr)*Time
         TT01(i,j,k) = (epU0*U1(i,j,k)+Pi01(i,j,k))*Time
         TT02(i,j,k) = (epU0*U2(i,j,k)+Pi02(i,j,k))*Time
         Rj(i,j,k)   = 0.0


2570  continue

       call dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT, Stotal,StotalBv,StotalSv,
     &  Ed,PL,Bd,Sd,Temp0,Temp,CMu, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,PNEW,NNEW)

    !EPS0 = 1.0d0
      DO 2600 J = NYPhy0-2,NYPhy+2
      DO 2600 I = NXPhy0-2,NXPhy+2
        EPS0(I,J) = Ed(I,J,NZ0)*HbarC
        V10(I,J)  = U1(I,J,NZ0)/U0(I,J,NZ0)
        V20(I,J)  = U2(I,J,NZ0)/U0(I,J,NZ0)
        TEM0(I,J) = Temp(I,J,NZ0)*HbarC

        If (ViscousC>1D-6) Then

          F0Pi00(I,J) = Pi00(I,J,NZ0)
          F0Pi01(I,J) = Pi01(I,J,NZ0)
          F0Pi02(I,J) = Pi02(I,J,NZ0)
          F0Pi11(I,J) = Pi11(I,J,NZ0)
          F0Pi12(I,J) = Pi12(I,J,NZ0)
          F0Pi22(I,J) = Pi22(I,J,NZ0)
          F0Pi33(I,J) = Pi33(I,J,NZ0)
        End If

2600  CONTINUE

            DO 2610 J = NYPhy0-2,NYPhy+2
            DO 2610 I = NXPhy0-2,NXPhy+2
                AEPS0(I,J) = Ed(I,J,NZ0)
                AV10(I,J)  = U1(I,J,NZ0)/U0(I,J,NZ0)
                AV20(I,J)  = U2(I,J,NZ0)/U0(I,J,NZ0)
                ATEM0(I,J) = Temp(I,J,NZ0)
2610        CONTINUE

         Print*,'after initializtion, 2d Time  ',Time


       Print*, 'Ed',Ed(0,0,NZ0), Ed(30,30,NZ0),
     &           Ed(50,50,NZ0),  Ed(60,60,NZ0),
     &           Ed(70,70,NZ0),  Ed(80,80,NZ0)
C     &           Ed(90,90,NZ0),  Ed(100,100,NZ0),
C     &           Ed(110,110,NZ0),  Ed(120,120,NZ0),
C     &           Ed(130,130,NZ0)

C       write(*,'(f6.3, 12f10.5)')Time,PPI(30,30,NZ0),PPI(40,40,NZ0),
C     &           PPI(50,50,NZ0),  PPI(60,60,NZ0),
C     &           PPI(70,70,NZ0),  PPI(80,80,NZ0),
C     &           PPI(90,90,NZ0),  PPI(100,100,NZ0)

      End Subroutine
