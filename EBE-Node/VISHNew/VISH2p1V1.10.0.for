! Modified according to InputFun-1.26test3 (1.27RC), call with Vx and Vy in regulatePi
! Multiple regulations

C*****************************************************************************
C                                                                            *
C         Viscous Israel Stewart Hydrodynamics in 2+1-dimention (VISH2+1)    *
C                                                                            *
C                           Huichao Song                                     *
C                                                                            *
C                    The Ohio State University                               *
C                                                                            *
C*****************************************************************************


C    REFERENCES

C   [1] H.Song and U.Heinz, Phys. Lett. B658, 279 (2008);
C   [2] H.Song and U.Heinz, Phys. Rev. C77, 064901 (2008);
C   [3] H.Song and U.Heinz, Phys. Rev. C78, 024902 (2008);
C   [4] H.Song and U.Heinz, arXiv:0909.1549 [nucl-th];
C   [5] H.Song, Ph.D thesis 2009, arXiv:0908.3656 [nucl-th].

!*******************************************************************************
!     This program is further modified by Zhi
!     Find changes made by Zhi by searching ---Changes-Zhi---
!                                   ---Last modified on Jan 9, 2010
!*******************************************************************************


C===============================================================================
C   Main Features includes:
C   --- Use "full" I-S equations [3].
C   --- Includes shear[1-2] and bulk viscosity[4].
C   --- Use Optical Glauber Initialization.
C   --- Has several available EOS's: EOS I, SM-EOS Q, and EOS L [5].
C   --- Choose to freeze-out by constant energy density/temperature
C
C   For details, please refer to Song's Ph.D thesis[5] Chap.3 as well as Ref[1-4].
C===============================================================================
!
! For list of changes since version 1.0, see ChangeLog.txt
!
C===============================================================================


!======== Pre-defined parameters =======================================
! compile with -cpp option
#define EOSDATALENGTH 155500
#define EOSMUDATALENGTH 501
#define EOSMUMAXPARTICLE 40
#define CONSTPI 3.14159265

! change to 1 to ignore checks
#define silent_checkPi 0
#define outputPiviolation .false.
#define outputMovie .false.

#define echo_level 5

!=======================================================================

       Program Main

       Implicit Double Precision (A-H, O-Z)
       Integer LS
       Integer NXPhy0, NYPhy0
       Integer NXPhy, NYPhy
       Integer NX, NY, NZ
       Integer NX0,NY0, NZ0         ! dimension

CSHEN===========================================================================
       Logical :: IOSCAR=.false.            ! trigger for output OSCAR format hydro results, False: not output, True: output
       Integer :: IVisflag=0          ! Flag for temperature dependent eta/s, 0: constant, 1: temperature dependent, which is defined in function ViscousCTemp(T)
       Integer :: IOSCARWrite       ! output OCCAR file path number
CSHEN===========================================================================

       Integer :: IhydroJetoutput   ! Output control for hydro evolution history
       Common /hydroJetoutput/ IhydroJetoutput

CSHEN===EOS from tables========================================================
      Integer, Parameter :: RegEOSdatasize = EOSDATALENGTH !converted EOS table size
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      
CSHEN======================================================================
C=============add chemical potential for EOS-PCE
      Integer, Parameter :: RegEOSMudatasize = EOSMUDATALENGTH   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = EOSMUMAXPARTICLE        !maximum allowed stable particles for partically chemical equilibrium EOS
      Integer :: Inumparticle       !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne            !total rows in mu table
      double precision :: EOS_Mu_e0   !lowest energy density in mu table
      double precision :: EOS_Mu_de   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      integer :: Muflag
CSHEN===EOS from tables end====================================================

       Common /LS/ LS
       Common /R0Bdry/ R0Bdry

       COMMON /EOSSEL/ IEOS   !Type of EOS
       COMMON /IEin/ IEin     !  type of initialization  entropy/enrgy
       COMMON /Initialization/ IInit     ! type of initialization CGC/Glauber

       Common /AWNBC/ A,Si0 !A Nuclei Number  Si0, Cross Section for NN
       Common /EK/ EK, HWN  !EK(T0) constant related to energy density,HWN percent of Wounded Nucleon
       Common /thick/ TRo0, TEta, TRA  !Para in Nuclear Thickness Function
       Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Shear Viscosity
       Common /ViscousBulk/ Visbulk, BulkTau,IRelaxBulk  ! Related to bulk Visousity

       Common /bb/ b  !impact parameter
       Common/dxdy/ ddx, ddy
       Common /TT0/ TT0   !T0

       Common /Timestep/ DT_1, DT_2

       common/Edec/Edec    !decoupling energy density
       common/Prefreezeout/Edec0, Ifreez
       Common/IEOS2dec/ IEOS2dec  ! IEOS=2 decouple by gluon/pion

       Common/OSCAR/ IOSCAR, IOSCARWrite  !CSHEN: for OSCAR output

       common /EOSdata/PEOSdata, SEOSdata, TEOSdata !CSHEN: for EOS from tables
       common /EOSdatastructure/ EOSe0, EOSde, EOSne
       common /EOSMudata/MuEOSdata, IMuflag
       common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne, 
     &                            Inumparticle

       parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm
       Character Cha

       Integer MaxT

      Double Precision :: cpu_start, cpu_end ! timing the application



***********************************************************************************
! ---Zhi-Changes---
      Common /RxyBlock/ Rx2,Ry2

      Common /DXY/ DX,DY

      Integer NDX, NDY, NDT ! used in Freeze-out subroutine
      Common /NXYTD/ NDX, NDY, NDT

      Double Precision T0 ! initial time tau_0
      Common /T0/ T0

C *******************************J.Liu changes*******************************

      Integer InitialURead
      Common/LDInitial/ InitialURead  ! IintURead =1 read initial velocity profile
C *******************************J.Liu changes end***************************

      call prepareInputFun() ! this is the initialization function in InputFun.for

      Print *, "Read parameters from Vishydro.inp file."

!----------Start of reading parameters from file------------------------
C========= Inputting Parameters ===========================================
      OPEN(1,FILE='Vishydro.inp',STATUS='OLD')

      READ(1,*) T0          !initial time (fm/c)
      READ(1,*) Edec        !Decoupling energy density   (GeV/fm^3)
      READ(1,*) IEOS        !type of EOS  (EOSQ:0  EOSI: 2, SM-EOSQ: 5, EOSL: 4(Katz05 data) )
      READ(1,*) IEOS2dec    !if IEOS=2 (EOSI),  0: decouple by gluon 1: decouple by pions
C----------------------------------------------------------------------
      Read(1,*) Cha
      Read(1,*) Cha

      READ(1,*) IInit       !initialization by Glauber or CGC  (0:Gaussian, 1:Glauber, 2:read from file )
      READ(1,*) IEin        !type of initialization  0:initialize by energy density 1: initialize by entropy density

C------- Parameter for viscous coeficients  ------------------------
      Read(1,*) Cha
      Read(1,*) Cha

      READ(1,*) ViscousC    ! eta/s  (constant)    (0: no shear vis )
      READ(1,*) VisBeta     !\tau_Pi=VisBeta*6.0\eta /(ST)
      READ(1,*) VisBulk     !VisBulk=C;  Xi/s= C* (Xi/S)_min  (C=0, no bulk vis; or C>1 )
      READ(1,*) IRelaxBulk  !type of bulk relaxation time (0: critical slowing down; 1: contant Relax Time
                            !2: \tau_PI=1.5/(2\piT))
      READ(1,*) BulkTau !constant relaxation time (fm/c) (require input IRelaxBulk=1)

C------- Parameter for freeze-out  ------------------------
      Read(1,*) Cha
      Read(1,*) Cha

      READ(1,*) NDX,NDY   ! freeze-out step in x and y dirction
      READ(1,*) NDT       ! freeze-out step in \tau dirction


C------ Lattice size and boundary R0 ----------------------------
      Read(1,*) Cha
      Read(1,*) Cha

      READ(1,*) DT_1      ! time step
      READ(1,*) LS        ! lattice size
      READ(1,*) R0Bdry    ! <x^2> and <y^2> for Gaussian initial condition

C------ Some uncommon parameters ----------------------------
      Read(1,*) Cha
      Read(1,*) Cha

      READ(1,*) Rx2,Ry2   ! <x^2> and <y^2> for Gaussian initial condition

C------ freeze out at first time below Edec  ----------------------------
      Read(1,*) Cha
      Read(1,*) Ifreez
      Read(1,*) Edec0   

C------ output hydro evolution file  ----------------------------
      Read(1,*) Cha
      Read(1,*) IhydroJetoutput

C ***************************J.Liu changes***************************
C------- Parameters for initial profile from Laudan matching-----------------------
      Read(1,*) Cha
      Read(1,*) Cha
      Read(1,*) InitialURead  
      CLOSE(1)
C ***************************J.Liu changes end***************************
C===========================================================================

      DX=0.1d0
      DY=0.1d0


!-----------End of reading parameters from file-------------------------

      Print *, "Now read parameters specified from CML"

      Call readInputFromCML2() ! check CML to see if there are any modifications on parameters
      Write (*,*) "Have:", "IEOS=", IEOS, "A=", A, ! write out parameter for a check
     &    "IInit=", IInit, "dT=", dT_1,
     &    "eta/s=",ViscousC,"b=",b,"Rx2=",Rx2,"Ry2=",Ry2,
     &    "EK=", EK, "tau0=", T0, "EDec=", EDec, ! energy unit: GeV/fm^3
     &    "LS=", LS, "R0Bdry", R0Bdry, "VisBeta=", VisBeta,
     &    "DX=", DX, "DY=", DY, "DT_1=", DT_1,
     &    "NDX=", NDX, "NDY=", NDY, "NDT=", NDT,
     &    "IhydroJetoutput=", IhydroJetoutput


      ddx=dx
      ddy=dy
      DZ=0.01d0


CSHEN======================================================================
C====Change to a smaller time step for small \tau_0 case ==================
      DT_2 = DT_1/10.0

      dT  = DT_1
CSHEN===END================================================================

      MaxT = 40.0/DT_1

      EK  = EK/Hbarc  !center energy density unit convert to fm^-4
      TT0 = T0

      CALL UNLINK ('results/surface.dat')
      CALL UNLINK ('results/decdat2.dat')
      CALL UNLINK ('results/decdat_mu.dat')   !CSHEN chemical potential for PCE

      OPEN(98,FILE='results/surface.dat',FORM='FORMATTED',
     &        STATUS='REPLACE')
      OPEN(99,FILE='results/decdat2.dat',FORM='FORMATTED',
     &        STATUS='REPLACE')

CSHEN======================================================================
C======output the chemical potential information at freeze out surface.====
      open(81,FILE='results/decdat_mu.dat',FORM='FORMATTED',
     &     STATUS='REPLACE')
      open(82,FILE='results/SurfaceX.dat',FORM='FORMATTED',
     &     STATUS='REPLACE')         !output the freeze-out surface along x-axis
      open(83,FILE='results/SurfaceY.dat',FORM='FORMATTED',
     &     STATUS='REPLACE')         !output the freeze-out surface along y-axis
      open(91,File='results/Temp.dat',status='REPLACE')
      open(93,File='results/Temp_evo.dat',status='REPLACE')
      open(90,File='results/APi.dat',status='REPLACE')
      open(89,File='results/AScource.dat',status='REPLACE')
      open(88,File='results/AScource2.dat',status='REPLACE')
      Open(377, FILE='results/ecc-evolution.dat',STATUS='REPLACE') ! anisotropy moments evolution

      !Open(3771,FILE='movie/DPc.dat',STATUS='REPLACE')
      !Close(3771)
      if (outputMovie) then 
         open(3773,FILE='movie/Evolution.dat',STATUS='REPLACE')
      endif
      !Open(3774,FILE='movie/U1.dat',STATUS='REPLACE')
      !Close(3774)
      !Open(3775,FILE='movie/Vy.dat',STATUS='REPLACE')
      !Close(3775)
      !Open(3776,FILE='movie/U2.dat',STATUS='REPLACE')
      !Close(3776)

      If (IOSCAR) Then
        IOSCARWrite = 43
        open(IOSCARWrite,File='results/OSCAR2008H.dat',Form='Formatted',
     &     status='REPLACE')   !output standard format VISH2+1 results for further hydro movie and jet quenching purpose
      EndIf
      if (outputPiviolation) then
         open(583, File='results/piViolation.dat',
     &        Form='Formatted', status='REPLACE')
         open(584, File='results/BulkpiViolation.dat',
     &        Form='Formatted', status='REPLACE')
      endif
CSHEN======================================================================

      open(92,File='results/anisotropy.dat',status='REPLACE')

      call InputEOSParameter !  EOS related parameter

      if(IEOS.eq.7) then      !CSHEN: input EOS from table
            call InputRegulatedEOS
      endif

      NXPhy0=-LS
      NYPhy0=-LS
      NXPhy=LS
      NYPhy=LS
      NX=NXPhy+5
      NY=NYPhy+5
      NZ=1
      NX0=NXPhy0-5
      NY0=NYPhy0-5
      NZ0=1

CSHEN======output OSCAR file Header=========================================
      if(IOSCAR) then
            call OSCARheaderoutput(IOSCARWrite, NXPhy0, NXPhy,
     &                             NYPhy0, NYPhy, T0, DX, DY)
      endif
CSHEN======output OSCAR file Header end=====================================

CSHEN======set up output file for hydro evolution history===================
      if(IhydroJetoutput .eq. 1) then
         Call setHydroFiles(NX0, NX, DX, 5, NY0, NY, DY, 5, T0, DT, 5)
      endif

      CALL CPU_TIME(cpu_start) !Tic
      Call Mainpro(NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NYPhy0,
     &          NXPhy,NYPhy,T0,DX,DY,DZ,DT,MaxT,NDX,NDY,NDT)   ! main program
      CALL CPU_TIME(cpu_end) ! Toc
      Print *, "Finished in", cpu_end-cpu_start, "seconds."

      Close(98)
      Close(99)
      Close(93)
      Close(81)
      Close(82)
      Close(83)
      close(IOSCARWrite)
      Close(377)
      if(outputMovie) then 
         close(3773)
      endif

      End
!-----------------------------------------------------------------------





C######################################################################################
      Subroutine Mainpro(NX0,NY0,NZ0, NX,NY,NZ,NXPhy0,NYPhy0,
     &         NXPhy,NYPhy,T0,DX,DY,DZ,DT,MaxT,NDX,NDY,NDT)

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

      !Dimension Vperp(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT00(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT01(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT02(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension PScT33(NX0:NX, NY0:NY, NZ0:NZ)

      Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi11A(NX0:NX, NY0:NY, NZ0:NZ)    !Pi11_attempt, used only when v_perp=0
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

      ! backups

      Double Precision Pi00Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi01Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi02Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi11Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi12Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi22Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Pi33Regulated(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PPIRegulated(NX0:NX,NY0:NY,NZ0:NZ)

CSHEN==========================================================================
C======output relaxation time for both shear and bulk viscosity================
      Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient relaxation time
      Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient relaxation time
CSHEN==========================================================================



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
      DIMENSION F0PPI(NX0:NX,NY0:NY),FPPI(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step


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
      Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling radius
      common/Edec/Edec
      common/Prefreezeout/Edec0, Ifreez
      common/Edec1/Edec1

      Common /Nsm/ Nsm
      Common /Accu/Accu

      COMMON /Initialization/ IInit     !
      Common/IEOS2dec/ IEOS2dec
      Common/R0Aeps/ R0,Aeps
      Common /R0Bdry/ R0Bdry

      Common /Timestep/ DT_1, DT_2

      COMMON /IEin/ IEin     !  type of initialization  entropy/energy
! ---Zhi-Changes---
      Common /RxyBlock/ Rx2, Ry2
      Common /EK/ EK, HWN
      Character (len=80) filename_buffer
      Integer r_power

      Integer NN ! For initial anisotropy calcualtion
      Double Precision XX, YY, TotalE, XC, YC, angle, RR
      Double Precision XN(0:9), YN(0:9), Weight ! For initial anisotropy calculation
      Double Precision XNP(0:9), YNP(0:9), WeightP ! XN' and YN', the one using r^n in the weight

      Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Shear Viscosity
      Common /ViscousBulk/ Visbulk, BulkTau,IRelaxBulk  ! Related to bulk Viscosity
      Common /sFactor/ sFactor

      Double Precision SEOSL7, PEOSL7, TEOSL7, SEOSL6
      Double Precision ss, ddt1, ddt2, ee1, ee2
      External SEOSL7
      Integer iRegulateCounter, iRegulateCounterBulkPi

!   ---Zhi-End---


CSHEN======================================================================
C==========OSCAR2008H related parameters===================================
      Integer :: ItimeOSCAR=800     !output time steps for OSCAR2008H file
      Integer :: IOSCARWrite
      Logical :: IOSCAR            ! trigger for output OSCAR format hydro results, False: not output, True: output
      Common/OSCAR/ IOSCAR,IOSCARWrite

      Integer :: Tau_idx
CSHEN=========end==========================================================

      Integer :: IhydroJetoutput
      Common/hydroJetoutput/ IhydroJetoutput

CSHEN===EOS from tables========================================================
      Integer, Parameter :: RegEOSdatasize = EOSDATALENGTH  !converted EOS table size
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                 !total rows of energy density

      Integer, Parameter :: RegEOSMudatasize = EOSMUDATALENGTH   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = EOSMUMAXPARTICLE        !maximum allowed stable particles for partially chemical equilibrium EOS
      Integer :: Inumparticle       !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne            !total rows in mu table
      double precision :: EOS_Mu_e0   !lowest energy density in mu table
      double precision :: EOS_Mu_de   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      Integer :: IMuflag
      double precision :: XMufreeze(1:IMax_mu)
CSHEN===EOS from tables end====================================================
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata !CSHEN: for EOS from tables
      common /EOSdatastructure/ EOSe0, EOSde, EOSne
      common /EOSMudata/MuEOSdata, IMuflag
      common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne, 
     &                            Inumparticle

      TFLAG = 0
      IMuflag = 0
      Edec1 = Edec

!=======================================================================
!============ Initialization ===========================================
      Call InitializeAll(NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0,NXPhy,NYPhy,T0,DX,DY,DZ,DT,MaxT,NDX,NDY,NDT,
     &  TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,PPI,PISc,XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Bd,Sd,Time,Temp0,Temp,CMu,T00,T01,T02,IAA,CofAA,PNEW,
     &  TEM0,ATEM0,Rj,EPS0,V10,V20,AEPS0,AV10,AV20,TFREEZ,TFLAG)

!=======================================================================

      !output hydro parameters for current run
      Open(111,File='./results/VISH2p1_tec.dat',status='Replace')
       Write(111,'(A,T10,A3,I1,T30,A)')"IEOS", " = ", IEOS, 
     &                                            "!IEOS   type for EOS"
       Write(111,'(A,T10,A3,F4.2,T30,A)')"Tau0", " = ", T0, 
     &                                "!tau0   initial time(unit: fm/c)"
       Write(111,'(A,T10,A3,F5.3,T30,A)')"Edec", " = ", EDec, 
     &              "!Edec   Decoupling energy density (unit: GeV/fm^3)"
       Write(111,'(A,T10,A3,F5.3,T30,A)')"Tdec", " = ", TFREEZ, 
     &                                  "!T_fo   freeze out temperature"
       Write(111,'(A1)')"*"
       Write(111,'(A,T10,A3,I1,T30,A)')"IInit", " = ", IInit, 
     &                                          "!Initialization method"
       Write(111,'(A,T10,A3,E13.7,T30,A)')"Norm"," = ", sFactor,
     &                       "!normalization factor for initial profile"
       Write(111,'(A1)')"*"
       Write(111,'(A,T10,A3,F5.3,T30,A)')"ViscousC", " = ", ViscousC, 
     &                                        "!ViscousC  eta/s(const.)"
       Write(111,'(A,T10,A3,F5.3,T30,A)')"VisBeta", " = ", VisBeta, 
     &                        "!relaxation constant for shear viscosity"
       Write(111,'(A,T10,A3,I1,T30,A)')"IVisflag", " = ", IVisflag, 
     &        "!flag for temperature dependent eta/s(T) (=0 for const.)"
       Write(111,'(A,T10,A3,F5.3,T30,A)')"VisBulk", " = ", VisBulk, 
     &                                  "!Bulk viscosity zeta/s(const.)"
       Write(111,'(A,T10,A3,I1,T30,A)')"IRelaxBulk", " = ", IRelaxBulk, 
     &                          "!relaxtion constant for bulk viscosity"
       Write(111,'(A1)')"*"
       Write(111,'(A,T10,A3,F4.2,T30,A)')"DT", " = ", DT_1, 
     &                                                 "!dT   time step"
       Write(111,'(A,T10,A3,F4.2,T30,A)')"DX", " = ", DX, 
     &                                   "!DX   lattice spacial spacing"
       Write(111,'(A,T10,A3,F4.2,T30,A)')"DY", " = ", DY, 
     &                                   "!DY   lattice spacial spacing"
       Write(111,'(A,T10,A3,I3,T30,A)')"ILS", " = ", NXPhy, 
     &                           "!lattice size (in positive direction)"
       Write(111,'(A,T10,A3,I1,T30,A)')"NDX", " = ", NDX,
     &                      "!NXD   freeze out skip step in x direction"
       Write(111,'(A,T10,A3,I1,T30,A)')"NDY", " = ", NDY,
     &                      "!NYD   freeze out skip step in y direction"
       Write(111,'(A,T10,A3,I1,T30,A)')"NDT", " = ", NDT, 
     &                              "!NTD   freeze out skip step in tau"
       Write(111,'(A,T10,A3,F4.1,T30,A)')"R0Bdry", " = ", R0Bdry, 
     &                                                    "!boundary R0"
      Close(111)

!     Output entropy and energy profiles:
      Open(111,File='./Initial/init-entropy.dat',status='Replace')
      Do I=NXPhy0,NXPhy
      Do J=NYPhy0,NYPhy
        Write (111,'(E20.8)',Advance="NO") Sd(I,J,1)
      End Do
        Write (111,*)
      End Do
      Close(111)

      Open(111,File='./Initial/init-energy.dat',status='Replace')
      Do I=NXPhy0,NXPhy
      Do J=NYPhy0,NYPhy
        Write (111,'(E20.8)',Advance="NO") Ed(I,J,1)
      End Do
        Write (111,*)
      End Do
      Close(111)
    
      ! output fluid cells outside the freeze out surface at initial time
      if(Ifreez .ne. 0) then
      DA0=(dx*NDX)*(dy*NDY)
      DA1=0.0
      DA2=0.0          

      DO J = NYPhy0, NYPhy
        IF (MOD(J, NDY) .NE. 0) cycle
        YY = J * DY
      DO I = NXPhy0, NXPHY
        IF (MOD(I, NDX) .NE. 0) cycle
        XX = I * DX
        If(Ed(I,J,1)*HbarC.lt.Edec .and. Ed(I,J,1)*HbarC.ge.Edec0) then
          VZCM = Vx(I,J,1)
          VRCM = Vy(I,J,1)
          BN = 0.0
          BAMU = 0.0                   
          SMU  = 0.0                  
          PDec2 = PL(I,J,1)

          CPi00 = Pi00(I,J,1)    
          CPi01 = Pi01(I,J,1)    
          CPi02 = Pi02(I,J,1)   
          CPi11 = Pi11(I,J,1)    
          CPi12 = Pi12(I,J,1)    
          CPi22 = Pi22(I,J,1)     
          CPi33 = Pi33(I,J,1)    

          WRITE(99,'(99E20.8E3)') Time,DA0,DA1,DA2,VZCM,VRCM,
     &                  Ed(I,J,1)*HbarC,BN,
     &                  Temp(I,J,1)*HbarC,BAMU,SMU, PDec2, 
     &                  CPi33*HbarC,
     &                  CPi00*HbarC,CPi01*HbarC,CPi02*HbarC,
     &                  CPi11*HbarC,CPi12*HbarC,CPi22*HbarC
          
          TM = Time
          XM = XX
          YM = YY

          WRITE(98,'(99E20.8)') Time, TM, XM,YM,
     &                  Sqrt(XM**2+YM**2), R0,Aeps    !surface    
        

          EDEC2 = Ed(I,J,1)*HbarC
          IF (IEOS.EQ.7) then
             do L=1,Inumparticle
                call interpCubic(MuEOSdata(:,L), RegEOSMudatasize, 
     &                   EOS_Mu_e0, EOS_Mu_de, EDEC2, XMufreezetemp)
                XMufreeze(L) = XMufreezetemp
             end do
             write(81,'(E15.6)', Advance='NO') EDEC2
             do L=1,Inumparticle
               write(81,'(E15.6)', Advance='NO') XMufreeze(L)
             enddo
             write(81,*)
          end if
        End If
       End Do
       End Do
       End If  !Ifreez
CSHEN=====================================================================
!=======================================================================
!=======================================================================
!----------- Initial profile related output starts here ----------------

      Print *, "-----------------------------------------------------"
      Print *, "Info for initial energy density profile..."

      ! First, check center:
      Print *
      Print *, "First, check where the center is:"
      XC = 0D0
      YC = 0D0
      TotalE = 0D0
      Do K = NZ0, NZ
      Do I = NXPhy0, NXPhy
      Do J = NYPhy0, NYPhy
        !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
        XX = I*DX
        YY = J*DY
        XC = XC + XX*Ed(I,J,K)
        YC = YC + YY*Ed(I,J,K)
        TotalE = TotalE + Ed(I,J,K)
      End Do
      End Do
      End Do
      XC = XC/TotalE
      YC = YC/TotalE
      Print *, "Centered at X=", XC, "Y=", YC

      ! Next, calculate eccentricity moments for different N:
      Print *
      Print *, "Next, moments that give minor axes (see 1007.5469)."
      Print *, "They are calcualted as -X(N)/X(0), -Y(N)/X(0)."
      Print *, "They are caluclated assuming the center is origin."
      Print *

      ! first, conventional eccentricities
      Open(375, FILE='results/ecc-init.dat',STATUS='REPLACE')
      Do NN = 1, 9
        Weight = 0D0
        WeightP = 0D0
        XN(NN) = 0D0
        YN(NN) = 0D0
        XNP(NN) = 0D0
        YNP(NN) = 0D0
        Do K = NZ0, NZ
        Do I = NXPhy0, NXPhy
        Do J = NYPhy0, NYPhy
          !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
          XX = I*DX - XC ! shift to the real center
          YY = J*DY - YC

          angle = NN*atan2(YY,XX)
          RR = sqrt(XX*XX+YY*YY)

          ! eccentricity, defined using r^2 as weight function:
          Weight=Weight+RR*RR*Ed(I,J,K) ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
          XN(NN)=XN(NN)+RR*RR*cos(angle)*Ed(I,J,K)
          YN(NN)=YN(NN)+RR*RR*sin(angle)*Ed(I,J,K)
          ! eccentricity, defined using r^n as weight function:
          WeightP = WeightP+RR**NN*Ed(I,J,K)
          XNP(NN) = XNP(NN)+RR**NN*cos(angle)*Ed(I,J,K)
          YNP(NN) = YNP(NN)+RR**NN*sin(angle)*Ed(I,J,K)
        End Do
        End Do
        End Do
        Print *, "The",NN,"-th eccentricity moment:"
        Print *, "X:", -XN(NN)/Weight, "Y:", -YN(NN)/Weight
        print *, "   Norm:", sqrt(XN(NN)*XN(NN)+YN(NN)*YN(NN))/Weight
        Print *, "X':", -XNP(NN)/WeightP, "Y':", -YNP(NN)/WeightP
        Print *, "   Norm':",
     &    sqrt(XNP(NN)*XNP(NN)+YNP(NN)*YNP(NN))/WeightP
        !Print *, "WeightP=", WeightP
        Write (375,'(6E20.8)') -XN(NN)/Weight, -YN(NN)/Weight, ! Note that I use minor axis to define eccentricity
     &            sqrt(XN(NN)*XN(NN)+YN(NN)*YN(NN))/Weight,
     &            -XNP(NN)/WeightP, -YNP(NN)/WeightP,
     &            sqrt(XNP(NN)*XNP(NN)+YNP(NN)*YNP(NN))/WeightP
      End Do ! NN=1, 9
      Close(375)

      ! next, repeat but using different weight-angle combinations
      Do r_power = 0, 9
      write(filename_buffer, '(a25,i1,a4)')
     &     'results/ecc-init-r_power-', r_power, '.dat'
      Open(375, FILE=filename_buffer,STATUS='REPLACE')
      Do NN = 1, 9
        Weight = 0D0
        WeightP = 0D0
        XN(NN) = 0D0
        YN(NN) = 0D0
        XNP(NN) = 0D0
        YNP(NN) = 0D0
        Do K = NZ0, NZ
        Do I = NXPhy0, NXPhy
        Do J = NYPhy0, NYPhy
          !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
          XX = I*DX - XC ! shift to the real center
          YY = J*DY - YC

          angle = NN*atan2(YY,XX)
          RR = sqrt(XX*XX+YY*YY)

          ! eccentricity, defined using r^2 as weight function:
          Weight=Weight+RR**r_power*Ed(I,J,K) ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
          XN(NN)=XN(NN)+RR**r_power*cos(angle)*Ed(I,J,K)
          YN(NN)=YN(NN)+RR**r_power*sin(angle)*Ed(I,J,K)
        End Do
        End Do
        End Do
        Write (375,'(4E20.8)') -XN(NN)/Weight, -YN(NN)/Weight, ! Note that I use minor axis to define eccentricity
     &            sqrt(XN(NN)*XN(NN)+YN(NN)*YN(NN))/Weight, Weight*DX*DY
      End Do ! NN=1, 9
      Close(375)
      End Do ! r_power

      Print *, "-----------------------------------------------------"

!-----------------------------------------------------------------------
      ! Then repeat above eccentricity calculation for entropy-defined eccentricities
      ! first center
      XC = 0D0
      YC = 0D0
      TotalE = 0D0
      Do K = NZ0, NZ
      Do I = NXPhy0, NXPhy
      Do J = NYPhy0, NYPhy
        !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
        XX = I*DX
        YY = J*DY
        XC = XC + XX*Sd(I,J,K)
        YC = YC + YY*Sd(I,J,K)
        TotalE = TotalE + Sd(I,J,K)
      End Do
      End Do
      End Do
      XC = XC/TotalE
      YC = YC/TotalE

      ! first conventional eccentricities
      Open(379, FILE='results/ecc-init-sd.dat',STATUS='REPLACE')
      Do NN = 1, 9
        Weight = 0D0
        WeightP = 0D0
        XN(NN) = 0D0
        YN(NN) = 0D0
        XNP(NN) = 0D0
        YNP(NN) = 0D0
        Do K = NZ0, NZ
        Do I = NXPhy0, NXPhy
        Do J = NYPhy0, NYPhy
          !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
          XX = I*DX - XC ! shift to the real center
          YY = J*DY - YC

          angle = NN*atan2(YY,XX)
          RR = sqrt(XX*XX+YY*YY)

          ! eccentricity, defined using r^2 as weight function:
          Weight=Weight+RR*RR*Sd(I,J,K) ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
          XN(NN)=XN(NN)+RR*RR*cos(angle)*Sd(I,J,K)
          YN(NN)=YN(NN)+RR*RR*sin(angle)*Sd(I,J,K)
          ! eccentricity, defined using r^n as weight function:
          WeightP = WeightP+RR**NN*Sd(I,J,K)
          XNP(NN) = XNP(NN)+RR**NN*cos(angle)*Sd(I,J,K)
          YNP(NN) = YNP(NN)+RR**NN*sin(angle)*Sd(I,J,K)
        End Do
        End Do
        End Do
        !Print *, "WeightP=", WeightP
        Write (379,381) -XN(NN)/Weight, -YN(NN)/Weight, ! Note that I use minor axis to define eccentricity
     &            sqrt(XN(NN)*XN(NN)+YN(NN)*YN(NN))/Weight,
     &            -XNP(NN)/WeightP, -YNP(NN)/WeightP,
     &            sqrt(XNP(NN)*XNP(NN)+YNP(NN)*YNP(NN))/WeightP
 381    Format(6(E20.8))
      End Do ! NN=1, 9
      Close(379)

      ! next, repeat but using different weight-angle combinations
      Do r_power = 0, 9
      write(filename_buffer, '(a28,i1,a4)')
     &     'results/ecc-init-sd-r_power-', r_power, '.dat'
      Open(379, FILE=filename_buffer,STATUS='REPLACE')
      Do NN = 1, 9
        Weight = 0D0
        WeightP = 0D0
        XN(NN) = 0D0
        YN(NN) = 0D0
        XNP(NN) = 0D0
        YNP(NN) = 0D0
        Do K = NZ0, NZ
        Do I = NXPhy0, NXPhy
        Do J = NYPhy0, NYPhy
          !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
          XX = I*DX - XC ! shift to the real center
          YY = J*DY - YC

          angle = NN*atan2(YY,XX)
          RR = sqrt(XX*XX+YY*YY)

          ! eccentricity, defined using r^2 as weight function:
          Weight=Weight+RR**r_power*Sd(I,J,K) ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
          XN(NN)=XN(NN)+RR**r_power*cos(angle)*Sd(I,J,K)
          YN(NN)=YN(NN)+RR**r_power*sin(angle)*Sd(I,J,K)
        End Do
        End Do
        End Do
        Write (379,'(4E20.8)') -XN(NN)/Weight, -YN(NN)/Weight, ! Note that I use minor axis to define eccentricity
     &            sqrt(XN(NN)*XN(NN)+YN(NN)*YN(NN))/Weight, Weight*DX*DY
      End Do ! NN=1, 9
      Close(379)
      End Do ! r_power

!-----------------------------------------------------------------------

      ! Finally, calculate the "overlap" area
      XC = 0D0
      YC = 0D0
      TotalE = 0D0
      Do K = NZ0, NZ
      Do I = NXPhy0, NXPhy
      Do J = NYPhy0, NYPhy
        !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
        XX = I*DX
        YY = J*DY
        XC = XC + XX*Ed(I,J,K)
        YC = YC + YY*Ed(I,J,K)
        TotalE = TotalE + Ed(I,J,K)
      End Do
      End Do
      End Do
      XC = XC/TotalE
      YC = YC/TotalE

      ! Next, calculate <x^2> and <y^2>
      XN(2) = 0D0
      YN(2) = 0D0
      Weight = 0D0
      Do K = NZ0, NZ
      Do I = NXPhy0, NXPhy
      Do J = NYPhy0, NYPhy
        !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
        XX = I*DX - XC ! shift to the real center
        YY = J*DY - YC
        XN(2)=XN(2)+XX*XX*Ed(I,J,K)
        YN(2)=YN(2)+YY*YY*Ed(I,J,K)
        Weight=Weight+Ed(I,J,K)
      End Do
      End Do
      End Do

      XN(2) = XN(2)/Weight ! <x^2>
      YN(2) = YN(2)/Weight ! <y^2>

      ! Output overlap area
      Open(383, FILE='results/overlap.dat',STATUS='REPLACE')
      Write (383, '(E20.12)') CONSTPI*sqrt(XN(2)*YN(2))
      Close (383)

!----------- Initial profile related output finishes here --------------
!=======================================================================
!=======================================================================


       IW = 0
       do 9999 ITime = 1,MaxT
!***********************  Begin  Time Loop ********************************
        Print *,ITime ,' time= ', Time
CSHEN===========================================================================
C======Using a smaller time step for short initialization time \tau_0
            if (Time .lt. 0.59) then
                  DT = DT_2
            else
                  DT = DT_1
            endif
CSHEN====END====================================================================


!   ---Zhi-Changes---
        Call determineR0(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,Sd,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33)  !fermi-dirac function to regulate boundary (by determine a suitable R0)

        DO 3006 K = NZ0,NZ !Store last step U0 U1 U2
        DO 3006 J = NY0,NY
        DO 3006 I = NX0,NX
            PU0(I,J,K)   = U0(I,J,K)
            PU1(I,J,K)   = U1(I,J,K)
            PU2(I,J,K)   = U2(I,J,K)
            Temp0(I,J,K)   = Temp(I,J,K)
            etaTtp0(I,J,K) = etaTtp(I,J,K)
            XiTtP0(I,J,K)  = XiTtP(I,J,K)
3006    Continue

      ! Initialize PiXXRegulated and backup oldTTXX
      if (ViscousC>1D-6 .or. VisBulk > 1D-6) then
        
        !shear
        If (ViscousC>1D-6) Then
          Pi00Regulated = Pi00
          Pi01Regulated = Pi01
          Pi02Regulated = Pi02
          Pi11Regulated = Pi11
          Pi22Regulated = Pi22
          Pi12Regulated = Pi12
          Pi33Regulated = Pi33
          iRegulateCounter = 0D0
          If (echo_level>=3) Print*, "Regulation counts:", 
     &                               iRegulateCounter
        else
          Pi00 = 0D0 
          Pi01 = 0D0 
          Pi02 = 0D0 
          Pi11 = 0D0 
          Pi22 = 0D0 
          Pi12 = 0D0 
          Pi33 = 0D0 
        endif
        
        !bulk
        if(VisBulk > 1D-6) then
          PPIRegulated = PPI
          iRegulateCounterBulkPi = 0D0
          If (echo_level>=3) Print*, "Regulation Bulk Pi counts:", 
     &                               iRegulateCounterBulkPi
        else
          PPI = 0D0
        endif
        

      call dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Bd,Sd,Temp0,Temp,CMu, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,PNEW,NNEW)  !PNEW NNEW  related to root finding

        DIFFC = 0.125D0
        !DIFFC = 0D0
        if(ViscousC > 1D-6) then 
         call UPShasta2 (Pi01,Vx,Vy,PScT01, NX0,NX,NY0,NY,NZ0,NZ, 0,0,  !Pi01
     &            DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,-1,1, DIFFC)
         call UPShasta2 (Pi02,Vx,Vy,PScT02, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &            DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,-1, DIFFC)
         call UPShasta2 (Pi00,Vx,Vy,PScT00, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi33,Vx,Vy,PScT33, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi11,Vx,Vy,PScT11, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi12,Vx,Vy,PScT12, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
         call UPShasta2 (Pi22,Vx,Vy,PScT22, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
        endif
        
        if(VisBulk > 1D-6) then
          call UPShasta2 (PPI,Vx,Vy,PISc, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &             DT,DX,DY, NXPhy0,NYPhy0,  NXPhy,NYPhy,1,1, DIFFC)
        endif

!CHANGES
! Boundary regulation immediately
        DO K=NZ0,NZ
        DO J=NYPhy0,NYPhy
        DO I=NXPhy0,NXPhy
          xx=DX*I
          yy=DY*J
          rr=sqrt(xx**2+yy**2)
          ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)
          Pi00(I,J,K) = Pi00(I,J,K)*ff
          Pi01(I,J,K) = Pi01(I,J,K)*ff
          Pi02(I,J,K) = Pi02(I,J,K)*ff
          Pi33(I,J,K) = Pi33(I,J,K)*ff
          Pi11(I,J,K) = Pi11(I,J,K)*ff
          Pi12(I,J,K) = Pi12(I,J,K)*ff
          Pi22(I,J,K) = Pi22(I,J,K)*ff
          PPI(I,J,K) = PPI(I,J,K)*ff
        End Do
        End Do
        End Do
        
        if(ViscousC > 1D-6) then
         if (outputPiviolation) then
            call checkPiandoutputViolation(Time, DX, DY, Vx, Vy, Ed, PL,
     &     NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &     Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
         endif
         Do ! pi evolution
           call checkPiAll(iFailed, II, JJ, Time, Vx, Vy, Ed, PL,
     &     NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &     Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &     Pi00Regulated,Pi01Regulated,Pi02Regulated,Pi33Regulated,
     &     Pi11Regulated,Pi12Regulated,Pi22Regulated,echo_level)
         If (iFailed==0) Exit
         call regulatePi(iRegulateCounter,Time,NX0,NY0,NZ0,NX,NY,NZ,
     &       NXPhy0,NXPhy,NYPhy0,NYPhy,
     &       Ed,PL,PPI,
     &       Pi00,Pi01,Pi02,Pi11,
     &       Pi12,Pi22,Pi33,Vx,Vy,II,JJ)
         If (echo_level>=7) Then
           Print*, "After: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &         Pi00(II,JJ,1),Pi11(II,JJ,1),
     &         Pi22(II,JJ,1),Pi33(II,JJ,1),
     &         Pi01(II,JJ,1),Pi02(II,JJ,1),
     &         Pi12(II,JJ,1)
         End If
         iRegulateCounter = iRegulateCounter + 1D0
         If (echo_level>=3) Then
           Print*, "Regulation counts:", iRegulateCounter
         EndIf
         End Do ! pi evolution
        endif
      
        if(VisBulk > 1D-6) then 
          if(outputPiviolation) then
            call checkBulkPiandoutputViolation(Time, Dx, Dy, Ed, PL,
     &        NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &        PPI)
          endif
          Do ! BulkPi evolution
            call checkBulkPi(iFailed, II, JJ, Time, Ed, PL,
     &      NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &      PPI, echo_level)
          If (iFailed==0) Exit
          call regulateBulkPi(iRegulateCounterBulkPi,Time,NX0,NY0,NZ0,
     &        NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy,Ed,PL,PPI,II,JJ)
          If (echo_level>=7) Then
            Print*, "After: PPI=", PPI(II,JJ,1)
          End If
          iRegulateCounterBulkPi = iRegulateCounterBulkPi + 1D0
          If (echo_level>=3) Then
            Print*, "Regulation BulkPi counts:", iRegulateCounterBulkPi
          EndIf
          End Do ! BulkPi evolution
        endif

      else ! ideal hydro

        Pi00 = 0D0
        Pi01 = 0D0
        Pi02 = 0D0
        Pi11 = 0D0
        Pi22 = 0D0
        Pi12 = 0D0
        Pi33 = 0D0
        PPI = 0D0

      call dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Bd,Sd,Temp0,Temp,CMu, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,PNEW,NNEW)  !PNEW NNEW  related to root finding

      End If ! If (ViscousC>1D-6 .or. VisBulk > 1D-6) Then


      ! T(mu,nu) evolution
      DIFFC = 0.125
      !DIFFC = 0D0
      call UPShasta2 (TT01,Vx,Vy,ScT01, NX0,NX,NY0,NY,NZ0,NZ, 0,0,   !T00
     &       DT,DX,DY,NXPhy0-3,NYPhy0-3,NXPhy+3,NYPhy+3,-1,1,DIFFC)
      call UPShasta2 (TT02,Vx,Vy,ScT02, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &       DT,DX,DY,NXPhy0-3,NYPhy0-3,NXPhy+3,NYPhy+3,1,-1,DIFFC)
      call UPShasta2 (TT00,Vx,Vy,ScT00, NX0,NX,NY0,NY,NZ0,NZ, 0,0,
     &       DT,DX,DY,NXPhy0-3,NYPhy0-3,NXPhy+3,NYPhy+3,1,1,DIFFC)

      DO K=NZ0,NZ
      DO J=NYPhy0,NYPhy
      DO I=NXPhy0,NXPhy
        xx=DX*I
        yy=DY*J
        rr=sqrt(xx**2+yy**2)
        ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)
        TT00(I,J,K) = TT00(I,J,K)*ff
        TT01(I,J,K) = TT01(I,J,K)*ff
        TT02(I,J,K) = TT02(I,J,K)*ff
      End Do
      End Do
      End Do

      if(outputMovie) then
         DO K=NZ0,NZ
         DO J=NYPhy0,NYPhy
         DO I=NXPhy0,NXPhy
           if(mod(I, 5) .eq. 0 .and. mod(J, 5) .eq. 0) then  
             write(3773, '(4F18.8)')Time, I*Dx, J*Dy,
     &             Temp(I, J, K)*0.19733D0
           endif
         enddo
         enddo
         enddo
      endif
      
      if(IhydroJetoutput .eq. 1) then
!        output hydro infos
!        Units: [ed]=GeV/fm^3, [sd]=fm^-3, [p]=GeV/fm^3, [T]=GeV, [Vx]=[Vy]=1
!        Units: [Pi]=GeV/fm^3, [PPi]=GeV/fm^3
         Call writeHydroBlock(ITime-1, Ed*HbarC, Sd, PL*HbarC,
     &      Temp*HbarC,Vx, Vy, Pi00*HbarC, Pi01*HbarC, Pi02*HbarC, 
     &      Pi02*HbarC*0.0d0, Pi11*HbarC, Pi12*HbarC, Pi12*HbarC*0.0d0,
     &      Pi22*HbarC, Pi22*HbarC*0.0d0, Pi33*HbarC, PPI*HbarC)
!        output hydro infos, end
      endif

CSHEN===========================================================================
C====output the OSCAR body file from hydro evolution============ ===============
      if(IOSCAR) then
       if(ITime .lt. ItimeOSCAR) then
        call OSCARbodyoutput(IOSCARWrite, NX0, NX, NY0, NY, NZ0, NZ,
     &                       NXPhy0, NXPhy, NYPhy0, NYPhy, ITime,
     &                       Ed, PL, Temp, Vx, Vy,
     &                       Pi00, Pi01, Pi02, Pi11, Pi12,
     &                       Pi22, Pi33, PPi, VRelaxT, VRelaxT0)
       endif
      endif
CSHEN===end=====================================================================


      Hc = HbarC

      Do J=0,NXPhy,NXPhy+1
      Do I=NYPhy0,NYPhy,10
        write(93, '(5e15.5)')Time, I*DX, J*DY, Temp(I,J,NZ0)*HBarC,
     &                       Ed(I,J,NZ0)*Hc
      enddo
      enddo

      DO 203 I=0,NXPhy,20
         Write(*,'(500e12.3)')(Temp(I,J,NZ0)*Hc,J=0,NYPhy,20 )
203   continue
         Write(*,*)


      Print*, 'Center Energy Density',  Ed(0,0,NZ0)*HBarC
      Print*, 'Center T',   Temp(0,0,NZ0)*HBarC
      Print*, '   '
      Print*, '   '

C      write(*,'(f6.3, 12f10.5)') PL(10,10,NZ0), PL(30,30,NZ0),
C     &           PL(50,50,NZ0),  PL(70,70,NZ0),
C     &           PL(80,80,NZ0),  PL(90,90,NZ0),
C     &           PL(100,100,NZ0),  PL(110,110,NZ0), PL(120,120,NZ0)


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C~~~~     Freezeout Procedure (rewritten from Petor's code azhydro0p2)  A   ~~~~
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C      NDT = 5
      N   = ITime
      T   = Time

      IF (MOD(N+1,NDT).EQ.0) THEN
         X0 = 0.0
         Y0 = 0.0
         DO 5300 J = NY0,NY
         DO 5300 I = NX0,NX
            EPS1(I,J) = Ed(I,J,NZ0)*HbarC
            V11(I,J)  = Vx(I,J,NZ0)
            V21(I,J)  = Vy(I,J,NZ0)
            TEM1(I,J) = Temp(I,J,NZ0)*HbarC

            If (ViscousC>1D-6) Then
              FPi00(I,J) = Pi00(I,J,NZ0)
              FPi01(I,J) = Pi01(I,J,NZ0)
              FPi02(I,J) = Pi02(I,J,NZ0)
              FPi11(I,J) = Pi11(I,J,NZ0)
              FPi12(I,J) = Pi12(I,J,NZ0)
              FPi22(I,J) = Pi22(I,J,NZ0)
              FPi33(I,J) = Pi33(I,J,NZ0)
            End If
            if(VisBulk > 1D-6) then
              FPPI(I,J) = PPI(I,J,NZ0)
            endif

5300  CONTINUE

!      Call FreezeoutPro9 (EDEC,TFREEZ, TFLAG, IEOS,NDX,NDY,NDT,
!     &     EPS0,EPS1,TEM0,TEM1, V10,V20,V11,V21, EINS,NINT, IW,
!     &     F0Pi00,F0Pi01,F0Pi02,F0Pi33, F0Pi11,F0Pi12,F0Pi22,
!     &     FPi00,FPi01,FPi02,FPi33,FPi11,FPi12,FPi22, N,T,X0,Y0,
!     &     DX,DY,DT,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ)

      Call FreezeoutPro10 (EDEC, IEOS, NDX, NDY, NDT,
     &     EPS0,EPS1,V10,V20,V11,V21, NINT, IW,
     &     F0Pi00,F0Pi01,F0Pi02,F0Pi33,F0Pi11,F0Pi12,F0Pi22,
     &     FPi00,FPi01,FPi02,FPi33,FPi11,FPi12,FPi22, 
     &     F0PPI, FPPI,
     &     N,T,DX,DY,DT,NXPhy,NYPhy,NX0,NX,NY0,NY)

      print*,'after freezeout'

      DO 5400 J=NY0,NY
      DO 5400 I=NX0,NX
c                IF (TFLAG .EQ. 1) THEN
          TEM0(I,J) = TEM1(I,J)
c                END IF
          EPS0(I,J) = EPS1(I,J)
          V10(I,J)  = V11(I,J)
          V20(I,J)  = V21(I,J)
          TEM0(I,J) = TEM1(I,J)

          If (ViscousC>1D-6) Then
            F0Pi00(I,J) = Pi00(I,J,NZ0)
            F0Pi01(I,J) = Pi01(I,J,NZ0)
            F0Pi02(I,J) = Pi02(I,J,NZ0)
            F0Pi11(I,J) = Pi11(I,J,NZ0)
            F0Pi12(I,J) = Pi12(I,J,NZ0)
            F0Pi22(I,J) = Pi22(I,J,NZ0)
            F0Pi33(I,J) = Pi33(I,J,NZ0)
          End If
          if(VisBulk > 1D-6) then
            F0PPI(I,J) = PPI(I,J,NZ0)
          endif
5400  CONTINUE

      Print*, 'NINT', NINT

      IF (NINT.EQ.0) THEN
        WRITE(*,*) 'Decoupling all done at proper time,  T=',T
        goto 10000
C        IF (EARTERM.EQ.'EARLY') THEN
C           WRITE(*,*) '    Early termination of program.'
C           GOTO 10000
C        End if

      END IF

      END IF    !  IF (MOD(N+1,NDT).EQ.0) THEN

C~~~~~~~~~~~~~~Freezeout Procedure (rewritten from Peters code azhydro0p2)     ~~~~~~~~~~~~~~~~
5555  continue

      Time=Time+DT

9999  Continue  !***********************  End Time Loop  ******************************************
10000 Continue

CSHEN===========================================================================
C============OSCAR2008H output format for hydro movie and jet quenching=========
C=========fill up rest redundant step for the output file=======================
      if(IOSCAR) then
5786     if(ITime .lt. ItimeOSCAR) then
         ITime =  ITime + 1
         call OSCARredundantoutput(IOSCARWrite, NZ0, NZ, ITime,
     &                             NXPhy0, NXPhy, NYPhy0, NYPhy)
         write(*,*) ITime
         goto 5786
         endif
      endif

CSHEN===End=====================================================================

      Return
      End

C##################################################################################
      Subroutine FreezeoutPro9 (EDEC,TFREEZ, TFLAG,IEOS, NDX,NDY,NDT,
     &     EPS0,EPS1,TEM0,TEM1, V10,V20,V11,V21, EINS,NINT, IW,
     &     F0Pi00,F0Pi01,F0Pi02,F0Pi33, F0Pi11,F0Pi12,F0Pi22,
     &     FPi00,FPi01,FPi02,FPi33,FPi11,FPi12,FPi22, N,T,X0,Y0,
     &     DX,DY,DT,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ)
*    a subroutine to calculate the freezeout surface which is changed from Peter  Azhydro0P2
*     TFLAG=0:      TFLAG=0 constant energy        (I)
*     EDEC: decoupling energy (TFLAG=0)            (I)
*     T=Time   N,Timestep for the largest Loop.     X0=0.0 Y0=0.0
      Implicit Double Precision (A-H, O-Z)

      DIMENSION EPS0(NX0:NX,NY0:NY),EPS1(NX0:NX,NY0:NY) ! Energy density in previous and current step
      DIMENSION TEM0(NX0:NX,NY0:NY),TEM1(NX0:NX,NY0:NY) ! Temperature density in previous and current step
      DIMENSION V10(NX0:NX,NY0:NY),V20(NX0:NX,NY0:NY)   !velocity in X Y in Previous step
      DIMENSION V11(NX0:NX,NY0:NY),V21(NX0:NX,NY0:NY)   !velocity in X Y in current step

      DIMENSION F0Pi00(NX0:NX,NY0:NY),FPi00(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi01(NX0:NX,NY0:NY),FPi01(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi02(NX0:NX,NY0:NY),FPi02(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi11(NX0:NX,NY0:NY),FPi11(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi12(NX0:NX,NY0:NY),FPi12(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi22(NX0:NX,NY0:NY),FPi22(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step
      DIMENSION F0Pi33(NX0:NX,NY0:NY),FPi33(NX0:NX,NY0:NY)   !Stress Tensor in previous and current step

      DIMENSION EPCUB(8),SURFS(12,0:2),VMID(0:2)
      Double Precision EpcubP(8), SSurfs(12,0:2), SVmid(0:2)! permuted EPCUB, to be averaged to get more accurate results for, e.g, VMID; summed unsigned Surfs and Vmid
      Double Precision signSurfs(12,0:2), signVmid(0:2) ! store signs, used in the averaging process
      DIMENSION IBIT(12,12)
C      Parameter(NDX=2,NDY=2, NDT=5)!
C      Common /Nfreeze/ NDX, NDY, NDT

      Parameter(MAXERR=150 )!Max Number erro in Cubcal
      PARAMETER (PI=3.141592653d0, HBARC=.19733d0)

      Logical INTERSECT, WRITING, WRITINGX, WRITINGY
      INTEGER TFLAG, EINS
      INTEGER NXPHY,NYPHY

!** Zhi ***
      Integer absI, absJ ! abs(I) and abs(J) used in the loop
      Integer tmpI

* ----------------------------------------------------------------------
*   ** A 'BITCHART' USED BY SUBROUINE 'CUBCAL' TO FIND AN APPROXIMATION
*   ** FOR A SEGMENT OF THE DECOUPLING SURFACE INTERSECTING A
*   ** ((dt)x(dx)x(dy)) 'CUBE'.
      DATA (IBIT(1,I),I=1,12)/  0,1,2,1,1,1,0,0,2,0,0,0/
      DATA (IBIT(2,I),I=1,12)/  1,0,1,2,0,1,1,0,0,2,0,0/
      DATA (IBIT(3,I),I=1,12)/  2,1,0,1,0,0,1,1,0,0,2,0/
      DATA (IBIT(4,I),I=1,12)/  1,2,1,0,1,0,0,1,0,0,0,2/
      DATA (IBIT(5,I),I=1,12)/  1,0,0,1,0,2,0,2,1,0,0,1/
      DATA (IBIT(6,I),I=1,12)/  1,1,0,0,2,0,2,0,1,1,0,0/
      DATA (IBIT(7,I),I=1,12)/  0,1,1,0,0,2,0,2,0,1,1,0/
      DATA (IBIT(8,I),I=1,12)/  0,0,1,1,2,0,2,0,0,0,1,1/
      DATA (IBIT(9,I),I=1,12)/  2,0,0,0,1,1,0,0,0,1,2,1/
      DATA (IBIT(10,I),I=1,12)/ 0,2,0,0,0,1,1,0,1,0,1,2/
      DATA (IBIT(11,I),I=1,12)/ 0,0,2,0,0,0,1,1,2,1,0,1/
      DATA (IBIT(12,I),I=1,12)/ 0,0,0,2,1,0,0,1,1,2,1,0/

      COMMON /TAREOS/ TBNP01,TEPP01,TDBNP1,TDEPP1,NTBNP1,NTEPP1,
     &                TBNP02,TEPP02,TDBNP2,TDEPP2,NTBNP2,NTEPP2,
     &                TMAT1(0:300,0:250),TMAT2(0:300,0:250)

      Integer IEOS2dec
      Common/IEOS2dec/ IEOS2dec  ! IEOS=2 decouple by gluon/pion

c      IF (MOD(N+1,NDT).EQ.0) THEN   !** N+1 *******************************************************************
       Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling redious
       Common/R0Aeps/ R0,Aeps

CSHEN======================================================================
C=============add chemical potential for EOS-PCE
      Integer, Parameter :: RegEOSMudatasize = EOSMUDATALENGTH   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = EOSMUMAXPARTICLE        !maximum allowed stable particles for partially chemical equilibrium EOS
      Integer :: Inumparticle       !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne            !total rows in mu table
      double precision :: EOS_Mu_e0   !lowest energy density in mu table
      double precision :: EOS_Mu_de   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      Integer :: IMuflag

      double precision :: XMufreeze(1:IMax_Mu)
      double precision :: XMufreezetemp

      common /EOSMudata/MuEOSdata, IMuflag
      common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne, 
     &                            Inumparticle
CSHEN======================================================================


      Rdec1=DX*NXPhy*1.414
      Rdec2=0.0

      DTD=NDT*DT

      NYFUL=NYPhy+2
      NXFUl=NXPhy+2

      NINT=0
      NINT2=0
      DO 1234 J=-NYFUL,NYFUL
        Y=Y0+J*DY
      DO 1235 I=-NXFUL,NXFUL
        IF (IEOS .EQ. 2) THEN
          if(IEOS2dec.eq.0) then ! decouple to gluons
            TEM1(I,J) =
     &          ((EPS1(I,J)-0.0)*120./(PI**2.)/169. *(HBARC**3.))**(.25) !changed by Evan  gluon
          else                   !decouple to poins
            TEM1(I,J)=SQRT(HBARC*SQRT(10.0*HBARC*EPS1(I,J))/PI)  !Ideal pion gas temperature
          end if
        ELSE IF (IEOS.EQ.0) THEN  !EOS=0 with phase transition
             BBBNN = 0.0  !  BNR1(I,J)
             EEENN = EPS1(I,J)
             IF (EEENN.GT.24) THEN
             TEM1(I,J) =((EEENN-.3642)*120./(PI**2.) /169.
     &              *(HBARC**3.))**(.25)
             ELSE
             TEM1(I,J) = EOUT(BBBNN,EEENN,0.0d0,
     &              TBNP01,TEPP01,TDBNP1,TDEPP1,NTBNP1,NTEPP1,
     &              TBNP02,TEPP02,TDBNP2,TDEPP2,NTBNP2,NTEPP2,
     &              TMAT1,TMAT2)
             END IF
        ELSE IF (IEOS.EQ.99) THEN
        ELSE
        END IF
 1235 CONTINUE
 1234 CONTINUE



      IF (TFLAG .EQ. 1) THEN
       FREEZ = TFREEZ
      ELSE
       FREEZ = EDEC
      END IF


      DO 500 absJ=0,NYFUL-NDY,NDY
      DO 510 absI=0,NXFUL-NDX,NDX

      Do 1911 tmpI=1,4 ! change this to tmpI=1,1 will output freeze-out data only in the 1st quadrant

      If (tmpI == 1) Then
        I = absI
        J = absJ
      ElseIf (tmpI == 2) Then
        I = -absI - NDX
        J = absJ
      ElseIf (tmpI == 3) Then
        I = -absI - NDX
        J = -absJ - NDY
      ElseIf (tmpI == 4) Then
        I = absI
        J = -absJ - NDY
      Else
        Print *, "You kidding. tmpI=", tmpI
        Stop
      End If

      Y = J*DY
      X = I*DX

      DXD = NDX*DX
      DYD = NDY*DY

      ! The cubic returned by SECTIO3 have the same orientation in all quadrants, determined by the order of x and y values of the vertices.
      IF (TFLAG .EQ. 1) THEN !freeze out by constant T
        CALL SECTIO2(FREEZ,INTERSECT,EPCUB,TEM0,TEM1,
     &        I,J,NDX,NDY,NX0,NY0,NX,NY)
      ELSE
        CALL SECTIO2(FREEZ,INTERSECT,EPCUB,EPS0,EPS1,
     &        I,J,NDX,NDY,NX0,NY0,NX,NY)
      ENDIF

      IF (INTERSECT) THEN   !***(INTERSECT)

      NINT    = NINT+1
      NINT2   = NINT2+1
      WRITING = .TRUE.
      if(absJ.eq.0) then
        WRITINGX = .TRUE.     !output freeze-out surface along x-axis
      else
        WRITINGX = .FALSE.
      endif
      if(absI.eq.0) then
        WRITINGY = .TRUE.     !output freeze-out surface along y-axis
      else
        WRITINGY = .FALSE.
      endif


      CALL CUBCAL(FREEZ,EPCUB,SURFS,NSURFS,VMID,IBIT,
     &             DTD,DXD,DYD,NERR,MAXERR)

      DAH0 = SURFS(1,0)
      DAH1 = SURFS(1,1)
      DAH2 = SURFS(1,2)

      DO 520 NS = 2,NSURFS
        DAH0 = DAH0+SURFS(NS,0)
        DAH1 = DAH1+SURFS(NS,1)
        DAH2 = DAH2+SURFS(NS,2)
520   CONTINUE

       TM = T-DTD+DT+VMID(0)
       XM = X+VMID(1)
       YM = Y+VMID(2)

      CALL P4(I,J,NDX,NDY,NDT,VMID,V10,V11,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,V1,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,V20,V21,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,V2,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi00,FPi00,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi00,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi01,FPi01,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi01,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi02,FPi02,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi02,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi11,FPi11,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi11,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi12,FPi12,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi12,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi22,FPi22,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi22,0)
      CALL P4(I,J,NDX,NDY,NDT,VMID,F0Pi33,FPi33,
     &        NX0,NY0,NX,NY,DTD,DXD,DYD,CPi33,0)


       BN   = 0.0
       PDec = PEPS(0.0d0,Edec)

       IF (TFLAG .EQ. 1) THEN   !freeze out by constant T
          CALL P4(I,J,NDX,NDY,NDT,VMID,EPS0,EPS1,
     &           NX0,NY0,NX,NY,DTD,DXD,DYD,EDEC,0)
          TDEC = FREEZ
       ELSE                    !freeze out by constant E
         IF (IEOS .EQ. 2) THEN  !freeze out by const.energy density EDec then find out TDec
             !Print *, "IEOS2dec=", IEOS2dec
             if(IEOS2dec.eq.0) then ! decouple to gluons
              TDEc=((EDEC-0.0)*120./(PI**2.)/169.
     &              *(HBARC**3.))**(.25)  !changed by Evan  !gluon
             else  ! decouple to pions
              TDEC=SQRT(HBARC*SQRT(10.0*HBARC*EDEC)/PI)
             end if
         ELSE IF (IEOS.EQ.0) THEN
               BN   = 0.0                 !FREEZ=EDec
               TDEC = EOUT(BN,FREEZ,0.0d0,
     &                     TBNP01,TEPP01,TDBNP1,TDEPP1,NTBNP1,NTEPP1,
     &                     TBNP02,TEPP02,TDBNP2,TDEPP2,NTBNP2,NTEPP2,
     &                     TMAT1,TMAT2)
         Else if (IEOS.eq.4) Then
               ee   = EDEC/HbarC          !chenge to fm
               TDEC = TEIEOS4(ee)*HBarC   !change to GEV
         ELSE IF (IEOS.eq.5) then         !CSHEN SM-EOS Q
               ee   = EDEC/HbarC          !1/fm^4
               TDEC = TEIEOS5pp(ee)*HbarC !GeV
         ELSE IF (IEOS.EQ.6) THEN         !CSHEN New EOS
               ee   = EDEC                !GeV/fm^3
               TDEC = TEOSL6(ee)          !GeV
         else if (IEOS.eq.7) then         !CSHEN EOS from table
               ee   = EDEC                !GeV/fm^3
               TDEC = TEOSL7(ee)          !GeV
         ELSE IF (IEOS.EQ.10) THEN
         ELSE IF (IEOS.EQ.99) THEN
         ELSE
         END IF
       END IF


c                 BAMU = EOUT(BN,EDEC,0.0,


       BAMU = 0.0
       SMU  = 0.0

       DA0  = DAH0*1.
       DA1  = DAH1*1.
       DA2  = DAH2*1.

       VZCM = V1
       VRCM = V2

       IF (WRITING) THEN   !---
         WRITE(98,2553) T-DTD+DT, TM, XM,YM,
     &                  Sqrt(XM**2+YM**2), R0,Aeps    !surface

         IF (WRITINGX) THEN
            write(82,2553) TM, Sqrt(XM**2+YM**2),
     &                   Sqrt(VZCM**2+VRCM**2),XM,VZCM,YM,VRCM
         end if
         if (WRITINGY) THEN
            write(83,2553) TM, Sqrt(XM**2+YM**2),
     &                   Sqrt(VZCM**2+VRCM**2),XM,VZCM,YM,VRCM
         end if
         IW = IW+1   !count the number in surface written file
         if(Sqrt(XM**2+YM**2).le.Rdec1) Rdec1=Sqrt(XM**2+YM**2)
         if(Sqrt(XM**2+YM**2).ge.Rdec1) Rdec2=Sqrt(XM**2+YM**2)

         IF (IEOS.EQ.99) THEN
         ELSE
C--------- ---with the same physical unit as energt dencity              EPS1(I,J)=Ed(I,J,NZ0)*HbarC
           WRITE(99,2555) TM,DA0,DA1,DA2,VZCM,VRCM,EDEC,BN, !
     &                  TDEC,BAMU,SMU, PDec, CPi33*HbarC,
     &                  CPi00*HbarC,CPi01*HbarC,CPi02*HbarC,
     &                  CPi11*HbarC,CPi12*HbarC,CPi22*HbarC
         END IF

CSHEN=====================================================================
C===============output chemcial potential for EOS-PCE=====================
         IF (((IEOS.EQ.6).or.(IEOS.eq.7)).and.(IMuflag.eq.0)) then
           if(Inumparticle.ne.0) then
            do L=1,Inumparticle
               call interpCubic(MuEOSdata(:,L), RegEOSMudatasize, 
     &                  EOS_Mu_e0, EOS_Mu_de, EDEC, XMufreezetemp)
               XMufreeze(L) = XMufreezetemp
            end do
            write(81,'(E15.6)', Advance='NO') EDEC
            do L=1,Inumparticle
              write(81,'(E15.6)', Advance='NO') XMufreeze(L)
            enddo
            write(81,*)
            IMuflag = 1
           endif
         endif
CSHEN==========end========================================================

         ep = Edec+Pdec
C         Write(79,2555)TM,XM,YM,V1,V2, CPi11/ep, CPi12/ep,CPi22/ep,
C     &                 CPi33/ep, CPi00/ep, CPi01/ep, CPi02/ep
C         Write(78,2555)TM,Sqrt(XM**2+YM**2),Sqrt(V1**2+V2**2),
C     &                 (CPi11+CPi22)/ep, CPi33/ep, CPi00/ep

       ELSE               !---
         NINT = NINT - 1
       END IF             !---

       END IF      !***(INTERSECT)

1911    Continue
 510    CONTINUE
 500    CONTINUE



c          END IF    !   IF (MOD(N+1,NDT).EQ.0) THEN   !** N+1 *******************************************************************
C          EINS=EINS+1



2553   FORMAT(8F20.8)
2555   FORMAT(19E20.8E3)
2565   FORMAT(14E14.6)


      Return
      End

C##################################################################################
      Subroutine FreezeoutPro10(Edec, IEOS, NDX, NDY, NDT,
     &     EPS0,EPS1,V10,V20,V11,V21, NINT, IW,
     &     F0Pi00,F0Pi01,F0Pi02,F0Pi33,F0Pi11,F0Pi12,F0Pi22,
     &     FPi00,FPi01,FPi02,FPi33,FPi11,FPi12,FPi22, 
     &     F0PPI, FPPI,
     &     N,T,DX,DY,DT,NXPhy,NYPhy,NX0,NX,NY0,NY)
*     a subroutine to calculate the freeze-out surface using cornelius from P. Huovinen
*     T=Time   N, Timestep for the largest Loop.
      Implicit none

      double precision :: Edec
      integer :: IEOS, NDX, NDY, NDT, NINT, IW
      integer :: N
      double precision :: T, X, Y
      double precision :: DX, DY, DT
      integer :: NXPHY, NYPHY, NX0, NX, NY0, NY

      double precision, Dimension(NX0:NX,NY0:NY) :: EPS0
      double precision, Dimension(NX0:NX,NY0:NY) :: EPS1 ! Energy density in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: V10
      double precision, Dimension(NX0:NX,NY0:NY) :: V20  !velocity in X Y in Previous step
      double precision, Dimension(NX0:NX,NY0:NY) :: V11
      double precision, Dimension(NX0:NX,NY0:NY) :: V21  !velocity in X Y in current step

      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi00
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi00   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi01
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi01   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi02
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi02   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi11
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi11   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi12
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi12   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi22
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi22   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0Pi33
      double precision, Dimension(NX0:NX,NY0:NY) :: FPi33   !Stress Tensor in previous and current step
      double precision, Dimension(NX0:NX,NY0:NY) :: F0PPI
      double precision, Dimension(NX0:NX,NY0:NY) :: FPPI   !Bulk Stress Tensor in previous and current step

      double precision, Dimension(0:1,0:1,0:1) :: Cube
      double precision, Dimension(0:2,4) :: dSigma
      double precision, Dimension(0:2,4) :: Vmid
      double precision, Dimension(0:2) :: Vmidpoint
      integer :: Nsurf, Nambi, Ndisc
      integer :: iSurf

      double precision, PARAMETER :: pi=3.141592653d0, HbarC=.19733d0

      Logical :: Intersect, WRITING, WRITINGX, WRITINGY
      double precision :: DTFreeze, DXFreeze, DYFreeze
      integer :: NXFUL, NYFUL
      double precision :: Tmid, Xmid, Ymid
      double precision :: v1mid, v2mid
      double precision :: CPi00, CPi01, CPi02
      double precision :: CPi11, CPi12, CPi22, CPi33
      double precision :: CPPI
      double precision :: BN, Pdec, ee, Tdec
      double precision :: BAMU, SMU
      double precision :: DA0, DA1, DA2
      double precision :: PEPS, TEIEOS4, TEIEOS5pp, TEOSL6, TEOSL7

!** Zhi ***
      Integer :: absI, absJ ! abs(I) and abs(J) used in the loop
      Integer :: I, J, L
      Integer :: tmpI

      double precision :: R0, Aeps
      Common/R0Aeps/ R0,Aeps

CSHEN======================================================================
C=============add chemical potential for EOS-PCE
      Integer, Parameter :: RegEOSMudatasize = EOSMUDATALENGTH   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = EOSMUMAXPARTICLE        !maximum allowed stable particles for partically chemical equilibrium EOS
      Integer :: Inumparticle       !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne            !total rows in mu table
      double precision :: EOS_Mu_e0   !lowest energy density in mu table
      double precision :: EOS_Mu_de   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      Integer :: IMuflag

      double precision :: XMufreeze(1:IMax_Mu)
      double precision :: XMufreezetemp

      common /EOSMudata/MuEOSdata, IMuflag
      common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne, 
     &                            Inumparticle
CSHEN======================================================================

      DTFreeze = NDT*DT
      DXFreeze = NDX*DX
      DYFreeze = NDY*DY

      NYFUL = NYPhy + 2
      NXFUl = NXPhy + 2

      NINT = 0

      DO 500 absJ=0,NYFUL-NDY,NDY
      DO 510 absI=0,NXFUL-NDX,NDX
      Do 1911 tmpI=1,4 ! change this to tmpI=1,1 will output freeze-out data only in the 1st quadrant
       If (tmpI == 1) Then
         I = absI
         J = absJ
       ElseIf (tmpI == 2) Then
         I = -absI - NDX
         J = absJ
       ElseIf (tmpI == 3) Then
         I = -absI - NDX
         J = -absJ - NDY
       ElseIf (tmpI == 4) Then
         I = absI
         J = -absJ - NDY
       Else
         Print *, "You kidding. tmpI=", tmpI
         Stop
       End If
       Y = J*DY
       X = I*DX

       CALL SECTIONCornelius(Edec, Intersect, Cube, EPS0, EPS1,
     &                       I,J,NDX,NDY,NX0,NY0,NX,NY)

       IF (INTERSECT) THEN   !***(INTERSECT)
         NINT = NINT+1
         WRITING = .TRUE.
         if(absJ.eq.0) then
           WRITINGX = .TRUE.     !output freeze-out surface along x-axis
         else
           WRITINGX = .FALSE.
         endif
         if(absI.eq.0) then
           WRITINGY = .TRUE.     !output freeze-out surface along y-axis
         else
           WRITINGY = .FALSE.
         endif
         dSigma = 0.0d0
         Nsurf = 0
         Vmid = 0.0d0
         Nambi = 0
         Ndisc = 0
         CALL Cornelius2(Edec, Cube, dSigma, Nsurf, Vmid,
     &                   DTFreeze, DXFreeze, DYFreeze, Nambi, Ndisc)
         Do iSurf = 1, Nsurf  ! loop over the freeze out surface in the cube
           Tmid = T - DTFreeze + DT + Vmid(0,iSurf)
           Xmid = X + Vmid(1,iSurf)
           Ymid = Y + Vmid(2,iSurf)
           Vmidpoint = Vmid(:,iSurf)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,V10,V11,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,v1mid,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,V20,V21,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,v2mid,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi00,FPi00,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi00,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi01,FPi01,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi01,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi02,FPi02,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi02,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi11,FPi11,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi11,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi12,FPi12,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi12,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi22,FPi22,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi22,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0Pi33,FPi33,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPi33,0)
           CALL P4(I,J,NDX,NDY,NDT,Vmidpoint,F0PPI,FPPI,
     &             NX0,NY0,NX,NY,DTFreeze,DXFreeze,DYFreeze,CPPI,0)
           BN   = 0.0
           Pdec = PEPS(0.0d0,Edec)

           if (IEOS.eq.4) Then
              ee   = Edec/HbarC          !change to fm
              Tdec = TEIEOS4(ee)*HbarC   !change to GEV
           ELSE IF (IEOS.eq.5) then         !CSHEN SM-EOS Q
              ee   = Edec/HbarC          !1/fm^4
              Tdec = TEIEOS5pp(ee)*HbarC !GeV
           ELSE IF (IEOS.EQ.6) THEN         !CSHEN New EOS
              ee   = Edec                !GeV/fm^3
              Tdec = TEOSL6(ee)          !GeV
           else if (IEOS.eq.7) then         !CSHEN EOS from table
              ee   = Edec                !GeV/fm^3
              Tdec = TEOSL7(ee)          !GeV
           ELSE
           END IF
           BAMU = 0.0
           SMU  = 0.0

           DA0  = dSigma(0, iSurf)
           DA1  = dSigma(1, iSurf)
           DA2  = dSigma(2, iSurf)

           IF (WRITING) THEN   !---
             ! write to surface.dat
             WRITE(98,2553) T-DTFreeze+DT, Tmid, Xmid, Ymid,
     &                      Sqrt(Xmid**2+Ymid**2), R0,Aeps
             ! write to decdat2.dat
             WRITE(99,2555) Tmid, DA0, DA1, DA2, v1mid, v2mid, Edec, BN,
     &                      Tdec, BAMU, SMU, Pdec, CPi33*HbarC,
     &                      CPi00*HbarC,CPi01*HbarC,CPi02*HbarC,
     &                      CPi11*HbarC,CPi12*HbarC,CPi22*HbarC,
     &                      CPPI*HbarC
             IF (WRITINGX) THEN  !write to SurfaceX.dat
                write(82,2553) Tmid, Sqrt(Xmid**2+Ymid**2),
     &                Sqrt(v1mid**2+v2mid**2), Xmid, v1mid, Ymid, v2mid
             end if
             if (WRITINGY) THEN  !write to SurfaceY.dat
                write(83,2553) Tmid, Sqrt(Xmid**2+Ymid**2),
     &                Sqrt(v1mid**2+v2mid**2), Xmid, v1mid, Ymid, v2mid
             end if
             IW = IW+1   !count the number in surface written file
             !output chemical potential for EOS-PCE
             IF (((IEOS.EQ.6).or.(IEOS.eq.7)).and.(IMuflag.eq.0)) then
               if(Inumparticle.ne.0) then
                do L=1,Inumparticle
                   call interpCubic(MuEOSdata(:,L), RegEOSMudatasize, 
     &                      EOS_Mu_e0, EOS_Mu_de, Edec, XMufreezetemp)
                   XMufreeze(L) = XMufreezetemp
                end do
                write(81,'(E15.6)', Advance='NO') Edec
                do L=1,Inumparticle
                  write(81,'(E15.6)', Advance='NO') XMufreeze(L)
                enddo
                write(81,*)
                IMuflag = 1
               endif
             endif
           ELSE               !write to file
             NINT = NINT - 1
           END IF             !write to file
         Enddo  !loop over Nsurf
       END IF      !***(INTERSECT)

1911    Continue
 510    CONTINUE
 500    CONTINUE

2553   FORMAT(8F20.8)
2555   FORMAT(20E20.8E3)
2565   FORMAT(14E14.6)

      Return
      End



C####################################################################
      Subroutine TransportPi6( Pi00,Pi01,Pi02,Pi33, Pi11,Pi12,Pi22,
     &  PPI,Ed, Sd, PL, Temp, Temp0, U0,U1,U2,PU0,PU1,PU2, DX,DY,DZ,DT,
     &  NX0,NY0,NZ0, NX,NY,NZ, Time, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &  VRelaxT, VRelaxT0)
C------- Transport Pi00,Pi01,Pi02,Pi03,Pi11, Pi12,Pi22 by first order theory
        Implicit Double Precision (A-H, O-Z)
        Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ) !used for checkPi
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ) !used for checkPi

        Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
        Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

       Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ) ! Stress Tensor
       Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)

       Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ) ! Stress Tensor
       Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)

       Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension DPc00(NX0:NX, NY0:NY, NZ0:NZ) ! Differential part of Pi source term
        Dimension DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DPc11(NX0:NX, NY0:NY, NZ0:NZ) ! Differential part of Pi source term
        Dimension DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DDU0(NX0:NX, NY0:NY, NZ0:NZ) ! Differential part of Pi source term
        Dimension DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU2(NX0:NX, NY0:NY, NZ0:NZ) !

       Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
       Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure density
       Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
       Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

       Dimension CMu(NX0:NX, NY0:NY, NZ0:NZ) !Local chemical potential

       Dimension VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient
       Dimension VCBeta(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient Beta2
       Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
       Dimension SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \sita
       Dimension DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms

        Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008
        Dimension etaTtp0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

       Dimension VBulk(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient-bulk vicosity \xi
       Dimension VCBeta0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient Beta0
       Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time \tau_PI
       Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008


       COMMON /EOSSEL/ IEOS !Type of EOS
       COMMON /IEin/ IEin  !parameter for initializtion by entropy 1 or energy 0
         Parameter(XC0=1.0d-5, AAC=1.0d-6)  !XC0 charge of root acuracy  !AAC discard overflow value


!   ---Changes-Zhi--- ***
      Double Precision R0,Aeps
      Common/R0Aeps/ R0,Aeps

      Vx = U1/U0
      Vy = U2/U0

       call getInitialR0(PU0,PU1,PU2,PU3,U0,U1,U2,U3,DX,DY,DZ,DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT, Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ, Ed,Sd,PL,VCoefi)
       Write (*,*) "After R0 initialization: R0=",R0
!   ---Zhi-End---

        call ViscousCoefi8(Ed,Sd,PL,Temp,
     &  VCoefi,VCBeta,VRelaxT,etaTtp, VBulk, VCBeta0,VRelaxT0,XiTtP,
     &  Time,DX,DY,DT,NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

         do 10 k=NZ0,NZ
         do 10 j=NY0,NY
         do 10 i=NX0,NX
              Temp0(i,j,k)=Temp(i,j,k)
10       continue

        call PiS4U5(PU0,PU1,PU2,PU3,U0,U1,U2,U3, DX,DY,DZ, DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT,  Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ)
C
      do 30 k=NZ0,NZ
      do 30 j=NYPhy0,NYPhy
      do 30 i=NXPhy0,NXPhy

        ConsP=VCoefi(i,j,k)*2


       Pi01(i,j,k)=ConsP*DPc01(i,j,k)
       Pi02(i,j,k)=ConsP*DPc02(i,j,k)
       Pi33(i,j,k)=ConsP*DPc33(i,j,k)

       !Pi00(i,j,k)=ConsP*DPc00(i,j,k)
       Pi11(i,j,k)=ConsP*DPc11(i,j,k)
       Pi12(i,j,k)=ConsP*DPc12(i,j,k)
       Pi22(i,j,k)=ConsP*DPc22(i,j,k)

       Pi00(i,j,k)=(Pi01(i,j,k)*U1(I,j,k)+Pi02(i,j,k)*U2(I,j,k))
     &      /DMax1(1.0*d-34,U0(i,j,k))

        PPI(I,J,K)=(-1.0)*VBulk(i,j,k)*SiLoc(i,j,k)
30    continue

       call TriSembdary3(Pi00,Pi01, Pi02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(Pi11,Pi22, Pi33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(Pi12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)


       Return
       End





C####################555555555555555555555555555555555555#########################
C--------------------------------------------------------------------------------
       Subroutine PiS4U5(PU0,PU1,PU2,PU3,U0,U1,U2,U3, DX,DY,DZ, DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT,  Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ)
        Implicit Double Precision (A-H, O-Z)

        Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
        Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension DPc00(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DPc11(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DDU0(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU2(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature  in last time step
        Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
        Dimension SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \sita
        Dimension DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms

        Common/R0Aeps/ R0,Aeps
        Common /Accu/Accu
        Common /R0Bdry/ R0Bdry

C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C       Dimension A00(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A01(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A02(NX0:NX, NY0:NY, NZ0:NZ)

C       Dimension A10(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A11(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A12(NX0:NX, NY0:NY, NZ0:NZ)

C       Dimension A20(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A21(NX0:NX, NY0:NY, NZ0:NZ)
C       Dimension A22(NX0:NX, NY0:NY, NZ0:NZ)
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c       goto 199

        DO 100 K=NZ0,NZ
        DO 100 J=NYPhy0,NYPhy
        DO 100 I=NXPhy0,NXPhy

          xx=DX*I
          yy=DY*J
          rr=sqrt(xx**2+yy**2)
          ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)


         D0U0=(U0(I,J,K)-PU0(I,J,K))/DT
         D0U1=(U1(I,J,K)-PU1(I,J,K))/DT
         D0U2=(U2(I,J,K)-PU2(I,J,K))/DT

       If(abs(Accu-3.0).le.0.00001) then  !3pt formula
         D1U0=(U0(I+1,J,K)-U0(I-1,J,K))/(2.0*DX)
         D1U1=(U1(I+1,J,K)-U1(I-1,J,K))/(2.0*DX)
         D1U2=(U2(I+1,J,K)-U2(I-1,J,K))/(2.0*DX)
         D2U0=(U0(I,J+1,K)-U0(I,J-1,K))/(2.0*DY)
         D2U1=(U1(I,J+1,K)-U1(I,J-1,K))/(2.0*DY)
         D2U2=(U2(I,J+1,K)-U2(I,J-1,K))/(2.0*DY)
       else if (abs(Accu-5.0).le.0.00001) then !5pt formula
         D1U0=(U0(I+1,J,K)*2.0d0/3.0d0-U0(I-1,J,K)*2.0d0/3.0d0
     &            -U0(I+2,J,K)/12.0d0+U0(I-2,J,K)/12.0d0)/DX
         D1U1=(U1(I+1,J,K)*2.0d0/3.0d0-U1(I-1,J,K)*2.0d0/3.0d0
     &            -U1(I+2,J,K)/12.0d0+U1(I-2,J,K)/12.0d0)/DX
         D1U2=(U2(I+1,J,K)*2.0d0/3.0d0-U2(I-1,J,K)*2.0d0/3.0d0
     &            -U2(I+2,J,K)/12.0d0+U2(I-2,J,K)/12.0d0)/DX

         D2U0=(U0(I,J+1,K)*2.0d0/3.0d0-U0(I,J-1,K)*2.0d0/3.0d0
     &           -U0(I,J+2,K)/12.0d0+U0(I,J-2,K)/12.0d0)/DY
         D2U1=(U1(I,J+1,K)*2.0d0/3.0d0-U1(I,J-1,K)*2.0d0/3.0d0
     &           -U1(I,J+2,K)/12.0d0+U1(I,J-2,K)/12.0d0)/DY
         D2U2=(U2(I,J+1,K)*2.0d0/3.0d0-U2(I,J-1,K)*2.0d0/3.0d0
     &           -U2(I,J+2,K)/12.0d0+U2(I,J-2,K)/12.0d0)/DY
        else
           Print*, "wrong input for Accu,
     &         Accu=3or5 for 3pt or 5pt cal of deriv. "
        end if

         CS=(D0U0+D1U1+D2U2+U0(I,J,K)/Time)/3.0

      If (Time > 0.8) Then
      If (CS .ne. CS) Then
        Print *, "Invalid CS!"
        Print *, "i,j,k=", i,j,k
        Print *, "D0U0=", D0U0, "D1U1=", D1U1, "D2U2=", D2U2
        Print *, "D2U1=", D2U1, "D2U0=", D2U0
        Print *, "U0=", U0(I,J,K), "Time=",Time
        Print *, "U0(I,J+1,K), U0(I,J-1,K)=", U0(I,J+1,K), U0(I,J-1,K)
        Print *, "U1(I,J+1,K), U1(I,J-1,K)=", U1(I,J+1,K), U1(I,J-1,K)
        Print *, "U2(I,J+1,K), U2(I,J-1,K)=", U2(I,J+1,K), U2(I,J-1,K)
      EndIf
      EndIf

      DU0=U0(I,J,K)*D0U0+U1(I,J,K)*D1U0+U2(I,J,K)*D2U0
      DU1=U0(I,J,K)*D0U1+U1(I,J,K)*D1U1+U2(I,J,K)*D2U1
      DU2=U0(I,J,K)*D0U2+U1(I,J,K)*D1U2+U2(I,J,K)*D2U2

      DDU0(i,j,k)=DU0
      DDU1(i,j,k)=DU1
      DDU2(i,j,k)=DU2

       DPc00(I,J,K)=D0U0-U0(I,J,K)*DU0+CS*(U0(I,J,K)**2-1.0)
       DPc01(I,J,K)=0.5*(D0U1-D1U0)-0.5*(U1(I,J,K)*DU0+U0(I,J,K)*DU1)
     &                      +CS*(U1(I,J,K)*U0(I,J,K))
       DPc02(I,J,K)=0.5*(D0U2-D2U0)-0.5*(U2(I,J,K)*DU0+U0(I,J,K)*DU2)
     &                      +CS*(U2(I,J,K)*U0(I,J,K))
       DPc33(I,J,K)=CS-U0(I,J,K)/Time

       DPc11(I,J,K)=(-1.0)*D1U1-U1(I,J,K)*DU1+CS*(U1(I,J,K)**2+1.0)
       DPc22(I,J,K)=(-1.0)*D2U2-U2(I,J,K)*DU2+CS*(U2(I,J,K)**2+1.0)
       DPc12(I,J,K)=(-0.5)*(D2U1+D1U2)-0.5*(U1(I,J,K)*DU2+U2(I,J,K)*DU1)
     &                      +CS*(U1(I,J,K)*U2(I,J,K))
!       call checkSigma(U0(I,J,K), U1(I,J,K), U2(I,J,K),
!     &    DPc00(I,J,K), DPc01(I,J,K), DPc02(I,J,K),
!     &    DPc11(I,J,K), DPc12(I,J,K), DPc22(I,J,K), DPc33(I,J,K),
!     &    Itrigger)
       !if(Itrigger .eq. 1) then
       !  print*, D0U2, D2U0
       !  print*, U0(I,J+1,K), U0(I,J-1,K)
       !  stop
       !endif

         SiLoc(I,J,K)=CS*3.0  !   CS=(D0U0+D1U1+D2U2+U0(I,J,K)/Time)/3.0

         D0LnT=(DLog(Temp(I,J,K))-Dlog(Temp0(I,J,K)))/DT

C       If(abs(Accu-3.0).le.0.00001) then  !3pt formula
          D1LnT=(DLog(Temp(I+1,J,K))-Dlog(Temp(I-1,J,K)))/(2.0*DX)
          D2LnT=(DLog(Temp(I,J+1,K))-Dlog(Temp(I,J-1,K)))/(2.0*DY)
C       else if (abs(Accu-5.0).le.0.00001) then  !5pt formula
C         D1LnT=(DLog(Temp(I+1,J,K))*2.0d0/3.0d0
C     &      -Dlog(Temp(I-1,J,K))*2.0d0/3.0d0
C     &      -DLog(Temp(I+2,J,K))/12.0d0+DLog(Temp(I-2,J,K))/12.0d0 )/DX
C         D2LnT=(DLog(Temp(I,J+1,K))*2.0d0/3.0d0
C     &      -Dlog(Temp(I,J-1,K))*2.0d0/3.0d0
C     &      -DLog(Temp(I,J+2,K))/12.0d0+DLog(Temp(I,J-2,K))/12.0d0 )/DY
C        else
C           Print*, "wrong input for Accu,
C     &         Accu=3or5 for 3pt or 5pt cal of deriv. "
C        end if
       DlnT(I,J,K)=U0(I,J,K)*D0LnT +U1(I,J,K)*D1LnT +U2(I,J,K)*D2LnT         !DLnT term, added 02/2008
       DlnT(I,J,K)= DlnT(I,J,K)*ff  !make boundary works
 100   Continue

      Return
      End





C##################################################################
C----------------        EOS      -----------------------------
***************************************************************

       Subroutine InputEOSParameter

       Implicit Double Precision (A-H, O-Z)
      COMMON /PAREOS/ BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
     &                BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,
     &                PMAT1(0:300,0:250),PMAT2(0:300,0:250)
      COMMON /TAREOS/ TBNP01,TEPP01,TDBNP1,TDEPP1,NTBNP1,NTEPP1,
     &                TBNP02,TEPP02,TDBNP2,TDEPP2,NTBNP2,NTEPP2,
     &                TMAT1(0:300,0:250),TMAT2(0:300,0:250)
      COMMON /BAREOS/ BABNP01,BAEPP01,BADBNP1,BADEPP1,NBABNP1,NBAEPP1,
     &                BABNP02,BAEPP02,BADBNP2,BADEPP2,NBABNP2,NBAEPP2,
     &                BAMAT1(0:300,0:250),BAMAT2(0:300,0:250)

      COMMON /EOSSEL/ IEOS


*     Pressure(e,n)
      OPEN(10,FILE='EOS/aa1_p.dat',STATUS='OLD')
      READ(10,*) BNP01,EPP01
      READ(10,*) DBNP1,DEPP1,NBNP1,NEPP1
       DO J=0,NEPP1
        READ(10,*) (PMAT1(I,J),I=0,NBNP1)
       END DO
      CLOSE(10)

      OPEN(10,FILE='EOS/aa2_p.dat',STATUS='OLD')
      READ(10,*) BNP02,EPP02
      READ(10,*) DBNP2,DEPP2,NBNP2,NEPP2
      DO J=0,NEPP2
        READ(10,*) (PMAT2(I,J),I=0,NBNP2)
      END DO
      CLOSE(10)

*     Temperature(e,n)
      OPEN(10,FILE='EOS/aa1_t.dat',STATUS='OLD')
      READ(10,*) TBNP01,TEPP01
      READ(10,*) TDBNP1,TDEPP1,NTBNP1,NTEPP1
      DO J=0,NTEPP1
        READ(10,*) (TMAT1(I,J),I=0,NTBNP1)
      END DO
      CLOSE(10)

      OPEN(10,FILE='EOS/aa2_t.dat',STATUS='OLD')
      READ(10,*) TBNP02,TEPP02
      READ(10,*) TDBNP2,TDEPP2,NTBNP2,NTEPP2
      DO J=0,NTEPP2
        READ(10,*) (TMAT2(I,J),I=0,NTBNP2)
      END DO
      CLOSE(10)

*      baryon chemical potential (e,n)
      OPEN(10,FILE='EOS/aa1_mb.dat',STATUS='OLD')
      READ(10,*) BABNP01,BAEPP01
      READ(10,*) BADBNP1,BADEPP1,NBABNP1,NBAEPP1
      DO J=0,NBAEPP1
        READ(10,*) (BAMAT1(I,J),I=0,NBABNP1)
      END DO
      CLOSE(10)

      OPEN(10,FILE='EOS/aa2_mb.dat',STATUS='OLD')
      READ(10,*) BABNP02,BAEPP02
      READ(10,*) BADBNP2,BADEPP2,NBABNP2,NBAEPP2
      DO J=0,NBAEPP2
        READ(10,*) (BAMAT2(I,J),I=0,NBABNP2)
      END DO
      CLOSE(10)


      Return
      End



C#############################################################################
       Subroutine StressTensorZero(Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &                         PPI, NX0,NY0,NZ0,  NX,NY,NZ) !Set Zero ST for Ideal Hydro
C-----------Used for Ideal Hydro related stress tensors are set zero-------------------
        Implicit Double Precision (A-H, O-Z)
        Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)
c        Dimension Pi03(NX0:NX, NY0:NY, NZ0:NZ)
c        Dimension Pi13(NX0:NX, NY0:NY, NZ0:NZ)

          do 10 I=NX0,NX
          do 10 J=NY0,NY
          do 10 K=NZ0,NZ
             Pi00(I,J,K)=0.0
             Pi01(I,J,K)=0.0
             Pi02(I,J,K)=0.0
             Pi33(I,J,K)=0.0
             Pi11(I,J,K)=0.0
             Pi12(I,J,K)=0.0
             Pi22(I,J,K)=0.0
             PPI(I,J,K)=0.0
c             Pi03(I,J,K)=0.0
c             Pi13(I,J,K)=0.0
 10      Continue

       Return
       End

C#####################################################
       Subroutine ViscousCoefi8(Ed,Sd,PL,Temp,
     &  VCoefi,VCBeta,VRelaxT,etaTtp, VBulk, VCBeta0,VRelaxT0,XiTtP,
     &  Time,DX,DY,DT,NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
!--------call visous coeficient from from entropy
!--------bulk viscosity are included
        Implicit Double Precision (A-H, O-Z)

       Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
       Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure density
       Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

       Dimension VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient shear viscosity eta
       Dimension VCBeta(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient Beta2
       Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time
       Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

       Dimension VBulk(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient-bulk vicosity \xi
       Dimension VCBeta0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient Beta0
       Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient relaxation time \tau_PI
       Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008


       Dimension AA(NX0:NX, NY0:NY, NZ0:NZ) !

       Integer :: IVisflag

       Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Viscous Coefficient eta and beta2
       Common /ViscousBulk/ Visbulk, BulkTau,IRelaxBulk  ! Related to bulk Viscous Coefficient Xi and beta0

       Common /Tde/ Tde, Rdec1, Rdec2,TempIni !Decoupling Temperature !decoupling radius
       Common/R0Aeps/ R0,Aeps
       Common/R0Bdry/ R0Bdry
       Parameter (HbarC=0.19733d0) !for changing between fm and GeV ! Hbarc=0.19733=GeV*fm

      do 10 k=1,1
      do 10 j=NYPhy0-2,NYPhy+2 ! -2,NYPhy+2
      do 10 i=NXPhy0-2,NXPhy+2
      xx=DX*I
      yy=DY*J
      rr=sqrt(xx**2+yy**2)
      ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)
CSHEN==========================================================================
C=======for temperature dependent \eta/s=======================================
      if(IVisflag.eq.1) then
        ViscousC = ViscousCTemp(Temp(i,j,k))      !CSHEN: for temperature dependent \eta/s
      endif
CSHEN======end=================================================================
!------------- for shear pressure -----------
      If (ViscousC.ge.0.00001) then
        VCoefi(i,j,k)=ViscousC*Sd(i,j,k)
        VCBeta(i,j,k)=VisBeta*3.0/dmax1(Sd(i,j,k)*Temp(i,j,k),1e-30)
        VRelaxT(i,j,k)=1.0/dmax1(2.0*VCoefi(i,j,k)*VCBeta(i,j,k),1e-30)
        etaTtp(i,j,k)=(VCoefi(i,j,k)*Temp(i,j,k))*VRelaxT(i,j,k)  ! A(e+p)T/6 !(eta T/tau_Pi) for extra term in full I-S
        etaTtp(i,j,k)=dmax1(etaTtp(i,j,k),1e-30)
        VCoefi(i,j,k)=VCoefi(i,j,k)*ff
        VRelaxT(i,j,k)=VRelaxT(i,j,k)*ff
        If (Time > 0.8) Then
        If (etaTtp(i,j,k) .ne. etaTtp(i,j,k)) Then
          Print *, "Invalid etaTtp!"
          Print *, "i,j,k=",i,j,k
          Print *, "VRelaxT=",VRelaxT(i,j,k)
          Print *, "VCoefi=", VCoefi(i,j,k), "Temp=", Temp(i,j,k)
          Print *, "VCBeta=", VCBeta(i,j,k)
          Print *, "VisBeta=", VisBeta
          Print *, "Sd=", Sd(i,j,k), "Ed=", Ed(i,j,k)
        EndIf
        EndIf
      else
        VCoefi(i,j,k)=0.0
        VCBeta(i,j,k)=0.0
        VRelaxT(i,j,k)=0.0
        etaTtp(i,j,k)=0.0
      end if

!------------ for bulk pressure------------
      ttemp= Temp(i,j,k)*HbarC ! input function as Temp
      e3pep=(Ed(i,j,k)-3.0*PL(i,j,k))/(Ed(i,j,k)+PL(i,j,k))

      If (VisBulk.ge.0.00001) then
        !eta=ViscousC*Sd(i,j,k)
        !VBulk(i,j,k)=VisBulk*BulkAdSH0(eta,ttemp)
        VBulk(i,j,k) = ViscousZetasTemp(Ed(i,j,k)*HbarC)*Sd(i,j,k)
        If (IRelaxBulk.eq.0) then
          TTpi=DMax1(0.1d0, 120* VBulk(i,j,k)/DMax1(Sd(i,j,k),0.1d0))
          VRelaxT0(i,j,k)=1.0/TTpi
        else if (IRelaxBulk.eq.1) then
          VRelaxT0(i,j,k)=1.0/BulkTau
        else if (IRelaxBulk.eq.2) then
          VRelaxT0(i,j,k)=2*3.1415926*Temp(i,j,k)/1.5
        else if (IRelaxBulk .eq. 3) then
          TauPi = 9.0*VBulk(i,j,k)/(Ed(i,j,k) - 3.*PL(i,j,k))
          VRelaxT0(i,j,k) = 1.0/DMax1(0.1d0, TauPi)
        else
          Print*,'This option is not supported by this version'
          Print*,'IRelaxBulk'
          Stop
        end if
        VCBeta0(i,j,k)=VisBulkBeta*6.0/(Sd(i,j,k)*Temp(i,j,k))
        XiTtP(i,j,k)=(VBulk(i,j,k)*Temp(i,j,k))*VRelaxT0(i,j,k)  !(Xi T/tau_Pi)  for extra term in full I-S

        VBulk(i,j,k)=VBulk(i,j,k)*ff
        VRelaxT0(i,j,k)=VRelaxT0(i,j,k)*ff
      else
        VBulk(i,j,k)=0.0
        VCBeta0(i,j,k)=0.0
        VRelaxT0(i,j,k)=0.0
        XiTtP(i,j,k)=0.0
      end if

10    continue

      Return
      End

CSHEN======================================================================
C====eta/s dependent on local temperature==================================

      double precision function ViscousCTemp(TT)
      Implicit double precision (A-H, O-Z)
      parameter (HbarC=0.19733d0)
      Parameter (pi=3.1415926d0)

      Integer :: IVisflag

      Common /ViscousC / ViscousC,VisBeta, IVisflag ! Related to Shear Viscosity

      ViscousCTemp = ViscousC
      return
      end

C====zeta/s dependent on local temperature==================================
      double precision function ViscousZetasTemp(Ed)
      ! Ed input should be in unit of GeV/fm^3
      Implicit double precision (A-H, O-Z)
      Parameter (pi=3.1415926d0)

      de = 0.05*Ed
      p1 = PEOSL7(Ed - de/2.)
      p2 = PEOSL7(Ed + de/2.)
      cs2 = (p2 - p1)/de   !cs^2 = dP/de

      ViscousZetasTemp = 1./(8*pi)*(1./3. - cs2)

      return
      end

CSHEN=====end==============================================================



C----------------------------------------------------
       Double Precision FUNCTION BulkAdSH0(eta,tt)  !Ads-CFT Bulk-Vis=0 for HRG
       Implicit Double Precision (A-H, O-Z)
          !eta shear viscosity  Cs speed of Sound:Katz05 Lattice EOS4 here
          !tt temperature !GeV

           Cper=1.0d0/(exp((0.17d0-tt)/0.001d0)+1.0d0)
           Bper=1.0d0/(exp((0.180d0-tt)/0.01d0)+1.0d0)

       Cs1=(0.3d0*(1.0d0/(exp((0.180d0-TT)/0.025d0)+1.0d0))*Bper
     &  +(0.150d0*(1.0-exp(-((0.185d0-TT)/0.037d0)**2))+0.035d0)
     &  *(1-Bper))*Cper    !QGP

CCcc       Width=0.02d0
C       Width=0.008d0      !sharp head
C       CC0=0.07458569d0
C       Cs2=((1.0d0/3.0d0-CC0)*(1.0-exp(-((0.172d0-TT)/Width)**2))+CC0)
C     &          *(1-Cper)  !HRG

       P2=0.07458569d0    !round head
       Cs2=((1.0d0/3.0d0-P2)/(exp((tt-0.157d0)/0.002d0)+1.0d0)+P2)
     &        *(1-Cper)

       Cs=Cs1+Cs2

           Zeta=eta*2.0*(1.0/3.0-Cs)
           BulkAdSH0=Zeta
       return
       end


C#######################################################
      Subroutine initialES2( Ed,Bd,Sd, DX,DY,DZ,DT,
     &    NX0,NY0,NZ0, NX,NY,NZ,NXPhy0,NYPhy0 ,NXPhy,NYPhy)
      Implicit Double Precision (A-H, O-Z)
      Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Dimension Bd(NX0:NX, NY0:NY, NZ0:NZ) !Baryon density
      Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !Pressure
      Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Dimension CMu(NX0:NX, NY0:NY, NZ0:NZ) !Local chemical potential

      Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm
      Parameter (pi=3.1415926d0)
      Parameter (gt=169.0d0/4.0d0)!total freedom of Quarks and Gluond  Nf=2.5 !change another in EntropyTEm

      Double Precision resultingEd

      COMMON /EOSSEL/ IEOS   !Type of EOS
      COMMON /IEin/ IEin  !parameter for initializtion by entropy 1 or energy 0
      Common /bb/ b  !impact parameter
      Common /EK/ EK, HWN !EK(T0) constant related to energy density
                           !HWN percent of Wounded Nucleon
      Common /RxyBlock/ Rx2, Ry2

      Integer, Parameter :: RegEOSdatasize = EOSDATALENGTH  !converted EOS table size
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                 !total rows of energy density

      Double Precision, External :: SEOSL7

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata !CSHEN: for EOS from tables
      common /EOSdatastructure/ EOSe0, EOSde, EOSne


      HBC=1.0-HWN       !HBC Percent of Binary Collision

      DN00=HWN*CWNBC (0.0d0,0.0d0,0.0d0,1.0d0)
     &           +HBC*CWNBC (0.0d0,0.0d0,0.0d0,2.0d0)
      Con1=EK/DN00  !factor K


      ConsG=gt*Pi**2/30.0
      ConsS0=(4.0/3.0)*ConsG**(0.25)
      ConsE1=(0.75)**(4.0/3.0)*(1.0/ConsG)*1.0/3.0

      If(IEin.eq.0 ) then  ! intialliztion by energy then find entropy
        Print*,'Con1,Ek,DN00', Con1,EK,DN00
        do 5 k=1,1
        do 5 j=NYPhy0,NYPhy
        do 5 i=NXPhy0,NXPhy
          xx=DX*i
          yy=DY*j
          zz=DZ*K
          Ed(i,j,k)=Con1*(HWN*CWNBC (xx,yy,b,1.0d0) !WN   !changed for 2d
     &                    +HBC*CWNBC (xx,yy,b,2.0d0))!BC
          If(IEOS.eq.2) then
            Sd(i,j,K)=ConsS0*ED(i,j,k)**(0.75)
          End if
 5      continue
        If(IEOS.ne.2) then
          call EntropyTemp3 (Ed,PL, Temp,CMu,Sd,
     &      NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
        End if
      End if ! IEin.eq.0


      If(IEin.eq.1 ) then !initializtion by entropy then find energy

        ! Ek is used as "Sk", where the conversion from GeV/fm^3 to fm^-4 is unneccessary.
        ! Here the unit for "Ek" from the input should be fm^-3.
        ! Note that "Ek" at the initialization stage is this "Ek/HbarC".
        Con1 = Con1*HbarC
        Print*,'Con1,Sk,DN00', Con1,EK*HbarC,DN00
        Do 15 k=1,1
        Do 15 j=NYPhy0,NYPhy
        Do 15 i=NXPhy0,NXPhy
          xx=DX*i
          yy=DY*j
          zz=DZ*K
          Sd(i,j,k)=Con1*(HWN*CWNBC (xx,yy,b,1.0d0) !WN   !changed for 2d
     &                    +HBC*CWNBC (xx,yy,b,2.0d0))!BC
          If(IEOS.eq.2) then
            Ed(i,j,k)=ConE1*Sd(i,j,k)**(4.0/3.0)
          End if
15      Continue
        If(IEOS==7) Then
          Do I = NXPhy0,NXPhy
          Do J = NYPhy0,NYPhy
            Call invertFunctionD(SEOSL7, 0D0, 1000D0, 1D-3, 0D0,
     &                        Sd(I,J,NZ0), resultingEd)
            Ed(I,J,NZ0) = resultingEd/HbarC ! to fm^(-4)
          End Do
          End Do
        End If
        Print *, Ed(0,0,NZ0)
      End if ! IEin.eq.7

       Return
       End



C#######################################################
       Subroutine initialES2Gaussian(Ed,Bd,Sd, DX,DY,DZ,DT,
     &    NX0,NY0,NZ0, NX,NY,NZ,NXPhy0,NYPhy0 ,NXPhy,NYPhy)
        Implicit Double Precision (A-H, O-Z)
       Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension Bd(NX0:NX, NY0:NY, NZ0:NZ) !Baryon density
       Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

       Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !Pressure
       Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
       Dimension CMu(NX0:NX, NY0:NY, NZ0:NZ) !Local chemical potential

       Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm
       Parameter (pi=3.1415926d0)
       Parameter (gt=169.0d0/4.0d0)!total freedom of Quarks and Gluond  Nf=2.5 !change another in EntropyTEm

       COMMON /EOSSEL/ IEOS   !Type of EOS
       COMMON /IEin/ IEin  !parameter for initializtion by entropy 1 or energy 0
       Common /bb/ b  !impact parameter
       Common /EK/ EK, HWN !EK(T0) constant related to energy density
                            !HWN percent of Wounded Nucleon
        Common /RxyBlock/ Rx2, Ry2
          HBC=1.0-HWN       !HBC Percent of Binary Collision

        DN00=HWN*CWNBC (0.0d0,0.0d0,0.0d0,1.0d0)
     &           +HBC*CWNBC (0.0d0,0.0d0,0.0d0,2.0d0)
        Con1=EK/DN00  !factor K
        Print*,'Con1,Ek,DN00', Con1,EK,DN00

        ConsG=gt*Pi**2/30.0
        ConsS0=(4.0/3.0)*ConsG**(0.25)
        ConsE1=(0.75)**(4.0/3.0)*(1.0/ConsG)*1.0/3.0

      If(IEin.eq.0 ) then  ! intialliztion by energy then find entropy
        do 5 k=1,1
         do 5 j=NYPhy0,NYPhy
          do 5 i=NXPhy0,NXPhy
                 xx=DX*i
                 yy=DY*j
                 zz=DZ*K
             Ed(i,j,k)=Ek*1.0
     &               *Exp(-xx*xx/(Rx2*2.0)-yy*yy/(Ry2*2.0))
            If(IEOS.eq.2) then
            Sd(i,j,K)=ConsS0*ED(i,j,k)**(0.75)
            End if
 5      continue
        If(IEOS.ne.2) then
            call EntropyTemp3 (Ed,PL, Temp,CMu,Sd,
     &      NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
        End if
      end if

      If(IEin.eq.1 ) then !initializtion by entropy then find energy
        Print *, "Gaussain initialization by entropy not surpported."
        Stop
      End if ! IEin.eq.7

       Return
       End
C###########################################################################
       Subroutine EntropyTemp3 (Ed,PL, Temp,CMu,Sd,
     &       NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
C------call Temperatuen mu  Entropy from ee pp-----------
        Implicit Double Precision (A-H, O-Z)
       Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension Bd(NX0:NX, NY0:NY, NZ0:NZ) !energy density
       Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !Pressure

       Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
       Dimension CMu(NX0:NX, NY0:NY, NZ0:NZ) !Local chemical potential
       Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

       Parameter (pi=3.1415926d0)
       Parameter (gt=169.0d0/4.0d0)!total freedom of Quarks and Gluond  Nf=2.5  !change another in InitialES
       Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm

      COMMON /PAREOS/ BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
     &                BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,
     &                PMAT1(0:300,0:250),PMAT2(0:300,0:250)
      COMMON /TAREOS/ TBNP01,TEPP01,TDBNP1,TDEPP1,NTBNP1,NTEPP1,
     &                TBNP02,TEPP02,TDBNP2,TDEPP2,NTBNP2,NTEPP2,
     &                TMAT1(0:300,0:250),TMAT2(0:300,0:250)
      COMMON /BAREOS/ BABNP01,BAEPP01,BADBNP1,BADEPP1,NBABNP1,NBAEPP1,
     &                BABNP02,BAEPP02,BADBNP2,BADEPP2,NBABNP2,NBAEPP2,
     &                BAMAT1(0:300,0:250),BAMAT2(0:300,0:250)

       COMMON /EOSSEL/ IEOS

      Do 1001 k=1,1                           !do 10 is in the unit of fm-1
      Do 1001 j=NYPhy0-2,NYPhy+2 ! -2,NYPhy
      Do 1001 i=NXPhy0-2,NXPhy+2
        Ed(i,j,k) = max(1e-30, Ed(i,j,k))
 1001 Continue




        If(IEOS.eq.2) then !Ideal QGP Gas (Ref. Wong Book)

        a=(gt*Pi**2)/90.0
        ConsT=(1.0/(3.0*a))**0.25    !ConsT=(30.0/(gt*Pi**2))**0.25
        ConsS=4.0d0*a                    !ConsS=4.0/3.0*gt*pi**2/30.0

           Print *, 'ConsT', ConsT
          do 10 k=1,1                           !do 10 is in the unit of fm-1
          do 10 j=NYPhy0-2,NYPhy+2 ! -2,NYPhy
          do 10 i=NXPhy0-2,NXPhy+2
             PL(i,j,k)=Ed(i,j,k)/3.0
             Temp(i,j,k)=ConsT*(Ed(i,j,k)**0.25)
             Sd(i,j,k)=ConsS*Temp(i,j,k)**3
          !If (Sd(i,j,k) .ne. Sd(i,j,k)) Then
          !  Print *, "Invalid Sd!"
          !  Print *, "i,j,k=", i,j,k
          !  Print *, "Ed=", Ed(i,j,k)
          !  Print *, "Sd=", Sd(i,j,k)
          !  Print *, "Temp=", Temp(i,j,k)
          !  Print *, "Pl=", PL(i,j,k)
          !EndIf
 10       continue
          Print*,'ConsS=',ConsS


       else if(IEOS.eq.0) then
       DO 98 I=NXPhy0-2,NXPhy+2
       DO 99 J=NYPhy0-2,NYPhy+2
         X=X0+I*DX
         Y=Y0+J*DY
         ee=Ed(I,J,NZ0)*HbarC      !fm-4 to GeV/fm3  peter unit
         ba=0.0 !Bd(I,J,NZ0)

           IF (ee.LT.24.) THEN
           TTEMP=EOUT(ba,ee,0.0d0,
     &          TBNP01,TEPP01,TDBNP1,TDEPP1,NTBNP1,NTEPP1,
     &          TBNP02,TEPP02,TDBNP2,TDEPP2,NTBNP2,NTEPP2,
     &          TMAT1,TMAT2)
           pp=PEPS(ba,ee)
           AMU=EOUT(ba,ee,0.0d0,
     &         BABNP01,BAEPP01,BADBNP1,BADEPP1,NBABNP1,NBAEPP1,
     &         BABNP02,BAEPP02,BADBNP2,BADEPP2,NBABNP2,NBAEPP2,
     &         BAMAT1,BAMAT2)
           ELSE
           TTEMP=( (ee-.3642)*120./(PI**2.) /169. *(HBARC**3.))**(.25)
           pp=(ee-4.*.3642)/3.
           AMU =  0.
*          with N_f=2.5    .3642 is our bag-constant in GeV/fm^3
           END IF
         IF(TTEMP.GT..00001) THEN
         ss=(ee+pp-AMU*ba)/TTEMP
         ELSE
         ss=0.
         END IF

       PL(I,J,NZ0)=pp/HbarC     !GeV/fm3 to fm-4
       Temp(I,J,NZ0)=TTEMP/HbarC !GeV/fm3 to fm-4
       CMu(I,J,NZ0)=AMU !HbarC?
       Sd(I,J,NZ0)=ss           !fm-3

 99   CONTINUE
 98   CONTINUE
      else if (IEOS.eq.4) then  ! the lattice based EOS fm

          do 110 k = 1,1                           !do 10 is in the unit of fm-1
          do 110 j = NYPhy0-2,NYPhy+2 ! -2,NYPhy
          do 110 i = NXPhy0-2,NXPhy+2
           ee = Ed(i,j,k)
           TT = TEIEOS4(ee)  !change the double precison fun for freezeourt corrispondingly
           pp = PEIEOS4(ee)
           ss = (pp+ee)/TT

           PL(i,j,k)   = pp
           Temp(i,j,k) = TT
           Sd(i,j,k)   = ss
110       continue

      else if (IEOS.eq.5) then      !CSHEN SM-EOS Q
          do 112 k = 1,1
          do 112 j = NYPhy0-2,NYPhy+2
          do 112 i = NXPhy0-2,NXPhy+2
            ee = Ed(i,j,k)
            TT = TEIEOS5pp(ee)
            pp = PEIEOS5pp(ee)
            ss = (pp+ee)/TT

            PL(i,j,k)   = pp
            Temp(i,j,k) = TT
            Sd(i,j,k)   = ss
112        continue

      else if (IEOS.eq.6) then  ! CSHEN New EOS from Lattice

          do 111 k = 1,1                           !do 10 is in the unit of fm-1
          do 111 j = NYPhy0-2,NYPhy+2 ! -2,NYPhy
          do 111 i = NXPhy0-2,NXPhy+2
           ee = Ed(i,j,k)*HbarC     !ee in Gev/fm^3
           TT = TEOSL6(ee)/HbarC    !TEOSL6 in GeV, TT in fm^-1
           pp = PEOSL6(ee)/HbarC          !PEOSL6 in GeV/fm^3, pp in fm^-4
           ss = SEOSL6(ee)                !SEOSL6 in fm^-3, ss in fm^-3

           PL(i,j,k)   = pp
           Temp(i,j,k) = TT
           Sd(i,j,k)   = ss
111       continue
      else if (IEOS.eq.7) then  ! CSHEN EOS from tables

          do k = 1,1
          do j = NYPhy0-2,NYPhy+2
          do i = NXPhy0-2,NXPhy+2
           ee = Ed(i,j,k)*HbarC     !ee in Gev/fm^3
           TT = TEOSL7(ee)/HbarC    !TEOSL7 in GeV, TT in fm^-1
           pp = PEOSL7(ee)/HbarC    !PEOSL7 in GeV/fm^3, pp in fm^-4
           ss = SEOSL7(ee)          !SEOSL7 in fm^-3, ss in fm^-3

           PL(i,j,k)   = pp
           Temp(i,j,k) = TT
           Sd(i,j,k)   = ss
          enddo   ! i loop
          enddo   ! j loop
          enddo   ! k loop


      else
      Print*,'alart EntropyTemp2 IEOS'
      end if

      Return
      end

***************************************************************************
      Double Precision FUNCTION PEPS(BN,EP) ! [EP]=GeV/fm^3
*****************************************************************
** THIS FUNCTION GIVES THE EQUATION OF STATE P=P(BN,EPS)        *
** BY INTERPOLATING THE TABLE OF VALUES OF PRESSURE GIVEN
** BY THE MATRICES 'PMAT1' AND 'PMAT2' AT COMMON /PAREOS/.
*****************************************************************
       Implicit Double Precision (A-H, O-Z)
      INTEGER RAPPSIZE
      PARAMETER (RAPPSIZE=76)
       Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm

      COMMON /PAREOS/ BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
     &                BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,
     &                PMAT1(0:300,0:250),PMAT2(0:300,0:250)

C      COMMON /RAPPEOS/ REOS(1:13,1:RAPPSIZE), pRAPP(1:RAPPSIZE),
C     &                 TRAPP(1:RAPPSIZE), sRAPP(1:RAPPSIZE),
C     &                 muNRAPP(1:RAPPSIZE), muNeRAPP(1:RAPPSIZE),
C     &                 mupiRAPP(1:RAPPSIZE), muKRAPP(1:RAPPSIZE),
C     &                 muetRAPP(1:RAPPSIZE),
C     &                 eRAPP(1:RAPPSIZE), esi2RAPP(1:RAPPSIZE),
C     &                 Tei2RAPP(1:RAPPSIZE), pei2RAPP(1:RAPPSIZE),
C     &                 mei2RAP(1:RAPPSIZE), meei2RAP(1:RAPPSIZE),
c     &                 mepi2RAP(1:RAPPSIZE), meKi2RAP(1:RAPPSIZE),
C     &                 meti2RAP(1:RAPPSIZE)

      COMMON /EOSSEL/ IEOS

      IF ((IEOS .EQ. 0)) THEN
       CALL POINT(BN,EP,F)        !Peter unit
       PEPS = F
      END IF

      IF ((IEOS .EQ. 2)) THEN
       PEPS = EP/3.0
      END IF

      IF ((IEOS .EQ. 4)) THEN       !peter Gev/fm-3
       ee   = Ep/HbarC              !my fm-4
       PEPS = PEIEOS4(ee)*HbarC     !GeV/fm^3
      END IF

      IF ((IEOS .eq. 5)) THEN       !CSHEN SM-EOS Q
       ee   = Ep/HbarC              !1/fm^4
       PEPS = PEIEOS5pp(ee)*HbarC   !GeV/fm^3
      end if

      IF ((IEOS .EQ. 6)) THEN       !CSHEN NEW EOS
       ee = Ep                      !ee in GeV/fm^3
       PEPS = PEOSL6(ee)            !GeV/fm^3
      end if

      IF ((IEOS .EQ. 7)) THEN       !CSHEN EOS from tables
       ee = Ep                      !ee in GeV/fm^3
       PEPS = PEOSL7(ee)            !GeV/fm^3
      end if

      RETURN
      END




C#####################################################################################
C###############      Lattice  EOS   Katz05   data   #################################
C#####################################################################################
       Double Precision Function PEIEOS4(ee)  !close to Lattice  Katz05 data
       Implicit Double Precision (A-H, O-Z)
       Aperp=1.0d0/(exp((4.0d0-ee)/1.0d0)+1.0d0)
       Cperp=1.0d0/(exp((1.8d0-ee)/1.5d0)+1.0d0)
       pp=(-2.2d0+2.655d0*LOG(1.18136d0+exp(0.1111d0*ee)))*Aperp
     &    +0.18d0*ee*(1.0d0-Cperp)+(0.40d0*Cperp)*(1-Aperp)
     &    -0.0886002645519078d0
         PEIEOS4=pp
       return
       end

C--------------------------------------------------------------------------------------
       Double Precision Function TEIEOS4(ee) !for Lattice
*       Double Precision Function TEIEOS4LL(ee) !for Latiice
       Implicit Double Precision (A-H, O-Z)
        Aper=1.0d0/(exp((60.0d0-ee)/60.0d0)+1.0d0)
        Bper=1.0d0/(exp((1.5d0-ee)/7.0d0)+1.0d0)
        Cper=1.0d0/(exp((9.0d0-ee)/3.5d0)+1.0d0)
        TT =((0.087d0*ee)**(0.25d0)*Aper
     &     +(0.102d0*ee)**(0.25d0)*(1.0d0-Aper)
     &              +0.15d0*(1.0d0-Bper))*Cper
     &     +(0.15d0*ee)**(0.15d0)*(1-Cper)
        TEIEOS4=TT
       return
       end
*     Temperature-energy density relation
*     T>165MeV modeling from Katz 05 Lattice (QM05)
*     T<165MeV free Massive resonance gas from Peter code
*     s=(e+p)/T


C####################################################################################
C##################      Peter EOS    SM-EOS Q      #################################
C####################################################################################
C----------------------------------------------------------------------
       Double Precision Function TEIEOS5pp(ee) !for peter EOS0 continous modeling
*       Double Precision Function TEIEOS4(ee) !for peter EOS0 continous modeling
       Implicit Double Precision (A-H, O-Z)
        Bper=1.0d0/(exp((8.0d0-ee)/3.0d0)+1.0d0)
        Cper=1.0d0/(exp((5.0d0-ee)/15.0d0)+1.0d0)
        TT=0.0+(0.07d0*ee)**(0.25d0)*Bper
     &     +(6.2d0*ee)**(0.16d0)*(1-Cper)*(1-Bper)
        TEIEOS5pp=TT
*      s=(e+p)/T
       return
       end
C-----------------------------------------------------------------------
       Double Precision Function PEIEOS5pp(ee)   !close to Peter---2nd order phase transition
*       Double Precision Function PEIEOS4(ee)   !close to Peter---2nd order phase transition
       Implicit Double Precision (A-H, O-Z)    !lower EOS E2-15fm close Peter IEOS0  other id the same
       Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm
C          EP=ee*HbarC      !fm-4 to GeV/fm3  peter unit
C          BN=0.0 !Bd(I,J,NZ0)
C          CALL POINT(BN,EP,F)
C          pp=F/HbarC     !GeV/fm3 to fm-4

         If(ee.gt.30.0d0)then        !Q (7.37 0.5)    H Exp(2.7 1.0) (2.45  0.1)
           pp=0.334d0*ee-2.46085238d0
         else if(ee.gt.5.0d0)then
           pp=(-2.46085238d0+0.167d0*LOG(2.696674d7+exp(2.0d0*ee)))!
         else
           pp=0.152000165d0*ee
     &       -0.0152d0*LOG(4.19314d10+Dexp(10.0d0*ee))
     &       -0.0248909*Dexp(-2.7d0*ee) +0.3966722716141761d0
         end if
**           Cs=1.0/3.0
**           Cs=0.333333d0*(1.0d0/(exp((7.37d0-ee)/0.5d0)+1.0d0))       !Q (7.37 0.5)
**           Cs=Dexp(-2.7d0*(ee+1.0d0))
**     &              +0.152d0*(1.0d0/(exp((ee-2.44593)/0.1d0)+1.0d0))  !H Exp(2.7 1.0) (2.45  0.1)
         PEIEOS5pp=pp
       return
       end


C#############################################################################################
C--------------------------Lattice  EOS   Pasi-09 ---------------------------------------
C------------------------- added by Shen Chun  ------------------------------------------
C-----------------------FIT FUNCTIONS GOT FROM TOM REILY---------------------------------
C#############################################################################################
       Double Precision Function PEOSL6(ee)  ! for lattice P(e)
C------  p(e)  p, e GeV/fm^3
       Implicit Double Precision (A-H, O-Z)     !h
      If (ee.lt.0.5028563305441270d0) Then
       pp=0.3299d0*(Exp(0.4346d0*ee) - 1)

C      Else If (ee.ge.0.5028563305441270d0 .AND.      !Tom original
C     & ee.lt.1.8408643375520979d0) Then
C	pp=0.0000001024d0*Exp(6.041d0*ee)
C     &  + 0.007273d0+0.14578d0*ee

CSHEN=================================================================
C====smooth cs======================================================
      Else If (ee.ge.0.5028563305441270d0 .AND.
     & ee.lt.1.62d0) Then
       pp=0.0000001024d0*Exp(6.041d0*ee)
     &  + 0.007273d0+0.14578d0*ee
      ELse if (ee.ge.1.62d0 .AND. ee.lt.1.86d0) then
       pp=0.30195*Exp(0.31309*ee)-0.256232
C      Else If (ee.ge.1.8408643375520979d0 .AND.      !Tome original
       Else if (ee.ge.1.86d0 .AND.
C===========end====================================================

     & ee.lt.9.9878355786273545d0) Then
       pp= 0.332d0*ee-0.3223d0*(ee**0.4585d0)
     &  - 0.003906d0*ee/Exp(0.05697d0*ee)
     &      + (0.1167d0/(ee**1.233d0)
     &  + 0.1436d0*ee/Exp(0.9131d0*ee))
      Else If (ee.ge.9.9878355786273545d0) Then
      pp = 0.3327d0*ee-0.3223d0*(ee**0.4585d0)
     &  - 0.003906d0*ee/Exp(0.05697d0*ee)
      End If
         PEOSL6=pp

       return
       end


C--------------------------------------------------------------------------------------
       Double Precision Function SEOSL6(ee) !for Lattice S(e)
       Implicit Double Precision (A-H, O-Z)
C------  s(e)   e GeV/fm^3  s fm^-3

      If (ee.lt.0.1270769021427449d0) Then
	ss=12.2304d0*(ee**1.16849d0)
      Else If (ee.ge.0.1270769021427449d0 .AND.
     &  ee.lt.0.4467079524674040d0) Then
	ss=11.9279d0*(ee**1.15635d0)
      Else If (ee.ge.0.4467079524674040d0 .AND.
     &  ee.lt.1.9402832534193788d0) Then
	ss=0.0580578d0 + 11.833d0*(ee**1.16187d0)
      Else If (ee.ge.1.9402832534193788d0 .AND.
     &  ee.lt.3.7292474570977285d0) Then
	ss=18.202d0*ee - 63.0218d0
     & - 4.85479d0*Exp(-0.0000000000272407d0*(ee**4.54886d0))
     & + 65.1272d0/(ee**0.128012d0)
     & /Exp(0.00369624d0*(ee**1.18735d0))
     & - ((4.75253d0/(ee**1.18423d0))-0.999986d0)
      Else If (ee.ge.3.7292474570977285d0) Then
	ss=18.202d0*ee - 63.0218d0
     & - 4.85479d0*Exp(-0.0000000000272407d0*(ee**4.54886d0))
     & + 65.1272d0/(ee**0.128012d0)
     & /Exp(0.00369624d0*(ee**1.18735d0))
      End If
        SEOSL6=ss**(3.0/4.0)

       return
       end

       Double Precision Function TEOSL6(ee) !for Lattice     T(e)
C------  T(e)   e GeV/fm^3  T GeV

      Implicit Double Precision (A-H, O-Z)

      If (ee.ge.0.5143939846236409d0) Then  ! 10^-16 accuracy
	tt=(ee+PEOSL6(ee))/(SEOSL6(ee))
      Else If (ee.lt.0.5143939846236409d0) Then
	tt=0.203054d0*(ee**0.30679d0)
      End If

      TEOSL6=tt

       return
       end



C######################################################################

       Double Precision Function CWNBC (xx,yy,b,WNBC)
C----- initial Wounded Nucleons density (WN) or Binary Collissions (BC)---
C       Glauber Model kjC       A Nuclei Number    b, impact parameter  Si0, Cross Section for NN
c       Character WNBC,WN,BC
        Implicit Double Precision (A-H, O-Z)
       Common /AWNBC/ A,Si0 !A Nuclei Number  Si0, Cross Section for NN

            TA=Thickness(xx+b/2.0,yy)
            TB=Thickness(xx-b/2.0,yy)

        If(abs(WNBC-1.0).le.0.001 ) then ! for Wounded Nucleons
          CWNBC=TA*(1.0-(1.0-(Si0*TB/A))**A)
     &             +TB*(1.0-(1.0-(Si0*TA/A))**A)
        end if

        If(abs(WNBC-2.0).le.0.001) then ! for Binary Collision
          CWNBC=Si0*TA*TB
        end if

       Return
       End


C####################################################################
       Function Thickness(xx,yy)
C----------Nuclear Thickness Function------------
       Implicit Double Precision (A-H, O-Z)
       integer iz, iz_min, iz_max
       Common /thick/ TRo0, TEta, TRA  !Para in Nuclear Thickness Function

       Thickness=0.0
       ddz=TRA/100.0
       iz_min = -5.0*TRA/ddz
       iz_max = 5.0*TRA/ddz
       do 10 iz = iz_min, iz_max, 1
            zz = iz*ddz
            r=sqrt(xx**2+yy**2+zz**2)
            RoAr=TRo0/(Exp((r-TRA)/TEta)+1.0)
            Thickness=Thickness+RoAr*ddz
 10    continue

       Return
       End
!C#################################################################
!       Subroutine CheckThickness
!C----------Check Thickness Func.
!       Implicit Double Precision (A-H, O-Z)
!       Common /thick/ TRo0, TEta, TRA  !Para in Nuclear Thickness Function
!         ARo=0.0
!         dx=TRA/10.0
!         dy=TRA/10.0
!         i=0
!         do 10 x=(-5.0)*TRA,(5.0)*TRA,dx
!            i=i+1
!           If(mod(i,10).eq.0)  Print*,x, 'check Thickness'
!          do 10 y=(-5.0)*TRA,(5.0)*TRA,dy
!              ARo=ARo+Thickness(x,y)*dx*dy
! 10      continue
!
!       Return
!       End



C###################################################################
      SUBROUTINE FYSCOR(T00,T01,T02,C)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
      Implicit Double Precision (A-H, O-Z)
      C = 1.0
      T00=dMAX1(0.0d0,T00)
      W=SQRT(T01*T01+T02*T02)
      IF (W.GT.T00) THEN
        C=0.999*T00/W
        T01=C*T01
        T02=C*T02
      END IF

      RETURN
      END

C###################################################################
      SUBROUTINE FYSCOR_velocity(T00,T01,T02,Ed,PL,BulkPi)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
      Implicit none
      double precision T00, T01, T02, Ed, PL, BulkPi
      double precision cstilde2
      double precision M0, M
      double precision scaleFactor
      double precision temp1

      cstilde2 = PL/Ed
      T00 = dMax1(0.0d0, T00)
      M0 = T00
      M = sqrt(T01**2 + T02**2)
      if (M0 .lt. M) then
        scaleFactor = 0.999*M0/M
        T01 = scaleFactor*T01
        T02 = scaleFactor*T02
      endif
      temp1 = 2.*sqrt(cstilde2)*M - M0*(1.+cstilde2)
      if (BulkPi .lt. temp1) then
        BulkPi = temp1 + 1D-3*abs(temp1)
      endif
      
      RETURN
      END

C###################################################################
      SUBROUTINE FYSCOR2(T00,T01,T02,PL,C)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
        Implicit Double Precision (A-H, O-Z)
      C = 1.0
      T00=dMAX1(0.0d0,T00)
      W=SQRT(T01*T01+T02*T02)
      IF (W.GT.(T00+PL)) THEN
        C=0.999*(T00+PL)/W
        T01=C*T01
        T02=C*T02
      END IF

      RETURN
      END

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C###################################################################
      SUBROUTINE FYSCOR_withBulkvis(T00,T01,T02,BulkPi)
*************************************************************
** THIS ROUTINE ASSERTS THAT THE VELOCITY OF THE ENERGY     *
** FLOW CAN'T EXCEED 1 AND THAT ENERGY DENSITY IS POSITIVE  *
** OR ZERO.                                                 *
*************************************************************
      Implicit Double Precision (A-H, O-Z)
      temp = dmax1(0.0d0, T00)
      temp = temp*(temp + BulkPi)
      temp = dmax1(0.0d0, temp)
      W = sqrt(T01*T01+T02*T02)
      IF (W .GT. sqrt(temp)) THEN ! regulate the violation
        C=0.999*sqrt(temp)/W
        T01=C*T01
        T02=C*T02
      END IF

      RETURN
      END

C###################################################################
      Subroutine dpSc8(TT00,TT01,TT02,ScT00,ScT01,ScT02,Vx,Vy,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22, PScT00,PScT01,PScT02,PScT33,
     &  PScT11,PScT12,PScT22,etaTtp0,etaTtp,  PPI,PISc, XiTtP0,XiTtP,
     &  U0,U1,U2, PU0,PU1,PU2,SxyT,Stotal,StotalBv,StotalSv,
     &  Ed,PL,Bd,Sd,Temp0,Temp,CMu, T00,T01,T02, IAA,CofAA,Time,DX,DY,
     &  DZ,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NX,NY0,NY,NZ0,NZ, PNEW,NNEW)
      !PNEW NNEW related to root finding

        Implicit Double Precision (A-H, O-Z)
        Dimension TT00(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension ScT00(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension TT01(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension ScT01(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension TT02(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension ScT02(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)
        Double Precision VxTmp(NX0:NX, NY0:NY, NZ0:NZ)
        Double Precision VyTmp(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension dummy(NX0:NX, NY0:NY, NZ0:NZ)

        Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
        Dimension PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
        Dimension PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

        Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT00(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT01(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT02(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT33(NX0:NX, NY0:NY, NZ0:NZ)  !

        Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT11(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT12(NX0:NX, NY0:NY, NZ0:NZ)  !
        Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
        Dimension PScT22(NX0:NX, NY0:NY, NZ0:NZ)  !

       Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)    !Bulk pressure
       Dimension PISc(NX0:NX, NY0:NY, NZ0:NZ)  !

C--------------------------------------------
        Dimension T00(NX0:NX, NY0:NY, NZ0:NZ)! ideal T00  energy momentum tensor
        Dimension T01(NX0:NX, NY0:NY, NZ0:NZ)! ideal T01
        Dimension T02(NX0:NX, NY0:NY, NZ0:NZ)! ideal T02
C-------------------------------------------
        Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
        Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
        Dimension Bd(NX0:NX, NY0:NY, NZ0:NZ) !net baryon density
        Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

        Dimension Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature  in last time step
        Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
        Dimension CMu(NX0:NX, NY0:NY, NZ0:NZ) !Local chemical potential

        Dimension VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient
        Dimension VCBeta(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient  Beta 2
        Dimension VRelaxT(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient relaxation Time
        Dimension SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \theta
        Dimension DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms  !added 02/2008

       Dimension VBulk(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient-bulk viscosity \xi
       Dimension VCBeta0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient Beta0
       Dimension VRelaxT0(NX0:NX, NY0:NY, NZ0:NZ) !viscous coefficient relaxation time \tau_PI

        Dimension etaTtp(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008
        Dimension etaTtp0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (eta T)/tau_pi terms in I-S eqn 02/2008

        Dimension XiTtP(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi terms in full I-S bulk eqn 08/2008
        Dimension XiTtP0(NX0:NX, NY0:NY, NZ0:NZ)  !extra (Xi T)/tau_Pi in last time step

        Dimension DDU0(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DDU2(NX0:NX, NY0:NY, NZ0:NZ) !


        Dimension DPc00(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension DPc11(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
        Dimension DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension CC(NX0:NX, NY0:NY, NZ0:NZ) !

        Dimension IAA(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

       Integer :: IVisflag

       Common /ViscousC / ViscousC,VisBeta, IVisflag
       Common /ViscousBulk/ Visbulk, BulkTau,IRelaxBulk

        DIMENSION PNEW(NNEW)!something related to root finding
       Parameter(XC0=1.0d-15, AAC=1.0d-16) !root finding accuracy

        LOGICAL V1FOUND,V2FOUND,V3FOUND,V4FOUND
        COMMON /EOSSEL/ IEOS
        Common /Newtonalart/ AI,AJ,AK,AE
        Common/R0Aeps/ R0,Aeps
        Common /R0Bdry/ R0Bdry

       Common/dxdy/ ddx, ddy
       Common /TT0/ TT0

       Parameter (HbarC=0.19733d0) !for changing between fm and GeV ! Hbarc=0.19733=GeV*fm

       Common /Nsm/ Nsm
       Common /Accu/Accu
       common/Edec1/Edec1

      ! ----- Use in root search -----
      Double Precision :: RSDM0, RSDM, RSPPI, RSee
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

         Nsm=5

        if(NZ0.ne.NZ)then
          Print *,' VSc2d , a 2-dimensinal subroutine '
          Stop
        end if


      PNEW(1)=XC0   !dimension related to newton root finding
      PNEW(2)=33335.5
      PNEW(3)=0.3
      PNEW(4)=0.3
      AAC0=AAC
      AAC10=AAC*0.1


!   ---Zhi-Changes---
!-------Regulate Pi(mu,nu) before adding it to T tensor
!      If (ViscousC>1D-6) Then
!        call regulatePi(Time,NX0,NY0,NZ0,NX,NY,NZ,
!     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
!     &  Ed,PL,PPI,
!     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33)
!      End If
!-------End of regulation---------
!       ---Zhi-End---

! Ver 1.6.29RC5: NXPhy0->NX0, etc

       DO 100 K=NZ0,NZ
       DO 100 J=NY0,NY
       DO 100 I=NX0,NX

        T00M=TT00(I,J,K)/Time-Pi00(I,J,K)   ! ---Zhi-Changes---
        T01M=TT01(I,J,K)/Time-Pi01(I,J,K)   ! ---Zhi-Changes---
        T02M=TT02(I,J,K)/Time-Pi02(I,J,K)   ! ---Zhi-Changes---
        PLM=PL(I,J,K)
        EdM = Ed(I,J,K)
        BulkPi = PPI(I,J,K)

        !*** ABORT SMALL NUMBERS TO AVOID UNDERFLOW **
        !CALL FYSCOR(T00M,T01M,T02M,CC)  !T00=AMIN1(T00,0.0)and  other treatment ???????
        !call FYSCOR_withBulkvis(T00M, T01M, T02M, BulkPi)
        call FYSCOR_velocity(T00M, T01M, T02M, EdM, PLM, BulkPi)

        if(T00M .gt. AAC) then
          T00(I,J,K) = T00M
          PPI(I,J,K) = BulkPi
          IF (ABS(T01M).GT.AAC) THEN
            T01(I,J,K)=T01M
          ELSE
            T01(I,J,K)=0.0
            Pi01(I,J,k)=0.0   !song viscous
          END IF
          IF (ABS(T02M).GT.AAC) THEN
            T02(I,J,K)=T02M
          ELSE
            T02(I,J,K)=0.0
            Pi02(I,J,K)=0.0
          END IF
        else
          T00(I,J,K) = 1D-30
          T01(I,J,K) = 0.0
          T02(I,J,K) = 0.0
          Pi01(I,J,k)=0.0   !song viscous
          Pi02(I,J,K)=0.0
          PPI(I,J,K) = 0.0
          PL(I,J,K) = 0.0
        endif
 100  CONTINUE

C---------------------------------------------------------------
      DO 300 K=NZ0,NZ
      DO 300 J=NYPhy0,NYPhy
!        EK=T00(0,J,K)
!        B0=0.0
!        B0=BN0(0,J,K)
      DO 310 I=NXPhy0,NXPhy !changed by song

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     The original root finding algorithm (generic one)
        !  T00IJ=T00(I,J,K)
        !  T01IJ=T01(I,J,K)
        !  T02IJ=T02(I,J,K)
        !  BN0IJ=0.0       !BN0(I,J) !with  Bayron current
        !  BulkPr=PPI(I,J,K)
        !
        !  AI=I
        !  AJ=J
        !  AK=K
        !  AE=Ed(I,J,K)
        !CALL NEWTON2(B0,EK,T00IJ,T01IJ,T02IJ,BN0IJ,BulkPr,PNEW,NNEW)  !bisection
        !
        !Ed(I,J,K)=DMAX1(0.0D0,EK)
        !Bd(I,J,K)=DMAX1(0.0D0,B0)
        !
        !ee=Ed(I,J,NZ0)*HbarC        !fm-4 to GeV/fm3  peter unit
        !Ba=Bd(I,J,K)                !unit
        !pp=DMAX1(0.0D0,PEPS(Ba,ee))
        !PL(I,J,K)=pp/HbarC           !GeV/fm3 to fm-4
C--------------------------------------

        !Vx(I,J,K)=T01(I,J,K)/Dmax1((T00(I,J,K)+PL(I,J,K)+Bulkpr),AAC)  !??
        !Vy(I,J,K)=T02(I,J,K)/Dmax1((T00(I,J,K)+PL(I,J,K)+Bulkpr),AAC)  !??



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     The following use quadESolver function to find root, which assumes Bd(I,J,K)=0
        Bd(I,J,K)=0
        T00IJ=T00(I,J,K)
        T01IJ=T01(I,J,K)
        T02IJ=T02(I,J,K)
        DM0 = T00IJ
        DM = sqrt(T01IJ*T01IJ + T02IJ*T02IJ)

        IDebug = 0

        RSDM0 = DM0
        RSDM = DM
        RSPPI = PPI(I,J,K)

        !eeH = findEdHook(1D0)
        !Call invertFunctionD(findEdHook,0D0,5D3,1D-3,ED(I,J,K),0D0,eeH)
        VP_local = findvHook(0.0D0)
        Call invertFunctionD(findvHook, 0.0D0, 1.0D0, 1D-6, 0.0, 0D0, 
     &                       VP_local)
        U0_local = findU0Hook(0.0D0)
        Call invertFunctionD(findU0Hook, 1D0, 5D3, 1D-6, 1.0, 0D0, 
     &                       U0_local)
        U0_critial = 1.21061
        if(U0_local .gt. U0_critial) then
          VP_local = sqrt(1. - 1./U0_local/U0_local)
        else
          U0_local = 1./sqrt(1. - VP_local*VP_local)
        endif
         
        eeH = DM0 - VP_local*DM

        !eeH = quadESolver(ED(I,J,K),DM0,DM,PPI(I,J,K),IDebug,I,J) ! use Ed from last time step as a starting point
        ED(I,J,K) = dmax1(eeH, 1D-10)
        DM = dmax1(DM, 1D-10)
        !VP = (DM0-eeH)/DM
        VP = VP_local

        ee=Ed(I,J,NZ0)*HbarC        !fm-4 to GeV/fm3  peter unit
        Ba=Bd(I,J,K)                !unit
        pp=DMAX1(0.0D0,PEPS(Ba,ee))
        PL(I,J,K)=pp/HbarC           !GeV/fm3 to fm-4

        Vx(I,J,K) = VP*T01IJ/DM
        Vy(I,J,K) = VP*T02IJ/DM
        if (isnan(Vx(I,J,K))) then
          print*, 'vx', VP_local, U0_local, VP, T01IJ, DM
          stop
        endif
        if (isnan(Vy(I,J,K))) then
          print*, 'vy', VP_local, U0_local, VP, T02IJ, DM
          stop
        endif
!----------------------------------------------------------------------
 310   CONTINUE
 300   CONTINUE


!---------- End of root finding ------------------------------------------------------
! CAREFUL

      TT00=(T00+Pi00)*Time ! ---Zhi-Changes---
      TT01=(T01+Pi01)*Time ! ---Zhi-Changes---
      TT02=(T02+Pi02)*Time ! ---Zhi-Changes---


!=============== Extrapolation Process ======================================

      DO K=NZ0,NZ
      DO J=NYPhy0,NYPhy
      DO I=NXPhy0,NXPhy
        U0(I,J,K)=1D0/sqrt(1D0-Vx(I,J,K)**2-Vy(I,J,k)**2+1D-10) ! also calculate U0 to save another separated loop
        U1(I,J,K)=U0(I,J,K)*Vx(I,J,K)
        U2(I,J,K)=U0(I,J,K)*Vy(I,J,K)
      End Do
      End Do
      End Do

      !Open(3773,FILE='movie/Vx.dat',STATUS='old',ACCESS='APPEND')
      !Open(3774,FILE='movie/U1.dat',STATUS='old',ACCESS='APPEND')
      !Open(3775,FILE='movie/Vy.dat',STATUS='old',ACCESS='APPEND')
      !Open(3776,FILE='movie/U2.dat',STATUS='old',ACCESS='APPEND')
      !Do I=NXPhy0, NXPhy
      !Do J=NYPhy0, NYPhy
      !Do K=NZ0, NZ
      !  Write(3773,'(f25.8)',ADVANCE='NO') Vx(I,J,K)
      !  Write(3774,'(f25.8)',ADVANCE='NO') U1(I,J,K)
      !  Write(3775,'(f25.8)',ADVANCE='NO') Vy(I,J,K)
      !  Write(3776,'(f25.8)',ADVANCE='NO') U2(I,J,K)
      !End Do
      !End Do
      !  Write(3773,*)
      !  Write(3774,*)
      !  Write(3775,*)
      !  Write(3776,*)
      !End Do
      !Close(3773)
      !Close(3774)
      !Close(3775)
      !Close(3776)

      call VSBdary3(U1,U2,NX0,NY0,NZ0,NX,NY,NZ, ! smearing and extrapolation
     & NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)

      DO K=NZ0,NZ
      DO J=NYPhy0-3,NYPhy+3
      DO I=NXPhy0-3,NXPhy+3
        U0(I,J,K)=sqrt(1D0+U1(I,J,K)**2+U2(I,J,k)**2) ! also calculate U0 to save another separated loop
        Vx(I,J,K)=U1(I,J,K)/U0(I,J,K)
        Vy(I,J,K)=U2(I,J,K)/U0(I,J,K)
        If (Vx(I,J,K)**2+Vy(I,J,K)**2>1D0) Then
          Print*, "Error:"
          Print*, "I,J=",I,J
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "norm=",Vx(I,J,K)**2+Vy(I,J,K)**2
          Print*, "U0,U1,U2=",U0(I,J,K),U1(I,J,K),U2(I,J,K)
          Stop
        End If
      End Do
      End Do
      End Do

!      call VSBdary3(Vx,Vy,NX0,NY0,NZ0,NX,NY,NZ, ! smearing and extrapolation
!     & NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)
      !Do K=NZ0,NZ
      !Do I=NXPhy0-3,NXPhy+3
      !Do J=NYPhy0-3,NYPhy+3
      !  dlength = sqrt(Vx(I,J,K)*Vx(I,J,K) + Vy(I,J,K)*Vy(I,J,K))
      !  If (dlength>1D0) Then
      !    Vx(I,J,K) = Vx(I,J,K)/dlength*0.99999999999D0
      !    Vy(I,J,K) = Vy(I,J,K)/dlength*0.99999999999D0
      !  End If
      !  U0(I,J,K)=1D0/max(sqrt(1D0-Vx(I,J,K)**2-Vy(I,J,k)**2),1D-30) ! also calculate U0 to save another separated loop
      !  U1(I,J,K)=U0(I,J,K)*Vx(I,J,K)
      !  U2(I,J,K)=U0(I,J,K)*Vy(I,J,K)
      !End Do
      !End Do
      !End Do

!=======================================================================

       call TriSembdary3(Bd,PL, Ed,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

      If (ViscousC>1D-6) Then

       call TriSembdary3(Pi00,Pi01, Pi02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(Pi11,Pi22, Pi33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(Pi12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)

      End If

       call Sembdary3(PPI, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)


C-----------------------------------------------------------------
      Call CoefASM2d(CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &               NX0,NY0,NZ0, NX,NY,NZ)
c------------------------------------------------------------------

      If (ViscousC>1D-6) Then

       DO 670 K=NZ0,NZ
       DO 670 J=NYPhy0,NYPhy
       DO 670 I=NXPhy0,NXPhy

        SDX1=(0.0 + Pi11(I+1,J,K) - Pi11(I-1,J,K)
     &   -Vx(I+1,J,K)*Pi01(I+1,J,K)
     &   +Vx(I-1,J,K)*Pi01(I-1,J,K))/(2.0*DX)
        SDY1=(Pi12(I,J+1,K)-Pi12(I,J-1,K)
     &   -Vy(I,J+1,K)*Pi01(I,J+1,K)
     &   +Vy(I,J-1,K)*Pi01(I,J-1,K))/(2.0*DY)
        ScT01(i,j,K)=(SDX1+SDY1)*Time  !ScT01

        SDX2=(Pi12(I+1,J,K)-Pi12(I-1,J,K)
     &   -Vx(I+1,J,K)*Pi02(I+1,J,K)
     &   +Vx(I-1,J,K)*Pi02(I-1,J,K))/(2.0*DX)
        SDY2=(0.0 + Pi22(I,J+1,K)-Pi22(I,J-1,K)
     &   -Vy(I,J+1,K)*Pi02(I,J+1,K)
     &   +Vy(I,J-1,K)*Pi02(I,J-1,K))/(2.0*DY)
        ScT02(i,j,K)=(SDX2+SDY2)*Time  !ScT02

        SDX0=(Pi01(I+1,J,K)-Pi01(I-1,J,K)+Vx(I+1,J,K)*(0.0
     & -Pi00(I+1,J,K))-Vx(I-1,J,K)*(0.0-Pi00(I-1,J,K)))/(2.0*DX)
        SDY0=(Pi02(I,J+1,K)-Pi02(I,J-1,K)+Vy(I,J+1,K)*(0.0
     & -Pi00(I,J+1,K))-Vy(I,J-1,K)*(0.0-Pi00(I,J-1,K)))/(2.0*DY)
        ScT00(i,j,K)=0.0+Time*(SDX0+SDY0)+0.0 !ScT00

 670   continue
      Else
        ScT01 = 0D0
        ScT02 = 0D0
        ScT00 = 0D0
      End If


      If (ViscousC>1D-6) Then

       call TriSembdary3(ScT00,ScT01,ScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ,NXPhy0,NYPhy0, NXPhy,NYPhy)
      End If


       DO 690 K=NZ0,NZ
       DO 690 J=NYPhy0,NYPhy
       DO 690 I=NXPhy0,NXPhy

        ScT01(i,j,K)=ScT01(i,j,k)
     &    +Time*(PL(I+1,J,K)-PL(I-1,J,K)
     &    +PPi(I+1,J,K)-PPi(I-1,J,K))/(2.0*DX)!ScT01
        ScT02(i,j,K)=ScT02(i,j,K)
     &    +Time*(PL(I,J+1,K)-PL(I,J-1,K)
     &    +PPi(I,J+1,K)-PPi(I,J-1,K))/(2.0*DY)

        ScT00(i,j,K)=ScT00(i,j,K)+PL(i,j,K)+PPi(i,j,K)+Pi33(I,J,K)
     & +Time*(Vx(I+1,J,K)*PL(I+1,J,K)-Vx(I-1,J,K)*PL(I-1,J,K))/(2.0*DX)
     &+Time*(Vx(I+1,J,K)*PPi(I+1,J,K)-Vx(I-1,J,K)*PPi(I-1,J,K))/(2.0*DX)
     & +Time*(Vy(I,J+1,K)*PL(I,J+1,K)-Vy(I,J-1,K)*PL(I,J-1,K))/(2.0*DY)
     &+Time*(Vy(I,J+1,K)*PPi(I,J+1,K)-Vy(I,J-1,K)*PPi(I,J-1,K))/(2.0*DY)
 690   continue
 666   continue

        call TriSembdary3(ScT00,ScT01, ScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

       If(Nsm.gt.0) then
       do MM=1,NSm
       call ASM2D6(ScT00,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(ScT01,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(ScT02,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       end do
       end if

C---------------------------------------------------------------------
        call EntropyTemp3 (Ed,PL, Temp,CMu,Sd,
     &         NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

        call ViscousCoefi8(Ed,Sd,PL,Temp,
     &  VCoefi,VCBeta,VRelaxT,etaTtp, VBulk, VCBeta0,VRelaxT0,XiTtP,
     &  Time,DX,DY,DT,NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
C--------------------

      If (ViscousC>1D-6) Then
        call PiS4U5(PU0,PU1,PU2,PU3,U0,U1,U2,U3, DX,DY,DZ, DT,
     & DPc00,DPc01,DPc02,DPc33, DPc11,DPc22,DPc12, DDU0,DDU1,DDU2,
     & Temp,Temp0,  SiLoc,DLnT,  Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     & NX0,NX,NY0,NY,NZ0,NZ)
      Else ! Ver 1.6.19RC5: set to 0
        DPc00=0D0
        DPc01=0D0
        DPc02=0D0
        DPc11=0D0
        DPc12=0D0
        DPc22=0D0
        DPc33=0D0
      End If


      If (ViscousC>1D-6) Then

       DO 700 K=NZ0,NZ
       DO 700 J=NYPhy0,NYPhy
       DO 700 I=NXPhy0,NXPhy

        Tr0=Pi00(i,j,k)*DDU0(i,j,k)-Pi01(i,j,k)*DDU1(i,j,k)
     &                           -Pi02(i,j,k)*DDU2(i,j,k)   !additonal term to keep transversality of pi
        Tr1=Pi01(i,j,k)*DDU0(i,j,k)-Pi11(i,j,k)*DDU1(i,j,k)
     &                           -Pi12(i,j,k)*DDU2(i,j,k)   !additonal term to keep transversality of pi
        Tr2=Pi02(i,j,k)*DDU0(i,j,k)-Pi12(i,j,k)*DDU1(i,j,k)
     &                           -Pi22(i,j,k)*DDU2(i,j,k)   !additonal term to keep transversality of pi

        TPi00=2.0*Tr0
        TPi01=Tr1+Vx(i,j,k)*Tr0
        TPi02=Tr2+Vy(i,j,k)*Tr0
        TPi11=2.0*Vx(i,j,k)*Tr1
        TPi12=Vx(i,j,k)*Tr2+Vy(i,j,k)*Tr1
        TPi22=2.0*Vy(i,j,k)*Tr2
        TPi33=0.0

        PScT00(i,j,K)=(0.0-TPi00)*(-1.0)
        PScT01(i,j,K)=(0.0-TPi01)*(-1.0)
        PScT02(i,j,K)=(0.0-TPi02)*(-1.0)
        PScT33(i,j,K)=(0.0-TPi33)*(-1.0)

        PScT11(i,j,K)=(0.0-TPi11)*(-1.0)
        PScT12(i,j,K)=(0.0-TPi12)*(-1.0)
        PScT22(i,j,K)=(0.0-TPi22)*(-1.0)
 700   continue
       Else ! Ver 1.6.19RC5: set to 0
        PScT00=0D0
        PScT01=0D0
        PScT02=0D0
        PScT11=0D0
        PScT12=0D0
        PScT22=0D0
        PScT33=0D0
      End If
C--------------------------------------------------------------------------------


      If (ViscousC>1D-6) Then

       call TriSembdary3(PScT00,PScT01, PScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(PScT11,PScT22, PScT33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(PScT12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)
      End If

      DO 710 K=NZ0,NZ
      DO 710 J=NYPhy0,NYPhy
      DO 710 I=NXPhy0,NXPhy
        xx=DX*I
        yy=DY*J
        rr=sqrt(xx**2+yy**2)
        ff=1.0/(Dexp((rr-R0Bdry)/Aeps)+1.0)

      Accu=3.0
       If(abs(Accu-3.0).le.0.00001) then  !3pt formula
          PA=(Vx(I+1,J,K)-Vx(I-1,J,K))/(2*DX)
     &      +(Vy(I,J+1,K)-Vy(I,J-1,K))/(2*DY)
       else if(abs(Accu-5.0).le.0.00001) then  !5pt formula
          PA=( Vx(I+1,J,K)*2.0d0/3.0d0-Vx(I-1,J,K)*2.0d0/3.0d0
     &       -Vx(I+2,J,K)/12.0d0+Vx(I-2,J,K)/12.0d0 )/DX
     &      +( Vy(I,J+1,K)*2.0d0/3.0d0-Vy(I,J-1,K)*2.0d0/3.0d0
     &       -Vy(I,J+2,K)/12.0d0+Vy(I,J-2,K)/12.0d0 )/DY
       else
          Print*, "wrong input for Accu"
       end if

          PT=VRelaxT(I,J,K)/U0(I,J,K)
          PS=2.0*VCoefi(I,J,K)

          PT0=VRelaxT0(I,J,K)/U0(I,J,K)
          PS0=VBulk(I,J,K)

       if(ViscousC.ge.0.00001) then
           D0Ln=(DLog(etaTtp0(I,J,K))-DLog(etaTtp(I,J,K)))/DT
           D1Ln=(DLog(etaTtp(I-1,J,K))-DLog(etaTtp(I+1,J,K)))/(2.0*DX)
           D2Ln=(DLog(etaTtp(I,J-1,K))-DLog(etaTtp(I,J+1,K)))/(2.0*DY)

           ADLnT=-0.5*(SiLoc(I,J,K)/U0(I,J,K)         !high precision version to reduce round-off error
     &           +(D0Ln+ Vx(I,J,K)*D1Ln +Vy(I,J,K)*D2Ln)*ff)

        else
           ADLnT=0.0
        end if

       if(VisBulk.ge.0.00001) then
           dB0=(U0(I,J,K)/XiTtp(i,j,k)-PU0(I,J,K)/XiTtp0(i,j,k))
           dB1=(U1(I+1,J,K)/XiTtp(I+1,j,k)-U1(I-1,J,K)/XiTtp(I-1,j,k))
           dB2=(U2(I,J+1,K)/XiTtp(i,j+1,k)-U2(I,J-1,K)/XiTtp(i,j-1,k))
          Badd2=-0.5*XiTtP(i,j,k)*(dB0/DT+dB1/(2.0*DX)+dB2/(2.0*DY)
     &       +U0(I,J,K)/XiTtP(i,j,k)/time)/U0(I,J,K) *ff

          DB0Ln=(DLog(XiTtp0(I,J,K))-DLog(XiTtp(I,J,K)))/DT
          DB1Ln=(DLog(XiTtp(I-1,J,K))-DLog(XiTtp(I+1,J,K)))/(2.0*DX)
          DB2Ln=(DLog(XiTtp(I,J-1,K))-DLog(XiTtp(I,J+1,K)))/(2.0*DX)

           Badd=-0.5*(SiLoc(I,J,K)/U0(I,J,K)    !high precision version to reduce round-off error
     &          +(DB0Ln +Vx(I,J,K)*DB1Ln +Vy(I,J,K)*DB2Ln)*ff)

       else
           Badd=0.0
       end if

          PScT00(i,j,K)=PScT00(i,j,K) +( Pi00(I,J,K)*ADLnT+
     &       PA*Pi00(I,J,K)-PT*(Pi00(I,J,K)-PS*DPc00(i,j,K)))*(-1.0)
          PScT01(i,j,K)=PScT01(i,j,K) +( Pi01(I,J,K)*ADLnT+
     &       PA*Pi01(I,J,K)-PT*(Pi01(I,J,K)-PS*DPc01(i,j,K)))*(-1.0)
          PScT02(i,j,K)=PScT02(i,j,K) +( Pi02(I,J,K)*ADLnT+
     &       PA*Pi02(I,J,K)-PT*(Pi02(I,J,K)-PS*DPc02(i,j,K)))*(-1.0)
          PScT33(i,j,K)=PScT33(i,j,K) +( Pi33(I,J,K)*ADLnT+
     &       PA*Pi33(I,J,K)-PT*(Pi33(I,J,K)-PS*DPc33(i,j,K)))*(-1.0)
          PScT11(i,j,K)=PScT11(i,j,K) +( Pi11(I,J,K)*ADLnT+
     &       PA*Pi11(I,J,K)-PT*(Pi11(I,J,K)-PS*DPc11(i,j,K)))*(-1.0)
          PScT22(i,j,K)=PScT22(i,j,K) +( Pi22(I,J,K)*ADLnT+
     &       PA*Pi22(I,J,K)-PT*(Pi22(I,J,K)-PS*DPc22(i,j,K)))*(-1.0)
          PScT12(i,j,K)=PScT12(i,j,K) +( Pi12(I,J,K)*ADLnT+
     &       PA*Pi12(I,J,K)-PT*(Pi12(I,J,K)-PS*DPc12(i,j,K)))*(-1.0)

          PScTSum = PScT00(i,j,k) + PScT01(i,j,k) + PScT02(i,j,k)
     &      + PScT11(i,j,k) + PScT12(i,j,k) + PScT22(i,j,k)
     &      + PScT33(i,j,k)

      If (TIME > TT0) Then
      If (PScTSum .ne. PScTSum) Then
        Print *, "Invalid PScT terms!"
        Print *, "i,j,k=", i,j,k
        Print *, "PScT00=", PScT00(i,j,K)
        Print *, "PScT01=", PScT01(i,j,K)
        Print *, "PScT02=", PScT02(i,j,K)
        Print *, "PScT11=", PScT11(i,j,K)
        Print *, "PScT12=", PScT12(i,j,K)
        Print *, "PScT22=", PScT22(i,j,K)
        Print *, "PScT33=", PScT33(i,j,K)
        Print *, "Pi00,ADLnT=", Pi00(I,J,K), ADLnT
        Print *, "PA", PA
        Print *, "Vx(I+1,J,K)", Vx(I+1,J,K), "Vx(I-1,J,K)", Vx(I-1,J,K)
        Print *, "Vy(I,J+1,K)", Vy(I,J+1,K), "Vy(I,J-1,K)", Vy(I,J-1,K)
        Print *, "PT=", PT, "Pi00=", Pi00(I,J,K)
        Print *, "PS=", PS, "DPc00=", DPc00(i,j,K)
        Print *, "SiLoc=", SiLoc(I,J,K), "U0=", U0(I,J,K)   !high precision version to reduce round-off error
        Print *, "D0Ln=", D0Ln, "Vx=", Vx(I,J,K), "Vy=", Vy(I,J,K)
        Print *, "D1Ln=", D1Ln, "D2Ln=", D2Ln, "ff=", ff
        Print *, "etaTtp(I,J-1,K)=", etaTtp(I,J-1,K)
        Print *, "etaTtp(I,J+1,K)=", etaTtp(I,J+1,K)

        Stop
      EndIf
      EndIf

          PISc(i,j,K)=0.0 +( PPI(I,J,K)*BAdd+
     &       PA*PPI(I,J,K)-PT0*(PPI(I,J,K)+PS0*SiLoc(i,j,K)))*(-1.0)

 710   continue

!-----------------------------------------------------------------------
!    Check if pi(mu,nu) goes away from sigma(mu,nu) (DPc)
      !dmaxRatio = 0D0
      !Iwhich = -1
      !Imax = -1000
      !Jmax = -1000
      !Kmax = -1000
      !
      !Do I=NXPhy0-3,NXPhy
      !Do J=NYPhy0-3,NYPhy
      !Do K=NZ0,NZ
      !
      !ratio = Pi00(I,J,K)/max(abs(DPc00(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 0
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !ratio = Pi01(I,J,K)/max(abs(DPc01(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 1
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !ratio = Pi02(I,J,K)/max(abs(DPc02(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 2
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !ratio = Pi11(I,J,K)/max(abs(DPc11(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 11
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !ratio = Pi12(I,J,K)/max(abs(DPc12(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 12
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !ratio = Pi22(I,J,K)/max(abs(DPc22(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 22
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !ratio = Pi33(I,J,K)/max(abs(DPc33(I,J,K)), 1D-30)
      !If (abs(ratio)>dmaxRatio) Then
      !  dmaxRatio = ratio
      !  Iwhich = 33
      !  Imax = I
      !  Jmax = J
      !  Kmax = K
      !End If
      !
      !End Do
      !End Do
      !End Do


!       If (Imax.ne.-1000 .OR. Jmax.ne.-1000) Then
!        Print*, "dmaxRatio=", dmaxRatio
!        Print*, "Iwhich=", Iwhich
!        Print*, "Imax,Jmax", Imax,Jmax
!        call printPi(id, Imax, Jmax,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!        Print*, "DPc00=", DPc00(Imax,Jmax,Kmax)
!        Print*, "DPc01=", DPc01(Imax,Jmax,Kmax)
!        Print*, "DPc02=", DPc02(Imax,Jmax,Kmax)
!        Print*, "DPc11=", DPc11(Imax,Jmax,Kmax)
!        Print*, "DPc12=", DPc12(Imax,Jmax,Kmax)
!        Print*, "DPc22=", DPc22(Imax,Jmax,Kmax)
!        Print*, "DPc33=", DPc33(Imax,Jmax,Kmax)
!      End If

!      Open(3771,FILE='movie/DPc.dat',STATUS='old',ACCESS='APPEND')
!      Do I=NXPhy0, NXPhy
!      Do J=NYPhy0, NYPhy
!      Do K=NZ0, NZ
!        Write(3771,'(f25.8)',ADVANCE='NO') abs(DPc00(I,J,K))
!     &    +abs(DPc01(I,J,K)) + abs(DPc02(I,J,K))
!     &    +abs(DPc11(I,J,K)) + abs(DPc12(I,J,K))
!     &    +abs(DPc22(I,J,K)) + abs(DPc33(I,J,K))
!      End Do
!      End Do
!        Write(3771,*)
!      End Do
!      Close(3771)
!
!      I0=13
!      J0=-77
!
!      I=I0-1
!      J=J0-1
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0-1
!      J=J0
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0-1
!      J=J0+1
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0
!      J=J0-1
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0
!      J=J0
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0
!      J=J0+1
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0+1
!      J=J0-1
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0+1
!      J=J0
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
!
!      I=I0+1
!      J=J0+1
!      K=NZ0
!      Print*, "At I=",I, "J=",J
!      call printPi(id, I, J,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
!      call printDPc(id, I, J, Sd,
!     &    NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
!     &    DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)












!-----------------------------------------------------------------------

! 777   continue

       call TriSembdary3(PScT00,PScT01, PScT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call TriSembdary3(PScT11,PScT22, PScT33,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
       call Sembdary3(PScT12, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)
       call Sembdary3(PISc, NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)

       If(Nsm.gt.0) then
       do MM=1,NSm
       call ASM2D6(PScT00,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT01,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT02,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT11,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT12,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT22,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                         NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PScT33,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
       call ASM2D6(PISc,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)

       end do
       end if

C*********************************************************************
       call TriSembdary3(T00, T01, T02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

       call TriSembdary3(TT00, TT01, TT02,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)

C-------------------------------------------------------------------------------------


       If(Nsm.gt.0) then
        do MM=1,Nsm
        Call ASM2D6bsm(TT00,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Call ASM2D6bsm(TT01,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Call ASM2D6bsm(TT02,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        end do
       end if

         If(Time.lt.TT0+dT) then
            StotalSv=0.0
            StotalBv=0.0
            Stotal=0.0
               K=1
            do 111 J=NYPhy0,NYPhy
            do 111 I=NXPhy0,NXPhy
             if(Ed(I-1,J-1,K).ge. Edec1/HbarC) Stotal=Stotal+Sd(i,j,K)
 111        continue
            Stotal=Stotal*time*dx*dy
            StotalSv=Stotal
            StotalBv=Stotal
         else
              DsP1=0.0
              DsP2=0.0
              K=1
            do 222 J=NYPhy0,NYPhy
            do 222 I=NXPhy0,NXPhy
            if(Ed(I-1,J-1,K).ge. Edec1/HbarC)  then
             if(ViscousC>1D-6) then
              Dsp1=Dsp1+(Pi00(I,J,K)**2 +Pi11(I,J,K)**2 +Pi22(I,J,K)**2
     &             +Pi33(I,J,K)**2 +2.0*Pi01(I,J,K)**2
     &             +2.0*Pi02(I,J,K)**2 +2.0*Pi12(I,J,K)**2)/
     &                  (2.0*VCoefi(I,J,K)*Temp(I,J,K))
             end if
             if(Visbulk.ge.0.00001.and.VBulk(I,J,K).ge.0.000001)then
             Dsp2=Dsp2+(PPI(I,J,K)*PPI(I,J,K))/
     &                    (VBulk(I,J,K)*Temp(I,J,K))
             end if
            end if
 222       continue
             StotalSv=StotalSv+ Dsp1*time*dx*dy*dt
             StotalBv=StotalBv+ Dsp2*time*dx*dy*dt

             Stotal=Stotal+(Dsp1+Dsp2)*time*dx*dy*dt
         end if

C            Print *, 'time',time,'Stotal', Stotal,StotalSv,StotalBv


         If(Time.le.TT0+dT) SxyT=0.0d0

        call anisotrop10 (Vx,Vy,U0,U1,U2, Ed,PL,Sd, PPI,
     &  Pi00,Pi01,Pi02, Pi33, Pi11,Pi12,Pi22, ScT00,ScT01,ScT02,
     &  DX,DY,DZ,DT,Time,EpsX,EpsP,TEpsP, Vaver,VavX,VavY, STotal22,
     &  AE,ATemp, APL, APLdx,APLdy, AScT00, AScT01, AScT02,
     &  AP11dx,AP22dy,  AP33,AP00, AP01,AP02, AP11, Ap12, AP22,SxyT,
     &  PaScT00,PaScT01,PaScT02, difScT0102,difPLxy, difPixy, difVxVy,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy, NX0,NX,NY0,NY,NZ0,NZ)  !PNEW NNEW something related to root finding



      AMV=Hbarc*1000.0 !fm-1 change to MEV
      GV=Hbarc !fm-1 change to GEV

      Write(92,'(2f8.3,9e14.4)') Time-TT0, VAver,
     &   EpsX,EpsP,TEpsP,Stotal

      Call snapshotEccentricity(Time-TT0,Vx,Vy,Ed,9,DX,DY,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ)

      Write(91,'(f8.3,5e14.4 )') Time-TT0,ATemp,APL,
     &   Temp(0,0,NZ0)*AMV, Temp(0,90,NZ0)*AMV, Temp(90,0,NZ0)*AMV

      Write(90,'(f8.3,7e12.3)') Time-TT0,
     &    AP33*GV,AP00*GV,AP01*GV,AP02*GV,AP11*GV,AP22*GV,AP12*GV

      Write(89,'(f8.3,10e14.4)') Time-TT0,-AScT00*GV,-AScT01*GV,
     &  -AScT02*GV, -PaScT00*GV, -PaScT01*GV, -PaScT02*GV

      Write(88,'(f8.3,10e14.4)') Time-TT0,-difScT0102*GV ,-difPLxy*GV,
     &  -difPixy*GV, -APL*GV, -APLdx*GV, -APLdy*GV,
     &   AP11dx*GV, AP22dy*GV


!     Checking codes deleted
!     See previous versions
!     Search "goto 888"

      Return
      End








C###########################################################################################

      subroutine anisotrop10 (Vx,Vy,U0,U1,U2, Ed,PL,Sd, PPI,
     &  Pi00,Pi01,Pi02, Pi33, Pi11,Pi12,Pi22, ScT00,ScT01,ScT02,
     &  DX,DY,DZ,DT,Time,EpsX,EpsP,TEpsP, Vaver,VavX,VavY, STotal,
     &  AE,ATemp, APL, APLdx,APLdy, AScT00, AScT01, AScT02,
     &  AP11dx,AP22dy,  AP33,AP00, AP01,AP02, AP11, Ap12, AP22,SxyT,
     &  PaScT00,PaScT01,PaScT02, difScT0102,difPLxy, difPixy, difVxVy,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy, NX0,NX,NY0,NY,NZ0,NZ)  !PNEW NNEW  related to root finding

      Implicit Double Precision (A-H, O-Z)
      Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ) !fluid velocity
      Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ) !fluid velocity

      Dimension U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Dimension U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Dimension Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Dimension PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Dimension Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
      Dimension Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature

      Dimension Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Dimension Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
       Dimension PPI(NX0:NX, NY0:NY, NZ0:NZ)    !Bulk pressure


      Dimension ScT00(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ
      Dimension ScT01(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ
      Dimension ScT02(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ

      common/Edec/Edec    !decoupling temperature

       EpsX1=0.0  !relate spacial ellipticity
       EpsX2=0.0  !relate spacial ellipticity
       EpsP1=0.0  !relate momentum ellipticity
       EpsP2=0.0  !relate momentum ellipticity
       TEpsP1=0.0  !relate momentum ellipticity
       TEpsP2=0.0  !relate momentum ellipticity
       STotal=0.0
       Vaver1=0.0    !average velocity
       Vaver2=0.0
       VavX1=0.0
       VavY1=0.0

      DO 100 K=NZ0,NZ
      DO 100 J=-NYPhy,NYPhy
      DO 100 I=-NXPhy,NXPhy
        !If (Ed(I,J,K) < EDec) Cycle
        gamma=1.0D0/sqrt(1D0-Vx(I,J,NZ0)**2-Vy(I,J,NZ0)**2)
        !gamma=1.0D0/max(sqrt(abs(1D0-Vx(I,J,K)**2-Vy(I,J,K)**2)),1D-15)

        xx=I*DX
        yy=J*DY

        x2=xx*xx
        y2=yy*yy

        EpsX1=EpsX1+(y2-x2)*Ed(I,J,K)*gamma   !*
        EpsX2=EpsX2+(y2+x2)*Ed(I,J,K)*gamma   !*

        !angle = 2*atan2(YY,XX)
        !RR = sqrt(XX*XX+YY*YY)
        !
        !gammaE=1.0/sqrt(1-Vx(I,J,K)**2-Vy(I,J,K)**2)*Ed(I,J,K)

        ! eccentricity, defined using r^2 as weight function:
        !Weight=Weight+RR*RR*gammaE ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
        !EpsX1=EpsX1+RR*RR*cos(angle)*gammaE
        !EpsX2=EpsX2+RR*RR*gammaE

        ep=Ed(I,J,K)+PL(I,J,K)
        TXX=ep*U1(I,J,K)*U1(I,J,K)+PL(I,J,K)
        TYY=ep*U2(I,J,K)*U2(I,J,K)+PL(I,J,K)

        EpsP1=EpsP1+(TXX-TYY)
        EpsP2=EpsP2+(TXX+TYY)

        TTXX=TXX+Pi11(I,J,K)+PPI(I,J,K)
     &   +PPI(I,J,K)*U1(I,J,K)*U1(I,J,K)
        TTYY=TYY+Pi22(I,J,K)+PPI(I,J,K)
     &   +PPI(I,J,K)*U2(I,J,K)*U2(I,J,K)

        TEpsP1=TEpsP1+(TTXX-TTYY)
        TEpsP2=TEpsP2+(TTXX+TTYY)

        Vr=sqrt(Vx(I,J,K)**2+Vy(I,J,K)**2)
        Vaver1=Vaver1+Ed(I,J,K)*Vr*gamma   !*
        VavX1=VavX1+Ed(I,J,K)*Vx(I,J,K)*gamma   !*
        VavY1=VavY1+Ed(I,J,K)*Vy(I,J,K)*gamma   !*
        Vaver2=Vaver2+Ed(I,J,K)*gamma   !*

        STotal=STotal+Sd(I,J,K)*U0(I,J,K)*Time*dx*dy
 100  continue

        EpsX=EpsX1/EpsX2  !spacial ellipticity
        EpsP=EpsP1/EpsP2  !Momentum  ellipticity
        TEpsP=TEpsP1/TEpsP2 ! total Momentum  ellipticity

        Vaver=Vaver1/Vaver2
        VavX=VaVX1/Vaver2
        VavY=VaVY1/Vaver2

      K=NZ0
      Sx=0.0d0
      I=NXPhy
      DO 110 J=NYPhy0,NYPhy
       Sx=Sx+1.0*Time*dy*dT*U1(I,J,K)*Sd(I,J,K)
 110  continue


      K=NZ0
      Sy=0.0d0
      J=NYPhy
      DO 120 I=NXPhy0,NXPhy
       Sy=Sy+1.0*Time*dx*dT*U2(I,J,K)*Sd(I,J,K)
 120  continue


         SxyT=SxyT+Sx+Sy
         STotal=(STotal+SxyT)*4.0
C        STotal=STotal*Time*dx*dy*DT


       AE=0.0
       ATemp=0.0
       AP00=0.0
       AP01=0.0
       AP02=0.0
       AP33=0.0
       AP11=0.0
       AP12=0.0
       AP22=0.0

       APLdx=0.0
       APLdy=0.0
       APL=0.0

       PaScT00=0.0
       PaSCT01=0.0
       PaSCT02=0.0

       AScT00=0.0
       ASCT01=0.0
       ASCT02=0.0


       difScT0102=0.0
       difPLxy=0.0
       difPixy=0.0
       difVxVy=0.0
C----------------------------------------------------------------
      DO 200 K=NZ0,NZ
      DO 200 J=0,NYPhy
      DO 200 I=0,NXPhy
      AE=AE+Ed(I,J,K)*gamma  !*
      ppe=PL(I,J,K)+Ed(I,J,K)
      ATemp=ATemp+Ed(I,J,K)*Temp(I,J,K)*gamma  !*

      AP00=AP00+Ed(I,J,K)*Pi00(I,J,K)*gamma/ppe  !*
      AP01=AP01+Ed(I,J,K)*Pi01(I,J,K)*gamma/ppe  !*
      AP02=AP02+Ed(I,J,K)*Pi02(I,J,K)*gamma/ppe  !*
      AP11=AP11+Ed(I,J,K)*Pi11(I,J,K)*gamma/ppe  !*
      AP12=AP12+Ed(I,J,K)*Pi12(I,J,K)*gamma/ppe  !*
      AP22=AP22+Ed(I,J,K)*Pi22(I,J,K)*gamma/ppe  !*
      AP33=AP33+Ed(I,J,K)*Pi33(I,J,K)*gamma/ppe  !*

      APL=APL+PL(I,J,K)*Ed(I,J,K)*gamma  !*
     & +Time*Ed(I,J,K)*gamma  !*
     & *((PL(I+1,J,K)*Vx(I+1,J,K) -PL(I-1,J,K)*Vx(I-1,J,K))/(2.0*DX)
     &    +(PL(I,J+1,K)*Vy(I,J+1,K)-PL(I,J-1,K)*Vy(I,J-1,K))/(2.0*DY))

      APLdx=APLdx+(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
      APLdy=APLdy+(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*

      AP11dx=AP11dx+(Pi11(I+1,J,K)-Pi11(I-1,J,K))
     &                                  *Ed(I,J,K)/(2.0*DX) *gamma  !*
      AP22dy=AP22dy+(Pi22(I,J+1,K)-Pi22(I,J-1,K))
     &                                  *Ed(I,J,K)/(2.0*DY) *gamma  !*

      AScT00=ASCT00+ScT00(I,J,K)*Ed(I,J,K)*gamma  !*
      AScT01=ASCT01+ScT01(I,J,K)*Ed(I,J,K)*gamma  !*
      AScT02=ASCT02+ScT02(I,J,K)*Ed(I,J,K)*gamma  !*

      difScT0102=difScT0102+(ScT01(I,J,K)-ScT02(I,J,K))*Ed(I,J,K)*gamma  !*
      difPLxy=difPLxy
     &     +Time*(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &     -Time*(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
      difPixy=difPixy
     &      +Time*(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &      +Time*(Pi11(I+1,J,K)-Pi11(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &      -Time*(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
     &      -Time*(Pi22(I,J+1,K)-Pi22(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*

      difVxVy=difVxVy+(Vx(I,J,K)-Vy(I,J,K))*Ed(I,J,K) *gamma  !*

      PaScT00=PaScT00+PL(I,J,K)*Ed(I,J,K)*gamma  !*
     &  +Time*Ed(I,J,K)*gamma  !*
     &  *((PL(I+1,J,K)*Vx(I+1,J,K) -PL(I-1,J,K)*Vx(I-1,J,K))/(2.0*DX)
     &    +(PL(I,J+1,K)*Vy(I,J+1,K)-PL(I,J-1,K)*Vy(I,J-1,K))/(2.0*DY))
     & +Ed(I,J,K)*Pi33(I,J,K)*gamma  !*

      PaSCT01=PaScT01
     &   +(PL(I+1,J,K)-PL(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
     &   +(Pi11(I+1,J,K)-Pi11(I-1,J,K))*Ed(I,J,K)/(2.0*DX)*gamma  !*
      PaSCT02=PaScT02
     &    +(PL(I,J+1,K)-PL(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
     &    +(Pi22(I,J+1,K)-Pi22(I,J-1,K))*Ed(I,J,K)/(2.0*DY)*gamma  !*
 200  continue

      ATemp=ATemp/AE

      AP00=AP00/AE
      AP01=AP01/AE
      AP02=AP02/AE
      AP11=AP11/AE
      AP12=AP12/AE
      AP22=AP22/AE
      AP33=AP33/AE

      APL=APL/AE
      APLdx=Time*APLdx/AE
      APLdy=Time*APLdy/AE
      AP11dx=Time*AP11dx/AE
      AP22dy=Time*AP22dy/AE

      AScT00=AScT00/AE
      AScT01=AScT01/AE
      AScT02=AScT02/AE

      PaScT00=PaScT00/AE
      PaScT01=Time*PaScT01/AE
      PaScT02=Time*PaScT02/AE

      difScT0102=difScT0102/AE
      difPLxy=difPLxy/AE
      difPixy=difPixy/AE
      difVxVy=difVxVy/AE

      return
      end


!======================================================================
      Double Precision Function quadESolver(ee0,DM0,DM,PPI,debug,I,J)
      ! Root finding by interating quadratic equation formula for the equation
      ! (M0+cs^2*e+Pi)(M0-e)=M^2

      Implicit None

      interface
        function PEPS(x,y)
        Double Precision :: x, y
        Double Precision :: PEPS
        end function PEPS
      end interface

      !Double Precision quadESolver

      Integer I,J
      Double Precision ED

      Double Precision ee, ee0, DM0, DM, PPI ! energy density, initial energy density, M0, M, Pi (see 0510014)
      Double Precision accuracy

      Double Precision ee1, cs2, pp ! energy denstiy from previous interation, p/e, presure
      ! Note that cs2 is not dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      Double Precision, Parameter :: HbarC=0.19733d0 !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm

      Integer, Parameter :: maxIter=300 ! maximum number of iterations
      Integer numIter ! number of iterations

      Integer debug

      Double Precision A,B ! intermedia step variables

      if (debug == 1) print *, "quadESolver started in debug moded..."

      ee1 = ee0
      cs2=PEPS(0.0d0, ee0*Hbarc)/dmax1(ee0,zero)/Hbarc ! check dimension?
      A=DM0*(1-cs2)+PPI
      B=DM0*(DM0+PPI)-DM*DM
      ee=2*B/dmax1((sqrt(A*A+4*cs2*B)+A), zero)

      If (debug.eq.1) Print *, "I,J=",I,J
      If (debug.eq.1) Print *, "ee,ee1,p/e=",ee,ee1,cs2

      numIter=1
      accuracy=ee*1e-6 ! relative accuracy: answer should be correct up to e*1e-10
      If (accuracy < 1e-15) accuracy=1e-15 ! No messing with machines precision (32bit double)
      !accuracy=1e-10
      Do While (abs(ee-ee1) > accuracy)
        ee1 = ee
        cs2=PEPS(0.0d0, ee1*Hbarc)/dmax1(ee1,zero)/Hbarc
        A=DM0*(1-cs2)+PPi
        B=DM0*(DM0+PPI)-DM*DM
        ee=2*B/dmax1((sqrt(A*A+4*cs2*B)+A), zero)
        If (debug.eq.1) Print *, "ee,ee1,p/e=",ee,ee1,cs2
        numIter=numIter+1
        If (numIter>maxIter) Then
          Print *, "quadESolver: too many iterations."
          Print *, "Number of iterations=", numIter
          Print *, "I,J=",I,J
          Print *, "M0,M=",DM0,DM
          Print *, "ee,ee1=",ee,ee1
          Print *, "cs2=", cs2
          !Stop
          Exit
        EndIf

      End Do

      quadESolver = ee
      If (debug.eq.1) Print *, "quadESolver=", quadESolver
      if (debug == 1) print *, "quadESolver finished."

      End Function
!----------------------------------------------------------------------






!======================================================================
      Double Precision Function findEdHook(ee)
      ! Root finding by iterating quadratic equation formula for the equation
      ! (M0+cs^2*e+Pi)(M0-e)=M^2

      Implicit None

      !interface
      !  function PEPS(x,y)
      !  Double Precision :: x, y
      !  Double Precision :: PEPS
      !  end function PEPS
      !end interface

      Double Precision PEPS

      Double Precision ee ! energy density

      Double Precision :: RSDM0, RSDM, RSPPI
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double Precision cs2 ! energy density from previous iteration, p/e, pressure
      ! Note that cs2 is not dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      Double Precision, Parameter :: HbarC=0.19733d0 !for changing between fm and GeV ! Hbarc=0.19733=GeV*fm

      Double Precision A,B ! intermedia step variables

      cs2=PEPS(0.0d0, ee*Hbarc)/dmax1(abs(ee),zero)/Hbarc ! check dimension?
      A=RSDM0*(1-cs2)+RSPPI
      B=RSDM0*(RSDM0+RSPPI)-RSDM*RSDM
      findEdHook = ee-2*B/dmax1((sqrt(A*A+4*cs2*B)+A), zero)

      End Function
!----------------------------------------------------------------------

!======================================================================
      Double Precision Function findvHook(v)
      ! Root finding by iterating quadratic equation formula for the equation
      ! (M0+cstilde^2*e+Pi)v=M

      Implicit None

      Double Precision PEPS
      Double Precision v

      Double Precision :: RSDM0, RSDM, RSPPI, RSee
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double Precision cstilde2 ! energy density from previous iteration, p/e, pressure
      ! Note that cstilde2 is NOT dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      Double Precision, Parameter :: HbarC=0.19733d0 !for changing between fm and GeV ! Hbarc=0.19733=GeV*fm

      Double Precision A ! temporary variables

      RSee = RSDM0 - v*RSDM
      cstilde2=PEPS(0.0d0, RSee*Hbarc)/dmax1(abs(RSee),zero)/Hbarc
      A=RSDM0*(1+cstilde2)+RSPPI

      findvHook = v - (2*RSDM)/(A + sqrt(dmax1(A*A 
     &                  - 4*cstilde2*RSDM*RSDM, 1D-30)) + 1D-30)
      if(isnan(findvHook)) then
        print*, cstilde2, A, v
        print*, A*A - 4*cstilde2*RSDM*RSDM
        print*, A + sqrt(A*A - 4*cstilde2*RSDM*RSDM)
        print*, RSDM0, RSDM, RSPPI, RSee
        stop
      endif
      End Function
!----------------------------------------------------------------------

!======================================================================
      Double Precision Function findU0Hook(U0)
      ! Root finding by iterating quadratic equation formula for the equation
      Implicit None

      Double Precision PEPS
      Double Precision v, U0

      Double Precision :: RSDM0, RSDM, RSPPI, RSee
      Common /findEdHookData/ RSDM0, RSDM, RSPPI ! M0, M, Pi (see 0510014)

      Double Precision cstilde2 ! energy density from previous iteration, p/e, pressure
      ! Note that cstilde2 is NOT dp/de, but rather p/e!!!

      Double Precision, Parameter :: zero=1e-30 ! a small number
      Double Precision, Parameter :: HbarC=0.19733d0 !for changing between fm and GeV ! Hbarc=0.19733=GeV*fm

      Double Precision A ! temporary variables

      RSee = RSDM0 - v*RSDM
      cstilde2=PEPS(0.0d0, RSee*Hbarc)/dmax1(abs(RSee),zero)/Hbarc
      A=RSDM0*(1+cstilde2)+RSPPI
      v = (2*RSDM)/(A+sqrt(dmax1(A*A-4*cstilde2*RSDM*RSDM, 1D-30)) 
     &              + 1D-30)

      findU0Hook = U0 - 1/sqrt(1-v*v)
      End Function
!----------------------------------------------------------------------


!======================================================================
      Subroutine snapshotEccentricity(relTime,Vx,Vy,Ed,nmoments,
     &  DX, DY,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ)
      Implicit None

      Integer NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ
      Integer nmoments
      Double Precision relTime, Edec, DX, DY
      Common /EDec/ EDec
      Double Precision Vx(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Vy(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision XN(1:nmoments), YN(1:nmoments), Weight
      Double Precision XNP(1:nmoments), YNP(1:nmoments), WeightP

      Double Precision XX, YY, gammaE, angle, RR, XC, YC, TotalE
      Integer I, J, K, NN

      Open(377,FILE='results/ecc-evolution.dat',
     &          STATUS='old',ACCESS='APPEND')

      Write (377,378,ADVANCE='NO') relTime

!     Calculate where the real center is:
      XC = 0D0
      YC = 0D0
      TotalE = 0D0
      Do K = NZ0, NZ
      Do I = NXPhy0, NXPhy
      Do J = NYPhy0, NYPhy
        !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
        XX = I*DX
        YY = J*DY
        XC = XC + XX*Ed(I,J,K)
        YC = YC + YY*Ed(I,J,K)
        TotalE = TotalE + Ed(I,J,K)
      End Do
      End Do
      End Do
      XC = XC/TotalE
      YC = YC/TotalE

      Do NN = 1, nmoments
        Weight = 0D0
        WeightP = 0D0
        XN(NN) = 0D0
        YN(NN) = 0D0
        XNP(NN) = 0D0
        YNP(NN) = 0D0
        Do K = NZ0, NZ
        Do I = NXPhy0, NXPhy
        Do J = NYPhy0, NYPhy
          !If (Ed(I,J,K)<EDec) Cycle ! Sum only inside the freeze-out surface
          XX = I*DX - XC ! shift to the real center
          YY = J*DY - YC

          angle = NN*atan2(YY,XX)
          RR = sqrt(XX*XX+YY*YY)

          gammaE=1.0/sqrt(1-Vx(I,J,K)**2-Vy(I,J,K)**2)*Ed(I,J,K)

          ! eccentricity, defined using r^2 as weight function:
          Weight=Weight+RR*RR*gammaE ! Note that this weight is repeatedly calculated. But since it is much cleaner written this way and it is very fast...
          XN(NN)=XN(NN)+RR*RR*cos(angle)*gammaE
          YN(NN)=YN(NN)+RR*RR*sin(angle)*gammaE
          ! eccentricity, defined using r^n as weight function:
          WeightP = WeightP+RR**NN*gammaE
          XNP(NN) = XNP(NN)+RR**NN*cos(angle)*gammaE
          YNP(NN) = YNP(NN)+RR**NN*sin(angle)*gammaE
        End Do
        End Do
        End Do
        Write (377,378,ADVANCE='NO')
     &          -XN(NN)/Weight, -YN(NN)/Weight, ! Note that I use minor axis to define eccentricity
     &          sqrt(XN(NN)*XN(NN)+YN(NN)*YN(NN))/Weight,
     &          -XNP(NN)/WeightP, -YNP(NN)/WeightP,
     &          sqrt(XNP(NN)*XNP(NN)+YNP(NN)*YNP(NN))/WeightP
 378    Format(7(E20.8))
      End Do ! NN=1, nmoments

      Close(377)

      End Subroutine
!----------------------------------------------------------------------


      Subroutine checkPi(id, Time, Vx, Vy, Ed, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
      ! Check traceless and transversality condition for pi(mu,nu)

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,id
      Double Precision Time

      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Vx(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Vy(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor


      Integer I,J,K
      Double Precision trace_pi, trans, toUse

      Double Precision :: accuracy=1D0

      If (silent_checkPi==1) Return

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        toUse = max(Ed(I,J,K)*1D0/(1D0-Vx(I,J,K)**2-Vy(I,J,K)**2),1D-1)

        trace_pi = Pi00(I,J,K)-Pi11(I,J,K)-Pi22(I,J,K)-Pi33(I,J,K)
        If (abs(trace_pi)>toUse*accuracy) Then
          Print*, "id=", id
          Print*, "Time=", Time
          Print*, "Pi should be traceless!"
          Print*, "I,J=", I,J
          Print*, "Pi00-Pi11-Pi22-Time**2*Pi33=", trace_pi
          Print*, "Pi00,Pi11,Pi22,Time**2*Pi33=",
     &      Pi00(I,J,K),Pi11(I,J,K),Pi22(I,J,K),Pi33(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Stop
        End If

        trans = Pi01(I,J,K)-Vx(I,J,K)*Pi11(I,J,K)-Vy(I,J,K)*Pi12(I,J,K)
        If (abs(trans)>toUse*accuracy) Then
          Print*, "id=", id
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi01-Vx*Pi11-Vy*Pi12=", trans
          Print*, "Pi01,Pi11,Pi12=",
     &      Pi01(I,J,K),Pi11(I,J,K),Pi12(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Stop
        End If

        trans = Pi02(I,J,K)-Vx(I,J,K)*Pi12(I,J,K)-Vy(I,J,K)*Pi22(I,J,K)
        If (abs(trans)>toUse*accuracy) Then
          Print*, "id=", id
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi02-Vx*Pi12-Vy*Pi22=", trans
          Print*, "Pi02,Pi12,Pi22=",
     &      Pi02(I,J,K),Pi12(I,J,K),Pi22(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Stop
        End If

        trans = Pi00(I,J,K)-Vx(I,J,K)*Pi01(I,J,K)-Vy(I,J,K)*Pi02(I,J,K)
        If (abs(trans)>toUse*accuracy) Then
          Print*, "id=", id
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi00-Vx*Pi01-Vy*Pi02=", trans
          Print*, "Pi00,Pi01,Pi02=",
     &      Pi00(I,J,K),Pi01(I,J,K),Pi02(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Stop
        End If

      End Do
      End Do
      End Do

      If (echo_level>=2) Print *, "checkPi passed; id=", id

      End Subroutine
!-----------------------------------------------------------------------


      Subroutine checkPiAll(failed, II, JJ, Time, Vx, Vy, Ed, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  Pi00Regulated,Pi01Regulated,Pi02Regulated,Pi33Regulated,
     &  Pi11Regulated,Pi12Regulated,Pi22Regulated,say_level)
      ! Check tracelessness + transversality + positiveness of Tr(pi^2)

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ
      Integer failed, II, JJ ! (II,JJ): where it fails
      Double Precision Time

      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Vx(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Vy(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision Pi00Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22Regulated(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Integer say_level

      Integer I,J,K
      Double Precision trace_pi, trans, TrPi2

      Double Precision :: Tideal_scale, pi_scale
      Double Precision :: absNumericalzero = 1D-2
      Double Precision :: relNumericalzero = 1D1  !Xsi_0 in Zhi's thesis

      Double Precision maxPiRatio
      Double Precision gamma_perp
      Common /maxPiRatio/ maxPiRatio

      failed = 0

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        Tideal_scale = sqrt(Ed(I,J,K)**2 + 3*PL(I,J,K)**2)
        
        TrPi2 = Pi00(I,J,K)**2+Pi11(I,J,K)**2+Pi22(I,J,K)**2
     &          +Pi33(I,J,K)**2
     &          -2*Pi01(I,J,K)**2-2*Pi02(I,J,K)**2+2*Pi12(I,J,K)**2

        pi_scale = sqrt(abs(TrPi2))
        !Positivity of Tr(pi^2)
        if(TrPi2 < -relNumericalzero*pi_scale) then
          If (say_level>=9) Then
            Print*, "Time=", Time
            Print*, "Trace Pi^2 is negative!"
            Print*, "I,J=", I,J
            Print*, "Trace pi^2=", TrPi2
            Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
            Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &        Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &        Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &        Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &        Pi12Regulated(I,J,K)
            Print*, "Before: Trace pi^2=",
     &        Pi00Regulated(I,J,K)**2+Pi11Regulated(I,J,K)**2
     &        +Pi22Regulated(I,J,K)**2+Pi33Regulated(I,J,K)**2
     &        -2*Pi01Regulated(I,J,K)**2-2*Pi02Regulated(I,J,K)**2
     &        +2*Pi12Regulated(I,J,K)**2
          End If
          II = I
          JJ = J
          failed = 1
          return
        endif

        !Largeness of pi tensor Tr(pi^2)
        If (pi_scale > max(maxPiRatio*Tideal_scale,
     &      absNumericalzero)) Then
          If (say_level>=9) Then
          Print*, "Time=", Time
          Print*, "pi is too large!"
          Print*, "I,J=", I,J
          Print*, "Trace pi^2=", TrPi2
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Trace pi^2=",
     &      Pi00Regulated(I,J,K)**2+Pi11Regulated(I,J,K)**2
     &      +Pi22Regulated(I,J,K)**2+Pi33Regulated(I,J,K)**2
     &      -2*Pi01Regulated(I,J,K)**2-2*Pi02Regulated(I,J,K)**2
     &      +2*Pi12Regulated(I,J,K)**2
          End If
          II = I
          JJ = J
          failed = 1
          return
        End If

        !trace of pi tensor
        trace_pi = Pi00(I,J,K)-Pi11(I,J,K)-Pi22(I,J,K)-Pi33(I,J,K)
        If(abs(trace_pi) > max(relNumericalzero*pi_scale,
     &     absNumericalzero) ) Then
          If (say_level>=9) Then
          Print*, "Time=", Time
          Print*, "Pi should be traceless!"
          Print*, "I,J=", I,J
          Print*, "Pi00-Pi11-Pi22-Time**2*Pi33=", trace_pi
          Print*, "Pi00,Pi11,Pi22,Time**2*Pi33=",
     &      Pi00(I,J,K),Pi11(I,J,K),Pi22(I,J,K),Pi33(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi00-Pi11-Pi22-Time**2*Pi33=",
     &      Pi00Regulated(I,J,K)-Pi11Regulated(I,J,K)
     &      -Pi22Regulated(I,J,K)-Pi33Regulated(I,J,K)
          End If
          II = I
          JJ = J
          failed = 1
          return
        End If

        !transversality of pi tensor  u_mu pi^{mu,x} = 0
        gamma_perp = sqrt(1./(1. - Vx(I,J,K)**2 - Vy(I,J,K)**2 + 1D-30))
        trans = gamma_perp*(Pi01(I,J,K) - Vx(I,J,K)*Pi11(I,J,K)
     &                      - Vy(I,J,K)*Pi12(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale, 
     &      absNumericalzero)) Then
          If (say_level>=9) Then
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi01-Vx*Pi11-Vy*Pi12=", trans
          Print*, "Pi01,Pi11,Pi12=",
     &      Pi01(I,J,K),Pi11(I,J,K),Pi12(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi01-Vx*Pi11-Vy*Pi12=",
     &      Pi01Regulated(I,J,K)-Vx(I,J,K)*Pi11Regulated(I,J,K)
     &      -Vy(I,J,K)*Pi12Regulated(I,J,K)
          End If
          II = I
          JJ = J
          failed = 1
          return
        End If

        !transversality of pi tensor  u_mu pi^{mu,y} = 0
        trans = gamma_perp*(Pi02(I,J,K) - Vx(I,J,K)*Pi12(I,J,K)
     &                      - Vy(I,J,K)*Pi22(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
          If (say_level>=9) Then
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi02-Vx*Pi12-Vy*Pi22=", trans
          Print*, "Pi02,Pi12,Pi22=",
     &      Pi02(I,J,K),Pi12(I,J,K),Pi22(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi02-Vx*Pi12-Vy*Pi22=",
     &      Pi02Regulated(I,J,K)-Vx(I,J,K)*Pi12Regulated(I,J,K)
     &      -Vy(I,J,K)*Pi22Regulated(I,J,K)
          End If
          II = I
          JJ = J
          failed = 1
          return
        End If

        !transversality of pi tensor  u_mu pi^{mu,tau} = 0
        trans = gamma_perp*(Pi00(I,J,K) - Vx(I,J,K)*Pi01(I,J,K)
     &                      - Vy(I,J,K)*Pi02(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
          If (say_level>=9) Then
          Print*, "Time=", Time
          Print*, "Tranversality violated!"
          Print*, "I,J=", I,J
          Print*, "Pi00-Vx*Pi01-Vy*Pi02=", trans
          Print*, "Pi00,Pi01,Pi02=",
     &      Pi00(I,J,K),Pi01(I,J,K),Pi02(I,J,K)
          Print*, "Vx,Vy=",Vx(I,J,K),Vy(I,J,K)
          Print*, "Ed,PL=", Ed(I,J,K),PL(I,J,K)
          Print*, "Before: Pi00,Pi11,Pi22,Pi33*T^2,Pi01,Pi02,Pi12=",
     &      Pi00Regulated(I,J,K),Pi11Regulated(I,J,K),
     &      Pi22Regulated(I,J,K),Pi33Regulated(I,J,K),
     &      Pi01Regulated(I,J,K),Pi02Regulated(I,J,K),
     &      Pi12Regulated(I,J,K)
          Print*, "Before: Pi00-Vx*Pi01-Vy*Pi02=",
     &      Pi00Regulated(I,J,K)-Vx(I,J,K)*Pi01Regulated(I,J,K)
     &      -Vy(I,J,K)*Pi02Regulated(I,J,K)
          End If
          II = I
          JJ = J
          failed = 1
          return
        End If


      End Do
      End Do
      End Do

      End Subroutine
!-----------------------------------------------------------------------

      Subroutine checkSigma(u0, u1, u2, 
     &           sigma00, sigma01, sigma02, sigma11, sigma12, sigma22,
     &           sigma33, trigger)
      !check the traceless + transversality of the velocity shear tensor
      Implicit None

      Double Precision u0, u1, u2

      Integer I,J,K
      Integer trigger
      Double Precision sigma00, sigma01, sigma02, 
     &                 sigma11, sigma12, sigma22, sigma33
      Double Precision trace, trans0, trans1, trans2
      
      Double Precision :: absNumericalzero = 1D-2
      Double Precision :: relNumericalzero = 1D-2  !Xsi_0 in Zhi's thesis
 
      trigger = 0

      trace = sigma00 - sigma11 - sigma22 - sigma33
      trans0 = u0*sigma00 - u1*sigma01 - u2*sigma02
      trans1 = u0*sigma01 - u1*sigma11 - u2*sigma12
      trans2 = u0*sigma02 - u1*sigma12 - u2*sigma22

      if(abs(trace) .gt. absNumericalzero) then
         trigger = 1
         print*, "velocity shear tensor: violate traceless"
         print*, "trace = ", trace
      endif
      if(abs(trans0) .gt. absNumericalzero) then
         trigger = 1
         print*, "velocity shear tensor: violate transversality"
         print*, "trans0 = ", trans0
      endif
      if(abs(trans1) .gt. absNumericalzero) then
         trigger = 1
         print*, "velocity shear tensor: violate transversality"
         print*, "trans1 = ", trans1
      endif
      if(abs(trans2) .gt. absNumericalzero) then
         trigger = 1
         print*, "velocity shear tensor: violate transversality"
         print*, "trans2 = ", trans2
      endif

      if(trigger .eq. 1) then
         print*, "velocity shear tensor is not right!"
         print*, u0, u1, u2
         print*, sigma00, sigma01, sigma02, 0.0d0
         print*, sigma01, sigma11, sigma12, 0.0d0
         print*, sigma02, sigma12, sigma22, 0.0d0
         print*, 0.0d0, 0.0d0, 0.0d0, sigma33
      endif

      end Subroutine

!-----------------------------------------------------------------------
      Subroutine checkPiandoutputViolation(Time, Dx, Dy, Vx, Vy, Ed, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22)
      ! Check size + tracelessness + transversality + positiveness of Tr(pi^2)
      ! and output all the violation points in the transverse plane

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ
      Double Precision Time, Dx, Dy

      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Vx(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Vy(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Integer I,J,K
      integer iFlag
      double precision violationType
      Double Precision trace_pi, trans, TrPi2

      Double Precision :: Tideal_scale, pi_scale
      Double Precision :: absNumericalzero = 1D-2
      Double Precision :: relNumericalzero = 1D-2  !Xsi_0 in Zhi's thesis

      Double Precision maxPiRatio
      Double Precision gamma_perp
      Common /maxPiRatio/ maxPiRatio

      iFlag = 0

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        violationType = 0D0
        Tideal_scale = sqrt(Ed(I,J,K)**2 + 3*PL(I,J,K)**2)
        
        TrPi2 = Pi00(I,J,K)**2+Pi11(I,J,K)**2+Pi22(I,J,K)**2
     &          +Pi33(I,J,K)**2
     &          -2*Pi01(I,J,K)**2-2*Pi02(I,J,K)**2+2*Pi12(I,J,K)**2

        pi_scale = sqrt(abs(TrPi2))
        !Positivity of Tr(pi^2)
        if(TrPi2 < -relNumericalzero*pi_scale) then
          violationType = violationType + 1.0D0
        endif

        !Largeness of pi tensor Tr(pi^2)
        If (pi_scale > max(maxPiRatio*Tideal_scale,
     &      absNumericalzero)) Then
          violationType = violationType + 0.1D0
        End If

        !trace of pi tensor
        trace_pi = Pi00(I,J,K)-Pi11(I,J,K)-Pi22(I,J,K)-Pi33(I,J,K)
        If(abs(trace_pi) > max(relNumericalzero*pi_scale,
     &     absNumericalzero) ) Then
          violationType = violationType + 0.01D0
        End If

        !transversality of pi tensor  u_mu pi^{mu,x} = 0
        gamma_perp = sqrt(1./(1. - Vx(I,J,K)**2 - Vy(I,J,K)**2 + 1D-30))
        trans = gamma_perp*(Pi01(I,J,K) - Vx(I,J,K)*Pi11(I,J,K) 
     &                      - Vy(I,J,K)*Pi12(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale, 
     &      absNumericalzero)) Then
          violationType = violationType + 0.001D0
        End If

        !transversality of pi tensor  u_mu pi^{mu,y} = 0
        trans = gamma_perp*(Pi02(I,J,K) - Vx(I,J,K)*Pi12(I,J,K) 
     &                      - Vy(I,J,K)*Pi22(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
          violationType = violationType + 0.001D0
        End If

        !transversality of pi tensor  u_mu pi^{mu,tau} = 0
        trans = gamma_perp*(Pi00(I,J,K) - Vx(I,J,K)*Pi01(I,J,K) 
     &                      - Vy(I,J,K)*Pi02(I,J,K))
        If (abs(trans) > max(relNumericalzero*pi_scale,
     &      absNumericalzero)) Then
          violationType = violationType + 0.001D0
        End If
        
        if(violationType > 0D0) then 
           iFlag = 1
           write(583, '(4F18.8)')Time, violationType, I*DX, J*DY
        endif

      End Do
      End Do
      End Do

      if(iFlag .eq. 0) then 
         write(583, '(4F18.8)')Time, 0D0, 10D0, 10D0
      endif
      End Subroutine
!-----------------------------------------------------------------------

      Subroutine checkBulkPi(failed, II, JJ, Time, Ed, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  PPI, say_level)
      ! Check the size of Bulk pressure

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ
      Integer failed, II, JJ ! (II,JJ): where it fails
      Double Precision Time

      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ)     !Bulk pressure

      Integer say_level  !Warning level

      Integer I,J,K
      Double Precision :: pressure_scale, bulkPi_scale
      Double Precision :: absNumericalzero = 1D-2
      Double Precision :: relNumericalzero = 1D-2  !Xsi_0 in Zhi's thesis

      Double Precision maxBulkPiRatio
      Common /maxBulkPiRatio/ maxBulkPiRatio

      failed = 0

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        !Largeness of bulk pressure
        pressure_scale = abs(PL(I,J,K))
        bulkPi_scale = abs(PPI(I,J,K))

        If (bulkPi_scale > max(maxBulkPiRatio*pressure_scale,
     &      absNumericalzero)) Then
          If (say_level>=9) Then
            Print*, "Time=", Time
            Print*, "Bulk Pi is larger than pressure!"
            Print*, "I,J=", I,J
            Print*, "Bulk Pi=", PPI(I,J,K)
            Print*, "Pressure=", PL(I,J,K)
          End If
          II = I
          JJ = J
          failed = 1
          return
        End If

      End Do
      End Do
      End Do

      End Subroutine
!-----------------------------------------------------------------------


      Subroutine checkBulkPiandoutputViolation(Time, Dx, Dy, Ed, PL,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  PPI)
      ! Check the size of Bulk pressure and output all violation 
      ! points in the transverse plane

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ
      Double Precision Time, Dx, Dy

      Double Precision Ed(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PL(NX0:NX,NY0:NY,NZ0:NZ)

      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ)     !Bulk pressure

      Integer violationType
      Integer I,J,K
      Double Precision :: pressure_scale, bulkPi_scale
      Double Precision :: absNumericalzero = 1D-2
      Double Precision :: relNumericalzero = 1D-2  !Xsi_0 in Zhi's thesis

      Double Precision maxBulkPiRatio
      Common /maxBulkPiRatio/ maxBulkPiRatio

      Do K=NZ0,NZ
      Do J=NYPhy0-3,NYPhy+3
      Do I=NXPhy0-3,NXPhy+3

        !Largeness of bulk pressure
        pressure_scale = abs(PL(I,J,K))
        bulkPi_scale = abs(PPI(I,J,K))

        If (bulkPi_scale > max(maxBulkPiRatio*pressure_scale,
     &      absNumericalzero)) Then
          write(584, '(3F18.8)')Time, I*DX, J*DY
        End If

      End Do
      End Do
      End Do

      End Subroutine
!-----------------------------------------------------------------------


      Subroutine printMore(id, I, J, Time,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  TT00,TT01,TT02, Ed, PL, Sd, Temp, Vx, Vy)
      ! Print values of Pi

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,id,I,J

      Double Precision Time

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision TT00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision TT01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision TT02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Temp(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Vx(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Vy(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision p00,p01,p02,p11,p12,p22,p33,TrPi

      Print*, "printPiAndTT: ID=", id
      Print*, "I,J=",I,J

      Print*, "Pi00=", Pi00(I,J,NZ0)
      Print*, "Pi01=", Pi01(I,J,NZ0)
      Print*, "Pi02=", Pi02(I,J,NZ0)
      Print*, "Pi11=", Pi11(I,J,NZ0)
      Print*, "Pi12=", Pi12(I,J,NZ0)
      Print*, "Pi22=", Pi22(I,J,NZ0)
      Print*, "Pi33*tau^2=", Pi33(I,J,NZ0)

      p00 = Pi00(I,J,NZ0)
      p01 = Pi01(I,J,NZ0)
      p02 = Pi02(I,J,NZ0)
      p11 = Pi11(I,J,NZ0)
      p12 = Pi12(I,J,NZ0)
      p22 = Pi22(I,J,NZ0)
      p33 = Pi33(I,J,NZ0)/(Time*Time) ! Pi33=pi^(3,3)*tau*tau

      TrPi = p00*p00+p11*p11+p22*p22+p33*Time*Time*p33*Time*Time
     &  -2*p01*p01-2*p02*p02+2*p12*p12

      Print*, "Tr(Pi^2)=",TrPi

      Print*, "TT00=", TT00(I,J,NZ0)
      Print*, "TT01=", TT01(I,J,NZ0)
      Print*, "TT02=", TT02(I,J,NZ0)
      Print*, "Ed=", Ed(I,J,NZ0)
      Print*, "PL=", PL(I,J,NZ0)
      Print*, "Sd=", Sd(I,J,NZ0)
      Print*, "Temp=", Temp(I,J,NZ0)
      Print*, "Vx=", Vx(I,J,NZ0)
      Print*, "Vy=", Vy(I,J,NZ0)
      Print*, "Time=", Time

      End Subroutine



!-----------------------------------------------------------------------
      Subroutine printDPc(id, I, J, Sd,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  DPc00,DPc01,DPc02,DPc33,DPc11,DPc12,DPc22)
      ! Print values of Pi

      Implicit None

      Integer NXPhy0,NXPhy,NYPhy0,NYPhy,NX0,NX,NY0,NY,NZ0,NZ,id,I,J

      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision DPc00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision DPc01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision DPc02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision DPc33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision DPc11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision DPc12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision DPc22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Print*, "DPc00,DPc*Sd=", DPc00(I,J,NZ), DPc00(I,J,NZ0)*Sd(I,J,NZ)
      Print*, "DPc01,DPc*Sd=", DPc01(I,J,NZ), DPc01(I,J,NZ0)*Sd(I,J,NZ)
      Print*, "DPc02,DPc*Sd=", DPc02(I,J,NZ), DPc02(I,J,NZ0)*Sd(I,J,NZ)
      Print*, "DPc11,DPc*Sd=", DPc11(I,J,NZ), DPc11(I,J,NZ0)*Sd(I,J,NZ)
      Print*, "DPc12,DPc*Sd=", DPc12(I,J,NZ), DPc12(I,J,NZ0)*Sd(I,J,NZ)
      Print*, "DPc22,DPc*Sd=", DPc22(I,J,NZ), DPc22(I,J,NZ0)*Sd(I,J,NZ)
      Print*, "DPc33,DPc*Sd=", DPc33(I,J,NZ), DPc33(I,J,NZ0)*Sd(I,J,NZ)

      End Subroutine
