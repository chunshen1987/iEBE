! Regualtion using tracelessness, transversality also

************************************************************************
!   Collections of Supportive Subroutines for VISH2+1 code
!
!   Author: Zhi Qiu (qiu@mps.ohio-state.edu)
!   -- For changes, see ChangeLogsInputFun.txt
!-----------------------------------------------------------------------
!   ToDos:
!   -- Add debug symbol (and if's)
!   -- Change many multiple-if's into single-line if'
!   -- Allow reading parameters outPiMovieDt and regMethod from CML & file
!   -- Add debug echo for parameters read from file
!
************************************************************************

#define say_level 3

      Subroutine prepareInputFun()
!     Purpose:
!       Do some common preparation for other functions in this file

      Implicit None

!----------Some global control symbols----------
      Integer debug ! used to control the amount of output
      Common /debug/ debug ! the larger the value, the more output on screen
!     Values for debug:
!     -- 0: no output
!     -- 1: minimum output (important values only)
!     -- 3: basic flow control output included (enter/leave subroutine)
!     -- 5: more flow control output included (if/do)
!     -- 9: variables value check output included

      Integer outPiTrace, outPiStrengthX
!     used to control whether the trace (radius only) of region (Pi>0.9*max_Pi) should be outputed to file
      Common /outPi/ outPiTrace, outPiStrengthX
!     outPiTrace: 0: no ouput; 1: output
!     -- Output PiTrace info. Format: (tau, minimal r where Pi is being regulated)
!     outPiStrengthX: 0: no ouput; 1: output
!     -- Output Pi info on X axes. Format (tau, list of Pi/Pi_max on X axes)

      Integer outPiComponents
!     if this is set to true, both the trace (see outPiTrace) and the strength on x axis (see outPiStrengthX) will be outputted for all components of Pi.
      Common /outPiComponents/ outPiComponents

      Double Precision outPiLastTime, outPiMovieDt, outPiMovieSize
      common /outPiMovie/ outPiLastTime, outPiMovieDt, outPiMovieSize
!     outPiCurTime stores the time when previous Pi output action happened
!     outPiMoveDt is the time interval by which a snapshot of Pi is outputted, 0.0 being no output
!     outPiMovieSize is the size of grid Pi is outputted: [0,outPiMovieSize]x[0,outPiMovieSize]

      Integer outEOnXY
      Common /outEOnXY/ outEOnXY

      Integer regMethod ! used to determine which method of regulation should be applied.
!     1: find R0 (see PiRatio) + Fermi-Dirac; 2: use tanh (see maxPiRatio)
      Common /regMethod/ regMethod

      Double Precision PiRatio ! used to determine R0; within r<R0, Pi/(e+p) < PiRatio
      Common /PiRatio/ PiRatio

      Double Precision maxPiRatio ! used in tanh regulation method: Pi=Pi_max*tanh(Pi/Pi_max), Pi_max=maxPiRatio*(e+p)
      Common /maxPiRatio/ maxPiRatio
      
      Double Precision maxBulkPiRatio ! used in tanh regulation method: Pi=Pi_max*tanh(Pi/Pi_max), Pi_max=maxPiRatio*(e+p)
      Common /maxBulkPiRatio/ maxBulkPiRatio

      Integer checkE
      Common /checkE/ checkE


!----------GLOBALS-----------------------------------------------
      debug = 0
      outPiTrace = 1
      outPiStrengthX = 0
      outPiComponents = 0
      outEOnXY = 0
      outPiLastTime = 0.0
      outPiMovieDt = 0.0
      outPiMovieSize = 130
      checkE = 0
!----------------------------------------------------------------


!     Read extraParas.inp file for these parameters
      OPEN(391,FILE='extraParas.inp',STATUS='OLD') ! a file contains some extra parameters

      READ(391,*) regMethod ! for the meaning of variables see their definitions
      Print *, "regMethod=", regMethod
      READ(391,*) PiRatio
      Print *, "PiRatio=", PiRatio
      READ(391,*) maxPiRatio
      Print *, "maxPiRatio=", maxPiRatio
      READ(391,*) maxBulkPiRatio
      Print *, "maxBulkPiRatio=", maxBulkPiRatio

      Close(391)

!     Establish output files
      If (outPiTrace .eq. 1) Then
        Open(392,FILE="movie/PiTrace.dat",STATUS='REPLACE')
        Close(392)
      EndIf
      If (outPiStrengthX .eq. 1) Then
        Open(394,FILE="movie/PiOnX.dat",STATUS='REPLACE')
        Close(394)
      EndIf
      If (outPiMovieDt >= 1e-6) Then
        Open(393,FILE="movie/PiMovie.dat",STATUS='REPLACE')
        Close(393)
      EndIf
      If (outPiComponents .eq. 1) Then
        Open(190,FILE="movie/PiComTrace.dat",STATUS='REPLACE')
        Close(190)
        Open(193,FILE="movie/Pi00X.dat",STATUS='REPLACE')
        Close(193)
        Open(194,FILE="movie/Pi01X.dat",STATUS='REPLACE')
        Close(194)
        Open(195,FILE="movie/Pi02X.dat",STATUS='REPLACE')
        Close(195)
        Open(196,FILE="movie/Pi11X.dat",STATUS='REPLACE')
        Close(196)
        Open(197,FILE="movie/Pi12X.dat",STATUS='REPLACE')
        Close(197)
        Open(198,FILE="movie/Pi22X.dat",STATUS='REPLACE')
        Close(198)
        Open(199,FILE="movie/Pi33X.dat",STATUS='REPLACE')
        Close(199)
      EndIf
      If (outEOnXY .eq. 1) Then
        Open(230,FILE="movie/eOnX.dat",STATUS='REPLACE')
        Close(230)
        Open(231,FILE="movie/eOnY.dat",STATUS='REPLACE')
        Close(231)
      EndIf
      If (checkE .eq. 1) Then
        Open(310,FILE="movie/checkE.dat",STATUS='REPLACE')
        Close(310)
      EndIf

      End Subroutine
!-----------------------------------------------------------------------



************************************************************************
      Subroutine readInputFromCML()
!     Purpose:
!     Read inputs from command line.

      Implicit None ! enforce explicit variable declaration

      Integer debug
      Common /debug/ debug

      Integer readFromCML ! variable for return value. 1: no CML argument; 0: succeeded in reading CML; stop otherwise
      Common /readFromCML/ readFromCML

      Integer IEOS, IEin, IInit
      Common /EOSSEL/ IEOS   !Type of EOS
      Common /IEin/ IEin     !  type of initialization  entropy/enrgy
      Common /Initialization/ IInit     ! type of initialization CGC/Glauber

      Double Precision A, Si0, EK, HWN, TRo0, TEta, TRA
      Common /AWNBC/ A,Si0 !A Nuclei Number  Si0, Cross Section for NN
      Common /EK/ EK, HWN  !EK(T0) constant related to energy density,HWN percent of Wounded Nucleon
      Common /thick/ TRo0, TEta, TRA  !Para in Nuclear Thickness Function

      Double Precision ViscousC, VisBeta, Visbulk, BulkTau, IRelaxBulk
      Common /ViscousC/ ViscousC, VIsBeta ! Related to Shear Viscosity
      Common /ViscousBulk/ Visbulk, BulkTau, IRelaxBulk  ! Related to bulk Visousity

      Double Precision ITeta, b, ddx, ddy, TT0
      Common /ITeta/ ITeta
      Common /bb/ b  !impact parameter
      Common/dxdy/ ddx, ddy
      Common /TT0/ TT0   ! T0, or tau_0

      Double Precision DT_1, DT_2 ! DT_1 is the standard time step, DT_2 is used as time step for early time (t<0.6 fm/c)
      Common /Timestep/ DT_1, DT_2

      Double Precision Edec
      Common/Edec/Edec    !decoupling temperature

      Integer IEOS2dec
      Common/IEOS2dec/ IEOS2dec  ! IEOS=2 decouple by gluon/pion

      Double Precision Rx2, Ry2 ! <x^2> and <y^2> used in Gaussian initial condition
      Common /RxyBlock/ Rx2, Ry2

      Integer NDX, NDY, NDT
      Common /NXYTD/ NDX, NDY, NDT

      Double Precision T0
      Common /T0/ T0

      Double Precision R0Bdry
      Common /R0Bdry/ R0Bdry
      Integer LS
      Common /LS/ LS

      Integer QNum, ArgIndex ! QNum is the total number of arguments, ArgIndex gives the index to the one currently reading

      Integer Qkind

!     Warning: when adding more parameters, be aware that fortran does not allow a line with more than (roughly) 73 Characters
!     Warning: float argument needs to have a period, so 1 should be written as 1.0 instead
      Character(len=40) :: ViscousCStr,bStr,T0Str,Rx2Str,Ry2Str
      Character(len=40) :: EKStr, EDecStr, dTStr
      Character(len=40) :: IInitStr,IEOSStr
      Character(len=40) :: MovieDtStr, QkindStr
      Character(len=40) :: LSStr, R0BdryStr, VisBetaStr

      If (debug>=3) Print *, "* readInputFromCML started"

      QNum = iargc ()

      readFromCML = 0 ! default to read from command line

      If (QNum .eq. 0) Then
        If (debug>=5) Print *, "-- No arguments."
        readFromCML = 1
      ElseIf (QNum .eq. 11) Then

        ArgIndex = 0

!       Read values from command line

        If (debug>=5) Print *,"-- Read command line parameters: 11"

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, IEOSStr)
        Read(unit=IEOSStr, fmt='(I5)') IEOS

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, QkindStr)
        Read(unit=QkindStr, fmt='(I5)') Qkind
        If (Qkind .eq. 1) Then
          A = 63.546
        ElseIf (Qkind .eq. 2) Then
          A = 196.96655
        EndIf

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, IInitStr)
        Read(unit=IInitStr, fmt='(I5)') IInit

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, dTStr)
        Read(unit=dTStr, fmt='(f15.8)') dT_1

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, ViscousCStr)
        Read(unit=ViscousCStr, fmt='(f15.8)') ViscousC

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, bStr)
        Read(unit=bStr, fmt='(f15.8)') b

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, Rx2Str)
        Read(unit=Rx2Str, fmt='(f15.8)') Rx2

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, Ry2Str)
        Read(unit=Ry2Str, fmt='(f15.8)') Ry2

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, EKStr)
        Read(unit=EKStr, fmt='(f15.8)') EK

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, T0Str)
        Read(unit=T0Str, fmt='(f15.8)') T0

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, EDecStr)
        Read(unit=EDecStr, fmt='(f15.8)') EDec

        ! give default values to LS and R0Bdry
        LS = 130
        R0Bdry = 12.0
        VisBeta = 0.5   !\tau_Pi=VisBeta*6.0\eta /(ST)


      ElseIf (QNum .eq. 13) Then ! read futher


        ArgIndex = 0

!       Read values from command line

        If (debug>=5) Print *,"-- Read command line parameters: 13"

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, IEOSStr)
        Read(unit=IEOSStr, fmt='(I5)') IEOS

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, QkindStr)
        Read(unit=QkindStr, fmt='(I5)') Qkind
        If (Qkind .eq. 1) Then
          A = 63.546
        ElseIf (Qkind .eq. 2) Then
          A = 196.96655
        EndIf

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, IInitStr)
        Read(unit=IInitStr, fmt='(I5)') IInit

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, dTStr)
        Read(unit=dTStr, fmt='(f15.8)') dT_1

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, ViscousCStr)
        Read(unit=ViscousCStr, fmt='(f15.8)') ViscousC

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, bStr)
        Read(unit=bStr, fmt='(f15.8)') b

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, Rx2Str)
        Read(unit=Rx2Str, fmt='(f15.8)') Rx2

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, Ry2Str)
        Read(unit=Ry2Str, fmt='(f15.8)') Ry2

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, EKStr)
        Read(unit=EKStr, fmt='(f15.8)') EK

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, T0Str)
        Read(unit=T0Str, fmt='(f15.8)') T0

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, EDecStr)
        Read(unit=EDecStr, fmt='(f15.8)') EDec

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, LSStr)
        Read(unit=LSStr, fmt='(I5)') LS

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, R0BdryStr)
        Read(unit=R0BdryStr, fmt='(f15.8)') R0Bdry

        VisBeta = 0.5   !\tau_Pi=VisBeta*6.0\eta /(ST)

      ElseIf (QNum .eq. 14) Then ! read futher, includes tau_pi (VisBeta)


        ArgIndex = 0

!       Read values from command line

        If (debug>=5) Print *,"-- Read command line parameters: 14"

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, IEOSStr)
        Read(unit=IEOSStr, fmt='(I5)') IEOS

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, QkindStr)
        Read(unit=QkindStr, fmt='(I5)') Qkind
        If (Qkind .eq. 1) Then
          A = 63.546
        ElseIf (Qkind .eq. 2) Then
          A = 196.96655
        EndIf

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, IInitStr)
        Read(unit=IInitStr, fmt='(I5)') IInit

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, dTStr)
        Read(unit=dTStr, fmt='(f15.8)') dT_1

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, ViscousCStr)
        Read(unit=ViscousCStr, fmt='(f15.8)') ViscousC

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, bStr)
        Read(unit=bStr, fmt='(f15.8)') b

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, Rx2Str)
        Read(unit=Rx2Str, fmt='(f15.8)') Rx2

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, Ry2Str)
        Read(unit=Ry2Str, fmt='(f15.8)') Ry2

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, EKStr)
        Read(unit=EKStr, fmt='(f15.8)') EK

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, T0Str)
        Read(unit=T0Str, fmt='(f15.8)') T0

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, EDecStr)
        Read(unit=EDecStr, fmt='(f15.8)') EDec

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, LSStr)
        Read(unit=LSStr, fmt='(I5)') LS

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, R0BdryStr)
        Read(unit=R0BdryStr, fmt='(f15.8)') R0Bdry

        ArgIndex = ArgIndex + 1
        Call getarg(ArgIndex, VisBetaStr)
        Read(unit=VisBetaStr, fmt='(f15.8)') VisBeta !\tau_Pi=VisBeta*6.0\eta /(ST)

      Else

        Print*,"Need: IEOS(Integer)",
     &    "(EOSQ:0  EOSI: 2, SM-EOSQ: 5, EOSL: 4(Katz05 data)), ",
     &    "Qkind(I)(1:Cu, 2:Au), ",
     &    "IInit(I)(0:Gaussian, 1:Glauber, 2:CGC), ",
     &    "dT(Double), ",
     &    "ViscousC(Double), b(D), Rx2(D), Ry2(D), EK(D), T0(D), ",
     &    "EDec(D), [LS(I), R0Bdry(D)], [VisBeta(D)]"
        Stop

      EndIf




!     Use default values for the rest of the parameters
        If (debug>=5) Print *, "-- Assign values to other parameters"
        IEOS2dec = 1    !If IEOS=2 (EOSI),  0: decouple by gluon 1: decouple by pions
        IEin = 0        !type of initialization  0:initialize by energy density 1: initialize by entropy density
        HWN = 1.0       !HWN  % of WN   (1-HWN) % of BC
        Si0 = 4.0       !cross section for NN
        TRo0 = 0.17     !nuclear density  (unit: fm^-3)
        TEta = 0.54     !dIffussion parameter of nuclear radii  (unit: fm)
        TRA = 6.37      !nuclear radii (unit: fm)

        VisBulk = 0.0   !VisBulk=C;  Xi/s= C* (Xi/S)_min  (C=0, no bulk vis; or C>1 )
        IRelaxBulk = 2  !type of bulk relaxation time (0: critical slowing down; 1: constant Relax Time
                            !2: \tau_PI=1.5/(2\piT))
        BulkTau = 5.0 !constant relaxation time (fm/c) (require input IRelaxBulk=1)

        NDX = 2
        NDY = 2     ! freeze-out step in x and y direction
        NDT = 5       ! freeze-out step in \tau direction

        Write (*,*) "Have:", "EOS=", IEOS, "Qkind=", Qkind, ! write out parameter for a check
     &    "Initialization=", IInit, "dT=", dT_1,
     &    "eta/s=",ViscousC,"b=",b,"Rx2=",Rx2,"Ry2=",Ry2,
     &    "EK=", EK, "tau0=", T0, "EDec=", EDec,
     &    "LS=", LS, "R0Bdry", R0Bdry, "VisBeta=", VisBeta


      If (debug>=3) Print *, "* readInputFromCML finished"

      End Subroutine
!-----------------------------------------------------------------------




************************************************************************
      Subroutine readInputFromCML2()
!     Purpose:
!     Read inputs from command line.

      Implicit None ! enforce explicit variable declaration

      Integer debug
      Common /debug/ debug

      Integer IEOS, IEin, IInit
      Common /EOSSEL/ IEOS   !Type of EOS
      Common /IEin/ IEin     !  type of initialization  entropy/enrgy
      Common /Initialization/ IInit     ! type of initialization CGC/Glauber

      Double Precision A, Si0, EK, HWN, TRo0, TEta, TRA
      Common /AWNBC/ A,Si0 !A Nuclei Number  Si0, Cross Section for NN
      Common /EK/ EK, HWN  !EK(T0) constant related to energy density,HWN percent of Wounded Nucleon
      Common /thick/ TRo0, TEta, TRA  !Para in Nuclear Thickness Function

      Double Precision ViscousC, VisBeta, Visbulk, BulkTau, IRelaxBulk
      Common /ViscousC/ ViscousC, VisBeta
      Common /ViscousBulk/ Visbulk, BulkTau, IRelaxBulk  ! Related to bulk Visousity

      Double Precision ITeta, b, ddx, ddy, TT0
      Common /ITeta/ ITeta
      Common /bb/ b  !impact parameter
      Common/dxdy/ ddx, ddy
      Common /TT0/ TT0   ! T0, or tau_0

      Double Precision DT_1, DT_2 ! DT_1 is the standard time step, DT_2 is used as time step for early time (t<0.6 fm/c)
      Common /Timestep/ DT_1, DT_2
      Double Precision DX, DY
      Common /DXY/ DX, DY

      Double Precision Edec
      Common/Edec/Edec    !decoupling temperature

      Integer IEOS2dec
      Common/IEOS2dec/ IEOS2dec  ! IEOS=2 decouple by gluon/pion

      Double Precision Rx2, Ry2 ! <x^2> and <y^2> used in Gaussian initial condition
      Common /RxyBlock/ Rx2, Ry2

      Double Precision sFactor ! multiplicity factor on entropy density
      Common /sFactor/ sFactor

      Integer NDX, NDY, NDT
      Common /NXYTD/ NDX, NDY, NDT

      Double Precision T0
      Common /T0/ T0

      Double Precision R0Bdry
      Common /R0Bdry/ R0Bdry
      Double Precision R0, Aeps
      Common /R0Aeps/ R0,Aeps
      Integer LS
      Common /LS/ LS

      Integer QNum, ArgIndex ! QNum is the total number of arguments, ArgIndex gives the index to the one currently reading

      Character*60 :: buffer
      Character*20 :: varName
      Integer IResult
      Double Precision DResult

      Aeps = 0.05D0

      If (debug>=3) Print *, "* readInputFromCML2 started"

      QNum = iargc ()

      sFactor = 1D0

      Do ArgIndex = 1, QNum
        Call getarg(ArgIndex, buffer)
        Call processAssignment(buffer, "=", varName, IResult, DResult)

        If (varName=="ieos") IEOS=IResult ! EOSQ: 0;  EOSI: 2; SM-EOSQ: 5; EOSL: 4(Katz05 data); s95p-PCE: 6

        If (varName=="iinit") IInit=IResult ! 0: Gaussian; 1: Optical Glauber; 2: From Initial/InitialEd.dat
        If (varName=="init") IInit=IResult
        If (varName=="ii") IInit=IResult

        If (varName=="iein") IEin=IResult ! 0: initialize by energy density; 1: initialize by entropy density
        If (varName=="iin") IEin=IResult

        If (varName=="dt") dT_1=DResult ! dT, DX, DY
        If (varName=="dx") dX=DResult
        If (varName=="dy") dy=DResult

        If (varName=="x2") Rx2=DResult ! <x^2> and <y^2>
        If (varName=="rx2") Rx2=DResult
        If (varName=="x^2") Rx2=DResult
        If (varName=="x_2") Rx2=DResult
        If (varName=="y2") Ry2=DResult
        If (varName=="ry2") Ry2=DResult
        If (varName=="y^2") Ry2=DResult
        If (varName=="y_2") Ry2=DResult

        If (varName=="ek") EK=DResult ! centeral energy density, in GeV/fm^3 (for Glauber and Gaussian)
        If (varName=="e0") EK=DResult
        If (varName=="ec") EK=DResult
        If (varName=="ecen") EK=DResult
        If (varName=="e_cen") EK=DResult

        If (varName=="sk") EK=DResult ! centeral entropy density, in fm^-3 (for Glauber)
        If (varName=="s0") EK=DResult
        If (varName=="sc") EK=DResult
        If (varName=="scen") EK=DResult
        If (varName=="s_cen") EK=DResult

        If (varName=="factor") sFactor=DResult ! VER-1.29RC3: final multiplicity factor on entropy density
        If (varName=="fac") sFactor=DResult
        If (varName=="ff") sFactor=DResult

        If (varName=="b") b=DResult ! impact parameter

        If (varName=="edec") EDec=DResult ! decouple energy density, in GeV/fm^3
        If (varName=="e_dec") EDec=DResult
        If (varName=="e_d") EDec=DResult

        If (varName=="t0") T0=DResult ! initial proper time tau_0, in fm/c
        If (varName=="viscousc") ViscousC=DResult ! variations for shear viscosities/entropy density ratio
        If (varName=="es") ViscousC=DResult
        If (varName=="e_s") ViscousC=DResult
        If (varName=="etas") ViscousC=DResult
        If (varName=="eta_s") ViscousC=DResult
        If (varName=="vis") ViscousC=DResult
        If (varName=="viscousc") ViscousC=DResult

        If (varName=="ils") LS=IResult ! Lattice size and R0Boudary
        If (varName=="r0") R0Bdry=DResult
        If (varName=="r0bdry") R0Bdry=DResult

        If (varName=="a") A=DResult ! Atom number A

        If (varName=="ndx") NDX=IResult ! freeze-out cell sizes
        If (varName=="ndy") NDY=IResult
        If (varName=="ndt") NDT=IResult

        If (varName=="visbeta") VisBeta=DResult ! VisBeta, used for proper time tau_pi

      End Do ! ArgIndex

      If (debug>=3) Print *, "* readInputFromCML finished"

      End Subroutine
!-----------------------------------------------------------------------




************************************************************************
      Subroutine processAssignment(string, separator,
     &                            varName, IResult, DResult)
!     This subroutine process a string assignment.
!     First it seprate string into LHS and RHS according to separator.
!     Then the LHS is converted into variable using only lower case
!     letters, and the RHS is converted into numerical values.
!     The variable IResult holds result for integer and DResult holds
!     one for double.
!     Convention: integer-valued variable should start with "I" or "N"
!     for its name.

      Implicit None

      Character (*) :: string, varName
      Character*60 :: LHS, RHS
      Character separator
      Integer IResult
      Double Precision DResult

      Integer break_here, I, cha

      varName = ""

      break_here = index(string, separator)
      LHS = adjustl(string(:break_here-1))
      RHS = adjustl(string(break_here+1:))

      ! convert LHS to lower case:
      Do I = 1,len_trim(LHS)
        cha = ichar(LHS(I:I))
        If (cha>=65 .and. cha<90) Then
          varName(I:I) = char(cha+32)
        Else
          varName(I:I) = LHS(I:I)
        EndIf
      EndDo

      ! convert RHS to numerics:
      If (varName(1:1)=="i" .or. varName(1:1)=="n") Then
        Read(RHS, fmt='(I5)') IResult
      Else
        Read(RHS, fmt='(f15.8)') DResult
      EndIf

      End Subroutine
!-----------------------------------------------------------------------




************************************************************************
      Subroutine getInitialR0(PU0,PU1,PU2,PU3,U0,U1,U2,U3,DX,DY,DZ,DT,
     &  DPc00,DPc01,DPc02,DPc33,DPc11,DPc22,DPc12,DDU0,DDU1,DDU2,
     &  Temp,Temp0,SiLoc,DLnT,Time, NXPhy0,NYPhy0,NXPhy,NYPhy,
     &  NX0,NX,NY0,NY,NZ0,NZ,Ed,Sd,PL,VCoefi)
!     Purpose:
!       Return a suitable initial (before the initialization of
!       Pi(mu,nu)) R0 (thru common) to regulate Pi(mu,nu)

      Implicit None

      Integer debug
      Common /debug/ debug

      Integer I,J,K
      Integer NX0,NY0,NZ0,NX,NY,NZ
      Integer NXPhy0,NYPhy0,NXPhy,NYPhy

      Common /dxdy/ ddx, ddy
      Double Precision ddx, ddy

      Double Precision PU0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity from last time step
      Double Precision PU1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision PU2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision PU3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Double Precision U0(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision U1(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision U2(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity
      Double Precision U3(NX0:NX, NY0:NY, NZ0:NZ) !Four velocity

      Double Precision DPc00(NX0:NX, NY0:NY, NZ0:NZ) ! DIfferential part of Pi source term
      Double Precision DPc01(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DPc02(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DPc33(NX0:NX, NY0:NY, NZ0:NZ) !

      Double Precision DPc11(NX0:NX, NY0:NY, NZ0:NZ) ! DIfferential part of Pi source term
      Double Precision DPc12(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DPc22(NX0:NX, NY0:NY, NZ0:NZ) !

      Double Precision DDU0(NX0:NX, NY0:NY, NZ0:NZ) ! DIfferential part of Pi source term
      Double Precision DDU1(NX0:NX, NY0:NY, NZ0:NZ) !
      Double Precision DDU2(NX0:NX, NY0:NY, NZ0:NZ) !

      Double Precision Temp0(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature  in last time step
      Double Precision Temp(NX0:NX, NY0:NY, NZ0:NZ) !Local Temperature
      Double Precision SiLoc(NX0:NX, NY0:NY, NZ0:NZ) ! Local expansion rate \sita
      Double Precision DLnT(NX0:NX, NY0:NY, NZ0:NZ) ! DlnT(x,y) terms

      Double Precision R0, Aeps, Accu
      Common /R0Aeps/ R0,Aeps
      Common /Accu/Accu  ! A parameter to determine the accuracy of Calculation
                         !Accu=3.0 used 3pt formula to cal derivative. Accu=5.0 use 5pt formula to cal Deriv.

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !pressure density
      Double Precision VCoefi(NX0:NX, NY0:NY, NZ0:NZ) !viscous coeficient shear viscosity eta
      Double Precision RMin, PiEPRatio, SigmaLargeness, EAndP

      Double Precision ViscousC, VisBeta
      Common /ViscousC / ViscousC, VisBeta

      Double Precision PiRatio ! used to determine R0; within r<R0, Pi/(e+p) < PiRatio
      Common /PiRatio/ PiRatio ! should already be setuped in prepareInputFun function

      Double Precision D0U0,D0U1,D0U2,D1U0,D1U1,D1U2,D2U0,D2U1,D2U2
      Double Precision CS,DT,DX,DY,DZ,Time,DU0,DU1,DU2

      Integer regMethod
      Common /regMethod/ regMethod

      If (debug>=3) Print *, "* Start GetInitialR0"

!      Aeps = 0.5

      If (regMethod .eq. 1) Then
        If (debug>=5) Print *,"-- To branch regMethod=1"

        DO 791 K=NZ0,NZ
        DO 791 J=NYPhy0,NYPhy
        DO 791 I=NXPhy0,NXPhy

        D0U0=(U0(I,J,K)-PU0(I,J,K))/DT
        D0U1=(U1(I,J,K)-PU1(I,J,K))/DT
        D0U2=(U2(I,J,K)-PU2(I,J,K))/DT

        If(abs(Accu-3.0).le.0.00001) Then  !3pt formula
          D1U0=(U0(I+1,J,K)-U0(I-1,J,K))/(2.0*DX)
          D1U1=(U1(I+1,J,K)-U1(I-1,J,K))/(2.0*DX)
          D1U2=(U2(I+1,J,K)-U2(I-1,J,K))/(2.0*DX)
          D2U0=(U0(I,J+1,K)-U0(I,J-1,K))/(2.0*DY)
          D2U1=(U1(I,J+1,K)-U1(I,J-1,K))/(2.0*DY)
          D2U2=(U2(I,J+1,K)-U2(I,J-1,K))/(2.0*DY)
        ElseIf (abs(Accu-5.0).le.0.00001) Then !5pt formula
          D1U0=(U0(I+1,J,K)*2.0d0/3.0d0-U0(I-1,J,K)*2.0d0/3.0d0
     &        -U0(I+2,J,K)/12.0d0+U0(I-2,J,K)/12.0d0)/DX
          D1U1=(U1(I+1,J,K)*2.0d0/3.0d0-U1(I-1,J,K)*2.0d0/3.0d0
     &        -U1(I+2,J,K)/12.0d0+U1(I-2,J,K)/12.0d0)/DX
          D1U2=(U2(I+1,J,K)*2.0d0/3.0d0-U2(I-1,J,K)*2.0d0/3.0d0
     &        -U2(I+2,J,K)/12.0d0+U2(I-2,J,K)/12.0d0)/DX
          D2U0=(U0(I,J+1,K)*2.0d0/3.0d0-U0(I,J-1,K)*2.0d0/3.0d0
     &        -U0(I,J+2,K)/12.0d0+U0(I,J-2,K)/12.0d0)/DY
          D2U1=(U1(I,J+1,K)*2.0d0/3.0d0-U1(I,J-1,K)*2.0d0/3.0d0
     &        -U1(I,J+2,K)/12.0d0+U1(I,J-2,K)/12.0d0)/DY
          D2U2=(U2(I,J+1,K)*2.0d0/3.0d0-U2(I,J-1,K)*2.0d0/3.0d0
     &        -U2(I,J+2,K)/12.0d0+U2(I,J-2,K)/12.0d0)/DY
        Else
          Print*, "Wrong input for Accu:",
     &    "Accu=3or5 for 3pt or 5pt cal of deriv."
        EndIf

        CS=(D0U0+D1U1+D2U2+U0(I,J,K)/Time)/3.0

        DU0=U0(I,J,K)*D0U0+U1(I,J,K)*D1U0+U2(I,J,K)*D2U0
        DU1=U0(I,J,K)*D0U1+U1(I,J,K)*D1U1+U2(I,J,K)*D2U1
        DU2=U0(I,J,K)*D0U2+U1(I,J,K)*D1U2+U2(I,J,K)*D2U2

        DPc00(I,J,K)=D0U0-U0(I,J,K)*DU0+CS*(U0(I,J,K)**2-1.0)
        DPc01(I,J,K)=0.5*(D0U1-D1U0)-0.5*(U1(I,J,K)*DU0+U0(I,J,K)*DU1)
     &              +CS*(U1(I,J,K)*U0(I,J,K))
        DPc02(I,J,K)=0.5*(D0U2-D2U0)-0.5*(U2(I,J,K)*DU0+U0(I,J,K)*DU2)
     &              +CS*(U2(I,J,K)*U0(I,J,K))
        DPc33(I,J,K)=CS-U0(I,J,K)/Time
        DPc11(I,J,K)=(-1.0)*D1U1-U1(I,J,K)*DU1+CS*(U1(I,J,K)**2+1.0)
        DPc22(I,J,K)=(-1.0)*D2U2-U2(I,J,K)*DU2+CS*(U2(I,J,K)**2+1.0)
        DPc12(I,J,K)=(-0.5)*(D2U1+D1U2)-0.5*(U1(I,J,K)*DU2
     &              +U2(I,J,K)*DU1)+CS*(U1(I,J,K)*U2(I,J,K))
 791    Continue

        RMin = NX*ddx+NY*ddy !---Upper-Bound-R0---

        DO 3007 K=NZ0,NZ !Check for Pi tensor
        DO 3007 J=NYPhy0,NYPhy
        DO 3007 I=NXPhy0,NXPhy
          EAndP = Abs(Ed(I,J,K)+PL(I,J,K))
          SigmaLargeness = 1/7.0*(Abs(DPc00(I,J,K))+
     &      Abs(DPc01(I,J,K))+Abs(DPc02(I,J,K))+Abs(DPc33(I,J,K))+
     &      Abs(DPc11(I,J,K))+Abs(DPc12(I,J,K))+Abs(DPc22(I,J,K)))
          PiEPRatio=2*ViscousC*Sd(I,J,K)*SigmaLargeness/EAndP

        If (PiEPRatio > PiRatio) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < RMin) Then
            RMin = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf

 3007   Continue
        R0 = RMin
      ElseIf (regMethod .eq. 2) Then ! use maximun possible R0
        If (debug>=5) Print *, "-- To branch regMethod=2"
        R0 = NX*ddx+NY*ddy
      Else  ! use R0=12
        If (debug>=5) Print *, "-- To other regMethod branch..."
        R0 = 12.0
      EndIf ! corresponding to the one on variable "regMethod"

      If (debug>=1) Print *,"R0=",R0 ! print out R0

      If (debug>=3) Print *, "& GetInitialR0 finished."

      End Subroutine
!-----------------------------------------------------------------------



************************************************************************
      Subroutine determineR0(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,Sd,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33)
!     Purpose:
!       Return a suitable R0 (thru common) to regulate Pi(mu,nu)

      Implicit None

      Integer debug
      Common /debug/ debug

      Integer NX0,NY0,NZ0,NX,NY,NZ
      Integer I,J,K

      Common /dxdy/ ddx, ddy
      Double Precision ddx, ddy

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor 
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision R0, Aeps
      Common /R0Aeps/ R0, Aeps

      Double Precision PiRatio ! used to determine R0; within r<R0, Pi/(e+p) < PiRatio
      Common /PiRatio/ PiRatio ! should already be setuped in prepareInputFun function


      Double Precision PiLargeness ! as a measurement of how large Pi is
      Double Precision PiEPRatio, EAndP
      Double Precision RMin ! an intermedia variable that trace the radius of the largest possible region

      Integer regMethod
      Common /regMethod/ regMethod

      If (debug>=3) Print *, "* Start DetermineR0"

!      Aeps = 0.4

      If (regMethod .eq. 1) Then
        If (debug>=5) Print *, "-- To branch regMethod=1"
        RMin = NX*ddx+NY*ddy !---Upper-Bound-R0---
        DO 3007 K=NZ0,NZ !Check for Pi tensor
        DO 3007 J=NY0,NY
        DO 3007 I=NX0,NX
          EAndP = Abs(Ed(I,J,K)+PL(I,J,K))
          PiLargeness = 1/7.0*(Abs(Pi00(I,J,K))+
     &      Abs(Pi01(I,J,K))+Abs(Pi02(I,J,K))+Abs(Pi33(I,J,K))+
     &      Abs(Pi11(I,J,K))+Abs(Pi12(I,J,K))+Abs(Pi22(I,J,K)))
          PiEPRatio=PiLargeness/EAndP

          If (PiEPRatio > PiRatio) Then
              If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < RMin) Then
                RMin = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
              EndIf
          EndIf
 3007   Continue
        R0 = RMin
      ElseIf (regMethod .eq. 2) Then ! use maximun possible R0
        If (debug>=5) Print *, "-- To regMethod=2 branch"
        R0 = (NX*ddx+NY*ddy)*2.0
      Else ! use R0=12.0
        If (debug>=5) Print *, "-- To other regMethod branch"
        R0 = 12.0
      EndIf

      If (debug>=1) Print *,"R0=",R0 ! output R0

      If (debug>=3) Print *, "* DetermineR0 finished"

      End Subroutine
!-----------------------------------------------------------------------


************************************************************************
      Subroutine regulateBulkPi(regStr,Time,NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
     &  Ed,PL,PPI,II,JJ)
!     Purpose:
!       Regulate Bulk pressure tensor by restrain it under a maximum
!       value using tanh function

      Implicit None

      Integer debug
      Common /debug/ debug

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Integer I,J,K,II,JJ,regStr

      Common/dxdy/ ddx, ddy ! lattice spacing
      Double Precision ddx, ddy

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ) ! Bulk Pressure Tensor

      Double Precision BulkPi   

      Integer regMethod
      Common /regMethod/ regMethod

      Double Precision :: Xsi0 = 1D0  !adaptive zero
      Double Precision :: pressure_scale, bulkPi_scale
      Double Precision regStrength

      Double Precision maxBulkPiRatio
      Common /maxBulkPiRatio/ maxBulkPiRatio

      Xsi0 = 1D-2/(regStr+1D0) ! VER-1.29RC: adaptive zero chooser VER-1.29RC4: bug fix: regStr -> regStr+1D0

      If (debug >= 3) Print *, "* Start RegulateBulkPi"

      If (regMethod == 2) Then ! do tanh regulation

        If (debug>=5) Print *, "-- To branch regMethod=2"
        If (debug>=5) Print *, "-- Pi regulation at time",Time

        DO 3019 K=NZ0,NZ
        DO 3019 J=NY0,NY
        DO 3018 I=NX0,NX

        regStrength = 1D-30
        
        pressure_scale = abs(PL(I,J,K))

        BulkPi = PPI(I,J,K)

        ! get Bulk pi scale
        bulkPi_scale = abs(BulkPi) + 1d-30
        if(bulkPi_scale .ne. bulkPi_scale) then
           print*, "Bulk Pi is NaN, I,J =", I, J
           stop
        endif

        ! find regulation strength using largeness comparison
        regStrength = max(bulkPi_scale/(maxBulkPiRatio*pressure_scale), 
     &                    regStrength)

        If ( say_level >=9 ) Then
          If (I==II.AND.J==JJ) Then
            Print*, "I,J=",I,J
            Print*, "regStrength=", regStrength
            Print*, "BulkPi = ", bulkPi_scale
            Print*, "maxPi=", maxBulkPiRatio*pressure_scale
            Print*, "PL=", PL(I,J,K)
            Print*, "Xsi0=", Xsi0
          endif
          If (I==0.AND.J==0) Then
            Print*, "I,J=",I,J
            Print*, "regStrength=", regStrength
            Print*, "BulkPi = ", bulkPi_scale
            Print*, "maxPi=", maxBulkPiRatio*pressure_scale
            Print*, "PL=", PL(I,J,K)
            Print*, "Xsi0=", Xsi0
          
          End If
        End If !If (debug>=9)

        PPI(I,J,K)=PPI(I,J,K)*(tanh(regStrength)/regStrength) ! Bulk pressure PPI is regulated here

3018    Continue
3019    Continue

      EndIf ! on regMethod

      If (debug>=3) Print *, "* RegulateBulkPi finished"

      End Subroutine
!-----------------------------------------------------------------------------

************************************************************************
      Subroutine regulatePi(regStr,Time,NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
     &  Ed,PL,PPI,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33,Vx,Vy,II,JJ)
!     Purpose:
!       Regulate Pi(mu,nu) tensor by restrain it under a maximum
!       value using tanh function

      Implicit None

      Integer debug
      Common /debug/ debug

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Integer I,J,K,II,JJ,regStr

      Common/dxdy/ ddx, ddy ! lattice spacing
      Double Precision ddx, ddy

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision Vx(NX0:NX, NY0:NY, NZ0:NZ)
      Double Precision Vy(NX0:NX, NY0:NY, NZ0:NZ)

      Double Precision PPI(NX0:NX, NY0:NY, NZ0:NZ) ! Bulk Pressure Tensor

      Double Precision Te00,Te01,Te02,Te11,Te12,Te22,Te33 ! These are T(mu,nu) in equilibrium

      Double Precision p00,p01,p02,p11,p12,p22,p33,vvx,vvy! These are just Pi at I,J,K (in a loop)

      Double Precision TrPi2 ! Tr(pi^2)

      Double Precision rTrPi2EAndP ! ratio between sqrt(Tr(pi^2)) and e+p

      Integer regMethod
      Common /regMethod/ regMethod

      Double Precision :: Xsi0 = 1D0  !adaptive zero
      Double Precision :: Tideal_scale, pi_scale
      Double Precision regStrength

      Double Precision maxPiRatio
      Common /maxPiRatio/ maxPiRatio
      Double Precision maxPi ! maxPi = maxPiRatio*(e+p)

      Double Precision PiPiMaxRatio1, PiPiMaxRatio2
      Double Precision PiPiMaxRatio3, PiPiMaxRatio4
      Double Precision rTrPi1, rTrPi2, rTrPi3, rTrPi4 ! radii determined by comparing PiPiMaxRatio to piRatioMax, piRatioAvg
      Double Precision gridInFz, inTP1, inTP2, inTP3, inTP4 ! count number of lattice points inside freezeout surface, TrPi1, TrPi2, TrPi3, TrPi4

      Double Precision PiAvg, PiRegAvg
      Integer PiCheckFlag, PiRegCheckFlag
      Double Precision gamma_perp

      Double Precision PiTr, PiTrSum, trans

      Xsi0 = 1D-2/(regStr+1D0) ! VER-1.29RC: adaptive zero chooser VER-1.29RC4: bug fix: regStr -> regStr+1D0

      If (debug >= 3) Print *, "* Start RegulatePi"

      If (regMethod == 2) Then ! do tanh regulation

        If (debug>=5) Print *, "-- To branch regMethod=2"
        If (debug>=5) Print *, "-- Pi regulation at time",Time

        DO 3009 K=NZ0,NZ
        DO 3009 J=NY0,NY
        DO 3008 I=NX0,NX

        regStrength = 1D-30

        vvx = Vx(I,J,K)
        vvy = Vy(I,J,K)
        gamma_perp = 1./sqrt(1. - vvx**2 - vvy**2 + 1D-30)
        Tideal_scale = sqrt(Ed(I,J,K)**2 + 3*PL(I,J,K)**2)

        p00 = Pi00(I,J,K)
        p01 = Pi01(I,J,K)
        p02 = Pi02(I,J,K)
        p11 = Pi11(I,J,K)
        p12 = Pi12(I,J,K)
        p22 = Pi22(I,J,K)
        p33 = Pi33(I,J,K)   ! pi are in t-xyz coordinate

        ! calculate Tr(pi^2)
        TrPi2 = p00*p00+p11*p11+p22*p22+p33*p33
     &    -2*p01*p01-2*p02*p02+2*p12*p12
        pi_scale = sqrt(abs(TrPi2)) + 1D-30

        ! find regulation strength

        ! first, tracelessness
        PiTr = p00-p11-p22-p33
        regStrength = max(abs(PiTr)/(Xsi0*MaxPiRatio*pi_scale), 
     &                    regStrength)

        !If (I==0.AND.J==0) Then
        !  Print*, "I,J=",I,J
        !  Print*, "regStrength1=", regStrength
        !End If

        ! next transversality
        trans = gamma_perp*(p01-vvx*p11-vvy*p12)
        regStrength = max(abs(trans)/(Xsi0*MaxPiRatio*pi_scale), 
     &                    regStrength)
        trans = gamma_perp*(p02-vvx*p12-vvy*p22)
        regStrength = max(abs(trans)/(Xsi0*MaxPiRatio*pi_scale),
     &                    regStrength)
        trans = gamma_perp*(p00-vvx*p01-vvy*p02)
        regStrength = max(abs(trans)/(Xsi0*MaxPiRatio*pi_scale), 
     &                    regStrength)

        !If (I==0.AND.J==0) Then
        !  Print*, "I,J=",I,J
        !  Print*, "regStrength2=", regStrength
        !End If


        ! largeness comparision
        rTrPi2EAndP = pi_scale/(MaxPiRatio*Tideal_scale) + 1e-30

        regStrength = max(rTrPi2EAndP, regStrength)


        !If (I==0.AND.J==0) Then
        !  Print*, "I,J=",I,J
        !  Print*, "regStrength3=", regStrength
        !End If

        !regStrength = exp(regStrength)-1D0 + 1D-30

        !If (I==0.AND.J==0) Then
        !  Print*, "I,J=",I,J
        !  Print*, "regStrength4=", regStrength
        !End If


        If ( say_level >=9 ) Then
        If (I==II.AND.J==JJ) Then
          Print*, "I,J=",I,J
          Print*, "regStrength=", regStrength
          Print*, "PiTr=", PiTr
          Print*, "numerical zero for pi=", Xsi0*MaxPiRatio*pi_scale
          Print*, "p01-vvx*p11-vvy*p12=", p01-vvx*p11-vvy*p12
          Print*, "regStrength1=",abs(p01-vvx*p11-vvy*p12)
     &                            /(Xsi0*MaxPiRatio*pi_scale)
          Print*, "p02-vvx*p12-vvy*p22=", p02-vvx*p12-vvy*p22
          Print*, "regStrength2=",abs(p02-vvx*p12-vvy*p22)
     &                            /(Xsi0*MaxPiRatio*pi_scale)
          Print*, "p00-vvx*p01-vvy*p02=", p00-vvx*p01-vvy*p02
          Print*, "regStrength3=",abs(p00-vvx*p01-vvy*p02)
     &                            /(Xsi0*MaxPiRatio*pi_scale)
          Print*, "sqrt(TrPi2)=", pi_scale
          Print*, "maxPi=", MaxPiRatio*Tideal_scale
          Print*, "Ed,PL=", Ed(I,J,K), PL(I,J,K)
          Print*, "Tideal_scale=", Tideal_scale
          Print*, "PiTr=", PiTr
          Print*, "regStrength0=", PiTr/(Xsi0*MaxPiRatio*pi_scale)
          Print*, "sqrt(abs(TrPi2)) / maxPi=", rTrPi2EAndP
          Print*, "Xsi0=", Xsi0
!
          call printMore(1, I, J, Time,
     &  NXPhy0, NXPhy, NYPhy0, NYPhy, NX0, NX, NY0, NY, NZ0, NZ,
     &  Pi00,Pi01,Pi02,Pi33,Pi11,Pi12,Pi22,
     &  0D0,0D0,0D0, Ed, PL, 0D0, 0D0, Vx, Vy)
        End If

        If (I==0.AND.J==0) Then
          Print*, "I,J=",I,J
          Print*, "regStrength=", regStrength
          Print*, "PiTr=", PiTr
          Print*, "numerical zero for pi=", Xsi0*MaxPiRatio*pi_scale
          Print*, "p01-vvx*p11-vvy*p12=", p01-vvx*p11-vvy*p12
          Print*, "regStrength1=",abs(p01-vvx*p11-vvy*p12)
     &                            /(Xsi0*MaxPiRatio*pi_scale)
          Print*, "p02-vvx*p12-vvy*p22=", p02-vvx*p12-vvy*p22
          Print*, "regStrength2=",abs(p02-vvx*p12-vvy*p22)
     &                            /(Xsi0*MaxPiRatio*pi_scale)
          Print*, "p00-vvx*p01-vvy*p02=", p00-vvx*p01-vvy*p02
          Print*, "regStrength3=",abs(p00-vvx*p01-vvy*p02)
     &                            /(Xsi0*MaxPiRatio*pi_scale)
          Print*, "sqrt(TrPi2)=", pi_scale
          Print*, "maxPi=", MaxPiRatio*Tideal_scale
          Print*, "Ed,PL=", Ed(I,J,K), PL(I,J,K)
          Print*, "Tideal_scale=", Tideal_scale
          Print*, "PiTr=", PiTr
          Print*, "regStrength0=", PiTr/(Xsi0*MaxPiRatio*pi_scale)
          Print*, "sqrt(abs(TrPi2)) / maxPi=", rTrPi2EAndP
          Print*, "Xsi0=", Xsi0

        End If
        End If !If (debug>=9)


        Pi00(I,J,K)=Pi00(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi01(I,J,K)=Pi01(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi02(I,J,K)=Pi02(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi11(I,J,K)=Pi11(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi12(I,J,K)=Pi12(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi22(I,J,K)=Pi22(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here
        Pi33(I,J,K)=Pi33(I,J,K)*(tanh(regStrength)/regStrength) ! Pi## is regulated here

        PiAvg = 1.0d0/7.0d0*
     &      (abs(Pi00(I,J,K))+abs(Pi01(I,J,K))
     &      +abs(Pi02(I,J,K))
     &      +abs(Pi11(I,J,K))+abs(Pi12(I,J,K))
     &      +abs(Pi22(I,J,K))
     &      +abs(Pi33(I,J,K)))

        If (PiAvg .ne. PiAvg) Then
          Print *, "Invalid PiAvg"
          Print *, "(I,J,K)=",I,J,K
          Print *, "e=", Ed(I,J,K)
          Print *, "p=", PL(I,J,K)
          Print *, "maxPi=", maxPi
          Print *, "TrPi2=",TrPi2
          Print *, "rTrPi2EAndP=",rTrPi2EAndP
          Print *, "Pi00=", Pi00(I,J,K)
          Print *, "Pi01=", Pi01(I,J,K)
          Print *, "Pi02=", Pi02(I,J,K)
          Print *, "Pi11=", Pi11(I,J,K)
          Print *, "Pi12=", Pi12(I,J,K)
          Print *, "Pi22=", Pi22(I,J,K)
          Print *, "Pi33=", Pi33(I,J,K)
          Stop
        EndIf

3008    Continue
3009    Continue

      EndIf ! on regMethod

      If (debug>=3) Print *, "* RegulatePi finished"

      End Subroutine
!-----------------------------------------------------------------------------




!*****************************************************************************
      Subroutine regulateAllPi(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,U0,U1,U2,Time,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,ratio)

      Implicit None

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Integer I,J,K
      Double Precision ratio

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision U0(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision U1(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision U2(NX0:NX, NY0:NY, NZ0:NZ) !local pressure

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision CC(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      CC = 1D0
      CC = min(CC, ratio*abs(Pi00)/max(abs(Ed*U0*U0), 1D-30))
      CC = min(CC, ratio*abs(Pi11)/max(abs((Ed+PL)*U1*U1+PL), 1D-30))
      CC = min(CC, ratio*abs(Pi22)/max(abs((Ed+PL)*U2*U2+PL), 1D-30))
      CC = min(CC, ratio*abs(Pi33)/max(abs(PL), 1D-30))
      CC = min(CC, ratio*abs(Pi01)/max(abs((Ed+PL)*U0*U1), 1D-30))
      CC = min(CC, ratio*abs(Pi02)/max(abs((Ed+PL)*U0*U2), 1D-30))
      CC = min(CC, ratio*abs(Pi12)/max(abs((Ed+PL)*U1*U2), 1D-30))


      Pi00 = CC*Pi00
      Pi11 = CC*Pi11
      Pi22 = CC*Pi22
      Pi33 = CC*Pi33
      Pi01 = CC*Pi01
      Pi02 = CC*Pi02
      Pi12 = CC*Pi12


      End Subroutine
















!----------------------------------------------------------------------------



!*****************************************************************************
      Subroutine doOtherOutputs(NX0,NY0,NZ0,NX,NY,NZ,Ed,PL,Time,
     &  NXPhy0,NXPhy,NYPhy0,NYPhy,
     &  Pi00,Pi01,Pi02,Pi11,Pi12,Pi22,Pi33,
     &  PiPiMaxRatio,maxPiRatio)

      Implicit None

      Integer debug
      Common /debug/ debug

      Integer outPiTrace, outPiStrengthX
      Common /outPi/ outPiTrace, outPiStrengthX

      Integer outPiComponents
      Common /outPiComponents/ outPiComponents

      Integer outEOnXY
      Common /outEOnXY/ outEOnXY

      Integer checkE
      Common /checkE/ checkE

      Double Precision maxE, maxEI, maxEJ, maxEK

      Double Precision outPiLastTime, outPiMovieDt, outPiMovieSize
      common /outPiMovie/ outPiLastTime, outPiMovieDt, outPiMovieSize

      Integer NX0,NY0,NZ0,NX,NY,NZ,NXPhy0,NXPhy,NYPhy0,NYPhy
      Integer I,J,K

      Common /dxdy/ ddx, ddy
      Double Precision ddx, ddy

      Double Precision Time

      Double Precision Ed(NX0:NX, NY0:NY, NZ0:NZ) !energy density
      Double Precision PL(NX0:NX, NY0:NY, NZ0:NZ) !local pressure
      Double Precision Sd(NX0:NX, NY0:NY, NZ0:NZ) !entropy density

      Double Precision Pi00(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi01(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi02(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi11(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi12(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi22(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor
      Double Precision Pi33(NX0:NX, NY0:NY, NZ0:NZ)    !Stress Tensor

      Double Precision EAndP, maxPiRatio, maxPi

      Double Precision PiPiMaxRatio, PiCut ! PiCut=(e+p)*maxPiRatio*PiPiMaxRatio

      Double Precision PiAvg

      Double Precision r00,r01,r02,r11,r12,r22,r33 !radius of region inside which PiXX > (e+p)*maxPiRatio*PiPiMaxRatio

!     Output outPiStrengthX
      If (outPiStrengthX .eq. 1) Then ! Output Pi value on X axis to file
        Open(394,FILE="movie/PiOnX.dat",STATUS='OLD',ACCESS='APPEND') ! file to output where Pi<0.9*max_pi
        Write(394,'(f15.8)',ADVANCE='NO') Time ! write column header
        DO 3108 K=NZ0,NZ0
        DO 3108 J=0,0
        DO 3108 I=0,NXPhy
          PiAvg = 1.0d0/7.0d0*
     &      (abs(Pi00(I,J,K))+abs(Pi01(I,J,K))+abs(Pi02(I,J,K))
     &      +abs(Pi11(I,J,K))+abs(Pi12(I,J,K))+abs(Pi22(I,J,K))
     &      +abs(Pi33(I,J,K)))
          Write(394,'(f15.8)',ADVANCE='NO') PiAvg/(Ed(I,J,K)+PL(I,J,K))
 3108   Continue
        Write(394,*) ! write a new-line-symbol to the file
        Close(394)
      EndIf

!     Output components of Pi
      If (outPiComponents .eq. 1) Then
!       Open files
        Open(190,FILE="movie/PiComTrace.dat",
     &                                  STATUS='OLD',ACCESS='APPEND')
        Open(193,FILE="movie/Pi00X.dat",STATUS='OLD',ACCESS='APPEND')
        Open(194,FILE="movie/Pi01X.dat",STATUS='OLD',ACCESS='APPEND')
        Open(195,FILE="movie/Pi02X.dat",STATUS='OLD',ACCESS='APPEND')
        Open(196,FILE="movie/Pi11X.dat",STATUS='OLD',ACCESS='APPEND')
        Open(197,FILE="movie/Pi12X.dat",STATUS='OLD',ACCESS='APPEND')
        Open(198,FILE="movie/Pi22X.dat",STATUS='OLD',ACCESS='APPEND')
        Open(199,FILE="movie/Pi33X.dat",STATUS='OLD',ACCESS='APPEND')


!       Initialize all the radius
        r00=NX*ddx+NY*ddy
        r01=NX*ddx+NY*ddy
        r02=NX*ddx+NY*ddy
        r11=NX*ddx+NY*ddy
        r12=NX*ddx+NY*ddy
        r22=NX*ddx+NY*ddy
        r33=NX*ddx+NY*ddy

        DO 300 K=NZ0,NZ
        DO 300 J=NY0,NY
        DO 300 I=NX0,NX

        PiCut = PiPiMaxRatio*maxPiRatio*(Ed(I,J,K)+PL(I,J,K)) ! max possible Pi value
        If (abs(Pi00(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r00) Then
            r00 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf
        If (abs(Pi01(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r01) Then
            r01 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf
        If (abs(Pi02(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r02) Then
            r02 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf
        If (abs(Pi11(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r11) Then
            r11 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf
        If (abs(Pi12(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r12) Then
            r12 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf
        If (abs(Pi22(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r22) Then
            r22 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf
        If (abs(Pi33(I,J,K)) > PiCut) Then
          If (sqrt(ddx*ddx*I*I+ddy*ddy*J*J) < r33) Then
            r33 = sqrt(ddx*ddx*I*I+ddy*ddy*J*J)
          EndIf
        EndIf

 300    Continue

        Write(190,"(f15.8,f15.8,f15.8,f15.8,f15.8,f15.8,f15.8,f15.8)")
     &              Time,r00,r01,r02,r11,r12,r22,r33


!       Next, write PiXX/(e+p) along x axis to file
        Write(193,'(f15.8)',ADVANCE='NO') Time ! write column header
        Write(194,'(f15.8)',ADVANCE='NO') Time ! write column header
        Write(195,'(f15.8)',ADVANCE='NO') Time ! write column header
        Write(196,'(f15.8)',ADVANCE='NO') Time ! write column header
        Write(197,'(f15.8)',ADVANCE='NO') Time ! write column header
        Write(198,'(f15.8)',ADVANCE='NO') Time ! write column header
        Write(199,'(f15.8)',ADVANCE='NO') Time ! write column header

        DO 1318 K=NZ0,NZ0
        DO 1318 J=0,0
        DO 1318 I=0,NXPhy
          EAndP = Ed(I,J,K)+PL(I,J,K)
          Write(193,'(f15.8)',ADVANCE='NO') Pi00(I,J,K)/EAndP
          Write(194,'(f15.8)',ADVANCE='NO') Pi01(I,J,K)/EAndP
          Write(195,'(f15.8)',ADVANCE='NO') Pi02(I,J,K)/EAndP
          Write(196,'(f15.8)',ADVANCE='NO') Pi11(I,J,K)/EAndP
          Write(197,'(f15.8)',ADVANCE='NO') Pi12(I,J,K)/EAndP
          Write(198,'(f15.8)',ADVANCE='NO') Pi22(I,J,K)/EAndP
          Write(199,'(f15.8)',ADVANCE='NO') Pi33(I,J,K)/EAndP
 1318   Continue
        Write(193,*) ! write a new-line-symbol to the file
        Write(194,*) ! write a new-line-symbol to the file
        Write(195,*) ! write a new-line-symbol to the file
        Write(196,*) ! write a new-line-symbol to the file
        Write(197,*) ! write a new-line-symbol to the file
        Write(198,*) ! write a new-line-symbol to the file
        Write(199,*) ! write a new-line-symbol to the file

!       Close files
        Close(190)
        Close(193)
        Close(194)
        Close(195)
        Close(196)
        Close(197)
        Close(198)
        Close(199)
      EndIf

!     Output e on x and y axes
      If (outEOnXY .eq. 1) Then
        Open(230,FILE="movie/eOnX.dat",STATUS='OLD',ACCESS='APPEND')
        Open(231,FILE="movie/eOnY.dat",STATUS='OLD',ACCESS='APPEND')
        Write(230,'(f15.8)',ADVANCE='NO') Time
        DO 310 I=0,NXPhy
          Write(230,'(f15.8)',ADVANCE='NO') Ed(I,0,1)
 310    Continue
        Write(231,'(f15.8)',ADVANCE='NO') Time
        DO 320 J=0,NYPhy
          Write(231,'(f15.8)',ADVANCE='NO') Ed(0,J,1)
 320    Continue
        Close(230)
        Close(231)
      EndIf

!     Output maximum value of Ed
      If (checkE .eq. 1) Then
        Open(230,FILE="movie/checkE.dat",STATUS='OLD',ACCESS='APPEND')
        maxE = 0.0;
        DO K=NZ0,NZ
        DO J=NY0,NY
        DO I=NX0,NX
          If (abs(Ed(I,J,K)) > abs(maxE)) Then
            maxE = Ed(I,J,K)
            maxEI = I; maxEJ = J; maxEK = k;
          EndIf
        End Do
        End Do
        End Do
        Write(230,'(f15.8,f15.8,I8,I8,I8)')
     &      Time, maxE, maxEI, maxEJ, maxEK
        Close(230)
        Print *, "Time, maxE, I, J, K",
     &    Time, maxE, maxEI, maxEJ, maxEK
      EndIf

      End Subroutine
