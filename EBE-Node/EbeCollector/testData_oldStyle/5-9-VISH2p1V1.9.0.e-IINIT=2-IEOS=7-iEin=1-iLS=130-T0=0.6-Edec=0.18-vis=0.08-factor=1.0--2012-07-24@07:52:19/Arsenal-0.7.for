************************************************************************
!   Collections of Subroutines related to Table operations
!
!   Author: Zhi Qiu (qiu@mps.ohio-state.edu)
!   -- For changes, see ChangeLogsInputFun.txt
************************************************************************

************************************************************************
      Subroutine interpCubic(table,tableSize,var0,dVar,varX,varResult)
!     Purpose:
!       Return the interpolated value at varX from table, using cubic polynomial interpolation.
!       -- table: 1d array, stores values corresponding to var0+index*dVar
!       -- tableSize: size of the talbe
!       -- var0: start of the independent variable
!       -- dVar: step of the independent variable
!       -- varX: the value of the varibale to be interpolated at
!       -- varResult: the return interpolated value

      Implicit None

!     declare input parameters
      Integer tableSize
      Double Precision table(0:tableSize-1)
      Double Precision var0, dVar, varX, varResult

!     declare local variables
      Integer idx ! varX is between var0+idx*dVar and var0+(idx+1)*dVar
      Double Precision A0,A1,A2,A3,varD ! varD=delta-var

      If ((varX-var0)<dVar*1e-15) Then
        varResult = table(0)
        Return
      End If

      idx = floor((varX-var0)/dVar)

      If (idx<0 .OR. idx>tableSize-1) Then
        Print *, "Subroutine interpCubic: varX out of bounds"
        Print *, "varX=",varX,"var0=",var0,"dVar=",dVar,"idx=",idx
        Stop
      End If

      If (idx==0) Then ! use quadratic interpolation
        A0 = table(0)
        A1 = table(1)
        A2 = table(2)
        varD = varX - var0
        varResult =
     &    (A0-2*A1+A2)/(2*dVar*dVar)*varD*varD
     &  - (3*A0-4*A1+A2)/(2*dVar)*varD
     &  + A0
        Return
      Else If (idx==tableSize-2) Then ! use quadratic interpolation
        A0 = table(tableSize-3)
        A1 = table(tableSize-2)
        A2 = table(tableSize-1)
        varD = varX - (var0 + (idx-1)*dVar)
        varResult =
     &    (A0-2*A1+A2)/(2*dVar*dVar)*varD*varD
     &  - (3*A0-4*A1+A2)/(2*dVar)*varD
     &  + A0
        Return
      Else ! use cubic interpolation
        A0 = table(idx-1)
        A1 = table(idx)
        A2 = table(idx+1)
        A3 = table(idx+2)
        varD = varX - (var0 + idx*dVar)
        varResult =
     &    (-A0+3*A1-3*A2+A3)/(6*dVar*dVar*dVar)*varD*varD*varD
     &  + (A0-2*A1+A2)/(2*dVar*dVar)*varD*varD
     &  - (2*A0+3*A1-6*A2+A3)/(6*dVar)*varD
     &  + A1
        Return
      End If

      End Subroutine
!-----------------------------------------------------------------------



************************************************************************
      Subroutine interpLinear(table,tableSize,var0,dVar,varX,varResult)
!     Purpose:
!       Return the interpolated value at varX from table, using linear interpolation.
!       -- table: 1d array, stores values corresponding to var0+index*dVar
!       -- tableSize: size of the talbe
!       -- var0: start of the independent variable
!       -- dVar: step of the independent variable
!       -- varX: the value of the varibale to be interpolated at
!       -- varResult: the return interpolated value

      Implicit None

!     declare input parameters
      Integer tableSize
      Double Precision table(0:tableSize-1)
      Double Precision var0, dVar, varX, varResult

!     declare local variables
      Integer idx ! varX is between var0+idx*dVar and var0+(idx+1)*dVar
      Double Precision A0,A1,A2,A3,varD ! varD=delta-var

      If ((varX-var0)<dVar*1e-15) Then
        varResult = table(0)
        Return
      End If

      idx = floor((varX-var0)/dVar)

      If (idx<0 .OR. idx>tableSize-1) Then
        Print *, "Subroutine interpLinear: varX out of bounds"
        Print *, "varX=",varX,"var0=",var0,"dVar=",dVar,"idx=",idx
        Stop
      End If

      varD = varX - (var0 + idx*dVar)
      varResult = (table(idx)*(dVar-varD) + table(idx+1)*varD)/dVar

      End Subroutine
!-----------------------------------------------------------------------



************************************************************************
      Subroutine invertNewtonL(table,tableSize,var0,dVar,varX,varResult)
!     Purpose:
!       Return the inverted value at varX from table using Newton method.
!       -- table: 1d array, stores values corresponding to var0+index*dVar
!       -- tableSize: size of the talbe
!       -- var0: start of the independent variable
!       -- dVar: step of the independent variable
!       -- varX: the value of the varibale to be inverted
!       -- varResult: the return inverted value
!
!   Solve: f(x)=0 with f(x)=table(x)-varX => f'(x)=table'(x)
!
      Implicit None

!     declare input parameters
      Integer tableSize
      Double Precision table(0:tableSize-1)
      Double Precision var0, dVar, varX, varResult

!     pre-fixed parameters
      Double Precision accuracy
      Double Precision ending
      Integer edge

!     declare local variables
      Double Precision XX1, XX2 ! used in iterations
      Double Precision F0, F1, F2, F3, F4, X1, X2, tmp ! intermedia variables
      Double Precision impatience ! number of iterations

!     initialize parameters
      accuracy = dVar*1e-3
      ending = var0+tableSize*dVar
      edge = 60

!     initial value, left point and middle point
      XX1 = var0 + tableSize*dVar/2D0
      XX2 = var0

      Do While (abs(XX2-XX1)>accuracy)
        XX1 = XX2 ! copy values

!       value of function at XX!
        Call interpCubic(table,tableSize,var0,dVar,XX1,tmp)
        F0 = tmp  - varX ! the value of the function at this point

!     decide X1 and X2 for differentiation
        If (XX1>var0+dVar) Then
          X1 = XX1 - dVar
        Else
          X1 = var0
        End If

        If (XX1<ending-dVar) Then
          X2 = XX1 + dVar
        Else
          X2 = ending
        End If

!     get values at X1 and X2
        Call interpCubic(table,tableSize,var0,dVar,X1,tmp)
        F1 = tmp
        Call interpCubic(table,tableSize,var0,dVar,X2,tmp)
        F2 = tmp

        F3 = (F1-F2)/(X1-X2) ! derivative at XX1

        XX2 = XX1 - F0/F3 ! Newton's mysterious method

        impatience = impatience + 1
        If (impatience>edge) Then
          Print *, "Subroutine invertNewton:",
     &                "max number of iterations reached!"
          Stop
        End If
      End Do ! <=> abs(XX2-XX1)>accuracy

      varResult = XX2

      End Subroutine
!-----------------------------------------------------------------------



************************************************************************
      Subroutine invertFunctionD(func,varL,varR,dd,varI,varX,varResult)
!     Purpose:
!       Return the varResult=func^-1(varX) using Newton method.
!       -- func: double precision 1-argument function to be inverted
!       -- varL: left boundary (for numeric derivative)
!       -- varR: right boundary (for numeric derivative)
!       -- dd: step (for numeric derivative)
!       -- varI: initial value
!       -- varX: the value of the varibale to be inverted
!       -- varResult: the return inverted value
!
!   Solve: f(x)=0 with f(x)=table(x)-varX => f'(x)=table'(x)
!
      Implicit None

!     declare input parameters
      Double Precision func
      Double Precision varL, varR, dd, varI, varX, varResult

!     pre-fixed parameters
      Double Precision accuracy
      Integer tolerance

!     declare local variables
      Double Precision XX1, XX2 ! used in iterations
      Double Precision F0, F1, F2, F3, X1, X2, tmp ! intermedia variables
      Integer impatience ! number of iterations

!     initialize parameters
      accuracy = dd*1e-3

      tolerance = 60
      impatience = 0

!     initial value, left point and middle point
      XX2 = varI
      XX1 = XX2-10*accuracy ! this value 10*accuracy is meanless, just to make sure the check in the while statement goes through

      Do While (abs(XX2-XX1)>accuracy)
        XX1 = XX2 ! copy values

!       value of function at XX!
        F0 = func(XX1)  - varX ! the value of the function at this point

!     decide X1 and X2 for differentiation
        If (XX1>varL+dd) Then
          X1 = XX1 - dd
        Else
          X1 = varL
        End If

        If (XX1<varR-dd) Then
          X2 = XX1 + dd
        Else
          X2 = varR
        End If

!     get values at X1 and X2
        F1 = func(X1)
        F2 = func(X2)
        F3 = (F1-F2)/(X1-X2) ! derivative at XX1

        XX2 = XX1 - F0/F3 ! Newton's mysterious method

        impatience = impatience + 1
        If (impatience>tolerance) Then
          Print *, "Subroutine invertFunctionD: ",
     &                "max number of iterations reached!"
          Stop
          !Print*, XX1, XX2
        End If
      End Do ! <=> abs(XX2-XX1)>accuracy

      varResult = XX2

      End Subroutine
!-----------------------------------------------------------------------


************************************************************************
      Subroutine QShasta2(U,vX,vY,S,NX0,NX,NY0,NY,NZ0,NZ,
     &  DT,DX,DY, NXPhy0,NYPhy0, NXPhy,NYPhy, AntiX,AntiY)
!     Purpose:
!       Return U(t=n+1) using 2D Shasta algorithm.
!       -- U: The quantitiy to be propogated
!       -- S: Source term

      Implicit None

!     declare input parameters
      Integer NX0,NX,NY0,NY,NZ0,NZ
      Integer NXPhy0,NYPhy0, NXPhy,NYPhy
      Double Precision DT,DX,DY, antiX, antiY
      Double Precision U(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision vX(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision vY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision S(NX0:NX,NY0:NY,NZ0:NZ)

!     declare local variables
      Double Precision UtX(NX0:NX,NY0:NY,NZ0:NZ) ! \tilde{U}^x
      Double Precision UtY(NX0:NX,NY0:NY,NZ0:NZ) ! \tilde{U}^y
      Double Precision Ut(NX0:NX,NY0:NY,NZ0:NZ) ! \tilde{U}
      Double Precision DelX(NX0:NX,NY0:NY,NZ0:NZ) ! \Delta^x
      Double Precision DelY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision QPX(NX0:NX,NY0:NY,NZ0:NZ) ! Q_{+}^x
      Double Precision QMX(NX0:NX,NY0:NY,NZ0:NZ) ! Q_{-}^x
      Double Precision QPY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision QMY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision ilX, ilY ! 1/\lambda^x
      Double Precision DelTX(NX0:NX,NY0:NY,NZ0:NZ) ! \tilde{\Delta}^x
      Double Precision DelTY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision AX(NX0:NX,NY0:NY,NZ0:NZ) ! A_{ad}^X
      Double Precision AY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision sigmaX(NX0:NX,NY0:NY,NZ0:NZ) ! \sigma^x
      Double Precision sigmaY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision ATX(NX0:NX,NY0:NY,NZ0:NZ) ! \tilde{A}^x
      Double Precision ATY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision UTMin(NX0:NX,NY0:NY,NZ0:NZ) ! \tilde{U}^{min}
      Double Precision UTMax(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision AIn(NX0:NX,NY0:NY,NZ0:NZ) ! A^{in}
      Double Precision AOut(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision FIn(NX0:NX,NY0:NY,NZ0:NZ) ! F^{in}
      Double Precision FOut(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision PFilter(NX0:NX,NY0:NY,NZ0:NZ) ! positive filter
      Double Precision NFilter(NX0:NX,NY0:NY,NZ0:NZ) ! negative filter
      Double Precision AHatX(NX0:NX,NY0:NY,NZ0:NZ) ! \hat{A}^x
      Double Precision AHatY(NX0:NX,NY0:NY,NZ0:NZ)
      Double Precision Q1(NX0:NX,NY0:NY,NZ0:NZ) ! tmp var

      Common /R0Bdry/ R0Bdry
      Double Precision R0Bdry

      Integer I,J,K
      Double Precision xx,yy,rr,ff

!     initialization to zeros
      UtX = 0D0
      UtY = 0D0
      Ut = 0D0
      DelX = 0D0
      DelY = 0D0
      QPX = 0.5D0
      QMX = 0.5D0
      QPY = 0.5D0
      QMY = 0.5D0
      ilX = 0D0
      ilY = 0D0
      DelTX = 0D0
      DelTY = 0D0
      UTMin = 0D0
      UTMax = 0D0
      AIn = 1D-30
      AOut = 1D-30
      FIn = 0D0
      FOut = 0D0
      AHatX = 0D0
      AHatY = 0D0

!     1/\lambda^x
      ilX = DX/DT ! should <=0.5
      ilY = DY/DT

!     prepare x direction
      DelX(NX0:NX-1,:,:) = U(NX0+1:NX,:,:)-U(NX0:NX-1,:,:)
      QPX(NX0:NX-1,:,:) = (ilX/2D0-vX(NX0:NX-1,:,:))
     &                  / (ilX+(vX(NX0+1:NX,:,:)-vX(NX0:NX-1,:,:)))
      QMX(NX0+1:NX,:,:) = (ilX/2D0+vX(NX0+1:NX,:,:))
     &                  / (ilX-(vX(NX0:NX-1,:,:)-vX(NX0+1:NX,:,:)))
      UtX(NX0+1:NX,:,:) = 0.5D0*(
     &   QPX(NX0+1:NX,:,:)**2*DelX(NX0+1:NX,:,:)
     & - QMX(NX0+1:NX,:,:)**2*DelX(NX0:NX-1,:,:))
     & + (QPX(NX0+1:NX,:,:)+QMX(NX0+1:NX,:,:))*U(NX0+1:NX,:,:)

!     prepare y direction
      DelY(:,NY0:NY-1,:) = U(:,NY0+1:NY,:)-U(:,NY0:NY-1,:)
      QPY(:,NY0:NY-1,:) = (ilY/2D0-vY(:,NY0:NY-1,:))
     &                  / (ilY+(vY(:,NY0+1:NY,:)-vY(:,NY0:NY-1,:)))
      QMY(:,NY0+1:NY,:) = (ilY/2D0+vY(:,NY0+1:NY,:))
     &                  / (ilY-(vY(:,NY0:NY-1,:)-vY(:,NY0+1:NY,:)))
      UtY(:,NY0+1:NY,:) = 0.5D0*(
     &   QPY(:,NY0+1:NY,:)**2*DelY(:,NY0+1:NY,:)
     & - QMY(:,NY0+1:NY,:)**2*DelY(:,NY0:NY-1,:))
     & + (QPY(:,NY0+1:NY,:)+QMY(:,NY0+1:NY,:))*U(:,NY0+1:NY,:)

!     1st order propagating
      Ut = UtX + UtY - U + DT*S

!     anti-diffusing preparation
      DelTX(NX0:NX-1,:,:) = Ut(NX0+1:NX,:,:) - Ut(NX0:NX-1,:,:)
      DelTY(:,NY0:NY-1,:) = Ut(:,NY0+1:NY,:) - Ut(:,NY0:NY-1,:)
      AX = antiX*DelTX/8D0
      AY = antiY*DelTY/8D0
      sigmaX = sign(1D0, AX)
      sigmaY = sign(1D0, AY)
      ATX(NX0+1:NX-1,:,:) = sigmaX(NX0+1:NX-1,:,:)*max(0D0, min(
     &      sigmaX(NX0+1:NX-1,:,:)*DelTX(NX0+2:NX,:,:),
     &      abs(AX(NX0+1:NX-1,:,:)),
     &      sigmaX(NX0+1:NX-1,:,:)*DelTX(NX0:NX-2,:,:)))
      ATY(:,NY0+1:NY-1,:) = sigmaY(:,NY0+1:NY-1,:)*max(0D0, min(
     &      sigmaY(:,NY0+1:NY-1,:)*DelTY(:,NY0+2:NY,:),
     &      abs(AY(:,NY0+1:NY-1,:)),
     &      sigmaY(:,NY0+1:NY-1,:)*DelTY(:,NY0:NY-2,:)))


!     more anti-diffusing preparation
      UTMin(NX0+1:NX-1,NY0+1:NY-1,:) = min(
     &    Ut(NX0+1:NX-1,NY0:NY-2,:), Ut(NX0:NX-2,NY0+1:NY-1,:),
     &    Ut(NX0+1:NX-1,NY0+1:NY-1,:), Ut(NX0+2:NX,NY0+1:NY-1,:),
     &    Ut(NX0+1:NX-1,NY0+2:NY,:))
      UTMax(NX0+1:NX-1,NY0+1:NY-1,:) = max(
     &    Ut(NX0+1:NX-1,NY0:NY-2,:), Ut(NX0:NX-2,NY0+1:NY-1,:),
     &    Ut(NX0+1:NX-1,NY0+1:NY-1,:), Ut(NX0+2:NX,NY0+1:NY-1,:),
     &    Ut(NX0+1:NX-1,NY0+2:NY,:))
      AIn(NX0+1:NX,NY0+1:NY,:) =
     &   max(0D0, ATX(NX0:NX-1,NY0+1:NY,:))
     & - min(0D0, ATX(NX0+1:NX,NY0+1:NY,:))
     & + max(0D0, ATY(NX0+1:NX,NY0:NY-1,:))
     & - min(0D0, ATY(NX0+1:NX,NY0+1:NY,:))
      AOut(NX0+1:NX,NY0+1:NY,:) =
     &   max(0D0, ATX(NX0+1:NX,NY0+1:NY,:))
     & - min(0D0, ATX(NX0:NX-1,NY0+1:NY,:))
     & + max(0D0, ATY(NX0+1:NX,NY0+1:NY,:))
     & - min(0D0, ATY(NX0+1:NX,NY0:NY-1,:))
      AIn = max(1D-30, AIn)
      AOut = max(1D-30, AOut)
      FIn = (UTMax - Ut) / AIn
      FOut = (Ut - UTMin) / AOut

      Q1 = sign(1D0,ATX) ! sign
      PFilter = (Q1 + abs(ATX))/2D0 ! positive part only
      NFilter = PFilter - Q1 ! negative part only
      AHatX(NX0:NX-1,:,:) = ATX(NX0:NX-1,:,:) * (
     &   min(1D0, FIn(NX0+1:NX,:,:), FOut(NX0:NX-1,:,:))*PFilter
     & + min(1D0, FIn(NX0:NX-1,:,:), FOut(NX0+1:NX,:,:))*NFilter )
      Q1 = sign(1D0,ATY) ! sign
      PFilter = (Q1 + abs(ATY))/2D0 ! positive part only
      NFilter = PFilter - Q1 ! negative part only
      AHatY(:,NY0:NY-1,:) = ATY(:,NY0:NY-1,:) * (
     &   min(1D0, FIn(:,NY0+1:NY,:), FOut(:,NY0:NY-1,:))*PFilter
     & + min(1D0, FIn(:,NY0:NY-1,:), FOut(:,NY0+1:NY,:))*NFilter )

!     anti-diffusing!!
      Ut(NX0+1:NX,NY0+1:NY,:) = Ut(NX0+1:NX,NY0+1:NY,:)
     &  - AHatX(NX0+1:NX,NY0+1:NY,:) - AHatY(NX0+1:NX,NY0+1:NY,:)
     &  + AHatX(NX0:NX-1,NY0+1:NY,:) + AHatY(NX0+1:NX,NY0:NY-1,:)

!     smearing
      !DO K=NZ0,NZ
      !DO J=NY0,NY
      !DO I=NX0,NX
      !
      !xx=DX*I
      !yy=DY*J
      !rr=sqrt(xx**2+yy**2)
      !ff=1.0/(Dexp((rr-R0Bdry)/0.4D0)+1.0)
      !
      !U(I,J,K) = ff*U(I,J,K)
      !
      !EndDo
      !EndDo
      !EndDo


!     boundary killer
      Ut(NX0:NXPhy0+2,:,:) = 0D0
      Ut(NXPhy+2:NX,:,:) = 0D0
      Ut(:,NY0:NYPhy0+2,:) = 0D0
      Ut(:,NYPhy+2:NY,:) = 0D0

!     return value
      U = Ut

      End Subroutine
!-----------------------------------------------------------------------
