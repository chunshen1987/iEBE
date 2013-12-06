!###########################################################################
! Version 1.02
! -- For changes, see ChangeLogsPhyBDary.txt
!---------------------------------------------------------------------------

C###########################################################################
      Subroutine UPShasta2(TT,Vx,Vy,ScT, NX0,NX,NY0,NY,NZ0,NZ,  IPX,IPY,
     &  DT,DX,DY,  NXPhy0,NYPhy0, NXPhy,NYPhy, ISYM1,ISYM2, DIFFC )
C------ using P.Colb 2+1 Shastala to cal the tansport part
C------ Transport diffusion part:
C-------limiter (antidiffusion part)
C-----  IPX,IPY ! differential x(or y) for source term
C------ NXPhy, NYphy :Real Physics demension < NX <NY
C------ ISYM1,ISYM2 symetric boundary condition for 2+1
        Implicit Double Precision (A-H, O-Z)

        Dimension TT(NX0:NX, NY0:NY, NZ0:NZ)  !Tansport density
        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)  !Velocity in X
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)  !Velocity in Y
        Dimension ScT(NX0:NX, NY0:NY, NZ0:NZ) !Source Term  ScT=ScX+ScY+ScZ

        Dimension PTT(NX0:NX, NY0:NY)
        Dimension PVx(NX0:NX, NY0:NY)
        Dimension PVy(NX0:NX, NY0:NY)
        Dimension PSc(NX0:NX, NY0:NY)
        DIMENSION WORK1(NX0:NX,NY0:NY),WORK2(NX0:NX,NY0:NY)

C        if(NX0.ne.-2.or.NY0.ne.-2.or.NX.ne.MAX .or. NY.ne.MAY ) then
C          print 'in UPShasta wrong dimension'
C          stop
C        end if

         do 100 I=NX0,NX   !(NX0,NX)
         do 100 J=NY0,NY   !(NY0,NY)
           PTT(I,J)=TT(I,J,NZ0)     !TT(I+3, J+3)
           PVx(I,J)=Vx(I,J,NZ0)
           PVy(I,J)=Vy(I,J,NZ0)
           PSc(I,J)=ScT(I,J,NZ0)
 100      continue

        CALL SHATRA(NXPhy0,NYPhy0,NXPhy,NYPhy,PTT,PVx,PVy,PSc, IPX, IPY,
     &         WORK1,  ISYM1,ISYM2, DT/DX,DT/DY,  NX0, NY0,  NX,NY)
        CALL SHACOR(NXPhy0,NYPhy0,NXPhy,NYPhy,PTT,PVx,PVy, WORK1,
     &    DT/DX,DT/DY, WORK2,WORK1,DIFFC, ISYM1,ISYM2, NX0, NY0,  NX,NY)

          If(IPX.eq.0.and.IPY.eq.0 ) then
C           DO 500 J=NYPhy0-2,NYPhy+2
C           DO 500 I=NXPhy0-2,NXPhy+2
           DO 500 J=NYPhy0,NYPhy
           DO 500 I=NXPhy0,NXPhy
           PTT(I,J)=PTT(I,J)-PSc(I,J)*DT
 500       Continue
          end if

         do 800 I=NX0,NX   !(NX0,NX)
         do 800 J=NY0,NY   !(NY0,NY)
           TT(I,J,NZ0)=PTT(I,J)     !TT(I+3, J+3)
 800      continue

      Return
      End

C##################################################################
*************************************************************

      SUBROUTINE SHATRA(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,B,ISL1,ISL2,AB,
     &        ISYM1,ISYM2,  TDX,TDY, NX0, NY0, NX,NY)
*     In order to get the symmetry relations right we need again
*     a routine to introduce the right conditions, especially for
*     the y-Component of the velocity, that is T02
        Implicit Double Precision (A-H, O-Z)
      DIMENSION A(NX0:NX,NY0:NY),U(NX0:NX,NY0:NY),
     &          V(NX0:NX,NY0:NY),B(NX0:NX,NY0:NY),
     &          AB(NX0:NX,NY0:NY)


      CALL SHAS1(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,B,ISL1,ISL2,AB,
     &                        TDX,TDY, NX0, NY0, NX,NY)

*      Ab  work1
*     Now AB is the transported and diffused field that will
*     replace A
*     Care must still be taken on the boundary conditions

      DO I=NXPhy0-1,NXPhy+1
cCC        AB(I,-1)=ISYM2*AB(I,1)
cCC        AB(I,-2)=ISYM2*AB(I,2)
!        AB(I,NYPhy+2) = 0.0
!        AB(I,NYPhy0-2) = 0.0
        AB(I,NYPhy+2) =2.0*AB(I,NYPhy+1)-AB(I,NYPhy)
        AB(I,NYPhy0-2)=2.0*AB(I,NYPhy0-1)-AB(I,NYPhy0)
      END DO

      DO J=NYPhy0-2,NYPhy+2
CCC        AB(-1,J)=ISYM1*AB(1,J)
CCC        AB(-2,J)=ISYM1*AB(2,J)
!        AB(NXPhy+2,J) = 0.0
!        AB(NXPhy0-2,J) = 0.0
        AB(NXPhy+2,J) =2.0*AB(NXPhy+1,J)-AB(NXPhy,J)
        AB(NXPhy0-2,J)=2.0*AB(NXPhy0-1,J)-AB(NXPhy0,J)
      END DO

      DO I=NXPhy0-2,NXPhy+2
        AB(I,NYPhy+2)=(2.0*AB(I,NYPhy+1)-AB(I,NYPhy))*0.5
     &                                      +0.5*AB(I,NYPhy+2)
        AB(I,NYPhy0-2)=(2.0*AB(I,NYPhy0-1)-AB(I,NYPhy0))*0.5
     &                                      +0.5*AB(I,NYPhy0-2)
      END DO

cCC        AB(-1,-1)=ISYM1*ISYM2*AB(1,1)
cCC        AB(-2,-1)=ISYM1*ISYM2*AB(2,1)
cCC        AB(-1,-2)=ISYM1*ISYM2*AB(1,2)
cCC        AB(-2,-2)=ISYM1*ISYM2*AB(2,2)

      RETURN
      END

*************************************************************

      SUBROUTINE SHACOR(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,AB,TDX,TDY,
     &                  AC1,GAB,DIFFC,ISYM1,ISYM2, NX0,NY0,NX,NY)
        Implicit Double Precision (A-H, O-Z)
*************************************************************
**  THIS SUBROUTINE CALLS THE ANTIDIFFUSION ALGORITHM TO    *
** DO THE CORRECTION ON THE TRANSPORTED AND DIFFUSED        *
** SOLUTION OF THE DIFFERENTIAL EQUATION                    *
**        DA/DT=-D(AU)/DX-D(AV)/DY-DB/DX                    *
**      USING 'SHASTA-MULTIDIMENSIONAL FCT' ALGORITHM       *
** ONE CALL OF THE ROUTINE ADVANCES THE DEPENDENT VARIABLE  *
** A BY ONE TIME STEP. THIS IS DONE IN TWO STEPS: THE FIRST *
** (SHAS1) IS THE 'TRANSPORT AND DIFFUSION', AND THE        *
** SECOND (SHAS2) IS THE 'ANTIDIFFUSION' STEP.              *
** AFTER EACH STEP THE BOUNDARIES MUST BE HANDLED           *
** SEPARATELY ACCORDING TO BOUNDARY CONDITIONS OR SYMMETRY. *
*************************************************************
** BOUNDARY IN THIS VERSION:                                *
** -SOLUTION IS ASSUMED TO BE SYMMETRIC OR ANTI-            *
**  SYMMETRIC ACCORDING TO PARAMETERS 'ISYM%' (=-1,0 OR 1)  *
**  WITH RESPECT TO X=X(MIN) (I=0).                         *
** -SOLUTION IS ASSUMED TO BE CONST. NEAR X=X(MAX).         *
*************************************************************

      DIMENSION A(NX0:NX,NY0:NY),U(NX0:NX,NY0:NY),
     &          V(NX0:NX,NY0:NY),
     &          AB(NX0:NX,NY0:NY),AC1(NX0:NX,NY0:NY),
     &          GAB(NX0:NX,NY0:NY)

      CALL SHAS2(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,AB,AC1,GAB,DIFFC,
     &                       TDX,TDY,NX0,NY0,NX,NY)

      DO 180 I=NXPhy0-2,NXPhy+2
        A(I,NYPhy+1)=0.0
        A(I,NYPhy+2)=0.0
        A(I,NYPhy0-1)=0.0
        A(I,NYPhy0-2)=0.0
CCC        A(I,-1)= ISYM2*A(I,1)
CCC        A(I,-2)= ISYM2*A(I,2)
180   CONTINUE

      DO 190 J=NYPhy0-2,NYPhy+2
        A(NXPhy+1,J)=0.0
        A(NXPhy+2,J)=0.0
        A(NXPhy0-1,J)=0.0
        A(NXPhy0-2,J)=0.0
CCC        A(-1,J)=ISYM1*A(1,J)
CCC        A(-2,J)=ISYM2*A(2,J)
190   CONTINUE

cCC      A(-1,-1)=ISYM1*ISYM2*A(1,1)
CCC      A(-2,-1)=ISYM1*ISYM2*A(2,1)
CCC      A(-1,-2)=ISYM1*ISYM2*A(1,2)
CCC      A(-2,-2)=ISYM1*ISYM2*A(2,2)

      RETURN
      END

*********************************************************

      SUBROUTINE SHAS1(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,B,ISL1,ISL2,AB,
     &                 TDX,TDY,NX0,NY0,NX,NY)
        Implicit Double Precision (A-H, O-Z)
*********************************************************
** THIS ROUTINE PRODUCES THE 'TRANSPORTED AND DIFFUSED' *
** SOLUTION AB FOR THE SHASTA ALGORITHM.                *
** NOTE: FOR INPUT: A(I,J),U(I,J),V(I,J),B(I,J)         *
**       FOR OUTPUT: AB(I,J)                            *
**            I=-1-N,...,N+1; J=-1-N,...,N+1            *
** NOTE: 'ISL' MUST BE EITHER 1 OR 0. (THIS IS NOT      *
**       CHECKED HERE!)                                 *
** NOTE: IF U(I,J)*DT/DX >0.5                           *
**  OR   IF V(I,J)*DT/DY >0.5 FOR SOME I, EXECUTION IS  *
**       TERMINATED AND AN ERROR MESSAGE PRINTED        *
*********************************************************

      DIMENSION A(NX0:NX,NY0:NY),U(NX0:NX,NY0:NY),
     &          V(NX0:NX,NY0:NY),B(NX0:NX,NY0:NY),
     &          AB(NX0:NX,NY0:NY)

          DO 100 I=NXPhy0-1,NXPhy+1
          DO 120 J=NYPhy0-1,NYPhy+1

          EIM1=U(I-1,J)*TDX
          EII=U(I,J)*TDX
          EIP1=U(I+1,J)*TDX
          EJM1=V(I,J-1)*TDY
          EJI=V(I,J)*TDY
          EJP1=V(I,J+1)*TDY

          IF (EJM1.GE.0.5) THEN
          WRITE(*,*) 'SHAS1: DT/DY TOO LARGE ! J=',J-1,'  I=',I
          WRITE(*,*) 'TDX=',TDX,'  V(I,J-1)=',V(I,J-1)
          STOP
          END IF
          IF (EIM1.GE.0.5) THEN
          WRITE(*,*) 'SHAS1: DT/DX TOO LARGE ! I=',I-1,'   J=',J
          WRITE(*,*) 'TDX=',TDX,'  U(I-1,J)=',U(I-1,J)
          STOP
          END IF


          QUP=(0.5-EII)/(1.0+(EIP1-EII))
          QUM=(0.5+EII)/(1.0-(EIM1-EII))
          QVP=(0.5-EJI)/(1.0+(EJP1-EJI))
          QVM=(0.5+EJI)/(1.0-(EJM1-EJI))
          AIM1=A(I,J)-A(I-1,J)
          AJM1=A(I,J)-A(I,J-1)
          AI=A(I+1,J)-A(I,J)
          AJ=A(I,J+1)-A(I,J)
          BIM1=ISL1*TDX*(B(I,J)-B(I-1,J)) !ISL1
          BJM1=ISL2*TDY*(B(I,J)-B(I,J-1)) !ISL2
          BI=ISL1*TDX*(B(I+1,J)-B(I,J))
          BJ=ISL2*TDY*(B(I,J+1)-B(I,J))


* This is the transportation step given by Boris and Book
* The twice appearing -0.5*A(I,J) is taking the 2 spacial directions
* into account

*          AB(I,J)=0.5*(QUP*QUP*AI-QUM*QUM*AIM1)
*     &             +QUP*(A(I,J)-BI)+QUM*(A(I,J)-BIM1)
*     &             -0.5*A(I,J)
*     &        +     0.5*(QVP*QVP*AJ-QVM*QVM*AJM1)
*     &             +QVP*(A(I,J)-BJ)+QVM*(A(I,J)-BJM1)
*     &             -0.5*A(I,J)

* The idea for the following transportation step is given
* in Rischkes Paper in Nucl Ph A595(1995)346


           AB(I,J)=0.5*(QUP*QUP*AI-QUM*QUM*AIM1)
     &             +QUP*A(I,J)+QUM*A(I,J)
     &             -0.5*ISL1*TDX*(B(I+1,J)-B(I-1,J))
     &             -0.5*A(I,J)
     &        +     0.5*(QVP*QVP*AJ-QVM*QVM*AJM1)
     &             +QVP*A(I,J)+QVM*A(I,J)
     &             -0.5*ISL2*TDY*(B(I,J+1)-B(I,J-1) )
     &             -0.5*A(I,J)


120    CONTINUE
100   CONTINUE


       IF (EJP1.GE.0.5) THEN
       WRITE(*,*) 'SHAS1: DT/DY TOO LARGE ! J=',J+1
       STOP
       END IF
       IF (EIP1.GE.0.5) THEN
       WRITE(*,*) 'SHAS1: DT/DX TOO LARGE ! I=',I+1
       STOP
       END IF

      RETURN
      END

******************************************************

      SUBROUTINE SHAS2(NXPhy0,NYPhy0,NXPhy,NYPhy,A,U,V,AB,AC1,
     &                    GAB,DIFFC,TDX,TDY,NX0,NY0,NX,NY)
        Implicit Double Precision (A-H, O-Z)
******************************************************
** THIS ROUTINE PRODUCES THE 'ANTIDIFFUSED' SOLUTION *
** A FROM THE 'TRANSPORTED AND DIFFUSED' SOLUTION AB *
** FOR THE SHASTA-FCT ALGORITHM USING ZALESAK:S      *
** METHOD                                            *
** NOTE: FOR INPUT: AB(I,J),U(I,J),V(I,J)            *
**       FOR OUTPUT: A(I,J)   I= 0,..,N; J= 0,..,N   *
** NOTE: IN THIS VERSION THE ANTIDIFFUSION CONSTANT  *
**       'DIFFC' IS FIXED NUMBER GIVEN TO ROUTINE    *
**       AS A PARAMETER.(SO U AND V ARE NOT NEEDED.) *
******************************************************
      DIMENSION AB(NX0:NX,NY0:NY),A(NX0:NX,NY0:NY),
     &           U(NX0:NX,NY0:NY),V(NX0:NX,NY0:NY),
     &         AC1(NX0:NX,NY0:NY),GAB(NX0:NX,NY0:NY)

      COMMON /PARAM/ DX, DY, DT, X0, Y0, ALPHA

      T1 = 0.0
      ACON = 0.0

      DO I=NXPhy0-2,NXPhy+2
       DO J=NYPhy0-2,NYPhy+2
       AC1(I,J)=A(I,J)
       END DO
      END DO



      DO I=NXPhy0,NXPhy
      IP=I-1
      IPP=I+1
      DO J=NYPhy0,NYPhy
        YY=0
        YPY1=0
        YMY1=0
      JP=J-1
      JPP=J+1

*     Here are some remnants to check the influence of the antidiffusion

      DIFF = DIFFC

*      EPSX = ABS(0.5*(U(IP+1,J)+U(IP,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPIP=DIFF*(AB(IP+1,J)-AB(IP,J))

*      EPSX = ABS(0.5*(U(IP,J)+U(IP-1,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPIM=DIFF*(AB(IP,J)-AB(IP-1,J))

*      EPSY = ABS(0.5*(V(IP,J+1)+V(IP,J)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPJP=DIFF*(AB(IP,J+1)-AB(IP,J))

*      EPSY = ABS(0.5*(V(IP,J)+V(IP,J-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPJM=DIFF*(AB(IP,J)-AB(IP,J-1))

*      EPSX = ABS(0.5*(U(I+1,J)+U(I,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIJIP=DIFF*(AB(I+1,J)-AB(I,J))

*      EPSX = ABS(0.5*(U(I,J)+U(I-1,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIJIM=DIFF*(AB(I,J)-AB(I-1,J))

*      EPSY = ABS(0.5*(V(I,J+1)+V(I,J)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIJJP=DIFF*(AB(I,J+1)-AB(I,J))

*      EPSY = ABS(0.5*(V(I,J)+V(I,J-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIJJM=DIFF*(AB(I,J)-AB(I,J-1))

*      EPSX = ABS(0.5*(U(IPP+1,J)+U(IPP,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPPIP=DIFF*(AB(IPP+1,J)-AB(IPP,J))

*      EPSX =ABS( 0.5*(U(IPP,J)+U(IPP-1,J)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AIPPIM=DIFF*(AB(IPP,J)-AB(IPP-1,J))

*      EPSY = ABS(0.5*(V(IPP,J+1)+V(IPP,J)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPPJP=DIFF*(AB(IPP,J+1)-AB(IPP,J))

*      EPSY = ABS(0.5*(V(IPP,J)+V(IPP,J-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AIPPJM=DIFF*(AB(IPP,J)-AB(IPP,J-1))

*      EPSX = ABS(0.5*(U(I+1,JP)+U(I,JP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPIP=DIFF*(AB(I+1,JP)-AB(I,JP))

*      EPSX = ABS(0.5*(U(I,JP)+U(I-1,JP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPIM=DIFF*(AB(I,JP)-AB(I-1,JP))

*      EPSY = ABS(0.5*(V(I,JP+1)+V(I,JP)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPJP=DIFF*(AB(I,JP+1)-AB(I,JP))

*     EPSY = ABS(0.5*(V(I,JP)+V(I,JP-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPJM=DIFF*(AB(I,JP)-AB(I,JP-1))

*      EPSX = ABS(0.5*(U(I+1,JPP)+U(I,JPP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPPIP=DIFF*(AB(I+1,JPP)-AB(I,JPP))

*      EPSX = ABS(0.5*(U(I,JPP)+U(I-1,JPP)))*TDX
*      DIFF = DIFFC + EPSX*EPSX*ACON
      AJPPIM=DIFF*(AB(I,JPP)-AB(I-1,JPP))

*      EPSY = ABS(0.5*(V(I,JPP+1)+V(I,JPP)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPPJP=DIFF*(AB(I,JPP+1)-AB(I,JPP))

*      EPSY = ABS(0.5*(V(I,JPP)+V(I,JPP-1)))*TDY
*      DIFF = DIFFC + EPSY*EPSY*ACON
      AJPPJM=DIFF*(AB(I,JPP)-AB(I,JPP-1))


C **** SEE PG.349 EQ-14 IN ZALESAK FOR THE BELOW CONDITIONS ********

      PRPI=AIJIP*(AB(I+1,J)-AB(I,J))
      PRPIPP=AIJIP*(AB(I+2,J)-AB(I+1,J))
      PRPIP=AIJIP*(AB(I,J)-AB(I-1,J))
      PRMI=AIJIM*(AB(I+1,J)-AB(I,J))
      PRMIP=AIJIM*(AB(I,J)-AB(I-1,J))
      PRMIP1=AIJIM*(AB(I-1,J)-AB(I-2,J))

      PRPJ=AIJJP*(AB(I,J+1)-AB(I,J))
      PRPJPP=AIJJP*(AB(I,J+2)-AB(I,J+1))
      PRPJP=AIJJP*(AB(I,J)-AB(I,J-1))
      PRMJ=AIJJM*(AB(I,J+1)-AB(I,J))
      PRMJP=AIJJM*(AB(I,J)-AB(I,J-1))
      PRMJP1=AIJJM*(AB(I,J-1)-AB(I,J-2))



      UIPJ=dMAX1(AC1(IP,J),GAB(IP,J))
      UIP1J=dMAX1(AC1(IP-1,J),GAB(IP-1,J))
      UIJ=dMAX1(AC1(I,J),GAB(I,J))
      UIPPJ=dMAX1(AC1(IPP,J),GAB(IPP,J))
      UIPP1J=dMAX1(AC1(IPP+1,J),GAB(IPP+1,J))
      UIJP=dMAX1(AC1(I,JP),GAB(I,JP))
      UIJP1=dMAX1(AC1(I,JP-1),GAB(I,JP-1))
      UIJPP=dMAX1(AC1(I,JPP),GAB(I,JPP))
      UIJPP1=dMAX1(AC1(I,JPP+1),GAB(I,JPP+1))
      UIPJP=dMAX1(AC1(IP,JP),GAB(IP,JP))
      UIPJPP=dMAX1(AC1(IP,JPP),GAB(IP,JPP))
      UIPPJP=dMAX1(AC1(IPP,JP),GAB(IPP,JP))
      UIPPJPP=dMAX1(AC1(IPP,JPP),GAB(IPP,JPP))

      VIPJ=dMIN1(AC1(IP,J),GAB(IP,J))
      VIP1J=dMIN1(AC1(IP-1,J),GAB(IP-1,J))
      VIJ=dMIN1(AC1(I,J),GAB(I,J))
      VIPPJ=dMIN1(AC1(IPP,J),GAB(IPP,J))
      VIPP1J=dMIN1(AC1(IPP+1,J),GAB(IPP+1,J))
      VIJP=dMIN1(AC1(I,JP),GAB(I,JP))
      VIJP1=dMIN1(AC1(I,JP-1),GAB(I,JP-1))
      VIJPP=dMIN1(AC1(I,JPP),GAB(I,JPP))
      VIJPP1=dMIN1(AC1(I,JPP+1),GAB(I,JPP+1))
      VIPJP=dMIN1(AC1(IP,JP),GAB(IP,JP))
      VIPJPP=dMIN1(AC1(IP,JPP),GAB(IP,JPP))
      VIPPJP=dMIN1(AC1(IPP,JP),GAB(IPP,JP))
      VIPPJPP=dMIN1(AC1(IPP,JPP),GAB(IPP,JPP))

      WMXIPJ=dMAX1(UIP1J,UIPJ,UIJ,UIPJP,UIPJPP)
      WMXIJ=dMAX1(UIPJ,UIJ,UIPPJ,UIJP,UIJPP)
      WMXIPPJ=dMAX1(UIJ,UIPPJ,UIPP1J,UIPPJP,UIPPJPP)
      WMXIJP=dMAX1(UIPJP,UIJP,UIPPJP,UIJP1,UIJ)
      WMXIJPP=dMAX1(UIPJPP,UIJPP,UIPPJPP,UIJ,UIJPP1)



      WMNIPJ=dMIN1(VIP1J,VIPJ,VIJ,VIPJP,VIPJPP)
      WMNIJ=dMIN1(VIPJ,VIJ,VIPPJ,VIJP,VIJPP)
      WMNIPPJ=dMIN1(VIJ,VIPPJ,VIPP1J,VIPPJP,VIPPJPP)
      WMNIJP=dMIN1(VIPJP,VIJP,VIPPJP,VIJP1,VIJ)
      WMNIJPP=dMIN1(VIPJPP,VIJPP,VIPPJPP,VIJ,VIJPP1)

      PIPJP=dMAX1(0.0d0,AIPIM)-dMIN1(0.0d0,AIPIP)
     &     +dMAX1(0.0d0,AIPJM)-dMIN1(0.0d0,AIPJP)
      QIPJP=(WMXIPJ-GAB(IP,J))
      IF (PIPJP.GT.0.0) THEN
        RIPJP=dMIN1(1.0d0,QIPJP/PIPJP)
      ELSE
        RIPJP=0.0
      END IF

      PIJP=dMAX1(0.0d0,AIJIM)-dMIN1(0.0d0,AIJIP)
     &    +dMAX1(0.0d0,AIJJM)-dMIN1(0.0d0,AIJJP)
      QIJP=(WMXIJ-GAB(I,J))
      IF (PIJP.GT.0.0) THEN
        RIJP=dMIN1(1.0d0,QIJP/PIJP)
      ELSE
        RIJP=0.0
      END IF

      PIPPJP=dMAX1(0.0d0,AIPPIM)-dMIN1(0.0d0,AIPPIP)
     &      +dMAX1(0.0d0,AIPPJM)-dMIN1(0.0d0,AIPPJP)
      QIPPJP=(WMXIPPJ-GAB(IPP,J))
      IF (PIPPJP.GT.0.0) THEN
        RIPPJP=dMIN1(1.0d0,QIPPJP/PIPPJP)
      ELSE
        RIPPJP=0.0
      END IF

      PIJPP=dMAX1(0.0d0,AJPIM)-dMIN1(0.0d0,AJPIP)
     &     +dMAX1(0.0d0,AJPJM)-dMIN1(0.0d0,AJPJP)
      QIJPP=(WMXIJP-GAB(I,JP))
      IF (PIJPP.GT.0.0) THEN
        RIJPP=dMIN1(1.0d0,QIJPP/PIJPP)
      ELSE
        RIJPP=0.0
      END IF

      PIJPPP=dMAX1(0.0d0,AJPPIM)-dMIN1(0.0d0,AJPPIP)
     &      +dMAX1(0.0d0,AJPPJM)-dMIN1(0.0d0,AJPPJP)
      QIJPPP=(WMXIJPP-GAB(I,JPP))
      IF (PIJPPP.GT.0.0) THEN
        RIJPPP=dMIN1(1.0d0,QIJPPP/PIJPPP)
      ELSE
        RIJPPP=0.0
      END IF

      PIPJM=dMAX1(0.0d0,AIPIP)-dMIN1(0.0d0,AIPIM)
     &     +dMAX1(0.0d0,AIPJP)-dMIN1(0.0d0,AIPJM)
      QIPJM=(GAB(IP,J)-WMNIPJ)
      IF (PIPJM.GT.0.0) THEN
        RIPJM=dMIN1(1.0d0,QIPJM/PIPJM)
      ELSE
        RIPJM=0.0
      END IF

      PIJM=DMAX1(0.0d0,AIJIP)-dMIN1(0.0d0,AIJIM)
     &    +DMAX1(0.0d0,AIJJP)-DMIN1(0.0d0,AIJJM)
      QIJM=(GAB(I,J)-WMNIJ)
      IF (PIJM.GT.0.0) THEN
        RIJM=DMIN1(1.0D0,QIJM/PIJM)
      ELSE
        RIJM=0.0
      END IF

      PIPPJM=DMAX1(0.0D0,AIPPIP)-DMIN1(0.0D0,AIPPIM)
     &      +DMAX1(0.0D0,AIPPJP)-dMIN1(0.0D0,AIPPJM)
      QIPPJM=(GAB(IPP,J)-WMNIPPJ)
      IF (PIPPJM.GT.0.0) THEN
        RIPPJM=DMIN1(1.0D0,QIPPJM/PIPPJM)
      ELSE
        RIPPJM=0.0
      END IF

      PIJPM=DMAX1(0.0D0,AJPIP)-DMIN1(0.0D0,AJPIM)
     &     +DMAX1(0.0D0,AJPJP)-DMIN1(0.0D0,AJPJM)
      QIJPM=(GAB(I,JP)-WMNIJP)
      IF (PIJPM.GT.0.0) THEN
        RIJPM=DMIN1(1.0D0,QIJPM/PIJPM)
      ELSE
        RIJPM=0.0
      END IF

      PIJPPM=DMAX1(0.0D0,AJPPIP)-DMIN1(0.0D0,AJPPIM)
     &      +DMAX1(0.0D0,AJPPJP)-DMIN1(0.0D0,AJPPJM)
      QIJPPM=(GAB(I,JPP)-WMNIJPP)
      IF (PIJPPM.GT.0.0) THEN
        RIJPPM=DMIN1(1.0D0,QIJPPM/PIJPPM)
      ELSE
        RIJPPM=0.0
      END IF

      IF (AIJIP.GE.0.0) THEN
       CIJIP=DMIN1(RIPPJP,RIJM)
      ELSE
       CIJIP=DMIN1(RIJP,RIPPJM)
      END IF

      IF (AIJIM.GE.0.0) THEN
       CIJIM=DMIN1(RIJP,RIPJM)
      ELSE
       CIJIM=DMIN1(RIPJP,RIJM)
      END IF

      IF (AIJJP.GE.0.0) THEN
       CIJJP=DMIN1(RIJPPP,RIJM)
      ELSE
       CIJJP=DMIN1(RIJP,RIJPPM)
      END IF

      IF (AIJJM.GE.0.0) THEN
       CIJJM=DMIN1(RIJP,RIJPM)
      ELSE
       CIJJM=DMIN1(RIJPP,RIJM)
      END IF


      ACIJIP=CIJIP*AIJIP
      ACIJJP=CIJJP*AIJJP
      test1= 0.5*(0.125-CIJJP*DIFFC)*YPY1*SIJJP
*      ACIJJP=CIJJP*AIJJP-0.5*(0.125-CIJJP*DIFFC)*YPY1*SIJJP
      ACIJIM=CIJIM*AIJIM
      ACIJJM=CIJJM*AIJJM
      test2= 0.5*(0.125-CIJJM*DIFFC)*YMY1*SIJJM
*      ACIJJM=CIJJM*AIJJM+0.5*(0.125-CIJJM*DIFFC)*YMY1*SIJJM

*      IF (test1.NE.0.) WRITE(*,*) test1
*      IF (test2.NE.0.) WRITE(*,*) test2

      A(I,J)=AB(I,J)-(ACIJIP-ACIJIM
     &         +ACIJJP-ACIJJM)


      END DO
      END DO

      RETURN
      END

C################################################################
       Double Precision FUNCTION  EOUT(BA, EN, STAN,
     &      BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
     &      BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,MAT1,MAT2)

       Implicit Double Precision (A-H, O-Z)
       double precision  MAT1(0:300,0:250), MAT2(0:300,0:250)
       double precision  MATRIX

      IF (EN .LT. EPP01 + DEPP1) THEN
        FX1 = STAN
        FY1 = MAT1(0,1)
        FY2 = MAT1(1,1)
        X1 = BNP01
        X2 = BNP01 + DBNP1
        Y1 = 0.0
        Y2 = EPP01 + DEPP1
        EOUT = RTER3D(X1, X2, Y1, Y2, FX1, FY1, FY2, BA, EN)
      ELSE
        IF (EN .LT. EPP01 + NEPP1*DEPP1) THEN
          EOUT = MATRIX(BNP01, DBNP1, NBNP1,
     &              EPP01, DEPP1, STAN, MAT1, BA, EN)
        ELSE
          IF (EN .LE. EPP02 + NEPP2*DEPP2) THEN
          EOUT = MATRIX(BNP02, DBNP2, NBNP2,
     &              EPP02, DEPP2, STAN, MAT2, BA, EN)
          ELSE
            FX1 = MAT2(NBNP2-2,NEPP2-2)
            FX2 = MAT2(NBNP2  ,NEPP2-2)
            FY1 = MAT2(NBNP2-2,NEPP2  )
            FY2 = MAT2(NBNP2  ,NEPP2  )
            X1 = BNP02 + (NBNP2-2)*DBNP2
            X2 = BNP02 + NBNP2*DBNP2
            Y1 = EPP02 + (NEPP2-2)*DEPP2
            Y2 = EPP02 + NEPP2*DEPP2
            EOUT = RTER2(X1, X2, Y1, Y2, FX1, FX2, FY1, FY2, BA, EN)

          END IF
        END IF
      END IF

      RETURN
      END
***************************************************************
      SUBROUTINE POINT(BA,EN,PRE)
******************************************************************
*  THIS ROUTINE INTERPOLATES PRESSURE VALUES (PRE) IN POINT      *
*  (BNP,EPP) USING MATRICES PMAT1 AND PMAT2.                     *
******************************************************************
      Implicit Double Precision (A-H, O-Z)
      COMMON /PAREOS/ BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
     &                BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,
     &                PMAT1(0:300,0:250),PMAT2(0:300,0:250)

      Double Precision  MATRIX

      STAN = 0.0
      IF (EN .LT. EPP01 + DEPP1) THEN
        FX1 = STAN
        FY1 = PMAT1(0,1)
        FY2 = PMAT1(1,1)
        X1 = BNP01
        X2 = BNP01 + DBNP1
        Y1 = 0.0
        Y2 = EPP01 + DEPP1
        PRE = RTER3D(X1, X2, Y1, Y2, FX1, FY1, FY2, BA, EN)
      ELSE
        IF (EN .LT. EPP01 + NEPP1*DEPP1) THEN
          PRE = MATRIX(BNP01,DBNP1,NBNP1,
     &              EPP01,DEPP1,STAN,PMAT1, BA, EN)
        ELSE
          IF (EN .LE. EPP02 + NEPP2*DEPP2) THEN
            PRE = MATRIX(BNP02,DBNP2,NBNP2,
     &              EPP02,DEPP2,STAN,PMAT2, BA, EN)
          ELSE
            FX1 = PMAT2(NBNP2-2,NEPP2-2)
            FX2 = PMAT2(NBNP2  ,NEPP2-2)
            FY1 = PMAT2(NBNP2-2,NEPP2  )
            FY2 = PMAT2(NBNP2  ,NEPP2  )
            X1 = BNP02 + (NBNP2-2)*DBNP2
            X2 = BNP02 + NBNP2*DBNP2
            Y1 = EPP02 + (NEPP2-2)*DEPP2
            Y2 = EPP02 + NEPP2*DEPP2
            PRE = RTER2(X1, X2, Y1, Y2, FX1, FX2, FY1, FY2, BA, EN)
          END IF
        END IF
      END IF
      RETURN
      END


***************************************************************

       Double Precision FUNCTION  MATRIX(BNP0, DBNP, NBNP,
     &     EPP0, DEPP, STAN, MAT, BA, EN)

       Implicit Double Precision (A-H, O-Z)
       Double Precision  MAT(0:300,0:250)

      IF(BA .GT.  BNP0 + NBNP*DBNP) THEN
        EPPY = (EN-EPP0)/DEPP
      J = INT(EPPY)
            FX1 = MAT(NBNP-1,J)
      FX2 = MAT(NBNP,J)
      FY1 = MAT(NBNP-1,J+1)
      FY2 = MAT(NBNP,J+1)
            X1 = BNP0 + (NBNP-1)*DBNP
      X2 = BNP0 + NBNP*DBNP
            Y1 = EPP0 + J * DEPP
      Y2 = EPP0 + (J+1) * DEPP
        MATRIX = RTER2(X1, X2, Y1, Y2, FX1, FX2, FY1, FY2, BA, EN)
        RETURN
      END IF

      BNPX = (BA-BNP0)/DBNP
      EPPY = (EN-EPP0)/DEPP
      I = INT(BNPX)
      J = INT(EPPY)
      FX1 = MAT(I,J)
      FX2 = MAT(I+1,J)
      FY1 = MAT(I,J+1)
      FY2 = MAT(I+1,J+1)

      IF ( (FX2 .LT. 0.0) .AND. (FY2 .LT. 0.0) ) THEN

       I = I - 1
        FX1 = MAT(I,J)
        FX2 = MAT(I+1,J)
        FY1 = MAT(I,J+1)
        FY2 = MAT(I+1,J+1)
      END IF

      X1 = BNP0 + I*DBNP
      X2 = BNP0 + (I+1)*DBNP
      Y1 = EPP0 + J * DEPP
      Y2 = EPP0 + (J+1) * DEPP
      IF ( (FX1 .LT. 0.0) .OR. (FY1 .LT. 0.0) ) THEN
       WRITE(*,*) 'BA ,EN not in physical region ! BA = ',BA,' EN=',EN
       MATRIX = STAN
      END IF

      IF ((FX2 .GE. 0.0) .AND. (FY2 .LT. 0.0) )THEN
        MATRIX = RTER3U(X1, X2, Y1, Y2, FX1, FX2, FY1, BA, EN)
      ELSE IF ((FX2 .LT. 0.0) .AND. (FY2 .GE. 0.0) )THEN
        MATRIX = RTER3D(X1, X2, Y1, Y2, FX1, FY1, FY2, BA, EN)
      ELSE
        MATRIX = RTER2(X1, X2, Y1, Y2, FX1, FX2, FY1, FY2, BA, EN)
      END IF

      RETURN
      END

*****************************************************************
       Double Precision FUNCTION RTER3D(X0, X1, Y0, Y1,
     &                           F00, F01, F11, X, Y)
       Implicit Double Precision (A-H, O-Z)
        DX = X1 - X0
        DY = Y1 - Y0
        P = (X-X0)/DX
        Q = (Y-Y0)/DY
        RTER3D = (1.0-Q)*F00 + (Q-P)*F01 + P*F11
        RETURN
      END
********************************************************************
       Double Precision  FUNCTION RTER2(X0, X1, Y0, Y1,
     &                        F00, F10, F01, F11, X, Y)
       Implicit Double Precision (A-H, O-Z)
        DX = X1 - X0
        DY = Y1 - Y0
        P = (X-X0)/DX
        Q = (Y-Y0)/DY
        RTER2 = (1.0-P)*(1.0-Q)*F00 + P*(1.0-Q)*F10 +
     &              Q*(1.0-P)*F01 + P*Q*F11
        RETURN
      END
*****************************************************************
       Double Precision  FUNCTION RTER3U(X0, X1, Y0, Y1,
     &                              F00, F10, F01, X, Y)
       Implicit Double Precision (A-H, O-Z)
       WRITE(*,*) 'RTER3U wird also doch noch gebraucht'
        DX = X1 - X0
        DY = Y1 - Y0
        P = (X-X0)/DX
        Q = (Y-Y0)/DY
        RTER3U = (1.0 - P - Q)*F00 + P*F10 + Q*F01
        RETURN
      END







C######################################################################
      SUBROUTINE SECTIO(E0,INTERSECT,EPCUB,EPS0,EPS1,
     &                                        I,J,NDX,NDY,MAX,MAY)

*****************************************************************
** THIS ROUTINE CHECKS WHETHER THE eps=E0 -SURFACE INTERSECTS   *
** WITH A 'CUBE' (I,J). IF THAT IS THE CASE, THE LOGICAL        *
** VARIABLE 'INTERSECT' IS RETURNED WITH A VALUE .TRUE. AND     *
** THE MATRIX EPCUB(8) WITH THE eps VALUES OF THE CUBE CORNERS. *
****************************************************************
       Implicit Double Precision (A-H, O-Z)
      DIMENSION EPS0(-2:MAX,-2:MAY),EPS1(-2:MAX,-2:MAY)
      DIMENSION EPCUB(8)
      LOGICAL INTERSECT

      IF ((E0-EPS0(I,J))*(EPS1(I+NDX,J+NDY)-E0).LT.0.0) THEN
        IF ((E0-EPS0(I+NDX,J))*(EPS1(I,J+NDY)-E0).LT.0.0) THEN
          IF ((E0-EPS0(I+NDX,J+NDY))*(EPS1(I,J)-E0).LT.0.0) THEN
            IF ((E0-EPS0(I,J+NDY))*(EPS1(I+NDX,J)-E0).LT.0.0) THEN
              INTERSECT=.FALSE.
              GO TO 200
            END IF
          END IF
        END IF
      END IF

100   INTERSECT=.TRUE.
      EPCUB(1)=EPS0(I,J)
      EPCUB(2)=EPS0(I+NDX,J)
      EPCUB(3)=EPS0(I+NDX,J+NDY)
      EPCUB(4)=EPS0(I,J+NDY)
      EPCUB(5)=EPS1(I,J)
      EPCUB(6)=EPS1(I+NDX,J)
      EPCUB(7)=EPS1(I+NDX,J+NDY)
      EPCUB(8)=EPS1(I,J+NDY)

200   RETURN
      END



*******************************************************************
      SUBROUTINE P4(I,J,NDX,NDY,NDT,VMID,F0,F1,
     &    NX0,NY0,NX,NY,DT,DX,DY,F,IDebug)

************************************************************
** THIS ROUTINE PERFORMS A THREE DIMENSIONAL INTERPOLATION
** USING DATA IN MATRICES F0 AND F1.
** Assumes NDX,NDY,NDT are all positive
************************************************************
      Implicit Double Precision (A-H, O-Z)
      Double Precision F0(NX0:NX,NY0:NY),F1(NX0:NX,NY0:NY)
      Double Precision V000,V100,V010,V001,V101,V011,V110,V111
      Double Precision VMID(0:2)
      Double Precision F

      Double Precision x,y,z

      !Print *, "P4 Started"
      !Print *, "I,J=",I,J
      !Print *, "NDX,NDY=",NDX,NDY

      x=VMID(1)/DX
      y=VMID(2)/DY
      z=VMID(0)/DT

      !Print *, "Division finished"

      V000=F0(I,J)
      V100=F0(I+NDX,J)
      V010=F0(I,J+NDY)
      V110=F0(I+NDX,J+NDY)

      !Print *, "Set variables halfway finished"

      V001=F1(I,J)
      V101=F1(I+NDX,J)
      V011=F1(I,J+NDY)
      V111=F1(I+NDX,J+NDY)

      !Print *, "Set variables finished"


      F = V000*(1-x)*(1-y)*(1-z) + V100*x*(1-y)*(1-z) +
     &    V010*(1-x)*y*(1-z) + V001*(1-x)*(1-y)*z +
     &    V101*x*(1-y)*z + V011*(1-x)*y*z +
     &    V110*x*y*(1-z) + V111*x*y*z

      If (IDebug == 1) Then
        Print *, "P4 debug"
        Print *, "x,y,z=",x,y,z
        Print *, "V000,V100,V010,V110=",V000,V100,V010,V110
        Print *, "V001,V101,V011,V111=",V001,V101,V011,V111
        Print *, "F=", F
        Print *, "VMID=", VMID
        Print *, "DX,DY,DT=", DX,DY,DT
        Print *, "NDX,NDY,NDT=", NDX,NDY,NDT
      End If

      RETURN
      END




*******************************************************************
       SUBROUTINE CUBCAL(E0,EPCUB,SURFS,NSURFS,VMID,IBIT,DT,DX,DY,
     &                  NERR,MAXERR)
        Implicit Double Precision (A-H, O-Z)
************************************************************
** THIS ROUTINE CALCULATES THE INTERSECTION OF AN 'EPS=E0'
** SURFACE WITH A CUBE WHICH IS DETERMINED BY GIVING ENERGY
** VALUES IN THE EIGHT CORNERS OF IT (VECTOR 'EPCUB')
** ROUTINE DIVIDES THIS INTERSECTION INTO TRIANGLES AND
** GIVES AS ITS OUTPUT THE COVARIANT SURFACE ELEMENT
** VECTORS ('SURFS') CORRESPONDING TO THESE TRIANGLES.
** NOTE: AN APPROXIMATION FOR THE TOTAL SURFACE ELEMENT
**       VECTOR CORRESPONDING TO THE PORTION OF THE SURFACE
**       INTERSECTING THE CUBE, IS GIVEN SIMPLY AS A SUM OF
**       THESE VECTORS!
************************************************************

      DIMENSION EPCUB(8),SURFS(12,0:2),CUTS(12,0:2)
      DIMENSION VMID(0:2)
      DIMENSION AD(0:2),BD(0:2)
      DIMENSION IBIT(12,12),IEDGE(12),ISID(6)


******************************************************
*** FIND THE INTERSECTION POINTS OF THE SURFACE WITH *
*** THE 12 EDGES OF THE CUBE. CALCULATE THE WEIGHTED *
*** AVERAGE OF THE VECTORS POINTING TO CUBE CORNERS  *
*** WITH EPS<E0 AND EPS>E0 SEPARATELY.               *
******************************************************
C       Print*,'Cubcal'
      NSURFS=0
      NLO=0
      NHI=0
      VL0=0.0
      VL1=0.0
      VL2=0.0
      VH0=0.0
      VH1=0.0
      VH2=0.0
      ELSUM=0.0
      EHSUM=0.0

      EK=EPCUB(1)
      EL=EPCUB(2)
      DEK=(E0-EK)
      ADEK=ABS(DEK)
      IF (DEK*(EL-E0).GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=1
        CUTS(NSURFS,0)=0.0
        CUTS(NSURFS,1)=(E0-EK)/(EL-EK)*DX
        CUTS(NSURFS,2)=0.0
      END IF
      IF (DEK.GT.0.0) THEN
        ELSUM=ELSUM+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
      END IF

      EK=EPCUB(2)
      EL=EPCUB(3)
      DEK=(E0-EK)
      ADEK=ABS(DEK)
      IF (DEK*(EL-E0).GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=2
        CUTS(NSURFS,0)=0.0
        CUTS(NSURFS,1)=DX
        CUTS(NSURFS,2)=(E0-EK)/(EL-EK)*DY
      END IF
      IF (DEK.GT.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL1=VL1+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH1=VH1+ADEK
      END IF

      EK=EPCUB(4)
      EL=EPCUB(3)
      DEK=(EL-E0)
      ADEK=ABS(DEK)
      IF ((E0-EK)*DEK.GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=3
        CUTS(NSURFS,0)=0.0
        CUTS(NSURFS,1)=(E0-EK)/(EL-EK)*DX
        CUTS(NSURFS,2)=DY
      END IF
      IF (DEK.LE.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL1=VL1+ADEK
        VL2=VL2+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH1=VH1+ADEK
        VH2=VH2+ADEK
      END IF

      EK=EPCUB(1)
      EL=EPCUB(4)
      DEK=(EL-E0)
      ADEK=ABS(DEK)
      IF ((E0-EK)*DEK.GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=4
        CUTS(NSURFS,0)=0.0
        CUTS(NSURFS,1)=0.0
        CUTS(NSURFS,2)=(E0-EK)/(EL-EK)*DY
      END IF
      IF (DEK.LE.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL2=VL2+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH2=VH2+ADEK
      END IF

      EK=EPCUB(1)
      EL=EPCUB(5)
      DEK=(EL-E0)
      ADEK=ABS(DEK)
      IF ((E0-EK)*DEK.GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=5
        CUTS(NSURFS,0)=(E0-EK)/(EL-EK)*DT
        CUTS(NSURFS,1)=0.0
        CUTS(NSURFS,2)=0.0
      END IF
      IF (DEK.LE.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL0=VL0+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH0=VH0+ADEK
      END IF

      EK=EPCUB(2)
      EL=EPCUB(6)
      DEK=(EL-E0)
      ADEK=ABS(DEK)
      IF ((E0-EK)*DEK.GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=6
        CUTS(NSURFS,0)=(E0-EK)/(EL-EK)*DT
        CUTS(NSURFS,1)=DX
        CUTS(NSURFS,2)=0.0
      END IF
      IF (DEK.LE.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL0=VL0+ADEK
        VL1=VL1+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH0=VH0+ADEK
        VH1=VH1+ADEK
      END IF

      EK=EPCUB(3)
      EL=EPCUB(7)
      DEK=(EL-E0)
      ADEK=ABS(DEK)
      IF ((E0-EK)*DEK.GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=7
        CUTS(NSURFS,0)=(E0-EK)/(EL-EK)*DT
        CUTS(NSURFS,1)=DX
        CUTS(NSURFS,2)=DY
      END IF
      IF (DEK.LE.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL0=VL0+ADEK
        VL1=VL1+ADEK
        VL2=VL2+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH0=VH0+ADEK
        VH1=VH1+ADEK
        VH2=VH2+ADEK
      END IF

      EK=EPCUB(4)
      EL=EPCUB(8)
      DEK=(EL-E0)
      ADEK=ABS(DEK)
      IF ((E0-EK)*DEK.GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=8
        CUTS(NSURFS,0)=(E0-EK)/(EL-EK)*DT
        CUTS(NSURFS,1)=0.0
        CUTS(NSURFS,2)=DY
      END IF
      IF (DEK.LE.0.0) THEN
        ELSUM=ELSUM+ADEK
        VL0=VL0+ADEK
        VL2=VL2+ADEK
      ELSE
        EHSUM=EHSUM+ADEK
        VH0=VH0+ADEK
        VH2=VH2+ADEK
      END IF

      EK=EPCUB(5)
      EL=EPCUB(6)
      IF ((E0-EK)*(EL-E0).GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=9
        CUTS(NSURFS,0)=DT
        CUTS(NSURFS,1)=(E0-EK)/(EL-EK)*DX
        CUTS(NSURFS,2)=0.0
      END IF

      EK=EPCUB(6)
      EL=EPCUB(7)
      IF ((E0-EK)*(EL-E0).GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=10
        CUTS(NSURFS,0)=DT
        CUTS(NSURFS,1)=DX
        CUTS(NSURFS,2)=(E0-EK)/(EL-EK)*DY
      END IF

      EK=EPCUB(8)
      EL=EPCUB(7)
      IF ((E0-EK)*(EL-E0).GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=11
        CUTS(NSURFS,0)=DT
        CUTS(NSURFS,1)=(E0-EK)/(EL-EK)*DX
        CUTS(NSURFS,2)=DY
      END IF

      EK=EPCUB(5)
      EL=EPCUB(8)
      IF ((E0-EK)*(EL-E0).GE.0.0) THEN
        NSURFS=NSURFS+1
        IEDGE(NSURFS)=12
        CUTS(NSURFS,0)=DT
        CUTS(NSURFS,1)=0.0
        CUTS(NSURFS,2)=(E0-EK)/(EL-EK)*DY
      END IF

*****************************************************************
** VDN:S CALCULATED BELOW ARE THE COMPONENTS OF THE DIFFERENCE **
** OF WEIGHTED AVERAGES OF POINTS WITH E<E0 AND E>E0.          **
** THIS VECTOR IS USED TO DEFINE THE POSITIVE SIDE OF SURFACE  **
** ELEMENTS IN THIS CUBE TO BE TOWARDS LOWER ENERGY.           **
*****************************************************************

      IF(ELSUM.EQ.0.0)GOTO 133
      VL0=VL0/ELSUM*DT
      VL1=VL1/ELSUM*DX
      VL2=VL2/ELSUM*DY

133   IF(EHSUM.EQ.0.0)GOTO 135
      VH0=VH0/EHSUM*DT
      VH1=VH1/EHSUM*DX
      VH2=VH2/EHSUM*DY

135   VD0=VL0-VH0
      VD1=VL1-VH1
      VD2=VL2-VH2

***   END INTERSECTION POINTS ***

*********************************************************
*** SORT THE INTERSECTION POINTS INTO A CIRCULAR ORDER **
*** USING A BITCHART IN MATRIX 'IBIT'.                 **
*********************************************************


      DO 100 I=1,NSURFS-2
        IE=IEDGE(I)
        IS=0
        DO 110 J=I+1,NSURFS
          JE=IEDGE(J)
          IF (IBIT(IE,JE).NE.0) THEN
            IS=IS+1
            ISID(IS)=J
          END IF
110     CONTINUE

        IF (IS.NE.1) THEN
          IF (IS.EQ.0) THEN
            GO TO 500
          END IF
          MINPTS=100000
          DO 120 J=1,IS
            JE=IEDGE(ISID(J))
            IPTS=0
            DO 130 K=1,I-1
              KE=IEDGE(K)
              IPTS=IPTS+10*MIN0(1,IBIT(KE,JE))
130         CONTINUE
            IPTS=IPTS+IBIT(JE,IE)-1
            IF (IPTS.LT.MINPTS) THEN
              JMIN=ISID(J)
              MINPTS=IPTS
            END IF
120       CONTINUE
        ELSE
          JMIN=ISID(1)
        END IF

        DO 140 K=0,2
          APU=CUTS(JMIN,K)
          CUTS(JMIN,K)=CUTS(I+1,K)
          CUTS(I+1,K)=APU
140     CONTINUE
        IE=IEDGE(JMIN)
        IEDGE(JMIN)=IEDGE(I+1)
        IEDGE(I+1)=IE
100   CONTINUE

      NSE=IEDGE(NSURFS)
      NSM=IEDGE(NSURFS-1)
      I1=IEDGE(1)
      IF (IBIT(NSE,NSM).EQ.0.OR.IBIT(NSE,I1).EQ.0) THEN
        GO TO 500
      END IF

****  END SORT ***

******************************************************
** CALCULATE THE MEAN VECTOR OF INTERSECTION POINTS  *
******************************************************

      V0=0.0
      V1=0.0
      V2=0.0

      DO 300 N=1,NSURFS
        V0=V0+CUTS(N,0)
        V1=V1+CUTS(N,1)
        V2=V2+CUTS(N,2)
300   CONTINUE

      VMID(0)=V0/NSURFS
      VMID(1)=V1/NSURFS
      VMID(2)=V2/NSURFS

      DO 400 N=1,NSURFS
        M=MOD(N,NSURFS)+1

        DO 410 I=0,2
          AD(I)=CUTS(M,I)-CUTS(N,I)
          BD(I)=VMID(I)-CUTS(N,I)
410     CONTINUE

****************************************************
*** COVARIANT COMPONENTS OF THE SURFACE VECTOR 'L'**
****************************************************

        SU0L=.5*(AD(2)*BD(1)-AD(1)*BD(2))
        SU1L=.5*(AD(0)*BD(2)-AD(2)*BD(0))
        SU2L=.5*(AD(1)*BD(0)-AD(0)*BD(1))

***************************************************************
*** CHOOSE THE POSITIVE DIRECTION TO BE TOWARDS LOWER ENERGY **
*** BY GIVING CORRECT SIGN TO SURFACE VECTOR.                **
***************************************************************

        VSUM=VD0*SU0L+VD1*SU1L+VD2*SU2L
        SIG=SIGN(1.0d0,VSUM)

        SURFS(N,0)=SIG*SU0L
        SURFS(N,1)=SIG*SU1L
        SURFS(N,2)=SIG*SU2L

400   CONTINUE

      RETURN

500   NERR=NERR+1
520   CONTINUE

      SURFS(1,0)=0.0
      SURFS(1,1)=0.0
      SURFS(1,2)=0.0
      NSURFS=0

      RETURN
      END





*******************************************************************
      SUBROUTINE PP3(I,J,NDX,NDY,VMID,F0,F1,NX0,NY0,NX,NY,DT,DX,DY,F)

************************************************************
** THIS ROUTINE PERFORMS A THREE DIMENSIONAL INTERPOLATION
** USING DATA IN MATRICES F0 AND F1.
************************************************************
        Implicit Double Precision (A-H, O-Z)
      DIMENSION F0(NX0:NX, NY0:NY),F1(NX0:NX, NY0:NY)
      DIMENSION VMID(0:2)

      A=VMID(1)/DX
      B=VMID(2)/DY
      C=VMID(0)/DT
      IP=I+NDX
      JP=J+NDY
      F0IPJP=F0(IP,JP)
      F0IPJ=F0(IP,J)
      F0IJP=F0(I,JP)
      F1IPJP=F1(IP,JP)
      F1IPJ=F1(IP,J)
      F1IJP=F1(I,JP)

      IF ((A+B).LT.1.0) THEN
        AB1=1-A-B
        F=(AB1*F0(I,J)+A*F0IPJ+B*F0IJP)*(1-C)
        F=F+(AB1*F1(I,J)+A*F1IPJ+B*F1IJP)*C
      ELSE
        A=A-1.0
        B=B-1.0
        AB1=1+A+B
        F=(AB1*F0IPJP-A*F0IJP-B*F0IPJ)*(1-C)
        F=F+(AB1*F1IPJP-A*F1IJP-B*F1IPJ)*C
      END IF

c      Print*,'pp3'
      RETURN
      END

*******************************************************************
      SUBROUTINE SECTIO2(E0,INTERSECT,EPCUB,EPS0,EPS1,
     &                        I,J,NDX,NDY,NX0,NY0,NX,NY)

*****************************************************************
** THIS ROUTINE CHECKS WHETHER THE eps=E0 -SURFACE INTERSECTS   *
** WITH A 'CUBE' (I,J). IF THAT IS THE CASE, THE LOGICAL        *
** VARIABLE 'INTERSECT' IS RETURNED WITH A VALUE .TRUE. AND     *
** THE MATRIX EPCUB(8) WITH THE eps VALUES OF THE CUBE CORNERS. *
****************************************************************
      Implicit Double Precision (A-H, O-Z)
      DIMENSION EPS0(NX0:NX, NY0:NY),EPS1(NX0:NX, NY0:NY)
      DIMENSION EPCUB(8)
      LOGICAL INTERSECT

      IF ((E0-EPS0(I,J))*(EPS1(I+NDX,J+NDY)-E0).LT.0.0) THEN
        IF ((E0-EPS0(I+NDX,J))*(EPS1(I,J+NDY)-E0).LT.0.0) THEN
          IF ((E0-EPS0(I+NDX,J+NDY))*(EPS1(I,J)-E0).LT.0.0) THEN
            IF ((E0-EPS0(I,J+NDY))*(EPS1(I+NDX,J)-E0).LT.0.0) THEN
              INTERSECT=.FALSE.
              GO TO 200
            END IF
          END IF
        END IF
      END IF

100   INTERSECT=.TRUE.
      EPCUB(1)=EPS0(I,J)
      EPCUB(2)=EPS0(I+NDX,J)
      EPCUB(3)=EPS0(I+NDX,J+NDY)
      EPCUB(4)=EPS0(I,J+NDY)
      EPCUB(5)=EPS1(I,J)
      EPCUB(6)=EPS1(I+NDX,J)
      EPCUB(7)=EPS1(I+NDX,J+NDY)
      EPCUB(8)=EPS1(I,J+NDY)
200   continue
c      Print*,'section2'
      RETURN
      END



*******************************************************************
      SUBROUTINE SECTIO3(E0,INTERSECT,EPCUB,EPS0,EPS1,
     &                        I,J,NDX,NDY,NX0,NY0,NX,NY)

** Zhi:
** The freeze-out cell is extended to the whole plane in a
** symmetric way, compare with SECTIO2 to see the difference.
*****************************************************************
** THIS ROUTINE CHECKS WHETHER THE eps=E0 -SURFACE INTERSECTS   *
** WITH A 'CUBE' (I,J). IF THAT IS THE CASE, THE LOGICAL        *
** VARIABLE 'INTERSECT' IS RETURNED WITH A VALUE .TRUE. AND     *
** THE MATRIX EPCUB(8) WITH THE eps VALUES OF THE CUBE CORNERS. *
****************************************************************
      Implicit Double Precision (A-H, O-Z)
      DIMENSION EPS0(NX0:NX, NY0:NY),EPS1(NX0:NX, NY0:NY)
      DIMENSION EPCUB(8)
      LOGICAL INTERSECT
      Integer IS,IL,JS,JL ! IS=smaller(I,I+NDX), used to give correct orientation in EPCUB, similar for others

      IF ((E0-EPS0(I,J))*(EPS1(I+NDX,J+NDY)-E0).LT.0.0) THEN
        IF ((E0-EPS0(I+NDX,J))*(EPS1(I,J+NDY)-E0).LT.0.0) THEN
          IF ((E0-EPS0(I+NDX,J+NDY))*(EPS1(I,J)-E0).LT.0.0) THEN
            IF ((E0-EPS0(I,J+NDY))*(EPS1(I+NDX,J)-E0).LT.0.0) THEN
              INTERSECT=.FALSE.
              GO TO 200
            END IF
          END IF
        END IF
      END IF

100   INTERSECT=.TRUE.

      If (NDX>0) Then
        IS = I
        IL = I + NDX
      Else
        IS = I + NDX
        IL = I
      End If

      If (NDY>0) Then
        JS = J
        JL = J + NDY
      Else
        JS = J + NDY
        JL = J
      End If

      EPCUB(1)=EPS0(IS,JS)
      EPCUB(2)=EPS0(IL,JS)
      EPCUB(3)=EPS0(IL,JL)
      EPCUB(4)=EPS0(IS,JL)
      EPCUB(5)=EPS1(IS,JS)
      EPCUB(6)=EPS1(IL,JS)
      EPCUB(7)=EPS1(IL,JL)
      EPCUB(8)=EPS1(IS,JL)
200   continue
      RETURN
      END

*******************************************************************
      SUBROUTINE SECTIONCornelius(E0,INTERSECT,ECUBE,EPS0,EPS1,
     &                            I,J,NDX,NDY,NX0,NY0,NX,NY)

*****************************************************************
** THIS ROUTINE CHECKS WHETHER THE eps=E0 -SURFACE INTERSECTS   *
** WITH A 'CUBE' (I,J). IF THAT IS THE CASE, THE LOGICAL        *
** VARIABLE 'INTERSECT' IS RETURNED WITH A VALUE .TRUE. AND     *
** THE MATRIX ECUBE(2,2,2) WITH THE eps VALUES OF THE CUBE      *
** CORNERS for Cornelius routine.                               *
*****************************************************************
      Implicit Double Precision (A-H, O-Z)
      DIMENSION EPS0(NX0:NX, NY0:NY),EPS1(NX0:NX, NY0:NY)
      DIMENSION ECUBE(0:1,0:1,0:1)
      LOGICAL INTERSECT

      IF ((E0-EPS0(I,J))*(EPS1(I+NDX,J+NDY)-E0).LT.0.0) THEN
        IF ((E0-EPS0(I+NDX,J))*(EPS1(I,J+NDY)-E0).LT.0.0) THEN
          IF ((E0-EPS0(I+NDX,J+NDY))*(EPS1(I,J)-E0).LT.0.0) THEN
            IF ((E0-EPS0(I,J+NDY))*(EPS1(I+NDX,J)-E0).LT.0.0) THEN
              INTERSECT=.FALSE.
              GO TO 200
            END IF
          END IF
        END IF
      END IF

      INTERSECT=.TRUE.
      ECUBE(0,0,0)=EPS0(I,J)
      ECUBE(0,1,0)=EPS0(I+NDX,J)
      ECUBE(0,1,1)=EPS0(I+NDX,J+NDY)
      ECUBE(0,0,1)=EPS0(I,J+NDY)
      ECUBE(1,0,0)=EPS1(I,J)
      ECUBE(1,1,0)=EPS1(I+NDX,J)
      ECUBE(1,1,1)=EPS1(I+NDX,J+NDY)
      ECUBE(1,0,1)=EPS1(I,J+NDY)
200   continue

      RETURN
      END

C###########################################################################
      Subroutine AFreezeoutPro3 (EDEC,TFREEZ, TFLAG,  IEOS,
     &     EPS0,EPS1,TEM0,TEM1,  V10,V20,V11,V21,  EINS,NINT, IW,
     &     N,T, X0,Y0,DX,DY,DT,NXPhy0,NYPhy0,NXPhy,NYPhy,NX0,NY0,NX,NY)
*    a subroutine to calculate the freezeout surface which is changed from Peter  Azhydro0P2
*     TFLAG=0:     TFLAG=1!freeze out by constant   TFLAG=0 constant energy        (I)
*     EDEC: decoupling energy (TFLAG=0)    TFREEZ: decoupling temperature (TFLAG=1)  (I)
*     T=Time   N,Timestep for the largest Loop.     X0=0.0 Y0=0.0
      Implicit Double Precision (A-H, O-Z)

      DIMENSION EPS0(NX0:NX, NY0:NY),EPS1(NX0:NX, NY0:NY) ! Energy density in previous and current step
      DIMENSION TEM0(NX0:NX, NY0:NY),TEM1(NX0:NX, NY0:NY) ! Temperature density in previous and current step
      DIMENSION V10(NX0:NX, NY0:NY),V20(NX0:NX, NY0:NY)   !velocity in X Y in Previous step
      DIMENSION V11(NX0:NX, NY0:NY),V21(NX0:NX, NY0:NY)   !velocity in X Y in current step

      DIMENSION EPCUB(8),SURFS(12,0:2),VMID(0:2)
      DIMENSION IBIT(12,12)

      Parameter(NXD=2,NYD=2, NTD=5)!
      Parameter(MAXERR=150 )!Max Number erro in Cubcal
      PARAMETER (PI=3.141592653d0, HBARC=.19733d0)

      Logical INTERSECT, WRITING
      INTEGER TFLAG, EINS
      INTEGER NXPHY,NYPHY

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
c      IF (MOD(N+1,NTD).EQ.0) THEN   !** N+1 *******************************************************************

        DTD=NTD*DT
        DXD=NXD*DX
        DYD=NYD*DY

         NYFUL=NYPhy+2
         NXFUL=NXPhy+2

         NYFUL0=NYPhy0-2
         NXFUL0=NXPhy0-2

            NINT=0
            NINT2=0
            DO 1234 J=NYFUL0,NYFUL
              Y=Y0 + J*DY
            DO 1235 I=0,NXFUL
                 IF (IEOS .EQ. 2) THEN
                  TEM1(I,J)=SQRT(HBARC*SQRT(10.0*HBARC*EPS1(I,J))/PI)
                 ELSE IF (IEOS.EQ.0) THEN
                 ELSE IF (IEOS.EQ.99) THEN
                 ELSE
                 END IF
 1235       CONTINUE
 1234       CONTINUE


c       TFLAG=0   ! TFLAG=1!freeze out by constant T =2 constant energy

        IF (TFLAG .EQ. 1) THEN
         FREEZ = TFREEZ
        ELSE
         FREEZ = EDEC
        END IF

            DO 500 J=NYFUL0,NYFUL
            IF (MOD(J,NYD) .NE. 0) GOTO 500
               Y=Y0 + J*DY
            DO 510 I=NXFUL0,NXFUL
            IF (MOD(I,NXD) .NE. 0) GOTO 510
               X=X0+I*DX
             Print *,I,J

                IF (TFLAG .EQ. 1) THEN !freeze out by constant T
                  CALL SECTIO2(FREEZ,INTERSECT,EPCUB,TEM0,TEM1,
     &                                I,J,NXD,NYD,NX0,NY0,NX,NY)
                ELSE
                  CALL SECTIO2(FREEZ,INTERSECT,EPCUB,EPS0,EPS1,
     &                                I,J,NXD,NYD,NX0,NY0,NX,NY)
                ENDIF

                IF (INTERSECT) THEN   !***(INTERSECT)
                          NINT=NINT+1
                          NINT2=NINT2+1
                          WRITING = .TRUE.     !--

                  CALL CUBCAL(FREEZ,EPCUB,SURFS,NSURFS,VMID,IBIT,
     &                        DTD,DXD,DYD,NERR,MAXERR)

                  DAH0=SURFS(1,0)
                  DAH1=SURFS(1,1)
                  DAH2=SURFS(1,2)

                  DO 520 NS=2,NSURFS
                    DAH0=DAH0+SURFS(NS,0)
                    DAH1=DAH1+SURFS(NS,1)
                    DAH2=DAH2+SURFS(NS,2)
520               CONTINUE

                  TM=T-DTD+DT+VMID(0)
                  XM=X+VMID(1)
                  YM=Y+VMID(2)

                  CALL PP3(I,J,NXD,NYD,VMID,V10,V11,
     &              NX0,NY0,NX,NY,DTD,DXD,DYD,V1)
                  CALL PP3(I,J,NXD,NYD,VMID,V20,V21,
     &              NX0,NY0,NX,NY,DTD,DXD,DYD,V2)
c                  CALL P3(I,J,NXD,NYD,VMID,BNR0,BNR1  !Bayon Number
c     &              ,MAX,MAY,DTD,DXD,DYD,BN)
                  BN=0.0

                  IF (TFLAG .EQ. 1) THEN   !freeze out by constant T
                    CALL PP3(I,J,NXD,NYD,VMID,EPS0,EPS1,
     &                NX0,NY0,NX,NY,DTD,DXD,DYD,EDEC)
                    TDEC=FREEZ
                  ELSE
                     IF (IEOS .EQ. 2) THEN
                       TDEC=SQRT(HBARC*SQRT(10.0*HBARC*EDEC)/PI)
                     ELSE IF (IEOS.EQ.10) THEN
                     ELSE IF (IEOS.EQ.99) THEN
                     ELSE
                     END IF
                  END IF

c                 BAMU = EOUT(BN,EDEC,0.0,
                  BAMU=0.0

                  SMU=0.0

                  DA0=DAH0*1.
                  DA1=DAH1*1.
                  DA2=DAH2*1.

                  VZCM=V1
                  VRCM=V2

                  IF (WRITING) THEN   !---
                         WRITE(98,2553) T-DTD+DT, TM, XM,YM,
     &                                               Sqrt(XM**2+YM**2)   !surface
                         IW=IW+1   !count the number in surface written file
                       IF (IEOS.EQ.99) THEN
                       ELSE IF (IEOS.eq.2)THEN !added by even
                         WRITE(99,2555) TM,DA0,DA1,DA2,VZCM,VRCM, !
     &                             EDEC,BN,TDEC,BAMU,SMU
                       ELSE ! IF (IEOS.eq0)THEN
                       END IF

                  ELSE               !---
                    NINT = NINT - 1
                  END IF             !---

                END IF      !***(INTERSECT)
510         CONTINUE
500         CONTINUE

c          END IF    !   IF (MOD(N+1,NTD).EQ.0) THEN   !** N+1 *******************************************************************
C          EINS=EINS+1

 2553             FORMAT(6F15.4)
 2555             FORMAT(12E14.6)
 2565             FORMAT(14E14.6)

      Return
      End

C################################################################
****************************************************************
      Subroutine  NEWTON2(B0,E0,T00IJ,T01IJ,T02IJ,BN0IJ,
     &                             BulkPr,PNEW,NNEW)  !bisection
****************************************************************
*     NEWTON used to call the Subroutine BISEC1, and if this Routine
*     failed to give a solution for the velocity, it called different
*     NEWTON algorithms. BISEC1 was never observed not to lead to a
*     solution the Newton algorithms were taken out of the program
*     (they still exist in the file NEWTON.f)
*     The program will now stop and display an error message in case
*     there was no solution found.
***************************************************************
C ---- --- updated |V| rooting finding with bulk pressure-BulkPr
      Implicit Double Precision (A-H, O-Z)
      DIMENSION PNEW(NNEW)
      Common /Newtonalart/ AI,AJ,AK,AE


      A=T00IJ
      B=SQRT(T01IJ*T01IJ+T02IJ*T02IJ)
      C=BN0IJ

      CALL BISEC2(A,B,C,BulkPr,EE0,BB0,IFOK,PNEW,NNEW)

      IF (IFOK .EQ. 2) THEN
        WRITE(*,*) 'Bisec did not yield a solution! IFOK 2 '
        WRITE(*,*) 'F1F2>0, F1F2<0'
        Write(*,*) 'I,J,K , Ed(I,J,K)', AI,AJ,AK,AE
!        STOP
      END IF

      IF (IFOK .EQ. 0) THEN
        WRITE(*,*) 'PNEW(2)is not enough for Bisec a solution IFork 0'
        Write(*,*) 'I,J,K , Ed(I,J,K)', AI,AJ,AK,AE
!        STOP
      ELSE
        E0 = EE0
        B0 = BB0
      END IF

      RETURN
      END


********************************************************************
      SUBROUTINE BISEC2(A,B,C,BulkPr,ED,BD,IFOK,PNEW,NNEW)
        Implicit Double Precision (A-H, O-Z)
      DIMENSION PNEW(NNEW)

      COMMON /PAREOS/ BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
     &                BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,
     &                PMAT1(0:300,0:250),PMAT2(0:300,0:250)


      IFOK = 1
      ERR=PNEW(1)
      ERR2=ERR*ERR
      ITMAX=INT(PNEW(2)+0.5)
      NIT = 0

      V1 = 0.0
      ED = A-B*V1
      GAMINV = SQRT(1.0-V1*V1)
      BD = C*GAMINV
      FU1 = FUNCS(V1,ED,BD,A,B,bulkPr)
      IF (ABS(FU1).LT. ERR) THEN     !==============================V=0.0  1
        V = 0.0
      ELSE
        V2 = 1.0-1e-15
        ED = A-B*V2
        GAMINV = SQRT(1.0-V2*V2)
        BD = C*GAMINV
        FU2 = FUNCS(V2,ED,BD,A,B,bulkPr)

        IF (ABS(FU2).LT. ERR) THEN   !============================ V=1.0 2
          V = V2
        ELSE
          IF (FU1*FU2 .GT. 0.0) THEN  !=========================== no solutin      $
*     * go and display error message
            IFOK = 2                   !IForK=2
            V = 0.0
            ED = A
            BD = C
          ELSE

            FU3 = 1.0
 199        CONTINUE
            IF ( (ABS(FU3).GT.ERR) .AND. (NIT .LT. ITMAX))THEN   !=====besiction 4

              NIT = NIT+1
              V3 = 0.5*(V1 + V2)
              ED = A-B*V3
              GAMINV = SQRT(1.0-V3*V3)
              BD = C*GAMINV
              FU3 = FUNCS(V3,ED,BD,A,B,bulkPr)  !include EOS
              IF (FU3*FU1 .GT. 0.0) THEN
                V1 = V3
                FU1 = FU3
              ELSE
                V2 = V3
                FU2 = FU3
              END IF
              GOTO 199
            ELSE
              IF (NIT.GT.(ITMAX-1)) THEN
                IFOK =0 ! 1   !IFOK=0? song
              END IF
              V = V3
            END IF
          END IF
        END IF
      END IF

      RETURN
      END

***********************************************************
      Double Precision FUNCTION FUNCS(V,ED,BD,A,B,bulkPr)
        Implicit Double Precision (A-H, O-Z)
***********************************************************
** THIS FUNCTION IS USED BY SUBROUTINE 'NEWTON4'          *
** IN SOLVING THE NEW BARYON NUMBER AND ENERGY DENSITIES  *
***********************************************************

c      COMMON /PAREOS/ BNP01,EPP01,DBNP1,DEPP1,NBNP1,NEPP1,
c     &                BNP02,EPP02,DBNP2,DEPP2,NBNP2,NEPP2,
c     &                PMAT1(0:300,0:250),PMAT2(0:300,0:250)
        Parameter (HbarC=0.19733d0) !for changcing between fm and GeV ! Hbarc=0.19733=GeV*fm

       ee=Ed*HbarC   ! fm^-4 to GeV/fm^3
       PR=PEPS(BD,ee)/HbarC  !GEV/fm^3 to fm^-4
       FUNCS = V*(A+PR+bulkPr) - B


      RETURN
      END

C###################################################################


C####################################################################
      Subroutine ASM2D6(AA,IAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Implicit Double Precision (A-H, O-Z)
       Dimension AA(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension IAA(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension FAA(NX0-3:NX+3, NY0-3:NY+3, NZ0:NZ)
        C0=0.486/1.514
        C1=0.343/1.514
        C2=-0.086/1.514

         do 100 k=NZ0,NZ
         do 100 i=NXPhy0-3,NXPhy+3
         do 100 j=NYPhy0-3,NYPhy+3
        FAA(i,j,k)=AA(i,j,k)
 100   continue
       do 300 i=NXPhy0-2,NXPhy+2
       do 300 j=NYPhy0-2,NYPhy+2
c       If(IAA(i,j,NZ0).eq.0) then
         AA(i,j,NZ0)=C0*FAA(i,j,NZ0) +C1*(FAA(i-1,j,NZ0)+FAA(i+1,j,NZ0)
     &     +FAA(i,j-1,NZ0)+FAA(i,j+1,NZ0)) +C2*(FAA(i-1,j-1,NZ0)
     &     +FAA(i+1,j+1,NZ0)+FAA(i+1,j-1,NZ0)+FAA(i-1,j+1,NZ0))
c       else
        AA(i,j,NZ0)=AA(i,j,NZ0)
c       end if
 300   continue
      Return
      End
C####################################################################
      Subroutine CoefASM2d(CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
      Implicit Double Precision (A-H, O-Z)
      Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)
      Common/dxdy/ ddx, ddy
      Common /R0Bdry/ R0Bdry
      Double Precision R0Bdry, AepsBdry

      AepsBdry = 0.5
      !Print *, "R0Bdry=", R0Bdry

        C0=0.486d0/1.514d0
        C1=0.343d0/1.514d0
        C2=-0.086d0/1.514d0

       do 300 i=NXPhy0-2,NXPhy+2
       do 300 j=NYPhy0-2,NYPhy+2
         xx=ddx*I
         yy=ddy*J
         rr=sqrt(xx**2+yy**2)
         ff=1.0/(Dexp((rr-R0Bdry)/AepsBdry)+1.0)

        C0a=(1.0-C0)*ff+C0
        C1a=C1*(1-ff)
        C2a=C2*(1-ff)

        CofAA(0,I,J,NZ0)=C0a
        CofAA(1,I,J,NZ0)=C1a
        CofAA(2,I,J,NZ0)=C2a

300   continue
      Return
      End



C####################################################################
      Subroutine ASM2D6bsm(AA,CofAA, NXPhy0,NYPhy0, NXPhy,NYPhy,
     &                        NX0,NY0,NZ0, NX,NY,NZ)
        Implicit Double Precision (A-H, O-Z)
       Dimension AA(NX0:NX, NY0:NY, NZ0:NZ)
       Dimension FAA(NX0-3:NX+3, NY0-3:NY+3, NZ0:NZ)
       Dimension CofAA(0:2,NX0:NX, NY0:NY, NZ0:NZ)

         do 100 k=NZ0,NZ
         do 100 i=NXPhy0-3,NXPhy+3
         do 100 j=NYPhy0-3,NYPhy+3
        FAA(i,j,k)=AA(i,j,k)
 100   continue

       do 300 i=NXPhy0-2,NXPhy+2
       do 300 j=NYPhy0-2,NYPhy+2

            C0=CofAA(0,I,J,NZ0)
            C1=CofAA(1,I,J,NZ0)
            C2=CofAA(2,I,J,NZ0)
          AA(i,j,NZ0)=C0*FAA(i,j,NZ0) +C1*(FAA(i-1,j,NZ0)+FAA(i+1,j,NZ0)
     &    +FAA(i,j-1,NZ0)+FAA(i,j+1,NZ0)) +C2*(FAA(i-1,j-1,NZ0)
     &    +FAA(i+1,j+1,NZ0)+FAA(i+1,j-1,NZ0)+FAA(i-1,j+1,NZ0))

300   continue
      Return
      End


C###############################################################################################
C###############################################################################################
      Subroutine VSBdary3(Vx,Vy, NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)
        Implicit Double Precision (A-H, O-Z)
        Dimension Vx(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension Vy(NX0:NX, NY0:NY, NZ0:NZ)
        LOGICAL V1FOUND,V2FOUND,V3FOUND,V4FOUND
C-----------Boundary Treatment-for Velocity ------------------

      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        Vx(I,NYPhy+kk,K)=2.0*Vx(I,NYPhy,K)-Vx(I,NYPhy-kk,K)
        Vx(I,NYPhy0-kk,K)=2.0*Vx(I,NYPhy0,K)-Vx(I,NYPhy0+kk,K)

        Vy(I,NYPhy+kk,K)=2.0*Vy(I,NYPhy,K)-Vy(I,NYPhy-kk,K)
        Vy(I,NYPhy0-kk,K)=2.0*Vy(I,NYPhy0,K)-Vy(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        Vx(NXPhy+kk,J,K)=2.0*Vx(NXPhy,J,K)-Vx(NXPhy-kk,J,K)
        Vx(NXPhy0-kk,J,K)=2.0*Vx(NXPhy0,J,K)-Vx(NXPhy0+kk,J,K)

        Vy(NXPhy+kk,J,K)=2.0*Vy(NXPhy,J,K)-Vy(NXPhy-kk,J,K)
        Vy(NXPhy0-kk,J,K)=2.0*Vy(NXPhy0,J,K)-Vy(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        Vx(I,NYPhy+kk,K)=(2.0*Vx(I,NYPhy,K)-Vx(I,NYPhy-kk,K))*0.5
     &                                      +0.5*Vx(I,NYPhy+kk,K)
        Vx(I,NYPhy0-kk,K)=(2.0*Vx(I,NYPhy0,K)-Vx(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*Vx(I,NYPhy0-kk,K)

        Vy(I,NYPhy+kk,K)=(2.0*Vy(I,NYPhy,K)-Vy(I,NYPhy-kk,K))*0.5
     &                                    +0.5*Vy(I,NYPhy+kk,K)
        Vy(I,NYPhy0-kk,K)=(2.0*Vy(I,NYPhy0,K)-Vy(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*Vy(I,NYPhy0-kk,K)
        end do
1820  CONTINUE



C          goto 999
C-------------Some other treatment of Velocity-------------------------
          DO 3310 K=NZ0,NZ
          DO 3310 I=NXPhy0-3,NXPhy+3
            V1FOUND=.TRUE.
            V2FOUND=.TRUE.
            V3FOUND=.TRUE.
            V4FOUND=.TRUE.

            DO 3311 J=NYPhy-1,0,-1
             AVy=Vy(I,J,K)
      IF (V1FOUND.AND.ABS(AVy).GT.AAC0) THEN
        V1FOUND=.FALSE.
          Vy(I,J+1,K)=2.0*Vy(I,J,K)-Vy(I,J-1,K)
        Vy(I,J+2,K)=2.0*Vy(I,J,K)-Vy(I,J-2,K)
        Vy(I,J+3,K)=2.0*Vy(I,J,K)-Vy(I,J-3,K)
        Vy(I,J+4,K)=2.0*Vy(I,J,K)-Vy(I,J-4,K)
      ENDIF
 3311     CONTINUE

            DO 3312 J=NYPhy0+1,0,+1
             AVy=Vy(I,J,K)
      IF (V2FOUND.AND.ABS(AVy).GT.AAC0) THEN
        V2FOUND=.FALSE.
          Vy(I,J-1,K)=2.0*Vy(I,J,K)-Vy(I,J+1,K)
        Vy(I,J-2,K)=2.0*Vy(I,J,K)-Vy(I,J+2,K)
        Vy(I,J-3,K)=2.0*Vy(I,J,K)-Vy(I,J+3,K)
        Vy(I,J-4,K)=2.0*Vy(I,J,K)-Vy(I,J+4,K)
      ENDIF
 3312     CONTINUE

 3310     CONTINUE

        DO 3320 K=NZ0,NZ
        DO 3320 J=NYPhy0-3,NYPhy+3

          V1FOUND=.TRUE.
          V2FOUND=.TRUE.
          V3FOUND=.TRUE.
          V4FOUND=.TRUE.

          DO 3321 I=NXPhy-1,0,-1
             AVx=Vx(I,J,K)
          IF (V3FOUND.AND.ABS(AVx).GT.AAC0) THEN
      V3FOUND=.FALSE.
C            Print *,I,J, AAC0,AVx
      Vx(I+1,J,K)=2.0*Vx(I,J,K)-Vx(I-1,J,K)
      Vx(I+2,J,K)=2.0*Vx(I,J,K)-Vx(I-2,J,K)
      Vx(I+3,J,K)=2.0*Vx(I,J,K)-Vx(I-3,J,K)
      Vx(I+4,J,K)=2.0*Vx(I,J,K)-Vx(I-4,J,K)
          END IF
 3321  CONTINUE

          DO 3322 I=NXPhy0+1,0,+1
             AVx=Vx(I,J,K)
          IF (V4FOUND.AND.ABS(AVx).GT.AAC0) THEN
      V4FOUND=.FALSE.
      Vx(I-1,J,K)=2.0*Vx(I,J,K)-Vx(I+1,J,K)
      Vx(I-2,J,K)=2.0*Vx(I,J,K)-Vx(I+2,J,K)
      Vx(I-3,J,K)=2.0*Vx(I,J,K)-Vx(I+3,J,K)
      Vx(I-4,J,K)=2.0*Vx(I,J,K)-Vx(I+4,J,K)
          END IF
 3322  CONTINUE

 3320  CONTINUE

 999   Continue
      Return
      End





C###############################################################################################
C###############################################################################################
      Subroutine VSBdary3Single(VV,NX0,NY0,NZ0,NX,NY,NZ,
     &  NXPhy0,NYPhy0, NXPhy,NYPhy,AAC0)
        Implicit Double Precision (A-H, O-Z)
        Dimension VV(NX0:NX, NY0:NY, NZ0:NZ)
        LOGICAL V1FOUND,V2FOUND,V3FOUND,V4FOUND
C-----------Boundary Treatment-for Velocity ------------------

      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        VV(I,NYPhy+kk,K)=2.0*VV(I,NYPhy,K)-VV(I,NYPhy-kk,K)
        VV(I,NYPhy0-kk,K)=2.0*VV(I,NYPhy0,K)-VV(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        VV(NXPhy+kk,J,K)=2.0*VV(NXPhy,J,K)-VV(NXPhy-kk,J,K)
        VV(NXPhy0-kk,J,K)=2.0*VV(NXPhy0,J,K)-VV(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        VV(I,NYPhy+kk,K)=(2.0*VV(I,NYPhy,K)-VV(I,NYPhy-kk,K))*0.5
     &                                      +0.5*VV(I,NYPhy+kk,K)
        VV(I,NYPhy0-kk,K)=(2.0*VV(I,NYPhy0,K)-VV(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*VV(I,NYPhy0-kk,K)
        end do
1820  CONTINUE



C          goto 999
C-------------Some other treatment of Velocity-------------------------
!          DO 3310 K=NZ0,NZ
!          DO 3310 I=NXPhy0-3,NXPhy+3
!            V1FOUND=.TRUE.
!            V2FOUND=.TRUE.
!            V3FOUND=.TRUE.
!            V4FOUND=.TRUE.
!
!            DO 3311 J=NYPhy-1,0,-1
!             AVy=Vy(I,J,K)
!      IF (V1FOUND.AND.ABS(AVy).GT.AAC0) THEN
!        V1FOUND=.FALSE.
!          Vy(I,J+1,K)=2.0*Vy(I,J,K)-Vy(I,J-1,K)
!        Vy(I,J+2,K)=2.0*Vy(I,J,K)-Vy(I,J-2,K)
!        Vy(I,J+3,K)=2.0*Vy(I,J,K)-Vy(I,J-3,K)
!        Vy(I,J+4,K)=2.0*Vy(I,J,K)-Vy(I,J-4,K)
!      ENDIF
! 3311     CONTINUE
!
!            DO 3312 J=NYPhy0+1,0,+1
!             AVy=Vy(I,J,K)
!      IF (V2FOUND.AND.ABS(AVy).GT.AAC0) THEN
!        V2FOUND=.FALSE.
!          Vy(I,J-1,K)=2.0*Vy(I,J,K)-Vy(I,J+1,K)
!        Vy(I,J-2,K)=2.0*Vy(I,J,K)-Vy(I,J+2,K)
!        Vy(I,J-3,K)=2.0*Vy(I,J,K)-Vy(I,J+3,K)
!        Vy(I,J-4,K)=2.0*Vy(I,J,K)-Vy(I,J+4,K)
!      ENDIF
! 3312     CONTINUE
!
! 3310     CONTINUE
!
!        DO 3320 K=NZ0,NZ
!        DO 3320 J=NYPhy0-3,NYPhy+3
!
!          V1FOUND=.TRUE.
!          V2FOUND=.TRUE.
!          V3FOUND=.TRUE.
!          V4FOUND=.TRUE.
!
!          DO 3321 I=NXPhy-1,0,-1
!             AVV=VV(I,J,K)
!          IF (V3FOUND.AND.ABS(AVV).GT.AAC0) THEN
!      V3FOUND=.FALSE.
!C            Print *,I,J, AAC0,AVV
!      VV(I+1,J,K)=2.0*VV(I,J,K)-VV(I-1,J,K)
!      VV(I+2,J,K)=2.0*VV(I,J,K)-VV(I-2,J,K)
!      VV(I+3,J,K)=2.0*VV(I,J,K)-VV(I-3,J,K)
!      VV(I+4,J,K)=2.0*VV(I,J,K)-VV(I-4,J,K)
!          END IF
! 3321  CONTINUE
!
!          DO 3322 I=NXPhy0+1,0,+1
!             AVV=VV(I,J,K)
!          IF (V4FOUND.AND.ABS(AVV).GT.AAC0) THEN
!      V4FOUND=.FALSE.
!      VV(I-1,J,K)=2.0*VV(I,J,K)-VV(I+1,J,K)
!      VV(I-2,J,K)=2.0*VV(I,J,K)-VV(I+2,J,K)
!      VV(I-3,J,K)=2.0*VV(I,J,K)-VV(I+3,J,K)
!      VV(I-4,J,K)=2.0*VV(I,J,K)-VV(I+4,J,K)
!          END IF
! 3322  CONTINUE
!
! 3320  CONTINUE

 999   Continue
      Return
      End




C########################################################################
      Subroutine Sembdary3(PL,NX0,NY0,NZ0,NX,NY,NZ,
     &                         NXPhy0,NYPhy0,NXPhy,NYPhy)
        Implicit Double Precision (A-H, O-Z)
        Dimension PL(NX0:NX, NY0:NY, NZ0:NZ)
C-----------Boundary Treatment-------------------
      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        PL(I,NYPhy+kk,K)=2.0*PL(I,NYPhy,K)-PL(I,NYPhy-kk,K)
        PL(I,NYPhy0-kk,K)=2.0*PL(I,NYPhy0,K)-PL(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        PL(NXPhy+kk,J,K)=2.0*PL(NXPhy,J,K)-PL(NXPhy-kk,J,K)
        PL(NXPhy0-kk,J,K)=2.0*PL(NXPhy0,J,K)-PL(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        PL(I,NYPhy+kk,K)=(2.0*PL(I,NYPhy,K)-PL(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL(I,NYPhy+kk,K)
        PL(I,NYPhy0-kk,K)=(2.0*PL(I,NYPhy0,K)-PL(I,NYPhy0+kk,K))*0.5
     &                                    +0.5*PL(I,NYPhy0-kk,K)
        end do
1820  CONTINUE

      Return
      End


C########################################################################
      Subroutine TriSembdary3(PL1, PL2, PL3,
     &        NX0,NY0,NZ0, NX,NY,NZ, NXPhy0,NYPhy0, NXPhy,NYPhy)
        Implicit Double Precision (A-H, O-Z)
        Dimension PL1(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension PL2(NX0:NX, NY0:NY, NZ0:NZ)
        Dimension PL3(NX0:NX, NY0:NY, NZ0:NZ)
C-----------Boundary Treatment-------------------

      DO 1800 K=NZ0,NZ
      DO 1800 I=NXPhy0,NXPhy
        do kk=1,3
        PL1(I,NYPhy+kk,K)=2.0*PL1(I,NYPhy,K)-PL1(I,NYPhy-kk,K)
        PL1(I,NYPhy0-kk,K)=2.0*PL1(I,NYPhy0,K)-PL1(I,NYPhy0+kk,K)

        PL2(I,NYPhy+kk,K)=2.0*PL2(I,NYPhy,K)-PL2(I,NYPhy-kk,K)
        PL2(I,NYPhy0-kk,K)=2.0*PL2(I,NYPhy0,K)-PL2(I,NYPhy0+kk,K)

        PL3(I,NYPhy+kk,K)=2.0*PL3(I,NYPhy,K)-PL3(I,NYPhy-kk,K)
        PL3(I,NYPhy0-kk,K)=2.0*PL3(I,NYPhy0,K)-PL3(I,NYPhy0+kk,K)
        end do
1800  CONTINUE

      DO 1810 K=NZ0,NZ
      DO 1810 J=NYPhy0-3,NYPhy+3
        do kk=1,3
        PL1(NXPhy+kk,J,K)=2.0*PL1(NXPhy,J,K)-PL1(NXPhy-kk,J,K)
        PL1(NXPhy0-kk,J,K)=2.0*PL1(NXPhy0,J,K)-PL1(NXPhy0+kk,J,K)

        PL2(NXPhy+kk,J,K)=2.0*PL2(NXPhy,J,K)-PL2(NXPhy-kk,J,K)
        PL2(NXPhy0-kk,J,K)=2.0*PL2(NXPhy0,J,K)-PL2(NXPhy0+kk,J,K)

        PL3(NXPhy+kk,J,K)=2.0*PL3(NXPhy,J,K)-PL3(NXPhy-kk,J,K)
        PL3(NXPhy0-kk,J,K)=2.0*PL3(NXPhy0,J,K)-PL3(NXPhy0+kk,J,K)
        end do
1810  CONTINUE

      DO 1820 K=NZ0,NZ
      DO 1820 I=NXPhy0-3,NXPhy-3
        do kk=1,3
        PL1(I,NYPhy+kk,K)=(2.0*PL1(I,NYPhy,K)-PL1(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL1(I,NYPhy+kk,K)
        PL1(I,NYPhy0-kk,K)=(2.0*PL1(I,NYPhy0,K)-PL1(I,NYPhy0+kk,K))*0.5
     &                                      +0.5*PL1(I,NYPhy0-kk,K)

        PL2(I,NYPhy+kk,K)=(2.0*PL2(I,NYPhy,K)-PL2(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL2(I,NYPhy+kk,K)
        PL2(I,NYPhy0-kk,K)=(2.0*PL2(I,NYPhy0,K)-PL2(I,NYPhy0+kk,K))*0.5
     &                                      +0.5*PL2(I,NYPhy0-kk,K)

        PL3(I,NYPhy+kk,K)=(2.0*PL3(I,NYPhy,K)-PL3(I,NYPhy-kk,K))*0.5
     &                                    +0.5*PL3(I,NYPhy+kk,K)
        PL3(I,NYPhy0-kk,K)=(2.0*PL3(I,NYPhy0,K)-PL3(I,NYPhy0+kk,K))*0.5
     &                                      +0.5*PL3(I,NYPhy0-kk,K)
        end do
1820  CONTINUE

      Return
      End
