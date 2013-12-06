!  cornelius2 version 1.3: Copyright 2012, Pasi Huovinen
!
! This subroutine is aimed to be used as a part of the fluid dynamical
! models of the heavy-ion physics community. Permission to use it for
! any purpose except for any commercial purpose is granted, provided
! that any publication cites the paper describing the algorithm: 
! P. Huovinen and H. Petersen, arXiv:1206.3371
!
! Permission to distribute this subroutine is granted, provided that no
! fee is charged, and that this copyright and permission notice appear
! in all the copies. Permission to modify this subroutine is granted
! provided that the modified subroutine is made publicly available
! latest when any results obtained using the modified subroutine are
! published, the modified subroutine is distributed under terms similar
! to this notice, and the modified code carries both the original
! copyright notice and notices stating that you modified it, and a
! relevant date/year.
!
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of FITNESS FOR
! A PARTICULAR PURPOSE.


   SUBROUTINE Cornelius2(E0,Cube,dSigma,Nsurf,Vmid,dt,dx,dy,Nambi,Ndisc)
!
! version 1.3 zeta
!
! This routine search for a 2-dimensional isosurface of constant X in
! a volume-element (cube) of 3-dimensional space when the values of X
! are known at the vertices (=corners) of the cube and X is interpolated 
! linearly between the vertices. I.e. the usual problem of finding the
! freeze-out surface.
! The routine devides this surface into triangles, evaluates their
! sizes and normal vectors, and provides an approximation (dSigma) of
! the normal vector of the surface as a sum of the normal vectors of
! the triangles. The length of the normal vector gives the size (area)
! of the surface.
!
!    Variables: E0:      the value of X on the surface (INPUT)
!               Cube:    3D cube (2x2x2) to store the values of X (INPUT)
!               dSigma:  3-vector table to store the normal vector(s) of
!                        the surface element(s), |dSigma| is the area
!                        of the surface element (OUTPUT)
!               Nsurf:   Number of separate surface elements within the cube
!                        (OUTPUT)
!               Vmid:    coordinates of the approximate centroid(s) of the
!                        element(s) if the origin is at (0,0,0) corner of
!                        the cube (OUTPUT)
!               dt, dx, dy:   the lengths of the edges of the cube (INPUT)
!               Nambi:   number of ambiguous faces on surfaces (INOUT)
!                        (do not set this to zero between successive calls)
!               Ndisc:   number of disconnected surface-elements so far (INOUT)
!                        (do not set this to zero between successive calls)
!
!        -- P. Huovinen, Jyvaskyla-Frankfurt, March 2005-March 2011 --
!
! Changes in version 1.2:
! - improves the handling of a case where a corner or the center of the face
!   is exactly at the freeze-out temperature and the face is ambiguous. In
!   previous version the behaviour dependend on which corner was at the FO
!   temperature, now the treatment is consistent depending only on the value
!   in the center of the face
!        -- PH, Heraklion, Sept 2011 --
!
! Changes in version 1.3:
! - added the license statement and a check that the cube really contains
!   a surface element
!        -- PH, Frankfurt, July 2012 --
!
! The ordering of values in Cube(t,i,j):
! first index time, second x, third y
!  Cube(0,0,0) <-> t=0,x=0,y=0
!  Cube(0,0,1) <-> t=0,x=0,y=1
!  Cube(0,1,0) <-> t=0,x=1,y=0
!  Cube(0,1,1) <-> t=0,x=1,y=1
!  Cube(1,0,0) <-> t=1,x=0,y=0
!  Cube(1,0,1) <-> t=1,x=0,y=1
!  Cube(1,1,0) <-> t=1,x=1,y=0
!  Cube(1,1,1) <-> t=1,x=1,y=1
!
     IMPLICIT NONE

     REAL(KIND(0D0)),INTENT(IN) :: E0
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1),INTENT(IN) :: Cube
     REAL(KIND(0D0)),DIMENSION(0:2,4),INTENT(OUT)      :: dSigma
     INTEGER,INTENT(OUT)        :: Nsurf
     REAL(KIND(0D0)),DIMENSION(0:2,4),INTENT(OUT)      :: Vmid
     REAL(KIND(0D0)),INTENT(IN) :: dt, dx, dy
     INTEGER,INTENT(INOUT)      :: Nambi, Ndisc

     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Edge ! Table for ends of edges  
                                                 ! i.e. corners of the polygons
     INTEGER :: Nedge, Ncorners
     LOGICAL :: Ambiguous
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut   ! Outside direction
     LOGICAL :: Pathological
     INTEGER,DIMENSION(5) :: EdgeSet ! which edges belong to different surfaces
     INTEGER :: j

     Ncorners = COUNT(Cube .ge. E0)
     dSigma = 0D0
! Check that the cube really contains a surface element
     IF ((Ncorners .gt. 0).and.(Ncorners .lt. 8)) THEN
        Nedge = 0
        Ut    = 0D0
        Ambiguous = .false.
        Pathological = .false.

! Find surface corners and edges:
        CALL Edges(E0,Cube,Edge,Ut,dt,dx,dy,Nedge,Ambiguous,Pathological)
        IF (Pathological) CALL DeadEnd(Cube,E0)

! Check disconnectedness:
        IF ((.not.(Ambiguous)).and.(Nedge .eq. 6))                           &
             Ambiguous = (((Cube(0,0,0)-E0)*(Cube(1,1,1)-E0).gt.0D0).and.    &
                          ((Cube(1,0,0)-E0)*(Cube(0,1,1)-E0).gt.0D0).and.    &
                          ((Cube(1,1,0)-E0)*(Cube(0,0,1)-E0).gt.0D0).and.    &
                          ((Cube(0,1,0)-E0)*(Cube(1,0,1)-E0).gt.0D0))
        IF (Ambiguous) THEN
           Nambi = Nambi + 1
           CALL Disconnected2D(Nedge,Edge,Ut,Nsurf,EdgeSet)
        ELSE
           Nsurf = 1
           EdgeSet(1) = 1
           EdgeSet(2) = Nedge+1
        END IF
        IF (Nsurf .gt. 1) Ndisc = Ndisc+1

! Evaluate the center and normal for each surface:
        DO j = 1,Nsurf
           Nedge = EdgeSet(j+1)-EdgeSet(j)
           CALL NormalVector(j,EdgeSet(j),EdgeSet(j+1)-1,Nedge,Edge,Ut,      &
                             dSigma,Vmid)
        END DO
     ELSE
        Nsurf = 0
        Vmid  = 0D0
     END IF

   END SUBROUTINE Cornelius2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE NormalVector(j,First,Last,Nedge,Edge,Ut,dSigma,Vmids)

     IMPLICIT NONE
     INTEGER :: i, j, First, Last, Nedge
     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Surfs, Ut
     REAL(KIND(0D0)),DIMENSION(0:2,4)    :: dSigma, Vmids
     REAL(KIND(0D0)),DIMENSION(0:2)      :: Vmid, Wmid, V, A, B, SuL, Vout
     REAL(KIND(0D0)) :: Area, AreaI, Vsum

! Calculate the mean vector of intersection points

     V = 0D0
     DO i = First,Last
        V = V + Edge(:,1,i) + Edge(:,2,i)
     END DO
     Vmid = V/(2*Nedge)

! Calculate the centroid, i.e., the center of gravity of the surface element
! (if Nedge=3, the surface element is a triangle and centroid is the mean of
!  corner coordinates)

     IF (Nedge .eq. 3) THEN
        Wmid = Vmid
     ELSE
        Area = 0D0
        V    = 0D0
        DO i = First,Last
           A = Edge(:,1,i) - Vmid
           B = Edge(:,2,i) - Vmid
           AreaI = Sqrt((A(2)*B(1) - A(1)*B(2))**2       &
                       +(A(0)*B(2) - A(2)*B(0))**2       &
                       +(A(1)*B(0) - A(0)*B(1))**2)/2
           Area = Area + AreaI
           V = V + (Edge(:,1,i) + Edge(:,2,i) + Vmid)*AreaI/3
        END DO
        Wmid = V/Area
     END IF

! Start calculating the surface vector

     DO i = First,Last
        A = Edge(:,1,i) - Wmid
        B = Edge(:,2,i) - Wmid

        SuL(0) = 0.5*(A(2)*B(1) - A(1)*B(2))     ! Covariant components
        SuL(1) = 0.5*(A(0)*B(2) - A(2)*B(0))     ! of the surface
        SuL(2) = 0.5*(A(1)*B(0) - A(0)*B(1))     ! normal vector

! Choose the direction towards lower energy

        Vout = Ut(:,i) - Wmid
        Vsum = Vout(0)*SuL(0) + Vout(1)*SuL(1) + Vout(2)*SuL(2)
        Surfs(:,i) = SIGN(1D0,Vsum)*SuL
     END DO

     dSigma(:,j) = SUM(Surfs(:,First:Last),DIM=2)  ! Surface normal
     Vmids(:,j) = Wmid

   END SUBROUTINE NormalVector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE DeadEnd(Cube,E0)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1) :: Cube

     WRITE(*,*) 'Error in CubeCut. Impossible surface.'
     WRITE(*,*) 'Freeze-out value:',E0
     WRITE(*,*) 'Values at Cube corners:'
     WRITE(*,*) 'E(1,0,1) =',Cube(1,0,1),'E(1,1,1) =',Cube(1,1,1)
     WRITE(*,*) 'E(1,0,0) =',Cube(1,0,0),'E(1,1,0) =',Cube(1,1,0)
     WRITE(*,*)
     WRITE(*,*) 'E(0,0,1) =',Cube(0,0,1),'E(0,1,1) =',Cube(0,1,1)
     WRITE(*,*) 'E(0,0,0) =',Cube(0,0,0),'E(0,1,0) =',Cube(0,1,0)
     STOP

   END SUBROUTINE DeadEnd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Subroutines below this line are identical to subroutines in     !
!!!! 3+1D cornelius, except for the size of the Edge array           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE Edges(E0,Cube,Edge,Ut,dt,dx,dy,Nedge,Outo,Patho)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dx,dy,dt
     REAL(KIND(0D0)),DIMENSION(0:1,0:1,0:1) :: Cube
     REAL(KIND(0D0)),DIMENSION(0:1,0:1)  :: Square
     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut
     REAL(KIND(0D0)),DIMENSION(2,4)      :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2)      :: Out
     INTEGER         :: Ncuts,Nedge, i
     LOGICAL         :: Outo, Patho

     DO i = 0,1  ! t = 0 or 1
        Square = Cube(i,:,:)
        CALL FindEdge(E0,Square,Cut,Out,dx,dy,Ncuts,Patho)
        IF (Patho) RETURN
        CALL StoreEdge(Ncuts,Cut,Out,i*dt,0,(/1,2/),Nedge,Edge,Ut,Outo)
     END DO

     DO i = 0,1  ! x = 0 or 1
        Square = Cube(:,i,:)
        CALL FindEdge(E0,Square,Cut,Out,dt,dy,Ncuts,Patho)
        IF (Patho) RETURN
        CALL StoreEdge(Ncuts,Cut,Out,i*dx,1,(/0,2/),Nedge,Edge,Ut,Outo)
     END DO

     DO i = 0,1  ! y = 0 or 1
        Square = Cube(:,:,i)
        CALL FindEdge(E0,Square,Cut,Out,dt,dx,Ncuts,Patho)
        IF (Patho) RETURN
        CALL StoreEdge(Ncuts,Cut,Out,i*dy,2,(/0,1/),Nedge,Edge,Ut,Outo)
     END DO

   END SUBROUTINE Edges


   SUBROUTINE FindEdge(E0,Square,Cut,Out,dx,dy,Ncuts,Patho)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dx,dy
     REAL(KIND(0D0)),DIMENSION(0:1,0:1) :: Square
     REAL(KIND(0D0)),DIMENSION(2,4)     :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2)     :: Out
     INTEGER :: Ncuts
     LOGICAL :: Patho

     CALL EndsOfEdge(E0,Square,Cut,dx,dy,Ncuts)

     IF (Ncuts .gt. 0) CALL FindOutside(Ncuts,E0,Square,Cut,Out,dx,dy)

     IF ((Ncuts .eq. 3).or.(Ncuts .eq. 1)) THEN
        WRITE(*,*) 'Error in FindEdge,',Ncuts,' cuts.'
        WRITE(*,*) 'Noncontinuous surface. E0 =',E0
        WRITE(*,*) 'Eps(0,0) =',Square(0,0),' Eps(1,0) =',Square(1,0)
        WRITE(*,*) 'Eps(0,1) =',Square(0,1),' Eps(1,1) =',Square(1,1)
        WRITE(*,*)
        Patho = .true.
     END IF

   END SUBROUTINE FindEdge


   SUBROUTINE EndsOfEdge(E0,Square,Cut,dx,dy,Ncuts)     

     IMPLICIT NONE
     REAL(KIND(0D0)) :: E0,dx,dy
     REAL(KIND(0D0)),DIMENSION(0:1,0:1) :: Square
     REAL(KIND(0D0)),DIMENSION(2,4)     :: Cut
     INTEGER :: Ncuts

     Ncuts = 0
     IF (((Square(0,0)-E0)*(Square(1,0)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts+1
        Cut(1,Ncuts) = (Square(0,0)-E0)/(Square(0,0)-Square(1,0))*dx
        Cut(2,Ncuts) = 0D0
     ELSE
        IF ((Square(0,0).eq.E0) .or. (Square(1,0).eq.E0))            &
             CALL EndsAtCorner(Square(0,0),Square(1,0),E0,Ncuts,     &
                               Cut(1,Ncuts+1),Cut(2,Ncuts+1),dx,0D0)
     END IF

     IF (((Square(0,0)-E0)*(Square(0,1)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts + 1
        Cut(1,Ncuts) = 0D0
        Cut(2,Ncuts) = (Square(0,0)-E0)/(Square(0,0)-Square(0,1))*dy
     ELSE
        IF ((Square(0,0).eq.E0) .or. (Square(0,1).eq.E0))            &
             CALL EndsAtCorner(Square(0,0),Square(0,1),E0,Ncuts,     &
                               Cut(2,Ncuts+1),Cut(1,Ncuts+1),dy,0D0)
     END IF

     IF (((Square(1,0)-E0)*(Square(1,1)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts+1
        Cut(1,Ncuts) = dx
        Cut(2,Ncuts) = (Square(1,0)-E0)/(Square(1,0)-Square(1,1))*dy
     ELSE
        IF ((Square(1,0).eq.E0) .or. (Square(1,1).eq.E0))            &
             CALL EndsAtCorner(Square(1,0),Square(1,1),E0,Ncuts,     &
                               Cut(2,Ncuts+1),Cut(1,Ncuts+1),dy,dx)
     END IF

     IF (((Square(0,1)-E0)*(Square(1,1)-E0)) .lt. 0D0) THEN
        Ncuts = Ncuts+1
        Cut(1,Ncuts) = (Square(0,1)-E0)/(Square(0,1)-Square(1,1))*dx
        Cut(2,Ncuts) = dy
      ELSE
        IF ((Square(0,1).eq.E0) .or. (Square(1,1).eq.E0))            &
             CALL EndsAtCorner(Square(0,1),Square(1,1),E0,Ncuts,     &
                               Cut(1,Ncuts+1),Cut(2,Ncuts+1),dx,dy)
    END IF

   END SUBROUTINE EndsOfEdge


   SUBROUTINE EndsAtCorner(A,B,E0,Ncuts,C1,C2,d1,d2)

     IMPLICIT NONE
     REAL(KIND(0D0)) :: A,B,E0,C1,C2,d1,d2
     INTEGER :: Ncuts

     IF ((A .eq. E0).and.(B .lt. E0)) THEN
        Ncuts = Ncuts+1
        C1 = 1D-9*d1
        C2 = d2
     END IF
     IF ((A .lt. E0).and.(B .eq. E0)) THEN
        Ncuts = Ncuts+1
        C1 = (1D0-1D-9)*d1
        C2 = d2
     END IF

   END SUBROUTINE EndsAtCorner


!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE FindOutside(Ncuts,E0,Square,Cut,Out,dx,dy)

! Finds a point outside the freeze-out surface and sorts ambiguous surfaces
! with two edges on one face of the cube. The rule is to interpolate the value
! at the center of the face, see if it is below or above freeze-out criterion
! and set the surface accordingly.

     IMPLICIT NONE
     INTEGER :: Ncuts
     REAL(KIND(0D0)) :: E0,dx,dy
     REAL(KIND(0D0)),DIMENSION(0:1,0:1) :: Square
     REAL(KIND(0D0)),DIMENSION(2,4)     :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2)     :: Out

     INTEGER :: i, j, Nout
     REAL(KIND(0D0)) :: Eave

     IF (Ncuts .eq. 4) THEN  ! Ambiguous surface check, interpolate the center
        Eave = 2.5D-1*SUM(Square)
        IF (    ((Square(0,0).lt.E0).and.(Eave.lt.E0)) &
            .or.((Square(0,0).ge.E0).and.(Eave.ge.E0))) THEN
           Out(:,1) = Cut(:,2)
           Cut(:,2) = Cut(:,3)
           Cut(:,3) = Out(:,1)
        END IF
        IF ((Eave-E0) .lt. 0D0) THEN ! Outward direction at the center
           Out(1,:) = 5D-1*dx
           Out(2,:) = 5D-1*dy
        ELSE
           IF ((Square(0,0)-E0) .lt. 0D0) THEN
              Out(:,1) = 0D0
              Out(1,2) = dx
              Out(2,2) = dy
           ELSE
              Out(1,1) = dx
              Out(2,1) = 0D0
              Out(1,2) = 0D0
              Out(2,2) = dy
           END IF
        END IF
     ELSE         ! Normal case, only one edge cutting the face of the cube
        Out = 0D0 ! Find the direction outwards (to lower value)
        Nout = 0
        DO i = 0,1
           DO j = 0,1
              IF (Square(i,j) .lt. E0) THEN
                 Out(1,1) = Out(1,1) + i*dx
                 Out(2,1) = Out(2,1) + j*dy
                 Nout = Nout + 1
              END IF
           END DO
        END DO
        IF (Nout .gt. 0) Out = Out/Nout
     END IF

   END SUBROUTINE FindOutside


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE StoreEdge(Ncuts,Cut,Out,dKnown,Trivial,nonTrivial,       &
                        Nedge,Edge,Ut,Outo)
     IMPLICIT NONE

     LOGICAL              :: Outo
     INTEGER              :: Ncuts, Nedge, Trivial
     INTEGER,DIMENSION(2) :: nonTrivial
     REAL(KIND(0D0)),DIMENSION(2,4) :: Cut
     REAL(KIND(0D0)),DIMENSION(2,2) :: Out
     REAL(KIND(0D0))                :: dKnown
     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut

     IF ((Ncuts .eq. 2).or.(Ncuts .eq. 4)) THEN
        Nedge = Nedge + 1
        Edge(Trivial,:,Nedge) = dKnown
        Edge(nonTrivial,:,Nedge) = Cut(:,1:2)
        Ut(Trivial,Nedge) = dKnown
        Ut(nonTrivial,Nedge) = Out(:,1)
     END IF
     IF (Ncuts .eq. 4) THEN
        Nedge = Nedge + 1
        Edge(Trivial,:,Nedge) = dKnown
        Edge(nonTrivial,:,Nedge) = Cut(:,3:4)
        Ut(Trivial,Nedge) = dKnown
        Ut(nonTrivial,Nedge) = Out(:,2)
        Outo = .true.
     END IF

   END SUBROUTINE StoreEdge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE Disconnected2D(N,Edge,Ut,Nsurfs,EdgeSet)

! Subroutine to check whether the surface cuts the cube once or several times

     IMPLICIT NONE
     INTEGER :: N, Nsurfs
     INTEGER,DIMENSION(5) :: EdgeSet
     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Edge
     REAL(KIND(0D0)),DIMENSION(0:2,2,12) :: Side
     REAL(KIND(0D0)),DIMENSION(0:2,12)   :: Ut, Ut2
     INTEGER :: i, j, Nside, Nedge

     Nsurfs = 1
     EdgeSet(1) = 1
     Side(:,:,1) = Edge(:,:,N)
     Ut2(:,1) = Ut(:,N)
     Nside = 1
     Nedge = N-1

     DO WHILE (Nedge .gt. 0)
        i = 1
        DO WHILE (ANY(Side(:,1,Nside).ne.Edge(:,1,i)).and.               &
                  ANY(Side(:,1,Nside).ne.Edge(:,2,i)).and.(i.le.Nedge))
           i = i+1
        END DO
        IF (i.le.Nedge) THEN
           IF (ALL(Side(:,1,Nside).eq.Edge(:,1,i))) THEN
              Side(:,2,Nside+1) = Edge(:,1,i)
              Side(:,1,Nside+1) = Edge(:,2,i)
           ELSE
              Side(:,:,Nside+1) = Edge(:,:,i)
           END IF
           Ut2(:,Nside+1) = Ut(:,i)
           Edge(:,:,i) = Edge(:,:,Nedge)
           Ut(:,i) = Ut(:,Nedge)
           EdgeSet(Nsurfs+1) = Nside+2
        ELSE
           Nsurfs = Nsurfs+1
           Side(:,:,Nside+1) = Edge(:,1:2,Nedge)
           Ut2(:,Nside+1) = Ut(:,Nedge)
        END IF
        Nside = Nside+1
        Nedge = Nedge-1
     END DO
     Edge(:,:,1:N) = Side(:,:,1:N)
     Ut(:,1:N) = Ut2(:,1:N)

   END SUBROUTINE Disconnected2D
