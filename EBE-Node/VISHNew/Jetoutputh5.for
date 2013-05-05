!***********************************************************************
      Subroutine setHydroFiles(XL_in, XH_in, DX_in, LSX_in, 
     &                         YL_in, YH_in, DY_in, LSY_in, 
     &                         Tau0_in, dTau_in, LST_in)
      
      Use HDF5
      Implicit none

      CHARACTER(LEN=20) :: H5hydroFilename = "results/JetData.h5" ! File name
      CHARACTER(LEN=8) :: groupname = "/Event" ! Group name
      Common /dataFile/ H5hydroFilename, groupname

      Integer :: XL_in, XH_in, YL_in, YH_in
      Double precision :: DX_in, DY_in, Tau0_in, dTau_in
      Integer :: LSX_in, LSY_in, LST_in

      Integer :: XL, XH, YL, YH
      Double precision :: DX, DY, Tau0, dTau
      Integer :: LSX, LSY, LST, LST_cur
      Integer:: XShift, YShift
      Common /hydroInfo/ XL, XH, DX, YL, YH, DY, Tau0, dTau
      Common /sparse/ XShift, LSX, YShift, LSY, LST, LST_cur
      
      Integer :: OutputViscousFlag = 1 ! Flag for whether to output shear stress tensor
      Common /OutputCtl/ OutputViscousFlag

      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER     ::   error ! Error flag

      XL = XL_in 
      XH = XH_in
      DX = DX_in
      YL = YL_in
      YH = YH_in
      DY = DY_in
      Tau0 = Tau0_in
      dTau = dTau_in

      LSX = LSX_in
      LSY = LSY_in
      LST = LST_in
      LST_cur = 0
      XShift = Abs(Mod(XL, LSX_in))
      YShift = Abs(Mod(YL, LSY_in))

      ! Initialize FORTRAN interface.
      CALL h5open_f(error)

      ! Create a new file using default properties.
      CALL h5fcreate_f(H5hydroFilename, H5F_ACC_TRUNC_F, file_id, error)

      ! Create group "Event" using absolute name.
      CALL h5gcreate_f(file_id, groupname, group_id, error)

      ! Write Attribute for group "Event"
      Call writeGroupattribute(group_id)

      ! Close the groups.
      CALL h5gclose_f(group_id, error)

      ! Terminate access to the file.
      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      CALL h5close_f(error)
      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine writeGroupattribute(group_id)
      Use HDF5
      Implicit none

      Integer :: XL, XH, YL, YH
      Double precision :: DX, DY, Tau0, dTau
      Integer :: LSX, LSY, LST, LST_cur
      Integer:: XShift, YShift
      Integer :: OutputViscousFlag

      Common /hydroInfo/ XL, XH, DX, YL, YH, DY, Tau0, dTau
      Common /sparse/ XShift, LSX, YShift, LSY, LST, LST_cur
      Common /OutputCtl/ OutputViscousFlag

      INTEGER(HID_T) :: group_id      ! Group identifier

      Call addGroupattributeInt(group_id, "XL",(XL+XShift)/LSX)
      Call addGroupattributeInt(group_id, "XH",(XH-XShift)/LSX)
      Call addGroupattributeInt(group_id, "YL",(YL+YShift)/LSY)
      Call addGroupattributeInt(group_id, "YH",(YH-YShift)/LSY)
      Call addGroupattributeDouble(group_id, "DX", DX*LSX)
      Call addGroupattributeDouble(group_id, "DY", DY*LSY)
      Call addGroupattributeDouble(group_id, "Tau0", Tau0)
      Call addGroupattributeDouble(group_id, "dTau", dTau*LST)
      Call addGroupattributeInt(group_id, "OutputViscousFlag", 
     &                             OutputViscousFlag)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine addGroupattributeInt(group_id, aname, avalue)
      Use HDF5
      Implicit none

      CHARACTER(LEN=*) :: aname       ! Attribute name
      Integer :: avalue      ! Attribute value
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      Integer(HID_T) :: atype_id
      
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      
      INTEGER     ::   error ! Error flag
     
      ! Create scalar data space for the attribute.
      CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create dataset attribute.
      CALL h5acreate_f(group_id, aname, H5T_NATIVE_INTEGER, aspace_id,
     &                  attr_id, error)
     
      ! Write the attribute data.
      CALL h5awrite_f(attr_id, H5T_NATIVE_INTEGER, avalue, adims, error)
     
      ! Close the attribute.
      CALL h5aclose_f(attr_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(aspace_id, error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine addGroupattributeDouble(group_id, aname, avalue)
      Use HDF5
      Implicit none

      CHARACTER(LEN=*) :: aname       ! Attribute name
      double precision :: avalue      ! Attribute value
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      
      INTEGER     ::   error ! Error flag
     
      ! Create scalar data space for the attribute.
      CALL h5screate_simple_f(arank, adims, aspace_id, error)

      ! Create dataset attribute.
      CALL h5acreate_f(group_id, aname, H5T_NATIVE_DOUBLE, aspace_id,
     &                  attr_id, error)
     
      ! Write the attribute data.
      CALL h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, avalue, adims, error)
     
      ! Close the attribute.
      CALL h5aclose_f(attr_id, error)

      ! Terminate access to the data space.
      CALL h5sclose_f(aspace_id, error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine writeHydroBlock(Time_id, Ed, Sd, P, Temp, Vx, Vy, 
     &   Pi00, Pi01, Pi02, Pi03, Pi11, Pi12, Pi13, Pi22, Pi23, Pi33, 
     &   BulkPi)

      Use HDF5
      Implicit none

      CHARACTER(LEN=20) :: H5hydroFilename ! File name
      CHARACTER(LEN=8) :: groupname        ! Group name
      CHARACTER(LEN=10) :: frameName       ! Group frame name
      Character(Len=4) :: frame_id_string
      Common /dataFile/ H5hydroFilename, groupname

      Integer :: XL, XH, YL, YH
      Double precision :: DX, DY, Tau0, dTau
      Integer :: LSX, LSY, LST, LST_cur
      Integer:: XShift, YShift
      Common /hydroInfo/ XL, XH, DX, YL, YH, DY, Tau0, dTau
      Common /sparse/ XShift, LSX, YShift, LSY, LST, LST_cur

      Integer :: OutputViscousFlag
      Common /OutputCtl/ OutputViscousFlag

      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: group_frame_id
      INTEGER(HID_T) :: dataset_id    ! Dataset identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier
      
      INTEGER     ::   error ! Error flag

      Integer :: Time_id, Frame_id
      Double precision, Dimension(XL:XH, YL:YH, 1:1) :: Ed, Sd, 
     &                  P, Temp, Vx, Vy
      Double precision, Dimension(XL:XH, YL:YH, 1:1):: Pi00, Pi01,
     &            Pi02, Pi03, Pi11, Pi12, Pi13, Pi22, Pi23, Pi33
      Double precision, Dimension(XL:XH, YL:YH, 1:1):: BulkPi
      
      INTEGER(HSIZE_T), DIMENSION(2) :: dims
      
      if(LST_cur /= 0) then  !no writing action
         LST_cur = LST_cur - 1
         return
      else   ! write to file and reset recurse idx
         LST_cur = LST - 1
      endif

      dims(1) = (XH - XL - 2*XShift)/LSX + 1
      dims(2) = (YH - YL - 2*YShift)/LSY + 1
      
      Frame_id = floor(DBLE(Time_id)/LST)
      write(unit=frame_id_string, fmt='(I4.4)') Frame_id
      frameName = "Frame_" // frame_id_string 

      ! Initialize FORTRAN interface.
      CALL h5open_f(error)

      ! Open an existing file.
      CALL h5fopen_f (H5hydroFilename, H5F_ACC_RDWR_F, file_id, error)
      
      ! Open an existing group in the specified file.
      CALL h5gopen_f(file_id, groupname, group_id, error)

      ! Create group "Frame_i" in group "Event" using relative name.
      CALL h5gcreate_f(group_id, frameName, group_frame_id, error)
      
      Call addGroupattributeDouble(group_frame_id, "Time", 
     &                             Tau0 + dTau*LST*Frame_id)

      ! Dump data into h5 file
      Call CSH5dumpBlockdata(group_frame_id, dims, "e", Ed)
      Call CSH5dumpBlockdata(group_frame_id, dims, "s", Sd)
      Call CSH5dumpBlockdata(group_frame_id, dims, "P", P)
      Call CSH5dumpBlockdata(group_frame_id, dims, "Temp", Temp)
      Call CSH5dumpBlockdata(group_frame_id, dims, "Vx", Vx)
      Call CSH5dumpBlockdata(group_frame_id, dims, "Vy", Vy)
      if(OutputViscousFlag .eq. 1) then
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi00", Pi00)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi01", Pi01)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi02", Pi02)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi03", Pi03)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi11", Pi11)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi12", Pi12)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi13", Pi13)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi22", Pi22)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi23", Pi23)
         Call CSH5dumpBlockdata(group_frame_id, dims, "Pi33", Pi33)
         Call CSH5dumpBlockdata(group_frame_id, dims, "BulkPi", BulkPi)
      endif

      ! Close the group.
      CALL h5gclose_f(group_frame_id, error)
      CALL h5gclose_f(group_id, error)

      ! Close the file.
      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      CALL h5close_f(error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      subroutine CSH5dumpBlockdata(group_id, dims, DatasetName, Dataset)
      Use HDF5
      Implicit none
      
      Character(Len=*) :: DatasetName
      Integer :: XL, XH, YL, YH
      Double precision :: DX, DY, Tau0, dTau
      Integer :: LSX, LSY, LST, LST_cur
      Integer:: XShift, YShift
      Common /hydroInfo/ XL, XH, DX, YL, YH, DY, Tau0, dTau
      Common /sparse/ XShift, LSX, YShift, LSY, LST, LST_cur

      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dataset_id    ! Dataset identifier
      INTEGER(HID_T) :: dataspace_id  ! Data space identifier

      INTEGER(HSIZE_T), DIMENSION(2) :: dims ! Datasets dimensions
      INTEGER(HSIZE_T), DIMENSION(2) :: dims_Cstyle ! Datasets dimensions
      Double precision, Dimension(XL:XH, YL:YH, 1:1) :: Dataset
      Double precision, Dimension(YL:YH, XL:XH, 1:1) :: Dataset_Cstyle

      Integer :: error
      Integer :: rank = 2

      Integer :: i, j
 
      ! convert shape of the matrix to C style for output
      dims_Cstyle(1) = dims(2)
      dims_Cstyle(2) = dims(1)
      do i = XL, XH, 1
         do j = YL, YH, 1
            Dataset_Cstyle(j, i, 1) = Dataset(i, j, 1)
         enddo
      enddo

      ! Create the data space for the first dataset.
      CALL h5screate_simple_f(rank, dims_Cstyle, dataspace_id, error)

      ! Create the dataset in group "Frame_i" with default properties.
      CALL h5dcreate_f(group_id, DatasetName, H5T_NATIVE_DOUBLE, 
     &                 dataspace_id, dataset_id, error)
      ! Write the first dataset.
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, 
     &  Dataset_Cstyle(YL+Yshift:YH-YShift:LSY,XL+XShift:XH-XShift:LSX,
     &                 1:1),dims, error)

      ! Close the dataspace for the dataset.
      CALL h5sclose_f(dataspace_id, error)

      ! Close the dataset.
      CALL h5dclose_f(dataset_id, error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydroFiles_initialEZ(H5hydroFilename_in)
      Implicit none
      CHARACTER(LEN=10) :: H5hydroFilename_in
      Call readHydroFiles_initial(H5hydroFilename_in, 0, 500)
      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydroFiles_initial(H5hydroFilename_in, 
     &                              InputViscousFlag_in, bufferSize_in)
      Use HDF5
      Implicit none

      CHARACTER(LEN=10) :: H5hydroFilename_in
      CHARACTER(LEN=10) :: hydroFileH5name ! File name
      CHARACTER(LEN=8) :: groupEventname = "/Event" ! Group name

      Common /fileInfo/ hydroFileH5name, groupEventname

      Integer :: InputViscousFlag_in, InputViscousFlag
      Common /InputCtl/ InputViscousFlag

      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: group_id      ! Group identifier

      Integer :: bufferSize_in
      INTEGER     ::   error ! Error flag

      hydroFileH5name = H5hydroFilename_in
      InputViscousFlag = InputViscousFlag_in

      ! Initialize FORTRAN interface.
      CALL h5open_f(error)

      ! open the hydro file
      CALL h5fopen_f (hydroFileH5name, H5F_ACC_RDWR_F, file_id, error)

      ! open a group
      CALL h5gopen_f(file_id, groupEventname, group_id, error)

      ! Read Attribute for group "Event"
      Call readHydrogridInfo(group_id)
      Call printHydrogridInfo()

      ! Read datasets from the file
      Call readHydroinfoBuffered_initialization(bufferSize_in)
      Call readHydroinfoBuffered_total(group_id)

      ! Close the groups.
      CALL h5gclose_f(group_id, error)

      ! Terminate access to the file.
      CALL h5fclose_f(file_id, error)

      ! Close FORTRAN interface.
      CALL h5close_f(error)
      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydrogridInfo(group_id)
      Use HDF5
      Implicit none

      CHARACTER(LEN=10) :: hydroFileH5name ! File name
      CHARACTER(LEN=8) :: groupEventname ! Group name

      Common /fileInfo/ hydroFileH5name, groupEventname

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes

      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes

      Integer :: IFlag
      Integer :: InputViscousFlag
      Common /InputCtl/ InputViscousFlag

      INTEGER(HID_T) :: group_id      ! Group identifier

      Integer :: error

      Call readH5Attribute_int(group_id, "XL", hydroGrid_XL)
      Call readH5Attribute_int(group_id, "XH", hydroGrid_XH)
      Call readH5Attribute_int(group_id, "YL", hydroGrid_YL)
      Call readH5Attribute_int(group_id, "YH", hydroGrid_YH)
      Call readH5Attribute_double(group_id, "DX", hydroGrid_DX)
      Call readH5Attribute_double(group_id, "DY", hydroGrid_DY)
      Call readH5Attribute_double(group_id, "Tau0", hydroGrid_Tau0)
      Call readH5Attribute_double(group_id, "dTau", hydroGrid_dTau)
      Call readH5Attribute_int(group_id, "OutputViscousFlag", IFlag)
      
      hydroGrid_X0 = hydroGrid_XL*hydroGrid_DX
      hydroGrid_Y0 = hydroGrid_YL*hydroGrid_DY
      
      Call h5gn_members_f(group_id, groupEventname, 
     &                    hydroGrid_numOfframes, error)
      hydroGrid_Taumax = hydroGrid_Tau0 
     &                   + (hydroGrid_numOfframes - 1)*hydroGrid_dTau

      InputViscousFlag = InputViscousFlag * IFlag

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine printHydrogridInfo()
      Implicit none
      
      CHARACTER(LEN=10) :: hydroFileH5name ! File name
      CHARACTER(LEN=8) :: groupEventname ! Group name

      Common /fileInfo/ hydroFileH5name, groupEventname

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes

      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes
      
      Integer :: InputViscousFlag
      Common /InputCtl/ InputViscousFlag

      write(*,'(A)')"--------------------------------------------------"
      write(*,'(A)')"--------------- hydro grid info ------------------"
      write(*,'(A)')"--------------------------------------------------"
      write(*,'(A, A)')"Filename : ", hydroFileH5name
      write(*,'(A, I5)') "XL = ", hydroGrid_XL
      write(*,'(A, I5)') "XH = ", hydroGrid_XH
      write(*,'(A, F5.3, A)') "DX = ", hydroGrid_DX, " fm"
      write(*,'(A, I5)') "YL = ", hydroGrid_YL
      write(*,'(A, I5)') "YH = ", hydroGrid_YH
      write(*,'(A, F5.3, A)') "DY = ", hydroGrid_DY, " fm"
      write(*,'(A, F5.3, A)') "Tau0 = ", hydroGrid_Tau0, " fm/c"
      write(*,'(A, F5.3, A)') "dTau = ", hydroGrid_dTau, " fm/c"
      write(*,'(A, I5)') "number of Frames = ", hydroGrid_numOfframes
      write(*,'(A, F7.3, A)') "Tau_max = ", hydroGrid_Taumax, " fm/c"
      write(*,"(A, I5)")"Readin viscous information : ",InputViscousFlag
      write(*,'(A)')"--------------------------------------------------"

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readH5Attribute_int(group_id, aname, avalue)
      Use HDF5
      Implicit none

      CHARACTER(LEN=*) :: aname       ! Attribute name
      Integer :: avalue      ! Attribute value
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      Integer(HID_T) :: atype_id
      
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      
      INTEGER     ::   error ! Error flag
      
      ! open an attribute
      Call h5aopen_name_f(group_id, aname, attr_id, error)
      ! read an attribute
      Call h5aread_f(attr_id, H5T_NATIVE_INTEGER, avalue, adims, error)
      ! close an attribute
      Call h5aclose_f(attr_id, error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readH5Attribute_double(group_id, aname, avalue)
      Use HDF5
      Implicit none

      CHARACTER(LEN=*) :: aname       ! Attribute name
      Double precision :: avalue      ! Attribute value
      INTEGER(HID_T) :: group_id      ! Group identifier

      INTEGER(HID_T) :: attr_id       ! Attribute identifier
      INTEGER(HID_T) :: aspace_id     ! Attribute Dataspace identifier
      Integer(HID_T) :: atype_id
      
      INTEGER(HSIZE_T), DIMENSION(1) :: adims = (/1/) ! Attribute dimension
      INTEGER     ::   arank = 1                      ! Attribure rank
      
      INTEGER     ::   error ! Error flag
      
      ! open an attribute
      Call h5aopen_name_f(group_id, aname, attr_id, error)
      ! read an attribute
      Call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, avalue, adims, error)
      ! close an attribute
      Call h5aclose_f(attr_id, error)

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydroinfoBuffered_initialization(bufferSize_in)
      Implicit None
      
      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes

      Integer:: bufferSize
      Integer:: bufferSize_in

      ! the last index is the buffer "layer" index that goes from 1 to bufferSize
      Double Precision, Pointer::
     & eM(:,:,:), PM(:,:,:), sM(:,:,:), TM(:,:,:), vxM(:,:,:), 
     & vyM(:,:,:),
     & pi00M(:,:,:), pi01M(:,:,:), pi02M(:,:,:), pi03M(:,:,:), 
     & pi11M(:,:,:), pi12M(:,:,:), pi13M(:,:,:), pi22M(:,:,:),
     & pi23M(:,:,:), pi33M(:,:,:), BulkPiM(:,:,:)

      Common /bufferedData/ bufferSize, 
     &  eM, PM, sM, TM, vxM, vyM, 
     &  pi00M, pi01M, pi02M, pi03M, pi11M, pi12M, pi13M, pi22M,
     &  pi23M, pi33M, BulkPiM

      bufferSize = bufferSize_in
      Allocate(eM(hydroGrid_XL:hydroGrid_XH, hydroGrid_YL:hydroGrid_YH, 
     &            1:bufferSize))
      Allocate(PM(hydroGrid_XL:hydroGrid_XH, hydroGrid_YL:hydroGrid_YH, 
     &            1:bufferSize))
      Allocate(sM(hydroGrid_XL:hydroGrid_XH, hydroGrid_YL:hydroGrid_YH, 
     &            1:bufferSize))
      Allocate(TM(hydroGrid_XL:hydroGrid_XH, hydroGrid_YL:hydroGrid_YH, 
     &            1:bufferSize))
      Allocate(vxM(hydroGrid_XL:hydroGrid_XH, hydroGrid_YL:hydroGrid_YH,
     &             1:bufferSize))
      Allocate(vyM(hydroGrid_XL:hydroGrid_XH, hydroGrid_YL:hydroGrid_YH,
     &             1:bufferSize))
      Allocate(pi00M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi01M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi02M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi03M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi11M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi12M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi13M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi22M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi23M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(pi33M(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))
      Allocate(BulkPiM(hydroGrid_XL:hydroGrid_XH, 
     &               hydroGrid_YL:hydroGrid_YH, 1:bufferSize))

      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydroinfoBuffered_total(group_id)
      Use HDF5
      Implicit None

      CHARACTER(LEN=10) :: hydroFileH5name ! File name
      CHARACTER(LEN=8) :: groupEventname ! Group name
      Common /fileInfo/ hydroFileH5name, groupEventname
      
      CHARACTER(LEN=10) :: frameName       ! Group frame name
      Character(Len=4) :: frame_id_string
      
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: groupFrame_id ! Group identifier
      Integer :: error

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes

      Integer :: InputViscousFlag
      Common /InputCtl/ InputViscousFlag

      Double Precision, Dimension(hydroGrid_XL:hydroGrid_XH,
     &  hydroGrid_YL:hydroGrid_YH, 1:1) :: Ed, P, Sd, Temp, Vx, Vy,
     &  pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33,
     &  BulkPi

      Double Precision :: Time

      Integer:: bufferSize

      ! the last index is the buffer "layer" index that goes from 1 to bufferSize
      Double Precision, Pointer::
     & eM(:,:,:), PM(:,:,:), sM(:,:,:), TM(:,:,:), vxM(:,:,:), 
     & vyM(:,:,:),
     & pi00M(:,:,:), pi01M(:,:,:), pi02M(:,:,:), pi03M(:,:,:), 
     & pi11M(:,:,:), pi12M(:,:,:), pi13M(:,:,:), pi22M(:,:,:),
     & pi23M(:,:,:), pi33M(:,:,:), BulkPiM(:,:,:)

      Common /bufferedData/ bufferSize, 
     &  eM, PM, sM, TM, vxM, vyM, 
     &  pi00M, pi01M, pi02M, pi03M, pi11M, pi12M, pi13M, pi22M,
     &  pi23M, pi33M, BulkPiM
       
      Integer :: J
      Integer :: xidx

      if(bufferSize .lt. hydroGrid_numOfframes) then
         write(*,*) "BufferSize is too small, increase it to at least", 
     &              hydroGrid_numOfframes
         stop
      endif
      Do J = 1, hydroGrid_numOfframes
        write(unit=frame_id_string, fmt='(I4.4)') J-1
        frameName = "Frame_" // frame_id_string 
        CALL h5gopen_f(group_id, frameName, groupFrame_id, error)
        
        Call readH5Dataset_double(groupFrame_id, "e", Ed)
        Call readH5Dataset_double(groupFrame_id, "s", Sd)
        Call readH5Dataset_double(groupFrame_id, "P", P)
        Call readH5Dataset_double(groupFrame_id, "Temp", Temp)
        Call readH5Dataset_double(groupFrame_id, "Vx", Vx)
        Call readH5Dataset_double(groupFrame_id, "Vy", Vy)
        if(InputViscousFlag .eq. 1) then
          Call readH5Dataset_double(groupFrame_id, "Pi00", pi00)
          Call readH5Dataset_double(groupFrame_id, "Pi01", pi01)
          Call readH5Dataset_double(groupFrame_id, "Pi02", pi02)
          Call readH5Dataset_double(groupFrame_id, "Pi03", pi03)
          Call readH5Dataset_double(groupFrame_id, "Pi11", pi11)
          Call readH5Dataset_double(groupFrame_id, "Pi12", pi12)
          Call readH5Dataset_double(groupFrame_id, "Pi13", pi13)
          Call readH5Dataset_double(groupFrame_id, "Pi22", pi22)
          Call readH5Dataset_double(groupFrame_id, "Pi23", pi23)
          Call readH5Dataset_double(groupFrame_id, "Pi33", pi33)
          Call readH5Dataset_double(groupFrame_id, "BulkPi", BulkPi)
        endif
        ! write to the buffer:
        eM(:,:,J) = Ed(:,:,1)
        PM(:,:,J) = P(:,:,1)
        sM(:,:,J) = Sd(:,:,1)
        TM(:,:,J) = Temp(:,:,1)
        vxM(:,:,J) = Vx(:,:,1)
        vyM(:,:,J) = Vy(:,:,1)
        if(InputViscousFlag .eq. 1) then
          pi00M(:,:,J) = pi00(:,:,1)
          pi01M(:,:,J) = pi01(:,:,1)
          pi02M(:,:,J) = pi02(:,:,1)
          pi03M(:,:,J) = pi03(:,:,1)
          pi11M(:,:,J) = pi11(:,:,1)
          pi12M(:,:,J) = pi12(:,:,1)
          pi13M(:,:,J) = pi13(:,:,1)
          pi22M(:,:,J) = pi22(:,:,1)
          pi23M(:,:,J) = pi23(:,:,1)
          pi33M(:,:,J) = pi33(:,:,1)
          BulkPiM(:,:,J) = BulkPi(:,:,1)
        endif
      
      ! Close the groups.
      CALL h5gclose_f(groupFrame_id, error)
      EndDo

      End
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readH5Dataset_double(group_id, datasetName, dset_data)
      Use HDF5
      Implicit None
      
      CHARACTER(LEN=*) :: datasetName       ! Dataset name
      INTEGER(HID_T) :: group_id      ! Group identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      
      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes
      
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims_Cstyle
      Double precision, Dimension(hydroGrid_XL:hydroGrid_XH, 
     &  hydroGrid_YL:hydroGrid_YH, 1:1) :: dset_data
      Double precision, Dimension(hydroGrid_YL:hydroGrid_YH, 
     &  hydroGrid_XL:hydroGrid_XH, 1:1) :: dset_data_Cstyle
      Integer :: error
      Integer :: i, j
      
      ! read in data matrix assuming in C style 
      ! need to perform transpose to convert it into fortran style
      data_dims(1) = hydroGrid_XH - hydroGrid_XL + 1
      data_dims(2) = hydroGrid_YH - hydroGrid_YL + 1
      data_dims_Cstyle(1) = data_dims(2)
      data_dims_Cstyle(2) = data_dims(1)

      ! Open an existing dataset.
      CALL h5dopen_f(group_id, datasetName, dset_id, error)

      ! Read the dataset.
      CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dset_data_Cstyle, 
     &               data_dims_Cstyle, error)
      do i = hydroGrid_XL, hydroGrid_XH, 1
        do j = hydroGrid_YL, hydroGrid_YH, 1
          dset_data(i, j, 1) = dset_data_Cstyle(j, i, 1)
        enddo
      enddo
     
      ! Close the dataset.
      CALL h5dclose_f(dset_id, error)
      
      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydroinfoBuffered_ideal(tau,x,y,e,p,s,T,vx,vy)
!     Return infos from hydro data file
!     -- tau,x,y: coordinates (in)
!     -- e,p,s,T,vx,vy: infos (out)
      Implicit None

      Double Precision:: tau,x,y,e,p,s,T,vx,vy

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes

      Double Precision, Dimension(1:2,1:2) :: Ed1,P1,Sd1,Temp1,VxB1,VyB1
      Double Precision, Dimension(1:2,1:2) :: Ed2,P2,Sd2,Temp2,VxB2,VyB2

      Integer:: tauI ! tau should lie between Tau0+dTau*tauI and Tau0+dTau*(tauI+1)
      Double Precision:: tauInc ! tau=Tau0+tauI1*dTau+Tau_inc*dTau, inc for increament (normalized to 0~1)

      Integer:: xi, yi ! similar to tau
      Double Precision:: xInc, yInc ! x=xi*DX+xInc*DX

      Double Precision:: var1 ! temporary variables

      ! first deal with tau
      if (tau < hydroGrid_Tau0 .AND. tau > hydroGrid_Tau0 - 1D-15) then
         tau = hydroGrid_Tau0
      endif
      var1 = (tau-hydroGrid_Tau0)/hydroGrid_dTau
      tauI = Floor(var1)
      tauInc = var1-tauI

      ! then x and y
      var1 = (x)/hydroGrid_DX
      xi = Floor(var1)
      xInc = var1-xi
      var1 = (y)/hydroGrid_DY
      yi = Floor(var1)
      yInc = var1-yi

!     For debug:
!      Print*, "tau, Tau0, dTau, tauI, tauInc=",
!     &        tau, Tau0, dTau, tauI, tauInc
!      Print*, "x,xi,xInc,y,yi,yInc=",
!     &        x,xi,xInc,y,yi,yInc

      Call readHydroBlockBufferedOrdered_ideal
     &    (tauI+1,xi,yi,Ed1,P1,Sd1,Temp1,VxB1,VyB1) ! tauI+1: tauI=0 <-> 1st block
      Call readHydroBlockBufferedOrdered_ideal
     &    (tauI+2,xi,yi,Ed2,P2,Sd2,Temp2,VxB2,VyB2)

      Call cubeInterp(xInc,yInc,tauInc,e,
     &  Ed1(1,1),Ed1(1+1,1),Ed1(1,1+1),Ed1(1+1,1+1),
     &  Ed2(1,1),Ed2(1+1,1),Ed2(1,1+1),Ed2(1+1,1+1))

      Call cubeInterp(xInc,yInc,tauInc,p,
     &  P1(1,1),P1(1+1,1),P1(1,1+1),P1(1+1,1+1),
     &  P2(1,1),P2(1+1,1),P2(1,1+1),P2(1+1,1+1))

      Call cubeInterp(xInc,yInc,tauInc,s,
     &  Sd1(1,1),Sd1(1+1,1),Sd1(1,1+1),Sd1(1+1,1+1),
     &  Sd2(1,1),Sd2(1+1,1),Sd2(1,1+1),Sd2(1+1,1+1))

      Call cubeInterp(xInc,yInc,tauInc,T,
     &  Temp1(1,1),Temp1(1+1,1),
     &  Temp1(1,1+1),Temp1(1+1,1+1),
     &  Temp2(1,1),Temp2(1+1,1),
     &  Temp2(1,1+1),Temp2(1+1,1+1))

      Call cubeInterp(xInc,yInc,tauInc,vx,
     &  VxB1(1,1),VxB1(1+1,1),VxB1(1,1+1),VxB1(1+1,1+1),
     &  VxB2(1,1),VxB2(1+1,1),VxB2(1,1+1),VxB2(1+1,1+1))

      Call cubeInterp(xInc,yInc,tauInc,vy,
     &  VyB1(1,1),VyB1(1+1,1),VyB1(1,1+1),VyB1(1+1,1+1),
     &  VyB2(1,1),VyB2(1+1,1),VyB2(1,1+1),VyB2(1+1,1+1))

      End Subroutine
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine readHydroBlockBufferedOrdered_ideal
     &    (idxTau,idxX,idxY,Ed22,P22,Sd22,Temp22,VxB22,VyB22)
!     Read only the 2x2 block froEdm the buffer or file at
!     tauI=idxTau, index_x=idxX, and index_y=idxY
!     This version assumes all the required hydro info are in the buffer
      Implicit None

      Integer :: idxTau, idxX, idxY

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes
      Double Precision, Dimension(hydroGrid_XL:hydroGrid_XH,
     & hydroGrid_YL:hydroGrid_YH, 1:1) :: Ed,P,Sd,Temp,VxB,VyB

      Double Precision, Dimension(1:2,1:2) ::
     &    Ed22,P22,Sd22,Temp22,VxB22,VyB22

      Integer:: bufferSize

      ! the last index is the buffer "layer" index that goes from 1 to bufferSize
      Double Precision, Pointer::
     & eM(:,:,:), PM(:,:,:), sM(:,:,:), TM(:,:,:), vxM(:,:,:), 
     & vyM(:,:,:),
     & pi00M(:,:,:), pi01M(:,:,:), pi02M(:,:,:), pi03M(:,:,:), 
     & pi11M(:,:,:), pi12M(:,:,:), pi13M(:,:,:), pi22M(:,:,:),
     & pi23M(:,:,:), pi33M(:,:,:), BulkPiM(:,:,:)

      Common /bufferedData/ bufferSize, 
     &  eM, PM, sM, TM, vxM, vyM, 
     &  pi00M, pi01M, pi02M, pi03M, pi11M, pi12M, pi13M, pi22M,
     &  pi23M, pi33M, BulkPiM

      If (idxTau < 0 .OR. idxTau > bufferSize .or. 
     &    idxTau > hydroGrid_numOfframes) Then
        Ed22(:,:) = 0D0
        P22(:,:) = 0D0
        Sd22(:,:) = 0D0
        Temp22(:,:) = 0D0
        VxB22(:,:) = 0D0
        VyB22(:,:) = 0D0
        Return
      End If

      If (idxX < hydroGrid_XL .OR. idxX > hydroGrid_XH) Then
        Ed22(:,:) = 0D0
        P22(:,:) = 0D0
        Sd22(:,:) = 0D0
        Temp22(:,:) = 0D0
        VxB22(:,:) = 0D0
        VyB22(:,:) = 0D0
        Return
      End If

      If (idxY < hydroGrid_YL .OR. idxY > hydroGrid_YH) Then
        Ed22(:,:) = 0D0
        P22(:,:) = 0D0
        Sd22(:,:) = 0D0
        Temp22(:,:) = 0D0
        VxB22(:,:) = 0D0
        VyB22(:,:) = 0D0
        Return
      End If

      ! read it from buffer:
      Ed22(:,:) = eM(idxX:idxX+1,idxY:idxY+1,idxTau)
      P22(:,:) = PM(idxX:idxX+1,idxY:idxY+1,idxTau)
      Sd22(:,:) = sM(idxX:idxX+1,idxY:idxY+1,idxTau)
      Temp22(:,:) = TM(idxX:idxX+1,idxY:idxY+1,idxTau)
      VxB22(:,:) = vxM(idxX:idxX+1,idxY:idxY+1,idxTau)
      VyB22(:,:) = vyM(idxX:idxX+1,idxY:idxY+1,idxTau)

      End Subroutine
!------------------------------------------------------------------------

!***********************************************************************
      Subroutine cubeInterp(x,y,z,Axyz,
     &                      A000,A100,A010,A110,A001,A101,A011,A111)
! Perform a 3d interpolation. The known data are A### located at the 8 corners,
! labels using the xyz order. Therefore A000 is value at the origin and A010
! is the value at (x=0,y=1,z=0). Note that the coordinate (x,y,z) must be
! constrained to the unit cube. Axyz is the return value.

      Implicit None
      Double Precision :: x,y,z,Axyz
      Double Precision :: A000,A100,A010,A110,A001,A101,A011,A111

      Axyz = A000*(1-x)*(1-y)*(1-z) + A100*x*(1-y)*(1-z) +
     &    A010*(1-x)*y*(1-z) + A001*(1-x)*(1-y)*z +
     &    A101*x*(1-y)*z + A011*(1-x)*y*z +
     &    A110*x*y*(1-z) + A111*x*y*z

      !If (debug == 1) Then
      !  Print *, "cubeInterp"
      !  Print *, "x,y,z=",x,y,z
      !  Print *, "A000,A100,A010,A110=",A000,A100,A010,A110
      !  Print *, "A001,A101,A011,A111=",A001,A101,A011,A111
      !  Print *, "Axyz=", Axyz
      !End If

      RETURN
      END
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine outputPlaintxtHuichaoFormat()
      Implicit None

      Double Precision:: tau,x,y,e,p,s,T,vx,vy

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes
      Integer :: Itau, Ix, Iy
      double precision :: outputDtau, outputDx, outputDy
      double precision :: outputXL, outputXH, outputYL, outputYH
      double precision :: outputTau0, outputTauMax
      integer :: outputNX, outputNY, outputNtau
      double precision :: local_Tau, local_x, local_y

      double precision :: hbarC

      hbarC = 0.19733

      outputDtau = 0.1d0
      outputDx = 0.5d0
      outputDy = 0.5d0
      outputXL = -11.0d0
      outputXH = 11.0d0
      outputYL = -11.0d0
      outputYH = 11.0d0
      outputTau0 = hydroGrid_Tau0
      outputTauMax = hydroGrid_Taumax


      outputNX = floor((outputXH - outputXL)/outputDx) + 1
      outputNY = floor((outputYH - outputYL)/outputDy) + 1
      outputNtau = floor((outputTauMax - outputTau0)/outputDtau) + 1

      open(20, FILE='hydroinfoPlaintxtHuichaoFormat.dat', 
     &     FORM='FORMATTED', STATUS='REPLACE')
      
      do Itau = 1, outputNtau, 1
        local_Tau = outputTau0 + (Itau - 1)*outputDtau
        do Ix = 1, outputNX, 1
          local_x = outputXL + (Ix - 1)*outputDx
          do Iy = 1, outputNY, 1
             local_y = outputYL + (Iy - 1)*outputDy
             call readHydroinfoBuffered_ideal(local_Tau, local_x, 
     &              local_y, e, p, s, T, vx, vy)
             write(20, '(3F10.2, 4F20.10)') local_x, local_y, local_Tau,
     &              e/hbarC, T/hbarC, vx, vy
          enddo
        enddo
      enddo
      close(20)

      return
      end
!-----------------------------------------------------------------------

!***********************************************************************
      Subroutine getJetDeltaTauMax(x0,y0,dirX,dirY,cutT,step,deltaTau)
!     Return the max possible deltaTau that determines the length of the
!     path of a jet positioned at (x0,y0) with direction (dirX,dirY). The
!     jet is assumed to travel at speed of light and the deltaTau value
!     is determined by the time the jet leaves the hydro data x-y-tau cube.
!     deltaTau is measured relative to tau0.
!     Must be called after getHydroDataRecLength.
!     After the crude maximum is determined, it is reduced by checking,
!     from the maximum possible value deltaTau, "step" by "step", if the
!     temperature is larger than cutT. It keeps shrinking if the readed
!     temperature is less than cutT (GeV)
!     Note that this method is slower than getJetDeltaTauMaxOld.

      Implicit None

      Integer :: hydroGrid_XL, hydroGrid_XH, hydroGrid_YL, hydroGrid_YH
      Double precision :: hydroGrid_X0, hydroGrid_Y0
      Double precision :: hydroGrid_DX, hydroGrid_DY
      Double precision :: hydroGrid_Tau0, hydroGrid_dTau
      Double precision :: hydroGrid_Taumax
      Integer :: hydroGrid_numOfframes
      Common /hydroGridinfo/ hydroGrid_XL, hydroGrid_XH, 
     &                       hydroGrid_X0, hydroGrid_DX, 
     &                       hydroGrid_YL, hydroGrid_YH, 
     &                       hydroGrid_Y0, hydroGrid_DY, 
     &                       hydroGrid_Tau0, hydroGrid_dTau,
     &                       hydroGrid_Taumax,
     &                       hydroGrid_numOfframes

      Double Precision:: x0,y0,dirX,dirY,deltaTau
      Double Precision:: rXL,rXH,rYL,rYH
      Double Precision:: dirNorm
      Double Precision:: cutT,step,jetLength,e,p,s,T,vx,vy

      ! first, calculate boundaries
      rXL = hydroGrid_XL*hydroGrid_DX
      rXH = hydroGrid_XH*hydroGrid_DX
      rYL = hydroGrid_YL*hydroGrid_DY
      rYH = hydroGrid_YH*hydroGrid_DY

      ! next, find intersections and length
      dirNorm = Sqrt(dirX*dirX+dirY*dirY)
      
      If (dirX>=0) Then ! point to the right
        If (dirY>=0) Then ! point to the top
          If (dirY*(rXH-x0)>dirX*(rYH-y0)) Then
            ! intercept with the top
            deltaTau = (rYH-y0)/dirY*dirNorm
          Else
            ! intercept with the right
            deltaTau = (rXH-x0)/dirX*dirNorm
          EndIf
        Else ! point to the bottom
          If ((-dirY)*(rXH-x0)>dirX*(y0-rYL)) Then
            ! intercept with the bottom
            deltaTau = (y0-rYL)/(-dirY)*dirNorm
          Else
            ! intercept with the right
            deltaTau = (rXH-x0)/dirX*dirNorm
          EndIf
        EndIf
      Else ! point to the left
        If (dirY>=0) Then ! point to the top
          If (dirY*(x0-rXL)>(-dirX)*(rYH-y0)) Then
            ! intercept with the top
            deltaTau = (rYH-y0)/dirY*dirNorm
          Else
            ! intercept with the left
            deltaTau = (x0-rXL)/(-dirX)*dirNorm
          EndIf
        Else ! point to the bottom
          If ((-dirY)*(x0-rXL)>(-dirX)*(y0-rYL)) Then
            ! intercept with the bottom
            deltaTau = (y0-rYL)/(-dirY)*dirNorm
          Else
            ! intercept with the left
            deltaTau = (x0-rXL)/(-dirX)*dirNorm
          EndIf
        EndIf      
      EndIf
      
      ! constrain from tau direction)
      deltaTau = Min(deltaTau, hydroGrid_dTau*(hydroGrid_numOfframes-1))

      ! try to shrink this length
      jetLength = deltaTau
      Do
        If (jetLength<0D0) Then
          deltaTau=0D0
          Exit
        EndIf
        if(jetLength < hydroGrid_Tau0) then
          call readHydroinfoBuffered_ideal(hydroGrid_Tau0,
     &                        x0+dirX/dirNorm*(jetLength), !
     &                        y0+dirY/dirNorm*(jetLength), !
     &                        e,p,s,T,vx,vy)
        else
          Call readHydroinfoBuffered_ideal(jetLength,
     &                        x0+dirX/dirNorm*(jetLength), !
     &                        y0+dirY/dirNorm*(jetLength), !
     &                        e,p,s,T,vx,vy)
        endif
        If (T>=cutT) Then
          deltaTau = jetLength
          Exit
        EndIf
        jetLength=jetLength-step
      End Do

      End Subroutine
!------------------------------------------------------------------------

!***********************************************************************
      Subroutine getJetavgLength_shell(cutT, Jet_maxlength, 
     &               Jet_avglength, Jet_avglength_in, Jet_avglength_out)
      Implicit none
      double precision :: cutT
      double precision :: Jet_maxlength, Jet_avglength, 
     &                    Jet_avglength_in, Jet_avglength_out

      call getJetavgLength_initial(cutT)
      call getJetavgLengthfromMaxlength(Jet_maxlength, Jet_avglength, 
     &                              Jet_avglength_in, Jet_avglength_out)
      return
      end
!------------------------------------------------------------------------

!***********************************************************************
      Subroutine getJetavgLength_initial(cutT)
!     For a given cut-off temperature, this subroutine performs the calculations 
!     for test particles' path lengths inside hydro medium and build up tables 
!     to store the averaged path lengths for test particles with some given 
!     maximum allowed jet path lengths.
!     Test particles are weighted by the local binary collision probabilities
!     that provided by hydro simulation.
      Implicit none
      double precision :: cutT
      Integer, Parameter :: NXPhy0=-130,NYPhy0=-130   !Physical initial grid size
      Integer, Parameter :: NXPhy=130,NYPhy=130       !Physical end grid size

      double precision, Parameter :: PI=3.1415926585d0
      Integer :: I=1, J=1, K=1, L=1     !loop variables

      double precision :: Binary_profile(NXPhy0:NXPhy, NYPhy0:NYPhy)

C=====variable for test particles====================================
      integer, parameter :: nphi = 100
      integer, parameter :: jetnum = 6812100     !Maximum number of testing jets
      double precision, parameter :: eps = 1D-4
      
      double precision :: x0, y0, phi0   !the origin of the jet position
      double precision :: dx, dy, dphi!the step length
      double precision :: dirx, diry
      double precision :: weight_j(jetnum)    !weight for jet production
      double precision :: phi_j(jetnum)
      double precision :: D_j(jetnum)    !record the distence the jets travels when it interact with medium

      Integer :: Idx                     !index of the testing jet
      Integer :: number_of_jets
      double precision :: maxTau

      double precision :: phi_boundary  !the maximum in plane jets emission angle
      integer, parameter :: n_max_pathlength = 100
      double precision :: max_pathlength0 = 1.0d0
      double precision :: dmax_pathlength0 = 0.2d0
      double precision :: max_pathlength(n_max_pathlength)
      double precision :: avg_pathlength(n_max_pathlength)
      double precision :: avg_pathlength_in(n_max_pathlength)
      double precision :: avg_pathlength_out(n_max_pathlength)
      double precision :: norm_pathlength(n_max_pathlength)
      double precision :: norm_pathlength_in(n_max_pathlength)
      double precision :: norm_pathlength_out(n_max_pathlength)

      common /testparticle_pathlength/avg_pathlength, avg_pathlength_in,
     &                                avg_pathlength_out, max_pathlength
      common /testparticle_pathlength_coeff/max_pathlength0, 
     &                                      dmax_pathlength0

C=====variable for jet quenching output end================================

      open(1, file='TATB_fromSd_order_2_block.dat', status='old')
      do I = NXPhy0, NXPhy, 1
         read(1, *) (Binary_profile(I, J), J = NYPhy0, NYPhy, 1)
      enddo
      close(1)
      
      dx = 0.1d0
      dy = 0.1d0
      dphi = 2.0d0*PI/nphi
      phi_boundary = PI/12.

      Idx = 1
      do 100 I = NXPhy0, NXPhy, 2
      do 100 J = NYPhy0, NYPhy, 2
         if(Binary_profile(I, J) .ge. eps) then
            x0 = I*dx
            y0 = J*dy
            do K = 1, nphi, 1 
               phi0 = 0. + (K - 1)*dphi
               dirx = cos(phi0)
               diry = sin(phi0)
               Call getJetDeltaTauMax(x0, y0, dirx, diry, cutT,
     &                                0.05d0,maxTau)
               D_j(Idx) = maxTau
               weight_j(Idx) = Binary_profile(I, J)
               phi_j(Idx) = phi0
               Idx = Idx + 1
            enddo
         endif
100   continue

      avg_pathlength = 0.0d0
      avg_pathlength_in = 0.0d0
      avg_pathlength_out = 0.0d0
      norm_pathlength = 0.0d0
      norm_pathlength_in = 0.0d0
      norm_pathlength_out = 0.0d0
      number_of_jets = Idx - 1
      do I = 1, number_of_jets
         do J = 1, n_max_pathlength
            max_pathlength(J) = max_pathlength0
     &                        + (J-1)*dmax_pathlength0
            if(D_j(I) .le. max_pathlength(J)) then
               avg_pathlength(J) = avg_pathlength(J) 
     &                           + D_j(I)*weight_j(I)
               norm_pathlength(J) = norm_pathlength(J) + weight_j(I)
               if((phi_j(I) .le. phi_boundary) .or.
     &            (phi_j(I) .ge. (2*PI - phi_boundary)) .or.
     &            (phi_j(I) .ge. (PI - phi_boundary) .and.
     &             phi_j(I) .le. (PI + phi_boundary)) ) then
                 avg_pathlength_in(J) = avg_pathlength_in(J)
     &                                + D_j(I)*weight_j(I)
                 norm_pathlength_in(J) = norm_pathlength_in(J) 
     &                                 + weight_j(I)
               else if(((phi_j(I) .ge. (PI/2. - phi_boundary)) .and.
     &                  (phi_j(I) .le. (PI/2. + phi_boundary))) .or.
     &                 (phi_j(I) .ge. (3.*PI/2. - phi_boundary) .and.
     &                  phi_j(I) .le. (3.*PI/2. + phi_boundary))) then
                 avg_pathlength_out(J) = avg_pathlength_out(J)
     &                                 + D_j(I)*weight_j(I)
                 norm_pathlength_out(J) = norm_pathlength_out(J) 
     &                                  + weight_j(I)
               endif
            endif
         enddo
      enddo
      do I = 1, n_max_pathlength
         avg_pathlength(I) = avg_pathlength(I)/norm_pathlength(I)
         avg_pathlength_in(I) = avg_pathlength_in(I)
     &                         /norm_pathlength_in(I)
         avg_pathlength_out(I) = avg_pathlength_out(I)
     &                         /norm_pathlength_out(I)
      enddo

      return
      end
!-----------------------------------------------------------------------

!***********************************************************************
      subroutine getJetavgLengthfromMaxlength(Maxlength, avglength, 
     &                  avglength_in, avglength_out)
!     Return the average test particles' path length (total mean, in-plane, 
!     out-of-plane) for a given maximum allowed path length. 
!     This subroutine must be called after the main program called 
!     getJetavgLength_initial to initialize the calculations.
!     This subroutine only perform a linear interpolation between the table
!     of averaged length that calculated in getJetavgLength_initial.
      Implicit none
      double precision :: Maxlength, avglength, avglength_in, 
     &                    avglength_out

      integer :: Idx
      double precision :: dx
      integer, parameter :: n_max_pathlength = 100
      double precision :: max_pathlength0
      double precision :: dmax_pathlength0
      double precision :: max_pathlength(n_max_pathlength)
      double precision :: avg_pathlength(n_max_pathlength)
      double precision :: avg_pathlength_in(n_max_pathlength)
      double precision :: avg_pathlength_out(n_max_pathlength)
      
      common /testparticle_pathlength/avg_pathlength, avg_pathlength_in,
     &                                avg_pathlength_out, max_pathlength
      common /testparticle_pathlength_coeff/max_pathlength0, 
     &                                      dmax_pathlength0
      if(Maxlength .lt. 0.0d0) then
         print*, "Max. allowed Jet length is < 0!"
         stop
      endif

      if(Maxlength .le. max_pathlength(1)) then
         dx = Maxlength
         avglength = avg_pathlength(1)/max_pathlength(1)*dx
         avglength_in = avg_pathlength_in(1)/max_pathlength(1)*dx
         avglength_out = avg_pathlength_out(1)/max_pathlength(1)*dx
      else if(Maxlength .ge. max_pathlength(n_max_pathlength)) then
         avglength = avg_pathlength(n_max_pathlength)
         avglength_in = avg_pathlength_in(n_max_pathlength)
         avglength_out = avg_pathlength_out(n_max_pathlength)
      else
         Idx = floor((Maxlength - max_pathlength0)/dmax_pathlength0)+1
         dx = Maxlength - max_pathlength(Idx)

         avglength = (avg_pathlength(Idx+1) - avg_pathlength(Idx))/
     &            dmax_pathlength0*dx + avg_pathlength(Idx)
         avglength_in = (avg_pathlength_in(Idx+1)
     &                   - avg_pathlength_in(Idx))/dmax_pathlength0*dx 
     &                  + avg_pathlength_in(Idx)
         avglength_out = (avg_pathlength_out(Idx+1)
     &                    - avg_pathlength_out(Idx))/dmax_pathlength0*dx
     &                   + avg_pathlength_out(Idx)
      endif

      return
      end
!-----------------------------------------------------------------------
