! Ver 1.2: EOS extended to infinite energy density by extrapolation
      Subroutine InputRegulatedEOS
      Implicit none

!=======list of parameters===================================================
      Integer, Parameter :: RegEOSdatasize = 155500  !converted EOS table data size
      Integer, Parameter :: RegEOSMudatasize = 501   !converted EOS Mu table data size
      Integer, Parameter :: IMax_Mu = 40        !maximum allowed stable particles for partically chemical equilibrium EOS
!=======list of parameters end===============================================

!=======declear variables for reading equation of state======================
      Integer :: I, J, K, L, M, N    !loop variables
      double precision :: EOSe0=0.0d0   !lowest energy density
      double precision :: EOSde=0.0d0   !spacing of energy density
      Integer :: EOSne=RegEOSdatasize   !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: ee
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2

      Integer :: Inumparticle       !number of stable particles in PCE, 0 for chemical equilibrium EOS
      Integer :: EOS_Mu_ne = RegEOSMudatasize   !total rows in mu table
      double precision :: EOS_Mu_e0 = 0.0d0   !lowest energy density in mu table
      double precision :: EOS_Mu_de = 0.0d0   !spacing of energy density in mu table
      double precision :: MuEOSdata(RegEOSMudatasize, IMax_Mu)
      Integer :: IMuflag
!
!=======declear variables for reading equation of state end==================

!=======common blocks========================================================
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne

      common /EOSMudata/MuEOSdata, IMuflag
      common /EOSMudatastructure/ EOS_Mu_e0, EOS_Mu_de, EOS_Mu_ne, 
     &                            Inumparticle
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2
!=======common blocks end====================================================
      

      open(5,FILE='EOS/EOS_tables/EOS_PST.dat',STATUS='OLD')
      do I=1,EOSne
        read(5,*) ee, PEOSdata(I), SEOSdata(I), TEOSdata(I)
        if(I.eq.1) EOSe0 = ee
        if(I.eq.2) EOSde = ee - EOSe0
      enddo 
      close(5)

      EOSEend = EOSe0 + EOSde*(EOSne-1)

      open(5,FILE='EOS/EOS_tables/EOS_particletable.dat', STATUS='OLD')
      read(5,*) Inumparticle
      close(5)
      if(Inumparticle.ne.0) then
         open(5,FILE='EOS/EOS_tables/EOS_Mu.dat', STATUS='OLD')
         do I=1,EOS_Mu_ne
            read(5,*) ee, (MuEOSdata(I,J), J=1, Inumparticle)
            if(I.eq.1) EOS_Mu_e0 = ee
            if(I.eq.2) EOS_Mu_de = ee - EOS_Mu_e0
         enddo
      endif

      open(5,FILE='EOS/EOS_tables/coeff.dat', STATUS='OLD')
      read(5,*) Pcoeff1, Pcoeff2
      read(5,*) Scoeff1, Scoeff2
      read(5,*) Tcoeff1, Tcoeff2
      close(5)


      end

C====EOS from table===================================================
      Double Precision Function PEOSL7(ee)  ! for lattice P(e)
      Implicit none
!=======list of parameters===================================================
      Integer, Parameter :: RegEOSdatasize = 155500  !converted EOS table data size
!=======list of parameters end===============================================
      double precision :: ee
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: p1
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2

      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2

      ee = abs(ee)

      EOSEend = EOSe0 + EOSde*(EOSne-1)
      if(ee.lt.0.0d0) then
       write(*,*)'Warning: out of the beginning of the table, PEOSL'
       write(*,*)'current e ', ee, ' smaller than Min e, ', EOSe0
       Stop
      else if (ee.lt.EOSe0) then
            p1 = PEOSdata(1)
            PEOSL7 = ee*p1/EOSe0
      else if (ee.lt.EOSEend) then 
       call interpCubic(PEOSdata, RegEOSdatasize, 
     &                  EOSe0, EOSde, ee, PEOSL7)
      else
       PEOSL7 = Pcoeff1*(ee**Pcoeff2)
       !write(*,*) 'Warning: out of the end of table, PEOSL'
       !write(*,*) 'current e ', ee, ' larger than  Max e, ', EOSEend
       !stop
      endif

      return
      end

      Double Precision Function SEOSL7(ee)  ! for lattice P(e)
      Implicit none
!=======list of parameters===================================================
      Integer, Parameter :: RegEOSdatasize = 155500  !converted EOS table data size
!=======list of parameters end===============================================
      double precision :: ee
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: S1
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2

      ee = abs(ee)

      EOSEend = EOSe0 + EOSde*(EOSne-1)
      if(ee.lt.0.0d0) then
       write(*,*)'Warning: out of the beginning of the table, PEOSL'
       write(*,*)'current e ', ee, ' smaller than Min e, ', EOSe0
       Stop
       else if (ee.lt.EOSe0) then
            s1 = SEOSdata(1)
            SEOSL7 = ee*S1/EOSe0
      else if (ee.lt.EOSEend) then 
       call interpCubic(SEOSdata, RegEOSdatasize, 
     &                  EOSe0, EOSde, ee, SEOSL7)
      else
       SEOSL7 = Scoeff1*(ee**Scoeff2)
       !write(*,*) 'Warning: out of the end of table, PEOSL'
       !write(*,*) 'current e ', ee, ' larger than  Max e, ', EOSEend
       !stop
      endif

      return
      end

      Double Precision Function TEOSL7(ee)  ! for lattice P(e)
      Implicit none
!=======list of parameters===================================================
      Integer, Parameter :: RegEOSdatasize = 155500  !converted EOS table data size
!=======list of parameters end===============================================
      double precision :: ee
      double precision :: PEOSdata(RegEOSdatasize),
     &                    SEOSdata(RegEOSdatasize),
     &                    TEOSdata(RegEOSdatasize)
      double precision :: EOSe0         !lowest energy density
      double precision :: EOSde         !spacing of energy density
      Integer :: EOSne                  !total rows of energy density
      double precision :: EOSEend   !the maximum energy density in the table
      double precision :: T1
      double precision :: Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                    Tcoeff1, Tcoeff2
      common /EOSdata/PEOSdata, SEOSdata, TEOSdata
      common /EOSdatastructure/ EOSe0, EOSde, EOSne
      common /EOScoeffs/ Pcoeff1, Pcoeff2, Scoeff1, Scoeff2,
     &                   Tcoeff1, Tcoeff2

      ee = abs(ee)

      EOSEend = EOSe0 + EOSde*(EOSne-1)
      if(ee.lt.0.0d0) then
       write(*,*)'Warning: out of the beginning of the table, PEOSL'
       write(*,*)'current e ', ee, ' smaller than Min e, ', EOSe0
       Stop
       else if (ee.lt.EOSe0) then
            T1 = TEOSdata(1)
            TEOSL7 = ee*T1/EOSe0
      else if (ee.lt.EOSEend) then 
       call interpCubic(TEOSdata, RegEOSdatasize, 
     &                  EOSe0, EOSde, ee, TEOSL7)
      else
       TEOSL7 = Tcoeff1*(ee**Tcoeff2)
       !write(*,*) 'Warning: out of the end of table, PEOSL'
       !write(*,*) 'current e ', ee, ' larger than  Max e, ', EOSEend
       !stop
      endif

      return
      end
 
