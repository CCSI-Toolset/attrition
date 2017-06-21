!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C 
!     Module name: DES_INIT_NAMELIST                                      C
!     Purpose: DES - initialize the des-namelist                          
!                                                                         C
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!     Reviewer:                                          Date:            C
!     Modified: Wesley Xu  (add desradiusnew array)                       C
!                                                                         C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
    
      SUBROUTINE DES_INIT_ARRAYS

      USE param
      USE param1
      USE discretelement
      USE indices
      USE geometry
      USE compar
      USE physprop
      USE des_bc
      USE run
      use desgrid 
      use desmpi 
      USE des_thermo
      USE des_rxns      
 
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I

!-----------------------------------------------

      DES_RADIUS(:) = ZERO
      DESRadiusNew(:)=ZERO !ATTRITION
      PMASS(:) = ZERO
      PVOL(:) = ZERO
      OMOI(:) = ZERO
      RO_Sol(:) = ZERO 

      DES_POS_OLD(:,:) = ZERO
      DES_POS_NEW(:,:) = ZERO
      DES_VEL_OLD(:,:) = ZERO
      DES_VEL_NEW(:,:) = ZERO

      DES_VEL_OOLD(:,:) = ZERO
      DES_ACC_OLD(:,:) = ZERO

      OMEGA_OLD(:,:) = ZERO
      OMEGA_NEW(:,:) = ZERO
      ROT_ACC_OLD(:,:) = ZERO

      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      SOLID_DRAG(:,:,:) = ZERO

      FC(:,:) = ZERO
      FN(:,:) = ZERO
      FT(:,:) = ZERO
      TOW(:,:) = ZERO

      PPOS(:,:) = ZERO
      GRAV(:) = ZERO
!      DES_WALL_POS(:,:) = UNDEFINED
!      DES_WALL_VEL(:,:) = UNDEFINED

      NEIGHBOURS(:,:) = -1
      NEIGHBOURS(:,1) = 0
      PN(:,:) = -1
      PN(:,1) = 0
      PV(:,:) = 1
      PFT(:,:,:) = ZERO

      IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. &
        DES_NEIGHBOR_SEARCH .EQ. 3) THEN
          LQUAD(:,:) = UNDEFINED_I
          CQUAD(:,:) = UNDEFINED
          PQUAD(:) = 0
      ENDIF      

! pradeep desgrid related routines are moved to desgrid module
!      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN
!         DESGRIDSEARCH_PIJK(:,:) = ZERO
!      ENDIF


      PINC(:) = ZERO
      PIJK(:,:) = ZERO

! pradeep removed this geometry is used to define des grid and setting zero here violates that 
      XE(:) = ZERO
      YN(:) = ZERO
      ZT(:) = ZERO

! J.Musser: DEM inlet/outlet
      PEA(:,:) = .FALSE.
! Pradeep not proper location for the following setting 
! If RESTART_1, PEA will be read in from the restart file
!      IF(RUN_TYPE == 'NEW') THEN
!         DO I=1, PARTICLES
!            PEA(I,1)=.TRUE.
!         ENDDO
!      ENDIF
      DES_BC_U_s(:) = ZERO
      DES_BC_V_s(:) = ZERO
      DES_BC_W_s(:) = ZERO

! T.Li : Hertzian collision model
      g_mod(:) = zero
      hert_kn(:,:) = zero
      hert_kwn(:) = zero
      hert_kt(:,:) = zero
      hert_kwt(:) = zero

!pradeep for parallel processin 
      iglobal_id = 0

! cohesion VDW forces
      IF(USE_COHESION) THEN
         Fcohesive(:,:) = ZERO
         PostCohesive (:) = ZERO
      ENDIF

! J.Musser: Energy and Species Equation Arrays
      IF(DES_ENERGY_EQ)THEN
         DES_T_s_OLD(:) = UNDEFINED
         DES_T_s_NEW(:) = UNDEFINED
         DES_C_PS(:) = UNDEFINED
         DES_X_s(:,:) = ZERO
         IF(DES_CONV_EQ) Qcv(:) = ZERO
         IF(DES_COND_EQ_PP) Qpp(:) = ZERO
         IF(DES_COND_EQ_PFP) Qpfp(:) = ZERO
         IF(DES_RADI_EQ) Qrd(:) = ZERO

         IF(FIND_THERMO_NBRHD) THEN
            THERMO_NBRHD(:,:) = -1
            THERMO_NBRHD(:,1) = 0
         ENDIF
         IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') &
            Qtotal_OLD(:) = ZERO
      ENDIF

      IF(ANY_DES_SPECIES_EQ)THEN
         DES_R_sp(:,:) = ZERO
         DES_R_sc(:,:) = ZERO
         Qint(:) = ZERO
         IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') THEN
            dMdt_OLD(:) = ZERO
            dXdt_OLD(:,:) = ZERO
            dRdt_OLD(:) = ZERO
         ENDIF
      ENDIF


      RETURN
      END SUBROUTINE DES_INIT_ARRAYS 


