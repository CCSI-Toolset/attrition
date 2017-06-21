!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!                                                                    
!  Module name: DES_ALLOCATE_ARRAYS                                     
!  Purpose: allocate arrays for DES

!  Reviewer: Wesley Xu, Dave Decroix                    Date: 10-May-12
!  Revision: Allocate attrition related arrays: desradiusnew and
!  attritionflag which is currently not used so commented
!
!  Revised: Wesley Xu                                   Date: 30-July-12
!  Comments: the commented line about attritionflag is removed
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      SUBROUTINE DES_ALLOCATE_ARRAYS 
                                                                   
!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE param 
      USE param1
      USE constant
      USE discretelement
      Use indices
      Use geometry
      Use compar
      Use physprop
      Use des_bc
      use funits
      use desgrid 
      use desmpi
      USE mfix_pic
      Use des_thermo
      Use des_rxns
      
      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------      
! indices      
      INTEGER I, J, K, IJK, M 
! domain volume      
!      DOUBLE PRECISION :: VOL_DOMAIN
! the maximum and minimum specified particle diameter 
!      DOUBLE PRECISION MAX_DIAM, MIN_DIAM
! for gener_part_config, the total solids volume fraction
!      DOUBLE PRECISION TOT_VOL_FRAC
! the number of particles in the system
      INTEGER NPARTICLES
!-----------------------------------------------
      INCLUDE 'function.inc'
      
      if(dmp_log.and.debug_des) WRITE(unit_log,'(1X,A)') &
         '---------- START DES_ALLOCATE_ARRAYS ---------->'
      if(dmp_log)write(unit_log,'(3X,A,I10)') &
         'Total number of particles = ', PARTICLES      
      if(dmp_log)write(unit_log,'(3X,A,I5)') 'Dimension = ', DIMN
      NWALLS = 2*DIMN

      IF(.NOT.NON_RECT_BC) THEN
! +nwalls is included since calc_force_des temporarily uses the variables 
! pos, vel, etc at elements beyond the array index given by particles. 
! unclear why additional array space is needed via particles_factor
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS
      ELSE
! Tingwen 19/01/2008 
! +2 to include the contact with edge and node
! only one edge contact or one node contact with wall is allowed for a particle
         NPARTICLES = PARTICLES * PARTICLES_FACTOR + NWALLS + 2 + NWALLS +1
      ENDIF

! J.Musser : Dynamic Particle Info
      IF(MAX_PIS /= UNDEFINED_I .AND. MAX_PIS .GT. NPARTICLES) THEN
         IF(DMP_LOG) WRITE(unit_log, 1002) MAX_PIS, NPARTICLES, MAX_PIS
         IF(myPE.EQ.PE_IO) WRITE(*, 1002) MAX_PIS, NPARTICLES, MAX_PIS
         NPARTICLES = MAX_PIS
      ENDIF


      
      MAXQUADS = 5*PARTICLES*MQUAD_FACTOR
      IF(MAXQUADS.LE.80000) MAXQUADS = 80000
      MAXNEIGHBORS = MN + 1 + NWALLS

      IF(DIMN.EQ.3) THEN
         NMQD = 11
      ELSE
         NMQD = 7
      ENDIF

!------------------------------------------------------------------------
! pradeep 
! For parallel processing the array size required should be either 
! specified by the user or could be determined from total particles 
! with some factor  
      if (numpes.gt.1) then    
      nparticles = (nparticles/numPEs)
      end if 
      if(nparticles .lt.1000) nparticles = 1000
!max_pip adjusted to accomodate temporary variables used for walls 
      
     
      max_pip = nparticles-2*nwalls -3  


      IF(DMP_LOG) WRITE(unit_log, 1003)  NPARTICLES, MAX_PIP
      IF(mype.eq.pe_IO) WRITE(*, 1003)  NPARTICLES, MAX_PIP 

      allocate (iglobal_id(nparticles))

      !check if max_pip is greater than maximum number of 
      !computational particles for the MPPIC case 
      !if not, then stop the code after printing the error message
      IF(GENER_PART_CONFIG.AND.MPPIC) THEN 
         IF(MAX_PIP.lt. SUM(PART_MPHASE(1:MMAX))) THEN 
            
            IF(DMP_LOG) WRITE(UNIT_LOG,1001) MAX_PIP, SUM(PART_MPHASE(1:MMAX))
            CALL MFIX_EXIT(myPE)
         ENDIF
      ENDIF

 1001 FORMAT(/1X,70('*')//' From: DES_ALLOCATE_ARRAYS',/,&
         ' Message: MAX_PIP',4x,i4,4x,'is smaller than number of particles ',4x,i4,4x, /, &
         & 'INCREASE PARTICLE FACTOR',/1X,70('*')/)

 1002 FORMAT(/2X & 
      & 'USER SUPPLIED MAX_PIS = ', i10,/2X & 
      & '> NPARTICLES (cumulative size of particle arrays) = ', i10, /2X, &
      & 'THEREFORE, SETTING NPARTICLES TO ', i10)

 1003 FORMAT(/2X,  & 
      & 'PARTICLE ARRAY SIZE ON EACH PROC  = ', i10,/2X & 
      & 'MAXIMUM PHYSICAL PARTICLES (MAX_PIP) ON EACH PROC', i10, /2X, &
      & 'NOTE THAT THIS VALUE OF MAX_PIP IS ONLY RELEVANT FOR A NEW RUN ', /2X &
      & 'FOR RESTARTS, MAX_PIP WILL BE SET LATER ON')


!------------------------------------------------------------------------

! DES Allocatable arrays
!-----------------------------------------------
! J.Musser : Dynamic Particle Info
! pradeep: for parallel processing added another index (from 3 to 4)for ghost  
      ALLOCATE( PEA (NPARTICLES, 4) )

! volume of nodes 
      ALLOCATE(DES_VOL_NODE(DIMENSION_3))

! ratio of actual volume of nodes to volume of nodes not corrected for
! on the wall or being outside the domain 
      ALLOCATE(DES_VOL_NODE_RATIO(DIMENSION_3))
! T. Li: Hertzian collision model
      allocate(hert_kn(MMAX,MMAX))
      allocate(hert_kt(MMAX,MMAX))
      allocate(hert_kwn(MMAX))
      allocate(hert_kwt(MMAX)) 
      allocate(g_mod(MMAX))
      
! COEFF OF RESITUTIONS 
      ALLOCATE(REAL_EN(MMAX,MMAX)) 
      ALLOCATE(REAL_EN_WALL(MMAX))
! for hertzian model need real_et, otherwise specify eta_t_fact 
      ALLOCATE(REAL_ET(MMAX,MMAX)) 
! for hertzian model need real_et_wall, otherwise specifiy eta_t_w_fact
      ALLOCATE(REAL_ET_WALL(MMAX)) 
      ALLOCATE(DES_ETAN(MMAX,MMAX))
      ALLOCATE(DES_ETAT(MMAX,MMAX))
      ALLOCATE(DES_ETAN_WALL(MMAX), DES_ETAT_WALL(MMAX))

      
! Particle attributes
! Radius, density, mass, moment of inertia           
      Allocate(  DES_RADIUS (NPARTICLES) )
      Allocate(  RO_Sol (NPARTICLES) )
      Allocate(  PVOL (NPARTICLES) )
      Allocate(  PMASS (NPARTICLES) )
      Allocate(  OMOI (NPARTICLES) )

! Attrition new radius and flag
      Allocate( DESRadiusNew (NPARTICLES) )
  
! Old and new particle positions, velocities (translational and
! rotational)       
      Allocate(  DES_POS_OLD (NPARTICLES,DIMN) )
      Allocate(  DES_POS_NEW (NPARTICLES,DIMN) )
      Allocate(  DES_VEL_OLD (NPARTICLES,DIMN) )
      Allocate(  DES_VEL_NEW (NPARTICLES,DIMN) )
      Allocate(  DES_VEL_OOLD(NPARTICLES,DIMN) )
      Allocate(  DES_ACC_OLD (NPARTICLES,DIMN) )

      IF(DIMN.GT.2) THEN
         Allocate(  OMEGA_OLD (NPARTICLES,DIMN) )
         Allocate(  OMEGA_NEW (NPARTICLES,DIMN) )
         ALLOCATE(  ROT_ACC_OLD (NPARTICLES,DIMN))
      ELSE
         Allocate(  OMEGA_OLD (NPARTICLES,1) )
         Allocate(  OMEGA_NEW (NPARTICLES,1) )
         ALLOCATE(  ROT_ACC_OLD (NPARTICLES,1))
      ENDIF        
     
! Total, normal and tangetial forces      
      Allocate(  FC (NPARTICLES,DIMN) )
      Allocate(  FN (NPARTICLES,DIMN) )
      Allocate(  FT (NPARTICLES,DIMN) )
      Allocate(  FTAN (DIMN) )
      Allocate(  FNORM (DIMN) )
      Allocate(  GRAV (DIMN) )

! Torque     
      IF(DIMN.EQ.3) THEN 
         Allocate(  TOW (NPARTICLES,DIMN) )
      ELSE
         Allocate(  TOW (NPARTICLES,1) )
      ENDIF
     
!      IF(.NOT.MPPIC) THEN 
         !particle positions at the last call neighbor search algorithm call 
         Allocate(  PPOS (NPARTICLES,DIMN) )

! Accumulated spring force      
         Allocate(  PFT (NPARTICLES,MAXNEIGHBORS,DIMN) )
! added by Tingwen to save the normal direction at previous time step
	 Allocate(  PFN (NPARTICLES,MAXNEIGHBORS,DIMN) )
! Tracking variables for particle contact history
         Allocate(  PN (NPARTICLES, MAXNEIGHBORS) )
         Allocate(  PV (NPARTICLES, MAXNEIGHBORS) )
         
         ! Neighbor search
         Allocate(  NEIGHBOURS (NPARTICLES, MAXNEIGHBORS) )

         IF (DES_NEIGHBOR_SEARCH .EQ. 2 .OR. DES_NEIGHBOR_SEARCH .EQ. 3) THEN
            Allocate(  LQUAD (MAXQUADS, NMQD) )
            Allocate(  PQUAD (NPARTICLES) )
            Allocate(  CQUAD (MAXQUADS, NWALLS) )
         ENDIF

!      ENDIF
    
! Temporary variables to store wall position, velocity and normal vector
!      Allocate(  DES_WALL_POS (NWALLS,DIMN) )
!      Allocate(  DES_WALL_VEL (NWALLS,DIMN) )
      Allocate(  WALL_NORMAL  (NWALLS,DIMN) )

! Variable that stores the particle in cell information (ID) on the
! computational grid defined by imax, jmax and kmax in mfix.dat
! pradeep- modified from 3d to 1d array 
      ALLOCATE(PIC(DIMENSION_3))
      DO IJK=1,DIMENSION_3
        NULLIFY(pic(ijk)%p) 
      ENDDO 

! Particles in a computational cell (for volume fraction)
      Allocate(  PINC (DIMENSION_3) )
! For each particle track its i,j,k location on computational grid
! defined by imax, jmax and kmax in mfix.dat and phase no.  
! plus to mark if the particle is close to wall by Tingwen
!      Allocate(  PIJK (NPARTICLES,6) )
      Allocate(  PIJK (NPARTICLES,5) )

! pradeep allocation of desgrid is moved to module desgrid
! For each particle track its i, j, k index according to the grid
! based search mesh when des_neighbor_search=4 
!      IF (DES_NEIGHBOR_SEARCH .EQ. 4) THEN      
!         ALLOCATE( DESGRIDSEARCH_PIJK (NPARTICLES,3) )
!      ENDIF

! Pradeep modifying the three dimensional arrays to one dimensional arrays
      IF(.NOT.MPPIC) THEN 
         IF(DES_INTERP_ON) THEN
            ALLOCATE(DRAG_AM(DIMENSION_3, MMAX))
            ALLOCATE(DRAG_BM(DIMENSION_3, DIMN, MMAX))
            ALLOCATE(WTBAR(DIMENSION_3,  MMAX))
            ALLOCATE(VEL_FP(NPARTICLES,3))
            ALLOCATE(F_gp(NPARTICLES ))  
            F_gp(1:NPARTICLES)  = ZERO
         ENDIF
      ELSE
         ALLOCATE(DRAG_AM(DIMENSION_3, MMAX))
         ALLOCATE(DRAG_BM(DIMENSION_3, DIMN, MMAX))
         ALLOCATE(WTBAR(DIMENSION_3,  MMAX))
         ALLOCATE(VEL_FP(NPARTICLES,3))
         ALLOCATE(F_gp(NPARTICLES ))  
         F_gp(1:NPARTICLES)  = ZERO
      ENDIF

! Drag exerted by the gas on solids
      Allocate(  SOLID_DRAG (DIMENSION_3, DIMENSION_M, DIMN) )
     

      ALLOCATE(MARK_PART(NPARTICLES))
      ALLOCATE(BED_HEIGHT(MMAX))

! Volume averaged solids volume in a cell      
      Allocate(  DES_U_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_V_s (DIMENSION_3, DIMENSION_M) )
      Allocate(  DES_W_s (DIMENSION_3, DIMENSION_M) )

! Averaged velocity obtained by avraging over all the particles
      ALLOCATE(  DES_VEL_AVG(DIMN) )

! Granular temperature
      Allocate(  DES_THETA (DIMENSION_3, DIMENSION_M) )

! Global Granular Energy
      ALLOCATE(  GLOBAL_GRAN_ENERGY(DIMN) )
      ALLOCATE(  GLOBAL_GRAN_TEMP(DIMN) )
    
! Cell faces
! Pradeep added 0 to store IMIN3 values 
      Allocate(  XE (0:DIMENSION_I) )
      Allocate(  YN (0:DIMENSION_J) )
      Allocate(  ZT (0:DIMENSION_K) )
! force due to gas-pressure gradient 
      ALLOCATE(P_FORCE(DIMENSION_3,DIMN))
! MP-PIC related 
      IF(MPPIC) THEN 
         Allocate( PS_FORCE_PIC(DIMENSION_3, DIMN))
         
         IF(.NOT.ALLOCATED(F_gp)) THEN 
            ALLOCATE(F_gp(NPARTICLES ))  
            F_gp(1:NPARTICLES)  = ZERO
            !allocate F_gp  for the mppic case. This will be used for the implicit treatment of drag force 
         ENDIF

         ALLOCATE(DES_STAT_WT(NPARTICLES))
         
         ALLOCATE(DES_VEL_MAX(DIMN))

         ALLOCATE(PS_GRAD(NPARTICLES, DIMN))
         
         ALLOCATE(AVGSOLVEL_P(NPARTICLES, DIMN))

         ALLOCATE(EPG_P(NPARTICLES))

      ENDIF 
     
! COHESION      
! Square-well potential parameters
      IF(USE_COHESION) THEN 
         IF(SQUARE_WELL)THEN ! these are used only with square-well method
           Allocate(  WELL_WIDTH (NPARTICLES) )
           Allocate(  WELL_DEPTH (NPARTICLES) )
! Matrix location of particle 
           Allocate(  PART_GRID (NPARTICLES,4) )
         ENDIF
         Allocate(  FCohesive (NPARTICLES,DIMN) )
         Allocate(  PostCohesive (NPARTICLES) )
! Does particle have at least one linked partner
!         Allocate(  IS_LINKED (NPARTICLES) )       ! array not used
! Does particle have at least one aggloerated partner
!         Allocate(  IS_AGGLOMERATED (NPARTICLES) )       ! array not used
! Array of linked partners
!         Allocate(  LINKS (NPARTICLES, MAXNEIGHBORS) )       ! array not used
! Array of agglomerated partners
!         Allocate(  AGGS (NPARTICLES, MAXNEIGHBORS) )       ! array not used
      ENDIF

! BEGIN Thermodynamic Allocation ---------------------------------------
      IF(DES_ENERGY_EQ)THEN

! Particle temperature
         Allocate( DES_T_s_OLD( NPARTICLES ) )
         Allocate( DES_T_s_NEW( NPARTICLES ) )
! Specific heat
         Allocate( DES_C_PS( NPARTICLES ) )
! Species mass fractions comprising a particle. This array may not be
! needed for all thermo problems.
         Allocate( DES_X_s( NPARTICLES, MAX_DES_NMAX) )
! Convection Specific arrays
         IF(DES_CONV_EQ) &
            Allocate( Qcv( NPARTICLES ) )
! Conduction Specific arrays
         IF(DES_COND_EQ_PP) &
            Allocate( Qpp( NPARTICLES ) )
! Conduction Specific arrays
         IF(DES_COND_EQ_PFP) &
            Allocate( Qpfp( NPARTICLES ) )
! Radiation Specific arrays
         IF(DES_RADI_EQ) &
            Allocate( Qrd( NPARTICLES ) )
! Stores number of neighbors based on neighbor search
         IF(FIND_THERMO_NBRHD) &
            Allocate( THERMO_NBRHD( NPARTICLES, MAXNEIGHBORS ) )
! Allocate the history variables for Adams-Bashforth integration
         IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') &
            Allocate( Qtotal_OLD( NPARTICLES ) )
      ENDIF
!------------------------------------------ End Thermodynamic Allocation


! BEGIN Species Allocation ---------------------------------------------
      IF(ANY_DES_SPECIES_EQ)THEN
! Rate of solids phase production for each species
         Allocate( DES_R_sp( NPARTICLES, MAX_DES_NMAX) )
! Rate of solids phase consumption for each species
         Allocate( DES_R_sc( NPARTICLES, MAX_DES_NMAX) )
! Allocate the history variables for Adams-Bashforth integration
         IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') THEN
! Rate of chnage of particle mass
            Allocate( dMdt_OLD( NPARTICLES ) )
! Rate of chnage of particle mass percent species
            Allocate( dXdt_OLD( NPARTICLES, MAX_DES_NMAX) )
         ENDIF

! Energy generation from reaction (cal/sec)
         Allocate( Qint( NPARTICLES ) )
         IF( TRIM(REACTION_MODEL) == 'SHRINKING_CORE')THEN
! Radius of unreacted core
            Allocate( CORE_RAD( NPARTICLES ) )
! Density of unreacted core
            Allocate( CORE_RHO( NPARTICLES ) )
! Allocate the history variables for Adams-Bashforth integration
            IF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') &
! Rate of chnage of radius of unreacted core
               Allocate( dRdt_OLD( NPARTICLES ) )
         ENDIF
      ENDIF
!------------------------------------------------ End Species Allocation

      
      if(dmp_log.and.debug_des)write(unit_log,'(1X,A)')&
         '<---------- END DES_ALLOCATE_ARRAYS ----------'

      RETURN
      END SUBROUTINE DES_ALLOCATE_ARRAYS 

      
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: ALLOCATE_DES_MIO                                       !
!                                                                      !
!  Purpose:                                                            !
!                                                                      !
!  Author: J.Musser                                   Date: 17-Aug-09  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

      SUBROUTINE ALLOCATE_DES_MIO

      USE des_bc
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local variables
!-----------------------------------------------
      INTEGER I     ! Loop counter for no. of DES_BCMI

!-----------------------------------------------

! Allocate/Initialize for inlets      
      IF(DES_BCMI /= 0)THEN

! Boundary condition ID array
         Allocate( DES_BC_MI_ID (DES_BCMI) )

! Distance offset of incoming particles in ghost cell
         Allocate( DES_BC_OFFSET (DES_BCMI) )

! Particle injection factor
         Allocate( PI_FACTOR (DES_BCMI) )

! Particle injection count (injection number)
         Allocate( PI_COUNT (DES_BCMI) )

! Particle injection time scale
         Allocate( DES_MI_TIME (DES_BCMI) )

! Boundary classification
         Allocate( DES_MI_CLASS (DES_BCMI) )
         Allocate( PARTICLE_PLCMNT (DES_BCMI) )

! Order inlet condition variables
! (only needed if particle_plcmt is assigned 'ordr')
         Allocate( MI_FACTOR (DES_BCMI) )
         Allocate( MI_WINDOW (DES_BCMI) )
         Allocate( MI_ORDER (DES_BCMI) )   ! type dmi
         Allocate( I_OF_MI ( DES_BCMI) )   ! type dmi
         Allocate( J_OF_MI ( DES_BCMI) )   ! type dmi

! Grid search loop counter array; 6 = no. of faces
         Allocate(  GS_ARRAY (DES_BCMI, 6) )

! Logical array stating if a bounday condition is polydisperse
         Allocate( DES_BC_POLY( DES_BCMI ) )

! Array used for polydisperse inlets: stores the particle number
! distribution of an inlet scaled with numfrac_limit
         Allocate( DES_BC_POLY_LAYOUT( DES_BCMI, NUMFRAC_LIMIT ) )

! Initializiation
! Logical for whether inlet is polydisperse         
         DES_BC_POLY(:) = .FALSE.
! Logical for inlet existance on IJK face         
         DES_MI_X = .FALSE.
         DES_MI_Y = .FALSE.
         DES_MI_Z = .FALSE.          
! Integer arrays
         DES_BC_MI_ID(:) = -1
         PI_FACTOR(:) = -1
         PI_COUNT(:) = -1
         MI_FACTOR(:) = -1
         MI_WINDOW(:) = -1
         GS_ARRAY(:,:) = -1
         DES_BC_POLY_LAYOUT(:,:) = -1
! Double precision arrays
         DES_MI_TIME(:) = UNDEFINED
! Character precision arrays
         DES_MI_CLASS(:) = UNDEFINED_C
         PARTICLE_PLCMNT(:) = UNDEFINED_C
! Derived data types
         DO I = 1,DES_BCMI
            NULLIFY( MI_ORDER(I)%VALUE )
            NULLIFY( I_OF_MI(I)%VALUE )
            NULLIFY( J_OF_MI(I)%VALUE )
         ENDDO

      ENDIF  ! end if des_bcmi /= 0


! Allocate/Initialize for outlets
      IF(DES_BCMO /= 0)THEN
 
! Boundary Condition ID array
         Allocate( DES_BC_MO_ID (DES_BCMO) )

! Boundary classification
         Allocate( DES_MO_CLASS (DES_BCMO) )

! Initializiation
! Integer arrays         
         DES_BC_MO_ID(:) = -1
! Character arrays         
         DES_MO_CLASS(:) = UNDEFINED_C         
! Logical for outlet existance on IJK face         
         DES_MO_X = .FALSE.
         DES_MO_Y = .FALSE.
         DES_MO_Z = .FALSE.         

      ENDIF   ! end if des_bcmo /= 0


      RETURN
      END SUBROUTINE ALLOCATE_DES_MIO

