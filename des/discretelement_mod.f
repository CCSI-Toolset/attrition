!                                                                       C
!   Module name: DISCRETELEMENT                                         C
!   Purpose: DES mod file 
!            Common Block containing DEM conditions 

!   PARALLEL PROCESSING: explanation of variables in parallel architecture  
!            particles - total number of particles in the whole system 
!            read from input or generated             
!            pip - particles in each processor (includes the ghost particles)
!            max_pis - Maximum allocated Particles in system  (input parameter) 
!            max_pip - maximum allocated particles in processor 
!  
!   Reviewer: Dave Decroix, Wesley Xu                     Date: 10-May-12
!   Revision: Added statements for attrition model
!
!   Revised: Wesley Xu                                    Date: 30-July-12
!   Comments: Added definition of the ABRAlpha for abrasive attrition
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   
      MODULE DISCRETELEMENT

      USE param
      USE param1

! Define interface - needed when passing arrays of assumed size      
      INTERFACE
         INTEGER FUNCTION DES_GETINDEXFROMPOS(LIM1,LIM2,PART_POS,&
            GRID_POS,AXIS,AXIS_INDEX)
            INTEGER :: LIM1, LIM2
            DOUBLE PRECISION :: PART_POS            
            DOUBLE PRECISION, DIMENSION(:) :: GRID_POS
            CHARACTER*1 :: AXIS,AXIS_INDEX
         END FUNCTION DES_GETINDEXFROMPOS
      END INTERFACE

! Total number of particles in simulation 
      INTEGER PARTICLES

! For parallel processing: global id of particles 
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGLOBAL_ID

! ghost cnt of particles on each processor
      INTEGER :: IGHOST_CNT
! Maximum global id, new particles global id will be assigned based on this value 
      Integer :: imax_global_id

!volume of each node. Used to obtain Eulerian fields
      double precision, allocatable, dimension(:) :: des_vol_node

!Ratio of actual volume of each node to volume of node not corrected for
!outside the domain or being in cut-cell 
      double precision, allocatable, dimension(:) :: des_vol_node_ratio
 

! Constant factor used to expand size of arrays beyond particle no.      
      DOUBLE PRECISION PARTICLES_FACTOR 

! Logic that controls whether to print data dem simulations (granular or
! coupled)
      LOGICAL PRINT_DES_DATA 

! logic that controls if des run time messages are printed on screen or not
      LOGICAL PRINT_DES_SCREEN  

! Usr specified time interval that controls frequency of writing DEM
! output and restart for pure granular flow; otherwise (when coupled)
! the frequency of writing output and restart is controlled by the
! value of spx_dt(1) and res_dt
      DOUBLE PRECISION DES_SPX_DT, DES_RES_DT

! This specifies the file type used for outputting DES data
! options are :      
!    TECPLOT - data is written in Tecplot format
!    undefined - data is written in ParaView format (default)
      CHARACTER*64 DES_OUTPUT_TYPE

! Used sporadically to control screen dumps (for debug purposes)
      LOGICAL :: DEBUG_DES

! Single particle no. index that is followed if debugging
      INTEGER FOCUS_PARTICLE

! Output file count for .vtp type files and for tecplot files;
! for vtp output used to label .vtp file names in sequential order
! and is saved so restarts begin at the correct count
      INTEGER VTP_FINDEX, TECPLOT_FINDEX

! DES - Continuum       
      LOGICAL DISCRETE_ELEMENT 
      LOGICAL DES_CONTINUUM_COUPLED
     
! DES - Invoke hybrid model where both the DEM and continuum model
! are employed to describe solids 
      LOGICAL DES_CONTINUUM_HYBRID 

! DES -
      LOGICAL DES_ONEWAY_COUPLED
!this logic is different from des_continuum coupled. Here the particles 
! see the fluid but the fluid does not see the particles 

! Only used when coupled and represents the number of times a pure
! granular flow simulation is run before the actual coupled simulation
! is started (i.e. for particle settling w/o fluid forces)
      INTEGER NFACTOR

! Switch to decide whether to call drag_gs or to call des_drag_gs via
! drag_fgs to calculate the gas-solids drag coefficient.  if false
! then drag_gs is used, otherwise des_drag_gs via drag_fgs is used
      LOGICAL DES_INTERP_ON

! Drag      
      LOGICAL TSUJI_DRAG

! Collision model, options are as follows
!   linear spring dashpot model (default/undefined)
!   'hertzian' model
      CHARACTER*64 DES_COLL_MODEL

! Integration method, options are as follows
!   'euler' first-order scheme (default)
!   'adams_bashforth' second-order scheme (by T.Li)
      CHARACTER*64 DES_INTG_METHOD 


! Value of solids time step based on particle properties
      DOUBLE PRECISION DTSOLID
! Currently obsolete quantity
      DOUBLE PRECISION DTSOLID_FACTOR 
! Run time value of simulation time used in dem simulation
      DOUBLE PRECISION S_TIME
   
! Neighbor search method, options are as follows
!   1= nsquare, 2=quadtree, 3=octree, 4=grid based search      
      INTEGER DES_NEIGHBOR_SEARCH
! Quantities used for reporting: max no. neighbors and max overlap
! that exists during last solid time step of dem simulation
      INTEGER NEIGH_MAX
      DOUBLE PRECISION OVERLAP_MAX
! Quantities used for neighbor search methods octree and quadtree
      INTEGER QLM, QLN, INIT_QUAD_COUNT, INQC
      INTEGER NQUAD, MAXQUADS, NMQD 
      DOUBLE PRECISION N2CT, QUADCT, OCTCT, MQUAD_FACTOR
      DOUBLE PRECISION RADIUS_EQ

! Related quantities used to determine whether neighbor search 
! should be called       
      INTEGER NEIGHBOR_SEARCH_N
      DOUBLE PRECISION NEIGHBOR_SEARCH_RAD_RATIO
      LOGICAL DO_NSEARCH


! Factor muliplied by sum of radii in grid based neighbor search and
! nsquare search method.  increases the effective radius of a particle
! for detecting particle contacts 
      DOUBLE PRECISION FACTOR_RLM

! User specified value for the max number of neighbors any particle is
! allowed
      INTEGER MN

! Slightly larger than MN: used to specify the max number of neighbors 
! any particle is allowed plus the number of walls in the simulation
! plus one
      INTEGER MAXNEIGHBORS

! Stores number of neighbors based on neighbor search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: NEIGHBOURS  !(PARTICLES, MAXNEIGHBORS)


! Neighbor search quantities for quad based search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LQUAD          !(MAXQUADS, NMQD)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PQUAD            !(PARTICLES) 
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CQUAD !(NWALLS,MAXQUADS)
      
! User specified dimension of the system (by default 2D, but if 3D system is
! desired then it must be explicitly specified)      
      INTEGER DIMN

! Variable that is set to the number of walls in the system (=2*DIMN)
      INTEGER NWALLS

! Position of domain boundaries generally given as 
!   (0, xlength, 0, ylength, 0, zlength)
      DOUBLE PRECISION WX1, EX2, BY1, TY2, SZ1, NZ2

! Gravity vector
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GRAV !(DIMN)

! X, Y, Z position of cell faces
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: XE  !(DIMENSION_I)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: YN  !(DIMENSION_J)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ZT  !(DIMENSION_K)

! Wall position, velocity and normal vector (used to temporarily store
! wall position and velocity)
!      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_POS !(NWALLS,DIMN)
!      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_WALL_VEL !(NWALLS,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: WALL_NORMAL  !(NWALLS,DIMN)

! Wall vibration parameters
      DOUBLE PRECISION DES_GAMMA, DES_F
!     if a value for des_f is given then this is used to move the 
!     bottom/top walls (xz plane) east/west      
      DOUBLE PRECISION LID_VEL 

! Particle treatment at the walls  
!     walldtsplit should be set to T for particle interaction with walls
      LOGICAL WALLDTSPLIT
      LOGICAL WALLREFLECT

! Currently obsolete variable
      LOGICAL NON_RECT_BC   
! Variable for cell_near_wall (obsolete)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_near_w

! Periodic wall BC
      LOGICAL DES_PERIODIC_WALLS
      LOGICAL DES_PERIODIC_WALLS_X
      LOGICAL DES_PERIODIC_WALLS_Y
      LOGICAL DES_PERIODIC_WALLS_Z

! Lees & Edwards wall BC
! Logic for Lees & Edwards BC (T = turn on LE BC)
      LOGICAL DES_LE_BC
! Relative velocity of LE boundaries (distance/time)      
      DOUBLE PRECISION DES_LE_REL_VEL
! Shear direction
!   2D options are DUDY or DVDX 
!   3D options are DUDY, DUDZ, DVDX, DVDZ, DWDX or DWDY
!   Note that all other directions are treated as periodic boundaries
      CHARACTER*4 DES_LE_SHEAR_DIR


! If gener_part_config is true, then the particle_input.dat file does
! not need to be supplied nor does the total number of particles as
! these are determined based on the specified volume fraction (vol_frac)
! in the specified domain (des_eps_xyzstart)     
      LOGICAL :: GENER_PART_CONFIG
      DOUBLE PRECISION ::  VOL_FRAC(DIM_M), DES_EPS_XSTART, &
                           DES_EPS_YSTART, DES_EPS_ZSTART
! The number of particles that belong to solid phase M according to the
! vol_frac and particle diameter. this information is used when
! gener_part_config is invoked for initialization 
      INTEGER PART_MPHASE(DIM_M)

! Assigns the initial particle velocity distribution based on user
! specified mean and standard deviation (regardless if already set
! within particle_input.dat)
      DOUBLE PRECISION pvel_mean, PVEL_StDev
! Constant input pressure gradient (currently unused?)
      DOUBLE PRECISION pgrad(3)

! Particle properties
! particle diameters
      DOUBLE PRECISION DES_D_P0 (DIM_M)
! particle densities      
      DOUBLE PRECISION DES_RO_s (DIM_M)
! number of solids phases
      INTEGER DES_MMAX
            
! Particle-particle and Particle-wall contact parameters
!     Spring contants      
      DOUBLE PRECISION KN, KN_W  !Normal
      DOUBLE PRECISION KT, KT_W, KT_FAC, KT_W_FAC  !Tangential factors = KT/KN and KT_w/KN_w, resp.
!     Damping coeffients      
      DOUBLE PRECISION ETA_DES_N, ETA_N_W  !Normal
      DOUBLE PRECISION ETA_DES_T, ETA_T_W  !Tangential
!     Tangential damping factors, eta_t = eta_t_factor * eta_N
      DOUBLE PRECISION DES_ETAT_FAC, DES_ETAT_W_FAC
!     Damping coeffients in array form 
      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE :: DES_ETAN, DES_ETAT   !(DES_MMAX,DES_MMAX)
      DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE :: DES_ETAN_WALL, DES_ETAT_WALL   !(DES_MMAX)
!     Friction coeficients
      DOUBLE PRECISION MEW, MEW_W
!     coeff of restituion input in one D array, solid solid
!     Tangential rest. coef. are used for hertzian collision model but not linear
      DOUBLE PRECISION DES_EN_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
      DOUBLE PRECISION DES_ET_INPUT(DIM_M+DIM_M*(DIM_M-1)/2)
!     coeff of restituion input in one D array, solid wall 
      DOUBLE PRECISION DES_EN_WALL_INPUT(DIM_M) 
      DOUBLE PRECISION DES_ET_WALL_INPUT(DIM_M)
!     actual coeff of rest.'s rearranged 
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::  REAL_EN, REAL_ET   !(DES_MMAX,DES_MMAX)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE ::  REAL_EN_WALL, REAL_ET_WALL   !(DES_MMAX)

! Hertzian model: T.Li
      DOUBLE PRECISION ew_young, vw_poisson
      DOUBLE PRECISION e_young(dim_m), v_poisson(dim_m)    
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hert_kn, hert_kt  !(DES_MMAX,DES_MMAX)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: hert_kwn, hert_kwt  !(DES_MMAX)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: g_mod               !(DES_MMAX)

! Particle attributes: radius, density, mass, moment of inertia      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_RADIUS !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RO_Sol     !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PVOL       !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PMASS      !(PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: OMOI       !(PARTICLES)
     
! Old and new particle positions, velocities (translational and
! rotational)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_OLD  !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_POS_NEW  !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_OLD  !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_NEW  !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_OLD    !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: OMEGA_NEW    !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPOS         !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_VEL_OOLD !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_ACC_OLD  !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ROT_ACC_OLD  !(PARTICLES,DIMN)

! Total, normal and tangetial forces      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FC    !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FN    !(PARTICLES,DIMN)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FT    !(PARTICLES,DIMN)

! There is no need to maintain FN and FT arrays. This can be accomplished by FTAN and FNORM. 
      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: FTAN, FNORM ! (DIMN)

! Drag exerted by the gas on solids
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: SOLID_DRAG ! (DIMENSION_3, DIMENSION_M, DIMN) 

! Torque      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TOW   !(PARTICLES,DIMN)
     
! Accumulated spring force
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFT !(PARTICLES,DIMN,MAXNEIGHBORS)
!     added by Tingwen to save the normal direction at previous time step
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: PFN ! (PARTICLES,DIMN,MAXNEIGHBORS)

! Variables used to track/store particle contact history      
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PN !(PARTICLES, MAXNEIGHBORS)
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PV !(PARTICLES, MAXNEIGHBORS)
 
! Run time logic.  if calc_fc is true, then the contact forces (FC) are 
! updated to include gas-solids drag and gas pressure in the call to 
! drag_fgs. calc_fc does not play a role in pure granular flow sims
      LOGICAL CALC_FC

! Run time logic. if callfromdes is true, then the pertinent mean fields
! (in this case ROP_S and F_GS) are not computed/updated in the call to
! drag_fgs.  it is done to speed up the simulation. callfromdes does
! not play a role in pure granular flow sims and is only relevant when
! des_interp_on is set to T
      LOGICAL CALLFROMDES


! Dynamic variable; for each cell stores the total number of particles
! and the id's of the particles in that cell
      TYPE iap1
         INTEGER, DIMENSION(:), POINTER:: p
      END TYPE iap1
!pradeep
! in order to facilitate the paralle processing the PIC is now defined as single array IJK
      TYPE(iap1), DIMENSION(:), ALLOCATABLE:: pic  ! (DIMENSION_3)

! Particles in a computational cell (for volume fraction)
      INTEGER, DIMENSION(:), ALLOCATABLE :: PINC  ! (DIMENSION_3)

! For each particle track its i,j,k location on grid and phase no.:      
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PIJK ! (PARTICLES,5)=>I,J,K,IJK,M 

! pradeep all variables related to desgrid is moved to module desgrid
! except the input variables

!!********************************************************************************
!! START when des_neighbor_search = 4 
!
!! Dynamic variable; for each cell stores the total number of particles
!! and the id's of the particles in that cell. This variable is used
!! exclusively for the des_neighbor_search option 4 (cell based search
!! method) and was created so that the mesh used to perform the cell
!! based search is independent of the standard mesh used to solve the
!! fluid equations in coupled dem cases
!! pradeep - for parallel processing
!! modifying from 3D array to 1D array
      TYPE(iap1), DIMENSION(:,:,:), ALLOCATABLE:: DESGRIDSEARCH_PIC  
!                                                       ! (desgs_imax2,desgs_jmax2,desgs_kmax2)
!
!! The number of i, j, k divisions in the grid used to perform the
!! cell based neighbor search
      INTEGER :: DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, &
                 DESGRIDSEARCH_KMAX
!! The number of i, j, k divisions in the grid used to perform the
!! cell based neighbor search plus 2 to account for ghost cells
      INTEGER :: DESGS_IMAX2, DESGS_JMAX2, &
                 DESGS_KMAX2
!
!! For each particle track its i,j,k location on the mesh used for the
!! cell based neighbor search option
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: DESGRIDSEARCH_PIJK  ! (PARTICLES,3)=>I,J,K
!
!! END when des_neighbor_search = 4       
!!********************************************************************************


! Volume averaged solids velocity in a fluid cell      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_U_s  
                        !(DIMENSION_3,DES_MMAX)
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_V_s
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_W_s

! Granular temperature in a fluid cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DES_THETA 
                        !(DIMENSION_3,DES_MMAX)

! Global average velocity: obtained by averaging over all the particles
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DES_VEL_AVG 
                        !(DIMN)
     
! Global granular energy & temp: obtained by averaging over all the particles
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GLOBAL_GRAN_ENERGY 
                        !(DIMN)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GLOBAL_GRAN_TEMP   
                        !(DIMN)

! Kinetic and potential energy of the system: obtained by averaging
! over all particles
      DOUBLE PRECISION DES_KE, DES_PE



! Additional quantities
      DOUBLE PRECISION :: MIN_RADIUS, MAX_RADIUS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MARK_PART

! Logic for bed height calculations (T = turn on bed height
! calculations)
      LOGICAL DES_CALC_BEDHEIGHT      
! Used to track bed height of solids phase M 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bed_height 

! START interpolation related data
! R.Garg 

! the coefficient add to gas momentum A matrix  at cell corners
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::drag_am 

! the coefficient add to gas momentum B matrix  at cell corners
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE ::drag_bm 

! fluid velocity at particle position
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::vel_fp 

      
      DOUBLE PRECISION, DIMENSION(:,:,:),POINTER :: weightp    
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_gp 
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: wtderivp
      double precision, dimension (:,:), allocatable ::  wtbar
      DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: sstencil
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE:: gstencil, vstencil, pgradstencil

!stencil for interpolation of solids pressure 
      DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE::  psgradstencil
!stencil for interpolation of solids velocity
      DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE::  VEL_SOL_STENCIL
      
! quantities are set in subroutine set_interpolation_scheme
! order = order of the interpolation method, ob2l = (order+1)/2,
! ob2r = order/2      
      CHARACTER(LEN=7):: scheme, interp_scheme
      INTEGER:: order, ob2l, ob2r
      
! END of interpolation related data


! START particle inlet/outlet related quantites 
! J.Musser

! Dynamic particle count elements:
! PEA(n,1) : This column identifies particle as 'existing' if true. 
! It is used with the inlet/outlet to skip indices that do not represent particles
! in the system or indices that represent particles that have exited the system.
! PEA(n,2) : This column identifies a particle as 'new' if true.
! Particles with a classification of 'new' do not react when in contact with a wall 
! or another particle, however existing particles do collide and interact with
! 'new' particles. The classification allows new particles to push particles 
! already in the system out of the way when entering to prevent overlap.
! PEA(n,3) : This column identifies a particle as 'exiting' if true. 
! If a particle initiates contact with a wall surface designated as a des outlet,
! this flag is set to true. With this classification the location of the particle 
! is checked to assess if the particle has fully exited the system.  At this point,
! the particle is removed from the list.

! pradeep: setting flags for ghost cell 
! PEA(n,4) : for ghost particles 
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: PEA! (MAX_PIS,3)

! Number of particles in the system (current)
      INTEGER PIP

! Maximum particles permitted in the system at once
      INTEGER MAX_PIS, MAX_PIP

! END particle inlet/outlet related quantities


!     Force due to gas-pressure gradient 
      DOUBLE PRECISION , DIMENSION(:,:), ALLOCATABLE :: P_FORCE
 
     


      INTEGER :: MAX_DES_BC_CELL 
      TYPE BDRY_TYPE
         INTEGER :: IJK_SCAL
         CHARACTER*100 :: DES_BC_TYPE
         DOUBLE PRECISION, DIMENSION (:) , ALLOCATABLE :: NORMAL 
         DOUBLE PRECISION, DIMENSION (:) , ALLOCATABLE :: CNOT, VEC
         !MI_BCID will contain the BC id for the mass inflow bc. It is made allocatable 
         !to reduce memory usage
         INTEGER, DIMENSION (:) , ALLOCATABLE :: MI_BCID
      END TYPE BDRY_TYPE 
      
      TYPE DES_BDR
      INTEGER :: COUNT_DES_BC

      TYPE(BDRY_TYPE), DIMENSION(:), ALLOCATABLE :: BDRY_LIST 
      END TYPE DES_BDR

      TYPE(DES_BDR), DIMENSION(:), ALLOCATABLE :: DES_CELLWISE_BCDATA 


      LOGICAL, DIMENSION(:), ALLOCATABLE :: TOBE_DELETED
      
!     MAX Velocity of particles in each direction
      DOUBLE PRECISION , DIMENSION(:), ALLOCATABLE :: DES_VEL_MAX




! Start Cohesion
! M.Weber      
     
!     Square-well potential parameters
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WELL_WIDTH ! (PARTICLES)
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: WELL_DEPTH ! (PARTICLES)
!     Store cohesive forces
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: FCohesive ! particles, dimn
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PostCohesive ! particles, dimn
      
!     Well depth read in from mfix.dat
      DOUBLE PRECISION MASTER_WELL_DEPTH
      DOUBLE PRECISION MASTER_WALL_WELL_DEPTH

!     Ratio of well width to inner radius
      DOUBLE PRECISION RADIUS_RATIO
      DOUBLE PRECISION WALL_RADIUS_RATIO

!     Array of linked partners
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LINKS ! (PARTICLES, MAXNEIGHBORS)

!     Array of agglomerated partners
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: AGGS ! (PARTICLES, MAXNEIGHBORS)

!     Does particle have at least one linked partner
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_LINKED ! (PARTICLES)

!     Switch to turn cohesion on and off
      LOGICAL USE_COHESION      

!     Switch to turn square well on and off
      LOGICAL SQUARE_WELL

!     Switch to run debuggin on and off
      INTEGER COHESION_DEBUG
      INTEGER COHESION_DEBUG_START

!     Specific particle to target in debugging
      INTEGER COH_DEBUG_PARTICLE

!     Specific time step to target in debugging
      INTEGER COH_DEBUG_STEP

!     1x2 array to hold the last two interacing particles
      INTEGER LAST_COLLISION(2)

!     Parameters to control writing of animation file
      INTEGER ANIMATION_WRITE_INTERVAL
      INTEGER ANIMATION_WRITE_COUNTER
      DOUBLE PRECISION TIME_FACTOR
      INTEGER MAX_ANIMATOR_STEP
      INTEGER MIN_ANIMATOR_STEP

!     Parameters to control writing of single particle log and collision log
      INTEGER MAX_LOG_STEP
      INTEGER MIN_LOG_STEP

!     Does particle have at least one aggloerated partner
      INTEGER, DIMENSION(:), ALLOCATABLE :: IS_AGGLOMERATED ! (PARTICLES)


!     Number of search grids in the x-direction
      INTEGER SEARCH_GRIDS(3)

!     Matrix of particles in each search grid
      INTEGER PART_IN_GRID(0:110, 0:110, 0:110, 0:55)

!     Matrix location of particle 
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: PART_GRID ! (PARTICLES,4)

!     Maximum number of particles in a search grid
      INTEGER MAX_PART_IN_GRID

!     Size of search grids
      DOUBLE PRECISION SEARCH_GRID_SIZE(3)

!     Use new algorithm for inelastic collisions
      LOGICAL USE_COL_MW

!     Van der Waals constants
      LOGICAL VAN_DER_WAALS
      DOUBLE PRECISION HAMAKER_CONSTANT
      DOUBLE PRECISION VDW_INNER_CUTOFF ! (in cm)
      DOUBLE PRECISION VDW_OUTER_CUTOFF
      DOUBLE PRECISION WALL_HAMAKER_CONSTANT
      DOUBLE PRECISION WALL_VDW_INNER_CUTOFF
      DOUBLE PRECISION WALL_VDW_OUTER_CUTOFF
      DOUBLE PRECISION SURFACE_ENERGY
      DOUBLE PRECISION WALL_SURFACE_ENERGY
      DOUBLE PRECISION Asperities ! average radius of asperities (default zero)      

!     Parameters to control Rhodes (2001) cohesion model (10/16/03)
      LOGICAL RHODES_COHESION
      DOUBLE PRECISION RHODES_COHESION_FACTOR
      DOUBLE PRECISION RHODES_COHESION_FACTOR_WALL
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE
      DOUBLE PRECISION RHODES_COHESION_LENGTH_SCALE_WALL

!     Variables used in net force measurements
      LOGICAL RECORD_NET_FORCES
      INTEGER VERTICAL_NET_FORCE_BINS
      DOUBLE PRECISION NET_FORCE_COUNTER
      DOUBLE PRECISION VERTICAL_NET_FORCE_INCREMENT
      DOUBLE PRECISION NET_FORCE_TIME_INCREMENT
      DOUBLE PRECISION TIME_AT_PREVIOUS_NET_FORCE
      DOUBLE PRECISION NET_PART_TAN_FORCE(20,2) !Contributions to net
      DOUBLE PRECISION NET_WALL_TAN_FORCE(20,2) !  force in x- and y-
      DOUBLE PRECISION NET_PART_NORM_FORCE(20,2) !  directions from up to
      DOUBLE PRECISION NET_WALL_NORM_FORCE(20,2) !  20 different heights
      DOUBLE PRECISION NET_GRAVITY_FORCE(20,2)
      DOUBLE PRECISION NET_PART_COH_FORCE(20,2)
      DOUBLE PRECISION NET_WALL_COH_FORCE(20,2)
      DOUBLE PRECISION NET_FLUID_DRAG_FORCE(20,2)
      DOUBLE PRECISION NET_PRESSURE_FORCE(20,2)
      INTEGER APP_COH_DIST_INT
      INTEGER CAP_COH_DIST_INT
      INTEGER ESC_COH_DIST_INT

! END Cohesion

! DES Attrition Model

      LOGICAL DES_Attrition
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DESRadiusNew ! (PARTICLES + more ???)
      DOUBLE PRECISION FractureHardness
      DOUBLE PRECISION FractureToughness
      DOUBLE PRECISION DESAlpha
      DOUBLE PRECISION ABRAlpha
      DOUBLE PRECISION DESAttritionThresh

! END DES Attrition

    END MODULE DISCRETELEMENT
