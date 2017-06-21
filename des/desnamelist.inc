       NAMELIST / DES_INPUT_DATA /&
!
! USE discretelement
	DISCRETE_ELEMENT, DES_CONTINUUM_COUPLED, DES_INTERP_ON, &
	 DES_ONEWAY_COUPLED, &
        TSUJI_DRAG, DES_INTG_METHOD, DES_CONTINUUM_HYBRID,&
	PARTICLES, DIMN, NFACTOR, &
	GENER_PART_CONFIG, VOL_FRAC, DES_EPS_XSTART,&
       	DES_EPS_YSTART, DES_EPS_ZSTART, &
        pgrad, pvel_mean, pvel_stdev, &
	MN, PARTICLES_FACTOR, &
	DES_NEIGHBOR_SEARCH, NEIGHBOR_SEARCH_N,&
        FACTOR_RLM, NEIGHBOR_SEARCH_RAD_RATIO,&
!	
! Used for grid/cell based search methods
        DESGRIDSEARCH_IMAX, DESGRIDSEARCH_JMAX, &
	DESGRIDSEARCH_KMAX, &
! 
! Used for quadtree or octree neighbor search methods
	MQUAD_FACTOR, QLM, QLN, INIT_QUAD_COUNT,&
!	
! Particle properties	
        KN, KN_W, KT_FAC, KT_W_FAC, MEW, MEW_W, & 
        DES_EN_INPUT, DES_ET_INPUT, DES_ETAT_FAC, &
        DES_EN_WALL_INPUT, DES_ET_WALL_INPUT, DES_ETAT_W_FAC, &
!      
! Additional properties for Hertzian model: T.Li
        DES_COLL_MODEL, ew_young, vw_poisson, e_young, v_poisson,&
!	
! Not presently used
	DTSOLID_FACTOR,&
!
! Wall conditions
        DES_PERIODIC_WALLS, DES_PERIODIC_WALLS_X,&
        DES_PERIODIC_WALLS_Y, DES_PERIODIC_WALLS_Z,&
        WALLDTSPLIT,&
! Lees Edwards BC (for DEM only)
        DES_LE_BC, DES_LE_REL_VEL, DES_LE_SHEAR_DIR,&
!	
! Wall conditions not presently used
        WALLREFLECT, NON_RECT_BC,& 
!
! DES wall BC: vibration parameters 
	DES_F, DES_GAMMA, LID_VEL,&
!
! DES wall BC: wall translational velocity
        DES_BC_Uw_s, DES_BC_Vw_s, DES_BC_Ww_s,&
!	
! DES flow boundary condition information: J.Musser
        MAX_PIS, &
        DES_BC_X_w, DES_BC_X_e, &
        DES_BC_Y_s, DES_BC_Y_n, &
        DES_BC_Z_b, DES_BC_Z_t, &
        DES_BC_VOLFLOW_s, DES_BC_MASSFLOW_s, &
        DES_BC_ROP_s, DES_BC_TYPE, &
        DES_BC_T_s, DES_BC_X_s, &
        DES_BC_CORE_Rho, &
       	FORCE_ORD_BC, &
!
! DES initial condition information: J.Musser
        DES_IC_X_w, DES_IC_X_e, &
        DES_IC_Y_s, DES_IC_Y_n, &
        DES_IC_Z_b, DES_IC_Z_t, &
        DES_IC_T_s, DES_IC_X_s, &
        DES_IC_CORE_Rho,        &
!
! DES output control information
        DEBUG_DES, FOCUS_PARTICLE,&
        PRINT_DES_DATA, DES_SPX_DT, DES_RES_DT, &
        DES_OUTPUT_TYPE, DES_CALC_BEDHEIGHT, &
!
! DES Cohesion Information: M.Weber
	MASTER_WELL_DEPTH, MASTER_WALL_WELL_DEPTH,&
        RADIUS_RATIO, WALL_RADIUS_RATIO, USE_COHESION, SQUARE_WELL,&
	MAX_LOG_STEP, VAN_DER_WAALS, HAMAKER_CONSTANT, VDW_INNER_CUTOFF,& 
        VDW_OUTER_CUTOFF, WALL_HAMAKER_CONSTANT,&
        WALL_VDW_INNER_CUTOFF, WALL_VDW_OUTER_CUTOFF,COHESION_DEBUG,&
        MAX_LOG_STEP, MIN_LOG_STEP, COHESION_DEBUG_START,Asperities,&
!
! DES Energy Equations: J.Musser
 DES_ENERGY_EQ, DES_CONV_EQ, DES_COND_EQ, DES_RADI_EQ,&
        DES_COND_EQ_PFP, DES_COND_EQ_PP, DES_CONV_CORR, FLPC, RDPC, & 
        DES_K_s0, DES_C_ps0, DES_Em, DES_MIN_COND_DIST, &
! DES Species Equations (reactive chemistry): J. Musser
       DES_NMAX,  DES_MW_s, DES_SPECIES_EQ, &
       DES_SPECIES_NAME, REACTION_MODEL, &
!
! Obsolete Controls
        ANIMATION_WRITE_INTERVAL, ANIMATION_WRITE_COUNTER, TIME_FACTOR,&
	MAX_ANIMATOR_STEP, MIN_ANIMATOR_STEP, & 
!


! MP-PIC related flags 

	 MPPIC, MPPIC_SOLID_STRESS_SNIDER, PSFAC_FRIC_PIC, FRIC_EXP_PIC, FRIC_NON_SING_FAC, CONSTANTNPC, &
	 CONSTANTWT, NPC_PIC, STATWT_PIC, CFL_PIC , MPPIC_COEFF_EN, MPPIC_PDRAG_IMPLICIT, &

! DES Attrition model flags and data

        DES_Attrition, FractureHardness, FractureToughness, DESAlpha, ABRAlpha, &
	 DESAttritionThresh, &

! des_cellwise_bc_data size 
	 MAX_DES_BC_CELL

