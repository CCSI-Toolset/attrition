!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CALC_FORCE_DES                                         C
!
!  Purpose: DES calculations of force acting on a particle, 
!           its velocity and its position                  
!
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 06-Dec-06  C
!  Reviewer: Rahul Garg                               Date: 02-Aug-07  C
!  Comments: Now includes particle-wall interaction history.           
!  Modified by Dave Decroix: add the attrition model
!           by Wesley Xu: add an temporary radius array into the 
!                         attrition model arguments                    C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CALC_FORCE_DES

      USE run      
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      IMPLICIT NONE

!-----------------------------------------------
! Local variables
!-----------------------------------------------      
      INTEGER I, J, LL, II, IW,Idim
      INTEGER NI, NLIM, N_NOCON, NEIGH_L
      INTEGER OVERLAP_MAXP
      INTEGER WALLCONTACT, WALLCHECK

      DOUBLE PRECISION OVERLAP_N, OVERLAP_T, SQRT_OVERLAP
      DOUBLE PRECISION FRAC_OVERLAP1, FRAC_OVERLAP2
! percent of particle radius when excess overlap will be flagged
      DOUBLE PRECISION, PARAMETER :: flag_overlap = 0.20d0

      DOUBLE PRECISION V_REL_TRANS_NORM, V_REL_TRANS_TANG, &
                       V_REL_NORM_OLD, VRN_OLD(DIMN)
      DOUBLE PRECISION FT_TMP(DIMN), PFT_TMP(DIMN)
! local normal and tangential forces      
      DOUBLE PRECISION FNS1(DIMN), FNS2(DIMN)
      DOUBLE PRECISION FTS1(DIMN), FTS2(DIMN)
! tmp variables to calculate magnitude of forces
      DOUBLE PRECISION FTMD, FNMD
      DOUBLE PRECISION NORMAL(DIMN), NORMAL_OLD(DIMN), TANGENT(DIMN)
      DOUBLE PRECISION DIST(DIMN), DIST_OLD(DIMN), DISTMOD, &
                       DISTMOD_OLD, R_LM
      DOUBLE PRECISION DTSOLID_TMP

! attrition temporary radius array for desradiusnew
      DOUBLE PRECISION desradius

! logic flag telling whether contact pair is old      
      LOGICAL ALREADY_NEIGHBOURS
! logic flag for local debug warnings
      LOGICAL DES_LOC_DEBUG
      LOGICAL ALREADY_EXISTS
! index to track accounted for particles      
      INTEGER PC
! store solids phase index of particle (i.e. pijk(np,5))
      INTEGER PHASEI, PHASELL
! local values used spring constants and damping coefficients
      DOUBLE PRECISION ETAN_DES, ETAN_DES_W, ETAT_DES, ETAT_DES_W,&
                       KN_DES, KN_DES_W, KT_DES, KT_DES_W
! local values used for calculating cohesive forces
      DOUBLE PRECISION FORCE_COH, EQ_RADIUS, DistApart, Norm_Dist, magGravity                       

! temp variable for tangential displacement calculation
! added by Tingwen on Feb 18, 2010
      DOUBLE PRECISION SIGMAT(DIMN),SIGMAT_OLD(DIMN)
      DOUBLE PRECISION NORM_OLD(DIMN),TANG_OLD(DIMN),TANG_NEW(DIMN),TMP_AX(DIMN)
      DOUBLE PRECISION TMP_MAG
      
! Temporary variables to store wall position, velocity and normal vector
      double precision w_pos_l(dimn),w_vel_l(dimn)   
     
! local value for particle radius
      DOUBLE PRECISION DES_R_L, DES_R_I
! local value for particle position
      DOUBLE PRECISION DES_POS_L(DIMN),DES_POS_I(DIMN)
      DOUBLE PRECISION DES_VEL_L(DIMN),DES_VEL_I(DIMN)
!      DOUBLE PRECISION DES_OME_L(DIMN),DES_OME_I(DIMN)
      
! Run time logic. is set to T when a sliding contact occurs
      LOGICAL PARTICLE_SLIDE  

      DOUBLE PRECISION TMP_FN(DIMN)

      DOUBLE PRECISION CROSSP(DIMN)
	  
!	  DOUBLE PRECISION omp_tow_tmp(DIMN)

!$      double precision omp_start, omp_end
!$      double precision omp_get_wtime     

 
!-----------------------------------------------      
! Functions
!-----------------------------------------------      
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT 

!-----------------------------------------------      

!!$      omp_start=omp_get_wtime()

! Calculate new values
!---------------------------------------------------------------------
      OVERLAP_MAXP = UNDEFINED_I
      DES_LOC_DEBUG = .FALSE.
      OVERLAP_MAX = ZERO
      NEIGH_MAX = -1
!!$   write(*,*) 'DEM call in calc_force_des.f'
	
      IF (S_TIME.LE.DTSOLID) THEN
         TANGENT(:) = ZERO
         NORMAL(:) = ZERO
         FC(:,:) = ZERO
         FN(:,:) = ZERO
         FT(:,:) = ZERO
!		 write(*,*)'calc_force_init_variables'
      ENDIF
! initialize cohesive forces
	IF(USE_COHESION) THEN
	 Fcohesive(:,:) = ZERO
	 PostCohesive (:) = ZERO
	ENDIF

!     Not sure why above is kept
	
!     Calculate contact force and torque
!---------------------------------------------------------------------
     
!      PC = 1      
!       by Tingwen      
!!$      omp_start=omp_get_wtime()

!$omp   parallel default(shared)                                  & 
!$omp   private(ll,fts1,fts2,fns1,fns2,ft_tmp,pft_tmp,            &
!$omp          PARTICLE_SLIDE,nlim,                               &
!$omp          des_r_l,des_r_i,des_pos_l,des_pos_i,               &
!$omp          des_vel_l,des_vel_i,                               &   !des_ome_l,des_ome_i,    
!$omp          n_nocon,ni,wallcheck,iw,wallcontact,i,             &
!$omp          already_neighbours,neigh_l,                        &
!$omp          w_pos_l,w_vel_l,                                   &
!$omp          r_lm,dist,distmod,frac_overlap1,normal,            &
!$omp          v_rel_trans_norm,tangent,                          &
!$omp          v_rel_trans_tang,                                  &
!$omp          overlap_n,overlap_t,dtsolid_tmp,phasell,           &
!$omp          sqrt_overlap,kn_des_w,kt_des_w,etan_des_w,         &
!$omp          etat_des_w,sigmat_old,norm_old,tmp_ax,tmp_mag,     &
!$omp          tang_old,tang_new,sigmat,                          &
!$omp          tmp_fn,ftmd,fnmd,                                  &
!$omp          crossp,ii,frac_overlap2,                           &
!$omp          dist_old,distmod_old,normal_old,                   &
!$omp          vrn_old,v_rel_norm_old,                            &
!$omp          phasei,kn_des,kt_des,etan_des,etat_des,            &
!$omp          Idim,force_coh,eq_radius,distapart,norm_dist,maggravity)          
!$omp do reduction(max:NEIGH_MAX,OVERLAP_MAX) schedule (dynamic,50)
      DO LL = 1, MAX_PIP
      
!!         IF(PC .GT. PIP) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE
!!         pc = pc+1
         IF(PEA(LL,4)) CYCLE
         
         IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
            IF (.NOT.DES_LOC_DEBUG) THEN
               DES_LOC_DEBUG = .TRUE.
               WRITE(*,1000)
            ENDIF
            WRITE(*,'(5X,A,I10)') 'On Particle ', LL
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y POS: ', DES_POS_NEW(LL,1), DES_POS_NEW(LL,2)
            WRITE(*,'(5X,A,2(ES15.7))') &
               'X,Y VEL: ', DES_VEL_NEW(LL,1), DES_VEL_NEW(LL,2)
         ENDIF

         FTS1(:) = ZERO
         FTS2(:) = ZERO
         FNS1(:) = ZERO
         FNS2(:) = ZERO

         FT_TMP(:) = ZERO
         PFT_TMP(:) = ZERO

         PARTICLE_SLIDE = .FALSE.

         IF(PN(LL,1).GE.1) THEN
            NLIM = PN(LL,1)+1
            N_NOCON = 0
! For each particle listed as in contact with particle LL in array PN,
! check the flag array PV to determine if particles remained in contact 
! after the previous call of CALC_FORCE_DES. 
            DO NI = 2, NLIM
               IF(PV(LL,NI-N_NOCON).EQ.0) THEN
! For each particle in PN(LL,2:MAXNEIGHBORS) that is no longer in 
! contact, shift the remaining particle contact information PN, PN,
! PFT left by one and reduce PN(LL,1) by one.
                  PN(LL,(NI-N_NOCON):MAXNEIGHBORS-1) = &
                     PN(LL,(NI-N_NOCON+1):MAXNEIGHBORS) 
                  PV(LL,(NI-N_NOCON):(MAXNEIGHBORS-1)) = &
                     PV(LL,(NI-N_NOCON+1):MAXNEIGHBORS)
                  PFT(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFT(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
!   added by Tingwen to save the normal direction at previous time step
                  PFN(LL,(NI-N_NOCON):(MAXNEIGHBORS-1),:) = &
                     PFN(LL,(NI-N_NOCON+1):MAXNEIGHBORS,:)
                  N_NOCON = N_NOCON + 1
                  PN(LL,1) = PN(LL,1) - 1
               ENDIF
            ENDDO
         ENDIF

! Initializing rest of the neighbor list which is not in contact and
! clean up after the above array left shifts
         NLIM = MAX(2,PN(LL,1) + 2) 
         PN(LL,NLIM:MAXNEIGHBORS) = -1
         PFT(LL,NLIM:MAXNEIGHBORS,:) = ZERO
!   added by Tingwen to save the normal direction at previous time step
         PFN(LL,NLIM:MAXNEIGHBORS,:) = ZERO   
         IF (PN(LL,1) .GT. NEIGH_MAX) NEIGH_MAX = PN(LL,1)


! Initializing the neighbor list contact information when particles are
! not in contact; i.e. when particle LL has no neighbors         
         IF (PN(LL,1).EQ.0) THEN
            PFT(LL,:,:) = ZERO
!   added by Tingwen to save the normal direction at previous time step
            PFN(LL,:,:) = ZERO
         ENDIF
         
! Reset the flag array PV; during each call to calc_force_des this
! variable tracks whether particle LL has any current neighbors
! the array is used in the next call to calc_force_des to update
! particle LL neighbor history above
         PV(LL,2:MAXNEIGHBORS) = 0

         DES_R_L=DES_RADIUS(LL)
         DES_POS_L(:)=DES_POS_NEW(LL,:)
         DES_VEL_L(:)=DES_VEL_NEW(LL,:)
!         DES_OME_L(:)=OMEGA_NEW(LL,:)
 
! Check particle LL for wall contacts
!---------------------------------------------------------------------
! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties; make sure the particle is not classified as
! a new 'entering' particle or is already marked as a potential exiting particle
         IF(WALLDTSPLIT .AND. .NOT.PEA(LL,2) .AND. .NOT.PEA(LL,3)  &
          )then
!!           .and.(pijk(ll,6).ne.0)) THEN	      

            WALLCHECK = 0
            DO IW = 1, NWALLS
               WALLCONTACT = 0
! Check to see if a particle is in contact with any of the walls
               CALL CFWALLCONTACT(IW, LL, WALLCONTACT)
               
               IF(WALLCONTACT.EQ.1) THEN
                  WALLCHECK = 1

! J.Musser : changed particles to max_pip                  
                  I = MAX_PIP + IW
                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF                        
                     ENDDO
                  ENDIF
                  
! Assign the wall particle a position and velocity
                  CALL CFWALLPOSVEL(LL, IW, w_pos_l, w_vel_l)                

                  R_LM = DES_R_L + DES_R_L 
                  DIST(:) = w_pos_l(:) - DES_POS_L(:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))
                  
! compute particle-wall VDW cohesive short-range forces
		  IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
		    DistApart = (DISTMOD-R_LM) ! distance between particle&wall surface
		    IF(DistApart < WALL_VDW_OUTER_CUTOFF)THEN
                      IF(DistApart > WALL_VDW_INNER_CUTOFF)THEN
                         FORCE_COH = WALL_HAMAKER_CONSTANT * DES_RADIUS(LL) / &
			             (6d0*DistApart**2) * &
                                     ( Asperities/(Asperities+DES_RADIUS(LL)) + &
			               ONE/(ONE+Asperities/DistApart)**2 )
                      ELSE

                         FORCE_COH = 4d0 * PI * WALL_SURFACE_ENERGY * DES_RADIUS(LL) * &
		  	            ( Asperities/(Asperities+DES_RADIUS(LL)) + &
			              ONE/(ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                      END IF                       
                      DO Idim=1,DIMN
                        Norm_Dist = DIST(Idim)/DISTMOD
		        Fcohesive(LL, Idim) = Fcohesive(LL, Idim) + Norm_Dist*FORCE_COH
                      END DO 
                    ENDIF     
		  ENDIF ! for using VDW cohesion model
                  

                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN 
                     
                     IF(((R_LM-DISTMOD)/R_LM*100.d0).GT.OVERLAP_MAX)THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
!                        OVERLAP_MAXP = LL
                     ENDIF

                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_R_L 
                     IF (FRAC_OVERLAP1 > flag_overlap) THEN
                        WRITE(*,'(5X,A,A,ES15.7)') &
                           'WARNING: excessive overlap detected ', &
                           'at time ', S_TIME
                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                           'between particle ', LL, 'and wall ',&
                           IW, 'with'
                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                          'overlap = ', (R_LM-DISTMOD), &
                           ' radius = ', DES_RADIUS(LL)
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF                             
                        WRITE(*,'(5X,A,I10,A,I5)') &
                           'DISTMOD is zero between particle ', LL, &
                           ' and wall ', IW
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair
                     CALL CFRELVEL_wall(LL, w_vel_l, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD)
                                         

! The normal overlap calculation was changed so that it no longer
! depends on the contact history (i.e., integration of incremental
! overlap found by velocity*dtsolid).  Now overlap is based purely on
! position of neighbors.  The change was made because the former
! method was found not to conserve energy 
                     OVERLAP_N =  R_LM-DISTMOD 

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                          DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                             'for first contact between particle', LL, &
                             'and wall ', IW, 'with'
                          WRITE(*,'(7X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF

                  ELSE
                     GOTO 200
                  ENDIF  !IF(R_LM - DISTMOD.GT.SMALL_NUMBER)

                  phaseLL = PIJK(LL,5) 

! T.Li : Hertz vs linear spring-dashpot contact model
                  IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
                     sqrt_overlap = SQRT(OVERLAP_N)
                     KN_DES_W = hert_kwn(phaseLL)*sqrt_overlap
                     KT_DES_W = hert_kwt(phaseLL)*sqrt_overlap
                     sqrt_overlap = SQRT(sqrt_overlap)
                     ETAN_DES_W = DES_ETAN_WALL(phaseLL)*sqrt_overlap
                     ETAT_DES_W = DES_ETAT_WALL(phaseLL)*sqrt_overlap
                  ELSE
                     KN_DES_W = KN_W
                     KT_DES_W = KT_W
                     ETAN_DES_W = DES_ETAN_WALL(phaseLL)
                     ETAT_DES_W = DES_ETAT_WALL(phaseLL)
                  ENDIF

                  FNS1(:) = -KN_DES_W * OVERLAP_N * NORMAL(:)
                  FNS2(:) = -ETAN_DES_W * V_REL_TRANS_NORM * NORMAL(:)
                  FN(LL,:) = FNS1(:) + FNS2(:) 

! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
! commented on Feb 18, 2010 by Tingwen Li
!                  PFT(LL,NI,:) = PFT(LL,NI,:)+OVERLAP_T*TANGENT(:)
!                  PFT_TMP(:) = PFT(LL,NI,:)   ! for an easy pass to des_dotprdct
!                  PFT_TMP(:) = PFT(LL,NI,:) - &
!                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
! -------------------------------------------------------------------------------------------------------
! Tingwen: New procedure to calculate the tangential displacement according to van der Hoef et al. (2006)
                  sigmat_old(:) = pft(ll,ni,:)                                        
                  norm_old(:) = pfn(ll,ni,:)
! calculate the unit vector for axies of rotation 
                if(dimn.eq.3)then
                  call des_crossprdct(tmp_ax,norm_old,normal)
                  tmp_mag=des_dotprdct(tmp_ax,tmp_ax)
                  if(tmp_mag .gt. zero)then
                    tmp_ax(:)=tmp_ax(:)/sqrt(tmp_mag)
! get the old tangential direction unit vector
                    call des_crossprdct(tang_old,tmp_ax,norm_old)
! get the new tangential direction unit vector due to rotation
                    call des_crossprdct(tang_new,tmp_ax,normal)
                    sigmat(:)=des_dotprdct(sigmat_old,tmp_ax)*tmp_ax(:) &
                              + des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                    sigmat(:)=sigmat(:)+overlap_t*tangent(:)                    
                  else
                    sigmat(:)=sigmat_old(:)+overlap_t*tangent(:)                                  
                  end if
                else
                  tang_old(1) =-norm_old(2)
                  tang_old(2) = norm_old(1)
                  tang_new(1) =-normal(2)
                  tang_new(2) = normal(1)
                  sigmat(:)=des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                  sigmat(:)=sigmat(:)+overlap_t*tangent(:)                    
                endif     
      
                  pft_tmp(:) =sigmat(:)
! end of new procedure
! -------------------------------------------------------------------------------------------------------

                  FTS1(:) = -KT_DES_W * PFT_TMP(:)
                  FTS2(:) = -ETAT_DES_W * V_REL_TRANS_TANG * TANGENT(:)
                  FT(LL,:) = FTS1(:) + FTS2(:) 

                  FT_TMP(:) = FT(LL,:)
                                                     
! Check for Coulombs friction law and limit the maximum value of the
! tangential force on a particle in contact with a wall
                  CALL CFSLIDEWALL(LL, TANGENT,PARTICLE_SLIDE)                                  
                  
! Calculate the total force FC and TOW on a particle in a particle-wall
! collision
                  CALL CFFCTOWALL(LL, NORMAL, DISTMOD)                                  

! Save the old normal direction
                  pfn(ll,ni,:) = normal(:)
                  
! Save the tangential displacement history with the correction of Coulomb's law
                  IF (PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = -( FT(LL,:) - FTS2(:) ) / KT_DES_W
                  ELSE
                     PFT(LL,NI,:) = PFT_TMP(:)
                  ENDIF
        
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000)
                     ENDIF                          
                     WRITE(*,*) '     STIME, DTSOLID = ', S_TIME, DTSOLID
                     WRITE(*,*) '     WALL CONTACT ON WALL =', IW
                     WRITE(*,*) '     ALREADY_NEIGHBOURS = ',&
                        ALREADY_NEIGHBOURS
                     WRITE(*,*) '     DES_VEL = ', DES_VEL_NEW(LL,1:DIMN),&
                        des_radius(LL)*OMEGA_NEW(LL,1)
                     WRITE(*,*) '     V-OMEGA R = ', &
                        DES_VEL_NEW(LL,1)+des_radius(LL)* OMEGA_NEW(LL,1),&
                        (DES_VEL_NEW(LL,1)+des_radius(LL)*OMEGA_NEW(LL,1))*DTSOLID
                     WRITE(*,*) '     M*g = ', PMASS(LL)*gravity
                     WRITE(*,*) '     KN_W, ETAN_W, KT_W, ETAT_W = ',&
                        KN_DES_W, ETAN_DES_W, KT_DES_W, ETAT_DES_W
                     WRITE(*,*) '     TANGENT= ', TANGENT
                     WRITE(*,*) '     HIST = ', PFT(LL,NI,1:2)
                     WRITE(*,*) '     PARTICLE_SLIDE ? ', PARTICLE_SLIDE
                     WRITE(*,*) '     FT and FN= ', FT( LL,:), FN(LL,:)
                     WRITE(*,*) '     KT_W*OT*TAN = ', &
                        KT_DES_W*((OVERLAP_T)) *TANGENT(:)
                     WRITE(*,*) '     OVERLAP_T = ', OVERLAP_T, TANGENT
                     FTMD = SQRT(DES_DOTPRDCT(FT_TMP,FT_TMP))
                     FNMD = SQRT(DES_DOTPRDCT(FN(LL,1:DIMN),FN(LL,1:DIMN)))
                     WRITE(*,*) '     FTMD, mu FNMD = ', FTMD, MEW_W*FNMD
                  ENDIF

                  PARTICLE_SLIDE = .FALSE.

               ENDIF   !wall contact
 200           CONTINUE
!
! Call Attrition model after completion of particle-wall calculation
!
                  IF(DES_ATTRITION) THEN
		    CALL CALC_ATTRITION_DES(LL,FN(LL,:),PFT_TMP(:),&
		    V_REL_TRANS_NORM,V_REL_TRANS_TANG,ALREADY_NEIGHBOURS,DESRadiusNew(LL))                 
                  ENDIF  
!
            ENDDO ! DO IW = 1, NWALLS
         ENDIF   !if(walldtsplit .and. .not.pea(LL,2))
!---------------------------------------------------------------------
! End check particle LL for wall contacts   
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  START PARTICLE-PARTICLE CONTACT 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Check particle LL neighbour contacts         
!---------------------------------------------------------------------
         IF (NEIGHBOURS(LL,1).GT.0) THEN
            DO II = 2, NEIGHBOURS(LL,1)+1
               I = NEIGHBOURS(LL,II)
!pradeep removed the check I.gt.LL
!               IF(I.GT.LL .AND. PEA(I,1)) THEN
!               IF( (PEA(I,4) .OR. I.GT.LL) .AND. PEA(I,1)) THEN
               IF(PEA(I,1)) THEN                  
                  ALREADY_NEIGHBOURS=.FALSE.
                  
                  IF(PN(LL,1).GT.0) THEN                     
                     DO NEIGH_L = 2, PN(LL,1)+1
                        IF(I.EQ. PN(LL,NEIGH_L)) THEN 
                           ALREADY_NEIGHBOURS=.TRUE.
                           NI = NEIGH_L
                           EXIT
                        ENDIF
                     ENDDO
                  ENDIF

         DES_R_I=DES_RADIUS(I)
         DES_POS_I(:)=DES_POS_NEW(I,:)
         DES_VEL_I(:)=DES_VEL_NEW(I,:)
!         DES_OME_I(:)=OMEGA_NEW(I,:)                  
!     Tingwen
!     Did not look into the treatment for peroidic wall
! Pradeep particles will be copied for periodic boundaries 
!                  IF(DES_PERIODIC_WALLS) THEN
!                     TEMPX = DES_POS_NEW(I,1)
!                     TEMPY = DES_POS_NEW(I,2)
!                     IF(DIMN.EQ.3) TEMPZ = DES_POS_NEW(I,3)        
!
!                     TEMPD = ABS(DES_POS_NEW(LL,1) - DES_POS_NEW(I,1))
!                     IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_X) THEN 
!                        IF(DEBUG_DES) THEN
!                           IF (.NOT.DES_LOC_DEBUG) THEN
!                              DES_LOC_DEBUG = .TRUE.
!                              WRITE(*,1000)
!                           ENDIF                                
!                           WRITE(*,'(5X,A,I,X,I)') &
!                              'PARTICLE LL & NEIGHBOR I: ', LL, I
!                           WRITE(*,'(5X,A,(ES15.7))') &
!                              'LL DES_POS = ', DES_POS_NEW(LL,:)
!                           WRITE(*,'(5X,A,(ES15.7))') &
!                              'I DES_POS = ', DES_POS_NEW(I,:)
!                        ENDIF
!                        IF(TEMPX.GT.DES_POS_NEW(LL,1)) THEN 
!                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) - (EX2-WX1)
!                           IF(DEBUG_DES) THEN
!                              IF (.NOT.DES_LOC_DEBUG) THEN
!                                 DES_LOC_DEBUG = .TRUE.
!                                 WRITE(*,1000) 
!                              ENDIF                                
!                              WRITE(*,*) '     NEW I POS WEST= ', DES_POS_NEW(I,1)
!                           ENDIF
!                        ELSE
!                           DES_POS_NEW(I,1) = DES_POS_NEW(I,1) + EX2 - WX1
!                           IF(DEBUG_DES) THEN
!                              IF (.NOT.DES_LOC_DEBUG) THEN
!                                 DES_LOC_DEBUG = .TRUE.
!                                 WRITE(*,1000)
!                              ENDIF
!                              WRITE(*,*) '     NEW I POS EAST = ', DES_POS_NEW(I,1)
!                           ENDIF
!                        ENDIF
!                     ENDIF
!                     
!                     TEMPD = ABS(DES_POS_NEW(LL,2) - DES_POS_NEW(I,2))
!                     IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Y) THEN
!                        IF(TEMPY.GT.DES_POS_NEW(LL,2)) THEN 
!                           DES_POS_NEW(I,2) = DES_POS_NEW(I,2) - (TY2-BY1)
!                        ELSE
!                           DES_POS_NEW(I,2) = DES_POS_NEW(I,2) + (TY2-BY1)
!                        ENDIF
!                     ENDIF
!                     
!                     IF(DIMN.EQ.3) THEN
!                        TEMPD = ABS(DES_POS_NEW(LL,3) - DES_POS_NEW(I,3))
!                        IF(TEMPD.GT.4.d0*MAX_RADIUS.AND.DES_PERIODIC_WALLS_Z) THEN 
!                           IF(TEMPZ.GT.DES_POS_NEW(LL,3)) THEN 
!                              DES_POS_NEW(I,3) = DES_POS_NEW(I,3) -(NZ2 - SZ1)
!                           ELSE
!                              DES_POS_NEW(I,3) = DES_POS_NEW(I,3) + (NZ2-SZ1)
!                           ENDIF
!                        ENDIF
!                     ENDIF
!                  ENDIF !   IF(DES_PERIODIC_WALLS) THEN
                  
                  R_LM = DES_R_L + DES_R_I
                  DIST(:) = DES_POS_I(:) - DES_POS_L(:)
                  DISTMOD = SQRT(DES_DOTPRDCT(DIST,DIST))

! compute particle-particle VDW cohesive short-range forces	
		  IF(USE_COHESION .AND. VAN_DER_WAALS) THEN
		    EQ_RADIUS = 2d0 * DES_RADIUS(LL)*DES_RADIUS(I) / &
		               (DES_RADIUS(LL)+DES_RADIUS(I))  ! for use in cohesive force
		    DistApart = (DISTMOD-R_LM) ! distance between particle surface
		    IF(DistApart < VDW_OUTER_CUTOFF)THEN
                      IF(DistApart > VDW_INNER_CUTOFF)THEN
                         FORCE_COH = HAMAKER_CONSTANT * EQ_RADIUS / (12d0*DistApart**2) * &
                                     ( Asperities/(Asperities+EQ_RADIUS) + &
			               ONE/(ONE+Asperities/DistApart)**2 )
                      ELSE

                         FORCE_COH = 2d0 * PI * SURFACE_ENERGY * EQ_RADIUS * &
		  	            ( Asperities/(Asperities+EQ_RADIUS) + &
			              ONE/(ONE+Asperities/VDW_INNER_CUTOFF)**2 )
                      END IF                       
                      DO Idim=1,DIMN
                        Norm_Dist = DIST(Idim)/DISTMOD
		        Fcohesive(LL, Idim) = Fcohesive(LL, Idim) + Norm_Dist*FORCE_COH
!		        Fcohesive(I, Idim) =  Fcohesive(I, Idim)  - Norm_Dist*FORCE_COH
                      END DO 
                    ENDIF     
		  ENDIF ! for using VDW cohesion model                  
                  
!                  IF(DES_PERIODIC_WALLS) THEN
!                     DES_POS_NEW(I,1) = TEMPX
!                     DES_POS_NEW(I,2) = TEMPY
!                     IF (DIMN.EQ.3) DES_POS_NEW(I,3) = TEMPZ  
!                  ENDIF

                  IF(R_LM - DISTMOD.GT.SMALL_NUMBER) THEN

                     IF(DEBUG_DES .AND. LL.EQ.FOCUS_PARTICLE) THEN
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000) 
                        ENDIF
                        WRITE(*,'(5X,A,(10I10))') 'NEIGHBORS: ', NEIGHBOURS(LL,:)
                     ENDIF

                     IF((((R_LM-DISTMOD)/R_LM)*100.d0).GT.OVERLAP_MAX)THEN
                        OVERLAP_MAX = (((R_LM-DISTMOD)/R_LM)*100.d0)
!                        OVERLAP_MAXP = LL
                     ENDIF
                     
                     FRAC_OVERLAP1 = (R_LM-DISTMOD)/DES_R_L !DES_RADIUS(LL)
                     FRAC_OVERLAP2 = (R_LM-DISTMOD)/DES_R_I !DES_RADIUS(I)
                     IF (FRAC_OVERLAP1 > flag_overlap .OR. &
                         FRAC_OVERLAP2 > flag_overlap) THEN
                        WRITE(*,'(5X,A,A,ES15.7)') &
                           'WARNING: excessive overlap detected ', &
                           'at time ', S_TIME
                        WRITE(*,'(7X,A,I10,2X,A,I5,2X,A)') &
                           'between particles ', LL, 'and ',&
                           I, 'with'
                        WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7,2X,ES15.7)') &
                           'overlap = ', (R_LM-DISTMOD), &
                           ' radii = ', DES_RADIUS(LL), DES_RADIUS(I)
                     ENDIF

                     IF(DISTMOD.NE.ZERO) THEN
                        NORMAL(:)= DIST(:)/DISTMOD
                     ELSE 
                        IF (.NOT.DES_LOC_DEBUG) THEN
                           DES_LOC_DEBUG = .TRUE.
                           WRITE(*,1000)
                        ENDIF
                        WRITE(*,'(5X,A,I10,I10)') &
                           'DISTMOD is zero between particle-pair ',&
                           LL, I
                        STOP
                     ENDIF

! Calculate the translational relative velocity for a contacting particle pair.
                     CALL CFRELVEL(LL, I, V_REL_TRANS_NORM, &
                        V_REL_TRANS_TANG, TANGENT, NORMAL, DISTMOD)                                                       

! Overlap calculation changed from history based to current position 
                     OVERLAP_N = R_LM-DISTMOD

                     IF(ALREADY_NEIGHBOURS) THEN 
                        PV(LL,NI) = 1
                        OVERLAP_T = V_REL_TRANS_TANG*DTSOLID
                     ELSE 
                        IF(DEBUG_DES) THEN
                           IF (.NOT.DES_LOC_DEBUG) THEN
                              DES_LOC_DEBUG = .TRUE.
                              WRITE(*,1000)
                           ENDIF
                           WRITE(*,'(5X,A,2(I10,X),A,ES15.7)') &
                              'Normal overlap for particle pair ',&
                              LL, I, ' : ', OVERLAP_N 
                        ENDIF
                        PN(LL,1) = PN(LL,1) + 1
                        NI = PN(LL,1) + 1
                        PN(LL,NI) = I
                        PV(LL,NI) = 1
                        IF (V_REL_TRANS_NORM .GT. ZERO) THEN
                           DTSOLID_TMP = OVERLAP_N/(V_REL_TRANS_NORM)
                        ELSEIF (V_REL_TRANS_NORM .LT. ZERO) THEN
                          DTSOLID_TMP = DTSOLID
                          DIST_OLD(:)=DES_POS_OLD(I,:)-DES_POS_OLD(LL,:)
                          DISTMOD_OLD=SQRT(DES_DOTPRDCT(DIST_OLD,DIST_OLD))
                          NORMAL_OLD(:)=DIST_OLD(:)/DISTMOD_OLD
                          VRN_OLD(:)=DES_VEL_OLD(LL,:)-DES_VEL_OLD(I,:)
                          V_REL_NORM_OLD=DES_DOTPRDCT(VRN_OLD,NORMAL_OLD)
                          WRITE(*,'(5X,A,A,ES15.7)') &
                             'WARNING: normal relative velocity less ',&
                             'than zero at time ', S_TIME
                          WRITE(*,'(7X,A,I10,2X,A,I10,2X,A)') &
                             'for first contact between particles', LL, &
                             'and ', I, 'with'
                          WRITE(*,'(7X,A,ES15.7,2X,A,ES15.7)') &
                             'V_REL_NORM = ', V_REL_TRANS_NORM, &
                             'and V_REL_NORM_OLD = ', V_REL_NORM_OLD
                        ELSE
                           DTSOLID_TMP = OVERLAP_N/&
                              (V_REL_TRANS_NORM+SMALL_NUMBER)
                        ENDIF
                        OVERLAP_T = V_REL_TRANS_TANG*&
                           MIN(DTSOLID,DTSOLID_TMP)
                     ENDIF
                  ELSE
                     GOTO 300
                  ENDIF

                  phaseLL = PIJK(LL,5)                  
                  phaseI = PIJK(I,5)

! T.Li : Hertz vs linear spring-dashpot contact model
                  IF (TRIM(DES_COLL_MODEL) .EQ. 'HERTZIAN') THEN
                     sqrt_overlap = SQRT(OVERLAP_N)
                     KN_DES = hert_kn(phaseLL,phaseI)*sqrt_overlap
                     KT_DES = hert_kt(phaseLL,phaseI)*sqrt_overlap
                     sqrt_overlap = SQRT(sqrt_overlap)
                     ETAN_DES = DES_ETAN(phaseLL,phaseI)*sqrt_overlap
                     ETAT_DES = DES_ETAT(phaseLL,phaseI)*sqrt_overlap
                  ELSE
                     KN_DES = KN
                     KT_DES = KT
                     ETAN_DES = DES_ETAN(phaseLL,phaseI)
                     ETAT_DES = DES_ETAT(phaseLL,phaseI)
                  ENDIF

                  FNS1(:) = -KN_DES * OVERLAP_N * NORMAL(:)
                  FNS2(:) = -ETAN_DES * V_REL_TRANS_NORM*NORMAL(:)
                  FN(LL,:) = FNS1(:) + FNS2(:)       

! Calculate the tangential displacement which is integration of tangential 
! relative velocity with respect to contact time. Correction in the tangential 
! direction is imposed                   
! commented on Feb 18, 2010 by Tingwen Li
!                  PFT(LL,NI,:) = PFT(LL,NI,:) + OVERLAP_T * TANGENT(:)
!                  PFT_TMP(:) = PFT(LL,NI,:)   ! for an easy pass to des_dotprdct
!                  PFT_TMP(:) = PFT(LL,NI,:) - &
!                     DES_DOTPRDCT(PFT_TMP,NORMAL)*NORMAL(:)
! -------------------------------------------------------------------------------------------------------
! Tingwen: New procedure to calculate the tangential displacement according to van der Hoef et al. (2006)
                  sigmat_old(:) = pft(ll,ni,:)                                        
                  norm_old(:) = pfn(ll,ni,:)
! calculate the unit vector for axies of rotation
                if(dimn.eq.3)then
                  call des_crossprdct(tmp_ax,norm_old,normal)
                  tmp_mag=des_dotprdct(tmp_ax,tmp_ax)
                  if(tmp_mag .gt. zero)then
                    tmp_ax(:)=tmp_ax(:)/sqrt(tmp_mag)
! get the old tangential direction unit vector
                    call des_crossprdct(tang_old,tmp_ax,norm_old)
! get the new tangential direction unit vector due to rotation
                    call des_crossprdct(tang_new,tmp_ax,normal)
                    sigmat(:)=des_dotprdct(sigmat_old,tmp_ax)*tmp_ax(:) &
                              + des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                    sigmat(:)=sigmat(:)+overlap_t*tangent(:)                    
                  else
                    sigmat(:)=sigmat_old(:)+overlap_t*tangent(:)                                  
                  end if
    else
                  tang_old(1) =-norm_old(2)
                  tang_old(2) = norm_old(1)
                  tang_new(1) =-normal(2)
                  tang_new(2) = normal(1)
                  sigmat(:)=des_dotprdct(sigmat_old,tang_old)*tang_new(:)
                  sigmat(:)=sigmat(:)+overlap_t*tangent(:)                    
                endif     
      
                  pft_tmp(:) =sigmat(:)
! end of new procedure
! -------------------------------------------------------------------------------------------------------

                  FTS1(:) = -KT_DES * PFT_TMP(:)
                  FTS2(:) = -ETAT_DES * V_REL_TRANS_TANG * TANGENT(:)
                  FT(LL,:) = FTS1(:) + FTS2(:) 

                  FT_TMP(:) = FT(LL,:)
                  call cfslide(ll,tangent,PARTICLE_SLIDE)
                  
! Calculate the total force FC and TOW on a particle in a particle-particle collision
                  CALL CFFCTOW(LL, I, NORMAL, DISTMOD)

! Save the old normal direction
                  pfn(ll,ni,:) =normal(:)
! Save the tangential displacement history with the correction of Coulomb's law
                  IF (PARTICLE_SLIDE) THEN
                     PFT(LL,NI,:) = -( FT(LL,:) - FTS2(:) ) / KT_DES
                  ELSE
                     PFT(LL,NI,:) = PFT_TMP(:)
                  ENDIF
                  
                  IF(DEBUG_DES.AND.LL.EQ.FOCUS_PARTICLE) THEN 
                     IF (.NOT.DES_LOC_DEBUG) THEN
                        DES_LOC_DEBUG = .TRUE.
                        WRITE(*,1000) 
                     ENDIF

                     PRINT*, '     EtaN, EtaT =  ', ETAN_DES, ETAT_DES
                     PRINT*, '     Percent overlap = ', (R_LM - DISTMOD)*100.d0/R_LM
                     PRINT*, '     rad ratio = ', DES_RADIUS(LL)/DES_RADIUS(I)
                     PRINT*, '     FNS1 and FNS2 = ', FNS1(:), FNS2(:)
                     PRINT*, '     PFT = ', PFT(LL,NI,:)
                     PRINT*, '     FORCEST = ', FT(LL,:)
                     PRINT*, '     FORCESN = ', FN(LL,:)
                     PRINT*, '     FORCEST = ', FT(LL,:)
                  ENDIF

                  IF(DEBUG_DES.AND.LL.eq.FOCUS_PARTICLE)THEN
                     INQUIRE(FILE='debug_file',EXIST=ALREADY_EXISTS)
                     IF(ALREADY_EXISTS)THEN
                        OPEN(UNIT=1,FILE='debug_file',STATUS='OLD',POSITION='APPEND')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1), 'FNy=',FN(LL,2)
                     ELSE
                        OPEN(UNIT=1,FILE='debug_file',STATUS='NEW')
                        WRITE(1,'(A,I5)')'CALC FORCE -- NEIGHBOR',II
                        WRITE(1,'(2(1x,A,E12.5))')&
                        'FNx=',FN(LL,1),'FNy=',FN(LL,2)
                     ENDIF
                     CLOSE (1)
                     PRINT*, 'PN', PN(LL,:)
                  ENDIF

                  PARTICLE_SLIDE = .FALSE.

               ENDIF         ! IF (I>LL .AND. PEA(I,1))

 300           CONTINUE
!
! Call Attrition model after completion of particle-particle calculation
!
                IF(DES_ATTRITION) THEN
		    CALL CALC_ATTRITION_DES(LL,FN(LL,:),PFT_TMP(:),&
		    V_REL_TRANS_NORM,V_REL_TRANS_TANG,ALREADY_NEIGHBOURS,DESRadiusNew(LL))                 
                ENDIF  

            ENDDO            ! DO II = 2, NEIGHBOURS(LL,1)+I
         ENDIF               ! IF(NEIGHBOURS(LL,1).GT.0)

!---------------------------------------------------------------------
! End check particle LL neighbour contacts         
!
      ENDDO   ! end loop over paticles LL
!$omp end parallel      
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'calc_force_loop:',omp_end - omp_start

! Treats wall interaction also as a two-particle interaction but accounting
! for the wall properties
      IF(DES_CONTINUUM_COUPLED) THEN
         CALL DRAG_FGS
      ENDIF
! cohesion part not called for VDW to use parallel and grid search capabilities 
! now available in dem.
! The square-well model is still available in the model/cohesion directory.
! COHESION
      IF(USE_COHESION .AND. .NOT.VAN_DER_WAALS)THEN
         CALL CALC_COHESIVE_FORCES
      ENDIF

! COHESION ! just for post-processing mag. of cohesive forces on each particle
      IF(USE_COHESION)THEN
      pc = 1
      magGravity = DSQRT(DES_DOTPRDCT(GRAV,GRAV))
      DO LL = 1, MAX_PIP
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(ll,1)) cycle 
         pc = pc+1
         if(pea(ll,4)) cycle 
          DO Idim=1,DIMN
            PostCohesive(LL) =  PostCohesive(LL) + Fcohesive(LL, Idim)**2
          ENDDO
          if(magGravity> ZERO) PostCohesive(LL) =  DSQRT(PostCohesive(LL)) / &
	                                           (PMASS(LL)*magGravity)
        ENDDO 
      ENDIF ! for cohesion model      
      
! Update the old values of particle position and velocity with the new values computed
!!$      omp_start=omp_get_wtime()
      CALL CFUPDATEOLD	  
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'cfupdateold:',omp_end - omp_start

      IF (DES_LOC_DEBUG) WRITE(*,1001)

 1000 FORMAT(5X,'---------- START CALC_FORCE_DES ---------->')
 1001 FORMAT(5X,'<---------- END CALC_FORCE_DES ----------') 

      RETURN
      END SUBROUTINE CALC_FORCE_DES

