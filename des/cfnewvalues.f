!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,       
!           position, angular velocity etc                            
!
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer:                                          Date:            C
!                                                                      C 
!
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing 
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!  
!  Reviewer: Wesley Xu, Dave Decroix                   Date: 10-May-12
!  Revision: Update particle radius according to attrition calculation
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
      SUBROUTINE CFNEWVALUES

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      
      USE mfix_pic 
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, I
      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST

! index to track accounted for particles  
      INTEGER PC 

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG


!-----------------------------------------------
! Functions 
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------      
      IF(MPPIC) THEN 
         
         if(MPPIC_SOLID_STRESS_SNIDER) THEN 
            
            !CALL MPPIC_COMPUTE_PS_GRAD            

            CALL CFNEWVALUES_MPPIC_SINDER 
         ELSE
            !CALL MPPIC_COMPUTE_PS_GRAD            
            
            CALL CFNEWVALUES_MPPIC
            !the above calls the coloring function like approach 
         ENDIF
         
         RETURN 
      ENDIF
	  
         DES_LOC_DEBUG = .FALSE.
!!      PC = 1
!$omp parallel do if(max_pip .ge. 10000) default(shared)        &
!$omp private(l,d,dist,neighbor_search_dist)                    &
!$omp reduction(.or.:do_nsearch) schedule (auto)                      	  
      DO L = 1, MAX_PIP
! pradeep skip ghost particles
!!         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
!!         pc = pc+1
         if(pea(l,4)) cycle 


! If a particle is classified as new, then forces are ignored. 
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.PEA(L,2))THEN 
            FC(L, :) = FC(L,:)/PMASS(L) + GRAV(:)
            IF(USE_COHESION .AND. VAN_DER_WAALS) FC(L, :) = FC(L,:) + Fcohesive(L, :)/PMASS(L)
         ELSE 
            FC(L,:) = ZERO
            TOW(L,:) = ZERO         
         ENDIF
         
         
! Advance particle position, velocity
!! comment out by Tingwen
        IF (TRIM(DES_INTG_METHOD) .EQ. 'EULER') THEN 
! first-order method              
            DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + FC(L,:)*DTSOLID
            DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
               DES_VEL_NEW(L,:)*DTSOLID 
! following is equivalent to x=xold + vold*dt + 1/2acc*dt^2
!         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + 0.5d0*&
!             (DES_VEL_NEW(L,:)+DES_VEL_OLD(L,:))*DTSOLID 
            OMEGA_NEW(L,:)   = OMEGA_OLD(L,:) + TOW(L,:)*OMOI(L)*DTSOLID

            ELSEIF (TRIM(DES_INTG_METHOD) .EQ. 'ADAMS_BASHFORTH') THEN
! T.Li:  second-order Adams-Bashforth scheme
           DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + 0.5d0*&
               ( 3.d0*DES_VEL_OLD(L,:)-DES_VEL_OOLD(L,:) )*DTSOLID
            DES_VEL_NEW(L,:) = DES_VEL_OLD(L,:) + 0.5d0*&
               ( 3.d0*FC(L,:)-DES_ACC_OLD(L,:) )*DTSOLID
            OMEGA_NEW(L,:)   =  OMEGA_OLD(L,:) + 0.5d0*&
               ( 3.d0*TOW(L,:)*OMOI(L)-ROT_ACC_OLD(L,:) )*DTSOLID
            DES_ACC_OLD(L,:) = FC(L,:)
            ROT_ACC_OLD(L,:) = TOW(L,:)*OMOI(L)
         ENDIF

! Update particle radius if using attrition model

         IF(DES_ATTRITION) THEN
! update radius
            DES_RADIUS(L)=DESRadiusNew(L)
         ENDIF

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so, 
! make sure that neighbor is called in des_time_march
         IF(.NOT.DO_NSEARCH) THEN
            D(:) = DES_POS_NEW(L,:) - PPOS(L,:)
            DIST = SQRT(DES_DOTPRDCT(D,D))
            NEIGHBOR_SEARCH_DIST = NEIGHBOR_SEARCH_RAD_RATIO*&
               DES_RADIUS(L)
            IF(DIST.GE.NEIGHBOR_SEARCH_DIST) DO_NSEARCH = .TRUE.
         ENDIF


! Check if the particle has moved a distance greater than or equal to 
! its radius during one solids time step. if so, call stop
         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         DIST = SQRT(DES_DOTPRDCT(D,D))
         IF(DIST.GE.DES_RADIUS(L)) THEN
            WRITE(*,1002) L, DIST, DES_RADIUS(L)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'old particle pos = ', DES_POS_OLD(L,:)
            WRITE(*,'(5X,A,3(ES17.9))') &
               'new particle pos = ', DES_POS_NEW(L,:)
            WRITE(*,'(5X,A,3(ES17.9))')&
               'new particle vel = ', DES_VEL_NEW(L,:) 
            WRITE(*,1003)
            STOP
         ENDIF

! ************************************************************************
! pradeep - Check has to be consolidated and placed in one place 
!
!! Warning message for particles moving into ghost cells:
!! Note that if this occurs then the particle_in_cell 
!! subroutine will call a stop
!         IF((DES_POS_NEW(L,1) < ZERO .OR. DES_POS_NEW(L,1) > XLENGTH) .AND.&
!         .NOT.DES_PERIODIC_WALLS_X .AND. &
!         .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
!! A new or exiting particle may exist in ghost cells
!            IF (.NOT.DES_LOC_DEBUG) THEN
!               DES_LOC_DEBUG = .TRUE.
!               WRITE(*,1000) 
!            ENDIF         
!            WRITE(*,'(5X,A,I10)') &
!               'X position outside domain for particle ', L
!            WRITE(*,'(7X,A,3(ES17.9))')&
!               'particle pos = ', DES_POS_NEW(L,:)
!            WRITE(*,'(7X,A,3(ES17.9))')&
!               'particle vel = ', DES_VEL_NEW(L,:)
!         ENDIF 
!
!         IF((DES_POS_NEW(L,2) < ZERO .OR. DES_POS_NEW(L,2) > YLENGTH) .AND.&
!         .NOT.DES_PERIODIC_WALLS_Y .AND. &
!         .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
!! A new or exiting particle may exist in ghost cells
!            IF (.NOT.DES_LOC_DEBUG) THEN
!               DES_LOC_DEBUG = .TRUE.
!               WRITE(*,1000) 
!            ENDIF         
!            WRITE(*,'(5X,A,I10)') &
!               'Y position outside domain for particle=: ', L
!            WRITE(*,'(7X,A,3(ES17.9))')&
!               'particle pos = ', DES_POS_NEW(L,:)
!            WRITE(*,'(7X,A,3(ES17.9))')&
!               'particle vel = ', DES_VEL_NEW(L,:)
!         ENDIF 
!
!         IF (DIMN > 2) THEN
!            IF((DES_POS_NEW(L,3) < ZERO .OR. &
!                DES_POS_NEW(L,3) > ZLENGTH) .AND.&
!            .NOT.DES_PERIODIC_WALLS_Z .AND. &
!            .NOT.PEA(L,2) .AND. .NOT.PEA(L,3)) THEN
!! A new or exiting particle may exist in ghost cells
!               IF (.NOT.DES_LOC_DEBUG) THEN
!                  DES_LOC_DEBUG = .TRUE.
!                  WRITE(*,1000) 
!               ENDIF         
!               WRITE(*,'(5X,A,I10)') &
!                  'Z position outside domain for particle ', L
!               WRITE(*,'(7X,A,3(ES17.9))')&
!                  'particle pos = ', DES_POS_NEW(L,:)
!               WRITE(*,'(7X,A,3(ES17.9))')&
!                  'particle vel = ', DES_VEL_NEW(L,:)
!            ENDIF
!         ENDIF 
! ************************************************************************

! Pradeep- no movement of particle for periodic boundaries
! it is taken care in desmpi mod

!! Periodic treatment
!         IF(DES_PERIODIC_WALLS) THEN
!            IF(DES_PERIODIC_WALLS_X) THEN
!               IF(DES_POS_NEW(L,1).GT.EX2) THEN
!                  DES_POS_NEW(L,1) = DES_POS_NEW(L,1) - (EX2 - WX1)
!                  PIJK(L,1) = 2
!               ELSEIF(DES_POS_NEW(L,1).LT.WX1) THEN
!                  DES_POS_NEW(L,1) = DES_POS_NEW(L,1) + (EX2 - WX1)
!                  PIJK(L,1) = IMAX1
!               ENDIF
!            ENDIF
!            IF(DES_PERIODIC_WALLS_Y) THEN
!               IF(DES_POS_NEW(L,2).GT.TY2) THEN
!                  DES_POS_NEW(L,2) = DES_POS_NEW(L,2) - (TY2 - BY1)
!                  PIJK(L,2) = 2
!               ELSEIF(DES_POS_NEW(L,2).LT.BY1) THEN
!                  DES_POS_NEW(L,2) = DES_POS_NEW(L,2) + (TY2 - BY1)
!                  PIJK(L,2) = JMAX1
!               ENDIF
!            ENDIF
!            IF(DES_PERIODIC_WALLS_Z) THEN
!               IF(DES_POS_NEW(L,3).GT.NZ2) THEN
!                  DES_POS_NEW(L,3) = DES_POS_NEW(L,3) - (NZ2 - SZ1)
!                  PIJK(L,3) = 2
!               ELSEIF(DES_POS_NEW(L,3).LT.SZ1) THEN
!                  DES_POS_NEW(L,3) = DES_POS_NEW(L,3) + (NZ2 - SZ1)
!                  PIJK(L,3) = KMAX1
!               ENDIF
!            ENDIF
!         ENDIF
! ************************************************************************

! Reset total contact force and torque      
         FC(L,:) = ZERO
         TOW(L,:) = ZERO


!         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO
!$omp end parallel do
	  
 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')  

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)         

      RETURN
      END SUBROUTINE CFNEWVALUES
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES_MPPIC_SINDER                               C
!
!                                                                      C 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C     
      SUBROUTINE CFNEWVALUES_MPPIC_SINDER 

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mppic_wallbc
      USE fldvar
      
      USE mfix_pic 
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM 
      INTEGER I, J, K, IJK, IJK_OLD, IJK2, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB

      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST, DP_BAR, COEFF_EN, MEANVEL(DIMN), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(DIMN), UPRIMETAU(DIMN), UPRIMETAU_INT(DIMN), MEAN_FREE_PATH, PS_FORCE(DIMN), VEL_ORIG(DIMN)
! index to track accounted for particles  
      INTEGER PC 

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC 
      DOUBLE PRECISION MAXDIST_PIC, UPRIMEMOD, UPRIMEMODNEW, signvel

! dt's in each direction  based on cfl_pic for the mppic case 
      
      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, MEANUS(DIMN, MMAX)
      DOUBLE PRECISION :: DPS_DXE, DPS_DXW, DPS_DYN, DPS_DYS, DPS_DZT, DPS_DZB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, epg_min2, velf_part(dimn)
      
      LOGICAL :: DELETE_PART
      INTEGER :: PIP_DEL_COUNT, count_bc 

      DOUBLE PRECISION MEANUS_e(DIMN, MMAX), MEANUS_w(DIMN, MMAX),MEANUS_n(DIMN, MMAX),MEANUS_s(DIMN, MMAX),MEANUS_t(DIMN, MMAX), MEANUS_b(DIMN, MMAX)      
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5), epg_min_loc(1)
!-----------------------------------------------
! Functions 
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------      
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      PC = 1
      FOCUS_PARTICLE = -1
      DTPIC_CFL = LARGE_NUMBER 
      
      if(dimn.eq.2) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(dimn.eq.3) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0 
      
      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)

      DO L = 1, MAX_PIP
         DELETE_PART = .false. 
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
         pc = pc+1
         if(pea(l,4)) cycle 

         DES_LOC_DEBUG = .FALSE.

         !if(mppic) then 
         !   IJK = PIJK(L, 4) 
            
         !   COUNT_BC = DES_CELLWISE_BCDATA(IJK)%COUNT_DES_BC 
         !   IF(COUNT_BC.GE.1) CALL MPPIC_ADD_FRIC_FORCE(L)
         !ENDIF
! If a particle is classified as new, then forces are ignored. 
! Classification from new to existing is performed in routine
! des_check_new_particle.f
         IF(.NOT.PEA(L,2))THEN 
            FC(L, :) = FC(L,:)/PMASS(L) + GRAV(:)
         ELSE 
            FC(L,:) = ZERO
         ENDIF

         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001) 
         !By comparing MFIX equations and equations in those papers, 
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         DP_BAR = F_gp(L)/(PMASS(L))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO

         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO 

         M = PIJK(L,5)
         IJK = PIJK(L,4)
               
         IJK_OLD = IJK
         
         PIJK_OLD(:) = PIJK(L,:)

         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         COEFF_EN =  MPPIC_COEFF_EN 
         UPRIMETAU(:) = ZERO

         VEL_ORIG(:) = DES_VEL_NEW(L,:)

         DES_VEL_NEW(L,:) = (DES_VEL_OLD(L,:) + & 
         & FC(L,:)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         IF(L.EQ.FOCUS_PARTICLE) THEN 
            
            WRITE(*,'(A20,2x,3(2x,i4))') 'CELL ID = ', PIJK(L,1:3)
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'EPS = ', 1.d0 - EP_g(PIJK(L,4))
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'DES_VEL ORIG = ', DES_VEL_NEW(L,:)
            
            WRITE(*,'(A20,2x,3(2x,g17.8))') 'FC = ', FC(L,:)
         ENDIF

         !MEANVEL(1) = DES_U_S(IJK_OLD,M)
         !MEANVEL(2) = DES_V_S(IJK_OLD,M)
         !IF(DIMN.EQ.3) MEANVEL(3) = DES_W_S(IJK_OLD,M)
               
         PS_FORCE(:) = PS_GRAD(L, :) 
         DELUP(:) = -( DTSOLID*PS_FORCE(:))/((1.d0+DP_BAR*DTSOLID))
         DELUP(:) = DELUP(:)/ROP_S(IJK_OLD,M)

         MEANVEL(:) = AVGSOLVEL_P(L,:)

         DO IDIM = 1, DIMN
            IF(PS_FORCE(IDIM).LE.ZERO) THEN 
               UPRIMETAU(IDIM) = MIN(DELUP(IDIM), (1+COEFF_EN)*(MEANVEL(IDIM)-DES_VEL_NEW(L,IDIM)))
               UPRIMETAU_INT(IDIM) = UPRIMETAU(IDIM)
               UPRIMETAU(IDIM) = MAX(UPRIMETAU(IDIM), ZERO) 
            ELSE 
               UPRIMETAU(IDIM) = MAX(DELUP(IDIM), (1+COEFF_EN)*(MEANVEL(IDIM)-DES_VEL_NEW(L,IDIM)))
               UPRIMETAU_INT(IDIM) = UPRIMETAU(IDIM)
               UPRIMETAU(IDIM) = MIN(UPRIMETAU(IDIM), ZERO) 
            END IF 

         ENDDO 

         DES_VEL_NEW(L,:) = DES_VEL_NEW(L,:) + UPRIMETAU(:)
         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
         DES_VEL_NEW(L,:)*DTSOLID 
               
         
         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(L,1:DIMN), DES_VEL_NEW(L, 1:DIMN)))
               
         RAD_EFF = DES_RADIUS(L) 
               !RAD_EFF = (DES_STAT_WT(L)**(1.d0/3.d0))*DES_RADIUS(L)
         MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
         MEAN_FREE_PATH = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2
               
         !IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then 
         !   DES_VEL_NEW(L,:) = (DES_VEL_NEW(L,:)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         !ENDIF

          D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
                              
          DIST = SQRT(DES_DOTPRDCT(D,D))

          IF(DIST.GT.MEAN_FREE_PATH) THEN 
          !WRITE(*,*) 'UPRIME OLD =  ', UPRIMETAU(:)
          !WRITE(*,*) 'DIST GT MEAN FREE PATH= ', DIST, MEAN_FREE_PATH

             DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
             DES_VEL_NEW(L,:)*DTSOLID*MEAN_FREE_PATH/DIST

             D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
             
             DIST = SQRT(DES_DOTPRDCT(D,D))
             !WRITE(*,*) 'new moved distance  = ', dist,1. -  ep_g(ijk)
          ENDIF


         CALL MPPIC_FIND_NEW_CELL(L)

         IJK = PIJK(L,4)
         
         IF(EP_G(IJK).LT.ep_star.and.fluid_at(ijk)) THEN !.and.(ijk.ne.ijk_old)) then 
            DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
            DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
            
            
            !DES_POS_NEW(L,:) = DES_POS_OLD(L,:)+ rand_vel(L, :)*dtsolid 
            
         ENDIF
         
         
         PIJK(L,:) = PIJK_OLD(:)
               

         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         D_GRIDUNITS(1) = ABS(D(1))/DX(PIJK(L,1))
         D_GRIDUNITS(2) = ABS(D(2))/DY(PIJK(L,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DIMN.EQ.3) D_GRIDUNITS(3) = ABS(D(3))/DZ(PIJK(L,3))

         DIST = SQRT(DES_DOTPRDCT(D,D))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER 
         IF(DIMN.EQ.3) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)
            
            
! Check if the particle has moved a distance greater than or equal to grid spacing 
! if so, then exit 
            
         DO IDIM = 1, DIMN 
            IF(D_GRIDUNITS(IDIM).GT.ONE) THEN 
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN 
                  
                  WRITE(UNIT_LOG, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  WRITE(*, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  DELETE_PART = .true. 
                     
               ENDIF
               !CALL mfix_exit(myPE)
            END IF
               
         END DO
         IF(.not.DELETE_PART) DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         FC(L,:) = ZERO

         IF(DELETE_PART) THEN 
            PEA(L,1) = .false. 
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO

      
      IF(MPPIC) THEN 
         CALL global_all_max(DTPIC_CFL)
         PIP = PIP - PIP_DEL_COUNT
         
         LPIP_DEL_COUNT_ALL(:) = 0
         LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT 
         CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL) 
         IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN 
            IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
                 & FREQUENTLY: MONITOR THIS MESSAGE'

            WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
                 & FREQUENTLY: MONITOR THIS MESSAGE'
            !DO IPROC = 0, NUMPES-1 
            !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
            !ENDDO
            
         ENDIF

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         IF(DTSOLID.GT.DTPIC_MAX) THEN 
            !IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2004) DTSOLID
         ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN 

            !IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
         ELSE
            !WRITE(*,'(A40,2x,g17.8)') 'DT
            !IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
         END IF


      ENDIF
 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')  

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)         

2001  FORMAT(/1X,70('*'),//,10X,  & 
           & 'MOVEMENT UNDESIRED IN CFNEWVALUES: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, & 
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, & 
           & 'DES_VEL_NEW = ',  3(2x, g17.8)) 

 2004 FORMAT(/10x, & 
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)
      
 2002 FORMAT(/10x, & 
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

 2003 FORMAT(/10x, & 
      & 'DTSOLID REMAINS UNCHANGED AT = ', g17.8)


      RETURN
      END SUBROUTINE CFNEWVALUES_MPPIC_SINDER 

      
      SUBROUTINE CFNEWVALUES_MPPIC

      USE param
      USE param1
      USE parallel
      USE physprop
      USE constant
      USE fldvar
      USE discretelement
      USE des_bc
      USE mpi_utility
      USE mppic_wallbc
      USE mfix_pic 
      USE randomno
      USE cutcell
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER L, M, IDIM 
      INTEGER I, J, K, IJK, IJK_OLD, IJK2, IJKE, IJKW, IJKN, IJKS, IJKT, IJKB

      DOUBLE PRECISION D(DIMN), DIST, &
                       NEIGHBOR_SEARCH_DIST, DP_BAR, COEFF_EN, MEANVEL(DIMN), D_GRIDUNITS(3)

      DOUBLE PRECISION DELUP(DIMN), UPRIMETAU(DIMN), UPRIMETAU_INT(DIMN), MEAN_FREE_PATH, PS_FORCE(DIMN), VEL_ORIG(DIMN)
! index to track accounted for particles  
      INTEGER PC 

! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

! maximum distance particles can move in MPPIC 
      DOUBLE PRECISION MAXDIST_PIC, UPRIMEMOD, UPRIMEMODNEW, signvel

! dt's in each direction  based on cfl_pic for the mppic case 
      
      DOUBLE PRECISION DTPIC_TMPX, DTPIC_TMPY , DTPIC_TMPZ, THREEINTOSQRT2, RAD_EFF, MEANUS(DIMN, MMAX), POS_Z
      DOUBLE PRECISION :: DPS_DXE, DPS_DXW, DPS_DYN, DPS_DYS, DPS_DZT, DPS_DZB
      DOUBLE PRECISION :: XI_EAST, XI_WEST, XI_NORTH, XI_SOUTH, XI_TOP, XI_BOTTOM, epg_min2, velf_part(dimn)
      
      LOGICAL :: DELETE_PART, INSIDE_DOMAIN
      INTEGER :: PIP_DEL_COUNT, count_bc 

      DOUBLE PRECISION MEANUS_e(DIMN, MMAX), MEANUS_w(DIMN, MMAX),MEANUS_n(DIMN, MMAX),MEANUS_s(DIMN, MMAX),MEANUS_t(DIMN, MMAX), MEANUS_b(DIMN, MMAX)      
      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1), PIJK_OLD(5), epg_min_loc(1)
      
      double precision  sig_u, mean_u,ymid 
      double precision, allocatable, dimension(:,:) ::  rand_vel

      double precision :: norm1, norm2, norm3 
!-----------------------------------------------
! Functions 
!-----------------------------------------------
      DOUBLE PRECISION, EXTERNAL :: DES_DOTPRDCT

!-----------------------------------------------      
      
      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

      PC = 1
      FOCUS_PARTICLE = -1
      DTPIC_CFL = LARGE_NUMBER 
      
      if(dimn.eq.2) THREEINTOSQRT2 = 2.D0*SQRT(2.D0)
      if(dimn.eq.3) THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      THREEINTOSQRT2 = 3.D0*SQRT(2.D0)
      DES_VEL_MAX(:) = ZERO
      PIP_DEL_COUNT = 0 
      
      !EPG_MIN2 = MINVAL(EP_G(:))
      !epg_min_loc = MINLOC(EP_G)
      !IJK = epg_min_loc(1)
      allocate(rand_vel(MAX_PIP, DIMN))
      do idim = 1, dimn 
         mean_u = zero
         sig_u = 1.d0 
         CALL NOR_RNO(RAND_VEL(1:MAX_PIP, IDIM), MEAN_U, SIG_U)
      enddo

      DO L = 1, MAX_PIP
         DELETE_PART = .false. 
! pradeep skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
         pc = pc+1
         if(pea(l,4)) cycle 

         DES_LOC_DEBUG = .FALSE.

         IF(.NOT.PEA(L,2))THEN 
            FC(L, :) = FC(L,:)/PMASS(L) + GRAV(:)
         ELSE 
            FC(L,:) = ZERO
         ENDIF
         
         !DP_BAR is the D_p in the Snider paper (JCP, vol 170, 523-549, 2001) 
         !By comparing the MFIX and equations in those papers, 
         !D_p = Beta/(EP_S*RHOP)
         !F_gp in drag_fgs.f  = Beta*PVOL/EP_S
         !Therefore, D_p = F_gp/(PVOL*RHOP) = F_gp/PMASS
         IF(DES_ONEWAY_COUPLED) F_gp(L) = ZERO 
         DP_BAR = F_gp(L)/(PMASS(L))
         IF(.NOT.DES_CONTINUUM_COUPLED) DP_BAR = ZERO

         
         if(.not.MPPIC_PDRAG_IMPLICIT) DP_BAR = ZERO 
         M = PIJK(L,5)
         IJK = PIJK(L,4)
               
         IJK_OLD = IJK
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         VEL_ORIG(:) = DES_VEL_NEW(L,:)

         DES_VEL_NEW(L,:) = (DES_VEL_OLD(L,:) + & 
         & FC(L,:)*DTSOLID)/(1.d0+DP_BAR*DTSOLID)

         MEANVEL(1) = des_u_s(ijk,m)
         MEANVEL(2) = des_v_s(ijk,m)
         IF(DIMN.EQ.3) MEANVEL(3) = DES_W_S(IJK,M)

         
                                !if(mod(L,400).eq.0) write(*,'(A,2x,10(2x,g17.8))') 'rand vel = ', rand_vel(1,:), sig_u
         do idim = 1, dimn 
            SIG_U = 0.005D0*MEANVEL(IDIM) 
            !DES_VEL_NEW(L, idim) = DES_VEL_NEW(L, idim) + sig_u*rand_vel(pc, idim ) 
            
            SIG_U = 0.005D0!*MEANVEL(IDIM) 
            rand_vel(L, idim)  = sig_u*DES_VEL_NEW(L, IDIM)*rand_vel(L, idim) 
            !DES_VEL_NEW(L, idim) = DES_VEL_NEW(L, idim) + rand_vel(pc, idim) 
         enddo
         
      ENDDO

      IF(.not.DES_ONEWAY_COUPLED) CALL MPPIC_APPLY_PS_GRAD

      
      PC = 1 
      DO L = 1, MAX_PIP
         DELETE_PART = .false. 
! pradeep: skip ghost particles
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
         pc = pc+1
         if(pea(l,4)) cycle 

         DES_LOC_DEBUG = .FALSE.

               
         M = PIJK(L,5)
         
         IJK = PIJK(L,4)
         IJK_OLD = IJK 
         PIJK_OLD(:) = PIJK(L,:)


         UPRIMEMOD = SQRT(DOT_PRODUCT(DES_VEL_NEW(L,1:DIMN), DES_VEL_NEW(L, 1:DIMN)))
               
         RAD_EFF = DES_RADIUS(L) 
         !RAD_EFF = (DES_STAT_WT(L)**(1.d0/3.d0))*DES_RADIUS(L)
         MEAN_FREE_PATH  = MAX(1.d0/(1.d0-EP_STAR), 1.d0/(1.D0-EP_G(IJK_OLD)))
         MEAN_FREE_PATH = MEAN_FREE_PATH*RAD_EFF/THREEINTOSQRT2
               
         IF(UPRIMEMOD*DTSOLID.GT.MEAN_FREE_PATH) then 
            DES_VEL_NEW(L,:) = (DES_VEL_NEW(L,:)/UPRIMEMOD)*MEAN_FREE_PATH/DTSOLID
         ENDIF

         DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
         DES_VEL_NEW(L,:)*DTSOLID + rand_vel(L, :)*dtsolid 
         
         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         
         IJK_OLD = PIJK(L,4)
         
         CALL MPPIC_FIND_NEW_CELL(L)

         
         IJK = PIJK(L,4)
         
         !IF((EP_G(IJK).LT.0.35.and.fluid_at(ijk)).or.(ep_g(ijk_old).lt.0.35)) then !.and.(ijk.ne.ijk_old)) then 
         !IF((EP_G(IJK).LT.EP_STAR.and.fluid_at(ijk)).and.(ijk.ne.ijk_old)) then 
         INSIDE_DOMAIN = .true. 

         INSIDE_DOMAIN = FLUID_AT(IJK)

         IF(CUT_CELL_AT(IJK)) THEN 
            POS_Z = zero 
            IF(DIMN.eq.3) POS_Z = DES_POS_NEW(L,3)
            CALL GET_DEL_H_DES(IJK,'SCALAR', & 
            & DES_POS_NEW(L,1),  DES_POS_NEW(L,2), &
            & POS_Z, & 
            & DIST, NORM1, NORM2, NORM3, .true.)
            
            IF(DIST.LE.ZERO) INSIDE_DOMAIN = .false. 
         ENDIF

         !IF((EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN)) then 
         !IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN) then  
         IF(EP_G(IJK).LT.EP_STAR.and.INSIDE_DOMAIN.and.IJK.NE.IJK_OLD) then 
            
            IF(CUT_CELL_AT(IJK)) then 
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
               DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + rand_vel(L, :)*dtsolid 
               DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
            ELSE
               !IF(IJK.NE.IJK_OLD) THEN 
                  DES_POS_NEW(L,:) = DES_POS_OLD(L,:)
                  DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + rand_vel(L, :)*dtsolid 
                  DES_VEL_NEW(L,:) = 0.8d0*DES_VEL_NEW(L,:)
               !ENDIF
            ENDIF
         ENDIF
         
         
         PIJK(L,:) = PIJK_OLD(:)
         DIST = SQRT(DES_DOTPRDCT(D,D))

          IF(DIST.GT.MEAN_FREE_PATH) THEN 
          !WRITE(*,*) 'UPRIME OLD = ', UPRIMETAU(:)
          !WRITE(*,*) 'DIST GT MEAN FREE PATH= ', DIST, MEAN_FREE_PATH

             !DES_POS_NEW(L,:) = DES_POS_OLD(L,:) + &
             !DES_VEL_NEW(L,:)*DTSOLID*MEAN_FREE_PATH/DIST

             !D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
             
                                !DIST = SQRT(DES_DOTPRDCT(D,D))
             !WRITE(*,*) 'new moved distance  = ', dist,1. -  ep_g(ijk)
          ENDIF



         D(:) = DES_POS_NEW(L,:) - DES_POS_OLD(L,:)
         D_GRIDUNITS(1) = ABS(D(1))/DX(PIJK(L,1))
         D_GRIDUNITS(2) = ABS(D(2))/DY(PIJK(L,2))
         D_GRIDUNITS(3) = 1.d0
         IF(DIMN.EQ.3) D_GRIDUNITS(3) = ABS(D(3))/DZ(PIJK(L,3))

         DIST = SQRT(DES_DOTPRDCT(D,D))
         DTPIC_TMPX = (CFL_PIC*DX(PIJK(L,1)))/(ABS(DES_VEL_NEW(L,1))+SMALL_NUMBER)
         DTPIC_TMPY = (CFL_PIC*DY(PIJK(L,2)))/(ABS(DES_VEL_NEW(L,2))+SMALL_NUMBER)
         DTPIC_TMPZ = LARGE_NUMBER 
         IF(DIMN.EQ.3) DTPIC_TMPZ = (CFL_PIC*DZ(PIJK(L,3)))/(ABS(DES_VEL_NEW(L,3))+SMALL_NUMBER)
            
            
! Check if the particle has moved a distance greater than or equal to grid spacing 
! if so, then exit 
            
         DO IDIM = 1, DIMN 
            IF(D_GRIDUNITS(IDIM).GT.ONE) THEN 
               IF(DMP_LOG.OR.myPE.eq.pe_IO) THEN 
                  
                  WRITE(UNIT_LOG, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'rand_vel = ', rand_vel(L,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_old = ', des_pos_old(l,:)
                  WRITE(UNIT_LOG, '(A,2x,3(g17.8))') 'des_pos_new = ', des_pos_new(L,:)
                  
                  WRITE(*, 2001) L, D_GRIDUNITS(:), DES_VEL_NEW(L,:)
                  DELETE_PART = .true. 
                     
               ENDIF
               !CALL mfix_exit(myPE)
            END IF
               
         END DO
         IF(.not.DELETE_PART) DTPIC_CFL = MIN(DTPIC_CFL, DTPIC_TMPX, DTPIC_TMPY, DTPIC_TMPZ)
         FC(L,:) = ZERO

         IF(DELETE_PART) THEN 
            PEA(L,1) = .false. 
            PIP_DEL_COUNT = PIP_DEL_COUNT + 1
         ENDIF
         IF (DES_LOC_DEBUG) WRITE(*,1001)
      ENDDO


      DEALLOCATE(RAND_VEL)      

      IF(MPPIC) THEN 
         CALL global_all_max(DTPIC_CFL)
         PIP = PIP - PIP_DEL_COUNT
         
         LPIP_DEL_COUNT_ALL(:) = 0
         LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT 
         CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL) 
         IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN 
            IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS STEPPING MORE THAN ONE GRID SPACE = ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
                 & FREQUENTLY: MONITOR THIS MESSAGE'

            WRITE(UNIT_LOG,'(/,2x,A,2x,i10,/,A)') 'TOTAL NUMBER OF PARTS  STEPPING MORE THAN ONE GRID SPACEC= ', SUM(LPIP_DEL_COUNT_ALL(:)), 'THIS SHOULD NOT HAPPEN &
                 & FREQUENTLY: MONITOR THIS MESSAGE'
            !DO IPROC = 0, NUMPES-1 
            !   WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
            !ENDDO
            
         ENDIF

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         IF(DTSOLID.GT.DTPIC_MAX) THEN 
            !IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2004) DTSOLID
         ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN 

            !IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
         ELSE
            !WRITE(*,'(A40,2x,g17.8)') 'DT
            !IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
         END IF


      ENDIF
 1000 FORMAT(3X,'---------- FROM CFNEWVALUES ---------->')
 1001 FORMAT(3X,'<---------- END CFNEWVALUES ----------')  

 1002 FORMAT(/1X,70('*')//&
         ' From: CFNEWVALUES -',/&
         ' Message: Particle ',I10, ' moved a distance ', ES17.9, &
         ' during a',/10X, 'single solids time step, which is ',&
         ' greater than',/10X,'its radius: ', ES17.9)
 1003 FORMAT(1X,70('*')/)         

2001  FORMAT(/1X,70('*'),//,10X,  & 
           & 'MOVEMENT UNDESIRED IN CFNEWVALUES: PARTICLE', i5, /,10X, &
           & ' MOVED MORE THAN A GRID SPACING IN ONE TIME STEP', /,10X, & 
           & 'MOVEMENT IN GRID UNITS = ', 3(2x, g17.8),/,10X,  &
           & 'TERMINAL ERROR: NOT STOPPING, BUT DELETING THE PARTICLE', &
           & /1X,70('*'), /10X, & 
           & 'DES_VEL_NEW = ',  3(2x, g17.8)) 

 2004 FORMAT(/10x, & 
      & 'DTSOLID SHUD BE REDUCED TO', g17.8)
      
 2002 FORMAT(/10x, & 
      & 'DTSOLID CAN BE INCREASED TO', g17.8)

 2003 FORMAT(/10x, & 
      & 'DTSOLID REMAINS UNCHANGED AT = ', g17.8)


      RETURN
      END SUBROUTINE CFNEWVALUES_MPPIC
    
      
!------------------------------------------------------------------------
! subroutine       : des_dbgpic 
! Author           : Pradeep G.
! Purpose          : For debugging the pic values 
! Parameters       : pstart - start indices of the particle 
!                    pend - end indices of the particle 
!                    pfreq - optional frequency (when the local count matches the 
!                    frequency the filw will be written)
!                    if not send then it prints the file 
!------------------------------------------------------------------------
      subroutine des_dbgpic (pstart,pend,pfreq)
      use discretelement 
      USE fldvar
      implicit none 
! dummy variables 
      integer pstart,pend
      integer,optional :: pfreq 
! local variables 
      integer lp,lijk
      integer, save :: lfcount = 0 ,lfreq =0 
      character(30) :: filename 
      if (present(pfreq)) then 
         lfreq = lfreq+1 
         if (lfreq .ne. pfreq) return
         lfreq =0 
      end if 
      lfcount = lfcount + 1 
      write(filename,'("debug",I3.3)'), lfcount
      open (unit=100,file=filename)
      do lp = pstart,pend 
         if (pea(lp,1) .and. .not.pea(lp,4)) then 
            lijk = pijk(lp,4) 
            write(100,*),"positon =",lijk,pijk(lp,1),pijk(lp,2),pijk(lp,3),ep_g(lijk),DES_U_s(lijk,1)
            write(100,*),"forces =", FC(lp,2),tow(lp,1)
         end if 
      end do  
      close (100)
      return 
      end 

!------------------------------------------------------------------------
! subroutine       : des_dbgtecplot 
! Author           : Pradeep G.
! Purpose          : prints the tecplot file for particle location 
! Parameters       : pstart - start indices of the particle 
!                    pend - end indices of the particle 

!                    pfreq - optional frequency (when the local count matches the 
!                    frequency the filw will be written)
!                    if not send then it prints the file 
!------------------------------------------------------------------------
      subroutine des_dbgtecplot (pstart,pend,pfreq)
      use discretelement 
      USE fldvar
      implicit none 
! dummy variables 
      integer pstart,pend
      integer,optional :: pfreq 
! local variables 
      integer lp,lijk
      integer, save :: lfcount = 0 ,lfreq =0 
      character(30) :: filename 
      if (present(pfreq)) then 
         lfreq = lfreq+1 
         if (lfreq .ne. pfreq) return
         lfreq =0 
      end if 
      lfcount = lfcount + 1 
      write(filename,'("new_tec",I3.3,".dat")'), lfcount
      open (unit=100,file=filename)
      write(100,*) 'VARIABLES = "ijk"   "x"   "y"   "vx"   "vy"   "ep_g", "FCX" ,"FCY", "TOW"'
      write(100,'(A,F14.7,A)') 'zone T = "' , s_time , '"' 
      do lp = pstart,pend 
         if (pea(lp,1)) then 
            lijk = pijk(lp,4) 
            write(100,*),lijk,des_pos_new(lp,1),des_pos_new(lp,2),des_vel_new(lp,1),des_vel_new(lp,2),ep_g(lijk),fc(lp,1),fc(lp,2),tow(lp,1)
         end if 
      end do  
      close (100)
      return 
      end 


      !
!  Module name: DES_PERIODIC_BC
!
!  Purpose: At this point the particles have already been advanced in
!     time according to the force balance.  Now, update the particle 
!     position/velocity according to simple periodic boundary conditions.
!     Only particles that have crossed a periodic BC will have their
!     positon modified accordingly. 
!
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE DES_PERIODIC_BC(L)

!-----------------------------------------------
! Modules
!-----------------------------------------------
     
      USE geometry
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! given particle ID number
      INTEGER, INTENT (IN) :: L

! local variables for system dimensions
      DOUBLE PRECISION LX, LXE, LXW, LY, LYN, LYS, LZ, LZT, LZB

! local variables for x, y, z position of the particle
      DOUBLE PRECISION XPOS, YPOS, ZPOS 

! local variables for i, j, k indices associated with particle position
      INTEGER I_INDEX, J_INDEX, K_INDEX, I, J, K      
!-----------------------------------------------
! Functions 
!-----------------------------------------------

!-----------------------------------------------

! assign temporary local variables for quick reference
      LXE = EX2
      LXW = WX1
      LX = LXE - LXW
      LYN = TY2
      LYS = BY1
      LY = LYN - LYS
      LZT = NZ2
      LZB = SZ1
      LZ = LZT - LZB

! initialize local variables 
      XPOS = DES_POS_NEW(L,1)
      YPOS = DES_POS_NEW(L,2)
      I_INDEX = PIJK(L,1)
      J_INDEX = PIJK(L,2)
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(L,3)
         K_INDEX = PIJK(L,3)
      ENDIF


      IF(DES_PERIODIC_WALLS_X) THEN
         IF(XPOS.GE.LXE) THEN
            XPOS = XPOS - LX
            I_INDEX = IMIN1
         ELSEIF(XPOS.LT.LXW) THEN
            XPOS = XPOS + LX
            I_INDEX = IMAX1
         ENDIF
      ENDIF
      IF(DES_PERIODIC_WALLS_Y) THEN
         IF(YPOS.GE.LYN) THEN
            YPOS = YPOS - LY
            J_INDEX = JMIN1
         ELSEIF(YPOS.LT.LYS) THEN
            YPOS = YPOS + LY
            J_INDEX = JMAX1
         ENDIF
      ENDIF
      IF(DIMN.EQ.3 .AND. DES_PERIODIC_WALLS_Z) THEN
         IF(ZPOS.GE.LZT) THEN
            ZPOS = ZPOS - LZ
            K_INDEX = KMIN1
         ELSEIF(ZPOS.LT.LZB) THEN
            ZPOS = ZPOS + LZ
            K_INDEX = KMAX1
         ENDIF
      ENDIF

! set particle position and index according to periodicity      
      DES_POS_NEW(L,1) = XPOS
      DES_POS_NEW(L,2) = YPOS
      PIJK(L,1) = I_INDEX
      PIJK(L,2) = J_INDEX
      IF (DIMN.EQ.3) THEN
         DES_POS_NEW(L,3) = ZPOS
         PIJK(L,3) = K_INDEX
      ENDIF

      RETURN

      END SUBROUTINE DES_PERIODIC_BC



!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!  Module name: DES_LEES_EDWARDS_BC
!
!  Purpose: At this point the particles have already been advanced in
!     time according to the force balance.  Now, update the particle 
!     position/velocity according to Lees & Edwards boundary conditions.
!     Only particles that have crossed a LE BC will have their positon/
!     velocity modified accordingly. 
!
!
!  Author: Janine Galvin                              Date: 
!  Reviewer:                                          Date: 
!
!  Comments: For further details recommend Computer Simulation of
!     Liquids, by M.P. Allen and D. J. Tildesley (section 8.2, 
!     pages 246-247). 
!
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      SUBROUTINE DES_LEES_EDWARDS_BC(L)

!-----------------------------------------------
! Modules
!-----------------------------------------------
     
      USE geometry
      USE discretelement

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! given particle ID number
      INTEGER, INTENT (IN) :: L

! x, y, z position of the particle
      DOUBLE PRECISION XPOS, YPOS, ZPOS

! x, y, z velocity of the particle      
      DOUBLE PRECISION XVEL, YVEL, ZVEL

! local variables for system dimensions
      DOUBLE PRECISION LX, LY, LZ, LXE, LXW, LYN, LYS, LZT, LZB
      
! i, j, k indices associated with particle position
      INTEGER I_INDEX, J_INDEX, K_INDEX, I, J, K

! the value of the i, j, k index minus two (i.e., I=I-2) or 1 (I=1)
! whichever is greater.  this will prevent errors in referencing the
! particle position by indices outside the array bounds.  For example,
! it is possible that a particle is in a cell with i = 2 to move to a
! cell with i = 4, wherein XE(I-2) may be used causing the code to fail.
      INTEGER IMINUS2, JMINUS2, KMINUS2

! local variable for shear direction
      CHARACTER*4 SHEAR_DIR      

! local variable for relative velocity of shear
      DOUBLE PRECISION REL_VEL

! quantity to indicate direction of particle movement resulting from
! crossing the LE BC
!   0 = no move (i.e. particle did not cross LE BC)      
!  -1 = west/south/bottom move
!   1 = east/north/top move
      INTEGER MOVE_DIR

! the distance between the periodic boundaries corresponding to the
! direction the shear is acting. for du/dy shear this corresponds to the
! x domain length
      DOUBLE PRECISION DOMAIN_SIZE

! determined by first calculating the distance the LE boundary (cell) 
! that was originally aligned with the center cell traveled in a given
! time step.  then integer multiples of the domain size are subtracted
! from this quantity until a distance less than the domain size remains
      DOUBLE PRECISION OFFSET_DISTANCE

! logical for local debug warnings
      LOGICAL DES_LOC_DEBUG

!-----------------------------------------------
! Functions 
!-----------------------------------------------

!-----------------------------------------------      

! Given that a particle's movement is currently restrained to be less 
! than its own radius during a single time step, the following
! manipulations only need to be conducted once in a given step.  If a
! particle could travel more than the domain length in a single time
! step then the 'IF' statements should be replaced with 'WHILE' loops.

! For example, in the case of DUDY shear the looping would need to be 
! performed until the new y position is smaller than ylength and/or 
! the distance between the periodic boundaries corresponding to the
! direction the shear is acting. For du/dy shear this corresponds to the
! x domain length greater than 0.  In such a case, looping would also be
! needed for adjusting the x position and additional measures would be 
! needed for locating the particle's i and j indices. Such calculations
! are not accounted for in the current code.

! If the particle does not cross the LE boundary nor the periodic
! boundary its indices will be correctly determined using the existing
! method found in particles_in_cell.f

! assign temporary local variables for quick reference
      LXE = EX2
      LXW = WX1
      LX = LXE - LXW
      LYN = TY2
      LYS = BY1
      LY = LYN - LYS
      LZT = NZ2
      LZB = SZ1
      LZ = LZT - LZB

      REL_VEL = DES_LE_REL_VEL
      SHEAR_DIR = TRIM(DES_LE_SHEAR_DIR)

! assign temporary local variables for manipulation/use
      XPOS = DES_POS_NEW(L,1)
      YPOS = DES_POS_NEW(L,2)
      XVEL = DES_VEL_NEW(L,1)
      YVEL = DES_VEL_NEW(L,2)
      I_INDEX = PIJK(L,1)
      J_INDEX = PIJK(L,2)      
      IF (DIMN .EQ. 3) THEN
         ZPOS = DES_POS_NEW(L,3)
         ZVEL = DES_VEL_NEW(L,3)
         K_INDEX = PIJK(L,3)
      ENDIF 


! initialize local quantities
      OFFSET_DISTANCE = 0
      IF (REL_VEL > 0 ) THEN
         MOVE_DIR = 1   ! particle is moved east, north or up
      ELSEIF (REL_VEL < 0 ) THEN
         MOVE_DIR = -1  ! particle is moved west, south, or down
      ELSE
         MOVE_DIR = 0
      ENDIF


      IF (DIMN .EQ. 2) THEN

! 2D shear : du/dy
! ----------------------------------------               
         IF(SHEAR_DIR.EQ.'DUDY') THEN
            DOMAIN_SIZE = LX
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF
                        
            IF (YPOS >= LYN) THEN
! particle crossed north Lees & Edwards boundary                    
               YPOS = YPOS - LY 
               J_INDEX = JMIN1
               XPOS = XPOS - OFFSET_DISTANCE
               XVEL = XVEL - REL_VEL
               MOVE_DIR = -1*MOVE_DIR
            ELSEIF (YPOS < LYS) THEN
! particle crossed south Lees & Edwards boundary            
               YPOS= YPOS +  LY
               J_INDEX = JMAX1
               XPOS = XPOS + OFFSET_DISTANCE
               XVEL = XVEL + REL_VEL
            ELSE
! particle did not cross the Lees & Edwards boundary                    
               MOVE_DIR = 0
            ENDIF

            IF (XPOS >= LXE) THEN
! particle crossed east periodic boundary 
               XPOS = XPOS - DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN 
! particle did not cross the LE boundary so the value of PIJK can
! readily be assigned as a particle is constrained to move less than 
! its own radius during a single time step
                  I_INDEX = IMIN1
               ELSE
! particle crossed the LE boundary and depending on the shear rate and
! particle position it is possible for the particle to have been moved
! more than its own radius and outside the domain.  the particle was
! moved east (move_dir = 1). 
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
! in this case the particle must have been on the east side of the
! domain for it to have crossed the east periodic boundary and the move
! will shift it to the west side of the domain. so to minimize the
! search start at the west most domain index and increment up toward the
! index corresponding to the old x position
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF

            ELSEIF (XPOS < LXW) THEN
! particle crossed west periodic boundary 
               XPOS = XPOS + DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN
! particle did not cross the LE boundary
                  I_INDEX = IMAX1
               ELSE
! particle crossed the LE boundary. the particle was moved left 
! (move_dir=-1).  
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
! in this case the particle must have been on the west side of the
! domain for it to have crossed the west periodic boundary and the move
! will shift it to the east side of the domain
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF

            ELSEIF (MOVE_DIR .NE. 0) THEN
! particle did not cross the periodic boundary but did cross the LE
! boundary (i.e. particle was moved but not across domain). identify the
! index by first searching if particle is in the 2 outer most cells then
! searching the the cells neighboring the old position

               I = I_INDEX   ! for shorthand/quick reference
               IMINUS2 = I - 2
               IF (IMINUS2 < 1) IMINUS2 = 1

               IF(XPOS >= XE(1) .AND. XPOS < XE(IMIN1)) THEN
                  I_INDEX = IMIN1
               ELSEIF (XPOS >= XE(IMAX1-1) .AND. XPOS < XE(IMAX1)) THEN
                  I_INDEX = IMAX1
               ELSEIF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN 
                  I_INDEX = I
               ELSEIF(XPOS >= XE(I) .AND. XPOS < XE(I+1)) THEN 
                  I_INDEX = I+1
               ELSEIF(XPOS >= XE(IMINUS2) .AND. XPOS < XE(I-1)) THEN
                     I_INDEX = I-1
               ELSEIF (MOVE_DIR .EQ. -1) THEN
! particle was moved west
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
! in this case the particle must have been on the east side of the domain 
! and the move will shift it to the west of the domain
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ELSEIF (MOVE_DIR .EQ. 1) THEN
! particle was moved east
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
! in this case the particle must have been on the LHS of the domain and
! the move will shift it to the RHS of the domain
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF  ! end final block to identify index

            ENDIF   ! the particle did not cross a LE or periodic boundary 


! 2D shear : dv/dx
! see 2D DUDY shear section for more details on code
! ---------------------------------------- 
         ELSEIF(SHEAR_DIR.EQ.'DVDX') THEN
            DOMAIN_SIZE = LY            
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF
          
            IF (XPOS >= LXE) THEN
! particle crossed east Lees & Edwards boundary
               XPOS = XPOS - LX 
               I_INDEX = IMIN1
               YPOS = YPOS - OFFSET_DISTANCE
               YVEL = YVEL - REL_VEL
               MOVE_DIR = -1*MOVE_DIR
            ELSEIF (XPOS < LXW) THEN
! particle crossed west Lees & Edwards boundary
               XPOS = XPOS + LX 
               I_INDEX = IMAX1
               YPOS = YPOS + OFFSET_DISTANCE
               YVEL = YVEL + REL_VEL
            ENDIF

            IF (YPOS >= LYN) THEN
! particle crossed north periodic boundary
               YPOS = YPOS - DOMAIN_SIZE
               IF (MOVE_DIR .EQ. 0) THEN
! particle did not cross the LE boundary
                  J_INDEX = JMIN1
               ELSE
! the particle crossed the LE boundary
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMIN1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMIN1,YPOS,&
                        YN,'y','j')
                  ENDIF
               ENDIF

            ELSEIF (YPOS < LYS) THEN
! particle crossed south periodic boundary
               YPOS = YPOS + DOMAIN_SIZE
               IF (MOVE_DIR .EQ. 0) THEN
! particle did not cross the LE boundary
                  J_INDEX = JMAX1
               ELSE
! particle crossed the LE boundary
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMAX1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMAX1,YPOS,&
                        YN,'y','j')
                  ENDIF
               ENDIF

            ELSEIF (MOVE_DIR .NE. 0) THEN
! particle did not cross the periodic boundary but did cross the LE
! boundary
               J = J_INDEX   ! for shorthand/quick reference
               JMINUS2 = J - 2
               IF (JMINUS2 < 1) JMINUS2 = 1

               IF(YPOS >= YN(1) .AND. YPOS < YN(JMIN1)) THEN
                  J_INDEX = JMIN1
               ELSEIF(YPOS >= YN(JMAX1-1) .AND. YPOS < YN(JMAX1)) THEN
                  J_INDEX = JMAX1
               ELSEIF(YPOS >= YN(J-1) .AND. YPOS < YN(J)) THEN
                  J_INDEX = J
               ELSEIF(YPOS >= YN(J) .AND. YPOS < YN(J+1)) THEN
                  J_INDEX = J+1
               ELSEIF(YPOS >= YN(JMINUS2) .AND. YPOS < YN(J-1)) THEN
                  J_INDEX = J-1
               ELSEIF (MOVE_DIR .EQ. -1) THEN
! particle was moved south
                  IF(DABS(OFFSET_DISTANCE) >= 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMIN1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMIN1,YPOS,&
                        YN,'y','j')
                  ENDIF
               ELSEIF(MOVE_DIR .EQ. 1) THEN
! the particle was moved north
                  IF (DABS(OFFSET_DISTANCE) >= 0.5d0*DOMAIN_SIZE) THEN
                     J_INDEX = DES_GETINDEXFROMPOS(JMAX1,J_INDEX,YPOS,&
                        YN,'y','j')
                  ELSE
                     J_INDEX = DES_GETINDEXFROMPOS(J_INDEX,JMAX1,YPOS,&
                        YN,'y','j')
                  ENDIF                 
               ENDIF   ! end final block to identify index

            ENDIF   ! the particle did not cross a LE or periodic boundary

         ENDIF   ! endif des_le_shear_dir == dudy or dvdx
      ENDIF  ! if dimn == 2


      IF (DIMN .EQ. 3) THEN
! 3D shear : du/dy
! ---------------------------------------- 
         IF(TRIM(DES_LE_SHEAR_DIR).EQ.'DUDY') THEN
            DOMAIN_SIZE = LX
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
               XPOS = XPOS - OFFSET_DISTANCE
               XVEL = XVEL - REL_VEL
               MOVE_DIR = -1*MOVE_DIR
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
               XPOS = XPOS + OFFSET_DISTANCE
               XVEL = XVEL + REL_VEL
            ELSE
               MOVE_DIR = 0
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN 
                  I_INDEX = IMIN1
               ELSE
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + DOMAIN_SIZE 
               IF (MOVE_DIR .EQ. 0) THEN
                  I_INDEX = IMAX1
               ELSE
                  IF (DABS(OFFSET_DISTANCE) < 0.5d0*DOMAIN_SIZE) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF
            ELSEIF (MOVE_DIR .NE. 0) THEN
               I = I_INDEX   ! for shorthand/quick reference
               IMINUS2 = I - 2
               IF (IMINUS2 < 1) IMINUS2 = 1

               IF(XPOS >= XE(1) .AND. XPOS < XE(IMIN1)) THEN
                  I_INDEX = IMIN1
               ELSEIF (XPOS >= XE(IMAX1-1) .AND. XPOS < XE(IMAX1)) THEN
                  I_INDEX = IMAX1
               ELSEIF(XPOS >= XE(I-1) .AND. XPOS < XE(I)) THEN 
                  I_INDEX = I
               ELSEIF(XPOS >= XE(I) .AND. XPOS < XE(I+1)) THEN 
                  I_INDEX = I+1
               ELSEIF(XPOS >= XE(IMINUS2) .AND. XPOS < XE(I-1)) THEN
                     I_INDEX = I-1
               ELSEIF (MOVE_DIR .EQ. -1) THEN
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMIN1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMIN1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ELSEIF (MOVE_DIR .EQ. 1) THEN
                  IF (DABS(OFFSET_DISTANCE) >= (0.5d0*DOMAIN_SIZE)) THEN
                     I_INDEX = DES_GETINDEXFROMPOS(IMAX1,I_INDEX,XPOS,&
                        XE,'x','i')
                  ELSE
                     I_INDEX = DES_GETINDEXFROMPOS(I_INDEX,IMAX1,XPOS,&
                        XE,'x','i')
                  ENDIF
               ENDIF  ! end final block to identify i index

            ENDIF   ! the particle did not cross a LE or periodic boundary 

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
            ENDIF

! 3D shear : du/dz
! ----------------------------------------             
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DUDZ') THEN
            DOMAIN_SIZE = LX
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
               XPOS = XPOS - OFFSET_DISTANCE
               XVEL = XVEL - REL_VEL
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
               XPOS = XPOS + OFFSET_DISTANCE
               XVEL = XVEL + REL_VEL
            ENDIF

! work out code to find i index 
            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
            ENDIF
            
            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
            ENDIF

! 3D shear : dv/dx
! ---------------------------------------- 
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DVDX') THEN
            DOMAIN_SIZE = LY
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
               I_INDEX = IMIN1
               YPOS = YPOS - OFFSET_DISTANCE
               YVEL = YVEL - REL_VEL
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
               I_INDEX = IMAX1
               YPOS = YPOS + OFFSET_DISTANCE
               YVEL = YVEL + REL_VEL
            ENDIF

! work out code to find j index            
            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
            ENDIF

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
            ENDIF

! 3D shear : dv/dz
! ----------------------------------------             
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DVDZ') THEN
            DOMAIN_SIZE = LY
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
               K_INDEX = KMIN1
               YPOS = YPOS - OFFSET_DISTANCE
               YVEL = YVEL - REL_VEL
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
               K_INDEX = KMAX1
               YPOS = YPOS + OFFSET_DISTANCE
               YVEL = YVEL + REL_VEL
            ENDIF

! work out code to find j index                   
            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
               I_INDEX = IMIN1
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
               I_INDEX = IMAX1
            ENDIF


! 3D shear : dw/dx
! ---------------------------------------- 
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DWDX') THEN
            DOMAIN_SIZE = LZ
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - (EX2 - LXW)
               I_INDEX = IMIN1
               ZPOS = ZPOS - OFFSET_DISTANCE
               ZVEL = ZVEL - REL_VEL
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + (EX2 - LXW)
               I_INDEX = IMAX1
               ZPOS = ZPOS + OFFSET_DISTANCE
               ZVEL = ZVEL + REL_VEL
            ENDIF

! work out code to find k index            
            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
            ENDIF

            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
            ENDIF


! 3D shear : dw/dy
! ----------------------------------------             
         ELSEIF(TRIM(DES_LE_SHEAR_DIR).EQ.'DWDY') THEN
            DOMAIN_SIZE = LZ
            IF (REL_VEL .NE. ZERO) THEN
               OFFSET_DISTANCE = REL_VEL*S_TIME - DOMAIN_SIZE *&
                  DBLE( FLOOR( REAL( (REL_VEL*S_TIME/DOMAIN_SIZE) )) )
            ENDIF

            IF (YPOS >= LYN) THEN
               YPOS = YPOS - LY
               J_INDEX = JMIN1
               ZPOS = ZPOS - OFFSET_DISTANCE
               ZVEL = ZVEL - REL_VEL
            ELSEIF (YPOS < LYS) THEN
               YPOS = YPOS + LY
               J_INDEX = JMAX1
               ZPOS = ZPOS + OFFSET_DISTANCE
               ZVEL = ZVEL + REL_VEL
            ENDIF

! work out code to find k index
            IF (ZPOS >= LZT) THEN
               ZPOS = ZPOS - LZ
            ELSEIF (ZPOS < LZB) THEN
               ZPOS = ZPOS + LZ
            ENDIF

            IF (XPOS >= LXE) THEN
               XPOS = XPOS - LX
               I_INDEX = IMIN1
            ELSEIF (XPOS < LXW) THEN
               XPOS = XPOS + LX
               I_INDEX = IMAX1
            ENDIF


         ENDIF   ! endif des_le_shear_dir == dudy, dudz, dvdx, dvdz, dwdx or dwdz

      ENDIF   ! endif dimn == 3

! set particle position and index according to periodicity/lees & edwards BC
      DES_POS_NEW(L,1) = XPOS  
      DES_POS_NEW(L,2) = YPOS
      DES_VEL_NEW(L,1) = XVEL
      DES_VEL_NEW(L,2) = YVEL
      PIJK(L,1) = I_INDEX
      PIJK(L,2) = J_INDEX
      IF (DIMN .EQ. 3) THEN
         DES_POS_NEW(L,3) = ZPOS  
         DES_VEL_NEW(L,3) = ZVEL
         PIJK(L,3) = K_INDEX
      ENDIF

      RETURN

 1002 FORMAT(/1X,70('*')//,' From: DES_DES_LEES_EDWARDS_BC',/,&
         ' Message: Calculation for the ', A, ' index associated with ',&
         'the particles new',/10X, A, '-position= ', ES12.5,&
         ' where the old ',A, '-position= ',/10X, ES12.5 ' and old ',&
         A, ' index= ', I5, /1X,70('*')/)



      END SUBROUTINE DES_LEES_EDWARDS_BC
      
