!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name:  PARTICLES_IN_CELL                                     C
!
!  Purpose: DES - Finding the fluid computational cell in which      
!           a particle lies, to calculte void fraction and also       
!           the volume averaged solids velocity of the cell            
!
!                                                                      C
!  Author: Jay Boyalakuntla                           Date: 12-Jun-04  C
!  Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C 
!  Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!  Comments: Removed the separate volume definitions and added pic     C
!            array formulation and bed height calculation.             C
!                                                                      C
!            For parallel processing indices are altered and changes   C
!            to variables related to desgridsearch are made; steps     C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE PARTICLES_IN_CELL

      USE param
      USE param1
      USE fldvar
      USE geometry
      USE indices
      USE physprop
      USE compar
      USE parallel
      USE sendrecv
      USE discretelement
      use desgrid 
      use desmpi
      USE cutcell 
      USE mfix_pic
      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
! particle no.
      INTEGER L
! accounted for particles
      INTEGER PC
! solids phase no.    
      INTEGER M
! ijk indices      
      INTEGER I, J, K, IJK, IPROC, IJK2, IJK_OLD
      INTEGER  IJKE, IJKW, IJKN, IJKS, IJKT, IJKB
! 1 over volume of fluid cell      
      DOUBLE PRECISION :: OVOL
! total volume of mth phase solids in cell and 1 over that value      
      DOUBLE PRECISION SOLVOLINC(DIMENSION_3,MMAX), OSOLVOL
! variables that count/store the number of particles in i, j, k cell
      INTEGER:: npic, pos
! pradeep,modified particle count from 3D to 1D 
      INTEGER, DIMENSION(DIMENSION_3):: particle_count
!      INTEGER, DIMENSION(DIMENSION_I,DIMENSION_J,DIMENSION_K):: particle_count
! particle x,y,z position
      DOUBLE PRECISION XPOS, YPOS, ZPOS
! Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.
! Logical for local debug warnings
      LOGICAL DES_LOC_DEBUG
! size of mesh for grid based search      
      DOUBLE PRECISION SIZEDX, SIZEDY, SIZEDZ 
! Pradeep -desgrid is moved to desgrid module 
!! variables that count/store the number of particles in i, j, k cell
!      INTEGER, DIMENSION(DESGS_IMAX2,DESGS_JMAX2,DESGS_KMAX2) :: DESGRIDSEARCH_NPIC
!      INTEGER, DIMENSION(DESGS_IMAX2,DESGS_JMAX2,DESGS_KMAX2) :: DESGS_particle_count
! Variables to calculate bed height of each solids phase
      DOUBLE PRECISION :: tmp_num(MMAX), tmp_den(MMAX), hcell, WTP, sum_eps, epg_old, VOL_AVG, epg_min, epg_min2, eps_max

! count for number of particles that were found in the ghost cell and hence
! removed       
      INTEGER :: PIP_DEL_COUNT, IMJK, ijmk, ijkm

      ! IER for error reporting
      INTEGER IER 

      double precision :: ugc, vgc 

      INTEGER :: LPIP_DEL_COUNT_ALL(0:numPEs-1)

      INTEGER :: COUNT_AVG
      INTEGER :: II, IPLUS1, IMINUS1 ! X-coordinate loop indices
      INTEGER :: JJ, JPLUS1, JMINUS1   ! Y-coordinate loop indices
      INTEGER :: KK, KPLUS1, KMINUS1   ! Z-coordinate loop indic

      INTEGER :: EPG_MIN_LOC(1), EPS_MAX_LOC(1)
      DOUBLE PRECISION :: MASS_SOL(MMAX), MASS_SOL2(MMAX)
!$      double precision omp_start1, omp_end1
!$      double precision omp_get_wtime	      
!-----------------------------------------------
      INCLUDE 'function.inc'
      INCLUDE 'ep_s1.inc'
      INCLUDE 'ep_s2.inc'


! following quantities are reset every call to particles_in_clel
      PIP_DEL_COUNT = 0 
      PINC(:) = 0
      SOLVOLINC(:,:) = ZERO
      DES_U_s(:,:) = ZERO
      DES_V_s(:,:) = ZERO
      DES_W_s(:,:) = ZERO
      DES_LOC_DEBUG = .FALSE.

      IF(FIRST_PASS) THEN
         if(dmp_log.and.debug_des) write(unit_log,'(3X,A)') &
            '---------- START FIRST PASS PARTICLES_IN_CELL ---------->'
      ENDIF

! pradeep : Call exchange particles - this will exchange particle 
! crossing boundaries as well as updates ghost particles information
!!$      omp_start1=omp_get_wtime()      
      call des_par_exchange 
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'des_par_exchange:',omp_end1 - omp_start1 
!!      PC = 1
!       by Tingwen      
!!$      omp_start1=omp_get_wtime()  
!$omp parallel do default(shared)                               &
!$omp private(l,m,xpos,ypos,zpos,i,j,k,ijk_old,ijk) schedule (guided,50)  
      DO L = 1, MAX_PIP
! pradeep - skip ghost particles
!!         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
!!         pc = pc+1
         if(pea(l,4)) cycle 
         
!!         WTP = ONE
         
!!         IF(MPPIC) WTP = DES_STAT_WT(L)
         IF(FIRST_PASS) THEN 
! Identify solid type based on diameter
!            DO M = 1, MMAX
!               IF(ABS(2.0d0*DES_RADIUS(L)-D_P0(M)).LT.SMALL_NUMBER.AND. &
!               ABS( RO_Sol(L)-RO_S(M)).LT.SMALL_NUMBER) THEN
!               PIJK(L,5) = M 
!               exit
!               ENDIF
!            ENDDO

!Only used for single solid phase attrition
	     PIJK(L,5) = 1

            IF(PIJK(L,5).EQ.0) THEN
               if(dmp_log) then 
                  write(unit_log,'(5X,A,A,I10)') &
                  'Problem determining the solids ',&
                  'association for particle: ',L
                  write(unit_log,'(7X,A,(ES15.9))') &
                  'Particle position = ', DES_POS_NEW(L,:)
                  write(unit_log,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                  'Particle diameter = ', 2.0*DES_RADIUS(L),&
                  'and D_P0(1:MMAX)= ', D_P0(1:MMAX)
                  write(unit_log,'(7X,A,ES15.9,/,7X,A,(ES15.9))')&
                  'Particle density = ', Ro_Sol(L), &
                  'and RO_s(1:MMAX) = ', RO_S(1:MMAX)
               end if 
            ENDIF
! Brute force technique to determine the particle locations in the Eulerian grid 
            XPOS = DES_POS_NEW(L,1)
            YPOS = DES_POS_NEW(L,2)
            IF (DIMN .EQ. 3) THEN
               ZPOS = DES_POS_NEW(L,3)
            ENDIF
! pradeep for parallel processing 
! in case of parallel processing the index of cell is from istart1 to iend1
! also particles can be at ghost cell for particles entering the domain 
! and exiting the domain (this might occur for restart) 
! XE will have values from istart2-1 to iend2 (see modified cfassign routine)
            do i = istart2, iend2
               if( xpos >= xe(i-1) .and. xpos < xe(i)) then
                  pijk(l,1) = i
                  exit 
               endif
            enddo
            do j = jstart2,jend2
               if(ypos >= yn(j-1) .and. ypos < yn(j)) then
                  pijk(l,2) = j
                  exit
               endif
            enddo
            if(dimn.eq.2) then
               k=1
               pijk(l,3) = 1
            else
               do k = kstart2,kend2
                  if(zpos >= zt(k-1) .and. zpos < zt(k)) then 
                     pijk(l,3) = k
                     exit
                  endif
               enddo
            endif
            IJK_OLD = funijk(i,j,k)
            
         ELSE                   ! if not first_pass
! pradeep - particle cannot jump by more than one cell distance or more
! than its radius in single time step as it will be captured in cfnewvalues
! hence first check if it is in previous location if not add +/- 1 based on location 
! the loop is modified b.c using I-2 index for particle entering which has I=1 
! might result in segmentation error
            i = pijk(l,1)
            j = pijk(l,2)
            k = pijk(l,3)
            
            IJK_OLD = pijk(l,4) 
            
            xpos = des_pos_new(l,1) 
            ypos = des_pos_new(l,2)
            if (dimn .eq. 3) then
               zpos = des_pos_new(l,3)
            endif
 
            if(xpos >= xe(i-1) .and. xpos < xe(i)) then 
               pijk(l,1) = i
            elseif(xpos >= xe(i)) then 
               pijk(l,1) = i+1
            else 
               pijk(l,1) = i-1 
            end if 

            if(ypos >= yn(j-1) .and. ypos < yn(j)) then 
               pijk(l,2) = j
            elseif(ypos >= yn(j))then 
               pijk(l,2) = j+1
            else
               pijk(l,2) = j-1
            end if 

            if(dimn.eq.2) then
               pijk(l,3) = 1
            else
               if(zpos >= zt(k-1) .and. zpos < zt(k)) then
                  pijk(l,3) = k
               elseif(zpos >= zt(k)) then 
                  pijk(l,3) = k+1
               else
                  pijk(l,3) = k-1
               end if 
            endif          
         endif   ! end of (if/else first_pass)

! sanity check for particles in the system 
         i = pijk(l,1)
         j = pijk(l,2)
         k = pijk(l,3)
         if (.not.pea(l,2) .and. .not.pea(l,3)) then 
            if (i.gt.iend1 .or. i.lt.istart1) then 
               If(MPPIC) THEN

                  !Rahul:
                  !this cud happen if the cell adjacent to the ghost cell 
                  !is a cut-cell and due to some tolerance issues, the particle
                  !is not detected outside the system. This will be a rare 
                  !occurence, but it can occur every now and then and there
                  !is no point in stalling the simulation here.
                  !Implementing an easy fix for now: delete this particle.
                  !To add more stuff later
                  !1. re distribute particle's weight among other particles 
                  !in the domain so that mass is conserved
                  !2. rather than deactivating the particle, reflect the particle
                  !inside the domain using the ghost cell bc's instead of cut-face bc 
                  IF(I.EQ.IEND1+1.AND.(XPOS >= XE(IEND1-1) .AND. XPOS <= XE(IEND1))) THEN 
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'I',I,'X',XPOS,DES_POS_OLD(L,1),'X',DES_VEL_NEW(L,1)
                     
                     PIJK(L,1) = IEND1

                  !in MPPIC a particle can lie on the surface of the wall as only the centers are tracked. 
                  else
                     write(*,*) 'piC: PIJK = ', PIJK(L,:)
                     if(dmp_log) WRITE(unit_log,1010) L,'I',I,'X',XPOS,DES_POS_OLD(L,1),'X',DES_VEL_NEW(L,1),DES_VEL_OLD(L,1), MPPIC, CARTESIAN_GRID, CUT_CELL_AT(IJK), FLUID_AT(IJK)

                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE 
                  endif
               ELSE
                  if(dmp_log) WRITE(unit_log,1007) L,'I',I,'X',XPOS,'X',DES_VEL_NEW(L,1)
                  call mfix_exit(mype)
               ENDIF
            end if 
            if (j.gt.jend1 .or. j.lt.jstart1) then 
               If(MPPIC) THEN
                  IF(J.EQ.JEND1+1.AND.(YPOS >= YN(JEND1-1) .AND. YPOS <= YN(JEND1))) THEN 
                     
                     IF(DMP_LOG) WRITE(UNIT_LOG,1011) L,'J',J,'Y',YPOS,DES_POS_OLD(L,2),'Y',DES_VEL_NEW(L,2)

                     PIJK(L,2) = JEND1
                  else
                     IF(DMP_LOG) WRITE(UNIT_LOG,1010) L,'J',J,'Y',YPOS,DES_POS_OLD(L,2),'Y',DES_VEL_NEW(L,2),DES_VEL_OLD(L,2), MPPIC, CARTESIAN_GRID, CUT_CELL_AT(IJK), FLUID_AT(IJK)

                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE
                  endif
               ELSE
                  if(dmp_log) WRITE(unit_log,1007) L,'J',J,'Y',YPOS,'Y',DES_VEL_NEW(L,2)
                  call mfix_exit(mype)
               ENDIF

            end if 

            if ((dimn.eq.3) .and. (k.gt.kend1 .or. k.lt.kstart1)) then 
               If(MPPIC) THEN
                  
                  IF(K.EQ.KEND1+1.AND.(ZPOS >= ZT(KEND1-1) .AND. ZPOS <= ZT(KEND1))) THEN 
                     
                     if(dmp_log) WRITE(unit_log,1011) L,'K',K,'Z',ZPOS, DES_POS_OLD(L,3),'Z',DES_VEL_NEW(L,3)
                     PIJK(L,3) = KEND1
                  else
                     
                     if(dmp_log) WRITE(unit_log,1010) L,'K',K,'Z',ZPOS,DES_POS_OLD(L,2),'Z',DES_VEL_NEW(L,3),DES_VEL_OLD(L,3), MPPIC, CARTESIAN_GRID, CUT_CELL_AT(IJK), FLUID_AT(IJK)
                     PIP_DEL_COUNT = PIP_DEL_COUNT + 1 
                     PEA(L,1) = .false. 
                     CYCLE 
                  ENDIF
               ELSE
                  if(dmp_log) WRITE(unit_log,1007) L,'K',K,'Z',ZPOS,'Z',DES_VEL_NEW(L,3)
                  call mfix_exit(mype)
               ENDIF
            end if 
         end if

         !set i, j, k again since pijk may have changed for the MPPIC case
!!         I = PIJK(L,1)
!!         J = PIJK(L,2)
!!         K = PIJK(L,3)

!!         IJK = FUNIJK(I,J,K)

         ! set particle in cell info and compute aggregates 
!!         PIJK(L,4) = IJK
!!         PINC(IJK) = PINC(IJK) + 1
!!         M = PIJK(L,5)
         
         !IF(M.Eq.0) WRITE(*,*) 'M EQUAL =', M, L
!!         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) +  PVOL(L)*WTP!*0.4
!!         DES_U_S(IJK,M) = DES_U_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,1)*WTP!*0.4
!!         DES_V_S(IJK,M) = DES_V_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,2)*WTP!*0.4
!!         IF(DIMN.EQ.3) DES_W_S(IJK,M) = DES_W_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,3)*WTP!*0.4

         !SOLVOLINC(IJK_OLD,M) = SOLVOLINC(IJK_OLD,M) +  PVOL(L)*WTP*0.6
         !DES_U_S(IJK_OLD,M) = DES_U_S(IJK_OLD,M) + PVOL(L)*DES_VEL_NEW(L,1)*WTP*0.6
         !DES_V_S(IJK_OLD,M) = DES_V_S(IJK_OLD,M) + PVOL(L)*DES_VEL_NEW(L,2)*WTP*0.6
         !IF(DIMN.EQ.3) DES_W_S(IJK_OLD,M) = DES_W_S(IJK_OLD,M) + PVOL(L)*DES_VEL_NEW(L,3)*WTP*0.6
      enddo      ! end loop over l = 1,particles
!$omp end parallel do 
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'first_loop:',omp_end1 - omp_start1 

!!$      omp_start1=omp_get_wtime()	  
!!$omp single private(l,wtp,i,j,k,ijk,m) !,omp_tp1,omp_tp2,omp_tp3)
!      if(first_pass)call des_par_exchange

         PC = 1
      DO L = 1, MAX_PIP
         if(pc.gt.pip) exit      
         IF(.NOT.PEA(L,1)) CYCLE
         PC = PC + 1
         if(pea(l,4)) cycle 
         WTP = ONE
         
         IF(MPPIC) WTP = DES_STAT_WT(L)
         
         I = PIJK(L,1)
         J = PIJK(L,2)
         K = PIJK(L,3)
		 
	 IJK = FUNIJK(I,J,K)
		 
         PIJK(L,4) = IJK
         PINC(IJK) = PINC(IJK) + 1
         M = PIJK(L,5)
!         pijk(l,6)=0
!         to mark if the particle is close to the walls	
!         if(i.eq.istart1 .or. i.eq.iend1 .or. j.eq.jstart1 .or. j.eq.jend1   &
!         .or. ((dimn.eq.3) .and. (k.eq.kend1 .or. k.eq.kstart1)))then
!         pijk(l,6)=1
!         endif
         
         !IF(M.Eq.0) WRITE(*,*) 'M EQUAL =', M, L
         SOLVOLINC(IJK,M) = SOLVOLINC(IJK,M) +  PVOL(L)*WTP!*0.4
         DES_U_S(IJK,M) = DES_U_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,1)*WTP!*0.4
         DES_V_S(IJK,M) = DES_V_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,2)*WTP!*0.4
         IF(DIMN.EQ.3) DES_W_S(IJK,M) = DES_W_S(IJK,M) + PVOL(L)*DES_VEL_NEW(L,3)*WTP!*0.4		 

      ENDDO      ! end loop over L = 1,particles
!!$omp end single   
!!$omp end parallel do 
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'second_loop:',omp_end1 - omp_start1 	  
      if(first_pass)call des_par_exchange

      if (FIRST_PASS)then 
          if(dmp_log.and.debug_des) write(unit_log,'(3x,a)') &
         '<---------- end first pass particles_in_cell ----------'
      end if 

      IF(MPPIC) THEN 
         PIP = PIP - PIP_DEL_COUNT

         LPIP_DEL_COUNT_ALL(:) = 0
         LPIP_DEL_COUNT_ALL(mype) = PIP_DEL_COUNT 
         CALL GLOBAL_ALL_SUM(LPIP_DEL_COUNT_ALL) 
         IF((DMP_LOG).AND.SUM(LPIP_DEL_COUNT_ALL(:)).GT.0) THEN 
            IF(PRINT_DES_SCREEN) WRITE(*,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARTICLES OUSIDE DOMAIN IN PIC = ', SUM(LPIP_DEL_COUNT_ALL(:))
            WRITE(UNIT_LOG,'(/,2x,A,2x,i10)') 'TOTAL NUMBER OF PARTICLES OUTSIDE DOMAIN IN PIC= ', SUM(LPIP_DEL_COUNT_ALL(:))
            DO IPROC = 0, NUMPES-1 
               WRITE(UNIT_LOG, '(/,A,i4,2x,A,2x,i5)') 'PARTICLES OUTSIDE DOMAIN (PIC)  ON PROC:', IPROC,' EQUAL TO', LPIP_DEL_COUNT_ALL(IPROC)
            ENDDO
            
         ENDIF
      END IF
      
! Assign/allocate the variable pic(i,j,k)%p(:). For each cell compare 
! the number of current particles in the cell to what was in the
! cell previously.  If different reallocate.  Store the particle ids
! ------------------------------------------------------------
! check all cells (including ghost cells); update entering/exiting 
! particle regions      
      MASS_SOL = ZERO 
	  
!!$      omp_start1=omp_get_wtime() 	  
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)           &
!$omp private(ijk,npic) !schedule (guided,50)     
      do ijk = ijkstart3, ijkend3
         npic = pinc(ijk)
         if (associated(pic(ijk)%p)) then
            if (npic.ne.size(pic(ijk)%p)) then
               deallocate(pic(ijk)%p)
               if (npic.gt.0) allocate(pic(ijk)%p(npic))
            endif
         else
            if (npic.gt.0) allocate(pic(ijk)%p(npic))
         endif
      enddo
!$omp end parallel do 
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'third_loop:',omp_end1 - omp_start1 

      particle_count(:) = 1
!!$      omp_start1=omp_get_wtime() 
      PC = 1
!!$omp parallel do default(shared)                               &
!!$omp private(l,ijk,pos) !schedule (guided,50)     
      DO L = 1, MAX_PIP
! pradeep skip ghost particles 
         if(pc.gt.pip) exit
         if(.not.pea(l,1)) cycle 
         pc = pc+1
         if(pea(l,4)) cycle 

         IJK = PIJK(L,4)
         pos = particle_count(IJK)
         pic(IJK)%p(pos) = L
         particle_count(IJK) = particle_count(IJK) + 1
      ENDDO
!!$omp end parallel do 
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'forth_loop:',omp_end1 - omp_start1 
	
! Calculate the cell average solids velocity, the bulk density (if not
! des_interp_on and not first_pass), the void fraction, and average
! height of each solids phase
      tmp_num(:) = ZERO 
      tmp_den(:) = ZERO 
!!$      omp_start1=omp_get_wtime() 
!$omp parallel do if(ijkend3 .ge. 2000) default(shared)        &
!$omp private(ijk,i,j,k,m,osolvol,ovol) !schedule (guided,50)     
      DO IJK = ijkstart3, ijkend3
         I = I_OF(IJK)
         J = J_OF(IJK)
         K = K_OF(IJK)
         IF(.NOT.FLUID_AT(IJK).OR..NOT.IS_ON_myPE_owns(I, J, K)) CYCLE 
         !for cut-cell, it is important to check both. FLUID_AT(IJK) 
         !alone is not enough as it is shared between procs and 
         !fluid_at(ijkend3) might be true when in fact it does 
         !not belong to that proc 
         ROP_SO(IJK,:) = ZERO 
         EP_G(IJK) = ONE   
         DO M = 1, MMAX
            IF(SOLVOLINC(IJK,M).GT.ZERO) THEN
               OSOLVOL = ONE/SOLVOLINC(IJK,M)   
               DES_U_s(IJK,M) = DES_U_s(IJK,M)*OSOLVOL
               DES_V_s(IJK,M) = DES_V_s(IJK,M)*OSOLVOL
               IF(DIMN.EQ.3) THEN
                  DES_W_s(IJK,M) = DES_W_s(IJK,M)*OSOLVOL
               ENDIF
            ENDIF
            IF(VOL(IJK).GT.0) THEN 
               OVOL = ONE/(VOL(IJK))
               IF(FIRST_PASS .OR. &
                 ((.NOT.FIRST_PASS).AND.(.NOT.DES_INTERP_ON))) THEN
                  ROP_S(IJK,M)  = RO_S(M)*SOLVOLINC(IJK,M)*OVOL
               ENDIF
            ENDIF
            IF(ROP_S(IJK,M) > ZERO) THEN
               ROP_SO(IJK,M)  = ROP_S(IJK,M) 
               MASS_SOL(M) = MASS_SOL(M) + ROP_S(IJK,M)*VOL(IJK)

               IF(.not.DES_ONEWAY_COUPLED) EP_G(IJK) = EP_G(IJK) - EP_S(IJK,M)
               IF(EP_G(IJK).LT.ZERO .AND.DES_CONTINUUM_COUPLED.and.(.not.MPPIC)) THEN 
! this does not matter if pure granular flow simulation (i.e. no fluid)
                  IF (.NOT.DES_LOC_DEBUG) THEN
                     DES_LOC_DEBUG = .TRUE.
                     if(dmp_log) write(unit_log,1000)
                  ENDIF
                  if(dmp_log) then 
                     write(unit_log,'(5X,A,I10,/,7X,A,I10,2X,I10,2X,A,ES15.9,/,7X,A,I10)') &
                     'WARNING EP_G LT zero at IJK: ', IJK,&
                     'I,J = ', I_OF(IJK), J, ' EP_S = ', EP_S(IJK,M), & 
                     'No. of particles in cell = ', PINC(IJK)
                     WRITE(unit_log, *) 'cut cell ?', cut_cell_at(ijk)
                  endif
                  
                  call mfix_exit(myPE)
               ENDIF
               ROP_G(IJK) = RO_G(IJK) * EP_G(IJK)
            ENDIF
         ENDDO   ! end loop over M=1,MMAX
      ENDDO     ! end loop over IJK=ijkstart3,ijkend3
!$omp end parallel do 
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'fifth_loop:',omp_end1 - omp_start1 

      IF(MPPIC) THEN 
         CALL MPPIC_COMPUTE_MEAN_FIELDS2
         !IF(CARTESIAN_GRID) CALL MPPIC_COMPUTE_MEAN_FIELDS_CG 
            
      end IF
      
      FIRST_PASS = .FALSE.
      IF (DES_LOC_DEBUG) then 
         if(dmp_log) write(unit_log,1001)
      end if 
      
!!$      omp_end1=omp_get_wtime()
!!$      write(*,*)'par_in_cell:',omp_end1 - omp_start1
!      open (unit=100,file='p_in_c.txt',status='unknown')
!      write(100,*)pijk,pinc,ROP_SO 
!      close(100)
      
      if(.not.mppic) then 
         RETURN 
      else
         if (nodesI*nodesJ*nodesK.gt.1) return
         !writing out some useful information on minimum epg
         !Currently, it is only valid for serial runs. RG 6/23/2011
      endif
      !IJK = funijk(11,2,kmax)
      !WRITE(*,'(A,10(2x,i5))') 'PINC in I, J, K = ', I_OF(IJK), J_OF(IJK), K_OF(IJK), PINC(IJK), PINC(IM_OF(IJK)), PINC(IP_OF(IJK))
      !WRITE(*,'(A,10(2x,g17.8))') 'EPS = ', EP_S(IJK,1), EP_S(IM_OF(IJK),1), EP_S(IP_OF(IJK),1)

!!      EPG_MIN2 = MINVAL(EP_G(:))
!!      epg_min_loc = MINLOC(EP_G(:))
!!      IJK = epg_min_loc(1)
!!      
!!      EPS_MAX = MAXVAL(ROP_S(:,1))
!!      eps_MAX = EPS_MAX/ro_s(1)
!!      eps_max_loc = MAXLOC(ROP_S(:,1))
!!      IJK2 = eps_max_loc(1)
!!      I = I_OF(IJK)
!!      J = J_OF(IJK)
!!      K = K_OF(IJK)
!!      IMJK = IM_OF(IJK)
!!      IJMK = JM_OF(IJK)
!!      IJKM = KM_OF(IJK)
!!
!!      UGC = HALF * (U_G(IJK) + U_G(IMJK))
!!      VGC = HALF * (V_G(IJK) + V_G(IJMK))
      !WRITE(*,'(10x,A,4(2x,g17.8))') 'EPGMIN1, EPS_MAX = ', epg_min2, eps_max
      
      !WRITE(*,'(10x,A,4(2x,g17.8))') 'MASS_SOL  =  ', pip*pmass(1)*des_stat_wt(1), sum(mass_sol), sum(mass_sol2)
      
      !WRITE(*,'(10x, A,4(2x,g17.8))') 'MEAN SOLID VEL = ',  DES_U_S(IJK,1),  DES_V_S(IJK,1) !,  DES_U_S(IJK,1)
      !WRITE(*,'(10x, A,4(2x,g17.8))') 'MEAN FLUID VEL = ',  UGC, VGC
      !WRITE(*,1014)  I_OF(IJK), j_of(ijk), k_of(ijk), &
      !& xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
      !& cut_cell_at(ijk),fluid_at(ijk)

      !WRITE(*,'(10x,A,4(2x,i5),2(2x,l2),/10x,A,i4)') 'EPS MAX LOCATION, CUT and FLUID CELL ? = ',IJK2, I_OF(IJK2), j_of(ijk2), k_of(ijk2), cut_cell_at(ijk2),fluid_at(ijk2)
      IF(DMP_LOG) THEN 
         !WRITE(UNIT_LOG,'(10x,A,4(2x,g17.8))') ' EPGMIN2, EPS_MAXS_SOL = ',  epg_min2
         
         !WRITE(UNIT_LOG,'(10x, A,4(2x,g17.8))') 'MEAN SOLID VEL = ',  DES_U_S(IJK,1),  DES_V_S(IJK,1) !,  DES_U_S(IJK,1)
         !WRITE(UNIT_LOG,'(10x, A,4(2x,g17.8))') 'MEAN FLUID VEL = ',  UGC, VGC
         WRITE(UNIT_LOG,1014) I_OF(IJK), j_of(ijk), k_of(ijk), &
         & xe(I) - 0.5*dx(i), yn(J)-0.5*DY(J), zt(K) - 0.5*DZ(K), & 
         & cut_cell_at(ijk),fluid_at(ijk)
      
         !WRITE(UNIT_LOG,'(10x,A,4(2x,i5),2(2x,l2),/10x,A,i4)') 'EPS MAX LOCATION, CUT and FLUID CELL ? = ',IJK2, I_OF(IJK2), j_of(ijk2), k_of(ijk2), xe(i_of(ijk)), yn(j_of(ijk)), zt(k_of(ijk)),  cut_cell_at(ijk2),fluid_at(ijk2)
      endif

      !WRITE(*,*) 'SOLID MASS: BEFORE AND AFTER: = ', mass_sol(1:mmax), mass_sol2(1:mmax)




 1000 FORMAT(3X,'---------- FROM PARTICLES_IN_CELL ---------->')
 1001 FORMAT(3X,'<---------- END PARTICLES_IN_CELL ----------') 

 1006 FORMAT(/1X,70('*')//&
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/&
         1X,70('*')/)

 1007 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/,&         
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/& 
         1X,70('*')/)

 1008 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/&
         ' Message: Check particle ',I8,' in cell ',A,/1X,A,&
         '-position: ',ES17.9,4X,A,'-velocity: ',ES17.9,/&
         1X,70('*')/)         

 1009 FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/&
         ' Message: ', A, ' index of particle ',I8,' is invalid.',/&
         ' index: ', I8,4X,A,'-position: ',ES17.9,/&
         1X,70('*')/)         



 1010     FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL -',/,&         
         ' Message: Particle ',I8,' moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,/& 
         '-velocity new and old: ',ES17.9,/& 
         ' MPPIC and Cartesian Grid ?', 2(L2,2x),/& 
         ' CUT_CELL and FLUID AT IJK ?', 2(L2,2x),/& 
         ' Marking this particle as inactive',/&          
          1X,70('*')/)

 1011     FORMAT(/1X,70('*')//&
         ' From: PARTICLES_IN_CELL: Particle recovered from ghost cell -',/,&         
         ' Message: Particle ',I8,' had moved into a',&
         ' ghost cell; from cell with ',A,' index : ',I8,/1X,A,&
         '-position new and old: ',2(ES17.9,4X),A,'-velocity: ',ES17.9,/& 
          1X,70('*')/)
            
 1012      FORMAT(/1X,70('*')//&
           & 'WARNING EP_G LT zero at IJK AFTER AVERAGING FOR CUT_CELL', /1x, &
           & 'CUT CELL ?', 1x, L2, /1x, &
           & 'IJK, I, J, K = ', 2x, i10, 3(2x,i5), /1x, &
           & 'EP_G = ', 2x, g17.8, /1x, & 
           & 'EP_S = ', 5(2x, g17.8))

 1014      FORMAT(10x,'EPG_MIN_LOCATION, I, J, K = ', 3(2x,i5),/, &
           &      10x,'XMID, YMID, ZMID FOR CELL = ', 3(2x,g17.8),/ & 
           &      10x,'CUT CELL, FLUID AT IJK ?    ', 2(2x, L2), /)

      RETURN
      END SUBROUTINE PARTICLES_IN_CELL

