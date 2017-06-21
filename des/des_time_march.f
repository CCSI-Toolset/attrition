!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: DES_TIME_MARCH                                         C
!
!     Purpose: Called in model/time_march.f to do DES calcs
!     Main DEM driver routine
!
!                                                                         C
!     Author: Jay Boyalakuntla                           Date: 21-Jun-04  C
!     Reviewer: Sreekanth Pannala                        Date: 09-Nov-06  C
!     Reviewer: Rahul Garg                               Date: 01-Aug-07  C
!     Comments: Changed the calling rules to neighbor search routines     C
!     Pradeep: do_nsearch has to be set for calling neighbour. Since this C 
!              this flag is used during exchange it has to be set before  C
!              calling particle in cell
!
!     Reviewer: Wesley Xu, Dave Decroix                  Date: 11-May-12
!     Revision: Added the initialization of desradiusnew of attrition model
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE DES_TIME_MARCH
     
!------------------------------------------------
! Modules
!------------------------------------------------
      USE param 
      USE param1 
      USE run
      USE output
      USE physprop
      USE fldvar
      USE geometry
      USE pgcor
      USE pscor
      USE cont
      USE coeff
      USE tau_g
      USE tau_s
      USE visc_g
      USE visc_s
      USE funits 
      USE vshear
      USE scalars
      USE drag
      USE rxns
      USE compar     
      USE time_cpu  
      USE discretelement   
      USE constant
      USE sendrecv
      USE des_bc
      USE cutcell 
      USE mppic_wallbc
      USE mfix_pic
      Use des_thermo
      Use des_rxns
      Use interpolation
      
      IMPLICIT NONE
!------------------------------------------------
! Local variables
!------------------------------------------------
      INTEGER NN, FACTOR, NP, IJK, I, J, K, BCV_I, LL, PC

!     Local variables to keep track of time when dem restart and des
!     write data need to be written when des_continuum_coupled is F
      DOUBLE PRECISION DES_RES_TIME, DES_SPX_TIME

!     Temporary variables when des_continuum_coupled is T to track
!     changes in solid time step 
      DOUBLE PRECISION TMP_DTS, DTSOLID_TMP 
      CHARACTER*5 FILENAME

!     Temporary variable used to track reporting frequency for the
!     maximum overlap and maximum no. neighbors for given NN loop
      DOUBLE PRECISION DES_TMP_TIME

!     Logical to see whether this is the first entry to this routine
      LOGICAL,SAVE:: FIRST_PASS = .TRUE.

!     temp
      DOUBLE PRECISION xpos,ypos, zpos, NORM_CF(3), DIST

      INTEGER INTERP_IJK(2**DIMN)

      DOUBLE PRECISION INTERP_WEIGHTS(2**DIMN)
      
      DOUBLE PRECISION :: DES_KE_VEC(DIMN)
! Identifies that the indicated particle is of interest for debugging
      LOGICAL FOCUS
      
      
!!	  integer cnt1,cnt2

!!$      double precision omp_start, omp_end
!!$      double precision omp_des_start,omp_des_end
!!$      double precision omp_get_wtime 
!!$      double precision omp_des_now,omp_des_last        
      
!-----------------------------------------------

      INCLUDE 'function.inc'
      INCLUDE 'fun_avg1.inc'
      INCLUDE 'fun_avg2.inc'

!!$      omp_des_start=omp_get_wtime()
! initialization      
      TMP_DTS = ZERO
      IF(FIRST_PASS.AND..NOT.MPPIC) THEN 

      if(dmp_log)write(unit_log,'(1X,A)')&
         '---------- FIRST PASS DES_TIME_MARCH ---------->'
         S_TIME = ZERO
         FIRST_PASS = .FALSE.
         INQC = INIT_QUAD_COUNT

! When no particles are present, skip the startup routines that loop over particles.
! That is, account for a setup that starts with no particles in the system.
         IF(PARTICLES /= 0) THEN
            if(do_nsearch) CALL NEIGHBOUR         

! Initialize attrition radius array desradiusnew, not suitable for MPPIC
         Do NP=1, MAX_PIP
            DESRadiusNew(NP)=DES_RADIUS(NP)
         ENDDO
         
! COHESION INITIALIZE ! for VDW this is now done in make_arrays_des.f
            IF(USE_COHESION .AND. .NOT.VAN_DER_WAALS)THEN
               CALL INITIALIZE_COHESION_PARAMETERS
               CALL INITIALIZE_COH_INT_SEARCH
            ENDIF

! To do only des in the 1st time step only for a new run so the particles settle down
! before the coupling is turned on

            IF(RUN_TYPE == 'NEW') THEN
               IF(DES_CONTINUUM_COUPLED.AND.(.NOT.USE_COHESION)) THEN
                  DES_CONTINUUM_COUPLED = .FALSE.
                  DO FACTOR = 1, NFACTOR
                     IF (FACTOR .EQ. 1) then 
                        if(dmp_log)write(unit_log,'(3X,A,/,5X,A,X,I10,X,A)') &
                        'FIRST PASS in DES_TIME_MARCH for new runs',&
                        'DEM settling period performed', NFACTOR, &
                        'times'
                      end if 
                     ! Force calculation         
                     CALL CALC_FORCE_DES

                     CALL CFNEWVALUES

! pradeep : set the flag do_nsearch before calling particle in cell
                     IF(MOD(FACTOR,NEIGHBOR_SEARCH_N).EQ.0) do_nsearch = .true. 
                     CALL PARTICLES_IN_CELL
                     if(do_nsearch) CALL NEIGHBOUR

                     ! Neighbor search                      
!                     IF(MOD(FACTOR,NEIGHBOR_SEARCH_N).EQ.0) THEN 
!                        CALL NEIGHBOUR
!                     ELSEIF(DO_NSEARCH) THEN 
!                        CALL NEIGHBOUR
!                        DO_NSEARCH = .FALSE.
!                     ENDIF
                  ENDDO
                  DES_CONTINUUM_COUPLED = .TRUE.
                  if(dmp_log)write(unit_log,'(3X,A)') 'END DEM settling period'
               ENDIF   ! end if coupled and no cohesion
!*************************************************************
! Pradeep following lines are not necessary 
! Setting flags calc_fc and callfromdes is not necessary
! also call to particle in cell is not necessary --- discussed over email 
               IF(DES_INTERP_ON) THEN 
                  CALC_FC = .FALSE.
                  CALLFROMDES = .FALSE.
               ENDIF
               CALL PARTICLES_IN_CELL
!*************************************************************
!Pradeep Remove ****** the comment 
!               CALL WRITE_DES_DATA
!**********************************************
               if(dmp_log)write(unit_log,'(3X,A,X,ES15.5)') &
                  'DES data file written at time =', S_TIME
            ENDIF   ! end if on new run type

            if(dmp_log)write(unit_log,'(1X,A)')&
               '<---------- END FIRST PASS DES_TIME_MARCH ----------'

         ENDIF   ! end if particles /= 0         
      ENDIF    ! end if first pass


! In case of restarts assign S_TIME from MFIX TIME 
      S_TIME = TIME
! pradeep following is added  to avoid error with check compile flag 
      DTSOLID_TMP = 0.0  
      IF(DES_CONTINUUM_COUPLED) THEN
         IF(DT.GE.DTSOLID) THEN
            FACTOR = CEILING(real(DT/DTSOLID))
         ELSE
            FACTOR = 1
            DTSOLID_TMP = DTSOLID
            DTSOLID = DT
         ENDIF
         
 1999    FORMAT(/1X,70('- '),//,10X,  & 
         & 'DEM SIMULATION CALLED ', 2x, i5, 2x, 'times this fluid step', /10x &
         & 'S_TIME, DT, DTSOLID and PIP = ', 3(2x,g17.8),2x, i10)
         
         if(dmp_log) WRITE(UNIT_LOG, 1999) factor, s_time, dt, dtsolid, pip

         if(mype.eq.pe_IO) WRITE(*, 1999) factor, s_time, dt, dtsolid, pip
         
      ELSE
         FACTOR = CEILING(real((TSTOP-TIME)/DTSOLID)) 
         if(myPE.eq.pe_IO) write(*,'(1X,A)')&
            '---------- START DES_TIME_MARCH ---------->'
         if(myPE.eq.pe_IO)write(*,'(3X,A,X,I10,X,A)') &
            'DEM SIMULATION will be called', FACTOR, 'times'
! Initialization for des_spx_time, des_res_time         
         IF(RUN_TYPE .EQ. 'NEW') THEN
            DES_SPX_TIME =  S_TIME
            DES_RES_TIME =  S_TIME
         ELSE
            DES_SPX_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) +&
               1 ) * DES_SPX_DT
            DES_RES_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_RES_DT) +&
               1 ) * DES_RES_DT
         ENDIF
      ENDIF   ! end if/else (des_continuum_coupled)

      IF (DES_CONTINUUM_COUPLED) DES_SPX_DT = SPX_DT(1)
      IF (RUN_TYPE .EQ. 'NEW') THEN
         DES_TMP_TIME = S_TIME
      ELSE
         DES_TMP_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) +1 ) *&
            DES_SPX_DT
      ENDIF

      IF(NEIGHBOR_SEARCH_N.GT.FACTOR) THEN 
         !NEIGHBOR_SEARCH_N = FACTOR
      ENDIF 

      IF(MPPIC) THEN 
         !compute the gas-phase pressure gradient at the beginning of the 
         !des loop as the gas-phase pressure field will not change during 
         !des calls
         IF(DES_CONTINUUM_COUPLED)   CALL COMPUTE_PG_GRAD
      ELSE
         IF(DES_CONTINUUM_COUPLED)   CALL COMPUTE_PG_GRAD
      ENDIF 

! JFD : Update radius 
         Do NP=1, MAX_PIP
            DESRadiusNew(NP)=DES_RADIUS(NP)
         ENDDO

!!$    omp_start = omp_get_wtime()
      DO NN = 1, FACTOR 
!!$      omp_des_now=omp_get_wtime()
!!        call system_clock( count=cnt1 )

         IF(DES_CONTINUUM_COUPLED) THEN
            IF(S_TIME.GE.(TIME+DT)) EXIT
! If the current time in the discrete loop exceeds the current time in
! the continuum simulation, exit the discrete loop
            IF((S_TIME+DTSOLID).GT.(TIME+DT)) THEN 
! If next time step in the discrete loop will exceed the current time 
! in the continuum simulation, modify the discrete time step so final
! time will match 
               TMP_DTS = DTSOLID
               DTSOLID = TIME + DT - S_TIME
            ENDIF 
            IF(DEBUG_DES) THEN
               !rAHUL: Pradeep, dont send this output to the log file
               !the log file becomes just to big. I think the idea
               !here is to proint this on screen 
               if(mype.eq.pe_IO) write(*,'(3X,A,X,I10,X,A,X,ES15.7)') &
                  'DES-COUPLED LOOP NO. =', NN, ' S_TIME =', S_TIME 
                if(mype.eq.pe_IO) write(*,'(3X,A,X,ES15.7)') &
                  'DTSOLID =', DTSOLID
            ENDIF
            

            CALC_FC = .TRUE.
            IF(DES_INTERP_ON .AND. NN.EQ.FACTOR) THEN 
! Toggle flag so the mean fields are calculated only at the last DEM
! time step
               CALLFROMDES = .FALSE.
            ELSE 
               CALLFROMDES = .TRUE.
            ENDIF
         ELSE   ! else if (des_continuum_coupled)
            IF(DEBUG_DES)then 
               if(mype.eq.pe_IO)write(*,'(3X,A,X,I10,X,A,X,ES15.7)') &
               'DEM LOOP NO. =', NN, ' S_TIME =', S_TIME 
             end if 
             
         ENDIF   ! end if/else (des_continuum_coupled) 
         
         
! If system is empty, skip force calculation calls
! pradeep removed the PIP check as communication between processor has 
! to be taken place all the time; further if pip =0 no time is spent
! in force computation 
!         IF (PIP /= 0) THEN
            IF(MPPIC) THEN 
               CALL MPPIC_COMPUTE_PS_GRAD            
               IF(DES_CONTINUUM_COUPLED)   CALL DRAG_FGS
               CALL CFUPDATEOLD
            ELSE 
!!$      omp_start=omp_get_wtime()
            	    CALL CALC_FORCE_DES
               !cfupdateold is called from inside clac_force_des
               !Update particle position, velocity
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'Calc_force_des:',omp_end - omp_start               
            ENDIF
    
!---------------------------------------------------------------------------->>>
! Loop over all particles
            PC = 1
            DO NP = 1, MAX_PIP
               IF(PC .GT. PIP) EXIT
               IF(.NOT.PEA(NP,1)) CYCLE
! Reset the debug flag
               FOCUS = .FALSE.
! Set the debugging flag
               IF(DEBUG_DES .AND. NP.EQ.FOCUS_PARTICLE) FOCUS = .TRUE.

! Calculate time dependent physical properties
               CALL DES_PHYSICAL_PROP(NP, FOCUS)

! Calculate cell-center interpolation weights and determine the 
! associated IJK values for the cells accounting for boundary
! conditions.
               IF(DES_INTERP_ON .AND. &
                  (ANY_DES_SPECIES_EQ .OR. DES_CONV_EQ)) THEN
                  INTERP_IJK(:) = -1
                  INTERP_WEIGHTS(:) = ZERO
                  CALL INTERPOLATE_CC(NP, INTERP_IJK, INTERP_WEIGHTS, &
                     FOCUS)
               ENDIF
! Calculate thermodynamic energy exchange
               IF(DES_ENERGY_EQ) CALL CALC_THERMO_DES(NP, &
                  INTERP_IJK, INTERP_WEIGHTS, FOCUS)

! Calculate reaction rates and interphase mass transfer
               IF(ANY_DES_SPECIES_EQ) CALL DES_RRATES(NP, &
                  INTERP_IJK, INTERP_WEIGHTS, FOCUS, 'SOLIDS')
! Increment the particle counter
               PC = PC + 1
            ENDDO
!----------------------------------------------------------------------------<<<    
            
!!$      omp_start=omp_get_wtime()
            CALL CFNEWVALUES
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'Cfnewvalues:',omp_end - omp_start    
    

!---------------------------------------------------------------------------->>>
! Loop over all particles
            PC = 1
            DO NP = 1, MAX_PIP
               IF(PC .GT. PIP) EXIT
               IF(.NOT.PEA(NP,1)) CYCLE
! Below line is added by Tingwe Not sure if it is necessary               
!               if(pea(np,4))cycle 
 
! Reset the debug flag
               FOCUS = .FALSE.
! Set the debugging flag
               IF(DEBUG_DES .AND. NP.EQ.FOCUS_PARTICLE) FOCUS = .TRUE.
! Update particle temperature
               IF(DES_ENERGY_EQ) &
                  CALL DES_THERMO_NEWVALUES(NP, FOCUS)
! Update particle from reactive chemistry process.
               IF(DES_SPECIES_EQ(PIJK(NP,5))) &
                  CALL DES_REACTION_MODEL(NP, FOCUS)
! Increment the particle counter
               PC = PC + 1
            ENDDO
!----------------------------------------------------------------------------<<<


	
! Impose the wall-particle boundary condition for mp-pic case 
            IF(MPPIC) CALL MPPIC_APPLY_WALLBC 
	 


! For systems with inlets/outlets check to determine if a particle has
! fully entered or exited the domain.  If the former, remove the status
! of 'new' and if the latter, remove the particle.
            IF (DES_MIO) CALL DES_CHECK_PARTICLE
! pradeep: set do_nsearch before calling particle_in_Cell
            IF(NN.EQ.1 .OR. MOD(NN,NEIGHBOR_SEARCH_N).EQ.0) do_nsearch =.true. 
!!$      omp_start=omp_get_wtime()
            CALL PARTICLES_IN_CELL
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'Particles_in_cell:',omp_end - omp_start    
            
            !IF(MPPIC.and.(.not.MPPIC_SOLID_STRESS_SNIDER)) THEN 
            !   CALL MPPIC_COMPUTE_PS_GRAD            
            !   CALL MPPIC_APPLY_PS_GRAD
            !endif
            !IF(MPPIC) CALL MPPIC_APPLY_PS_GRAD

!!$      omp_start=omp_get_wtime()            
            if (do_nsearch.AND.(.NOT.MPPIC)) then 
               IF(DEBUG_DES) then 
                  if(dmp_log)write(unit_log,'(3X,A,I10,/,5X,A,I10)') &
                     'Calling NEIGHBOUR: during iteration NN =', NN
               end if 
               call neighbour 
            end if 
!!$      omp_end=omp_get_wtime()
!!$      write(*,*)'Neighbor_search:',omp_end - omp_start	

!            IF(NN.EQ.1 .OR. MOD(NN,NEIGHBOR_SEARCH_N).EQ.0) THEN 
!               IF(DEBUG_DES) if(dmp_log)write(unit_log,'(3X,A,A,/,5X,A,I10)') &
!                  'Calling NEIGHBOUR: NN=1 or ',&
!                  'MOD(NN,NEIGHBOR_SEARCH_N)=0',&
!                  'NEIGHBOR_SEARCH_N = ', NEIGHBOR_SEARCH_N
!               CALL NEIGHBOUR
!            ELSEIF(DO_NSEARCH) THEN 
!               IF(DEBUG_DES) if(dmp_log)write(unit_log,'(3X,A,A,/,5X,A,A,L)') &
!                  'Calling NEIGHBOUR: a particle moved ',&
!                  'more than its radius since', 'last time NEIGHBOUR ',&
!                  'was called; DO_NSEARCH = ', DO_NSEARCH
!               CALL NEIGHBOUR
!               DO_NSEARCH = .FALSE.
!            ENDIF

!         ENDIF   ! end if particles /= 0


! Update time to reflect changes 
         S_TIME = S_TIME + DTSOLID

! When coupled the granular temperature subroutine is only calculated at end 
! of the current DEM simulation 
         IF(DES_CONTINUUM_COUPLED .AND. NN.EQ.FACTOR) &
            CALL DES_GRANULAR_TEMPERATURE

! When coupled, all write calls are made in time_march (the continuum 
! portion) according to user settings for spx_time and res_time.
! The following section targets data writes for DEM only cases:
         IF(.NOT.DES_CONTINUUM_COUPLED) THEN    
! Keep track of TIME for DEM simulations
            TIME = S_TIME
 
! Write data using des_spx_time and des_res_time; note the time will
! reflect current position of particles  
            IF(PRINT_DES_DATA) THEN
               IF ( (S_TIME+0.1d0*DTSOLID >= DES_SPX_TIME) .OR. &
                    (S_TIME+0.1d0*DTSOLID >= TSTOP) .OR. &
                    (NN == FACTOR) ) THEN
                  DES_SPX_TIME = &
                     ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) &
                     + 1 )*DES_SPX_DT
! Granular temperature subroutine should be called/calculated when
! writing DES data 
                  CALL DES_GRANULAR_TEMPERATURE
!pradeep remove ***************
                  CALL WRITE_DES_DATA
!pradeep remove ***********
                  if(dmp_log)write(unit_log,'(3X,A,X,ES15.5)') &
                     'DES data file written at time =', S_TIME
               ENDIF
            ENDIF

            IF ( (S_TIME+0.1d0*DTSOLID >= DES_RES_TIME) .OR. &
                 (S_TIME+0.1d0*DTSOLID >= TSTOP) .OR. &
                 (NN == FACTOR) ) THEN
               DES_RES_TIME = &
                  ( INT((S_TIME+0.1d0*DTSOLID)/DES_RES_DT) &
                  + 1 )*DES_RES_DT
! Pradeep new subroutine for writing des_restart file 
!                  CALL WRITE_DES_RESTART
                  call des_write_restart
! Write RES1 here since it won't be called in time_march.  This will
! also keep track of TIME
               CALL WRITE_RES1 
               if(dmp_log)write(unit_log,'(3X,A,X,ES15.5)') &
               'DES.RES and .RES files written at time =', S_TIME
            ENDIF
         ENDIF  ! end if (.not.des_continuum_coupled)


! J.Musser : mass inlet/outlet -> particles entering the system
         IF(DES_MI)THEN 
            DO BCV_I = 1, DES_BCMI
               IF(PI_FACTOR(BCV_I) .GT. 1)THEN
                  IF(DES_MI_TIME(BCV_I) .LE. S_TIME) THEN   !Verify start time
                     CALL DES_MASS_INLET(BCV_I)
                     DES_MI_TIME(BCV_I) = S_TIME + PI_FACTOR(BCV_I) *&
                        DTSOLID
                  ENDIF
               ELSE
                  CALL DES_MASS_INLET(BCV_I)
               ENDIF
            ENDDO
         ENDIF

! Report some statistics on overlap and neighbors to screen log
         IF ( (S_TIME+0.1d0*DTSOLID >= DES_TMP_TIME) .OR. &
              ( (S_TIME+0.1d0*DTSOLID >= TSTOP) .AND. &
               (.NOT.DES_CONTINUUM_COUPLED) ) .OR. &          
              (NN .EQ. FACTOR) ) THEN
            DES_TMP_TIME = ( INT((S_TIME+0.1d0*DTSOLID)/DES_SPX_DT) &
               + 1 )*DES_SPX_DT
               !if(dmp_log)write(unit_log,'(3X,A,I10,A,/,5X,A,X,I5,2X,A,X,ES15.7)') &
               !'For loop ', NN, ' :',&
               !'MAX no. of neighbors =',NEIGH_MAX,&
               !'and MAX % overlap =', OVERLAP_MAX
         ENDIF

!!         call system_clock( count=cnt2 )
!!         write(*,*)'One_time_march:',cnt2-cnt1,'microseconds'	 
!!$      omp_des_last=omp_get_wtime()
!!$      write(*,*)'<--------DES_TIME_MARCH-',omp_des_last-omp_des_now,'SECONDS -------->'
      ENDDO     ! end do NN = 1, FACTOR
!!$      omp_end = omp_get_wtime()
!!$      write(*,*)'Total wall time:  ', omp_end-omp_start, '  Seconds' 

! When coupled, and if needed, reset the discrete time step accordingly
      IF(DT.LT.DTSOLID_TMP) THEN
         DTSOLID = DTSOLID_TMP
      ENDIF

      IF(TMP_DTS.NE.ZERO) THEN
         DTSOLID = TMP_DTS
         TMP_DTS = ZERO
      ENDIF

      IF(MPPIC) THEN 
         IF(DMP_LOG) WRITE(UNIT_LOG, 2000) DTSOLID, DTPIC_CFL, DTPIC_TAUP, MIN(DTPIC_CFL, DTPIC_TAUP)
         IF(myPE.eq.pe_IO) WRITE(*, 2000) DTSOLID, DTPIC_CFL, DTPIC_TAUP, MIN(DTPIC_CFL, DTPIC_TAUP)

         DTPIC_MAX = MIN(DTPIC_CFL, DTPIC_TAUP)

         IF(DTSOLID.GT.DTPIC_MAX) THEN 
            IF(DMP_LOG) WRITE(UNIT_LOG, 2001) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2001) DTPIC_MAX
            DTSOLID = DTPIC_MAX
         ELSEIF(DTSOLID.LT.DTPIC_MAX) THEN 

            IF(DMP_LOG) WRITE(UNIT_LOG, 2002) DTSOLID
            IF(myPE.eq.pe_IO) WRITE(*, 2002) DTPIC_MAX
            DTSOLID = DTPIC_MAX
         ELSE
            !WRITE(*,'(A40,2x,g17.8)') 'DT
            IF(DMP_LOG) WRITE(UNIT_LOG, 2003) DTSOLID 
            IF(mype.eq.pe_IO) WRITE(*,2003) DTSOLID
         END IF
      END IF

 2000 FORMAT(/10x, & 
      & 'ADAPTING DTSOLID FOR MPPIC', /10x, &
      & 'DTSOLID CURRENT  = ', g17.8, /10x,  &
      & 'DTPIC_CFL :  = ', g17.8, /10x,  &
      & 'DTPIC TAUP:  = ', g17.8, /10x, &
      & 'DTPIC_MAX :  = ', g17.8)

 2001 FORMAT(/10x, & 
      & 'REDUCING CURRENT DTSOLID TO', g17.8)
      
 2002 FORMAT(/10x, & 
      & 'INCREASING CURRENT DTSOLID TO', g17.8)

 2003 FORMAT(/10x, & 
      & 'DTSOLID REMAINS UNCHANGED AT = ', g17.8)
      
      IF(.NOT.DES_CONTINUUM_COUPLED)then 
         if(dmp_log)write(unit_log,'(1X,A)')&
         '<---------- END DES_TIME_MARCH ----------'
!Pradeep call send recv for variables 
      else 
         call send_recv(ep_g,2)
         call send_recv(rop_g,2)
         call send_recv(des_u_s,2)
         call send_recv(des_v_s,2) 
         if(dimn.eq.3) call send_recv(des_w_s,2) 
         call send_recv(rop_s,2)
      end if 
!!$      omp_des_end=omp_get_wtime()
!!$      write(*,*)'          DES_TIME_MARCH: ',omp_des_end - omp_des_start,' SECONDS'

      CONTAINS
      
      SUBROUTINE CALCULATE_GLOB_ENERGY 
      Implicit None 

      INTEGER PC, LL 
! squared particle velocity v.v
      DOUBLE PRECISION SQR_VEL, THETA
            DES_KE = ZERO
      DES_PE = ZERO 
      DES_VEL_AVG(:) = ZERO
      DES_KE_VEC(:) = zero 
! Calculate global average velocity, kinetic energy &
! potential energy
      PC = 1
      DO LL = 1, MAX_PIP
         IF(PC .GT. PIP) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE
         
         PC = PC+1
         IF(PEA(LL,4)) CYCLE 

         SQR_VEL = ZERO
         DO I = 1, DIMN
            SQR_VEL = SQR_VEL + DES_VEL_NEW(LL,I)**2
            DES_KE_VEC(I) = DES_KE_VEC(I) + 0.5d0*PMASS(LL)* (DES_VEL_NEW(LL,I)**2.d0)
                     
         ENDDO

         DES_KE = DES_KE +  SQR_VEL 
         DES_PE = DES_PE + PMASS(LL)*DBLE(ABS(GRAV(2)))*&
            DES_POS_NEW(LL,2)
         DES_VEL_AVG(:) =  DES_VEL_AVG(:) + DES_VEL_NEW(LL,:)

      ENDDO

      DES_KE = DES_KE/DBLE(PIP)

!J.Musser changed PARTICLES TO PIS 
      DES_VEL_AVG(:) = DES_VEL_AVG(:)/DBLE(PIP)
      
! The following quantities are primarily used for debugging/developing
! and allow a quick check of the energy conservation in the system.
! In their current form they are best applied to monodisperse cases. 
! Calculate x,y,z components of global energy & granular temperature
      GLOBAL_GRAN_ENERGY = ZERO
      GLOBAL_GRAN_TEMP  = ZERO
      DES_KE_VEC = ZERO 
      PC = 1
      DO LL = 1, MAX_PIP
         IF(PC.GT.PIP) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE 
         PC = PC+1
         IF(PEA(LL,4)) CYCLE 


         IF(PC .GT. PIP) EXIT
         IF(.NOT.PEA(LL,1)) CYCLE

         GLOBAL_GRAN_ENERGY(:) = GLOBAL_GRAN_ENERGY(:) + &
         & 0.5d0*PMASS(LL)*(DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2
         GLOBAL_GRAN_TEMP(:) = GLOBAL_GRAN_TEMP(:) + &
         &  (DES_VEL_NEW(LL,:)-DES_VEL_AVG(:))**2
         
         DES_KE_VEC(:) = DES_KE_VEC(:) + &
         & (DES_VEL_NEW(LL,:))**2

      ENDDO
      DES_KE = SUM(DES_KE_VEC(:))/DBLE(PIP)

      THETA = SUM(GLOBAL_GRAN_TEMP(:))/DBLE(PIP)
      
      WRITE(*,*) 'KE, THETA = ', DES_KE, THETA, NN !, DOT_PRODUCT(DES_VEL_AVG(:), DES_VEL_AVG(:))
      GLOBAL_GRAN_ENERGY(:) =  GLOBAL_GRAN_ENERGY(:)/DBLE(PIP)
      GLOBAL_GRAN_TEMP(:) =  GLOBAL_GRAN_TEMP(:)/DBLE(PI)
      RETURN

      END SUBROUTINE CALCULATE_GLOB_ENERGY 

      END SUBROUTINE DES_TIME_MARCH


!!$
!!$      LL = 51
!!$      WRITE(*,*) 'INIT*******************************'
!!$      WRITE(*,'(A,3(2x,g17.8))') 'DES_POS = ', DES_POS_NEW(LL,:)
!!$      WRITE(*,'(A,3(2x,g17.8))') 'DES_VEL = ', DES_VEL_NEW(LL,:)
!!$      WRITE(*,'(A,3(2x,i10))') 'DES_CELL = ', PIJK(LL,1:4)
!!$      WRITE(*,*) 'ACTIVE, FLUID, SMALL, CUT = ', PEA(LL, 1), FLUID_AT(PIJK(LL,4)), CUT_CELL_AT(PIJK(LL,4)), SMALL_CELL_AT(PIJK(LL,4))
!!$      WRITE(*,*) 'IJK = ', PIJK(LL,4), DES_CELLWISE_BCDATA(PIJK(LL,4))%COUNT_DES_BC 
!!$      WRITE(*,*) 'IMJK = ', IM_OF(PIJK(LL,4)), DES_CELLWISE_BCDATA(IM_OF(PIJK(LL,4)))%COUNT_DES_BC 
!!$      
!!$      WRITE(*,*) 'IPJK = ', IP_OF(PIJK(LL,4)), DES_CELLWISE_BCDATA(IP_OF(PIJK(LL,4)))%COUNT_DES_BC 
!!$      XPOS = DES_POS_NEW(LL,1) 
!!$      YPOS = DES_POS_NEW(LL,2)
!!$      ZPOS = ZERO 
!!$      IF (DIMN .EQ. 3) THEN
!!$         ZPOS = DES_POS_NEW(LL,3)
!!$      ENDIF
            
!!$
!!$      CALL GET_DEL_H_DES(PIJK(LL,4),'SCALAR',XPOS , YPOS, ZPOS,& 
!!$           & DIST, NORM_CF(1), NORM_CF(2), NORM_CF(3), .true.)
!!$      
!!$      WRITE(*,*) 'DIST FROM CUTCELL = ', DIST
!!$      
!!$      WRITE(*,*) 'INIT*******************************'
!!$      
