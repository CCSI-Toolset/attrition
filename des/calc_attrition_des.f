!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                         C  
!     Module name: CALC_ATTRITION_DES                                     C
!
!     Purpose: Called in model/des/des_time_march to do DES attrition calcs
!
!                                                                         C
!     Author: David DeCroix, Wesley Xu                   Date: 10-Apr-12  C
!     Reviewer:                                                           C
!     Comments: This subroutine implements a particle fracture or attrition
!               model based on the following papers:
!     Modified by Wesley Xu: calculate velocity based on impact force
!
!     Revised: Wesley Xu                                 Date: 30-July-12
!     Comments: Abrasion induced attrition mechanism is included. Impact 
!               induced attrition is also modified to correct some errors.
!                                                                         C
!
!     Revised: Wesley Xu				 Date: 5-Nov-12   C
!     Comments: Model changed back to velocity based. 						
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C

      SUBROUTINE calc_attrition_des(LL,force,slide,vel_rel_n,&
      vel_rel_t,isincontact, desradius)

      USE run
      USE param1
      USE discretelement
      USE geometry
      USE compar
      USE constant
      IMPLICIT NONE

!------------------------------------------------------------------------------
!------------------  Local Variables ------------------------------------------
!------------------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION force(DIMN)
      DOUBLE PRECISION slide(DIMN)
      DOUBLE PRECISION vel_rel_n
      double precision vel_rel_t
      Logical isincontact
      INTEGER LL
!
      double precision velm
      DOUBLE PRECISION Vloss_imp !Volume loss due to impact attrition
!
      DOUBLE PRECISION Vloss_abr !Volume loss due to abrasive attrition
!
      DOUBLE PRECISION desradius
!
      DOUBLE PRECISION desradius_new !temporarily stored desradius
!--------------------------------------------------------------------------------------!      
!     The following was used for interaction-force based impact attrition calculation. !
!--------------------------------------------------------------------------------------!
!      DOUBLE PRECISION Force_mag !magnitude of the normal force                       !
!      DOUBLE PRECISION nindex !index value for different models: Ghadiri or Evans     !
!      DOUBLE PRECISION ImpactVel                                                      !
!                                                                                      !
!Calculate the magnitude of the force                                                  !
!      Force_mag=0.0                                                                   !
!      DO I=1, DIMN                                                                    !
!         Force_mag=Force(I)**2.0+Force_mag                                            !
!      ENDDO                                                                           !
!      Force_mag=Force_mag**0.5                                                        !
!--------------------------------------------------------------------------------------! 
!Calculate impact velocity. Here DES_Radius(LL) is used because the interaction        !
!between the target particle and any of the neighboring particles occurs simultaneously!
!                                                                                      !
!      ImpactVel=0.0                                                                   !
!      ImpactVel=Force_mag/(SQRT(RO_SOL(LL)*FractureHardness)&                         !
!                   *DES_Radius(LL)**2.0)                                              !
!---------------------------------------------------------------------------------------
!Calculate the volume loss due to impact attrition                                     !
!      nindex=1.0 !Ghadiri (For Evans, it will be 2/3)                                 ! 
!      Vloss_imp=0.0                                                                   !
!      IF(ImpactVel.GE.DESAttritionThresh) THEN                                        !
!         Vloss_imp=DESAlpha/FractureToughness**(2.0*nindex)&                          !
!         *Force_mag**((3.0+nindex)/2.0)&                                              !
!         /FractureHardness**(3.0*(1.0-nindex)/2.0)                                    !
!         Vloss_imp=DESAlpha*Force_mag**2.0/FractureToughness**2.0                     ! 
!         ELSE                                                                         !
!         Vloss_imp=0.0                                                                !
!      ENDIF                                                                           !   
!--------------------------------------------------------------------------------------!
!Calculate the volume loss due to impact attrition based on relative velocity
      velm=sqrt(vel_rel_n**2.0+vel_rel_t**2.0)
    
      if (isincontact) then
         vloss_imp=0.0
      elseif (velm.ge.desattritionthresh) then
         vloss_imp=desalpha*ro_sol(LL)*velm**2.0*des_radius(ll)&
         *fracturehardness/fracturetoughness**2.0
      else
         vloss_imp=0.0
      endif

!Calculate the volume loss due to abrasive attrition
      Vloss_abr=0.0
      DO I=1, DIMN
         Vloss_abr=Vloss_abr+(Force(I)*Slide(I))**2.0
      ENDDO
      Vloss_abr=ABRAlpha*SQRT(Vloss_abr)/FractureHardness
!--------------------------------------------------------------------------------
!Calculate the new radius
      desradius_new=0.0
!
!The following equation uses desradius to add up the attrition loss by all neighboring particles
!Note here when the particle is completely worn out, the particle radius is set to a tolerance value.
!However, to be more accurate, we should remove the particle from the fluid domain.
!
      desradius_new=desradius*(1.0-(Vloss_imp+Vloss_abr))**(1.0/3.0)
      
      if(desradius_new.le.0.005) then
	desradius_new=0.005
	else
	endif

!--------------------------------------------------------------------------------
!Update the desradius. 
!
      desradius=desradius_new

      RETURN               
      END SUBROUTINE calc_attrition_des
      
