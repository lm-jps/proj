MODULE CHANGE_VAR
!
! July 20, 2011
! R. Centeno
! Module that contains routines to compute specific variable changes in the 
! model atmosphere parameters. These routines are called only from the main
! program (invert.f90), and only affect the derivatives of the chi2 function
! (DO_CHANGE_DER takes care of this).
! The gradient and Hessian are computed with the new variables, hence the 
! perturbations to the model are obtained for the new combination of parameters.
! The variable change has to be undone in order to correct the model atmosphere
! for the next iteration (UNDO_CHANGE_DMODEL).
!
! If we want to add any other variable change, this will be the only module that
! needs modifying.


CONTAINS
!
!-----------------------------------------------------
!
! Changes the variables of MODEL and stores the new variables in MODELC
!
  SUBROUTINE DO_CHANGE_VAR(MODEL, MODELC)

    USE CONS_PARAM
    IMPLICIT NONE
    REAL(DP),              DIMENSION(10)          :: MODEL, MODELC

    MODELC = MODEL
    ! We compute sqrt(eta0)*DopplerWidth instead DopplerWidth
    MODELC(5) = MODEL(5) * dsqrt(MODEL(1)) 
    ! We compute the sum of the Source Function and its gradient instead of its gradient.
    MODELC(9) = MODEL(9) + MODEL(8) 

  END SUBROUTINE DO_CHANGE_VAR
!
!-----------------------------------------------------
!
! Routine that undoes the variable change of the perturbations to the model
!
  SUBROUTINE UNDO_CHANGE_DMODEL(INMODEL, DMODEL)

    USE CONS_PARAM
    USE FILT_PARAM
    IMPLICIT NONE
    
    REAL(DP),   DIMENSION(10)          :: INMODEL, DMODEL 
    
    ! INMODEL is the original model (no variable change)
    ! DMODEL is the perturbation obtained from the inversion, with variable change. 
    ! We overwrite it to obtained the perturbation to the original variables

    DMODEL(5) = (DMODEL(5) - 0.5/DSQRT(INMODEL(1))*INMODEL(5)*DMODEL(1)) / &
         (DSQRT(INMODEL(1)))
!         (0.5/DSQRT(INMODEL(1))*DMODEL(1) + DSQRT(INMODEL(1)))
    DMODEL(9) = DMODEL(9) - DMODEL(8)

  END SUBROUTINE UNDO_CHANGE_DMODEL

!
! -------------------------------------------------------
!
! Changing the variables in the derivatives
!
 SUBROUTINE DO_CHANGE_DER(INMODEL, DSYN)

    USE CONS_PARAM
    USE FILT_PARAM
    IMPLICIT NONE
    
    REAL(DP),   DIMENSION(10)          :: INMODEL, CHMODEL
    REAL(DP), INTENT(OUT),  DIMENSION(10,NBINS,4) :: DSYN

    ! INMODEL is the original model with no variable change
    ! CHMODEL is the model with the variable change.
    ! DSYN are the derivatives without the variable change. They are overwritten with 
    ! those that do have the change in variable.


    CALL DO_CHANGE_VAR(INMODEL, CHMODEL)

    ! For eta0*DopplerWidth 
    DSYN(1,:,:)=DSYN(1,:,:)-0.5*CHMODEL(5)/CHMODEL(1)/dsqrt(CHMODEL(1))*DSYN(5,:,:)
    DSYN(5,:,:)=DSYN(5,:,:)/dsqrt(CHMODEL(1)) 
    ! For Source Function + Gradient
    DSYN(8,:,:)=DSYN(8,:,:)-DSYN(9,:,:)

  END SUBROUTINE DO_CHANGE_DER


END MODULE CHANGE_VAR
!CVSVERSIONINFO "$Id: change_var.f90,v 1.1 2011/10/14 17:24:47 keiji Exp $"
