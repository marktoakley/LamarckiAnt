MODULE MODCUDALBFGS

#ifndef DUMMY_CUDA

USE COMMONS, ONLY : RMS, DEBUG, CUDATIMET, MAXBFGS, MAXERISE, CUDAPOT, NPCALL, COLDFUSION, COLDFUSIONLIMIT, MYUNIT, DGUESS, &
                    BQMAX, FREEZE, SEEDT, FREEZECORE, FROZEN, NFREEZE, QUENCHDOS, CENT, CALCQT, INTMINT, DUMPT, COMPRESSRIGIDT, & 
                    AMBER12T

USE GENRIGID, ONLY : ATOMRIGIDCOORDT, DEGFREEDOMS, NRIGIDBODY, NSITEPERBODY, RIGIDGROUPS, MAXSITE, SITESRIGIDBODY, & 
                     RIGIDSINGLES, IINVERSE, AACONVERGENCET

USE MODAMBER9, ONLY : STEEREDMINT

#endif

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_BOOL, C_CHAR

IMPLICIT NONE

#ifndef DUMMY_CUDA
INTERFACE
    SUBROUTINE CUDA_LBFGS(N, C_XCOORDS, EPS, C_MFLAG, C_ENERGY, ITMAX, C_ITDONE, C_MAXBFGS, C_MAXERISE, C_RMS, C_CUDAPOT, & 
                         C_DEBUG, C_CUDATIMET, ECALLS, C_COLDFUSION, C_COLDFUSIONLIMIT, &
                         C_DGUESS, MUPDATE, C_ATOMRIGIDCOORDT, C_DEGFREEDOMS, C_NRIGIDBODY, C_NSITEPERBODY, &
                         C_RIGIDGROUPS, C_MAXSITE, C_SITESRIGIDBODY, C_RIGIDSINGLES, C_BQMAX, C_IINVERSE, &
                         PROJECT, C_FREEZE, C_FROZEN, C_NFREEZE, POTENTIALTIME, C_AACONVERGENCET) BIND(C,NAME="setup_lbfgs")

        IMPORT :: C_INT, C_DOUBLE, C_BOOL, C_CHAR

        INTEGER(KIND=C_INT), INTENT(IN) :: N, & ! 3*no. of atoms
                                           ITMAX, & ! Max. no. of steps allowed in minmization
                                           C_DEGFREEDOMS, & ! Rigid Body Framework (RBF): no. of degrees of freedom
                                           C_NRIGIDBODY, & ! RBF: no. of rigid bodies
                                           C_MAXSITE, & ! RBF: max. no. of sites in a rigid body
                                           MUPDATE, & ! History size
                                           C_NFREEZE ! No. of frozen atoms

        INTEGER(KIND=C_INT), DIMENSION(C_NRIGIDBODY), INTENT(IN) :: C_NSITEPERBODY ! RBF: no. of rigid body sites

        INTEGER(KIND=C_INT), DIMENSION(C_MAXSITE*C_NRIGIDBODY), INTENT(IN) :: C_RIGIDGROUPS ! RBF: list of atoms in rigid bodies

        INTEGER(KIND=C_INT), DIMENSION(C_DEGFREEDOMS/3 - 2*C_NRIGIDBODY), INTENT(IN) :: C_RIGIDSINGLES ! RBF: list of atoms not in rigid bodies

        INTEGER(KIND=C_INT), INTENT(OUT) :: C_ITDONE, & ! No. of LBFGS iterations done
                                            ECALLS ! Number of potential calls made during this call to LBFGS

        REAL(KIND=C_DOUBLE), INTENT(IN) :: EPS, & ! Convergence tolerance for RMS force
                                           C_MAXBFGS, & ! Max. step size allowed
                                           C_MAXERISE, & ! Max. energy rise allowed
                                           C_COLDFUSIONLIMIT, & ! Limit below which cold fusion is diagnosed and minimization terminated
                                           C_DGUESS, & ! Initial guess for inverse Hessian diagonal elements
                                           C_BQMAX ! Sloppy quench tolerance for RMS gradient

        REAL(KIND=C_DOUBLE), DIMENSION(C_MAXSITE*3*C_NRIGIDBODY), INTENT(IN) :: C_SITESRIGIDBODY ! RBF: coordinates of the rigid body sites

        REAL(KIND=C_DOUBLE), DIMENSION(C_NRIGIDBODY*3*3), INTENT(IN) :: C_IINVERSE ! RBF: inverse eigenvalues of the unweighted tensor of gyration

        REAL(KIND=C_DOUBLE), DIMENSION(N), INTENT(INOUT) :: C_XCOORDS ! Coordinates

        REAL(KIND=C_DOUBLE), INTENT(OUT) :: C_ENERGY, & ! Energy
                                            C_RMS, & ! RMS force
                                            POTENTIALTIME ! Time taken in calculating potential - not used in GMIN

        LOGICAL(KIND=C_BOOL), INTENT(IN) :: C_DEBUG, & ! If true, print debug info.
                                            C_CUDATIMET, & ! If true, print timing info. 
                                            C_ATOMRIGIDCOORDT, & ! If false, use rigid body coordinates
                                            C_FREEZE, & ! If true, freeze some specified atoms
                                            PROJECT, & ! PROJECT is OPTIM only, always false
                                            C_AACONVERGENCET ! If true, use more accurate method of calculating rigid RMS force

        LOGICAL(KIND=C_BOOL), DIMENSION(N/3), INTENT(IN) :: C_FROZEN ! Logical array specifying frozen atoms

        LOGICAL(KIND=C_BOOL), INTENT(OUT) :: C_MFLAG, & ! True if quench converged
                                             C_COLDFUSION ! Set to true during minimization if cold fusion diagnosed

        CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: C_CUDAPOT ! Character specifying the CUDA potential to be used

    END SUBROUTINE CUDA_LBFGS
END INTERFACE

INTERFACE
    SUBROUTINE CUDA_ENEGRAD_CPUTOGPU(NATOMS, COORDS, C_TOTENERGY, C_GRADIENTS) BIND(C,NAME="gminoptim_enegrad_cputogpu")

        IMPORT :: C_INT, C_DOUBLE

        INTEGER(KIND=C_INT), INTENT(IN) :: NATOMS ! No. of atoms

        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS), INTENT(IN) :: COORDS ! Atomic coordinates
        REAL(KIND=C_DOUBLE), INTENT(OUT) :: C_TOTENERGY ! Total energy of the system
        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS), INTENT(OUT) :: C_GRADIENTS ! Gradient of the energy w.r.t. each atomic coordinate

    END SUBROUTINE CUDA_ENEGRAD_CPUTOGPU    
END INTERFACE
#endif /* DUMMY_CUDA */

CONTAINS

    SUBROUTINE CUDA_LBFGS_WRAPPER(N, MUPDATE, XCOORDS, EPS, MFLAG, ENERGY, ITMAX, ITDONE, RESET)

        ! Variables passed as *arguments through this wrapper* (not common) with intent in for CUDA_LBFGS are converted directly
#ifndef DUMMY_CUDA
        INTEGER(KIND=C_INT) :: N, MUPDATE, ITMAX, C_ITDONE, ECALLS, C_DEGFREEDOMS, C_NRIGIDBODY, C_MAXSITE, C_NFREEZE
        INTEGER(KIND=C_INT), DIMENSION(NRIGIDBODY) :: C_NSITEPERBODY
        INTEGER(KIND=C_INT), DIMENSION(MAXSITE*NRIGIDBODY) :: C_RIGIDGROUPS
        INTEGER(KIND=C_INT), DIMENSION(DEGFREEDOMS/3 - 2*NRIGIDBODY) :: C_RIGIDSINGLES

        REAL(KIND=C_DOUBLE) :: EPS, C_MAXBFGS, C_MAXERISE, C_ENERGY, C_RMS, C_COLDFUSIONLIMIT, C_DGUESS, C_BQMAX, POTENTIALTIME
        REAL(KIND=C_DOUBLE), DIMENSION(MAXSITE*3*NRIGIDBODY) :: C_SITESRIGIDBODY
        REAL(KIND=C_DOUBLE), DIMENSION(NRIGIDBODY*3*3) :: C_IINVERSE
        REAL(KIND=C_DOUBLE), DIMENSION(N) :: C_XCOORDS

        LOGICAL(KIND=C_BOOL) :: C_MFLAG, C_DEBUG, C_CUDATIMET, C_ATOMRIGIDCOORDT, C_COLDFUSION, C_FREEZE, PROJECT, C_AACONVERGENCET
        LOGICAL(KIND=C_BOOL), DIMENSION(N/3) :: C_FROZEN

        CHARACTER(LEN=1, KIND=C_CHAR) :: C_CUDAPOT
#endif
        ! Same as above, but now dimension(1).
#ifdef DUMMY_CUDA
        INTEGER(KIND=C_INT) :: N, MUPDATE, ITMAX, C_ITDONE, ECALLS, C_DEGFREEDOMS, C_NRIGIDBODY, C_MAXSITE, C_NFREEZE
        INTEGER(KIND=C_INT), DIMENSION(1) :: C_NSITEPERBODY
        INTEGER(KIND=C_INT), DIMENSION(1) :: C_RIGIDGROUPS
        INTEGER(KIND=C_INT), DIMENSION(1) :: C_RIGIDSINGLES

        REAL(KIND=C_DOUBLE) :: EPS, C_MAXBFGS, C_MAXERISE, C_ENERGY, C_RMS, C_COLDFUSIONLIMIT, C_DGUESS, C_BQMAX, POTENTIALTIME
        REAL(KIND=C_DOUBLE), DIMENSION(1) :: C_SITESRIGIDBODY
        REAL(KIND=C_DOUBLE), DIMENSION(1) :: C_IINVERSE
        REAL(KIND=C_DOUBLE), DIMENSION(1) :: C_XCOORDS

        LOGICAL(KIND=C_BOOL) :: C_MFLAG, C_DEBUG, C_CUDATIMET, C_ATOMRIGIDCOORDT, C_COLDFUSION, C_FREEZE, PROJECT, C_AACONVERGENCET
        LOGICAL(KIND=C_BOOL), DIMENSION(N/3) :: C_FROZEN

        CHARACTER(LEN=1, KIND=C_CHAR) :: C_CUDAPOT
#endif

        ! Variables passed as *arguments through this wrapper* (not common) with intent out for CUDA_LBFGS are not passed into it
        ! Therefore uninitialised C types are passed in and converted types are copied back after the call

        INTEGER :: I, J, K, ITDONE
        DOUBLE PRECISION :: ENERGY, POTEL
        DOUBLE PRECISION, DIMENSION(N) :: XCOORDS
        LOGICAL :: MFLAG, RESET
#ifndef DUMMY_CUDA
        COMMON /MYPOT/ POTEL

        IF (.NOT. RESET) THEN
            WRITE(MYUNIT,'(A)') "modcudalbfgs> Warning: LBFGS resetting, though RESET is false. "
        END IF

        IF (DUMPT .AND. DEBUG) THEN
            WRITE(MYUNIT,'(A)') "modcudalbfgs> Warning: printing behaviour of DUMP and DEBUG during minimization is not implemented. "
        END IF

        IF ((SEEDT.AND.FREEZECORE) .OR. STEEREDMINT .OR. QUENCHDOS .OR. CENT .OR. CALCQT .OR. INTMINT .OR. COMPRESSRIGIDT) THEN
            WRITE(MYUNIT,'(A)') "modcudalbfgs> Keyword SEED with FREEZECORE is not yet supported. SEED can be used with NOFREEZE. &
                                 Keywords STEEREDMIN, QUENCHDOS, CENTRE, CALCQ, INTMIN, COMPRESSRIGID are not yet supported. & 
                                 Contact rgm38 if you would like a feature to be added. "
            WRITE(MYUNIT,'(A)') "modcudalbfgs> Disclaimer: this list might not be exhaustive! "
            STOP
        END IF

        ! Variables from common blocks or modules with intent in or inout are copied into C types

        DO K = 1,NRIGIDBODY
            DO J = 1,3
                DO I = 1,MAXSITE
                     C_SITESRIGIDBODY((K - 1)*3*MAXSITE + (J - 1)*MAXSITE + I) = SITESRIGIDBODY(I,J,K)
                END DO
            END DO
        END DO

        DO J = 1,NRIGIDBODY
            DO I = 1,MAXSITE
                C_RIGIDGROUPS((J - 1)*MAXSITE + I) = RIGIDGROUPS(I,J)
            END DO
        END DO

        DO I = 1,NRIGIDBODY
            C_NSITEPERBODY(I) = NSITEPERBODY(I)
        END DO

        DO I = 1,(DEGFREEDOMS/3 - 2*NRIGIDBODY)
            C_RIGIDSINGLES(I) = RIGIDSINGLES(I)
        END DO

        DO K = 1,3
            DO J = 1,3
                DO I = 1,NRIGIDBODY
                   C_IINVERSE((K - 1)*3*NRIGIDBODY + (J - 1)*NRIGIDBODY + I) = IINVERSE(I,J,K)
                END DO
            END DO
        END DO

        DO I = 1,N
            C_XCOORDS(I) = XCOORDS(I)
        END DO

        DO I = 1,(N/3)
            C_FROZEN(I) = FROZEN(I)
        END DO

        C_CUDAPOT = CUDAPOT
        C_DEBUG = DEBUG
        C_CUDATIMET = CUDATIMET
        C_MAXBFGS = MAXBFGS
        C_MAXERISE = MAXERISE
        C_ATOMRIGIDCOORDT = ATOMRIGIDCOORDT
        C_DEGFREEDOMS = DEGFREEDOMS
        C_NRIGIDBODY = NRIGIDBODY
        C_MAXSITE = MAXSITE
        C_COLDFUSIONLIMIT = COLDFUSIONLIMIT
        C_DGUESS = DGUESS
        C_BQMAX = BQMAX
        C_FREEZE = FREEZE
        C_NFREEZE = NFREEZE
        PROJECT = .FALSE.
        C_AACONVERGENCET = AACONVERGENCET

        ! 'C_' prefix denotes those variables which have intent out or inout or are copies of those from common blocks/modules 
        CALL CUDA_LBFGS(N, C_XCOORDS, EPS, C_MFLAG, C_ENERGY, ITMAX, C_ITDONE, C_MAXBFGS, C_MAXERISE, C_RMS, C_CUDAPOT, & 
                       C_DEBUG, C_CUDATIMET, ECALLS, C_COLDFUSION, C_COLDFUSIONLIMIT, C_DGUESS, MUPDATE, & 
                       C_ATOMRIGIDCOORDT, C_DEGFREEDOMS, C_NRIGIDBODY, C_NSITEPERBODY, C_RIGIDGROUPS, C_MAXSITE, & 
                       C_SITESRIGIDBODY, C_RIGIDSINGLES, C_BQMAX, C_IINVERSE, PROJECT, C_FREEZE, C_FROZEN, C_NFREEZE, & 
                       POTENTIALTIME, C_AACONVERGENCET)

        ! Make sure C types with intent out or inout are coverted back to Fortran ones

        DO I = 1,N
            XCOORDS(I) = DBLE(C_XCOORDS(I))
        END DO

        ENERGY = DBLE(C_ENERGY)
        RMS = DBLE(C_RMS)
        MFLAG = LOGICAL(C_MFLAG)
        ITDONE = INT(C_ITDONE)
        COLDFUSION = LOGICAL(C_COLDFUSION)

        IF (COLDFUSION) THEN
            WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
            WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENERGY,COLDFUSIONLIMIT
            ENERGY = 1.0D6
            POTEL = 1.0D6
            RMS = 1.0D0
        END IF

        NPCALL = NPCALL + INT(ECALLS)
#endif /* DUMMY_CUDA */

    END SUBROUTINE CUDA_LBFGS_WRAPPER

    SUBROUTINE CUDA_ENEGRAD_WRAPPER(NATOMS, COORDS, TOTENERGY, GRADIENTS)

        INTEGER(KIND=C_INT) :: NATOMS

        REAL(KIND=C_DOUBLE) :: C_TOTENERGY
        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS) :: COORDS, C_GRADIENTS

        INTEGER :: X

        DOUBLE PRECISION :: TOTENERGY
        DOUBLE PRECISION, DIMENSION(3*NATOMS) :: GRADIENTS

#ifndef DUMMY_CUDA
        ! Calculates the energy and gradients on the GPU using the GB potential
        CALL CUDA_ENEGRAD_CPUTOGPU(NATOMS, COORDS, C_TOTENERGY, C_GRADIENTS)

        TOTENERGY = DBLE(C_TOTENERGY)
        
        DO X = 1,(3*NATOMS)
            GRADIENTS(X) = DBLE(C_GRADIENTS(X))
        END DO
#endif
    END SUBROUTINE CUDA_ENEGRAD_WRAPPER

    SUBROUTINE CUDA_NUMERICAL_HESS(NATOMS, COORDS, HESSIAN, DELTA)
        IMPLICIT NONE

        INTEGER(KIND=C_INT) :: NATOMS
        REAL(KIND=C_DOUBLE) :: C_ENERGY
        REAL(KIND=C_DOUBLE) :: COORDS(3*NATOMS), C_GRADIENTS(3*NATOMS)
        DOUBLE PRECISION    :: HESSIAN(3*NATOMS, 3*NATOMS)
        DOUBLE PRECISION    :: DELTA
        DOUBLE PRECISION    :: GRAD_PLUS(3*NATOMS), GRAD_MINUS(3*NATOMS)
        INTEGER             :: I

#ifndef DUMMY_CUDA
        DO I = 1, 3*NATOMS
            ! Plus
            COORDS(I) = COORDS(I) + DELTA
            CALL CUDA_ENEGRAD_CPUTOGPU(NATOMS, COORDS, C_ENERGY, C_GRADIENTS)
            GRAD_PLUS(:) = DBLE(C_GRADIENTS(:))
            ! Minus
            COORDS(I) = COORDS(I) - 2.0D0 * DELTA
            CALL CUDA_ENEGRAD_CPUTOGPU(NATOMS, COORDS, C_ENERGY, C_GRADIENTS)
            GRAD_MINUS(:) = DBLE(C_GRADIENTS(:))
            ! Reset coords
            COORDS(I) = COORDS(I) + DELTA
            ! Calculate hessian
            HESSIAN(I, :) = (GRAD_PLUS(:) - GRAD_MINUS(:)) / (2.0D0 * DELTA)
        END DO
#endif
    END SUBROUTINE CUDA_NUMERICAL_HESS

END MODULE MODCUDALBFGS
