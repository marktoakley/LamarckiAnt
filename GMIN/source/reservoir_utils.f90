SUBROUTINE RESERVOIR_EXCHANGE(X,Y,Z,XO,YO,ZO,VOLD, &                                           
                              VNEW,GRAD,BETA, &
                              LBFGS_ITERATIONS,MCSTEP,ITRAJ)

   USE COMMONS 

   IMPLICIT NONE

   DOUBLE PRECISION :: X(NATOMS),Y(NATOMS),Z(NATOMS),Xhsa(NATOMS),Yhsa(NATOMS),Zhsa(NATOMS),XO(NATOMS),YO(NATOMS),ZO(NATOMS),E
   DOUBLE PRECISION :: POTEL,POTEL_HSA,BETA(0:NPAR-1),GRAD(3*NATOMS),VOLD,VNEW
   DOUBLE PRECISION :: CONFIG(3*NATOMS),WORKF,WORKR,WORK,W,DPRAND,RANDOM,MCSTEP
   INTEGER LBFGS_ITERATIONS,WELLLOW,WELLphys,J1,ITRAJ

   DOUBLE PRECISION Xpre(NATOMS),Ypre(NATOMS),Zpre(NATOMS),Xhpre(NATOMS),Yhpre(NATOMS),Zhpre(NATOMS)

   INTEGER tau

   ! calc WorkF
   DO J1=1,NATOMS
      CONFIG(3*(J1-1)+1) = X(J1)
      CONFIG(3*(J1-1)+2) = Y(J1)
      CONFIG(3*(J1-1)+3) = Z(J1)
   ENDDO

   CALL POTENTIAL(CONFIG,GRAD,POTEL,.FALSE.,.FALSE.)

   Xpre(:) = X(:)
   Ypre(:) = Y(:)
   Zpre(:) = Z(:)

   
   CALL RESERVOIR_CALC_EHSA(X,Y,Z,WELLphys,POTEL,POTEL_HSA,BETA,LBFGS_ITERATIONS)
   WORKR = BETA_RES*POTEL_HSA - BETA(MYNODE)*POTEL
   
   !tau = 1 
   !CALL switch(X,Y,Z,WORKR,WELLphys,tau,1,BETA,LBFGS_ITERATIONS)

   ! sample config from HSA
   CALL RESERVOIR_SAMPLE_HSA(Xhsa,Yhsa,Zhsa,WELLLOW,POTEL,POTEL_HSA,BETA_RES)
   Xhpre(:) = Xhsa(:)
   Yhpre(:) = Yhsa(:)
   Zhpre(:) = Zhsa(:)
   !CALL switch(Xhsa,Yhsa,Zhsa,WORKF,WELLLOW,tau,-1,BETA,LBFGS_ITERATIONS)
   WORKF = BETA(USERES)*POTEL - BETA_RES*POTEL_HSA

   WORK = WORKF + WORKR
   
   W=MIN(1.0D0,DEXP(-WORK))
   ! check if NaN
   IF(WORK.EQ.(WORK+1.0)) W=0.0

   WRITE(MYUNIT,'(A,4(G14.7,2x),2I4,F20.1,I10.1)') 'reservoir> hsawork:  ',WORKF,WORKR,WORK,W,WELLLOW,WELLphys,MCSTEP,ITRAJ
   !WRITE(MYUNIT,'(A,4(G14.7,xx),2I4,F20.1)') 'reservoir> hsawork: ',WORKF,WORKR,WORK,W,WELLLOW,WELLphys,MCSTEP
   !WRITE(MYUNIT,*) 'reservoir> hsawork: ',WORKF,WORKR,WORK,W,WELLLOW,WELLphys,MCSTEP

   RANDOM=DPRAND()

   IF (W.GT.RANDOM) THEN
      ! move is accepted

      XO(:) = X(:)
      YO(:) = Y(:)
      ZO(:) = Z(:)
      X(:) = Xhsa(:)
      Y(:) = Yhsa(:)
      Z(:) = Zhsa(:)

      DO J1=1,NATOMS
         COORDS(3*(J1-1)+1,MYNODE+1)=X(J1)
         COORDS(3*(J1-1)+2,MYNODE+1)=Y(J1)
         COORDS(3*(J1-1)+3,MYNODE+1)=Z(J1)
      ENDDO
      
      VOLD = POTEL
      E = POTEL 

      VNEW = POTEL
      !VOLD = E
      !E = POTEL
   ENDIF

   RETURN

END SUBROUTINE RESERVOIR_EXCHANGE

SUBROUTINE switch(Xs,Ys,Zs,WORK,WELLLOW,tau,par,BETA,LBFGS_ITERATIONS)

   USE COMMONS
  
   ! par =-1 for forward : hsa  -> phys direction
   ! par = 1 for reverse : phys -> hsa  direction

   IMPLICIT NONE

   INTEGER time,tau,exab_count,WELLLOW,J1,LBFGS_ITERATIONS,par

   DOUBLE PRECISION Xso(NATOMS),Yso(NATOMS),Zso(NATOMS),Xs(NATOMS),Ys(NATOMS),Zs(NATOMS),PERTCOORD(3*NATOMS)
   DOUBLE PRECISION w,WORK,Ephys,Ephysnew,GRAD(3*NATOMS),vnew,Ehsa,Ehsapert,Ephyspert,Eps,Epspert
   logical :: jumpt
   double precision peint, peqv(NENRPER),BETA(0:NPAR-1)

   ! bug dirty fix
   JUMPT = .FALSE.

   Xso(:) = Xs(:)
   Yso(:) = Ys(:)
   Zso(:) = Zs(:)

   Ephys=0.0
   Ehsa=0.0
   CALL RESERVOIR_CALC_EHSA(Xs,Ys,Zs,WELLLOW,Ephys,Ehsa,BETA,LBFGS_ITERATIONS)

   w = 0
   DO time=1,tau-1
      w = w + par * (1./(tau)) * (Ehsa-Ephys) 

      ! RELAXATION at constant lambda !!!!
         ! MC trial move
      CALL BSPT_TAKESTEP(Xs,Ys,Zs,peint,peqv,vnew,exab_count)
      IF (DEBUG) THEN
         IF (par.eq.1) THEN
            WRITE (MYUNIT,"(A,2G20.10,I3.1)") "switch> HSA: Ehsa, Ephys ",Ehsa,Ephys,time
         ENDIF
      ENDIF
         ! alignment 
      CALL RESERVOIR_CALC_EHSA(Xs,Ys,Zs,WELLLOW,Ephyspert,Ehsapert,BETA,LBFGS_ITERATIONS)
         ! acceptance step
      CALL myaccrej(Xs,Ys,Zs,Xso,Yso,Zso,Ephys,Ephyspert,Ehsa,Ehsapert,par,time,tau,BETA,LBFGS_ITERATIONS)
      !!!!
   ENDDO
   CALL RESERVOIR_CALC_EHSA(Xs,Ys,Zs,WELLLOW,Ephys,Ehsa,BETA,LBFGS_ITERATIONS)
   w = w + par*(1./(tau)) * (Ehsa-Ephys) 

   !IF (par.eq.-1) THEN
   IF (DEBUG) THEN
      WRITE (MYUNIT,"(A,I2.0,2G20.10,2I3.0,G20.10)") "switch> HSA: Ehsa, Ephys ",USERES,Ehsa,Ephys,time,par&
        & ,Xs(1)
   ENDIF
   WORK = BETA(USERES) * w
   !IF (Ehsa.gt.100.0.and.par.eq.-1) THEN
   !   DO J1=1,NATOMS
   !      WRITE (MYUNIT,"(3G20.10)") Xs(J1),Ys(J1),Zs(J1)
   !   ENDDO
   !   !STOP
   !ENDIF


   RETURN

END SUBROUTINE switch

SUBROUTINE myaccrej(Xa,Ya,Za,Xao,Yao,Zao,Ephys,Ephyspert,Ehsa,Ehsapert,par,time,tau,BETA,LBFGS_ITERATIONS)

   ! currently only coded for BETA_RES = BETA(USERES) !

   USE COMMONS

   IMPLICIT NONE

   INTEGER LBFGS_ITERATIONS,par,time,tau
   DOUBLE PRECISION Xa(NATOMS),Ya(NATOMS),Za(NATOMS),Xao(NATOMS),Yao(NATOMS),Zao(NATOMS)
   DOUBLE PRECISION RANDOM,Weight,BETA(0:NPAR-1),DPRAND,Ephys,Ephyspert,Ehsa,Ehsapert,Eps,Epspert

   RANDOM = DPRAND()

   ! calculation of E_lambda = H_lambda(x) := H_i(x) + lambda * (H_f(x) -H_i(x))
   Eps = Ephys + (0.5*(1.-par) + par*time*1./(tau)) * (Ehsa- Ephys)
   Epspert = Ephyspert + (0.5*(1.-par) + par*time*1./(tau)) * (Ehsapert - Ephyspert)
   !Epspert = Ephyspert + (time*1./(tau)) * (Ehsapert - Ephyspert)

   Weight=MIN(1.0D0,EXP(-BETA(USERES)*(Epspert-Eps)))

   !IF (par.eq.-1) THEN
   !   WRITE (MYUNIT,"(A,2G20.10)") "acc> HSA: Xs, Xso ",Xa(1),Xao(1)
   !ENDIF
   !WRITE(MYUNIT,'(A,5G20.10,I6,4G20.10)') "switch> xo x",Xao(1),Xa(1),Weight,Epspert,Eps,par,Ephys,Ephyspert,Ehsa,Ehsapert

   IF (Weight.LT.RANDOM) THEN
      ! rej
      Xa(:) = Xao(:)
      Ya(:) = Yao(:)
      Za(:) = Zao(:)
   ELSE
      ! accept
      Xao(:) = Xa(:)
      Yao(:) = Ya(:)
      Zao(:) = Za(:)
      Ehsa = Ehsapert
      Ephys = Ephyspert
   ENDIF

   RETURN
END SUBROUTINE myaccrej


SUBROUTINE RESERVOIR_CALC_EHSA(Xl,Yl,Zl,WELLLOW,POTEL,POTEL_HSA,BETA,LBFGS_ITERATIONS)

   USE COMMONS 

   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(IN) :: Xl(NATOMS),Yl(NATOMS),Zl(NATOMS)
   !DOUBLE PRECISION, INTENT(INOUT) :: POTEL
   DOUBLE PRECISION POTEL,DOSSTATS(3*NATOMS)
   DOUBLE PRECISION POTELOLD,DISTANCE,VNEWLOW,VNEW
   DOUBLE PRECISION, INTENT(OUT) :: POTEL_HSA
   DOUBLE PRECISION, INTENT(IN) :: BETA(0:NPAR-1)
   INTEGER, INTENT(OUT) :: WELLLOW

   LOGICAL CONVERGEDab
   DOUBLE PRECISION CONFIGQUENCHED(3*NATOMS),CONFIGALIGNED(3*NATOMS),POTELquenched,DVWELL,DUMMY,NDUMMY,CONFIGQ(3*NATOMS)
   DOUBLE PRECISION CONFIG(3*NATOMS),GRAD(3*NATOMS),DISTANCE2
   DOUBLE PRECISION RMATBEST(3,3),DIST2,DISTTOL,DISTANCEq
   INTEGER J1,J2,J3,LBFGS_ITERATIONS

   !ab2111>
   DISTTOL=0.01 ! for aligned structure tolerance for minpermdist

   !ab2111 todo: test this
   DO J1=1,NATOMS
      CONFIGQUENCHED(3*(J1-1)+1) = Xl(J1)
      CONFIGQUENCHED(3*(J1-1)+2) = Yl(J1)
      CONFIGQUENCHED(3*(J1-1)+3) = Zl(J1)
   ENDDO

   CONFIGALIGNED(:) = CONFIGQUENCHED(:)

   CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGEDab,CONFIGQUENCHED(:),DOSSTATS)
   CALL POTENTIAL(CONFIGQUENCHED(:),GRAD,POTELquenched,.FALSE.,.FALSE.)

   WELLLOW=0
   DVWELL=10000000.0
   DO J1=1,NRESMIN
      !WRITE(*,'(A,G20.10,I6,G20.10)') "emin dev ", POTELquenched,J1, ABS(EMIN(J1)-POTELquenched)
      IF ((DVWELL.GT.ABS(EMIN(J1)-POTELquenched)).AND.(ABS((EMIN(J1)-POTELquenched)/POTELquenched).LT.0.01)) THEN
         DVWELL=ABS(EMIN(J1)-POTELquenched)
         WELLLOW=J1
      ENDIF
   ENDDO

   !WRITE(MYUNIT,'(A,G20.10,I2.0)') "reservoir_utils>: potential quenched ", POTELquenched, WELLLOW
   !ab2111> debug: only well 1
   !WELLLOW=1

   !DO J2=1,3*NATOMS,3
   !   WRITE(*,'(3G20.10)') CONFIGALIGNED(J2),CONFIGALIGNED(J2+1),CONFIGALIGNED(J2+2)
   !ENDDO
   CALL MINPERMDIST(RESPOINTS(:,WELLLOW),CONFIGQUENCHED(:),NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,&
&       TWOD,DISTANCEq,DIST2,RIGID,RMATBEST)
   !WRITE(*,'(A,3G20.10)') "reservoir_utils> rmatbest",CONFIGALIGNED(1),CONFIGALIGNED(2),CONFIGALIGNED(3)
!   CALL MINPERMDIST(RESPOINTS(:,WELLLOW),CONFIGALIGNED(:),NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,&
!&       TWOD,DISTANCEq,DIST2,RIGID,RMATBEST)

   !DO J1=1,3
   !   DO J2=1,3
   !      WRITE(*,'(A,2G20.10)') "r", RMATBEST(J1,J2),DISTANCEq
   !   ENDDO
   !ENDDO
   !DO J2=1,3*NATOMS,3
   !   WRITE(*,'(3G20.10)') CONFIGALIGNED(J2),CONFIGALIGNED(J2+1),CONFIGALIGNED(J2+2)
   !ENDDO
   !STOP
   

   IF (DISTANCE2.GT.DISTTOL) THEN
      WRITE(*,'(A,G10.10)') "reservoir_utils> Distance match for quenched configs not satisfied",DISTANCE
      !STOP
   ENDIF 

   !WRITE(MYUNIT,'(A,3G20.10)') "config before alignment> ",CONFIGALIGNED(1),CONFIGALIGNED(3),CONFIGALIGNED(4)
   !WRITE(*,'(A,6G20.10)') "config before alignment> ",CONFIGALIGNED(1),CONFIGALIGNED(3),CONFIGALIGNED(4)&
   !  &,DISTANCE,RMATBEST(1,1),RESPOINTS(1,1)

   !ab2111> debug
   !CALL MINPERMDIST(RESPOINTS(:,WELLLOW),CONFIGALIGNED(:),NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,&
   !  &       TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

   CALL MINPERMDIST(CONFIGQUENCHED(:),CONFIGALIGNED(:),NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,.FALSE.,&
     &       TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

   !WRITE(*,'(A,6G20.10)') "config after alignment> ",CONFIGALIGNED(1),CONFIGALIGNED(3),CONFIGALIGNED(4)&
   !  &,DISTANCE,RMATBEST(1,2),RESPOINTS(1,1)
   !STOP
     !ab2111 
   !POTEL_HSA = EMIN(WELLLOW) - (1./BETA(USERES)) * LOG(PW(WELLLOW))
   POTEL_HSA = EMIN(WELLLOW) 
   VNEW = EMIN(WELLLOW) 

   DO J1=1,3*NATOMS
      DO J2=1,3*NATOMS
         DO J3=1,3*NATOMS
            POTEL_HSA = POTEL_HSA + 0.5 * (CONFIGALIGNED(J1) - RESPOINTS(J1,WELLLOW)) * HESSEIGVC(J1,J2,WELLLOW) &
   &          * HESSEIGVA(J2,WELLLOW) * HESSEIGVC(J3,J2,WELLLOW) * (CONFIGALIGNED(J3) - RESPOINTS(J3,WELLLOW))
            VNEW = VNEW + 0.5 * (CONFIGQUENCHED(J1) - RESPOINTS(J1,WELLLOW)) * HESSEIGVC(J1,J2,WELLLOW) &
     &       * HESSEIGVA(J2,WELLLOW) * HESSEIGVC(J3,J2,WELLLOW) * (CONFIGQUENCHED(J3) - RESPOINTS(J3,WELLLOW))
         ENDDO
      ENDDO
   ENDDO
   !ab2111> debug
   IF (.FALSE.) THEN
   !IF ((BETA(MYNODE)*POTEL_HSA).GE.10.0) THEN
       ! find new minimum 
       CONVERGEDab=.FALSE.
       CONFIGQ(:)=CONFIG(:)
       POTELOLD=POTEL
       CALL QUENCH(.FALSE.,MYNODE+1,LBFGS_ITERATIONS,DUMMY,NDUMMY,CONVERGEDab,CONFIG(:),DOSSTATS)
       !CALL QUENCH(.FALSE.,MYNODE+1,100,DUMMY,NDUMMY,CONVERGEDab,CONFIGQ(:),DOSSTATS)
       IF (CONVERGEDab.NEQV..TRUE.) WRITE(MYUNIT, '(A)') 'bspt> WARNING - quench did not converge' 
       CALL POTENTIAL(CONFIG(:),GRAD,POTEL,.TRUE.,.FALSE.)
       !WRITE(MYUNIT,'(A,G14.7,A,G14.7)') 'new quenched energy: ',POTELOLD,' --->',POTEL
       !WRITE(MYUNIT,'(A,I4,A,I4)') 'WELLLOW VWELLLOW ',WELLLOW,VNEWLOW
   ENDIF

   ! Calculate work associated with switching to other replica
   CALL POTENTIAL(CONFIGALIGNED(:),GRAD,POTEL,.FALSE.,.FALSE.)
   !ab2111> debug
   !WRITE(MYUNIT,'(A,G14.7,A,2G14.7,I4,A,3G14.7)') 'Final Work elements JLOW+1',POTEL,'--->',POTEL_HSA,VNEW,WELLLOW," Ephys(quen)",POTELquenched&
   !  & ,DISTANCEq,DISTANCE

   IF (POTEL_HSA.GT.0.0) THEN
   !IF (POTEL_HSA> 100.0) THEN
      !DO J1=1,NATOMS
      !   WRITE(MYUNIT,'(3F20.10)') CONFIGALIGNED(3*J1-2),CONFIGALIGNED(3*J1-1),CONFIGALIGNED(3*J1)
      !   !WRITE(MYUNIT,'(3F20.10)') CONFIGQUENCHED(3*J1-2),CONFIGQUENCHED(3*J1-1),CONFIGQUENCHED(3*J1)
      !ENDDO
      !STOP
   ENDIF

   RETURN 
   END SUBROUTINE RESERVOIR_CALC_EHSA

subroutine RESERVOIR_SAMPLE_HSA(Xhsa,Yhsa,Zhsa,WELLLOW,POTEL,POTEL_HSA,MYBETA)

   USE COMMONS 
   use random_normal_module

   IMPLICIT NONE

   DOUBLE PRECISION, INTENT(INOUT) :: Xhsa(NATOMS),Yhsa(NATOMS),Zhsa(NATOMS),POTEL
   DOUBLE PRECISION, INTENT(OUT) :: POTEL_HSA
   DOUBLE PRECISION, INTENT(IN) :: MYBETA 
   INTEGER, INTENT(OUT) :: WELLLOW

   LOGICAL CONVERGEDab
   DOUBLE PRECISION CONFIGQUENCHED(3*NATOMS),CONFIGALIGNED(3*NATOMS),POTELquenched,DVWELL,DUMMY,NDUMMY,GRAD(3*NATOMS)
   DOUBLE PRECISION RMATBEST(3,3),DIST2,DISTTOL,RANDOM,RANDVEC(3*NATOMS),CONFIG(3*NATOMS),DPRAND,PNEW,VNEW
   INTEGER J1,J2,J3,WELL

   ! choose well number
   J1=0
   DUMMY=0.0
   RANDOM = DPRAND()
   DO WHILE (DUMMY < RANDOM)
      J1=J1+1
      DUMMY = DUMMY + PW(J1)
   ENDDO
   WELL=J1

   !ab2111> debug: sample only from well 1
   !WELL=1
   !WRITE(*,'(A,G20.10,I6)') "reservoir> random WELL", RANDOM,WELL

   ! generate random vector
   DO J1=1,1000
      DO J2=1,3*NATOMS-6
         RANDOM = RANDOM_NORMAL()
         RANDVEC(J2) = RANDOM
      ENDDO
   ENDDO

   ! set last 6 elements to 0 to ensure no sampling occurs in translation /
   ! rotation directions
   DO J2=3*NATOMS-5,3*NATOMS
      RANDVEC(J2) = 0.0
   ENDDO

   ! draw sample from reservoir
   DO J2=1,3*NATOMS
      CONFIG(J2) = RESPOINTS(J2,WELL)
      DO J1=1,3*NATOMS-6 ! setting summation bound of NATOMS-2 sets 1/HESSEIGVA = 0
         IF (DEBUG) THEN
            IF (HESSEIGVC(J2,J1,WELL).EQ.0.0) THEN
               WRITE(MYUNIT,'(A,2I2)') 'error: zeros!',J2,J1
            ENDIF
         ENDIF
         CONFIG(J2) = CONFIG(J2) + HESSEIGVC(J2,J1,WELL) * DSQRT(1./MYBETA)&
     &    * (1./DSQRT(HESSEIGVA(J1,WELL))) * RANDVEC(J1)
      ENDDO
   ENDDO

   ! simple calculation of Energy
   WELLLOW=WELL
   !ab2111
   !POTEL_HSA = EMIN(WELLLOW) - (1./MYBETA) * LOG(PW(WELLLOW))
   POTEL_HSA = EMIN(WELLLOW) 
   DO J1=1,3*NATOMS
      POTEL_HSA = POTEL_HSA + 0.5 * (1./MYBETA)  * RANDVEC(J1) * RANDVEC(J1)
   ENDDO

   ! calculate energy of selected minimum
   PNEW=0.0
   !ab2111 debug>
   !DO WELL=1,NRESMIN
   DO WELL=1,WELLLOW
      VNEW = EMIN(WELL) 
      DO J1=1,3*NATOMS
         ! consistency check for energy / sampling
         DO J2=1,3*NATOMS
            DO J3=1,3*NATOMS
   !            VNEW = VNEW + 0.5 * (CONFIG(J1) - RESPOINTS(J1,1)) * HESSEIGVC(J1,J2) &
   !     &       * HESSEIGVA(J2) * HESSEIGVC(J3,J2) * (CONFIG(J3) - RESPOINTS(J3,1))
               VNEW = VNEW + 0.5 * (CONFIG(J1) - RESPOINTS(J1,WELL)) * HESSEIGVC(J1,J2,WELL) &
        &       * HESSEIGVA(J2,WELL) * HESSEIGVC(J3,J2,WELL) * (CONFIG(J3) - RESPOINTS(J3,WELL))
            ENDDO
         ENDDO

      ENDDO
      !PNEW = PNEW + EXP(-MYBETA * VNEW)
   ENDDO
   !VNEW = -(1./MYBETA) * LOG(PNEW)


   CALL POTENTIAL(CONFIG(:),GRAD,POTEL,.FALSE.,.FALSE.)

!   WRITE(*,'(A,I6,3F20.10,A,F20.10)') 'reservoir> Energy of minimum, sampled config:',&
!     &     WELLLOW,EMIN(WELLLOW),BETA(USERES)*POTEL_HSA,BETA(USERES)*POTEL," wF(x)",BETA(USERES)*(POTEL-POTEL_HSA)

   ! Assign trial configuration to SA-sampled config
   DO J2=1,NATOMS
      Xhsa(J2) = CONFIG(3*(J2-1)+1)
      Yhsa(J2) = CONFIG(3*(J2-1)+2)
      Zhsa(J2) = CONFIG(3*(J2-1)+3)
   ENDDO

   RETURN 
   end subroutine RESERVOIR_SAMPLE_HSA
