! 
! Two-dimensional periodic XY model. ch558
!
SUBROUTINE Energy_2d_PBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE COMMONS, ONLY : NONEDAPBC, NATOMS, XYPHI ! XYPHI has dim 2N*N
  IMPLICIT NONE
  INTEGER N, i1, i2
  DOUBLE PRECISION, dimension((3*NATOMS)*(3*NATOMS)) :: theta, GRAD
  DOUBLE PRECISION :: Energy
  LOGICAL GTEST,SECT
    
  n=NONEDAPBC
  theta(N*N)=0
  Energy=0
  i2=0
  DO WHILE ( i2.LT.(N*N))
     i1=0
     DO WHILE(i1.LT.(N))
        IF(i1.EQ.(N-1)) THEN
           Energy = Energy +  (cos(xyphi(i1+i2+1) + theta(i2+1) - theta(i1+i2+1)));
        ELSE
           Energy = Energy + ( cos(xyphi(i1+i2+1) + theta(i1+i2+2) - theta(i1+i2+1)));
        ENDIF
        
        IF(i2.EQ.((N*N) - N)) THEN
           Energy = Energy + ( cos(xyphi(i1+i2+(N*N)+1) + theta(i1+1) - theta(i1+i2+1)));
        ELSE
           Energy = Energy + (cos(xyphi(i1+i2+(N*N)+1) + theta(i1+i2+N+1)- theta(i1+i2+1)));
        ENDIF
        i1=i1+1
     ENDDO
     i2 = i2 + N
  END DO


  Energy = 2 - (Energy/(N*N))

  IF (.NOT.GTEST) RETURN
  i2=0
  
  DO WHILE(i2.LT.(N*N))
     i1=0
     DO WHILE(i1.LT.N)
        
        IF(i1.EQ.(N-1) .AND. i2.GT.0 .AND. i2.LT.(N*N-N)) THEN
           grad(i1+i2+1) = -sin(xyphi(N+i2)+theta(i2+1)-theta(N+i2)) + sin(xyphi(N-1+i2)+theta(N+i2)-theta(N-1+i2)) &
- sin( xyphi( N + i2 + N**2) + theta(N+i2+N) - theta(N+i2) ) + sin(xyphi(N+i2-N+N*N)+theta(N+i2)-theta(N+i2-N))
        ELSE IF(i1.EQ.0 .AND. i2.EQ.0) THEN
           grad(i1+i2+1)= -sin(xyphi(1)+theta(2)-theta(1)) + sin(xyphi(N)+theta(1)-theta(N)) &
- sin(xyphi(N*N+1)+theta(N+1)-theta(1)) + sin(xyphi(N*N-N+N*N+1)+theta(1)-theta(N*N-N+1))
        ELSE IF(i1.EQ.0 .AND. i2.GT.0 .AND. i2.LT.(N*N-N)) THEN
           grad(i1+i2+1) = -sin(xyphi(i2+1)+theta(2+i2)-theta(i2+1)) + sin(xyphi(N+i2)+theta(i2+1)-theta(N+i2)) &
+sin(xyphi(i2-N+N*N+1)+theta(i2+1)-theta(i2-N+1)) - sin(xyphi(i2+N*N+1)+theta(i2+N+1) - theta(i2+1))
        ELSE IF(i2.EQ.(N*N-N) .AND. i1.GT.0 .AND. i1.LT.(N-1)) THEN
           grad(i1+i2+1) = -sin(xyphi(i1+i2+1) + theta(i1+2+i2)-theta(i1+i2+1)) + sin(xyphi(i1+i2) + & 
theta(i1+i2+1) - theta(i1+i2)) + sin(xyphi(i1+i2-N +N*N+1)+theta(i1+i2+1)- theta(i1+i2-N+1)) &
-sin(xyphi(i1+i2+N*N+1) + theta(i1+1)-theta(i1+i2+1))
        ELSE IF(i2.EQ.0 .AND. i1.GT.0 .AND. i1.LT.(N-1)) THEN
           grad(i1+i2+1)=  -sin(xyphi(i1+1)+theta(i1+2)-theta(i1+1)) + sin(xyphi(i1)+theta(i1+1)-theta(i1)) &
-sin(xyphi(i1+N*N+1)+theta(i1+N+1)-theta(i1+1))+sin(xyphi(i1+N*N-N +N*N+1)+theta(i1+1)-theta(i1+N*N-N+1))
        ELSE IF(i1.EQ.(N-1) .AND. i2.EQ.0) THEN
           grad(i1+i2+1) =  -sin(xyphi(N)+theta(1)-theta(N)) + sin(xyphi(N-1) + theta(N)-theta(N-1)) & 
- sin(xyphi(N +N*N)+theta(N+N) - theta(N))+sin(xyphi(N*N +(N*N))+theta(N)-theta(N*N))
        ELSE IF(i1.EQ.0 .AND. i2.EQ.(N*N-N)) THEN
           grad(i1+i2+1) = -sin(xyphi(N*N-N+1)+theta(2+N*N-N)-theta(N*N-N+1)) & 
+ sin(xyphi(N + N*N-N)+theta(N*N-N+1)-theta(N+N*N-N) ) &
-sin(xyphi(N*N-N+N*N+1)+theta(1)-theta(N*N-N+1)) &
+ sin(xyphi(N*N-N-N+N*N+1)+theta(N*N-N+1) - theta(N*N-N-N+1))
        ELSE
           grad(i1+i2+1)= -sin(xyphi(i1+i2+1)+theta(i1+2+i2)-theta(i1+i2+1))+ sin(xyphi(i1+i2)+theta(i1+i2+1)-theta(i1+i2)) &
- sin( xyphi(i1+i2+N*N+1) + theta(i1+i2+N+1)-theta(i1+i2+1)) + sin(xyphi(i1+i2-N +N*N+1)+theta(i1+i2+1)-theta(i1+i2-N+1))
        ENDIF
        i1=i1+1
     ENDDO
     i2=i2+N
  ENDDO

grad(N*N)=0

END SUBROUTINE Energy_2d_PBC



