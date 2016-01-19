! 
! Three-dimensional Anti-periodic XY model. ch558
!
SUBROUTINE Energy_3d_APBC(THETA,GRAD,ENERGY,GTEST,SECT)
  USE COMMONS, ONLY : NONEDAPBC, NATOMS, XYPHI ! XYPHI has dim 3N*N*N
  IMPLICIT NONE
  INTEGER N, i1, i2,i3
  DOUBLE PRECISION, dimension((3*NATOMS)*(3*NATOMS)*(3*NATOMS)) :: theta, GRAD
  DOUBLE PRECISION :: Energy
  LOGICAL GTEST,SECT
    
  N=NONEDAPBC
  Energy=0
  i3=0

  DO WHILE ( i3.LT.(N*N*N))
     i2=0
     DO WHILE ( i2.LT.(N*N))
        i1=0
        DO WHILE(i1.LT.(N))
           
           IF(i1.EQ.(N-1)) THEN
              Energy = Energy +  ( cos(xyphi(i1+i2+i3+1) - theta(i2+i3+1) - theta(i1+i2+13+1)));
           ELSE
              Energy = Energy + (cos(xyphi(i1+i2+i3+1) + theta(i1+i2+i3+2) - theta(i1+i2+i3+1)));
           ENDIF
           
           IF(i2.EQ.((N*N) - N)) THEN
              Energy = Energy + ( cos(xyphi(i1+i2+i3+(N*N*N)+1) - theta(i1+i3+1) - theta(i1+i2+i3+1)));
           ELSE
              Energy = Energy + (cos(xyphi(i1+i2+i3+(N*N*N)+1) + theta(i1+i2+i3+N+1)- theta(i1+i2+i3+1)));
           ENDIF
           
           IF (i3.EQ.(N*N*N-N*N)) THEN
              Energy = Energy + ( cos(xyphi(i1+i2+i3+2*(N*N*N)+1) - theta(i1+i2+1) - theta(i1+i2+i3+1)))
           ELSE
              Energy = Energy + ( cos(xyphi(i1+i2+i3+2*N*N*N+1) + theta(i1+i2+i3+N*N+1) - theta(i1+i2+i3+1)))
           ENDIF

           i1=i1+1
        ENDDO
        i2 = i2 + N
     END DO
     i3 = i3 + (N*N)
  ENDDO
  
  Energy= 3 - (Energy/(N**3)) 


  IF (.NOT.GTEST) RETURN
  i3=0

  DO WHILE (i3.LT.(N*N*N))
     i2=0
     DO WHILE (i2.LT.(N*N))
        i1=0
        DO WHILE (i1.LT.N)
           


           if (i1.EQ.0 .AND. i2.EQ.0 .AND. i3.EQ.0) THEN 
     
              grad(1+0)= -sin(xyphi(1+0)+theta(1+1)-theta(1+0))-sin(xyphi(1+N-1+i2+i3)-theta(1+i1+i2+i3)-theta(1+N-1+i2+i3)) &
-sin(xyphi(1+i1+i2+i3+(N*N*N))+theta(1+i1+i2+i3+N)-theta(1+i1+i2+i3))& 
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))+theta(1+i1+i2+N*N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+i2+N*N*N-N*N+(2*N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+i2+N*N*N-N*N))
              
              
              
           else if (i1.EQ.0 .AND. i2.EQ.0 .AND. i3.EQ.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3)= -sin(xyphi(1+i1+i2+i3)+theta(1+i1+i2+i3+1)-theta(1+i1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i1+i2+i3)-theta(1+N-1+i2+i3)) &
-sin(xyphi(1+i1+i2+i3+(N*N*N))+theta(1+i1+i2+i3+N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-N*N+(2*N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-N*N))
              
           else if (i1.EQ.0 .AND. i2.EQ.(N*N-N) .AND. i3.EQ.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3)= -sin(xyphi(1+i1+i2+i3)+theta(1+i1+i2+i3+1)-theta(1+i1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i1+i2+i3)-theta(1+N-1+i2+i3))&
+sin(xyphi(1+i1+i2-N+i3+(N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2-N+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-N*N+(2*N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-N*N))
              
           else if (i1.EQ.0 .AND. i2.EQ.(N*N-N) .AND. i3.EQ.0) THEN 
              
              grad(1+i1+i2+i3)= -sin(xyphi(1+i1+i2+i3)+theta(1+i1+i2+i3+1)-theta(1+i1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i1+i2+i3)-theta(1+N-1+i2+i3))&
+sin(xyphi(1+i1+i2-N+i3+(N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2-N+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))+theta(1+i1+i2+i3+N*N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+i2+N*N*N-N*N+(2*N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+i2+N*N*N-N*N))
              
           else if (i1.EQ.(N-1) .AND. i2.EQ.0 .AND. i3.EQ.0) THEN 
              
              grad(1+i1+i2+i3)= sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3))&
-sin(xyphi(1+i1+i2+i3+(N*N*N))+theta(1+i1+i2+i3+N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))+theta(1+i1+i2+N*N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+i2+N*N*N-N*N+(2*N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+i2+N*N*N-N*N))
              
           else if (i1.EQ.(N-1) .AND. i2.EQ.0 .AND. i3.EQ.(N*N*N-N*N)) THEN 
     
              grad(1+i1+i2+i3)=     sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3)) &
-sin(xyphi(1+i1+i2+i3+(N*N*N))+theta(1+i1+i2+i3+N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-N*N+(2*N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-N*N))
              
           else if (i1.EQ.(N-1) .AND. i2.EQ.(N*N-N) .AND. i3.EQ.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3)=    sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3)) &
+sin(xyphi(1+i1+i2-N+i3+(N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2-N+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-N*N+(2*N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-N*N))
              
           else if (i1.EQ.(N-1) .AND. i2.EQ.(N*N-N) .AND. i3.EQ.0) THEN 
              
              grad(1+i1+i2+i3)=   sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3))&
+sin(xyphi(1+i1+i2-N+i3+(N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2-N+i3))&
-sin(xyphi(1+i1+(N*N)-N+i3+(N*N*N))-theta(1+i1+i3)-theta(1+i1+N*N-N+i3))&
-sin(xyphi(1+i1+i2+i3+(2*N*N*N))+theta(1+i1+i2+i3+N*N)-theta(1+i1+i2+i3))&
-sin(xyphi(1+i1+i2+N*N*N-N*N+(2*N*N*N))-theta(1+i1+i2+i3)-theta(1+i1+i2+N*N*N-N*N))
              
           else if(i3.EQ.0 .AND. i1.EQ.(N-1) .AND. i2.gt.0 .AND. i2.lt.(N*N-N)) THEN 
    
              grad(1+i1+i2) = -sin(xyphi(1+N-1+i2) - theta(1+i2) - theta(1+N-1+i2))&
+sin(xyphi(1+N-2+i2)+theta(1+N-1+i2)-theta(1+N-2+i2)) &
- sin(xyphi(1+N-1+i2+(N*N*N))+theta(1+N-1+i2+N)-theta(1+N-1+i2)) &
+ sin(xyphi(1+N-1+i2-N+(N*N*N))+theta(1+N-1+i2)-theta(1+N-1+i2-N)) &
- sin(xyphi(1+i1+i2+(2*N*N*N))+theta(1+i1+i2+(N*N))-theta(1+i1+i2)) &
- sin(xyphi(1+i1+i2+(N*N*N)-(N*N)+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+(N*N*N)-(N*N)))
              

              
           else if(i3.EQ.0 .AND. i1.EQ.0 .AND. i2.gt.0 .AND. i2.lt.(N*N-N)) THEN 
              
              grad(1+i1+i2) = -sin(xyphi(1+i2)+theta(1+1+i2)-theta(1+i2)) &
- sin(xyphi(1+N-1+i2)-theta(1+i2)-theta(1+N-1+i2))&
+sin(xyphi(1+i2-N+(N*N*N))+theta(1+i2)-theta(1+i2-N)) &
- sin(xyphi(1+i2+(N*N*N))+theta(1+i2+N) - theta(1+i2)) &
- sin(xyphi(1+i1+i2+(2*N*N*N))+theta(1+i1+i2+(N*N))-theta(1+i1+i2)) &
- sin(xyphi(1+i1+i2+(N*N*N)-(N*N)+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+(N*N*N)-(N*N)))
              
           else if(i3.EQ.0 .AND. i2.EQ.(N*N-N) .AND. i1.gt.0 .AND. i1.lt.(N-1)) THEN 
              
              grad(1+i1+i2) = -sin(xyphi(1+i1+i2) + theta(1+i1+1+i2)-theta(1+i1+i2)) &
+ sin(xyphi(1+i1-1+i2) + theta(1+i1+i2) - theta(1+i1-1+i2)) &
+sin(xyphi(1+i1+i2-N + (N*N*N))+theta(1+i1+i2) - theta(1+i1+i2-N)) &
- sin(xyphi(1+i1+i2+(N*N*N)) - theta(1+i1)-theta(1+i1+i2)) &
- sin(xyphi(1+i1+i2+(2*N*N*N))+theta(1+i1+i2+(N*N))-theta(1+i1+i2)) &
- sin(xyphi(1+i1+i2+(N*N*N)-(N*N)+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+(N*N*N)-(N*N)))
              
           else if(i3.EQ.0 .AND. i2.EQ.0 .AND. i1.gt.0 .AND. i1.lt.(N-1)) THEN 
              
              grad(1+i1+i2)=  -sin(xyphi(1+i1)+theta(1+i1+1)-theta(1+i1)) &
+ sin(xyphi(1+i1-1)+theta(1+i1)-theta(1+i1-1))&
-sin(xyphi(1+i1+(N*N*N))+theta(1+i1+N)-theta(1+i1))&
-sin(xyphi(1+i1+(N*N)-N +(N*N*N))-theta(1+i1)-theta(1+i1+(N*N)-N)) &
- sin(xyphi(1+i1+i2+(2*N*N*N))+theta(1+i1+i2+(N*N))-theta(1+i1+i2)) &
- sin(xyphi(1+i1+i2+(N*N*N)-(N*N)+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+(N*N*N)-(N*N)))
              

              
           else if(i3.EQ.0 .AND. i1.lt.(N-1) .AND. i1.gt.0 .AND. i2.lt.(N*N-N) .AND. 0.lt.i2) THEN 
              grad(1+i1+i2)= -sin(xyphi(1+i1+i2)+theta(1+i1+1+i2)-theta(1+i1+i2))&
+ sin(xyphi(1+i1-1+i2)+theta(1+i1+i2)-theta(1+i1-1+i2)) &
- sin(xyphi(1+i1+i2+(N*N*N))+theta(1+i1+i2+N)-theta(1+i1+i2)) &
+ sin(xyphi(1+i1+i2-N +(N*N*N))+theta(1+i1+i2)-theta(1+i1+i2-N)) &
- sin(xyphi(1+i1+i2+(2*N*N*N))+theta(1+i1+i2+(N*N))-theta(1+i1+i2)) &
- sin(xyphi(1+i1+i2+(N*N*N)-(N*N)+(2*N*N*N))-theta(1+i1+i2)-theta(1+i1+i2+(N*N*N)-(N*N)))
              
           else if(i3.EQ.(N*N*N-N*N) .AND. i1.EQ.(N-1) .AND. i2.gt.0 .AND. i2.lt.(N*N-N)) THEN 
              
              grad(1+i1+i2+i3) = -sin(xyphi(1+N-1+i2+i3) - theta(1+i2+i3) - theta(1+N-1+i2+i3))&
+sin(xyphi(1+N-2+i2+i3)+theta(1+N-1+i2+i3)-theta(1+N-2+i2+i3)) &
- sin(xyphi(1+N-1+i2+(N*N*N)+i3)+theta(1+N-1+i2+N+i3)-theta(1+N-1+i2+i3)) &
+ sin(xyphi(1+N-1+i2-N+(N*N*N)+i3)+theta(1+N-1+i2+i3)-theta(1+N-1+i2-N+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) - theta(1+i1+i2) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
    
              
           else if(i3.EQ.(N*N*N-N*N) .AND. i1.EQ.0 .AND. i2.gt.0 .AND. i2.lt.(N*N-N)) THEN 
              grad(1+i1+i2+i3) = -sin(xyphi(1+i2+i3)+theta(1+1+i2+i3)-theta(1+i2+i3)) &
- sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3))&
+sin(xyphi(1+i2-N+(N*N*N)+i3)+theta(1+i2+i3)-theta(1+i2-N+i3)) &
- sin(xyphi(1+i2+(N*N*N)+i3)+theta(1+i2+N+i3) - theta(1+i2+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) - theta(1+i1+i2) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
              
           else if(i3.EQ.(N*N*N-N*N) .AND. i2.EQ.(N*N-N) .AND. i1.gt.0 .AND. i1.lt.(N-1)) THEN 
              
              grad(1+i1+i2+i3) = -sin(xyphi(1+i1+i2+i3) + theta(1+i1+1+i2+i3)-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1-1+i2+i3) + theta(1+i1+i2+i3) - theta(1+i1-1+i2+i3)) &
+sin(xyphi(1+i1+i2-N + (N*N*N)+i3)+theta(1+i1+i2+i3) - theta(1+i1+i2-N+i3)) &
- sin(xyphi(1+i1+i2+(N*N*N)+i3) - theta(1+i1+i3)-theta(1+i1+i2+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) - theta(1+i1+i2) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
              
              
           else if(i3.EQ.(N*N*N-N*N) .AND. i2.EQ.0 .AND. i1.gt.0 .AND. i1.lt.(N-1)) THEN 
    
              grad(1+i1+i2+i3)=  -sin(xyphi(1+i1+i3)+theta(1+i1+1+i3)-theta(1+i1+i3)) &
+ sin(xyphi(1+i1-1+i3)+theta(1+i1+i3)-theta(1+i1-1+i3))&
-sin(xyphi(1+i1+(N*N*N)+i3)+theta(1+i1+N+i3)-theta(1+i1+i3))&
-sin(xyphi(1+i1+(N*N)-N +(N*N*N)+i3)-theta(1+i1+i3)-theta(1+i1+(N*N)-N+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) - theta(1+i1+i2) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
      
              
              
           else if(i3.EQ.(N*N*N-N*N) .AND. i1.lt.(N-1) .AND. i1.gt.0 .AND. i2.lt.(N*N-N) .AND. 0.lt.i2) THEN 
              
              grad(1+i1+i2+i3)= -sin(xyphi(1+i1+i2+i3)+theta(1+i1+1+i2+i3)-theta(1+i1+i2+i3))+ &
sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3)) &
- sin(xyphi(1+i1+i2+(N*N*N)+i3)+theta(1+i1+i2+N+i3)-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1+i2-N +(N*N*N)+i3)+theta(1+i1+i2+i3)-theta(1+i1+i2-N+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) - theta(1+i1+i2) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
              
    
           else if(i2.EQ.0 .AND. i1.EQ.(N-1) .AND. i3.gt.0 .AND. i3.lt.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3) = -sin(xyphi(1+N-1+i3) - theta(1+i3) - theta(1+N-1+i3))&
+sin(xyphi(1+N-2+i3)+theta(1+N-1+i3)-theta(1+N-2+i3))&
- sin(xyphi(1+N-1+i3+(N*N*N))+theta(1+N-1+i3+N)-theta(1+N-1+i3)) &
- sin(xyphi(1+N-1+i3+N*N-N+(N*N*N))-theta(1+N-1+i3)-theta(1+N-1+i3+N*N-N)) &
- sin(xyphi(1+i1+i3+(2*N*N*N))+theta(1+i1+i3+(N*N))-theta(1+i1+i3)) &
+ sin(xyphi(1+i1+i3-(N*N)+(2*N*N*N))+theta(1+i1+i3)-theta(1+i1+i3-(N*N)))
              
           else if(i2.EQ.0 .AND. i1.EQ.0 .AND. i3.gt.0 .AND. i3.lt.(N*N*N-N*N)) THEN 
    
              grad(1+i1+i2+i3) = -sin(xyphi(1+i3)+theta(1+1+i3)-theta(1+i3)) &
- sin(xyphi(1+N-1+i3)-theta(1+i3)-theta(1+N-1+i3))&
- sin(xyphi(1+i1+i3+(N*N*N))+theta(1+i1+i3+N)-theta(1+i1+i3)) &
- sin(xyphi(1+i1+i3+N*N-N+(N*N*N))-theta(1+i1+i3)-theta(1+i1+i3+N*N-N)) &
 - sin(xyphi(1+i1+i3+(2*N*N*N))+theta(1+i1+i3+(N*N))-theta(1+i1+i3)) &
+ sin(xyphi(1+i1+i3-(N*N)+(2*N*N*N))+theta(1+i1+i3)-theta(1+i1+i3-(N*N)))
      
           else if(i2.EQ.0 .AND. i1.lt.(N-1) .AND. i1.gt.0 .AND. i3.lt.(N*N*N-N*N) .AND. 0.lt.i3) THEN 
              
              grad(1+i1+i2+i3)= -sin(xyphi(1+i1+i3)+theta(1+i1+1+i3)-theta(1+i1+i3))&
+ sin(xyphi(1+i1-1+i3)+theta(1+i1+i3)-theta(1+i1-1+i3)) &
 - sin(xyphi(1+i1+i3+(N*N*N))+theta(1+i1+i3+N)-theta(1+i1+i3))&
 - sin(xyphi(1+i1+i3+N*N-N+(N*N*N))-theta(1+i1+i3)-theta(1+i1+i3+N*N-N)) &
 - sin(xyphi(1+i1+i3+(2*N*N*N))+theta(1+i1+i3+(N*N))-theta(1+i1+i3)) &
+ sin(xyphi(1+i1+i3-(N*N)+(2*N*N*N))+theta(1+i1+i3)-theta(1+i1+i3-(N*N)))
              
              
           else if (i2.EQ.(N*N-N) .AND. i1.EQ.(N-1) .AND. i3.gt.0 .AND. i3.lt.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3) = -sin(xyphi(1+N-1+i2+i3) - theta(1+i2+i3) - theta(1+N-1+i2+i3))&
+sin(xyphi(1+N-2+i2+i3)+theta(1+N-1+i2+i3)-theta(1+N-2+i2+i3)) &
- sin(xyphi(1+N-1+i2+(N*N*N)+i3)-theta(1+N-1+i3)-theta(1+N-1+i2+i3))&
 + sin(xyphi(1+N-1+i2-N+(N*N*N)+i3)+theta(1+N-1+i2+i3)-theta(1+N-1+i2-N+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) + theta(1+i1+i2+i3+(N*N)) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
      
           else if(i2.EQ.(N*N-N) .AND. i1.EQ.0 .AND. i3.gt.0 .AND. i3.lt.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3) = -sin(xyphi(1+i2+i3)+theta(1+1+i2+i3)-theta(1+i2+i3)) &
- sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3))&
-sin(xyphi(1+i1+i2+(N*N*N)+i3)-theta(1+i1+i3)-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1+i2-N+(N*N*N)+i3)+theta(1+i1+i2+i3) - theta(1+i1+i2-N+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) + theta(1+i1+i2+i3+(N*N)) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
              
           else if(i2.EQ.(N*N-N) .AND. i1.lt.(N-1) .AND. i1.gt.0 .AND. i3.lt.(N*N*N-N*N) .AND. 0.lt.i3) THEN 
              grad(1+i1+i2+i3)= -sin(xyphi(1+i1+i2+i3)+theta(1+i1+1+i2+i3)-theta(1+i1+i2+i3))&
+ sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3)) &
-sin(xyphi(1+i1+i2+(N*N*N)+i3)-theta(1+i1+i3)-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1+i2-N+(N*N*N)+i3)+theta(1+i1+i2+i3) - theta(1+i1+i2-N+i3)) &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N)) + theta(1+i1+i2+i3+(N*N)) - theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N)) + theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
              
              
           else if(i1.EQ.0 .AND. i2.lt.(N*N-N) .AND. 0.lt.i2 .AND. 0.lt.i3 .AND. i3.lt.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3)=  -sin(xyphi(1+i2+i3) + theta(1+1+i2+i3) - theta(1+i2+i3))&
-sin(xyphi(1+N-1+i2+i3)-theta(1+i2+i3)-theta(1+N-1+i2+i3)) & 
- sin(xyphi(1+i2+i3+(N*N*N))+theta(1+i2+i3+N)-theta(1+i2+i3))&
 + sin(xyphi(1+i2+i3-N+(N*N*N))+theta(1+i2+i3)-theta(1+i2+i3-N))  &
- sin(xyphi(1+i2+i3+(2*N*N*N))+theta(1+i2+i3+(N*N))-theta(1+i2+i3))&
 + sin(xyphi(1+i2+i3-(N*N)+(2*N*N*N))+theta(1+i2+i3)-theta(1+i2+i3-(N*N)))
              
           else if(i1.EQ.(N-1) .AND. i2.lt.(N*N-N) .AND. 0.lt.i2 .AND. 0.lt.i3 .AND. i3.lt.(N*N*N-N*N)) THEN 
              
              grad(1+i1+i2+i3)=  -sin(xyphi(1+i1+i2+i3) - theta(1+i2+i3) - theta(1+i1+i2+i3))&
+sin(xyphi(1+N-2+i2+i3)+theta(1+i1+i2+i3)-theta(1+N-2+i2+i3)) &
- sin(xyphi(1+i1+i2+i3+(N*N*N))+theta(1+i1+i2+i3+N)-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1+i2+i3-N+(N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-N))  &
- sin(xyphi(1+i1+i2+i3+(2*N*N*N))+theta(1+i1+i2+i3+(N*N))-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1+i2+i3-(N*N)+(2*N*N*N))+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-(N*N)))
              
           else 
              
              grad(1+i1+i2+i3) = -sin(xyphi(1+i1+i2+i3)+theta(1+i1+1+i2+i3)-theta(1+i1+i2+i3)) &
+ sin(xyphi(1+i1-1+i2+i3)+theta(1+i1+i2+i3)-theta(1+i1-1+i2+i3))&
-sin(xyphi(1+i1+i2+i3+N*N*N)+theta(1+i1+i2+N+i3)-theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2-N+i3+N*N*N)+theta(1+i1+i2+i3)-theta(1+i1+i2-N+i3))&
-sin(xyphi(1+i1+i2+i3+2*N*N*N)+theta(1+i1+i2+i3+N*N)-theta(1+i1+i2+i3))&
+sin(xyphi(1+i1+i2+i3-N*N+2*N*N*N)+theta(1+i1+i2+i3)-theta(1+i1+i2+i3-N*N))
      
              
              
           END IF
           i1=i1+1
        ENDDO
        i2=i2+N
     ENDDO
     i3=i3+N*N
  ENDDO

END SUBROUTINE ENERGY_3D_APBC
