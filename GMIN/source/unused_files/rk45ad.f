C
C PAGE 326-328: NUMERICAL MATHEMATICS AND COMPUTING, CHENEY/KINCAID, 1985
C
C FILE: RK45AD.FOR
C
C ADAPTIVE SCHEME BASED ON RUNGE-KUTTA-FEHLBERG METHOD (RK45AD,RK45,F)
C
      EXTERNAL F  
      TA = 1.0
      XA = 2.0
      H = 7.8125E-3 
      TB = 1.5625 
      ITMAX = 100 
      EMIN = 1.0E-8 
      EMAX = 1.0E-4 
      HMIN = 1.0E-6 
      HMAX = 1.0
      CALL RK45AD(F,TA,XA,H,TB,ITMAX,EMIN,EMAX,HMIN,HMAX,IFLAG)     
      PRINT 3,H,TA,XA
    3 FORMAT(5X,E10.3,2(5X,E22.14))
      PRINT 5,IFLAG 
    5 FORMAT('0 IFLAG =',I2)
      STOP
      END 
  
      FUNCTION F(T,X) 
      F = 2.0+(X-T-1.0)**2  
      RETURN
      END 
  
      SUBROUTINE RK45AD(F,T,X,H,TB,ITMAX,EMIN,EMAX,HMIN,HMAX,IFLAG) 
      EXTERNAL  F 
      DATA  EPSI/1.0E-6/   
      IFLAG = 2
      NSTEP = 0   
      PRINT 4     
      PRINT 3,H,T,X 
   2  IF(ABS(H) .LT. HMIN)  H = HMIN
      IF(ABS(H) .GT. HMAX)  H = HMAX
      DT = ABS(TB - T)      
      IF(DT .LE. ABS(H))  THEN
        IFLAG = 0 
        IF(DT .LE. EPSI*AMAX1(ABS(TB),ABS(T)))  RETURN    
        H = SIGN(DT,H)      
      END IF      
      XSAVE = X   
      TSAVE = T   
      CALL RK45(F,T,X,H,EST)
      NSTEP = NSTEP + 1     
      PRINT 3,H,T,X,EST     
      IF(IFLAG .EQ. 0)  RETURN
      IF(NSTEP .GT. ITMAX)  THEN
        IFLAG = 1 
        RETURN
      END IF      
      IF((EMIN .LE. EST) .AND. (EST .LE. EMAX))  GO TO 2  
      IF(EST .LT. EMIN)  THEN 
        H = 2.0*H 
      ELSE
        H = 0.5*H 
      END IF      
      X = XSAVE   
      T = TSAVE   
      GO TO 2 
   3  FORMAT(5X,E10.3,2(5X,E22.14),5X,E10.3)    
   4  FORMAT(9X,'H',22X,'T',26X,'X',18X,'EST')  
      END 
  
      SUBROUTINE RK45(F,T,X,H,EST)    
      DATA  C21,C31,C32, C41,C42,C43, C51,C52,C53,C54,    
     A      C61,C62,C63,C64,C65, A1,A3,A4,A5, B1,B3,B4,B5, B6, C40  
     B  /0.25,0.09375,0.28125,
     C  0.87938097405553,-3.2771961766045,3.3208921256258,
     D  2.0324074074074,-8.0,7.1734892787524,-0.20589668615984,     
     E  -0.2962962962963,2.0,-1.3816764132554,0.45297270955166,-0.275,
     F  0.11574074074074,0.54892787524366,0.5353313840156,-0.2,     
     G  0.11851851851852,0.51898635477583,0.50613149034201,-0.18,   
     H  0.036363636363636, 0.92307692307692/    
      F1 = H*F(T,X) 
      F2 = H*F(T+ 0.25*H,X + C21*F1)  
      F3 = H*F(T+0.375*H,X + C31*F1 + C32*F2)   
      F4 = H*F(T+C40*H,X + C41*F1 + C42*F2 + C43*F3)  
      F5 = H*F(T+H    ,X + C51*F1 + C52*F2 + C53*F3 + C54*F4)       
      F6 = H*F(T+0.5*H,X + C61*F1 + C62*F2 + C63*F3 + C64*F4 + C65*F5)
      X5 = X + B1*F1 + B3*F3 + B4*F4 + B5*F5 + B6*F6      
      X  = X + A1*F1 + A3*F3 + A4*F4 + A5*F5    
      T = T + H   
      EST = ABS(X - X5)     
      RETURN      
      END 
