C   Copyright (C) 1992  N.M. Maclaren
C   Copyright (C) 1992  The University of Cambridge

C   This software may be reproduced and used freely, provided that all
C   users of it agree that the copyright holders are not liable for any
C   damage or injury caused by use of this software and that this
C   condition is passed onto all subsequent recipients of the software,
C   whether modified or not.



        SUBROUTINE SDPRND_UNIVERSAL (ISEED1)
        DOUBLE PRECISION XMOD, YMOD, POLY1(101), OTHER1, OFFSET1, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED1, INDEX1, IX, IY, IZ, I
        LOGICAL INITAL1
        SAVE INITAL1
        COMMON /RANDDP1/ POLY1, OTHER1, OFFSET1, INDEX1
        DATA INITAL1/.TRUE./
C
C   ISEED should be set to an integer between 0 and 9999 inclusive;
C   a value of 0 will initialise the generator only if it has not
C   already been done.
C
        IF (INITAL1 .OR. ISEED1 .NE. 0) THEN
            INITAL1 = .FALSE.
        ELSE
            RETURN
        END IF
C
C   INDEX must be initialised to an integer between 1 and 101
C   inclusive, POLY(1...N) to integers between 0 and 1000009710
C   inclusive (not all 0), and OTHER to a non-negative proper fraction
C   with denominator 33554432.  It uses the Wichmann-Hill generator to
C   do this.
C
        IX = MOD(ABS(ISEED1),10000)+1
        IY = 2*IX+1
        IZ = 3*IX+1
        DO 10 I = -10,101
            IF (I .GE. 1) POLY1(I) = AINT(XMOD*X)
            IX = MOD(171*IX,30269)
            IY = MOD(172*IY,30307)
            IZ = MOD(170*IZ,30323)
            X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+
     1        DBLE(IZ)/30323.0D0,1.0D0)
  10    CONTINUE
        OTHER1 = AINT(YMOD*X)/YMOD
        OFFSET1 = 1.0D0/YMOD
        INDEX1 = 1
        END

        DOUBLE PRECISION FUNCTION DPRAND_UNIVERSAL()
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY1(101),
     1    OTHER1, OFFSET1, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0,
     1    XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0,
     2    TINY = 1.0D-17)
        INTEGER INDEX1, N
        LOGICAL INITAL1
        SAVE INITAL1
        COMMON /RANDDP1/ POLY1, OTHER1, OFFSET1, INDEX1
        DATA INITAL1/.TRUE./
C
C   This returns a uniform (0,1) random number, with extremely good
C   uniformity properties.  It assumes that double precision provides
C   at least 33 bits of accuracy, and uses a power of two base.
C
        IF (INITAL1) THEN
            CALL SDPRND_UNIVERSAL (0)
            INITAL1 = .FALSE.
        END IF
C
C   See [Knuth] for why this implements the algorithm described in
C   the paper.  Note that this code is tuned for machines with fast
C   double precision, but slow multiply and divide; many, many other
C   options are possible.
C
        N = INDEX1-64
        IF (N .LE. 0) N = N+101
        X = POLY1(INDEX1)+POLY1(INDEX1)
        X = XMOD4-POLY1(N)-POLY1(N)-X-X-POLY1(INDEX1)
        IF (X .LT. 0.0D0) THEN
            IF (X .LT. -XMOD) X = X+XMOD2
            IF (X .LT. 0.0D0) X = X+XMOD
        ELSE
            IF (X .GE. XMOD2) THEN
                X = X-XMOD2
                IF (X .GE. XMOD) X = X-XMOD
            END IF
            IF (X .GE. XMOD) X = X-XMOD
        END IF
        POLY1(INDEX1) = X
        INDEX1 = INDEX1+1
        IF (INDEX1 .GT. 101) INDEX1 = INDEX1-101
C
C   Add in the second generator modulo 1, and force to be non-zero.
C   The restricted ranges largely cancel themselves out.
C
   10   Y = 37.0D0*OTHER1+OFFSET1
        OTHER1 = Y-AINT(Y)
        IF (OTHER1 .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER1
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND_UNIVERSAL = X+TINY
        END
