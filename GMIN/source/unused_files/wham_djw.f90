!       WHAM program for histograms P(T,V,Q1,Q2)

	PARAMETER(OUT=1E35)
	PARAMETER(NTMAX=32)          ! NTMAX is the maximum number of replicas
	PARAMETER(NVMAX=1000,NINI=5) ! NVMAX is the maximum number of bins
	PARAMETER(NQ1MAX=100,NQ2MAX=100)
	PARAMETER(KB=1.D0)

	INTEGER HISTO(NTMAX,NVMAX,NQ1MAX,NQ2MAX)
	INTEGER HISTV(NVMAX)

	DOUBLE PRECISION  HIST_T(NVMAX,NQ1MAX,NQ2MAX)

	DOUBLE PRECISION V,VMIN,VMAX,DV,S,XP,PIVOT,U,T(NTMAX)
	DIMENSION N(NTMAX,NVMAX),NPA(NVMAX),NMAX(NTMAX)
	DIMENSION JPMAX(NTMAX)
	DOUBLE PRECISION B(NTMAX,NVMAX),ETOT(NTMAX)
	DOUBLE PRECISION AM(NTMAX,NTMAX),BM(NTMAX),X(NTMAX),Y(NVMAX)
	DOUBLE PRECISION EMIN(NTMAX),EMAX(NTMAX),ENUL
	DOUBLE PRECISION XNPA(NVMAX),TMIN,TMAX,TL,DT
	DOUBLE PRECISION Z,CV,VMINNEW,VMAXNEW,KADM

	DOUBLE PRECISION KDIB(NVMAX)
	DOUBLE PRECISION ETOTMIN,ETOTMAX,DETOT,EETOT

	CHARACTER*50 FILENAME, ISTR

	OPEN(10,FILE='histo.para')
	READ(10,*) VMIN,VMAX,NV
	READ(10,*) Q1MIN,Q1MAX,NQ1
	READ(10,*) Q2MIN,Q2MAX,NQ2
	READ(10,*) NE
	DO I=1,NE
	   READ(10,*) T(I)
	ENDDO
	CLOSE(10)

	DV=(VMAX-VMIN)/(NV-1)
	DQ1=(Q1MAX-Q1MIN)/(NQ1-1)
	DQ2=(Q2MAX-Q2MIN)/(NQ2-1)

	DO I=1,NE     ! I is the number of temperatures
	   DO J=1,NV  ! J is the number of bins
	      N(I,J)=0
	      B(I,J)=0.0D0
	   ENDDO
	ENDDO

	DO J=1,NV
	   HISTV(J)=0
	ENDDO

	DO I=1,NE
           WRITE (ISTR, '(i10)') i-1 
           FILENAME="FreeEnStatHistE.HistO1.HistO2."//trim(adjustl(istr))
	   OPEN(UNIT=20,FILE=FILENAME,STATUS='old')
	   NMAX(I)=0
	   PRINT *,'T=',T(i)
	   
 0001	   READ(20,*,END=0009) J,IQ1,IQ2,NHIST
	   HISTO(I,J,IQ1,IQ2)=NHIST
	   N(I,J)=N(I,J)+HISTO(I,J,IQ1,IQ2)
	   HISTV(J)=HISTV(J)+HISTO(I,J,IQ1,IQ2)
	   IF (NMAX(I).LT.N(I,J)) THEN
	      NMAX(I)=N(I,J)
	      JPMAX(I)=J
	   ENDIF
	   GOTO 0001
 0009	   CONTINUE
	   PRINT *,'n(',jpmax(i),')=',nmax(i)
	   CLOSE(20)
	ENDDO
	
	DO I=1,NE
	   
	   PRINT *,I
	   
	   DO J=JPMAX(I)+1,NV
	      IF (N(I,J).LE.NINI) THEN
		 DO K=J,NV
		    N(I,K)=0
		 ENDDO
		EMAX(I)=VMIN+(J-1)*DV
		GOTO 475
	     ENDIF
	  ENDDO
 475	  DO J=JPMAX(I),1,-1
	     IF (N(I,J).LE.NINI) THEN
		DO K=J,1,-1
		   N(I,K)=0
		ENDDO
		EMIN(I)=VMIN+J*DV
		GOTO 476
	     ENDIF
	  ENDDO
	  
 476	  DO J=1,NV
	     B(I,J)=0.
	     V=VMIN+(J-1)*DV
	     IF ((V.GE.EMIN(I)).AND.(V.LE.EMAX(I))) THEN
		IF (N(I,J).NE.0) THEN
		   B(I,J)=LOG(1.D0*N(I,J))+V/T(I)
		ENDIF
	     ENDIF
	  ENDDO
	ENDDO

	PRINT *,'n(i,j) calculated.'

	OPEN(50,FILE='nij',status='unknown')
	DO I=1,NE
	   DO J=1,NV
	      WRITE(50,*) REAL(VMIN+(J-1)*DV),N(I,J)
	   ENDDO
	ENDDO
	CLOSE(50)

 1985	DO J=1,NV
	   NPA(J)=0
	   DO I=1,NE
	      NPA(J)=NPA(J)+N(I,J)
	   ENDDO
	   IF (NPA(J).EQ.0) NPA(J)=1
	   XNPA(J)=DFLOAT(NPA(J))
	ENDDO

	DO I=1,NE
	   BM(I)=0.
	   DO J=1,NV
	      IF (N(I,J).NE.0) BM(I)=BM(I)+B(I,J)*DFLOAT(N(I,J))
	      
	      XP=0.
	      DO K=1,NE
		 IF (N(K,J).NE.0) XP=XP+DFLOAT(N(K,J))*B(K,J)
	      ENDDO
	      BM(I)=BM(I)-XP*DFLOAT(N(I,J))/XNPA(J)
	   ENDDO
	ENDDO

	S=0.
	DO I=1,NE
	   S=S+BM(I)
	ENDDO
	PRINT *,'B calculated. Sum of terms is ',s

	DO I=1,NE

	   DO IP=1,NE

	      AM(I,IP)=0.
	      DO K=1,NV
		 AM(I,IP)=AM(I,IP)-DFLOAT(N(I,K))*DFLOAT(N(IP,K))/XNPA(K)
	      ENDDO
	      
	      IF (I.EQ.IP) THEN
		 DO K=1,NV
		    AM(I,IP)=AM(I,IP)+DFLOAT(N(I,K))
		 ENDDO
	      ENDIF
	      
	   ENDDO
	   
	ENDDO
	
	PRINT *,'A calculated'

	DO J=1,NE
	   S=0.
	   DO I=1,NE
	      S=S+AM(I,J)
	   ENDDO
	   PRINT *,'Column ',j,': ',s
	ENDDO

	DO I=1,NE-1

	   JPIVOT=0
	   PIVOT=0.
	   DO J=I,NE
	      U=AM(J,I)
	      IF (PIVOT.LT.ABS(U)) THEN
		 PIVOT=ABS(U)
		 JPIVOT=J
	      ENDIF
	   ENDDO
	   J=JPIVOT
	   U=AM(J,I)
	   
	   DO K=1,NE
	      S=AM(I,K)
	      AM(I,K)=AM(J,K)
	      AM(J,K)=S
	   ENDDO
	   S=BM(I)
	   BM(I)=BM(J)
	   BM(J)=S

	   DO K=I+1,NE
	      S=AM(K,I)	
	      AM(K,I)=0.
	      DO NK=I+1,NE
		 AM(K,NK)=AM(K,NK)-S*AM(I,NK)/U
	      ENDDO
	      BM(K)=BM(K)-S*BM(I)/U
	   ENDDO

	   DO K=1,I-1
	      S=AM(K,I)
	      AM(K,I)=0.
	      DO NK=I+1,NE
		 AM(K,NK)=AM(K,NK)-S*AM(I,NK)/U
	      ENDDO
	      BM(K)=BM(K)-S*BM(I)/U
	   ENDDO

	ENDDO

	AM(NE,NE)=1.
	BM(NE)=0.
	
	PRINT *,'System solved'

	X(NE)=BM(NE)/AM(NE,NE)

	DO I=NE-1,1,-1
	   S=BM(I)
	   X(I)=S/AM(I,I)
	ENDDO

	OPEN(17,FILE='solx',status='unknown')
	DO J=1,NE
	   WRITE(17,*) J,X(J)
	ENDDO
	CLOSE(17)

	DO J=1,NV
	   Y(J)=0.
	   DO I=1,NE
	      IF (N(I,J).NE.0) Y(J)=Y(J)+N(I,J)*(B(I,J)-X(I))
	   ENDDO
	   Y(J)=Y(J)/XNPA(J)
	ENDDO

	OPEN(20,FILE='S_wham',status='unknown')
	DO I=1,NV
	   IF (Y(I).NE.0) WRITE(20,*) VMIN+(I-1)*DV,Y(I),I
	ENDDO
	CLOSE(20)

	PRINT *,'Entropy calculated'

        IF (NQ1.EQ.1) RETURN 

	TMIN=T(1)
	TMAX=T(NE)
	NEXP=-1000000000
	NPOINTS=500
	DT=(TMAX-TMIN)/(NPOINTS-1)
	
	OPEN(11,FILE='U_wham',status='unknown')
	OPEN(12,FILE='Cv_wham',status='unknown')
	OPEN(13,FILE='Q1_wham',status='unknown')
	OPEN(14,FILE='Q2_wham',status='unknown')

	DO I=1,NV
	   IF (Y(I).NE.0) THEN
	      NPOPO=INT(Y(I)-(VMIN+(I-1)*DV)/TMIN)
	      IF (NPOPO.GT.NEXP) NEXP=NPOPO
	   ENDIF
	ENDDO

	DO NT=1,NPOINTS
	   TL=TMIN+(NT-1)*DT

	   PRINT *,'T=',Tl

	   DO J=1,NV
	      DO IQ1=1,NQ1
		 DO IQ2=1,NQ2
		    HIST_T(J,IQ1,IQ2)=0.D0
		 ENDDO
	      ENDDO
	   ENDDO
 2001	   Z=0.
	   DO J=1,NV
	      IF (Y(J).NE.0) THEN
		 BOLTZ=EXP(Y(J)-(VMIN+(J-1)*DV)/TL-NEXP*1.)
		 Z=Z+BOLTZ
	      ENDIF
	   ENDDO
           IF (Z.LT.1) THEN
              NEXP=NEXP-2
              GOTO 2001
           ENDIF
           IF (Z.GT.100) THEN
              NEXP=NEXP+2
              GOTO 2001
           ENDIF

           DO J=1,NV
              IF (Y(J).NE.0) THEN
                 BOLTZ=EXP(Y(J)-(VMIN+(J-1)*DV)/TL-NEXP*1.)
                 Z=Z+BOLTZ
		 DO IT=1,NE
		    DO IQ1=1,NQ1
		       DO IQ2=1,NQ2
			  HIST_T(J,IQ1,IQ2)=HIST_T(J,IQ1,IQ2) &
     &		    + HISTO(IT,J,IQ1,IQ2)*BOLTZ/HISTV(J)
		       ENDDO
		    ENDDO
		 ENDDO
	      ENDIF
	   ENDDO

	   XP=0
	   Z=0.
	   U=0.
	   S=0.
	   Q1AV=0.D0
	   Q2AV=0.D0
	   Q1AV2=0.D0
	   Q2AV2=0.D0
	   AVNUM=0.D0
	   DO J=1,NV
	      V=VMIN+(J-1)*DV
	      IF (Y(J).NE.0) THEN
		 XP=EXP(Y(J)-V/TL-NEXP*1.)
		 DO IQ1=1,NQ1
		    Q1=Q1MIN+(IQ1-0.5)*DQ1
		    DO IQ2=1,NQ2
		       Q2=Q2MIN+(IQ2-0.5)*DQ2
		       Q1AV=Q1AV+Q1*HIST_T(J,IQ1,IQ2)
		       Q2AV=Q2AV+Q2*HIST_T(J,IQ1,IQ2)
		       Q1AV2=Q1AV2+Q1*Q1*HIST_T(J,IQ1,IQ2)
		       Q2AV2=Q2AV2+Q2*Q2*HIST_T(J,IQ1,IQ2)
		       AVNUM=AVNUM+HIST_T(J,IQ1,IQ2)
		    ENDDO
		 ENDDO
	      ELSE
		 XP=0.0
	      ENDIF
	      Z=Z+XP
	      U=U+V*XP
	      S=S+V*V*XP
	   ENDDO
	   U=U/Z
	   S=S/Z
	   CV=(S-U*U)/(TL*TL)
	   Q1AV=Q1AV/AVNUM
	   Q2AV=Q2AV/AVNUM
	   Q1AV2=Q1AV2/AVNUM
	   Q2AV2=Q2AV2/AVNUM
	   FQ1=SQRT(Q1AV2-Q1AV**2)
	   FQ2=SQRT(Q2AV2-Q2AV**2)

	   WRITE(11,*) TL,U
	   WRITE(12,*) TL,CV
	   WRITE(13,*) TL,Q1AV,FQ1
	   WRITE(14,*) TL,Q2AV,FQ2

	PRINT *,'Cv,Q1,Q2=',cv,Q1av,Q2av
	   
	ENDDO

	CLOSE(11)
	CLOSE(12)
	CLOSE(13)
	CLOSE(14)

	PRINT *,'Files U, Cv, Q1, Q2 calculated'

	STOP 'Normal end.'

	END

