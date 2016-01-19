!
! Calculate NGT rates from pathsample database for different values of AMH
! epsilon parameter and regroupfree values.
!
PROGRAM RATES

DOUBLE PRECISION EPSILON, PLANCK, RT, RATEFAC, RGF, RATE

WRITE(*,'(A)') '     epsilon         regroupfree/epsilon       PLANCK    room temperature/epsilon     rate factor         rate/Hz'

DO J1=1,20    ! epsilon
   EPSILON=(J1*1.0D0)*0.1D0 ! kcal/mol
!
!  prefactor is sqrt(1000*4.184/6.022D23*12*1.661D-27*1.0D-20)
!
   RATEFAC=5.8785D12*SQRT(EPSILON)
   DO J2=1,20 ! regroupfree values
      RGF=(1.0D0*(J2-1)) ! REGROUPFREE parameter is in EPSILON units
      CALL SYSTEM('cp pathdata.template pathdata')
      OPEN(UNIT=1,FILE='temp',STATUS='UNKNOWN')
!
!  prefactor is 6.626D-34/sqrt(12*1.661D-27*1.0D-20*1000*4.184/6.022D23)
!
      PLANCK=0.563D0/SQRT(EPSILON)
      RT=0.59D0/EPSILON
      WRITE(1,'(A,G20.10)') 'PLANCK ',PLANCK
      WRITE(1,'(A,G20.10)') 'REGROUPFREE ',RGF
      WRITE(1,'(A,G20.10)') 'TEMPERATURE ',RT
      CLOSE(1)
      CALL SYSTEM('cat temp >> pathdata')
      CALL SYSTEM('~wales/svn/PATHSAMPLE/bin/pathsample.2.1 > output.NGT')
      CALL SYSTEM('grep kNSS output.NGT | tail -1 | sed -e "s/.*B)=//" > temp2')
      OPEN(UNIT=1,FILE='temp2',STATUS='OLD')
      READ(1,*) RATE
      CLOSE(1)
      WRITE(*,'(6G20.10)') EPSILON, RGF, PLANCK, RT, RATEFAC, RATE*RATEFAC
      CALL FLUSH(6)
   ENDDO
ENDDO

END
