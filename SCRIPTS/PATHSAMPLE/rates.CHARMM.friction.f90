!
! Calculate NGT rates from pathsample database for different values of 
! regroupfree and friction parameters.
!
! For CHARMM and AMBER diagonalise the reciprocal mass-weighted Hessian in OPTIM, where the various masses are known. 
! For convenience the frequency unit conversion is done in OPTIM as well, so that the rate constants calculated in 
! PATHSAMPLE are in s and do not need to be converted. The value required for the Planck constant is therefore different 
! because the frquency is not in reduced units. 
! Instead, we need to convert  to kcal/mol, since these are the units of kT. 
! Hence we need (h/Js)/(epsilon/J), where epsilon is one kcal/mol. Since 1kcal/mol is 6.948D-21 J the required value 
! for the PLANCK keyword in regrouping calculations is 9.536D-14..
!
PROGRAM RATES

DOUBLE PRECISION EPSILON, PLANCK, RT, RGF, RATE, FRIC

WRITE(*,'(A)') &
&'  epsilon/(kcal/mol) regroupfree/epsilon       PLANCK*   room temperature/epsilon     rate factor         rate/Hz    ' &
& // '        friction/Hz'

EPSILON=1.0D0 ! energy unit is one kcal/mol
DO J2=1,20 ! regroupfree values
   RGF=(1.0D0*(J2-1)) ! REGROUPFREE parameter is in EPSILON units, i.e. kcal/mol
!
! Since the frequencies are obtained in Hz in pathsample, the friction coefficient should
! be in Hz.
!
   DO J3=1,20
      FRIC=(J3-1)*0.1D0*1.0D12 ! convert to per picosecond
      CALL SYSTEM('cp pathdata.template pathdata')
      OPEN(UNIT=1,FILE='temp',STATUS='UNKNOWN')
      PLANCK=9.536D-14
      RT=0.592D0 ! this is kT in kcal/mol at T=298 K
      WRITE(1,'(A,G20.10)') 'PLANCK ',PLANCK
      WRITE(1,'(A,G20.10)') 'REGROUPFREE ',RGF
      WRITE(1,'(A,G20.10)') 'TEMPERATURE ',RT
      WRITE(1,'(A,G20.10)') 'FRICTION ',FRIC
      WRITE(1,'(A,I6)') 'CYCLES ',0
      WRITE(1,'(A,I6)') 'NGT ',0
      CLOSE(1)
      CALL SYSTEM('cat temp >> pathdata')
      CALL SYSTEM('~wales/svn/PATHSAMPLE/bin/pathsample.2.1 > output.NGT')
      CALL SYSTEM('grep kNSS output.NGT | tail -1 | sed -e "s/.*B)=//" > temp2')
      OPEN(UNIT=1,FILE='temp2',STATUS='OLD')
      READ(1,*) RATE
      CLOSE(1)
      WRITE(*,'(6G20.10)') EPSILON, RGF, PLANCK, RT, RATE, FRIC
      CALL FLUSH(6)
   ENDDO
      WRITE(*,'(A)') ' '
ENDDO

END
