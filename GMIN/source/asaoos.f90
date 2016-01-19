!|gd351>

SUBROUTINE ASAOOSPOT (X, G, ENERGY, GTEST)

  USE COMMONS, ONLY: NATOMS, SIGMAP

  IMPLICIT NONE
  
  INTEGER          :: J, J1, J2, J3
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
  DOUBLE PRECISION :: ENERGY, R, RSQ, ARG, V, DV
  DOUBLE PRECISION :: RI(NATOMS,3), RIJ(3)
  LOGICAL          :: GTEST
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
  DOUBLE PRECISION :: AT1, UEL, UET, UECL, UECL1, UECL2, UECH, UESM1, UESM2, UESM
  DOUBLE PRECISION :: UELSQ, UETSQ, UECL1SQ, UECL2SQ, UECLSQ, UECHSQ, UESM1SQ, UESM2SQ, UESMSQ
  DOUBLE PRECISION :: VLV0, VLK0, V_1, PA1, PA2, PA3, PA4
  DOUBLE PRECISION :: VS1A0, VS1A1, VS1A2, VS1A3, VS1A4, VS2A0, VS2A1, VS2A2, VS2A3, VS2A4
  DOUBLE PRECISION :: SIGMAC, Q, ZP, BETA, SIGMAPSQ, SIGMACSQ, RANGESQ, ONEPLUSQCB, PREFACTOR, VMIN, PFATH, PFATHD
  

  !hard core approximation:
  !             r < 0.997    : ~ linear
  ! 0.997     < r < 0.998    : ~ TANH[ AT1*(r-0.9975) ]
  ! 0.998     < r < 1.000    : hard core and attractive potential smoothly merged by polynomial of grade 3
  !attractive potential
  ! 1.00000   < r < 1+SIGMAP : Asakura-Oosawa potential
  ! 1.+SIGMAP < r            : 0 

  !parameters change when ZP, BETA, SIGMAC or SIGMAP are changed

  SIGMAC=1.D0
!  SIGMAP=2.D-1

  AT1  = 100.D0 !!higher AT1 -> steeper hardcore approx
  UEL  = 0.985D0
  UET  = 0.99D0
  UESM = 1.02D0

  UELSQ  = (SIGMAC*UEL)**2
  UETSQ  = (SIGMAC*UET)**2
  UESMSQ = (SIGMAC*UESM)**2

  Q=SIGMAP/SIGMAC

  ZP=1.D2
  BETA=1.D0

  SIGMACSQ=SIGMAC**2
  SIGMAPSQ=SIGMAP**2
  RANGESQ=(SIGMAC+SIGMAP)**2

  ONEPLUSQCB=(1.D0+Q)**3
  PREFACTOR=-PI*ZP*ONEPLUSQCB*(SIGMAP**3)/(6.D0*BETA*Q**3)

  VMIN = PREFACTOR * (1.D0 - (3.D0)/(2.D0*SIGMAC*(1.D0+Q)) + (1.D0)/(2.D0*ONEPLUSQCB*SIGMAC**3))
  PREFACTOR = PREFACTOR / ABS(VMIN)

  V_1=PREFACTOR*(1.D0 - (3.D0)/(2.D0*SIGMAC*(1.D0+Q)) + (1.D0)/(2.D0*ONEPLUSQCB*SIGMAC**3))

  !parameters have to be recalculated when changing AT1
  !parameters for smoothing function VS (polynomial of grade 3)
  !VS(UET)  = VTANH(UET)
  !VS'(UET) = VTANH'(UET)
  !VS(1.00)   = VAO(1.00)
  !VS'(1.00)  = VAO'(1.00)

  VS1A0 = 5.06351526092385D6
  VS1A1 = -2.0396594216841813D7
  VS1A2 = 3.081184200916632D7
  VS1A3 = -2.0687966411502566D7
  VS1A4 = 5.209202358254224D6

  PA1  = 1.9764998955815453D7
  PA2  = 1.976499995581546D7
  PA3  = 1.0081475977469167D7
  PA4  = 2.9846474933284618D7
  VS2A0 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = -7.802999982557824D7
  PA2  = -7.802999982557827D7
  PA3  = -3.979920141105787D7
  PA4  = -1.178291997366361D8
  VS2A1 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = 1.1551499974180764D8
  PA2  = 1.1551499974180768D8
  PA3  = 5.8916474868341096D7
  PA4  = 1.7443147461014873D8
  VS2A2 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = -7.599999983014236D7
  PA2  = -7.59999998301424D7
  PA3  = -3.876124941338524D7
  PA4  = -1.147612497435276D8
  VS2A3 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = 1.8749999958097514D7
  PA2  = 1.874999995809752D7
  PA3  = 9.562499978632849D6
  PA4  = 2.831249993673036D7
  VS2A4 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))


!!$  VS1A0 = 5.06351526092385D6
!!$  VS1A1 = -2.0396594216841813D7
!!$  VS1A2 = 3.081184200916632D7
!!$  VS1A3 = -2.0687966411502566D7
!!$  VS1A4 = 5.209202358254224D6
!!$
!!$  PA1  = 3.080599642740617D8
!!$  PA2  = 3.0805996527406156D8
!!$  PA3  = 1.555754334632798D8
!!$  PA4  = 4.636353977373414D8
!!$  VS2A0 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = -1.2241198620131054D9
!!$  PA2  = -1.2241198620131047D9
!!$  PA3  = -6.18195831816136D8
!!$  PA4  = -1.8423156923292408D9
!!$  VS2A1 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = 1.8240597943881383D9
!!$  PA2  = 1.8240597943881376D9
!!$  PA3  = 9.211653461652914D8
!!$  PA4  = 2.745225140553429D9
!!$  VS2A2 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = -1.2079998638332074D9
!!$  PA2  = -1.207999863833207D9
!!$  PA3  = -6.10044930735294D8
!!$  PA4  = -1.8180447950685012D9
!!$  VS2A3 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = 2.999999661841128D8
!!$  PA2  = 2.9999996618411267D8
!!$  PA3  = 1.514999829228588D8
!!$  PA4  = 4.514999491069715D8
!!$  VS2A4 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))



  !parameter for core linear
  !VL(UEL)  = VTANH(UEL)
  !VL'(UEL) = VTANH'(UEL)

  VLV0 =  48.7D0
  VLK0 = -50.0D0


  ENERGY = 0.D0
  IF (GTEST) G(:) = 0.D0
  
  DO J1 = 1, NATOMS
    
    J3 = 3*J1
    RI(J1,1:3) = X(J3-2:J3)
       
  END DO


  DO J1 = 1, NATOMS

    DO J2 = J1+1, NATOMS

      RIJ(:) = RI(J1,:) - RI(J2,:)
      RSQ = DOT_PRODUCT(RIJ(:),RIJ(:))


      IF (RSQ.LE.UELSQ) THEN
      !core - linear
        R = SQRT(RSQ)
        
        ENERGY = ENERGY + (VLV0 + VLK0 * R)  

        IF (GTEST) THEN
          !calculate derivatives
          DV = VLK0 / R   !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 


      ELSEIF ((RSQ.GT.UELSQ).AND.(RSQ.LE.UETSQ)) THEN
      !core - ArcTanh
        R = SQRT(RSQ)
        ARG = AT1*(R-UEL)
        V =  0.5D0*(1.D0-((EXP(ARG)-EXP(-ARG))/(EXP(ARG)+EXP(-ARG)))) - 1.05D0 ! = 0.5*(1-TANH(ARG))           
 
        ENERGY = ENERGY + V

        IF (GTEST) THEN
          !calculate derivatives
          DV = -2.D0*AT1 / ( (EXP(ARG)+EXP(-ARG))**2 * R )  !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 


      ELSEIF ((RSQ.GT.UETSQ).AND.(RSQ.LE.SIGMACSQ)) THEN
      !smooth merging
        R = SQRT(RSQ)
        V = VS1A0 + VS1A1 * R + VS1A2 * RSQ + VS1A3 * RSQ*R + VS1A4 * RSQ*RSQ

        ENERGY = ENERGY + V

        IF (GTEST) THEN
          !calculate derivatives
          DV = VS1A1 / R + 2.D0 * VS1A2 + 3.D0 * VS1A3 * R + 4.D0 * VS1A4 * RSQ!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 

      ELSEIF ((RSQ.GT.SIGMACSQ).AND.(RSQ.LE.UESMSQ)) THEN
      !smooth merging
        R = SQRT(RSQ)
        V = VS2A0 + VS2A1 * R + VS2A2 * RSQ + VS2A3 * RSQ*R + VS2A4 * RSQ*RSQ

        ENERGY = ENERGY + V

        IF (GTEST) THEN
          !calculate derivatives
          DV = VS2A1 / R + 2.D0 * VS2A2 + 3.D0 * VS2A3 * R + 4.D0 * VS2A4 * RSQ!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 


      ELSEIF ((RSQ.GT.UESMSQ).AND.(RSQ.LE.RANGESQ)) THEN
      !Asakura-Oosawa attraction
        R = SQRT(RSQ)
        V = PREFACTOR * (1.D0 - (3.D0*R)/(2.D0*SIGMAC*(1.D0+Q)) + (RSQ*R)/(2.D0*ONEPLUSQCB*SIGMAC**3))

        ENERGY = ENERGY + V   

        IF (GTEST) THEN
          !calculate derivatives
          DV = PREFACTOR * ( - 3.D0/(R*2.D0*SIGMAC*(1.D0+Q)) + 3.D0*R/(2.D0*ONEPLUSQCB*SIGMAC**3))!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
        END IF 

      END IF

    END DO
  END DO

END SUBROUTINE ASAOOSPOT


SUBROUTINE ASAOOSPRINT()

  USE COMMONS, ONLY: NATOMS, SIGMAP

  IMPLICIT NONE
  
  INTEGER          :: J, J1, J2, J3
  DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS)
  DOUBLE PRECISION :: ENERGY, R, RSQ, ARG, V, DV
  DOUBLE PRECISION :: RI(NATOMS,3), RIJ(3)
  LOGICAL          :: GTEST
  DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
  DOUBLE PRECISION :: AT1, UEL, UET, UECL, UECL1, UECL2, UECH, UESM1, UESM2, UESM
  DOUBLE PRECISION :: UELSQ, UETSQ, UECL1SQ, UECL2SQ, UECLSQ, UECHSQ, UESM1SQ, UESM2SQ, UESMSQ
  DOUBLE PRECISION :: VLV0, VLK0, V_1, PA1, PA2, PA3, PA4
  DOUBLE PRECISION :: VS1A0, VS1A1, VS1A2, VS1A3, VS1A4, VS2A0, VS2A1, VS2A2, VS2A3, VS2A4
  DOUBLE PRECISION :: SIGMAC, Q, ZP, BETA, SIGMAPSQ, SIGMACSQ, RANGESQ, ONEPLUSQCB, PREFACTOR, VMIN, PFATH, PFATHD
  
  
  GTEST=.TRUE.  
  OPEN(84,file='asaoospot.dat')
  OPEN(85,file='asaoosder.dat')
  
  !hard core approximation:
  !             r < 0.997    : ~ linear
  ! 0.997     < r < 0.998    : ~ TANH[ AT1*(r-0.9975) ]
  ! 0.998     < r < 1.000    : hard core and attractive potential smoothly merged by polynomial of grade 3
  !attractive potential
  ! 1.00000   < r < 1+SIGMAP : Asakura-Oosawa potential
  ! 1.+SIGMAP < r            : 0 

  !parameters change when ZP, BETA, SIGMAC or SIGMAP are changed

  SIGMAC=1.D0
!  SIGMAP=2.D-1

  AT1  = 100.D0 !!higher AT1 -> steeper hardcore approx
  UEL  = 0.985D0
  UET  = 0.99D0
  UESM = 1.02D0

  UELSQ  = (SIGMAC*UEL)**2
  UETSQ  = (SIGMAC*UET)**2
  UESMSQ = (SIGMAC*UESM)**2

  Q=SIGMAP/SIGMAC

  ZP=1.D2
  BETA=1.D0

  SIGMACSQ=SIGMAC**2
  SIGMAPSQ=SIGMAP**2
  RANGESQ=(SIGMAC+SIGMAP)**2

  ONEPLUSQCB=(1.D0+Q)**3
  PREFACTOR=-PI*ZP*ONEPLUSQCB*(SIGMAP**3)/(6.D0*BETA*Q**3)

  VMIN = PREFACTOR * (1.D0 - (3.D0)/(2.D0*SIGMAC*(1.D0+Q)) + (1.D0)/(2.D0*ONEPLUSQCB*SIGMAC**3))
  PREFACTOR = PREFACTOR / ABS(VMIN)

  V_1=PREFACTOR*(1.D0 - (3.D0)/(2.D0*SIGMAC*(1.D0+Q)) + (1.D0)/(2.D0*ONEPLUSQCB*SIGMAC**3))

  !parameters have to be recalculated when changing AT1
  !parameters for smoothing function VS (polynomial of grade 3)
  !VS(UET)  = VTANH(UET)
  !VS'(UET) = VTANH'(UET)
  !VS(1.00)   = VAO(1.00)
  !VS'(1.00)  = VAO'(1.00)

  VS1A0 = 5.06351526092385D6
  VS1A1 = -2.0396594216841813D7
  VS1A2 = 3.081184200916632D7
  VS1A3 = -2.0687966411502566D7
  VS1A4 = 5.209202358254224D6

  PA1  = 1.9764998955815453D7
  PA2  = 1.976499995581546D7
  PA3  = 1.0081475977469167D7
  PA4  = 2.9846474933284618D7
  VS2A0 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = -7.802999982557824D7
  PA2  = -7.802999982557827D7
  PA3  = -3.979920141105787D7
  PA4  = -1.178291997366361D8
  VS2A1 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = 1.1551499974180764D8
  PA2  = 1.1551499974180768D8
  PA3  = 5.8916474868341096D7
  PA4  = 1.7443147461014873D8
  VS2A2 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = -7.599999983014236D7
  PA2  = -7.59999998301424D7
  PA3  = -3.876124941338524D7
  PA4  = -1.147612497435276D8
  VS2A3 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))

  PA1  = 1.8749999958097514D7
  PA2  = 1.874999995809752D7
  PA3  = 9.562499978632849D6
  PA4  = 2.831249993673036D7
  VS2A4 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))


!!$  VS1A0 = 5.06351526092385D6
!!$  VS1A1 = -2.0396594216841813D7
!!$  VS1A2 = 3.081184200916632D7
!!$  VS1A3 = -2.0687966411502566D7
!!$  VS1A4 = 5.209202358254224D6
!!$
!!$  PA1  = 3.080599642740617D8
!!$  PA2  = 3.0805996527406156D8
!!$  PA3  = 1.555754334632798D8
!!$  PA4  = 4.636353977373414D8
!!$  VS2A0 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = -1.2241198620131054D9
!!$  PA2  = -1.2241198620131047D9
!!$  PA3  = -6.18195831816136D8
!!$  PA4  = -1.8423156923292408D9
!!$  VS2A1 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = 1.8240597943881383D9
!!$  PA2  = 1.8240597943881376D9
!!$  PA3  = 9.211653461652914D8
!!$  PA4  = 2.745225140553429D9
!!$  VS2A2 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = -1.2079998638332074D9
!!$  PA2  = -1.207999863833207D9
!!$  PA3  = -6.10044930735294D8
!!$  PA4  = -1.8180447950685012D9
!!$  VS2A3 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))
!!$
!!$  PA1  = 2.999999661841128D8
!!$  PA2  = 2.9999996618411267D8
!!$  PA3  = 1.514999829228588D8
!!$  PA4  = 4.514999491069715D8
!!$  VS2A4 = PA1 - PA2/(((1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))) - &
!!$  PA3/(((1.D0+SIGMAP)**3)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP))) + &
!!$  PA4/((1.D0+SIGMAP)*(1.D0+0.5D0/(1.D0+SIGMAP)**3-1.5D0/(1.D0+SIGMAP)))



  !parameter for core linear
  !VL(UEL)  = VTANH(UEL)
  !VL'(UEL) = VTANH'(UEL)

  VLV0 =  48.7D0
  VLK0 = -50.0D0


  ENERGY = 0.D0
  IF (GTEST) G(:) = 0.D0

  
    DO J1 = 0, 2000

      RSQ=(0.985D0+(J1/10000.D0))**2


      IF (RSQ.LE.UELSQ) THEN
      !core - linear
        R = SQRT(RSQ)
        
        ENERGY = ENERGY + (VLV0 + VLK0 * R)  
        WRITE(84,*) 'lin',R, VLV0 + VLK0 * R

        IF (GTEST) THEN
          !calculate derivatives
          DV = VLK0 / R   !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
          WRITE(85,*) 'lin',R, DV*R
        END IF 


      ELSEIF ((RSQ.GT.UELSQ).AND.(RSQ.LE.UETSQ)) THEN
      !core - ArcTanh
        R = SQRT(RSQ)
        ARG = AT1*(R-UEL)
        V =  0.5D0*(1.D0-((EXP(ARG)-EXP(-ARG))/(EXP(ARG)+EXP(-ARG)))) -1.05D0 ! = 0.5*(1-TANH(ARG))           
 
        ENERGY = ENERGY + V
        WRITE(84,*) 'ath',R, V

        IF (GTEST) THEN
          !calculate derivatives
          DV = -2.D0*AT1 / ( (EXP(ARG)+EXP(-ARG))**2 * R )  !factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
          WRITE(85,*) 'ath',R, DV*R
        END IF 


      ELSEIF ((RSQ.GT.UETSQ).AND.(RSQ.LE.SIGMACSQ)) THEN
      !smooth merging
        R = SQRT(RSQ)
        V = VS1A0 + VS1A1 * R + VS1A2 * RSQ + VS1A3 * RSQ*R + VS1A4 * RSQ*RSQ

        ENERGY = ENERGY + V
        WRITE(84,*) 'vs1',R, V

        IF (GTEST) THEN
          !calculate derivatives
          DV = VS1A1 / R + 2.D0 * VS1A2 + 3.D0 * VS1A3 * R + 4.D0 * VS1A4 * RSQ!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
          WRITE(85,*) 'vs1',R, DV*R
        END IF 

      ELSEIF ((RSQ.GT.SIGMACSQ).AND.(RSQ.LE.UESMSQ)) THEN
      !smooth merging
        R = SQRT(RSQ)
        V = VS2A0 + VS2A1 * R + VS2A2 * RSQ + VS2A3 * RSQ*R + VS2A4 * RSQ*RSQ

        ENERGY = ENERGY + V
        WRITE(84,*) 'vs2',R, V

        IF (GTEST) THEN
          !calculate derivatives
          DV = VS2A1 / R + 2.D0 * VS2A2 + 3.D0 * VS2A3 * R + 4.D0 * VS2A4 * RSQ!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
          WRITE(85,*) 'vs2',R, DV*R
        END IF 


      ELSEIF ((RSQ.GT.UESMSQ).AND.(RSQ.LE.RANGESQ)) THEN
      !Asakura-Oosawa attraction
        R = SQRT(RSQ)
        V = PREFACTOR * (1.D0 - (3.D0*R)/(2.D0*SIGMAC*(1.D0+Q)) + (RSQ*R)/(2.D0*ONEPLUSQCB*SIGMAC**3))

        ENERGY = ENERGY + V 
        WRITE(84,*) 'aso',R, V

        IF (GTEST) THEN
          !calculate derivatives
          DV = PREFACTOR * ( - 3.D0/(R*2.D0*SIGMAC*(1.D0+Q)) + 3.D0*R/(2.D0*ONEPLUSQCB*SIGMAC**3))!factor 1/R from dR/dXi = Xi/R already included
          G(J1*3-2:J1*3) = G(J1*3-2:J1*3) + DV * RIJ(:)
          G(J2*3-2:J2*3) = G(J2*3-2:J2*3) - DV * RIJ(:)
          WRITE(85,*) 'aso',R, DV*R
        END IF 

      END IF


  END DO

CLOSE(84)
CLOSE(85)

END SUBROUTINE ASAOOSPRINT



!<gd351|
