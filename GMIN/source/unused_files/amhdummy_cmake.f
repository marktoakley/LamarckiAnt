      SUBROUTINE WALESAMH_INITIAL()

      RETURN
      END

      SUBROUTINE NUM_TO_CHAR(COUNT,CCOUNT)

      INTEGER COUNT
      CHARACTER CCOUNT*1

      RETURN
      END

!       SUBROUTINE E_WRITE(AVEP,T,NUMPRO,I_TEMP)

!       DOUBLE PRECISION T,AVEP(:,:,:)
!       INTEGER NUMPRO,I_TEMP

!       RETURN
!       END

      SUBROUTINE WALESAMH_INTERFACE(COORD_MCP,GRAD_FOR_WALES,E_FOR_WALES)
      IMPLICIT NONE

      DOUBLE PRECISION GRAD_FOR_WALES(*),E_FOR_WALES
      DOUBLE PRECISION COORD_MCP(*)

      RETURN
      END

!      MODULE AMHGLOBALS

!      INTEGER NMRES
!      INTEGER OMOVI
!      INTEGER NUMPRO
!      INTEGER IRES(500)
!      DOUBLE PRECISION AVEP(10,2,50)

 !     end MODULE AMHGLOBALS

!      MODULE AMH_INTERFACES
!      contains
!	SUBROUTINE E_WRITE(DUMMY1,DUMMY2,DUMMY3,DUMMY4)
!	END SUBROUTINE E_WRITE
!      END MODULE AMH_INTERFACES
