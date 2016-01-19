!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales 
!   This file is part of GMIN.                                      
!
!   GMIN is free software; you can redistribute it and/or modIFy
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version. 
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; IF not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-130 US
!
!=============================================================
!   All routines in this file were implemented by
!   Brooke Husic (beh35) in 2015. 
!=============================================================
!
!   What the mlowest routines do and where they are:
!
!   In commons.f90, three new parameters are defined:
!
!     1. Logical MLOWEST: If MLOWEST is true, GMIN diverts to special
!        behavior where it will not stop until it reaches mutliple
!        specified target energies.
!        
!        Logicals MLOWEST and TARGET are mutually exclusive.
!        Specifying MLOWEST instead of TARGET results in similar
!        behavior with a few changes to account for needing to hit
!        multiple energy targets before quitting GMIN.
!
!     2. Integer MTARGETS intidcates how many energies are specified
!
!     3. Type MTARGETS_ARRAY: This replaces the type TARGETS which is
!        used when TARGET is specified. Each MTARGETS_ARRAY type is
!        an ordered pair of a target energy (MTARGETS_ARRAY%SOUGHT) and
!        a logicial initialized as FALSE (MTARGETS_ARRAY%FOUND), the
!        latter of which identifies whether the corresponding energy
!        has been hit yet or not.
!
!   In keywords.f, more framework:
!
!     1. MLOWEST logical defaults to FALSE
!          If MLOWEST keyword is specified, MLOWEST logical changes to
!          TRUE. This is mutually exclusive with TARGET=TRUE
!          
!          A type MTARGETS_ARRAY called OBJ is allocated dimensions
!          equal to MTARGETS and the energies are inserted into
!          the array using create_mtargets_array in mlowest.f90 
!
!   In io1.F, if MLOWEST is TRUE, target energies are assigned to
!   the SOUGHT position of the  MTARGETS_ARRAY called OBJ.
!
!   In saveit.F, if MLOWEST is TRUE, when a target energy is hit,  
!   the OBJ%FOUND(i) logical corresponding to the target energy value
!   OBJ%SOUGHT(i) is changed to TRUE. Then, decide_to_quit from
!   mlowest.f90 is called to check and see if ALL(OBJ%SOUGHT) is TRUE.
!   If so, HIT is set to TRUE which ends GMIN searching just like it
!   does when TARGET is used.
!
!
!
!!

SUBROUTINE create_mtargets_array(targ) 

! This subroutine takes every target and makes it an ordered pair
! containing the target energy value and the logical "false"

USE COMMONS, ONLY : OBJ, MTARGETS 

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: targ(MTARGETS)
    INTEGER                       :: i

    ! Make MTARGETS_ARRAY type with the same number of dimensions
    ! as targets listed where the first element of each entry
    ! is the target energy and the second is a logical "false"
    ! meaning the target hasn't been hit yet

DO i=1,MTARGETS

    obj(i)%sought = targ(i)
    obj(i)%found  = .FALSE.

ENDDO

END SUBROUTINE create_mtargets_array
    

SUBROUTINE decide_to_quit()

! This subroutine inputs the working array of mtargets and
! if all targets have been found it sets HIT to true 

USE COMMONS, ONLY: OBJ, HIT 

    IMPLICIT NONE

!    IF(MLOWEST) THEN

        IF  (ALL(OBJ%found)) THEN
            HIT = .TRUE.
        END IF
    
!    ELSE
!        IF (ANY(ntargets%found)) THEN
!            ! space for something else to happen here
!        END IF
!
!    END IF

END SUBROUTINE decide_to_quit



