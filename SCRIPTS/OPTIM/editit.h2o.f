C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      PROGRAM EDIT
      LOGICAL TEST

      PRINT*,'punch coordinates energy gradient secder transform'
      PRINT*,'title'
      PRINT*,'(H2O) minimum optimization - RHF/aug-ccpVTZ'
      
      PRINT*,'geometry'
      READ(*,*) X,Y,Z
      WRITE(*,'(3F20.10,A)') X,Y,Z,' 8.0 o'
      READ(*,*) X,Y,Z
      WRITE(*,'(3F20.10,A)') X,Y,Z,' 1.0 h1'
      READ(*,*) X,Y,Z
      WRITE(*,'(3F20.10,A)') X,Y,Z,' 1.0 h2'
      PRINT*,'end'
      PRINT*,'basis'
      PRINT*,'S   H1'
      PRINT*,'0.01968500        13.01000000'
      PRINT*,'0.13797700         1.96200000'
      PRINT*,'0.47814800         0.44460000'
      PRINT*,'S   H1 '
      PRINT*,'1.00000000         0.12200000'
      PRINT*,'P   H1 '
      PRINT*,'1.00000000         0.72700000'
      PRINT*,'S   H1 '
      PRINT*,'1.00000000         0.02974000'
      PRINT*,'P   H1 '
      PRINT*,'1.00000000         0.14100000'
      PRINT*,'S   H2'
      PRINT*,'0.01968500        13.01000000'
      PRINT*,'0.13797700         1.96200000'
      PRINT*,'0.47814800         0.44460000'
      PRINT*,'S   H2 '
      PRINT*,'1.00000000         0.12200000'
      PRINT*,'P   H2 '
      PRINT*,'1.00000000         0.72700000'
      PRINT*,'S   H2 '
      PRINT*,'1.00000000         0.02974000'
      PRINT*,'P   H2 '
      PRINT*,'1.00000000         0.14100000'
      PRINT*,'S   O '
      PRINT*,'0.00071000     11720.00000000'
      PRINT*,'0.00547000      1759.00000000'
      PRINT*,'0.02783700       400.80000000'
      PRINT*,'0.10480000       113.70000000'
      PRINT*,'0.28306200        37.03000000'
      PRINT*,'0.44871900        13.27000000'
      PRINT*,'0.27095200         5.02500000'
      PRINT*,'0.01545800         1.01300000'
      PRINT*,'S   O '
      PRINT*,'-0.00016000     11720.00000000'
      PRINT*,'-0.00126300      1759.00000000'
      PRINT*,'-0.00626700       400.80000000'
      PRINT*,'-0.02571600       113.70000000'
      PRINT*,'-0.07092400        37.03000000'
      PRINT*,'-0.16541100        13.27000000'
      PRINT*,'-0.11695500         5.02500000'
      PRINT*,'0.55736800         1.01300000'
      PRINT*,'S   O '
      PRINT*,'1.00000000         0.30230000'
      PRINT*,'P   O '
      PRINT*,'0.04301800        17.70000000'
      PRINT*,'0.22891300         3.85400000'
      PRINT*,'0.50872800         1.04600000'
      PRINT*,'P   O '
      PRINT*,'1.00000000         0.27530000'
      PRINT*,'D   O '
      PRINT*,'1.00000000         1.18500000'
      PRINT*,'S   O '
      PRINT*,'1.00000000         0.07896000'
      PRINT*,'P   O '
      PRINT*,'1.00000000         0.06856000'
      PRINT*,'D   O '
      PRINT*,'1.00000000         0.33200000'
      PRINT*,'end'
      PRINT*,'scftype mp2'

      STOP
      END
