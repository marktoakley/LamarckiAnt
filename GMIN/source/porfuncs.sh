#!/bin/bash

# GPL License Info  {{{
#
#   Portability functions module generator
#   Copyright (C) 2003-2014 Semen A. Trygubenko and David J. Wales
#   This file is part of GMIN.
#
#   GMIN is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   GMIN is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#}}}

if test "x$1" == "x"; then :
  FC="pgi"
else
  FC=$1
  if test "x${FC}" != "xnag" -a "x${FC}" != "xifort" -a "x${FC}" != "xpgi" -a "x${FC}" != "xifc" -a "x${FC}" != "xg95" -a "x${FC}" != "xgfortran" -a "x${FC}" != "xpathscale"; then :
    echo "Unknown compiler!"
    exit 1
  fi
fi

echo "MODULE PORFUNCS"
if test "x${FC}" == "xnag"; then :
  echo "     use f90_unix, only: getarg"
  echo "     use f90_unix_proc, only: system, exit"
fi

echo "     implicit none"
echo "     contains"

if test "x${FC}" == "xnag"; then :
  echo "          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be"
  echo "          use f90_unix, NAGflush => flush ! connected for formatted sequential output; ISTAT is ignored"
  echo "               implicit none"
  echo "               integer,intent(in) :: UNIT"
  echo "               integer,intent(out),optional :: ISTAT"
  echo "               call NAGflush(UNIT)"
  echo "          end subroutine flush"
  echo " "
elif test "x${FC}" == "xg95"; then :
  echo "          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be"
  echo "               implicit none ! connected for formatted sequential output; ISTAT is ignored"
  echo "               integer,intent(in) :: UNIT"
  echo "               integer,intent(out),optional :: ISTAT"
  echo "               call flush(UNIT)"
  echo "          end subroutine flush"
  echo " "
elif test "x${FC}" == "xgfortran"; then :
  echo "          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be"
  echo "               implicit none ! connected for formatted sequential output; ISTAT is ignored"
  echo "               integer,intent(in) :: UNIT"
  echo "               integer,intent(out),optional :: ISTAT"
  echo "!              call flush(UNIT)"
  echo "          end subroutine flush"
  echo " "
fi

echo "          subroutine getarg_subr(position,value) ! wraps getarg function so it can be use-associated"

if test "x${FC}" == "xnag"; then :
  echo "               use f90_unix, only: getarg"
fi

echo "               implicit none"
echo " "
echo "               integer,intent(in) :: position"
echo "               character(len=*),intent(out) :: value"
echo " "
echo " "
echo "               call getarg(position,value)"
echo "          end subroutine getarg_subr"
echo " "
echo "          subroutine iargc_subr(n) ! wraps iargc function so it can be use-associated"

if test "x${FC}" == "xnag"; then :
  echo "               use f90_unix, only: iargc"
fi

echo "               implicit none"
echo "               integer,intent(out) :: n"
echo " "

if test "x${FC}" != "xnag"; then :
  echo "               integer iargc"
  echo " "
fi

echo "               n = iargc()"
echo "          end subroutine iargc_subr"
echo " "
echo "          subroutine fork_subr(pid) ! returns zero in the child process, PID of child in parent process"

if test "x${FC}" == "xnag"; then :
  echo "               use f90_unix_proc, only: fork"
fi

echo "               implicit none"
echo "               integer, intent(inout) :: pid"
echo " "

if test "x${FC}" == "xpgi"; then :
     echo "               integer fork"
fi

echo " "

if test "x${FC}" == "xpgi"; then :
     echo "               pid=fork()"
elif test "x${FC}" == "xifort"; then :
     echo "               integer ierror"
     echo "               call pxffork(pid,ierror)"
elif test "x${FC}" == "xnag"; then :
     echo "               call fork(pid)"
elif test "x${FC}" == "xpathscale"; then :
     echo "               integer fork"
     echo "               pid=fork()"
fi

echo "          end subroutine fork_subr"
echo " "
echo "          subroutine system_subr(JobString,ExitStatus)"

if test "x${FC}" == "xnag"; then :
  echo "               use f90_unix_proc, only: system"
fi

echo "               implicit none"
echo " "
echo "               character(len=*),intent(in) :: JobString"
echo "               integer,intent(out) :: ExitStatus"
echo " "

if test "x${FC}" == "xifort"; then :
  echo "               integer shiftr,system"
  echo " "
  echo "               ExitStatus=system(JobString)"
elif test "x${FC}" == "xpgi"; then :
  echo "               integer system"
  echo " "
  echo "               ExitStatus=system(JobString)"
  echo "               ExitStatus=ishft(ExitStatus,-8)"
elif test "x${FC}" == "xnag"; then :
  echo "               call system(JobString,ExitStatus)"
  echo "               ExitStatus=ishft(ExitStatus,-8)"
elif test "x${FC}" == "xpathscale"; then :
  echo "               integer system"
  echo "               ExitStatus=system(JobString)"
fi

echo "          end subroutine system_subr"
echo " "
echo "          subroutine wait_subr(pid,ExitStatus)"

if test "x${FC}" == "xnag"; then :
  echo "               use f90_unix_proc, only: wait"
fi

echo "               implicit none"
echo " "
echo "               integer,intent(inout) :: pid,ExitStatus"
echo " "

if test "x${FC}" == "xifort"; then :
   echo "               integer shiftr,ierror"
   echo " "
   echo "               call pxfwait(ExitStatus,pid,ierror)"
elif test "x${FC}" == "xpgi"; then :
   echo "               integer wait"
   echo " "
   echo "               pid=wait(ExitStatus)"
   echo "               ExitStatus=ishft(ExitStatus,-8)"
elif test "x${FC}" == "xnag"; then :
  echo "               call wait(ExitStatus,pid)"
  echo "               ExitStatus=ishft(ExitStatus,-8)"
elif test "x${FC}" == "xpathscale"; then :
  echo "               INTEGER WAIT"
  echo "               pid=wait(ExitStatus)"
fi

echo "          end subroutine wait_subr"
echo " "
echo "          subroutine getpid_subr(pid)"

if test "x${FC}" == "xnag"; then :
  echo "               use f90_unix, only: getpid"
fi

echo "               implicit none"
echo " "
echo "               integer,intent(out) :: pid"
echo " "

if test "x${FC}" == "xifort" -o "x${FC}" == "xpgi" -o "x${FC}" == "xifc" -o "x${FC}" == "xpathscale"; then :
  echo "               integer getpid"
  echo " "
fi

echo "               pid=getpid()"
echo "          end subroutine getpid_subr"
echo "END MODULE PORFUNCS"
