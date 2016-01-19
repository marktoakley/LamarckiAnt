#!/bin/csh -f

# GPL License Info  {{{
#
#   Portability functions module generator
#   Copyright (C) 2003-2005 Semen A. Trygubenko and David J. Wales
#   This file is part of OPTIM.
#
#   OPTIM is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   OPTIM is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#}}}

# read in compiler info from command-line input {{{
if ( $#argv == 1) then
  set compiler = $argv[1]
  if ( $compiler != 'nag' && $compiler != 'ifort' && $compiler != 'pgi' && $compiler != 'ifc' && $compiler != 'g95'  && $compiler != 'gfortran' && $compiler != 'pathscale' ) then
    echo 'Unknown compiler!'
    exit 1
  endif
else
  set compiler = 'pgi'
endif
#}}}

echo "MODULE PORFUNCS"

if ( $compiler == 'nag' ) then
  echo "     use f90_unix, only: getarg"
  echo "     use f90_unix_proc, only: system, exit"
endif

echo "     implicit none"
echo "     contains"
# FLUSH {{{
if ( $compiler == 'nag' ) then
  echo "          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be"
  echo "          use f90_unix, NAGflush => flush ! connected for formatted sequential output; ISTAT is ignored"
  echo "               implicit none"
  echo "               integer,intent(in) :: UNIT"
  echo "               integer,intent(out),optional :: ISTAT"
  echo "               call NAGflush(UNIT)"
  echo "          end subroutine flush"
  echo " "
else if ( $compiler == 'g95' ) then
  echo "          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be"
  echo "               implicit none ! connected for formatted sequential output; ISTAT is ignored"
  echo "               integer,intent(in) :: UNIT"
  echo "               integer,intent(out),optional :: ISTAT"
  echo "               call flush(UNIT)"
  echo "          end subroutine flush"
  echo " "
else if ( $compiler == 'gfortran' ) then
  echo "          subroutine flush(UNIT,ISTAT) ! flushes the output buffer of logical unit UNIT which must be"
  echo "               implicit none ! connected for formatted sequential output; ISTAT is ignored"
  echo "               integer,intent(in) :: UNIT"
  echo "               integer,intent(out),optional :: ISTAT"
  echo "!              call flush(UNIT)"
  echo "          end subroutine flush"
  echo " "
endif
#}}}
#op226> getarg  {{{
echo "          subroutine getarg_subr(position,value) ! wraps getarg function so it can be use-associated"

#nag{{{
if ( $compiler == 'nag' ) then
  echo "               use f90_unix, only: getarg"
endif
#}}}

echo "               implicit none"
echo " "
echo "               integer,intent(in) :: position"
echo "               character(len=*),intent(out) :: value"
echo " "

# not nag{{{
if ( $compiler != 'nag' ) then
  #echo "               integer getarg"
  echo " "
endif 
#}}}

echo "               call getarg(position,value)"
echo "          end subroutine getarg_subr"
echo " "
#}}}
# IARGC {{{
echo "          subroutine iargc_subr(n) ! wraps iargc function so it can be use-associated"

if ( $compiler == 'nag' ) then
  echo "               use f90_unix, only: iargc"
endif

echo "               implicit none"
echo "               integer,intent(out) :: n"
echo " "

if ( $compiler != 'nag' ) then
  echo "               integer iargc"
  echo " "
endif 

echo "               n = iargc()"
echo "          end subroutine iargc_subr"
echo " "
#}}}
# FORK {{{
echo "          subroutine fork_subr(pid)" ! returns zero in the child process, PID of child in parent process
if ( $compiler == 'nag' ) then
  echo "               use f90_unix_proc, only: fork"
endif
echo "               implicit none"
echo "               integer, intent(inout) :: pid"
echo " "
if ( $compiler == 'pgi' ) then
     echo "               integer fork"
endif
echo " "
if ( $compiler == 'pgi' ) then
     echo "               pid=fork()"
else if ( $compiler == 'ifort' ) then
     echo "               integer ierror"
     echo "               call pxffork(pid,ierror)"
else if ( $compiler == 'nag' ) then
     echo "               call fork(pid)"
else if ( $compiler == 'pathscale' ) then
     echo "               integer fork"
     echo "               pid=fork()"
else if ( $compiler == 'gfortran' ) then
endif
echo "          end subroutine fork_subr"
echo " "
#}}}
# SYSTEM {{{
echo "          subroutine system_subr(JobString,ExitStatus)"
if ( $compiler == 'nag' ) then
  echo "               use f90_unix_proc, only: system"
endif
echo "               implicit none"
echo " "
echo "               character(len=*),intent(in) :: JobString"
echo "               integer,intent(out) :: ExitStatus"
echo " "
if ( $compiler == 'ifort' ) then
  echo "               integer shiftr,system"
  echo " "
  echo "               ExitStatus=system(JobString)"
#  echo "               ExitStatus=shiftr(ExitStatus,-8)"
else if ( $compiler == 'pgi' ) then
  echo "               integer system"
  echo " "
     echo "               ExitStatus=system(JobString)"
     echo "               ExitStatus=ishft(ExitStatus,-8)"
else if ( $compiler == 'nag' ) then
  echo "               call system(JobString,ExitStatus)"
  echo "               ExitStatus=ishft(ExitStatus,-8)"
else if ( $compiler == 'pathscale' ) then
  echo "               integer system"
  echo "               ExitStatus=system(JobString)"
endif
echo "          end subroutine system_subr"
echo " "
#}}}
# WAIT {{{
echo "          subroutine wait_subr(pid,ExitStatus)"
if ( $compiler == 'nag' ) then
  echo "               use f90_unix_proc, only: wait"
endif
echo "               implicit none"
echo " "
echo "               integer,intent(inout) :: pid,ExitStatus"
echo " "
if ( $compiler == 'ifort' ) then
   echo "               integer shiftr,ierror"
   echo " "
   echo "               call pxfwait(ExitStatus,pid,ierror)"
else if ( $compiler == 'pgi' ) then
   echo "               integer wait"
   echo " "
   echo "               pid=wait(ExitStatus)"
   echo "               ExitStatus=ishft(ExitStatus,-8)"
else if ( $compiler == 'nag' ) then
  echo "               call wait(ExitStatus,pid)"
  echo "               ExitStatus=ishft(ExitStatus,-8)"
else if ( $compiler == 'pathscale' ) then
  echo "               INTEGER WAIT"
  echo "               pid=wait(ExitStatus)"
else if ( $compiler == 'gfortran' ) then
endif
echo "          end subroutine wait_subr"
echo " "
#}}}
# GETPID {{{
echo "          subroutine getpid_subr(pid)"
if ( $compiler == 'nag' ) then
  echo "               use f90_unix, only: getpid"
endif
echo "               implicit none"
echo " "
echo "               integer,intent(out) :: pid"
echo " "
if ( $compiler == 'ifort' || $compiler == 'pgi' || $compiler == 'ifc' || $compiler == 'pathscale' ) then
  echo "               integer getpid"
  echo " "
endif
echo "               pid=getpid()"
echo "          end subroutine getpid_subr"
#}}}
echo "END MODULE PORFUNCS"
