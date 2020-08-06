MODULE comv
!
!**********************************************************************************
!*  COMV                                                                          *
!**********************************************************************************
!* This module contains the global variables used by Atomsk.                      *
!* Global variables should be limited to the strict minimum, and NEVER            *
!* include variables or arrays about atomic systems (atom positions etc.).        *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: by Linus Schönström on 2020-08-06                           *
!**********************************************************************************
!* This program is free software: you can redistribute it and/or modify           *
!* it under the terms of the GNU General Public License as published by           *
!* the Free Software Foundation, either version 3 of the License, or              *
!* (at your option) any later version.                                            *
!*                                                                                *
!* This program is distributed in the hope that it will be useful,                *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
!* GNU General Public License for more details.                                   *
!*                                                                                *
!* You should have received a copy of the GNU General Public License              *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
!**********************************************************************************
!
CHARACTER(LEN=24),PARAMETER:: version = 'master'
INTEGER:: nwarn, nerr  !number of warnings/errors encountered during run
INTEGER,PARAMETER:: il = SELECTED_INT_KIND(9)        !integers up to 10^9
INTEGER,PARAMETER:: dp = SELECTED_REAL_KIND(15,307)  !reals with 64-bits precision
INTEGER(il),PARAMETER:: NATOMS_MAX = HUGE(INT(0,il)) !maximum number of atoms that Atomsk can handle
!
!**********************************
!*  ENVIRONMENT-DEPENDENT VARIABLES
!**********************************
CHARACTER(LEN=3):: lang !language in which the program will run (should be 2 letters)
CHARACTER(LEN=1):: langyes, langBigYes, langno !one-letter shortcuts for "yes" and "no", e.g. "y", "n"
#if defined(WINDOWS)
  CHARACTER(LEN=1),PARAMETER:: pathsep='\'     !path separator for Windows
  CHARACTER(LEN=3),PARAMETER:: system_ls='dir' !command to list current directory
  CHARACTER(LEN=16),PARAMETER:: pathnull='2>nul'  !redirection to NULL for Windows
#else
  CHARACTER(LEN=1),PARAMETER:: pathsep='/'     !path separator for UNIX/Linux
  CHARACTER(LEN=3),PARAMETER:: system_ls='ls ' !command to list current directory
  CHARACTER(LEN=16),PARAMETER:: pathnull='2>/dev/null'  !redirection to NULL for UNIX/Linux
#endif
!
!**********************************
!*  PROGRAM BEHAVIOR
!**********************************
CHARACTER(LEN=128):: logfile !name of logfile for the program
LOGICAL:: overw, ignore      !automatically overwrite/ignore existing files?
INTEGER:: verbosity          !level of verbosity of the program
!
!**********************************
!*  OUTPUT
!**********************************
!The following array contains a list of formats available *FOR OUTPUT* only.
!It should be updated when new formats are added to Atomsk
!Note that each entry must be *exactly* 5 characters long (add spaces if necessary)
!Let's keep to alphabetical order in this list as well (except for exyz and sxyz)
CHARACTER(LEN=5),DIMENSION(28),PARAMETER:: listofformats =                             &
& (/'abin ','atsk ','bop  ','bx   ','cel  ','cfg  ','cif  ','coo  ','csv  ','d12  ',   &
&   'dlp  ','fdf  ','gin  ','imd  ','jems ','lmp  ','mol  ','pos  ','pw   ','stru ',   &
&   'trc  ','vesta','xmd  ','xsf  ','xv   ','xyz  ','exyz ','sxyz '    &
& /)
!
END MODULE comv
