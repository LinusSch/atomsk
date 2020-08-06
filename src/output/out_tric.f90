MODULE out_tric
!
!**********************************************************************************
!*  OUT_TRIC
!**********************************************************************************
!* This module reads in an arrays containing atomic positions, and                *
!* writes a partial TRIC input file with the ending .TRC.                         *
!*                                                                                *
!* This file format is not canonically described anywhere, as far as I'm aware,   *
!* but can be created by the GUI scripts of TRIC and read by the Fortran core of  *
!* the same. The TRIC program is written by Vasily Khodyrev and described in      *
!*     Khodyrev et.al.: The shower approach in the simulation of ion scattering   *
!*                      from solids                                               *
!* available at                                                                   *
!*     researchgate.net/publication/51466825                                      *
!*                                                                                *
!**********************************************************************************
!* Based on OUT_XMD:                                                              *
!* (C) Nov. 2013 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 March 2014                                    *
!*                                                                                *
!* Conversion into OUT_TRIC:                                                      *
!* (C) 2020-08 - Linus Schönström                                                 *
!*     Uppsala University                                                         *
!*     linus.schonstrom@physics.uu.se                                             *
!* Last modified by Linus Schönström on 2020-08-06                                *
!*                                                                                *
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
!
!
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_TRIC(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j, Nspecies
INTEGER:: typecol, masscol !index of charges, types, mass in AUX
REAL(dp):: smass
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
masscol = 0
typecol = 0
!
msg = 'entering WRITE_TRIC'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Determine if atom positions are in reduced coordinates
CALL FIND_IF_REDUCED(P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Determine how many different species are present
CALL FIND_NSP(P(:,4),atypes)
!
!
100 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=1000)
!Write comment(s)
DO i=1,SIZE(comment)
  WRITE(40,'(a)') '#'//TRIM(comment(1))
ENDDO
!
!
!If box is not orthogonal, warn the user
IF( DABS(H(1,2))>1.d-6 .OR. DABS(H(1,3))>1.d-6 .OR.  &
  & DABS(H(2,1))>1.d-6 .OR. DABS(H(2,3))>1.d-6 .OR.  &
  & DABS(H(3,1))>1.d-6 .OR. DABS(H(3,2))>1.d-6      ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(3709,(/""/),(/0.d0/))
ENDIF
!
!Write box size
IF( isreduced ) THEN
  !All atom positions will be scaled to box size
  WRITE(40,'(a10,3(f12.6,1X))') 'BOX SCALE ', H(1,1), H(2,2), H(3,3)
ELSE
  WRITE(40,'(a4,3(f12.6,1X))') 'BOX ', H(1,1), H(2,2), H(3,3)
ENDIF
!
!
!Write atom positions
WRITE(msg,*) SIZE(P,1)
WRITE(40,'(a)') "POSITION "//TRIM(ADJUSTL(msg))
!
DO i=1,SIZE(P,1)
  IF(typecol.NE.0) THEN
    !use the defined types
    Nspecies = NINT(AUX(i,typecol))
  ELSE
    !Replace atomic number by atom type
    DO j=1,SIZE(atypes,1)
      IF( atypes(j,1)==INT(P(i,4)) ) Nspecies = j
    ENDDO
  ENDIF
  !
  WRITE(temp,150) Nspecies, P(i,1), P(i,2), P(i,3)
  !Write line to file
  WRITE(40,'(a)') TRIM(ADJUSTL(temp))
ENDDO
150 FORMAT(i2,1X,6(f14.8,1X))
!
!
!Write atom masses and species
DO i=1,SIZE(atypes,1)
  WRITE(40,'(a12,i2)') "SELECT TYPE ", i
  CALL ATOMSPECIES(atypes(i,1),species)
  IF(masscol>0) THEN
    !Atom mass is defined as an auxiliary property
    !Find an atom of this type and use its mass
    smass=0.d0
    j=0
    DO WHILE(smass<=0.d0)
      j=j+1
      IF( NINT(P(j,4)) == NINT(atypes(i,1)) ) THEN
        smass = AUX(j,masscol)
      ENDIF
    ENDDO
  ELSE
    !Compute atom mass
    CALL ATOMMASS(species,smass)
  ENDIF
  WRITE(40,'(a5,f12.3)') "MASS ", smass
  WRITE(40,'(a9,i2,1X,a2)') "TYPENAME ", i, species
ENDDO
!
!
!
200 CONTINUE
CLOSE(40)
!
msg = "TRIC"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
1000 CONTINUE
!
!
END SUBROUTINE WRITE_TRIC
!
END MODULE out_tric
