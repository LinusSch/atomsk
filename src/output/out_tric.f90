MODULE out_tric
!
!**********************************************************************************
!*  OUT_TRIC                                                                      *
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
!* This source file is liberally commented, as the author intends to maintain it  *
!* as needed despite not expecting to maintain his (limited) Fortran knowledge.   *
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



USE atoms        !provides subroutines ATOMMASS, ATOMNUMBER and ATOMSPECIES
USE comv         !global variables, e.g. dp (64 bit reals)
USE messages     !provides subroutine ATOMSK_MSG
USE subroutines  !provides subroutines FIND_IF_REDUCED and FIND_NSP
                 !this file also calls INVMAT(), a subroutine defined in module math, I guess chain loading?
IMPLICIT NONE


CONTAINS   !just the one
    SUBROUTINE WRITE_TRIC(H,P,comment,AUXNAMES,AUX,outputfile)

        !Input variables
        CHARACTER(LEN=*),INTENT(IN):: outputfile                          !Filename to write to
        REAL(dp),DIMENSION(3,3),INTENT(IN):: H                            !Base vectors of the supercell
        REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P                !Positions
        REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX              !Auxiliary properties
        CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
        CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment  !unused as this format does not support a comment



        !Variables for reduction of coordinates, if needed
        LOGICAL:: isreduced
        REAL(dp),DIMENSION(3,3):: G      !inverse of the cell (H)
        REAL(dp):: P1, P2, P3

        !Variables I'm not quite sure how they're used yet
        REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
        INTEGER:: typecol, masscol       !index of type and mass column in AUX property array

        !Variables holding information on the current atom type or the current atom
        INTEGER:: i, j, atype            !loop indices, atom type
        REAL(dp):: mass                  !atom mass
        CHARACTER(LEN=2):: species       !atom symbol

        !Variables used in formatting and writing
        CHARACTER(LEN=4096):: msg, temp  !string variables: mutable, fixed size, gets filled with blanks often
        CHARACTER(LEN=128):: header, header_template, temp2, note
        INTEGER:: left, right, center    !used in formatting header lines
        INTEGER:: align_note_r           !used so that notes can be inserted anywhere in the line



        !Formatting parameters, but caution: if notes or comments are set too far left they may overwrite data!
        !
        !The line with lattice parameters is the pain point, the number of
        !digits before decimal point in a moves the end of that data.
        !
        !At some point I should parameterize alignment of notes in the atom
        !list as well, that is buried in a format specifier for now.
        !
        INTEGER,PARAMETER:: indent= 4                 ! indent of header lines
        INTEGER,PARAMETER:: align_c= 3*10 +1          ! place all comments this far into the line
        INTEGER,PARAMETER:: align_note= align_c - 5   ! place type notes this far into the line

        !Debug message to keep track of execution
        msg = 'entering WRITE_TRIC'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !
        !Initialize variables
        masscol = 0
        typecol = 0
        isreduced = .FALSE.

        !If the cell is not orthogonal, this file cannot be written
        IF( DABS(H(1,2))>1.d-6 .OR. DABS(H(1,3))>1.d-6 .OR.  &
            & DABS(H(2,1))>1.d-6 .OR. DABS(H(2,3))>1.d-6 .OR.  &
            & DABS(H(3,1))>1.d-6 .OR. DABS(H(3,2))>1.d-6      ) THEN
            nwarn=nwarn+1
            CALL ATOMSK_MSG(3719,(/"TRIC"/),(/0.d0/))
            GOTO 1000  !end of subroutine
        ENDIF

        !Determine if atom positions are in reduced coordinates
        CALL FIND_IF_REDUCED(P,isreduced)
        WRITE(msg,*) 'isreduced:', isreduced
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !if not, calculate inverse of matrix H for use in reducing when writing
        IF( .NOT.isreduced ) THEN
            msg = 'inverting matrix H'
            CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
            CALL INVMAT(H,G)
        ENDIF

        !Determine how many different species are present
        CALL FIND_NSP(P(:,4),atypes)

        !Prepare some template strings
        header_template= repeat(' ',indent) // repeat('=',64)
        center= len_trim(adjustL(header_template)) / 2 + indent


        !********************************
        !* Opening the file for writing *
        !********************************
        OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=1000)

        !Write section header line
        header= 'ATOMIC COMPOSITION'
        left= center - len_trim(header) / 2
        right= left + len_trim(header) + 1
        temp= header_template
        temp(left:right)= ' '//header//' ' !replaces the indexed characters with header, concatenated with surrounding spaces
        write(40,'(a)') trim(temp)         !a format specifier with a single string avoids a leading space created by the
                                           !right-adjusted list-oriented write(unit,*)

        !Write number of types
        write(temp2,*) SIZE(atypes,1)      !writing to a string variable is how you convert from number to string
        temp= adjustL(temp2)
        temp(align_c:)= 'Number of different atomic sites'  !replacing characters with open-ended indexing fills with blanks
        write(40,'(a)') trim(temp)

        !Write atom info for each type
        DO i=1,SIZE(atypes,1)              !looping over the length of the atypes array

            !create note
            write(temp2,*) i
            note= '('//trim(adjustL(temp2))//')'
            align_note_r= align_note + len_trim(note)

            !get symbol
            CALL ATOMSPECIES(atypes(i,1),species)

            !write symbol
            write(temp2,*) species
            temp= adjustL(temp2)
            temp(align_c:)= 'Atom symbol'
            temp(align_note:align_note_r)= note
            write(40,'(a)') trim(temp)

            !write number
            write(temp2,*) int(atypes(i,1))
            temp= adjustL(temp2)
            temp(align_c:)= 'Atom number'
            temp(align_note:align_note_r)= note
            write(40,'(a)') trim(temp)

            !get mass
            IF(masscol>0) THEN
                !Atom mass is defined as an auxiliary property
                !Find an atom of this type and use its mass
                mass=0.d0
                j=0
                DO WHILE(mass<=0.d0)
                    j=j+1
                    IF( NINT(P(j,4)) == NINT(atypes(i,1)) ) THEN
                        mass = AUX(j,masscol)
                    ENDIF
                ENDDO
            ELSE
                !get atom mass from atoms.f90
                CALL ATOMMASS(species,mass)
            ENDIF

            !write mass
            write(temp2,'(f10.6)') mass
            temp= adjustL(temp2)
            temp(align_c:)= 'Atom mass'
            temp(align_note:align_note_r)= note
            write(40,'(a)') trim(temp)

            !write thermal vibrations placeholder
            write(temp2,'(3(f5.3,1X))') 0.0, 0.0, 0.0   ! the '1X' in the format string inserts one space
            temp= adjustL(temp2)
            temp(align_c:)= 'Amplitudes of thermal vib., Angstrom'
            temp(align_note:align_note_r)= note
            write(40,'(a)') trim(temp)

        ENDDO  !done with per-type info



        !Write section header line
        header= 'CRYSTAL STRUCTURE'
        left= center - len_trim(header) / 2
        right= left + len_trim(header) + 1
        temp= header_template
        temp(left:right)= ' '//header//' '
        write(40,'(a)') trim(temp)

        !Write lattice type line
        temp= '5'//repeat(' ',5)//'(Custom)'
        temp(align_c:)= 'Lattice type'
        write(40,'(a)') trim(temp)

        !Write (super) cell size
        write(temp2,'(3(f9.4,1X))') H(1,1), H(2,2), H(3,3)
        temp= adjustL(temp2)
        temp(align_c:)= 'Cell size (lattice parameters), Angstrom'
        write(40,'(a)') trim(temp)

        !Write number of atoms
        write(temp2,*) SIZE(P,1)
        temp= adjustL(temp2)
        temp(align_c:)= 'Number of atoms in cell'
        write(40,'(a)') trim(temp)

        !Write list of atom positions
        DO i=1,SIZE(P,1)

            !get type identifier
            IF(typecol.NE.0) THEN
                !use the defined types
                atype = NINT(AUX(i,typecol))
            ELSE
                !Replace atomic number by atom type
                DO j=1,SIZE(atypes,1)
                    IF( atypes(j,1)==INT(P(i,4)) ) atype = j
                ENDDO
            ENDIF

            !get symbol
            CALL ATOMSPECIES(P(i,4),species)

            !write identifier line
            WRITE(temp,'(i2,5X,a4)') atype, '('//species//')'
            WRITE(40,'(a)') TRIM(ADJUSTL(temp))

            !create position line with reduced coordinates
            !yanked from out_crystal.f90 by Pierre Hirel, with slight modification
            IF(isreduced) THEN
                WRITE(temp,'( 3(f11.8,2X) )') P(i,1), P(i,2), P(i,3)
            ELSE
                P1 = P(i,1)
                P2 = P(i,2)
                P3 = P(i,3)
                WRITE(temp,'( 3(f11.8,2X) )')  &
                      &  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                      &  P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                      &  P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
            ENDIF

            !write position line
            write(40,'(a)') trim(temp)
        ENDDO

        !********************
        !* Closing the file *
        !********************
        CLOSE(40)


        !Success message
        msg = "TRIC"
        temp = outputfile
        CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))

        !The end
        1000 CONTINUE
    END SUBROUTINE WRITE_TRIC
END MODULE out_tric
