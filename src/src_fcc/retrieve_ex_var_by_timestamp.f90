!***************************************************************************
! retrieve_ex_var_by_timestamp.f90
! --------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2015, LI-COR Biosciences
!
! This file is part of EddyPro (TM).
!
! EddyPro (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyPro (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!***************************************************************************
!
! \brief       Retrieve "essentials" information from file, based on
!              timestamp provided on input
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveExVarsByTimestamp(unt, Timestamp, lEx, endReached, skip)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    type(DateType), intent(in) :: Timestamp
    logical, intent(out) :: endReached
    logical, intent(out) :: skip
    type(ExType), intent(out) :: lEx
    !> Local variables
    type(DateType) :: ExTimestamp
    logical :: EndOfFileReached
    logical :: ValidRecord

    skip = .false.
    endReached = .false.
    do
        call ReadExRecord('', unt, -1, lEx, ValidRecord, EndOfFileReached)

        !> If end of files was reached, exit routine with error flag
        if (EndOfFileReached) then
            endReached = .true.
            skip = .true.
            return
        end if

        !> If timestamp matches, exit routine (with error flag if the case)
        call DateTimeToDateType(lEx%date, lEx%time, ExTimestamp)
        if (ExTimestamp == Timestamp) then
            if (.not. ValidRecord) skip = .true.
            return
        end if

        !> If timestamp exceeds the one looked for, backwards essentials unit
        !> and exit with error code
        if (ExTimestamp >= Timestamp) then
            backspace(unt)
            skip = .true.
            return
        end if
    end do
end subroutine RetrieveExVarsByTimestamp
