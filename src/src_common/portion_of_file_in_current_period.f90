!***************************************************************************
! portion_of_file_in_current_period.f90
! -------------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
! \brief       Calculate first and last record of current file that are to be
!              associated to current period
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PortionOfFileInCurrentPeriod(InitialTimestamp, FinalTimestamp, FileInitialTimestamp, FirstRecord, LastRecord)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(DateType), intent(in) :: InitialTimestamp
    type(DateType), intent(in) :: FinalTimestamp
    type(DateType), intent(in) :: FileInitialTimestamp
    integer, intent(out) :: FirstRecord
    integer, intent(out) :: LastRecord


    !> Note:  LastRecord = -1 means: read to the end of the file.

    !> Initialization
    FirstRecord = nint(error)
    LastRecord = nint(error)

    !> Ideal case: File matches flux period
    if (FileInitialTimestamp == InitialTimestamp .and. &
        FileInitialTimestamp + DatafileDateStep == FinalTimestamp) then
        FirstRecord = 1
        LastRecord = -1
        return
    end if

    !> If beginning of file is same or after InitialTimestamp..
    if (FileInitialTimestamp >= InitialTimestamp) then
        !> ..FirstRecord is 1
        FirstRecord = 1
    else
        !> Function TimeLag calculates the time (in fractions of a day) between the beginning of the file
        !> and the BEGINNING of the current period, then converted to number of records
        FirstRecord = nint(TimeLag(InitialTimestamp, FileInitialTimestamp) * 1440d0 * 60d0 * Metadata%ac_freq) + 1
    end if

    !> If end of file is before (or equal to) FinalTimestamp, LastRecord is -1
    if (FileInitialTimestamp + DatafileDateStep <= FinalTimestamp) then
        LastRecord = -1
    else
        !> Function TimeLag calculates the time (in fractions of a day) between the beginning of the file
        !> and the END of the current period, then converted to number of records
        LastRecord = nint(TimeLag(FinalTimestamp, FileInitialTimestamp) * 1440d0 * 60d0 * Metadata%ac_freq)
    end if

end subroutine PortionOfFileInCurrentPeriod
