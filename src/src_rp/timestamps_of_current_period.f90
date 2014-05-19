!***************************************************************************
! timestamps_of_current_period.f90
! --------------------------------
!Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
!Copyright (C) 2011, LI-COR Biosciences
!
!This file is part of EddyPro (TM).
!
!EddyPro (TM) is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!EddyPro (TM) is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!***************************************************************************

!***************************************************************************
! \file        src/timestamps.f90
! \brief       Defines initial and final timestamps of current period
! \version     4.1.0
! \date        2012-10-11
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TimestampsOfCurrentPeriod(InitialFile, InitialTimestamp, FinalTimestamp)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: InitialFile
    type (DateType), intent(out) :: InitialTimestamp
    type (DateType), intent(out) :: FinalTimestamp
    !> local variables
    character(10) :: file_date
    character(5) :: file_time
    character(10) :: date
    character(5) :: time

    !> Grab timestamp written on file
    call FilenameToDateTime(InitialFile, EddyProProj%fproto, file_date, file_time)

    !> Takes it to the beginning of the period if timestamp on file refers to end of period
    if (EddyProLog%tstamp_end) then
        call SubtractDateStep(file_date, file_time, date, time, DatafileDateStep)
    else
        date = file_date
        time = file_time
    end if

    !> Convert into initial timestamp
    call DateTimeToDateType(date, time, InitialTimestamp)

    !> Determine end of current averaging interval
    call AddDateStep(date, time, date, time, DateStep)
    call DateTimeToDateType(date, time, FinalTimestamp)

end subroutine TimestampsOfCurrentPeriod
