!***************************************************************************
! ts_extract_subperiod_indexes.f90
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
! \brief       From a time series, retrieve indexes of a selected subperiod
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine tsExtractSubperiodIndexes(TimeSeries, nrow, StartTimestamp, &
    EndTimestamp, StartIndex, EndIndex)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    type (DateType), intent(in) :: TimeSeries(nrow)
    type (DateType), intent(in) :: StartTimestamp
    type (DateType), intent(in) :: EndTimestamp
    integer, intent(out) :: StartIndex
    integer, intent(out) :: EndIndex
    !> Local variables
    integer :: i

    StartIndex = nint(error)
    EndIndex = nint(error)

    !> Returns errors if periods do not overlap
    if (EndTimestamp < TimeSeries(1) &
        .or. StartTimestamp > TimeSeries(nrow)) return

    !> Search StartTimestamp
    if (StartTimestamp < TimeSeries(1))  then
        StartIndex = 1
    else
        do i = 1, nrow
            if (TimeSeries(i) >= StartTimestamp) then
                StartIndex = i
                exit
            end if
        end do
    end if

    !> Search EndTimestamp
    if (EndTimestamp >= TimeSeries(nrow)) then
        EndIndex = nrow - 1
    else
        do i = StartIndex, nrow
            if (TimeSeries(i) > EndTimestamp) then
                EndIndex = i - 1
                exit
            end if
        end do
    end if
end subroutine tsExtractSubperiodIndexes
