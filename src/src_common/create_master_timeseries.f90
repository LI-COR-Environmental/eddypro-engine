!***************************************************************************
! create_master_timeseries.f90
! ----------------------------
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
! \brief       Create timestamp array based on start/end timestamp and step
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateTimeSeries(StartTimestamp, EndTimestamp, &
    Step, RawTimeSeries, nrow)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(inout) :: EndTimestamp
    type(DateType), intent(in) :: Step
    type(DateType), intent(out) :: RawTimeSeries(nrow)
    !> in/out variables
    integer :: cnt


    write(*, '(a)', advance = 'no') ' Creating time series &
        &for the time period covered by raw data files.. '

    !> create master timestamps array
    RawTimeSeries(1) = StartTimestamp
    cnt = 1
    do
        if (RawTimeSeries(cnt) >= EndTimestamp) exit
        cnt = cnt + 1
        if (cnt > nrow) exit
        RawTimeSeries(cnt) = RawTimeSeries(cnt - 1) + Step
    end do

    write(*, '(a)') ' done.'
end subroutine CreateTimeSeries

!***************************************************************************
!
! \brief       Calculate number of periods in time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function NumOfPeriods(StartTimestamp, EndTimestamp, Step)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(inout) :: EndTimestamp
    type(DateType), intent(in) :: Step
    !> Local variables
    integer :: cnt
    type (DateType) :: Timestamp


    Timestamp = StartTimestamp
    cnt = 1
    do
        if (Timestamp >= EndTimestamp) exit
        cnt = cnt + 1
        Timestamp = Timestamp + Step
    end do
    NumOfPeriods = cnt - 1
end function NumOfPeriods
