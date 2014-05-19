!***************************************************************************
! write_out_stats.f90
! -------------------
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
! \file        src/write_out_stats.f90
! \brief       Create date array based on start/end timestamp and step
! \version     3.0.0
! \date        2011-03-24
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateDateArray(start_date, start_time, end_date, end_time, &
                           Step, DateArray, ndates)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: start_date
    character(*), intent(in) :: end_date
    character(*), intent(in) :: start_time
    character(*), intent(in) :: end_time
    type(DateType), intent(in) :: Step
    type(DateArrayType), intent(out) :: DateArray(MaxNumDates)
    integer, intent(out) :: ndates
    !> in/out variables
    character(10) :: curr_date
    character(5):: curr_time
    type(DateType) :: DataSetCurrTime
    type(DateType) :: DataSetEndTime

	
    write(*, '(a)', advance = 'no') ' Creating timestamp array for processed time period.. '
    call log_msg( ' inf=creating timestamp array for processed time period')
    LogString = ' start_date=' // start_date
    call log_msg(LogString)
    LogString = ' end_date=' // end_date
    call log_msg(LogString)

    !> Initialize dataset end time
    call DateTimeToDateType(end_date, end_time, DataSetEndTime)
    DataSetEndTime = DataSetEndTime + Step

    !> create dataset timestamp array
    call DateTimeToDateType(start_date, start_time, DataSetCurrTime)
    ndates = 0
    do
        if (DataSetCurrTime > DataSetEndTime) exit
        ndates = ndates + 1
        call DateTypeToDateTime(DataSetCurrTime, curr_date, curr_time)
        DateArray(ndates)%text = curr_date(1:10) // ',' // curr_time(1:5)
        DataSetCurrTime = DataSetCurrTime + Step
    end do
    ndates = ndates - 1
    DateArray(1:ndates)%present = .false.
    write(*, '(a)') ' done.'
    call log_msg( ' inf=timestamp array created correctly.')
end subroutine CreateDateArray
