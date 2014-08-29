!***************************************************************************
! extract_timeperiod_stats.f90
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
! \brief       extract a subset of statistics from a larger set, based on \n
!              timestamps
! \author      Gerardo Fratini
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExtractTimePeriodStats(DateArray, StatsIn, StatsOut, Nin, ndates, Nout)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: Nin
    integer, intent(in) :: ndates
    Type(StatsType), intent(in) :: StatsIn(MaxNumFile)
    type (DateArrayType), intent(in) :: DateArray(MaxNumFile)
    integer, intent(out) :: Nout
    Type(StatsType), intent(out) :: Statsout(MaxNumFile)
    !local variables
    integer :: i
    integer :: j
    character(16) :: InDate

    write(*, '(a)') ' Extracting statistics for selected time period.. '
    call log_msg( ' inf=extracting statistics for selected time period')

    Nout = 0
    ol: do i =1, Nin
        InDate = StatsIn(i)%Date(1:10) // ',' // StatsIn(i)%time(1:5)
        il: do j = 1, ndates
            if (index(InDate(1:len_trim(InDate)), &
                DateArray(j)%text(1:len_trim(DateArray(j)%text))) /= 0) then
                Nout = Nout + 1
                StatsOut(Nout) = StatsIn(i)
            exit il
            end if
        end do il
    end do ol
    if (Nout == 0) then
        call log_msg(' stop=no statistics data found for the selected time period. terminating execution')
        call ErrorHandle(0, 0, 11)
    end if
    write(LogInteger, '(i6)') Nout
    call SchrinkString(LogInteger)
    write(*, '(a)') '  ' // LogInteger(1:len_trim(LogInteger)) // ' stats available for the selected time period.'
    write(*, '(a)') ' done.'
    LogString = ' n_aval_stats=' // LogInteger(1:len_trim(LogInteger))
    call log_msg(LogString)
    call log_msg( ' inf=statistics for the selected time period extracted correctly')
end subroutine ExtractTimePeriodStats
