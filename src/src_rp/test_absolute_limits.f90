!***************************************************************************
! test_absolute_limits.f90
! ------------------------
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
! \brief       Checks for data outside realistic ranges, hard-flag \n
!              file accordingly and eliminate those values if requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestAbsoluteLimits(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    logical :: printout
    !> local variables
    integer :: j = 0
    integer :: hflags(gas4)
    real(kind = dbl) :: HorVel(N)

    if (printout) write(*, '(a)', advance = 'no') '   Absolute limits test..'

    !> initializations
    hflags = 0

    !> Flag wind components
    HorVel(:) = sqrt((Set(:, u)**2) + (Set(:, v)**2))
    if (any(HorVel > al%u_max)) then
        hflags(u) = 1
        hflags(v) = 1
    end if
    if (any(abs(Set(:, w)) > al%w_max)) hflags(w) = 1
    !> Flag sonic temperature
    if (any(Set(:, ts) < al%t_min + 273.15d0) .or. &
         any(Set(:, ts) > al%t_max + 273.15d0)) hflags(ts) = 1

    !> Flag co2, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(co2)%present) then
        if(E2Col(co2)%measure_type == 'molar_density') then
            if (any(Set(:, co2) * StdVair * 1d3 < al%co2_min) &
                .or. any(Set(:, co2) * StdVair * 1d3 > al%co2_max)) hflags(co2) = 1
        else
            if (any(Set(:, co2) < al%co2_min) &
                .or. any(Set(:, co2) > al%co2_max)) hflags(co2) = 1
        end if
    end if

    !> Flag h2o, expressed as [mmol m-3] if molar_density, [mmol mol-1] otherwise
    if (E2Col(h2o)%present) then
        if(E2Col(h2o)%measure_type == 'molar_density') then
            if (any(Set(:, h2o) * StdVair * 1d3 < al%h2o_min) &
                .or. any(Set(:, h2o) * StdVair * 1d3 > al%h2o_max)) hflags(h2o) = 1
        else
            if (any(Set(:, h2o) < al%h2o_min) &
                .or. any(Set(:, h2o) > al%h2o_max)) hflags(h2o) = 1
        end if
    end if

    !> Flag ch4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(ch4)%present) then
        if(E2Col(ch4)%measure_type == 'molar_density') then
            if (any(Set(:, ch4) * StdVair * 1d3 < al%ch4_min) &
                .or. any(Set(:, ch4) * StdVair * 1d3 > al%ch4_max)) hflags(ch4) = 1
        else
            if (any(Set(:, ch4) < al%ch4_min) &
                .or. any(Set(:, ch4) > al%ch4_max)) hflags(ch4) = 1
        end if
    end if

    !> Flag gas4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(gas4)%present) then
        if(E2Col(gas4)%measure_type == 'molar_density') then
            if (any(Set(:, gas4) * StdVair * 1d3 < al%gas4_min) &
                .or. any(Set(:, gas4) * StdVair * 1d3 > al%gas4_max)) hflags(gas4) = 1
        else
            if (any(Set(:, gas4) < al%gas4_min) &
                .or. any(Set(:, gas4) > al%gas4_max)) hflags(gas4) = 1
        end if
    end if

    !> Create an 8-digits number containing the values of the hflags
    IntHF%al = 900000000
    do j = 1, gas4
        IntHF%al = IntHF%al + hflags(j) * (10**(gas4 - j))
    end do

    !> Filter sonic measurements for absolute limits
    if(RPsetup%filter_al) then
        where (HorVel > al%u_max)
            Set(:, u) = error
            Set(:, v) = error
        end where
        where (abs(Set(:, w)) > al%w_max)
            Set(:, w) = error
        end where
        where ( Set(:, ts) < al%t_min + 273.15d0 .or. &
            Set(:, ts) > al%t_max + 273.15d0)
            Set(:, ts) = error
        end where
    end if

    if (printout) write(*,'(a)') ' Done.'
end subroutine TestAbsoluteLimits
