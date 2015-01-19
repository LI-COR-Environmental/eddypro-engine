!***************************************************************************
! test_absolute_limits.f90
! -------------------------
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
    integer :: i = 0
    integer :: j = 0
    integer :: hflags(gas4)
    real(kind = dbl) :: HorVel(N)

    if (printout) write(*, '(a)', advance = 'no') '   Absolute limits test..'

    !> initializations
    hflags = 9

    !> Checks for values out of bounds in the whole dataset
    do i = 1, N
        !>  Wind components
        HorVel(i) = sqrt((Set(i, u)**2) + (Set(i, v)**2))
        if (HorVel(i) > al%u_max) then
            hflags(u) = 1
            hflags(v) = 1
            if(RPsetup%filter_al) Set(i, u) = error
            if(RPsetup%filter_al) Set(i, v) = error
        else
            hflags(u) = 0
            hflags(v) = 0
        end if
        if(abs(Set(i, w)) > al%w_max) then
            hflags(w) = 1
            if(RPsetup%filter_al) Set(i, w) = error
        else
            hflags(w) = 0
        end if
        !>  Sonic temperature
        if ((Set(i, ts) < al%t_min + 273.15d0) .or. &
            (Set(i, ts) > al%t_max + 273.15d0) ) then
            hflags(ts) = 1
            if(RPsetup%filter_al) Set(i, ts) = error
        else
            hflags(ts) = 0
        end if

        !> co2, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
        if (E2Col(co2)%present) then
            if(E2Col(co2)%measure_type == 'molar_density') then
                if( (Set(i, co2) * StdVair * 1d3 < al%co2_min) .or. &
                    (Set(i, co2) * StdVair * 1d3 > al%co2_max) ) then
                    hflags(co2) = 1
                else
                    hflags(co2) = 0
                end if
            else
                if( (Set(i, co2) < al%co2_min) .or. &
                    (Set(i, co2) > al%co2_max) ) then
                    hflags(co2) = 1
                else
                    hflags(co2) = 0
                end if
            end if
        end if

        !> h2o, expressed as [mmol m-3] if molar_density, [mmol mol-1] otherwise
        if (E2Col(h2o)%present) then
            if(E2Col(h2o)%measure_type == 'molar_density') then
                if( (Set(i, h2o) * StdVair < al%h2o_min) .or. &
                    (Set(i, h2o) * StdVair > al%h2o_max) ) then
                    hflags(h2o) = 1
                else
                    hflags(h2o) = 0
                end if
            else
                if( (Set(i, h2o) < al%h2o_min) .or. &
                    (Set(i, h2o) > al%h2o_max) ) then
                    hflags(h2o) = 1
                else
                    hflags(h2o) = 0
                end if
            end if
        end if

        !> ch4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
        if (E2Col(ch4)%present) then
            if(E2Col(ch4)%measure_type == 'molar_density') then
                if( (Set(i, ch4) * StdVair * 1d3 < al%ch4_min) .or. &
                    (Set(i, ch4) * StdVair * 1d3 > al%ch4_max) ) then
                    hflags(ch4) = 1
                else
                    hflags(ch4) = 0
                end if
            else
                if( (Set(i, ch4) < al%ch4_min) .or. &
                    (Set(i, ch4) > al%ch4_max) ) then
                    hflags(ch4) = 1
                else
                    hflags(ch4) = 0
                end if
            end if
        end if

        !> gas4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
        if (E2Col(gas4)%present) then
            if(E2Col(gas4)%measure_type == 'molar_density') then
                if( (Set(i, gas4) * StdVair * 1d3 < al%gas4_min) .or. &
                    (Set(i, gas4) * StdVair * 1d3 > al%gas4_max) ) then
                    hflags(gas4) = 1
                else
                    hflags(gas4) = 0
                end if
            else
                if( (Set(i, gas4) < al%gas4_min) .or. &
                    (Set(i, gas4) > al%gas4_max) ) then
                    hflags(gas4) = 1
                else
                    hflags(gas4) = 0
                end if
            end if
        end if
    end do
    !> Create an 8-digits number containing the values of the hflags
    IntHF%al = 900000000
    do j = 1, gas4
        IntHF%al = IntHF%al + hflags(j) * (10**(gas4 - j))
    end do
     if (printout) write(*,'(a)') ' Done.'
end subroutine TestAbsoluteLimits
