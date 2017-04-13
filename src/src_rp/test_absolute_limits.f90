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
    integer :: i = 0
    integer :: j = 0
    integer :: cnt1 = 0
    integer :: cnt2 = 0
    integer :: hflags(gas4)
    real(kind = dbl) :: HorVel


    if (printout) write(*, '(a)', advance = 'no') '   Absolute limits test..'

    !> initializations
    hflags = 0

    !> Flag and filter wind components
    cnt1 = 0
    cnt2 = 0
    do i = 1, N
        if (all(Set(i, u:w) /= error)) then
            !> Horizontal wind
            HorVel = sqrt((Set(i, u)**2) + (Set(i, v)**2))
            if (HorVel > al%u_max) then
                cnt1 = cnt1 + 1
                hflags(u) = 1
                hflags(v) = 1
                if(RPsetup%filter_al) then
                    Set(i, u) = error
                    Set(i, v) = error
                end if
            end if
            !> Vertical wind
            if (abs(Set(i, w)) > al%w_max) then
                cnt2 = cnt2 + 1
                hflags(w) = 1
                if(RPsetup%filter_al) Set(i, w) = error
            end if
        end if
    end do
    Essentials%al_s(u) = cnt1
    Essentials%al_s(v) = cnt1
    Essentials%al_s(w) = cnt2

    !> Flag and filter sonic temperature
    Essentials%al_s(ts) = count(Set(:, ts) /= error .and. &
                                (Set(:, ts) < al%t_min + 273.15d0 .or. &
                                 Set(:, ts) > al%t_max + 273.15d0))
    if (Essentials%al_s(ts) > 0) hflags(ts) = 1
    if(RPsetup%filter_al) then
        where (Set(:, ts) < al%t_min + 273.15d0 .or. &
            Set(:, ts) > al%t_max + 273.15d0)
            Set(:, ts) = error
        end where
    end if

    !> Flag co2, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(co2)%present) then
        if(E2Col(co2)%measure_type == 'molar_density') then
            Essentials%al_s(co2) = count(Set(:, co2) /= error .and. &
                                         (Set(:, co2) * StdVair * 1d3 < al%co2_min .or. &
                                          Set(:, co2) * StdVair * 1d3 > al%co2_max))
        else
            Essentials%al_s(co2) = count(Set(:, co2) /= error .and. &
                                         (Set(:, co2) < al%co2_min .or. &
                                          Set(:, co2) > al%co2_max))
        end if
        if (Essentials%al_s(co2) > 0) hflags(co2) = 1
    else
        Essentials%al_s(co2) = ierror
        hflags(co2) = 9
    end if

    !> Flag h2o, expressed as [mmol m-3] if molar_density, [mmol mol-1] otherwise
    if (E2Col(h2o)%present) then
        if(E2Col(h2o)%measure_type == 'molar_density') then
            Essentials%al_s(h2o) = count(Set(:, h2o) /= error .and. &
                                         (Set(:, h2o) * StdVair < al%h2o_min .or. &
                                          Set(:, h2o) * StdVair > al%h2o_max))
        else
            Essentials%al_s(h2o) = count(Set(:, h2o) /= error .and. &
                                         (Set(:, h2o) < al%h2o_min .or. &
                                          Set(:, h2o) > al%h2o_max))
        end if
        if (Essentials%al_s(h2o) > 0) hflags(h2o) = 1
    else
        Essentials%al_s(h2o) = ierror
        hflags(h2o) = 9
    end if

    !> Flag ch4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(ch4)%present) then
        if(E2Col(ch4)%measure_type == 'molar_density') then
            Essentials%al_s(ch4) = count(Set(:, ch4) /= error .and. &
                                         (Set(:, ch4) * StdVair * 1d3 < al%ch4_min .or. &
                                          Set(:, ch4) * StdVair * 1d3 > al%ch4_max))
        else
            Essentials%al_s(ch4) = count(Set(:, ch4) /= error .and. &
                                         (Set(:, ch4) < al%ch4_min .or. &
                                          Set(:, ch4) > al%ch4_max))
        end if
        if (Essentials%al_s(ch4) > 0) hflags(ch4) = 1
    else
        Essentials%al_s(ch4) = ierror
        hflags(ch4) = 9
    end if

    !> Flag gas4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
    if (E2Col(gas4)%present) then
        if(E2Col(gas4)%measure_type == 'molar_density') then
            Essentials%al_s(gas4) = count(Set(:, gas4) /= error .and. &
                                         (Set(:, gas4) * StdVair * 1d3 < al%gas4_min .or. &
                                          Set(:, gas4) * StdVair * 1d3 > al%gas4_max))
        else
            Essentials%al_s(gas4) = count(Set(:, gas4) /= error .and. &
                                         (Set(:, gas4) < al%gas4_min .or. &
                                          Set(:, gas4) > al%gas4_max))
        end if
        if (Essentials%al_s(gas4) > 0) hflags(gas4) = 1
    else
        Essentials%al_s(gas4) = ierror
        hflags(gas4) = 9
    end if

    !> Create an 8-digits number containing the values of the hflags
    IntHF%al = 900000000
    do j = 1, gas4
        IntHF%al = IntHF%al + hflags(j) * (10**(gas4 - j))
    end do

    if (printout) write(*,'(a)') ' Done.'
end subroutine TestAbsoluteLimits
