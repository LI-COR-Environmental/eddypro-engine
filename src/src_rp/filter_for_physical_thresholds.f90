!***************************************************************************
! filter_for_physical_thresholds.f90
! ----------------------------------
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
! \brief       Apply "absolute limits" filter to variables selected in SelectWhat
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilterForPhysicalThresholds(Set, N, M, FilterWhat)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    logical, intent(in) :: FilterWhat(M)
    real(kind = dbl), intent(inout) :: Set(N, M)


    !> Eliminate for values out of bounds in the whole dataset
    do i = 1, N
        !> co2, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
        if (E2Col(co2)%present .and. FilterWhat(co2)) then
            if(E2Col(co2)%measure_type == 'molar_density') then
                if( (Set(i, co2) * Ambient%Va * 1d3 < al%co2_min) .or. &
                    (Set(i, co2) * Ambient%Va * 1d3 > al%co2_max)) Set(i, co2) = error
            else
                if( (Set(i, co2) < al%co2_min) .or. &
                    (Set(i, co2) > al%co2_max)) Set(i, co2) = error
            end if
        end if

        !> h2o, expressed as [mmol m-3] if molar_density, [mmol mol-1] otherwise
        if (E2Col(h2o)%present .and. FilterWhat(h2o)) then
            if(E2Col(h2o)%measure_type == 'molar_density') then
                if( (Set(i, h2o) * Ambient%Va < al%h2o_min) .or. &
                    (Set(i, h2o) * Ambient%Va > al%h2o_max)) Set(i, h2o) = error
            else
                if( (Set(i, h2o) < al%h2o_min) .or. &
                    (Set(i, h2o) > al%h2o_max)) Set(i, h2o) = error
            end if
        end if

        !> ch4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
        if (E2Col(ch4)%present .and. FilterWhat(ch4)) then
            if(E2Col(ch4)%measure_type == 'molar_density') then
                if( (Set(i, ch4) * Ambient%Va * 1d3 < al%ch4_min) .or. &
                    (Set(i, ch4) * Ambient%Va * 1d3 > al%ch4_max)) Set(i, ch4) = error
            else
                if( (Set(i, ch4) < al%ch4_min) .or. &
                    (Set(i, ch4) > al%ch4_max)) Set(i, ch4) = error
            end if
        end if

        !> gas4, expressed as [mmol m-3] if molar_density, [umol mol-1] otherwise
        if (E2Col(gas4)%present .and. FilterWhat(gas4)) then
            if(E2Col(gas4)%measure_type == 'molar_density') then
                if( (Set(i, gas4) * Ambient%Va * 1d3 < al%gas4_min) .or. &
                    (Set(i, gas4) * Ambient%Va * 1d3 > al%gas4_max)) Set(i, gas4) = error
            else
                if( (Set(i, gas4) < al%gas4_min) .or. &
                    (Set(i, gas4) > al%gas4_max)) Set(i, gas4) = error
            end if
        end if
    end do
end subroutine FilterForPhysicalThresholds
