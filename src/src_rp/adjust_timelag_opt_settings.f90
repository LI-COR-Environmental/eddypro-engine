!***************************************************************************
! adjust_timelag_opt_settings.f90
! -------------------------------
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
! \brief       Adjust time lag opt settings if user did not set or set improperly
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjustTimelagOptSettings()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    real(kind = dbl) :: nominal
    real(kind = dbl) :: mult(GHGNumVar)
    real(kind = dbl) :: gui_tlag_threshold = -1000d0
    real(kind = dbl) :: tube_time(GHGNumVar)
    real(kind = dbl) :: tube_volume(GHGNumVar)
    real(kind = dbl) :: cell_time(GHGNumVar)
    real(kind = dbl) :: cell_volume(GHGNumVar)
    real(kind = dbl) :: safety


    !> Initialization to zero of all timelags
    E2Col(:)%min_tl = 0d0
    E2Col(:)%max_tl = 0d0

    !> Initialize multiplier
    mult(:) = 2d0     !< For passive gases
    mult(h2o) = 10d0  !< For active gases
    safety = 0.3d0    !< Safety margin for min/max setting

    !> Transit time in cell and sampling lines of closed path instruments
    where (E2Col(co2:gas4)%instr%path_type == 'closed')
        tube_volume(co2:gas4) = &
            (p * (E2Col(co2:gas4)%instr%tube_d / 2d0)**2 * &
            E2Col(co2:gas4)%instr%tube_l)
        tube_time(co2:gas4) =  tube_volume(co2:gas4) &
            / E2Col(co2:gas4)%instr%tube_f

        cell_volume(co2:gas4) = &
            (p * (E2Col(co2:gas4)%instr%hpath_length / 2d0)**2 * &
                                E2Col(co2:gas4)%instr%vpath_length)
        cell_time(co2:gas4) = cell_volume(co2:gas4) &
            / E2Col(co2:gas4)%instr%tube_f
    elsewhere
        tube_time(co2:gas4) = 0d0
        cell_time(co2:gas4) = 0d0
    end where

    !> If user didn't set min and max time lags, does so by using tube properties for closed path
    !> and distances for open path
    do gas = co2, gas4
        if (E2Col(gas)%present) then
            if (TOSetup%min_lag(gas) < gui_tlag_threshold &
                .or. TOSetup%max_lag(gas) < gui_tlag_threshold) then
                if (E2Col(gas)%instr%path_type == 'closed') then
                    !> Closed path
                    nominal = tube_time(gas) + cell_time(gas)
                    E2Col(gas)%min_tl = max(0d0, nominal - 2d0)
                    E2Col(gas)%max_tl = min(nominal + mult(gas) * nominal, &
                        RPsetup%avrg_len * 60d0) + safety
                else
                    !> Open path
                    E2Col(gas)%min_tl = - dsqrt(E2Col(gas)%instr%hsep**2 &
                        + E2Col(gas)%instr%vsep**2) * 2d0 - safety
                    E2Col(gas)%max_tl = + dsqrt(E2Col(gas)%instr%hsep**2 &
                        + E2Col(gas)%instr%vsep**2) * 2d0 + safety
                end if
            else
                E2Col(gas)%min_tl = TOSetup%min_lag(gas)
                E2Col(gas)%max_tl = TOSetup%max_lag(gas)
            end if
        end if
    end do
end subroutine AdjustTimelagOptSettings
