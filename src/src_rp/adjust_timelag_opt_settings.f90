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
    real(kind = dbl) :: nominal
    real(kind = dbl) :: pg_mult = 2d0    !< multiplier for Passive Gases
    real(kind = dbl) :: ag_mult = 10d0   !< multiplier for Active Gases
    real(kind = dbl) :: gui_tlag_threshold = -1000d0


    !> If user didn't set min and max time lags, does so by using tube properties for closed path
    !> and distances for open path
    !> CO2
    if (TOSetup%co2_min_lag < gui_tlag_threshold .or. TOSetup%co2_max_lag < gui_tlag_threshold) then
        if (E2Col(co2)%instr%path_type == 'closed') then
            nominal = (p * (E2Col(co2)%instr%tube_d /2d0)**2 * E2Col(co2)%instr%tube_l) &
                    / E2Col(co2)%instr%tube_f
            E2Col(co2)%min_tl = max(0d0, nominal - 2d0)
            E2Col(co2)%max_tl = min(nominal + pg_mult * nominal, RPsetup%avrg_len * 60d0)
        else
            E2Col(co2)%min_tl = - dsqrt(E2Col(co2)%instr%hsep**2 + E2Col(co2)%instr%vsep**2) * 2d0 - 0.3d0
            E2Col(co2)%max_tl =   dsqrt(E2Col(co2)%instr%hsep**2 + E2Col(co2)%instr%vsep**2) * 2d0 + 0.3d0
        end if
    else
        E2Col(co2)%min_tl = TOSetup%co2_min_lag
        E2Col(co2)%max_tl = TOSetup%co2_max_lag
    end if

    !> H2O
    if (TOSetup%h2o_min_lag < gui_tlag_threshold .or. TOSetup%h2o_max_lag < gui_tlag_threshold) then
        if (E2Col(h2o)%instr%path_type == 'closed') then
            nominal = (p * (E2Col(h2o)%instr%tube_d /2d0)**2 * E2Col(h2o)%instr%tube_l) &
                    / E2Col(h2o)%instr%tube_f
            E2Col(h2o)%min_tl = max(0d0, nominal - 2d0)
            E2Col(h2o)%max_tl = min(nominal + ag_mult * nominal, RPsetup%avrg_len * 60d0)
        else
            E2Col(h2o)%min_tl = - dsqrt(E2Col(h2o)%instr%hsep**2 + E2Col(h2o)%instr%vsep**2) * 2d0 - 0.3d0
            E2Col(h2o)%max_tl =   dsqrt(E2Col(h2o)%instr%hsep**2 + E2Col(h2o)%instr%vsep**2) * 2d0 + 0.3d0
        end if
    else
        E2Col(h2o)%min_tl = TOSetup%h2o_min_lag
        E2Col(h2o)%max_tl = TOSetup%h2o_max_lag
    end if

    !> CH4
    if (TOSetup%ch4_min_lag < gui_tlag_threshold .or. TOSetup%ch4_max_lag < gui_tlag_threshold) then
        if (E2Col(ch4)%instr%path_type == 'closed') then
            nominal = (p * (E2Col(ch4)%instr%tube_d /2d0)**2 * E2Col(ch4)%instr%tube_l) &
                    / E2Col(ch4)%instr%tube_f
            E2Col(ch4)%min_tl = max(0d0, nominal - 2d0)
            E2Col(ch4)%max_tl = min(nominal + pg_mult * nominal, RPsetup%avrg_len * 60d0)
        else
            E2Col(ch4)%min_tl = - dsqrt(E2Col(ch4)%instr%hsep**2 + E2Col(ch4)%instr%vsep**2) * 2d0 - 0.3d0
            E2Col(ch4)%max_tl =   dsqrt(E2Col(ch4)%instr%hsep**2 + E2Col(ch4)%instr%vsep**2) * 2d0 + 0.3d0
        end if
    else
        E2Col(ch4)%min_tl = TOSetup%ch4_min_lag
        E2Col(ch4)%max_tl = TOSetup%ch4_max_lag
    end if

    !> GAS4
    if (TOSetup%gas4_min_lag < gui_tlag_threshold .or. TOSetup%gas4_max_lag < gui_tlag_threshold) then
        if (E2Col(gas4)%instr%path_type == 'closed') then
            nominal = (p * (E2Col(gas4)%instr%tube_d /2d0)**2 * E2Col(gas4)%instr%tube_l) &
                    / E2Col(gas4)%instr%tube_f
            E2Col(gas4)%min_tl = max(0d0, nominal - 2d0)
            E2Col(gas4)%max_tl = min(nominal + pg_mult * nominal, RPsetup%avrg_len * 60d0)
        else
            E2Col(gas4)%min_tl = - dsqrt(E2Col(gas4)%instr%hsep**2 + E2Col(gas4)%instr%vsep**2) * 2d0 - 0.3d0
            E2Col(gas4)%max_tl =   dsqrt(E2Col(gas4)%instr%hsep**2 + E2Col(gas4)%instr%vsep**2) * 2d0 + 0.3d0
        end if
    else
        E2Col(gas4)%min_tl = TOSetup%gas4_min_lag
        E2Col(gas4)%max_tl = TOSetup%gas4_max_lag
    end if
end subroutine AdjustTimelagOptSettings
