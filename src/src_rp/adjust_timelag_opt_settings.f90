!***************************************************************************
! adjust_timelag_opt_settings.f90
! -------------------------------
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!***************************************************************************
!
! \brief       Adjust time-lag opt settings if user did not set or set improperly
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

    !> If user didn't set min and max time-lags, does so by using tube properties for closed path
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
