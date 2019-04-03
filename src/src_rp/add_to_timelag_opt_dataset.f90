!***************************************************************************
! add_to_timelag_opt_dataset.f90
! ------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Store calculated time-lags and other variables used for the
!              time-lag optimization, if all conditions are met
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddToTimelagOptDataset(TimelagOpt, nrow, n)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: n
    type(TimeLagOptType), intent(inout):: TimelagOpt(nrow)


    !> Passive gases
    if (E2Col(co2)%present &
        .and. dabs(Flux0%co2) > TOSetup%co2_min_flux &
        .and. Essentials%used_timelag(co2) /= E2Col(co2)%max_tl &
        .and. Essentials%used_timelag(co2) /= E2Col(co2)%min_tl) then
            TimelagOpt(n)%tlag(co2) = Essentials%used_timelag(co2)
    else
        TimelagOpt(n)%tlag(co2) = error
    end if

    if (E2Col(ch4)%present &
        .and. Flux0%ch4 > TOSetup%ch4_min_flux &
        .and. Essentials%used_timelag(ch4) /= E2Col(ch4)%max_tl &
        .and. Essentials%used_timelag(ch4) /= E2Col(ch4)%min_tl) then
        TimelagOpt(n)%tlag(ch4) = Essentials%used_timelag(ch4)
    else
        TimelagOpt(n)%tlag(ch4) = error
    end if

    if (E2Col(gas4)%present &
        .and. Flux0%gas4 > TOSetup%gas4_min_flux &
        .and. Essentials%used_timelag(gas4) /= E2Col(gas4)%max_tl &
        .and. Essentials%used_timelag(gas4) /= E2Col(gas4)%min_tl) then
        TimelagOpt(n)%tlag(gas4) = Essentials%used_timelag(gas4)
    else
        TimelagOpt(n)%tlag(gas4) = error
    end if

    !> Water vapor and RH
    if (E2Col(h2o)%present) then
        if (Flux0%LE > TOSetup%le_min_flux &
            .and. Essentials%used_timelag(h2o) /= E2Col(h2o)%max_tl &
            .and. Essentials%used_timelag(h2o) /= E2Col(h2o)%min_tl) then
            TimelagOpt(n)%tlag(h2o) = Essentials%used_timelag(h2o)
        else
            TimelagOpt(n)%tlag(h2o) = error
        end if
        if (Stats%RH >= 0d0 .and. Stats%RH <= 100d0) then
            TimelagOpt(n)%RH = Stats%RH
        else
            TimelagOpt(n)%RH = error
        end if
    else
        TimelagOpt(n)%tlag(h2o) = error
        TimelagOpt(n)%RH = error
    end if
end subroutine AddToTimelagOptDataset
