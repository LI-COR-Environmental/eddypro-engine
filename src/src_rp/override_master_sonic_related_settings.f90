!***************************************************************************
! override_master_sonic_related_settings.f90
! ------------------------------------------
! Copyright (C) 2016-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief       Override or adjust processing settings that are related to the
!              specific sonic anemometer model
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OverrideMasterSonicRelatedSettings()
    use m_rp_global_var
    implicit none


    !> If AoA selection is set to 'automatic', it's time to retrieve it
    if (RPsetup%calib_aoa == 'automatic') call InferAoaMethod()

    !> If running with EXP settings, override any previous setting and: (1) do not
    !> apply AoA correction; (2) apply the w-boost if MasterSonic is WM with
    !> appropriate firmware version
    if (EddyProProj%run_mode == 'express') then
        RPsetup%calib_aoa = 'none'
        RPsetup%calib_wboost = .true.
    end if

    !> If MasterSonic is other than a WM/WMPro with specific firmware versions,
    !> the w-boost cannot be applied
    if (.not. SonicDataHasWBug) then
        if (RPsetup%calib_wboost .and. EddyProProj%run_mode /= 'express') &
            call ExceptionHandler(95)
        RPsetup%calib_wboost = .false.
    end if

    !> If w-boost is enabled (either because selected and applicable or
    !> because running in express mode), do not apply AoA
    if (RPSetup%calib_wboost) RPSetup%calib_aoa = 'none'

    ! !> If Nakai and Shimoyama 2012 was selected or inferred as an AoA method,
    ! !> regardless of the sonic model, the w-boost correction is not to be applied.
    ! if (RPSetup%calib_aoa == 'nakai_12') RPsetup%calib_wboost = .false.

    !> Handle inappropriate AoA selections
    if (RPsetup%calib_aoa == 'nakai_06') then
        select case (MasterSonic%model(1:len_trim(MasterSonic%model) - 2))
            case ('r2','r3_50','r3_100','r3a_100')
                continue
            case default
                call ExceptionHandler(94)
                RPsetup%calib_aoa = 'none'
        end select
    end if
    if (RPsetup%calib_aoa == 'nakai_12') then
        select case (MasterSonic%model(1:len_trim(MasterSonic%model) - 2))
            case ('wm', 'wmpro')
                continue
            case default
                call ExceptionHandler(94)
                RPsetup%calib_aoa = 'none'
        end select
    end if

    !> Cross wind correction must be applied for R2
    if (MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'r2') &
        RPsetup%calib_cw = .true.

    !> Cross wind correction should not be applied for CSAT3
    if (MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'csat3') &
        RPsetup%calib_cw = .false.

    !> Cross wind correction should not be applied for CSAT3B
    if (MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'csat3b') &
        RPsetup%calib_cw = .false.

end subroutine OverrideMasterSonicRelatedSettings
