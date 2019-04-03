!***************************************************************************
! interpret_diagnostics.f90
! -------------------------
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
! \brief       Interpret binary diagnostic flags and accumulate info
!              on bad flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InterpretLicorDiagnostics(DiagSet, nrow, ncol)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: DiagSet(nrow, ncol)
    !> local variables
    integer :: i
    integer :: n72
    integer :: n75
    integer :: n77
    include '../src_common/interfaces.inc'


    n72 = 0
    n75 = 0
    n77 = 0
    Diag7200 = diag7200type(Diag7200%present, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0d0)
    Diag7500 = diag7500type(Diag7500%present, 0, 0, 0, 0, 0d0)
    Diag7700 = diag7700type(Diag7700%present, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    do i = 1, nrow
        !> LI-7200
        if (DiagSet(i, diag72) /= error) then
            n72 = n72 + 1
            Diag7200%head_detect = Diag7200%head_detect + ibits(nint(DiagSet(i, diag72)),  12, 1)
            Diag7200%t_out       = Diag7200%t_out       + ibits(nint(DiagSet(i, diag72)),  11, 1)
            Diag7200%t_in        = Diag7200%t_in        + ibits(nint(DiagSet(i, diag72)),  10, 1)
            Diag7200%aux_in      = Diag7200%aux_in      + ibits(nint(DiagSet(i, diag72)),  9,  1)
            Diag7200%delta_p     = Diag7200%delta_p     + ibits(nint(DiagSet(i, diag72)),  8,  1)
            Diag7200%chopper     = Diag7200%chopper     + ibits(nint(DiagSet(i, diag72)),  7,  1)
            Diag7200%detector    = Diag7200%detector    + ibits(nint(DiagSet(i, diag72)),  6,  1)
            Diag7200%pll         = Diag7200%pll         + ibits(nint(DiagSet(i, diag72)),  5,  1)
            Diag7200%sync        = Diag7200%sync        + ibits(nint(DiagSet(i, diag72)),  4,  1)
        end if
        !> LI-7500
        if (DiagSet(i, diag75) /= error) then
            n75 = n75 + 1
            Diag7500%chopper     = Diag7500%chopper     + ibits(nint(DiagSet(i, diag75)), 7, 1)
            Diag7500%detector    = Diag7500%detector    + ibits(nint(DiagSet(i, diag75)), 6, 1)
            Diag7500%pll         = Diag7500%pll         + ibits(nint(DiagSet(i, diag75)), 5, 1)
            Diag7500%sync        = Diag7500%sync        + ibits(nint(DiagSet(i, diag75)), 4, 1)
        end if
        !> LI-7700
        if (DiagSet(i, diag77) /= error) then
            n77 = n77 + 1
            Diag7700%not_ready              = Diag7700%not_ready              + &
                ibits(nint(DiagSet(i, diag77)),  15, 1)
            Diag7700%no_signal              = Diag7700%no_signal              + &
                ibits(nint(DiagSet(i, diag77)),  14, 1)
            Diag7700%re_unlocked            = Diag7700%re_unlocked            + &
                ibits(nint(DiagSet(i, diag77)),  13, 1)
            Diag7700%bad_temp               = Diag7700%bad_temp               + &
                ibits(nint(DiagSet(i, diag77)),  12, 1)
            Diag7700%laser_temp_unregulated = Diag7700%laser_temp_unregulated + &
                ibits(nint(DiagSet(i, diag77)),  11, 1)
            Diag7700%block_temp_unregulated = Diag7700%block_temp_unregulated + &
                ibits(nint(DiagSet(i, diag77)),  10, 1)
            Diag7700%motor_spinning         = Diag7700%motor_spinning         + &
                 ibits(nint(DiagSet(i, diag77)),   9, 1)
            Diag7700%pump_on                = Diag7700%pump_on                + &
                 ibits(nint(DiagSet(i, diag77)),   8, 1)
            Diag7700%top_heater_on          = Diag7700%top_heater_on          + &
                 ibits(nint(DiagSet(i, diag77)),   7, 1)
            Diag7700%bottom_heater_on       = Diag7700%bottom_heater_on       + &
                 ibits(nint(DiagSet(i, diag77)),   6, 1)
            Diag7700%calibrating            = Diag7700%calibrating            + &
                 ibits(nint(DiagSet(i, diag77)),   5, 1)
            Diag7700%motor_failure          = Diag7700%motor_failure          + &
                 ibits(nint(DiagSet(i, diag77)),   4, 1)
            Diag7700%bad_aux_tc1            = Diag7700%bad_aux_tc1            + &
                 ibits(nint(DiagSet(i, diag77)),   3, 1)
            Diag7700%bad_aux_tc2            = Diag7700%bad_aux_tc2            + &
                 ibits(nint(DiagSet(i, diag77)),   2, 1)
            Diag7700%bad_aux_tc3            = Diag7700%bad_aux_tc3            + &
                 ibits(nint(DiagSet(i, diag77)),   1, 1)
            Diag7700%box_connected          = Diag7700%box_connected          + &
                 ibits(nint(DiagSet(i, diag77)),   0, 1)
        end if
    end do

    !> AGC or RSSI (depending on sw version) stored anyway in variable "AGC".
    if (.not. CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('5.3.0'))) then
        do i = 1, nrow
            !> LI-7200
            if (DiagSet(i, diag72) /= error) Diag7200%AGC = Diag7200%AGC + &
                    (ibits(nint(DiagSet(i, diag72)), 0, 4) * 6.25d0) + 6.25d0
            !> LI-7500
            if (DiagSet(i, diag75) /= error) Diag7500%AGC = Diag7500%AGC + &
                (ibits(nint(DiagSet(i, diag75)), 0, 4) * 6.25d0) + 6.25d0
        end do
    else
        do i = 1, nrow
            !> LI-7200
            if (DiagSet(i, diag72) /= error) Diag7200%AGC = Diag7200%AGC + &
                ibits(nint(DiagSet(i, diag72)), 0, 4) * 6.6667d0
            !> LI-7500
            if (DiagSet(i, diag75) /= error) Diag7500%AGC = Diag7500%AGC + &
                ibits(nint(DiagSet(i, diag75)), 0, 4) * 6.6667d0
        end do
    end if

    !> Adjusts values that report "1" for "good" and averages AGC
    if (n72 > 0) then
        Diag7200%head_detect = n72 - Diag7200%head_detect
        Diag7200%t_out       = n72 - Diag7200%t_out
        Diag7200%t_in        = n72 - Diag7200%t_in
        Diag7200%aux_in      = n72 - Diag7200%aux_in
        Diag7200%delta_p     = n72 - Diag7200%delta_p
        Diag7200%chopper     = n72 - Diag7200%chopper
        Diag7200%detector    = n72 - Diag7200%detector
        Diag7200%pll         = n72 - Diag7200%pll
        Diag7200%sync        = n72 - Diag7200%sync
        Diag7200%AGC = Diag7200%AGC / dfloat(n72)
    end if

    if (n75 > 0) then
        Diag7500%chopper     = n75 - Diag7500%chopper
        Diag7500%detector    = n75 - Diag7500%detector
        Diag7500%pll         = n75 - Diag7500%pll
        Diag7500%sync        = n75 - Diag7500%sync
        Diag7500%AGC = Diag7500%AGC / dfloat(n75)
    end if
end subroutine InterpretLicorDiagnostics
