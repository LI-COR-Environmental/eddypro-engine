!***************************************************************************
! set_licor_diagnostics.f90
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
! \brief       Set AGC and RSSI as available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SetLicorDiagnostics(M)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: M
    !> local variables
    integer :: i

    Essentials%AGC72 = error
    Essentials%AGC75 = error
    Essentials%RSSI77 = error

    !> First, uses user-columns labeled as either AGC or RSSI. Note that AGC72(75) holds either
    !> AGC or RSSI from LI-7200(LI-7500), because no instrument will provide both
    do i = 1, M
        !> LI-7200
        if (UserCol(i)%var == 'AGC' .and. index(UserCol(i)%instr%model, 'li7200') /= 0) then
            Essentials%AGC72 = UserStats%Mean(i)
            exit
        end if
        if (UserCol(i)%var == 'RSSI' .and. index(UserCol(i)%instr%model, 'li7200') /= 0) then
            Essentials%AGC72 = UserStats%Mean(i)
            exit
        end if
        !> LI-7500
        if (UserCol(i)%var == 'AGC' .and. index(UserCol(i)%instr%model, 'li7500') /= 0) then
            Essentials%AGC75 = UserStats%Mean(i)
            exit
        end if
        if (UserCol(i)%var == 'RSSI' .and. index(UserCol(i)%instr%model, 'li7500') /= 0) then
            Essentials%AGC75 = UserStats%Mean(i)
            exit
        end if
        !> LI-7700
        if (UserCol(i)%var == 'RSSI' .and. index(UserCol(i)%instr%model, 'li7700') /= 0) then
            Essentials%RSSI77 = UserStats%Mean(i)
            exit
        end if
    end do
    !> Then, if AGC from LI-COR's flags are available, override any previous value with it
    !> This is because the user may call AGC a column that is not the actual AGC, while the
    !> value from LI-COR's flags are surely correct
    if (Diag7200%AGC  /= 0d0) Essentials%AGC72 = Diag7200%AGC
    if (Diag7500%AGC /= 0d0) Essentials%AGC75 = Diag7500%AGC
end subroutine SetLicorDiagnostics
