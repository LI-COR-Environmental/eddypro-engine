!***************************************************************************
! set_licor_diagnostics.f90
! -------------------------
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
