!***************************************************************************
! detect_master_sonic.f90
! -----------------------
! Copyright (C) 2015, LI-COR Biosciences
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
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DetectMasterSonic(LocCol, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    integer, intent(in) :: ncol
    !> local variables
    integer :: i

    MasterSonic = NullInstrument
    do i = 1, ncol
        !> u-component of wind vector
        if (trim(adjustl(LocCol(i)%var)) == 'u' &
            .and. LocCol(i)%Instr%master_sonic) then
            MasterSonic = LocCol(i)%Instr
            exit
        end if
    end do

    !> Convenient variable that tells if sonic data is biased by the w-boost bug
    !> Note: If SwVer is not available for the sonic, and the sonic is a WM/WMP, then
    !> assume the bug is present
    if ((MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'wm' .or. &
         MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'wmpro') .and. &
         (MasterSonic%sw_ver%major == 2329 .and. MasterSonic%sw_ver%minor < 700)) then
         SonicDataHasWBug = .true.
    else
         SonicDataHasWBug = .false.
    end if
end subroutine DetectMasterSonic
