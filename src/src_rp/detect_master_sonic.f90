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
            return
        end if
    end do

end subroutine DetectMasterSonic
