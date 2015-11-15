!***************************************************************************
! infer_aoa_method.f90
! --------------------
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
! \brief       Infer most appropriate AoA method based on , \n
!              sonic anemometer model and serial number
! \author      Gerardo Fratini
! \notes
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InferAoaMethod(mSonic)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(InstrumentType), intent(inout) :: mSonic
    include '../src_common/interfaces_1.inc'
    character(len(mSonic%sw_ver_string)), external :: replace
    type(SwVerType) :: SwVer

    !> AoA selection based only on sonic model
    select case(mSonic%model)
        case ('r3_50','r3_100', 'r2')
            RPsetup%calib_aoa = 'nakai_06'
        case ('wm','wmpro')
            !>Replace dash with dot if the case
            mSonic%sw_ver_string = replace(mSonic%sw_ver_string, &
                '-', '.', len(mSonic%sw_ver_string))
            !>Build SwVer object from string
            SwVer = SwVerFromString(mSonic%sw_ver_string)
            if (SwVer%minor < 700) then
                RPsetup%calib_aoa = 'nakai_12'
            else
                RPsetup%calib_aoa = 'none'
            end if
         case default
            RPsetup%calib_aoa = 'none'
    end select
end subroutine InferAoaMethod
