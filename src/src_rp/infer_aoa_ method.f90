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
subroutine InferAoaMethod(SonicModel, SwVer)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: SonicModel
    type(SwVerType), intent(in) :: SwVer
    include '../src_common/interfaces_1.inc'


    if (EqualSwVer(SwVer, errSwVer)) then
        !> AoA selection based only on sonic model
        select case(SonicModel)
            case ('r3_50','r3_100', 'r2')
                RPsetup%calib_aoa = 'nakai_06'
            case ('wm','wmpro')
                RPsetup%calib_aoa = 'nakai_12'
            case default
                RPsetup%calib_aoa = 'none'
        end select

    else
        !> Insert here handling of AoA selection based
        !> on sonic's firmware version

    end if
end subroutine InferAoaMethod
