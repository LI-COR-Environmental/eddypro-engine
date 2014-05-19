!***************************************************************************
! tests_not_performed.f90
! -----------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
! \brief       Sets flags to 99999... for tests not performed
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestsNotPerformed()
    use m_rp_global_var
    implicit none

    if(.not.Test%sr) IntHF%sr = 99999999
    if(.not.Test%ar) IntHF%ar = 99999999
    if(.not.Test%do) IntHF%do = 99999999
    if(.not.Test%al) IntHF%al = 99999999
    if(.not.Test%sk) IntHF%sk = 99999999
    if(.not.Test%sk) IntSF%sk = 99999999
    if(.not.Test%ds) IntHF%ds = 99999999
    if(.not.Test%ds) IntSF%ds = 99999999
    if(.not.Test%tl) IntHF%tl = 9999
    if(.not.Test%tl) IntSF%tl = 9999
    if(.not.Test%aa) IntHF%aa = 9
    if(.not.Test%ns) IntHF%ns = 9
end subroutine TestsNotPerformed
