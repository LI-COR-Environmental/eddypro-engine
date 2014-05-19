!***************************************************************************
! check_file_prototype.f90
! ------------------------
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
! \brief       Weak check of compatibility of provided raw file prototype
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CheckFilePrototype()
    use m_common_global_var
    implicit none
    !> local variables
    character(64) :: Pattern

    !date patterns: yyyy, yy, ddd, dd, mm
    !time patterns: HH MM
    Pattern = EddyProProj%fproto(1:len_trim(EddyProProj%fproto))

    !> Weak test
    if ( index(Pattern, 'yy') == 0 &
    .or. index(Pattern, 'dd') == 0 &
    .or. index(Pattern, 'HH') == 0 &
    .or. index(Pattern, 'MM') == 0) call ErrorHandle(0, 0, 20)
end subroutine CheckFilePrototype
