!***************************************************************************
! validate_filename_template.f90
! ------------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Weak check of compatibility of provided raw file prototype
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ValidateFilenameTemplate()
    use m_common_global_var
    implicit none
    !> local variables
    character(FilenameLen) :: Template

    !date patterns: yyyy, yy, ddd, dd, mm
    !time patterns: HH MM
    Template = trim(adjustl(EddyProProj%fname_template))

    !> Weak test
    if ( index(Template, 'yy') == 0 &
    .or. index(Template, 'dd') == 0 &
    .or. index(Template, 'HH') == 0 &
    .or. index(Template, 'MM') == 0) call ExceptionHandler(20)
end subroutine ValidateFilenameTemplate
