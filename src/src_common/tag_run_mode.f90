!***************************************************************************
! tag_run_mode.f90
! ----------------
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
! \brief       Add token to Timestamp_FilePadding to tag run mode
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TagRunMode()
    use m_common_global_var
    implicit none


    select case(EddyProProj%run_mode)
        case('express')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_exp'
        case('advanced')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_adv'
        case('md_retrieval')
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_mdr'
        case default
            Timestamp_FilePadding = trim(Timestamp_FilePadding) // '_nul'
    end select
end subroutine TagRunMode
