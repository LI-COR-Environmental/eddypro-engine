!***************************************************************************
! m_numeric_kinds.f90
! -------------------
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
! \brief       Define numeric kinds
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_numeric_kinds
    implicit none
    save


    !> Declare parameters
    integer, parameter :: short = kind(2)
    !integer, parameter :: int   = kind(4)
    integer, parameter :: long  = kind(8)
    integer, parameter :: sgl   = kind(0.0)
    integer, parameter :: dbl   = kind(0.0d0)
    integer, parameter :: utf8 = selected_char_kind('ISO_10646')
end module m_numeric_kinds
