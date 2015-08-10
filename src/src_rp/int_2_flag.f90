!***************************************************************************
! int_2_flag.f90
! --------------
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
! \brief       Converts integer flags into characher flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Int2Flags(len)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: len

    call int2char(IntHF%sr, CharHF%sr, len)
    call int2char(IntHF%ar, CharHF%ar, len)
    call int2char(IntHF%do, CharHF%do, len)
    call int2char(IntHF%al, CharHF%al, len)
    call int2char(IntHF%sk, CharHF%sk, len)
    call int2char(IntSF%sk, CharSF%sk, len)
    call int2char(IntHF%ds, CharHF%ds, len)
    call int2char(IntSF%ds, CharSF%ds, len)
    call int2char(IntHF%tl, CharHF%tl, len)
    call int2char(IntSF%tl, CharSF%tl, len)
    call int2char(IntHF%aa, CharHF%aa, len)
    call int2char(IntHF%ns, CharHF%ns, len)
end subroutine Int2Flags
