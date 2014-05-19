!***************************************************************************
! replace_seonic_temperature.f90
! ------------------------------
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
! \brief       If sonic (or fast) temperature is out-ranged, tries to replace it
!              with any other available fast temperature reading
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReplaceSonicTemperature(Set, nrow, ncol, UserSet, unrow, uncol)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: unrow, uncol
    real(kind = dbl), intent(in) :: UserSet(unrow,uncol)
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Local variables
    integer :: i
    integer :: ord


    !> Look for a suitable temperature reading
    ord = nint(error)
    do i = 1, uncol
        if (UserCol(i)%var == 'ts' .and. &
            (UserStats%Mean(i) > 220d0 .and. UserStats%Mean(i) < 340d0)) ord = i
    end do

    do i = 1, uncol
        if (UserCol(i)%var == 'fast_t' .and. &
            (UserStats%Mean(i) > 220d0 .and. UserStats%Mean(i) < 340d0)) ord = i
    end do

    !> If a suitable temperature was found, replace the current one with the good one
    if (ord /= nint(error)) Set(1:nrow, ts) = UserSet(1:nrow, ord)
end subroutine ReplaceSonicTemperature
