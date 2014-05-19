!***************************************************************************
! power_of_two.f90
! ----------------
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
! \brief       Calculates closest power of two to N
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PowerOfTwo(N, n2)
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(out) :: n2
    !> local variables
    integer :: i = 0

    loop: do i = 1, 25
        if((2**i) > N) then
            n2 = 2**(i - 1)
            exit loop
        end if
    end do loop
end subroutine PowerOfTwo
