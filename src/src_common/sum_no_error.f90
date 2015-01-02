!***************************************************************************
! sum_no_error.f90
! --------------------
! Copyright (C) 2015, LI-COR Biosciences
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
! \brief       Calculate column-wise sums on a 2d array \n
!              ignoring specified error values \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SumNoError(Set, nrow, ncol, Summ, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Summ(ncol)
    !> local variables
    integer :: i
    integer :: j
    logical :: data_exist

    Summ = 0d0
    do j = 1, ncol
        data_exist = .false.
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                data_exist = .true.
                Summ(j) = Summ(j) + Set(i, j)
            end if
        end do
        if (.not. data_exist) Summ(j) = err_float
    end do
end subroutine SumNoError

