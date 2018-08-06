!***************************************************************************
! fisher.f90
! ----------
! Copyright (C) 2018, LI-COR Biosciences
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
! \brief       Compute Fisher test on covariances with and without repeated 
!              values
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine fisher(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: var
    integer :: i
    integer :: j
    real(kind = dbl) :: Corr(ncol, ncol)
    real(kind = dbl) :: mSet(nrow, ncol)
    real(kind = dbl) :: mCorr(ncol, ncol)


    !> Compute mSet by eliminating repeated values in Set
    write(*, '(a)', advance = 'no') "  Evaluating correlation differences with and without repeated values.."
    mSet = Set
    do var = u, gas4
        do i = 2, nrow
            if ((Set(i, var) - Set(i-1, var)) < 1d-6) mSet(i, var) = error
        end do
    end do
    call CorrelationMatrixNoError(Set, size(Set, 1), size(Set, 2), Corr, error)
    call CorrelationMatrixNoError(mSet, size(mSet, 1), size(mSet, 2), mCorr, error)
    do i = u, gas4
        do j = u, gas4
            if (Corr(i, j) /= error .and. mCorr(i,j) /= error) then
                Essentials%CorrDiff(i, j) = dabs(Corr(i, j) - mCorr(i, j))
            else
                Essentials%CorrDiff(i, j) = error
            end if
        end do
    end do
    write(*, *) " Done."

end subroutine fisher