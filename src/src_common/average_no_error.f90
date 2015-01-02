!***************************************************************************
! average_no_error.f90
! --------------------
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
! \brief       Calculate column-wise averages on a 2d array \n
!              ignoring specified error values \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AverageNoError(Set, nrow, ncol, Mean, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Mean(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: Nact = 0
    real(kind = dbl) :: RawMean(ncol)


    RawMean = 0d0
    do j = 1, ncol
        Nact = 0
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                Nact = Nact + 1
                RawMean(j) = RawMean(j) + Set(i, j)
            end if
        end do
        if (Nact /= 0) then
            RawMean(j) = RawMean(j) / dble(Nact)
        else
            RawMean(j) = err_float
        end if
    end do
    Mean = 0.d0
    do j = 1, ncol
        if (RawMean(j) /= err_float) then
            Nact = 0
            do i = 1, nrow
                if (Set(i, j) /= err_float) then
                    Nact = Nact + 1
                    Mean(j) = Mean(j) + Set(i, j) - RawMean(j)
                end if
            end do
            if (Nact /= 0) then
                Mean(j) = Mean(j) / dble(Nact)
            else
                Mean(j) = err_float
            end if
        else
            Mean(j) = err_float
        end if
    end do

    where (Mean(:) /= err_float)
        Mean(:) = Mean(:) + RawMean(:)
    elsewhere
        Mean(:) = err_float
    end where
end subroutine AverageNoError
