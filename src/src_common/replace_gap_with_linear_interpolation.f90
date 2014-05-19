!***************************************************************************
! replace_gap_with_linear_interpolation.f90
! -----------------------------------------
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
! \brief       replace gaps with linear interpolation of neighbouring data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReplaceGapWithLinearInterpolation(vec, N, N2, err_code)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: err_code
    real(kind = dbl), intent(inout) :: vec(N)
    integer, intent(inout) :: N2
    !> local variables
    integer :: cnt
    integer :: i
    integer :: ii
    integer :: k
    integer :: NN
    real(kind = dbl) :: slope
    real(kind = dbl) :: q


    N2 = N
    !> First, eliminate leading and trailing error codes
    do i = 1, N
        if (vec(1) == err_code) then
            vec(1: N2-1) = vec(2: N2)
            N2 = N2 - 1
            cycle
        end if
        exit
    end do
    NN = N2

    do i = NN, 1, -1
        if (vec(i) == err_code) then
            N2 = N2 - 1
            cycle
        end if
        exit
    end do

    !> Now eliminate internal ones
    i = 0
    do
        i = i + 1
        if (i >= N2) exit
        !> If an error code is found, look for consecutive ones
        cnt = 0
        if (vec(i) == err_code) then
            cnt = cnt + 1
            il: do ii = i + 1, N2
                if (vec(ii) == err_code) then
                    cnt = cnt + 1
                else
                    exit il
                end if
            end do il
            !> Now a sequence of error codes as long as "cnt" has been
            !> detected, and can be replaced with linear interpolation
            slope = (vec(i + cnt) - vec(i - 1)) / (dble(cnt + 1))
            q = vec(i - 1)
            do k = i, i + cnt - 1
                vec(k) = (slope * (dble(k - (i - 1))) + q)
            end do
            i = i + cnt - 1
        end if
    end do
end subroutine ReplaceGapWithLinearInterpolation
