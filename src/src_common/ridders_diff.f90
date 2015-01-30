!***************************************************************************
! ridders_diff.f90
! ----------------
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
! \brief       Returns the derivative of a function f at a point x by Ridder's implementation
!              of Neville's ploynomial extrapolation to zero, applied to the central difference
!              approximation of the derivative: f' = (f(x + dh) - f(x - dh)) / (2 * dh).
! \author      Gerardo Fratini, Antonio Forgione
! \notes       Ridder's algorithm works by evaluating derivatives with smaller and smaller step sizes,
!              fitting a polynomial and extrapolating its value to zero. This algorithm is generally
!              more accurate and much more reliable, but much slower than using simple finite differences.
!              It is robust to awkward functions and does not require careful tuning of the step size,
!              furthermore it provides an estimate of the errors.
!              http://www.edwardrosten.com/cvd/toon/html-user/group__gFunctions.html
!              https://github.com/edrosten/TooN/blob/master/functions/derivatives.h
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function ridders_diff(f, x_var, p1, p2, p3, partial_var_index, init_step_size, abs_err)
    use m_common_global_var
    implicit none
    !> In/out variables
    real(kind = dbl), external :: f
    real(kind = dbl), intent(in) :: x_var
    real(kind = dbl), intent(in) :: p1, p2, p3
    integer, intent(in) :: partial_var_index  !> index of the variable about which to differentiate between x_var, p1, p2 and p3
    real(kind = dbl), intent(in) ::  init_step_size   !> initial step size, input
    real(kind = dbl), intent(out) :: abs_err
    !> Local variables
    integer :: i
    integer :: j
    integer, parameter :: MAX_ITER = 10   !> 400 to tune for accuracy instead of speed (see above references)
    real(kind = dbl) :: current_err
    real(kind = dbl) :: err_1
    real(kind = dbl) :: err_2
    real(kind = dbl) :: factor
    real(kind = dbl) :: dh
    real(kind = dbl) :: table(MAX_ITER, MAX_ITER) ! Naville's table
    real(kind = dbl), parameter :: CONV_FACTOR = 1.4d0
    real(kind = dbl), parameter :: SQUARED_CONV_FACTOR = CONV_FACTOR ** 2.
    real(kind = dbl), parameter :: ERR_GUARD = 2d0  ! stop when abs_err is greater than ERR_GUARD times the previous iteration


    ridders_diff = error
    abs_err = 1d50
    dh = init_step_size

    select case (partial_var_index)
        case (1)
            table(1, 1) = (f(x_var + dh, p1, p2, p3) - f(x_var - dh, p1, p2, p3)) / (2.d0 * dh)
        case (2)
            table(1, 1) = (f(x_var, p1 + dh, p2, p3) - f(x_var, p1 - dh, p2, p3)) / (2.d0 * dh)
        case (3)
            table(1, 1) = (f(x_var, p1, p2 + dh, p3) - f(x_var, p1, p2 - dh, p3)) / (2.d0 * dh)
        case (4)
            table(1, 1) = (f(x_var, p1, p2, p3 + dh) - f(x_var, p1, p2, p3 - dh)) / (2.d0 * dh)
    end select

    do i = 2, MAX_ITER
        dh = dh / CONV_FACTOR
        factor = SQUARED_CONV_FACTOR

        select case (partial_var_index)
            case (1)
                table(1, i) = (f(x_var + dh, p1, p2, p3) - f(x_var - dh, p1, p2, p3)) / (2.d0 * dh)
            case (2)
                table(1, i) = (f(x_var, p1 + dh, p2, p3) - f(x_var, p1 - dh, p2, p3)) / (2.d0 * dh)
            case (3)
                table(1, i) = (f(x_var, p1, p2 + dh, p3) - f(x_var, p1, p2 - dh, p3)) / (2.d0 * dh)
            case (4)
                table(1, i) = (f(x_var, p1, p2, p3 + dh) - f(x_var, p1, p2, p3 - dh)) / (2.d0 * dh)
        end select

        do j = 2, i
            table(j, i) = (table(j - 1, i) * factor - table(j - 1, i - 1)) / (factor - 1.)
            factor = SQUARED_CONV_FACTOR * factor
            err_1 = dabs(table(j, i) - table(j - 1, i)) ! difference between the current point (high order, i-th column of the
                                                        ! Naville's table) and the corresponding lower order point at the current step size
            err_2 = dabs(table(j, i) - table(j - 1, i - 1)) ! difference between two consecutive points at the
                                                            ! corresponding lower order point at the larger stepsize
            current_err = max(err_1, err_2)
            if (current_err <= abs_err) then
                abs_err = current_err
                ridders_diff = table(j, i)
            endif
        enddo
        if (dabs(table(i, i) - table(i - 1, i - 1)) >= ERR_GUARD * abs_err) return
    enddo
end function ridders_diff
!
!**********************************************************************************************
!**********************************************************************************************
!
double precision function func(indvar, p1, p2, p3)
    use m_common_global_var
    implicit none
    real(kind = dbl), intent(in) :: p1
    real(kind = dbl), intent(in) :: p2
    real(kind = dbl), intent(in) :: p3
    real(kind = dbl), intent(in) :: indvar

    if (p1 == error .or. p2 == error .or. p3 == error) then
        func = error
    else
        func = dlog(p1 * indvar/p2 / (1.d0 + (indvar/p2)**(2.d0*p3))**(1.1666667d0/p3))
    end if
end function func
