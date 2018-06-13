!***************************************************************************
! kid.f90
! -------
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
! \brief       Compute Kurtosis index on differenced variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine KID(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: var
    real(kind = dbl) :: Primes(nrow, ncol)

    do var = u, gas4
        call VariableStochasticDetrending(Set(:, var), Primes(:, var), nrow)
        call KurtosisNoError(Primes(:, var), nrow, 1, Essentials%KID(var), error)
        Essentials%ZCD(var) = count(abs(Primes(:, var)) < 1d-6)
    end do
end subroutine KID