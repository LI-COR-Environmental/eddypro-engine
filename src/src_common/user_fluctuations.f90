!***************************************************************************
! user_fluctuations.f90
! ---------------------
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
! \brief       Calculate fluctuations around a trend, defined by the \n
!              chosen detrending method
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UserFluctuations(Set, Primes, nrow, ncol, Tconst, LocStats, LocCol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: Tconst
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    type(UserStatsType), intent(in) :: LocStats
    type(ColType), intent(in) :: LocCol(ncol)
    real(kind = dbl), intent(out) :: Primes(nrow, ncol)
    !> local variables
    integer :: var


    select case(Meth%det(1:len_trim(Meth%det)))
        case('ld')
            call LinDetrend(Set, Primes, Tconst, LocCol, nrow, ncol)
        case('rm')
            call RunningMean(Set, Primes, Tconst, LocCol, nrow, ncol)
        case('ewa')
            call ExpWeightAvrg(Set, Primes, Tconst, LocCol, nrow, ncol)
        case('ba')
            do var = 1, ncol
                if (LocCol(var)%present) then
                    where(Set(1:nrow, var) /= error)
                        Primes(1:nrow, var) = Set(1:nrow, var) - LocStats%Mean(var)
                    elsewhere
                        Primes(1:nrow, var) = error
                    end where
                end if
            end do
    end select
end subroutine UserFluctuations
