!***************************************************************************
! test_higher_moments.f90
! -----------------------
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
! \brief       Checks for skewness and kurtosis of the whole file outside \n
!              realistic ranges and hard-flag file accordingly
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestHigherMoments(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: j = 0
    integer :: hflags(E2NumVar)
    integer :: sflags(E2NumVar)
    real(kind = dbl) :: LocSet(N, E2NumVar)
    real(kind = dbl) :: Primes(N, E2NumVar)
    real(kind = dbl) :: Skw(E2NumVar)
    real(kind = dbl) :: Kur(E2NumVar)

    write(*, '(a)', advance = 'no') '   Skewness & kurtosis test..'

    !> Initializations
    hflags = 9
    sflags = 9

    !> Define LocSet, limited to variables u to gas4
    LocSet(1:N, u:E2NumVar) = Set(1:N, u:E2NumVar)

    !> Linear detrending
    call LinDetrend(LocSet, Primes, RPsetup%avrg_len * 60, E2Col, N, E2NumVar)

    !> Hihgher moments 
    call KurtosisNoError(Primes, N, E2NumVar, Kur, error)
    call SkewnessNoError(Primes, N, E2NumVar, Skw, error)

    !> Hard/soft flags for skewness or kurtorsis out of bounds
    do j = u, GHGNumVar
        if (E2Col(j)%present) then
            if ((Skw(j) <= sk%hf_skmin) .or. (Skw(j) >= sk%hf_skmax) .or. &
                (Kur(j) <= sk%hf_kumin) .or. (Kur(j) >= sk%hf_kumax)) then
                hflags(j) = 1
            else
                hflags(j) = 0
            end if
            if ((Skw(j) <= sk%sf_skmin) .or. (Skw(j) >= sk%sf_skmax) .or. &
                (Kur(j) <= sk%sf_kumin) .or. (Kur(j) >= sk%sf_kumax)) then
                sflags(j) = 1
            else
                sflags(j) = 0
            end if
            Essentials%sk_s_skw(j) = Skw(j)
            Essentials%sk_s_kur(j) = Kur(j)
        else
            Essentials%sk_s_skw(j) = error
            Essentials%sk_s_kur(j) = error
        end if
    end do

    !>  Create 8-digits numbers containing hflag/sflag values
    IntHF%sk = 900000000
    IntSF%sk = 900000000
    do j = u, GHGNumVar
        IntHF%sk = IntHF%sk + hflags(j) * 10 **(GHGNumVar - j)
        IntSF%sk = IntSF%sk + sflags(j) * 10 **(GHGNumVar - j)
    end do
    write(*,'(a)') ' Done.'
end subroutine TestHigherMoments
