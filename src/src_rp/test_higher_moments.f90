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
    integer :: i = 0
    integer :: j = 0
    integer :: hflags(GHGNumVar)
    integer :: sflags(GHGNumVar)
    real(kind = dbl) :: LocSet(N, GHGNumVar)
    real(kind = dbl) :: Trend(N, GHGNumVar)
    real(kind = dbl) :: Primes(N, GHGNumVar)
    real(kind = dbl) :: Mean(GHGNumVar)
    real(kind = dbl) :: Var(GHGNumVar)
    real(kind = dbl) :: Skw(GHGNumVar)
    real(kind = dbl) :: Kur(GHGNumVar)
    real(kind = dbl) :: sumx1(GHGNumVar)
    real(kind = dbl) :: sumx2(GHGNumVar)
    real(kind = dbl) :: sumtime
    real(kind = dbl) :: sumtime2
    real(kind = dbl) :: b(GHGNumVar)

    write(*, '(a)', advance = 'no') '   Skewness & kurtosis test..'

    !> Initializations
    hflags = 9
    sflags = 9

    !> Define LocSet, limited to variables u to gas4
    LocSet(1:N, u:GHGNumVar) = Set(1:N, u:GHGNumVar)

    !> Linear detrending
    !> mean values
    Mean(:) = sum(LocSet(1:N, :), dim = 1)
    Mean(:) = Mean(:) / dfloat(N)
    sumx1 = 0.d0
    sumx2 = 0.d0
    sumtime = 0.d0
    sumtime2 = 0.d0
    do i = 1, N
        sumx1(:) = sumx1(:) + (LocSet(i, :)*(dble(i - 1)))
        sumx2(:) = sumx2(:) + LocSet(i, :)
        sumtime = sumtime + (dble(i - 1))
        sumtime2 = sumtime2 + (dble(i - 1))**2
    end do
    b(:) = (sumx1(:) - (sumx2(:)*sumtime) / dble(N)) / &
          (sumtime2 - (sumtime*sumtime) / dble(N))
    !> Trend
    do i = 1, N
        Trend(i, :) = Mean(:) + b(:) * (dble(i - 1) - sumtime / dble(N))
    end do
    !> Fluctuations
    do i = 1, N
        Primes(i, :) = LocSet(i, :) - Trend(i, :)
    end do

    !> Standard deviations
    Var = 0.d0
    do i = 1, N
        Var(:) = Var(:) + Primes(i, :) **2
    end do
    Var = Var / dble(N - 1)

    !> Skewness & Kurtosis
    Skw = 0.d0
    Kur = 0.d0
    do i = 1, N
        Skw(:) = Skw(:) + (Primes(i, :) / dsqrt(Var(:))) **3
        Kur(:) = Kur(:) + (Primes(i, :) / dsqrt(Var(:))) **4
    end do
    Skw(:) = Skw(:) / dble(N)
    Kur(:) = Kur(:) / dble(N)

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
