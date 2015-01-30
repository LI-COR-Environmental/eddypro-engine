!***************************************************************************
! fluctuations.f90
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
subroutine Fluctuations(Set, Primes, nrow, ncol, Tconst, LocStats, LocCol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: Tconst
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    type(StatsType), intent(in) :: LocStats
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
                else
                    Primes(1:nrow, var) = error
                end if
            end do
    end select
end subroutine Fluctuations

!***************************************************************************
!
! \brief       Linear de-trending
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LinDetrend(Set, Primes, Tconst, LocCol, N, M)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: Tconst
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    type(ColType), intent(in) :: LocCol(M)
    real(kind = dbl), intent(out) :: Primes(N, M)
    !> local variables
    integer :: i
    integer :: chuncks
    integer :: j
    integer :: nn(M)
    integer :: mm(M)
    integer :: Rconst
    integer :: imin
    integer :: imax
    real(kind = dbl) :: sumx1(M)
    real(kind = dbl) :: sumx2(M)
    real(kind = dbl) :: mean(M)
    real(kind = dbl) :: sumtime(M)
    real(kind = dbl) :: sumtime2(M)
    real(kind = dbl) :: b(M)
    real(kind = dbl) :: Trend(N, M)


    !> Time constant in terms of rows
    Rconst = Tconst * idint(Metadata%ac_freq)
    chuncks = 0
    do
        if (N < Rconst) then
            imin = 1
            imax = N
        else
            chuncks = chuncks + 1
            imin = 1 + (chuncks-1) * Rconst
            if (N - chuncks * Rconst >= Rconst) then
                imax = chuncks * Rconst
            else
                imax = N
            end if
        end if

        !> Linear regression
        sumx1 = 0d0
        sumx2 = 0d0
        sumtime = 0d0
        sumtime2 = 0d0
        nn = 0
        do j = 1, M
            if(LocCol(j)%present) then
                do i = imin, imax
                    if (Set(i, j) /= error) then
                        nn(j) = nn(j) + 1
                        sumx1(j) = sumx1(j) + (Set(i, j)*(dble(nn(j) - 1)))
                        sumx2(j) = sumx2(j) + Set(i, j)
                        sumtime(j) = sumtime(j) + (dble(nn(j) - 1))
                        sumtime2(j) = sumtime2(j) + (dble(nn(j) - 1))**2
                    end if
                end do
                where (nn(:) /= 0)
                    mean(:) = sumx2(:) / dble(nn(:))
                end where
            end if
        end do

        !> Trend
        mm = 0
        do j = 1, M
            if (LocCol(j)%present) then
                b(j) = (sumx1(j) - (sumx2(j)*sumtime(j)) / dble(nn(j))) / &
                      (sumtime2(j) - (sumtime(j)*sumtime(j)) / dble(nn(j)))
                do i = imin, imax
                    mm(j) = mm(j) + 1
                    if(Set(i, j) /= error) then
                        Trend(i, j) = mean(j) + b(j) * &
                            (dble(mm(j) - 1) - sumtime(j) / dble(nn(j)))
                    else
                        Trend(i, j) = error
                    end if
                end do
            end if
        end do

        !> Fluctuations
        do j = 1, M
            if (LocCol(j)%present) then
                do i = imin, imax
                    if(Set(i, j) /= error .and. Trend(i, j) /= error) then
                        Primes(i, j) = Set(i, j) - Trend(i, j)
                    else
                        Primes(i, j) = error
                    end if
                end do
            end if
        end do
        if (imax == N) exit
    end do
end subroutine LinDetrend

!***************************************************************************
!
! \brief       Running mean
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RunningMean(Set, Primes, Tconst, LocCol, N, M)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: Tconst
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    type(ColType), intent(in) :: LocCol(M)
    real(kind = dbl), intent(out) :: Primes(N, M)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: back = 0
    integer :: cnt = 0
    integer :: Rconst
    real(kind = dbl) :: Trend(N, M)


    !> Time constant in terms of rows
    Rconst = Tconst * idint(Metadata%ac_freq)

    !> Trend
    do k = 1, M
        if (LocCol(k)%present) then
            !> First Rconst elements in the time series have a special treatment
            do i = 1, Rconst
                Trend(i, k) = 0.d0
                cnt = 0
                do j = 1, i
                    if (Set(j, k) /= error) then
                        cnt = cnt + 1
                        Trend(i, k) = Trend(i, k) + Set(j, k)
                    end if
                end do
                if (cnt /= 0) then
                    Trend(i, k) = Trend(i, k) / dble(cnt)
                else
                    Trend(i, k) = error
                end if
            end do
            !> Then, normal running mean
            do i = Rconst + 1, N
!                if (Set(i, k) /= error .and. Set(i - Rconst, k) /= error) then
!                    Trend(i, k) = Trend(i - 1, k) + &
!                        (Set(i, k) - Set(i - RConst, k)) / dble(Rconst)
!                else
!                    Trend(i, k) = Trend(i - 1, k)
!                end if

                if (Set(i, k) /= error .and. Set(i - Rconst, k) /= error) then
                    do back = 1, Rconst
                        if (Trend(i - back, k) /= error) then
                            Trend(i, k) = Trend(i - back, k) + &
                                (Set(i, k) - Set(i - RConst, k)) / dble(Rconst)
                            exit
                        end if
                    end do
                else
                    Trend(i, k) = Trend(i - 1, k)
                end if
            end do
            !> Fluctuations
            where (Trend(1:N, k) /= error .and. Set(1:N, k) /= error)
                Primes(1:N, k) = Set(1:N, k) - Trend(1:N, k)
            elsewhere
                Primes(1:N, k) = error
            end where
        else
            Primes(:, k) = error
        end if
    end do
end subroutine RunningMean

!***************************************************************************
!
! \brief       Exponentially-weighted running average
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExpWeightAvrg(Set, Primes, Tconst, LocCol, N, M)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: Tconst
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    type(ColType), intent(in) :: LocCol(M)
    real(kind = dbl), intent(out) :: Primes(N, M)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: back = 0
    integer :: cnt = 0
    integer :: Rconst
    real(kind = dbl) :: b
    real(kind = dbl) :: dt
    real(kind = dbl) :: Trend(N, M)


    !> Time consant in terms of rows
    Rconst = Tconst * idint(Metadata%ac_freq)

    !> Weighing factor b
    dt = 1.d0 / Metadata%ac_freq
    b = dexp(-dt / dfloat(Tconst))

    !> Trend
    do k = 1, M
        if (LocCol(k)%present) then
            !> First Rconst elements in the time series have a special treatment
            do i = 1, Rconst
                Trend(i, k) = 0.d0
                cnt = 0
                do j = 1, i
                    if (Set(j, k) /= error) then
                        cnt = cnt + 1
                        Trend(i, k) = Trend(i, k) + Set(j, k)
                    end if
                end do
                if (cnt /= 0) then
                    Trend(i, k) = Trend(i, k) / dble(cnt)
                else
                    Trend(i, k) = error
                end if
            end do
            !> Then, normal running mean
            do i = Rconst + 1, N
                if (Set(i, k) /= error) then
                    do back = 1, Rconst
                        if (Trend(i - back, k) /= error) then
                            Trend(i, :) = b * Trend(i - back, :) + (1.d0 - b) * Set(i, :)
                            exit
                        end if
                    end do
                else
                    Trend(i, k) = Trend(i - 1, k)
                end if
            end do
            !> Fluctuations
            where (Trend(1:N, k) /= error .and. Set(1:N, k) /= error)
                Primes(1:N, k) = Set(1:N, k) - Trend(1:N, k)
            elsewhere
                Primes(1:N, k) = error
            end where
        else
            Primes(:, k) = error
        end if
    end do
end subroutine ExpWeightAvrg
