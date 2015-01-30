!***************************************************************************
! user_basic_stats.f90
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
! \brief       Calculates main statistics on non-sensitive variables, using \n
!              different algorithms at different processing levels
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UserBasicStats(UserSet, N, M, nfold)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    integer, intent(in) :: nfold
    real(kind = dbl), intent(in) :: UserSet(N, M)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: Nact = 0
    real(kind = dbl) :: RawMean(M)
    real(kind = dbl) :: Prime(N, M)
    real(kind = dbl) :: sumi
    real(kind = dbl) :: sumj

    if (nfold == 1) then
        write(*, '(a)', advance = 'no') '  Calculating user statistics..'
    else
        write(*, '(a,i1,a)', advance = 'no') '  Re-calculating user statistics (',nfold,')..'
    end if

    !> covariances
    do i = 1, M
        do j = 1, M
            sumi = 0d0
            sumj = 0d0
            UserStats%Cov(i, j) = 0d0
            Nact = 0
            do k = 1, N
                if (UserSet(k, i) /= error .and. UserSet(k, j) /= error) then
                    Nact = Nact + 1
                    UserStats%Cov(i, j) = UserStats%Cov(i, j) + UserSet(k, i) * UserSet(k, j)
                    sumi = sumi + UserSet(k, i)
                    sumj = sumj + UserSet(k, j)
                end if
            end do
            if (Nact /= 0) then
                sumi = sumi / dble(Nact)
                sumj = sumj / dble(Nact)
                UserStats%Cov(i, j) = UserStats%Cov(i, j) / dble(Nact)
                UserStats%Cov(i, j) = UserStats%Cov(i, j) - sumi * sumj
            else
                UserStats%Cov(i, j) = error
            end if
        end do
    end do

    if (nfold <= 6) then
        !> mean values (only before detrending, after is deleterious)
        RawMean = 0d0
        do j = 1, M
            Nact = 0
            do i = 1, N
                if (UserSet(i, j) /= error) then
                    Nact = Nact + 1
                    RawMean(j) = RawMean(j) + UserSet(i, j)
                end if
            end do
            if (Nact /= 0) then
                RawMean(j) = RawMean(j) / dble(Nact)
            else
                RawMean(j) = error
            end if
        end do
        UserStats%Mean = 0.d0
        do j = 1, M
            if (RawMean(j) /= error) then
                Nact = 0
                do i = 1, N
                    if (UserSet(i, j) /= error) then
                        Nact = Nact + 1
                        UserStats%Mean(j) = UserStats%Mean(j) + UserSet(i, j) - RawMean(j)
                    end if
                end do
                if (Nact /= 0) then
                    UserStats%Mean(j) = UserStats%Mean(j) / dble(Nact)
                else
                    UserStats%Mean(j) = error
                end if
            else
                UserStats%Mean(j) = error
            end if
        end do

        where (UserStats%Mean(1:M) /= error)
            UserStats%Mean(1:M) = UserStats%Mean(1:M) + RawMean(:)
        elsewhere
            UserStats%Mean(1:M) = error
        end where

        !> fluctuations (only before detrending, after is deleterious)
        do j = 1, M
            do i = 1, N
                if (UserSet(i, j) /= error) then
                    Prime(i, j) = UserSet(i, j) - UserStats%Mean(j)
                else
                    Prime(i, j) = error
                end if
            end do
        end do
    else
        !> after detrending, UserSet contains fluctuations
        Prime = UserSet
    end if

    !> Standard deviations
    do j = 1, M
        UserStats%StDev(j) = dsqrt(UserStats%Cov(j, j))
    end do

    !> skewness and kurtosis
    UserStats%Skw = 0.d0
    UserStats%Kur = 0.d0
    do j = 1, M
        Nact = 0
        do i = 1, N
            if (Prime(i, j) /= error) then
                Nact = Nact + 1
                UserStats%Skw(j) = UserStats%Skw(j) + (Prime(i, j))**3
                UserStats%Kur(j) = UserStats%Kur(j) + (Prime(i, j))**4
            end if
        end do
        if (Nact /= 0) then
            UserStats%Skw(j) = UserStats%Skw(j) / (UserStats%StDev(j)**3) / dble(Nact - 1)
            UserStats%Kur(j) = UserStats%Kur(j) / (UserStats%StDev(j)**4) / dble(Nact - 1)
        else
            UserStats%Skw(j) = error
            UserStats%Kur(j) = error
        end if
    end do
    write(*,'(a)') ' Done.'
end subroutine UserBasicStats
