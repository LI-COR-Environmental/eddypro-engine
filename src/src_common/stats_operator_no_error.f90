!***************************************************************************
! stats_operator_no_error.f90
! ---------------------------
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
! \brief       Calculate column-wise sums on a 2d array \n
!              ignoring specified error values \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SumNoError(Set, nrow, ncol, Summ, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Summ(ncol)
    !> local variables
    integer :: i
    integer :: j
    logical :: data_exist

    Summ = 0d0
    do j = 1, ncol
        data_exist = .false.
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                data_exist = .true.
                Summ(j) = Summ(j) + Set(i, j)
            end if
        end do
        if (.not. data_exist) Summ(j) = err_float
    end do
end subroutine SumNoError

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

!***************************************************************************
!
! \brief       Calculate column-wise angular averages on a 2d array \n
!              ignoring specified error values. In EddyPro, mainly meant for \n
!              calculation of mean wind direction given a set of wind direction \n
!              measurements.
!
!              Implementation reference:
!              "Circular Statistics in R"
!              by A. Pewsey, M. Neuhaeuser and G. D. Ruxton.
!
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AngularAverageNoError(Set, nrow, ncol, Mean, err_float)
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
    real(kind = dbl) :: CosSum
    real(kind = dbl) :: SinSum


    do j = 1, ncol

        !> Calculate a (CosSum) and b (SinSum)
        CosSum = 0d0
        SinSum = 0d0
        Nact = 0
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                Nact = Nact + 1
                CosSum = CosSum + dcos(Set(i, j)/180d0*p)
                SinSum = SinSum + dsin(Set(i, j)/180d0*p)
            end if
        end do
        if (Nact /= 0) then
            CosSum = CosSum / dble(Nact)
            SinSum = SinSum / dble(Nact)
        else
            Mean(j) = err_float
            cycle
        end if

        !> Angular average is atan2 of b and a (the expressed in degrees)
        Mean(j) = datan2(SinSum, CosSum) * 180d0 / p

        !> Take angle form (-180, 180) to (0, 360)
        if (Mean(j) < 0d0) Mean(j) = 360d0 + Mean(j)
    end do
end subroutine AngularAverageNoError


!***************************************************************************
!
! \brief       Calculates standard deviation of array (column-wise) ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StDevNoError(Set, nrow, ncol, StDev, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: StDev(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: Nact = 0
    real(kind = dbl) :: Mean(ncol)


    !> Initializations
    StDev = 0d0

    !> Calculate mean values
    call AverageNoError(Set, nrow, ncol, Mean, err_float)

    !> Sum of squared residuals
    do j = 1, ncol
        if (Mean(j) == err_float) then
            StDev(j) = err_float
        else
            Nact = 0
            do i = 1, nrow
                if (Set(i, j) /= err_float) then
                    Nact = Nact + 1
                    StDev(j) = StDev(j) + (Set(i, j) - Mean(j)) **2
                end if
            end do
            if (Nact /= 0 .and. StDev(j) >= 0d0) then
                StDev(j) = dsqrt(StDev(j) / dble(Nact-1))
            else
                StDev = err_float
            end if
        end if
    end do

end subroutine StDevNoError

!***************************************************************************
!
! \brief       Calculates covariance matrix of given array, ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CovarianceMatrixNoError(Set, nrow, ncol, Cov, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Cov(ncol, ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: Nact = 0
    real(kind = dbl) :: sumi
    real(kind = dbl) :: sumj

    do i = 1, ncol
        do j = 1, ncol
            sumi = 0d0
            sumj = 0d0
            Cov(i, j) = 0d0
            Nact = 0
            do k = 1, nrow
                if (Set(k, i) /= err_float .and. Set(k, j) /= err_float) then
                    Nact = Nact + 1
                    Cov(i, j) = Cov(i, j) + Set(k, i) * Set(k, j)
                    sumi = sumi + Set(k, i)
                    sumj = sumj + Set(k, j)
                end if
            end do
            if (Nact /= 0) then
                sumi = sumi / dble(Nact)
                sumj = sumj / dble(Nact)
                Cov(i, j) = Cov(i, j) / dble(Nact)
                Cov(i, j) = Cov(i, j) - sumi * sumj
            else
                Cov(i, j) = err_float
            end if
        end do
    end do
end subroutine CovarianceMatrixNoError

!***************************************************************************
!
! \brief       Calculates covariance matrix of given arrays applying \n
!              given lag, ignoring provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function LaggedCovarianceNoError(col1, col2, nrow, rlag, err_float)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: rlag
    real(kind = dbl), intent(in) :: col1(nrow)
    real(kind = dbl), intent(in) :: col2(nrow)
    real(kind = dbl), intent(in) :: err_float
    !> Local variables
    integer :: lag
    integer :: i
    integer :: n
    real(kind = dbl) :: cov
    real(kind = dbl) :: sumi
    real(kind = dbl) :: sumj


    sumi = 0d0
    sumj = 0d0
    cov = 0d0
    n = 0
    if (rlag >= 0) then
        !> Positive lags are interpreted as col2 being "late"
        do i = 1, nrow - rlag
            if (col1(i) /= err_float .and. col2(i+rlag) /= err_float) then
                n = n + 1
                cov = cov + col1(i) * col2(i+rlag)
                sumi = sumi + col1(i)
                sumj = sumj + col2(i+rlag)
            end if
        end do
    else
        !> Positive lags are interpreted as col1 being "late"
        lag = -rlag
        do i = 1, nrow - lag
            if (col2(i) /= err_float .and. col1(i+lag) /= err_float) then
                n = n + 1
                cov = cov + col2(i) * col1(i+lag)
                sumi = sumi + col2(i)
                sumj = sumj + col1(i+lag)
            end if
        end do
    end if

    !> Finish up
    if (n /= 0) then
        sumi = sumi / dble(n)
        sumj = sumj / dble(n)
        cov = cov / dble(n)
        LaggedCovarianceNoError = cov - sumi * sumj
    else
        LaggedCovarianceNoError = err_float
    end if
end function LaggedCovarianceNoError
