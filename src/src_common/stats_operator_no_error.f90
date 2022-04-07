!***************************************************************************
! stats_operator_no_error.f90
! ---------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
!              Yamartino, 1984:
!              https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281984%29023%3C1362%3AACOSPE%3E2.0.CO%3B2
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
                CosSum = CosSum - dcos(Set(i, j) / 180d0 * p)
                SinSum = SinSum - dsin(Set(i, j) / 180d0 * p)
            end if
        end do
        if (Nact /= 0) then
            CosSum = CosSum / dble(Nact)
            SinSum = SinSum / dble(Nact)
        else
            Mean(j) = err_float
            cycle
        end if

        !> Angular average is atan2 of b and a
        !> "+p" adjust quadrant, then express in degrees
        Mean(j) = (datan2(SinSum, CosSum) + p) * 180d0 / p
    end do
end subroutine AngularAverageNoError

!***************************************************************************
!
! \brief       Calculate column-wise angular stdev on a 2d array \n
!              ignoring specified error values. In EddyPro, mainly meant for \n
!              calculation of wind direction standard deviation given a set of wind direction. \n
!
!              Implementation reference:
!              Yamartino, 1984:
!              https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281984%29023%3C1362%3AACOSPE%3E2.0.CO%3B2
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AngularStDevApproxNoError(Set, nrow, ncol, AngStDev, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: AngStDev(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: Nact = 0
    real(kind = dbl) :: CosSum
    real(kind = dbl) :: SinSum
    real(kind = dbl) :: eps


    do j = 1, ncol

        !> Calculate a (CosSum) and b (SinSum)
        CosSum = 0d0
        SinSum = 0d0
        Nact = 0
        do i = 1, nrow
            if (Set(i, j) /= err_float) then
                Nact = Nact + 1
                CosSum = CosSum - dcos(Set(i, j) / 180d0 * p)
                SinSum = SinSum - dsin(Set(i, j) / 180d0 * p)
            end if
        end do
        if (Nact /= 0) then
            CosSum = CosSum / dble(Nact)
            SinSum = SinSum / dble(Nact)
        else
            AngStDev(j) = err_float
            cycle
        end if

        !> Approximate standard deviation of wind direction (see Yamartino, 1984:
        !> https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281984%29023%3C1362%3AACOSPE%3E2.0.CO%3B2)
        eps = dsqrt(1. - (CosSum**2 + SinSum**2))
        AngStDev = (dasin(eps) * (1. + (2 / sqrt(3.) - 1) * eps**3)) * 180d0 / p
    end do
end subroutine AngularStDevApproxNoError

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
! \brief       Calculates correlation matrix of given array, ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CorrelationMatrixNoError(Set, nrow, ncol, Corr, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Corr(ncol, ncol)
    !> local variables
    real(kind = dbl) :: Cov(ncol, ncol)
    real(kind = dbl) :: StDev(ncol)
    integer :: i
    integer :: j
    integer :: var

    
    call CovarianceMatrixNoError(Set, size(Set, 1), size(Set, 2), Cov, err_float)

    do var = u, gas4
        call StDevNoError(Set, size(Set, 1), size(Set, 2), StDev, err_float)
    end do

    do i = u, gas4
        do j = u, gas4
            if (Cov(i, j) /= err_float .and. StDev(i) > 0d0 .and. StDev(j) > 0d0) then
                Corr(i, j) = Cov(i, j) / (StDev(i) * StDev(j))
            else
                Corr(i, j) = err_float
            end if
        end do 
    end do
end subroutine CorrelationMatrixNoError

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

!***************************************************************************
!
! \brief       Calculates Skewness of timeserie ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SkewnessNoError(Set, nrow, ncol, Skw, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(out) :: Skw(ncol)
    !> Local variables
    integer :: i
    integer :: j
    integer :: Nact
    real(kind = dbl) :: StDev(ncol)


    !> Compute StDev, needed for Skewness
    call StDevNoError(Set, nrow, ncol, StDev, err_float)

    Skw = 0.d0
    do j = u, ncol
        if (E2Col(j)%present) then
            Nact = 0
            do i = 1, nrow
                if (Set(i, j) /= error) then
                    Nact = Nact + 1
                    Skw(j) = Skw(j) + (Set(i, j))**3
                end if
            end do
            if (Nact /= 0) then
                Skw(j) = Skw(j) / (StDev(j)**3) / dble(Nact - 1)
            else
                Skw(j) = error
            end if
        else
            Skw(j) = error
        end if
    end do
end subroutine SkewnessNoError

!***************************************************************************
!
! \brief       Calculates Kurtosis index  of timeserie ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine KurtosisNoError(Set, nrow, ncol, Kur, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(out) :: Kur(ncol)
    !> Local variables
    integer :: i
    integer :: j
    integer :: Nact
    real(kind = dbl) :: StDev(ncol)


    !> Compute StDev, needed for Kurtosis
    call StDevNoError(Set, nrow, ncol, StDev, err_float)

    Kur = 0.d0
    do j = u, ncol
        if (E2Col(j)%present) then
            Nact = 0
            do i = 1, nrow
                if (Set(i, j) /= error) then
                    Nact = Nact + 1
                    Kur(j) = Kur(j) + (Set(i, j))**4
                end if
            end do
            if (Nact /= 0) then
                Kur(j) = Kur(j) / (StDev(j)**4) / dble(Nact - 1)
            else
                Kur(j) = error
            end if
        else
            Kur(j) = error
        end if
    end do
end subroutine KurtosisNoError

!***************************************************************************
!
! \brief       Calculates quantiles of array (column-wise) ignoring \n
!              provided error code
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine QuantileNoError(Set, nrow, ncol, Quantile, qin, err_float)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: err_float
    real(kind = dbl), intent(in) :: qin
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: Quantile(ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: cnt = 0
    integer :: M = 0
    logical :: mask(nrow)
    real(kind = dbl), allocatable :: x(:)
    real(kind = dbl), external :: quantile_sas5


    !> Initializations
    Quantile = err_float

    !> Sum of squared residuals
    do j = 1, ncol
        mask(:) = Set(:, j) /= err_float
        M = count(mask)
        if (M > 0) then
            allocate(x(M))
            cnt = 0
            do i = 1, nrow
                if (Set(i, j) /= err_float) then 
                    cnt = cnt + 1
                    x(cnt) = Set(i, j)
                end if
            end do
            if (size(x) > 1) then
                Quantile(j) = quantile_sas5(x, size(x), qin)
            end if
            deallocate(x)
        end if
    end do
end subroutine QuantileNoError


subroutine unbiased_correlation(arr1, arr2, n, err_float, lag, r, t, m)
    ! Computes the correlation of two series A and B of length
    ! N as a function of lag. The correlation is unbiased
    ! because the sums in the covariance and standard devia-
    ! tions are divided by m, not m-1, where m is the number
    ! of overlapping grid points.

    ! The correlation r is defined as 
    !
    !           cov(A,B)
    !     r = -------------
    !          s(A) * s(B)
    !
    ! where cov() is covariance and s() is standard deviation.

    ! The series A and B must be prepared such that they
    ! are sent to this routine with zero lag. If necessary, this
    ! is accomplished by padding the beginning or ends (or both)
    ! of each time series with missing values (err_float) so that their
    ! elements correspond one to one. The resultant time series
    ! will be of equal length N. Also, missing values (err_float)
    ! distributed thoughout each time series is perfectly
    ! acceptable. Time series B will be shifted by the amount
    ! lag relative to time series A:

    ! RETURNS: r, the unbiased correlation, t, the significance
    ! of the unbiased correlation ( t is set to err_float if B = A
    ! at lag 0, i.e. an autocorrelation at lag 0, or if the
    ! correlation is 1.0 in general), and m, the number of
    ! overlapping grid points.

    ! A:      t1  t2  t3  t4        ...       tN
    ! B:              t1  t2  t3  t4        ...       tN 

    !         \_ _/
    !           V
    !           lag = 2 (Defined to be > 0.) 

    implicit none
    !> In/out variables
    integer, intent(in) :: n
    real, dimension(n), intent(in) :: arr1 
    real, dimension(n), intent(in) :: arr2           
    real, intent(in) :: err_float 
    integer, intent(in) :: lag 
    real, intent(out) :: r
    real, intent(out) :: t
    integer, intent(out) :: m
    !> Local variables
    real, dimension(:), allocatable :: x
    real, dimension(:), allocatable :: y
    real, dimension(:), allocatable :: xdev
    real, dimension(:), allocatable :: ydev
    real, dimension(:), allocatable :: xdevydev
    real, dimension(:), allocatable :: xdevxdev
    real, dimension(:), allocatable :: ydevydev
    real :: xmn
    real :: ymn
    real :: covxy
    real :: sx
    real :: sy


    allocate( x(1:3*n) )
    allocate( y(1:3*n) )

    x(:) = err_float
    y(:) = err_float

    !> Extend arrays with copy for easier handling of shifting
    x(n+1:2*n) = arr1(:)
    y(n+1:2*n) = arr2(:)

    !> Shift y
    y = eoshift(y, shift = -lag, boundary = err_float)
    where (x == err_float) y = err_float
    where (y == err_float) x = err_float
    m = count(x /= err_float)

    !> Mean values
    xmn = sum(x, dim=1, mask=x /= err_float) / float(m)
    ymn = sum(y, dim=1, mask=y /= err_float) / float(m)

    !> Fluctuations around mean
    allocate(xdev(1:3*n))
    allocate(ydev(1:3*n))
    xdev(:) = x(:) - xmn
    where (x == err_float) xdev(:) = err_float
    ydev(:) = y(:) - ymn
    where (y == err_float) ydev(:) = err_float

    !> Covariance
    allocate( xdevydev(1:3*n) )
    xdevydev(:) = xdev(:) * ydev(:)
    where (x == err_float) xdevydev(:) = err_float
    covxy = sum(xdevydev, dim=1, mask=xdevydev /= err_float) / float(m)

    !> Standard Deviations
    allocate(xdevxdev(1:3*n))
    xdevxdev(:) = xdev(:) * xdev(:)
    where (x == err_float) xdevxdev(:) = err_float
    sx = sqrt(sum(xdevxdev, dim=1, mask=xdevxdev /= err_float) / float(m))

    allocate(ydevydev(1:3*n))
    ydevydev(:) = ydev(:) * ydev(:)
    where (y == err_float) ydevydev(:) = err_float
    sy = sqrt(sum(ydevydev, dim=1, mask=ydevydev /= err_float) / float(m))

    !> Correlation coefficient 
    r = covxy / ( sx * sy )

    if (r /= 1.0) then
        t = r * sqrt( (m - 2) / ( 1 - r*r ) )
    else
        t = err_float
    end if
    deallocate( x, y, xdev, ydev, xdevydev, xdevxdev, ydevydev )
end subroutine unbiased_correlation 


!***************************************************************************
!
! \brief       Cross-correlation function for passed arrays and specified
!              lag-time boundaries
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
! function CrossCorrelation(arr1, arr2, nrow, lagmin, lagmax) result(CCF)
!     use m_common_global_var
!     implicit none
!     !> in/out variables
!     integer, intent(in) :: nrow
!     real(kind = dbl), intent(in) :: arr1(nrow)
!     real(kind = dbl), intent(in) :: arr2(nrow)
!     integer, intent(in) :: lagmin
!     integer, intent(in) :: lagmax
!     real(kind = dbl), dimension(51) :: CCF
!     !> local variables
!     integer :: lag
!     real(kind = dbl), external :: LaggedCovarianceNoError

!     print*, arr1(1), arr2(1), nrow, lagmin, lagmax
!     stop
!     ! do lag = lagmin, lagmax
!     !     CCF(lag) = LaggedCovarianceNoError(arr1, arr2, nrow, lag, error)
!     ! end do

! end function CrossCorrelation


subroutine CrossCorrelation(arr1, arr2, nrow, lagmin, lagmax, CCF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: arr1(nrow)
    real(kind = dbl), intent(in) :: arr2(nrow)
    integer, intent(in) :: lagmin
    integer, intent(in) :: lagmax
    real(kind = dbl), intent(out) :: CCF(lagmin: lagmax)
    !> local variables
    integer :: lag
    real(kind = dbl) :: sig1(1), sig2(1)
    real(kind = dbl), external :: LaggedCovarianceNoError


    !> Cross-covariance function
    do lag = lagmin, lagmax
        CCF(lag) = LaggedCovarianceNoError(arr1, arr2, nrow, lag, error)
    end do

    !> Normalize by standard deviations
    call StDevNoError(arr1, nrow, 1, sig1, error)
    call StDevNoError(arr2, nrow, 1, sig2, error)
    CCF = CCF / (sig1(1) * sig2(1))

end subroutine CrossCorrelation