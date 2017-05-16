!***************************************************************************
! timelag_handle.f90
! ------------------
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
! \brief       Calculates time lags (in terms of data rows) for all scalars \n
!              not measured by the anemometer. Also calculates covariances \n
!              of H2O and Cell T with time-lags of other scalars (from the \n
!              same instrument) for proper WPL of closed path systems.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TimeLagHandle(TlagMeth, Set, nrow, ncol, ActTLag, TLag, &
    DefTlagUsed, InTimelagOpt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    character(*), intent(in) :: TlagMeth
    logical, intent(in) :: InTimelagOpt
    logical, intent(out) :: DefTlagUsed(ncol)
    real(kind = dbl), intent(out) :: ActTLag(ncol)
    real(kind = dbl), intent(out) :: TLag(ncol)
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: def_rl(ncol)
    integer :: min_rl(ncol)
    integer :: max_rl(ncol)
    real(kind = dbl) :: ColW(nrow)
    real(kind = dbl) :: ColH2O(nrow)
    real(kind = dbl) :: ColTC(nrow)
    real(kind = dbl) :: FirstCol(nrow)
    real(kind = dbl) :: SecondCol(nrow)
    real(kind = dbl) :: TmpSet(nrow, ncol)

    if  (.not. InTimelagOpt) write(*, '(a)', advance = 'no') &
        '  Compensating time-lags..'

    !> for E2Set scalars, initialise auxiliary vars to zero
    def_rl(:) = 0
    min_rl(:) = 0
    min_rl(:) = 0
    !> Define "row-lags" for scalars, using time-lags
    !> retrieved from metadata file
    where (E2Col(ts:pe)%present)
        def_rl(ts:pe) = nint(E2Col(ts:pe)%def_tl * Metadata%ac_freq)
        min_rl(ts:pe) = nint(E2Col(ts:pe)%min_tl * Metadata%ac_freq)
        max_rl(ts:pe) = nint(E2Col(ts:pe)%max_tl * Metadata%ac_freq)
    end where

    !> calculate actual time-lags according to the chosen method
    select case(TlagMeth)
        case ('constant')
            !> constant timelags are set equal to default values (user selected)
            RowLags(ts:pe) = def_rl(ts:pe)
            TLag(ts:pe)    = E2Col(ts:pe)%def_tl
            ActTLag(ts:pe) = E2Col(ts:pe)%def_tl
            DefTlagUsed(ts:pe) = .true.
        case ('maxcov', 'maxcov&default')
            !> covariance maximization method, with or without default
            do j = ts, pe
                !> Only for present variables,
                !> with both min and max "row lags" /= 0
                if (E2Col(j)%present &
                    .and. (min_rl(j) /= 0 .or. max_rl(j) /= 0)) then
                    FirstCol(:)  = Set(:, w)
                    SecondCol(:) = Set(:, j)
                    call CovMax(min_rl(j), max_rl(j), &
                        FirstCol, SecondCol, size(FirstCol), &
                        TLag(j), RowLags(j))
                    ActTLag(j) = TLag(j)
                    !> If no max cov has been detected within the interval, \n
                    !> sets the time lag to the suggested values
                    if (TlagMeth == 'maxcov&default') then
                        if ( (RowLags(j) == min_rl(j)) .or. (RowLags(j) == max_rl(j)) ) then
                            DefTlagUsed(j) = .true.
                            TLag(j) = dble(def_rl(j)) / Metadata%ac_freq
                            RowLags(j) = def_rl(j)
                        end if
                    end if
                else
                    RowLags(j) = ierror
                    TLag(j) = error
                    ActTLag(j) = error
               end if
            end do
        case ('none')
            !> not compensating for timelags
            RowLags(ts:pe) = 0
            TLag(ts:pe) = 0d0
    end select

    if  (.not. InTimelagOpt) then
        !> For closed path instruments, calculate H2O covariances
        !> for time-lags of other scalars from the same instrument
        Stats%h2ocov_tl_co2 = error
        Stats%h2ocov_tl_ch4 = error
        Stats%h2ocov_tl_gas4 = error
        if (E2Col(h2o)%present &
            .and. E2Col(h2o)%instr%path_type == 'closed') then
            ColW(1:nrow) = Set(1:nrow, w)
            ColH2O(1:nrow) = Set(1:nrow, h2o)
            if (E2Col(co2)%present &
                .and. E2Col(co2)%instr%model == E2Col(h2o)%instr%model &
                .and. RowLags(co2) > 0) &
                call CovarianceW(ColW, ColH2O, size(ColW), &
                    RowLags(co2), Stats%h2ocov_tl_co2)
            if (E2Col(ch4)%present &
                .and. E2Col(ch4)%instr%model == E2Col(h2o)%instr%model &
                .and. RowLags(ch4) > 0) &
                call CovarianceW(ColW, ColH2O, size(ColW), &
                    RowLags(ch4), Stats%h2ocov_tl_ch4)
            if (E2Col(gas4)%present &
                .and. E2Col(gas4)%instr%model == E2Col(h2o)%instr%model &
                .and. RowLags(gas4) > 0) &
                call CovarianceW(ColW, ColH2O, size(ColW), &
                RowLags(gas4), Stats%h2ocov_tl_gas4)
        end if

        !> Calculate cell temperature covariances with
        !> time-lags of scalars from the same instrument
        Stats%tc_cov_tl_co2 = error
        Stats%tc_cov_tl_h2o = error
        Stats%tc_cov_tl_ch4 = error
        Stats%tc_cov_tl_gas4 = error
        if (E2Col(tc)%present) then
            !> Store vertical wind component and tc in ad-hoc arrays
            ColW(1:nrow) = Set(1:nrow, w)
            ColTC(1:nrow) = Set(1:nrow, tc)
            if (E2Col(co2)%present &
                .and. E2Col(co2)%instr%model == E2Col(tc)%instr%model &
                .and. RowLags(co2) > 0) &
                call CovarianceW(ColW, ColTC, size(ColTC), &
                    RowLags(co2), Stats%tc_cov_tl_co2)
            if (E2Col(h2o)%present &
                .and. E2Col(h2o)%instr%model == E2Col(tc)%instr%model &
                .and. RowLags(h2o) > 0) &
                call CovarianceW(ColW, ColTC, size(ColTC), &
                    RowLags(h2o), Stats%tc_cov_tl_h2o)
            if (E2Col(ch4)%present &
                .and. E2Col(ch4)%instr%model == E2Col(tc)%instr%model &
                .and. RowLags(ch4) > 0) &
                call CovarianceW(ColW, ColTC, size(ColTC), &
                    RowLags(ch4), Stats%tc_cov_tl_ch4)
            if (E2Col(gas4)%present &
                .and. E2Col(gas4)%instr%model == E2Col(tc)%instr%model &
                .and. RowLags(gas4) > 0) &
                call CovarianceW(ColW, ColTC, size(ColTC), &
                    RowLags(gas4), Stats%tc_cov_tl_gas4)
        end if
    end if

    !> Align data according to relevant time-lags,
    !> filling remaining with error code.
    do j = u, pe
        if (E2Col(j)%present) then
            if (RowLags(j) >= 0) then
                !> For positive lags
                do i = 1, nrow - RowLags(j)
                    TmpSet(i, j) = Set(i + RowLags(j), j)
                end do
                do i = nrow - Rowlags(j) + 1, nrow
                    TmpSet(i, j) = error
                end do
            else
                !> For negative lags
                do i = 1, abs(RowLags(j))
                    TmpSet(i, j) = error
                end do
                do i = abs(RowLags(j)) + 1, nrow
                    TmpSet(i, j) = Set(i + RowLags(j), j)
                end do
            end if
        else
            TmpSet(1:nrow, j) = error
        end if
    end do
    Set = TmpSet
    if  (.not. InTimelagOpt) write(*,'(a)') ' Done.'
end subroutine TimeLagHandle

!*******************************************************************************
!
! \brief       Performs covariance analysis for determining the "optimal" \n
!              time lag, the one that maximizes the covariance.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine CovMax(lagmin, lagmax, Col1, Col2, nrow, TLag, RLag)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: lagmin
    integer, intent(in) :: lagmax
    real(kind = dbl), intent(in) :: Col1(nrow)
    real(kind = dbl), intent(in) :: Col2(nrow)
    integer, intent(out) :: RLag
    real(kind = dbl), intent(out) :: TLag
    !> local variables
    integer :: i = 0
    integer :: ii = 0
    integer :: N2
    real(kind = dbl), allocatable :: ShSet(:, :)
    real(kind = dbl), allocatable :: ShPrimes(:, :)
    real(kind = dbl) :: CovMat(2,2)
    real(kind = dbl) :: Cov
    real(kind = dbl) :: MaxCov

    Cov = 0.d0
    MaxCov = 0.d0
    TLag = 0.d0
    do i = lagmin, lagmax
        N2 = nrow - abs(i)
        allocate(ShSet(N2, 2))
        allocate(ShPrimes(N2, 2))

        !> Align the two timeseries at the current time-lag 
        do ii = 1, N2
            if (i < 0) then
                ShSet(ii, 1) = Col1(ii - i)
                ShSet(ii, 2) = Col2(ii)
            else
                ShSet(ii, 1) = Col1(ii)
                ShSet(ii, 2) = Col2(ii + i)
            end if
        end do

        !> Block average
        ShPrimes = ShSet

        !> Linear detrending
        ! call VariableLinearDetrending(ShSet(:, 1), ShPrimes(:, 1), N2)
        ! call VariableLinearDetrending(ShSet(:, 2), ShPrimes(:, 2), N2)

        !> Stochastic detrending
        ! call VariableStochasticDetrending(ShSet(:, 1), ShPrimes(:, 1), N2)
        ! call VariableStochasticDetrending(ShSet(:, 2), ShPrimes(:, 2), N2)

        call CovarianceMatrixNoError(ShPrimes, size(ShPrimes, 1), size(ShPrimes, 2), CovMat, error)
        Cov = CovMat(1, 2)

        !> Max cov and actual time lag
        if (abs(Cov) > MaxCov) then
            MaxCov = abs(Cov)
            TLag = dble(i) / Metadata%ac_freq
            RLag = i
        end if
        deallocate(ShSet)
        deallocate(ShPrimes)
    end do
end subroutine CovMax


!***************************************************************************
!
! \brief       Calculate covariance between two arrays using an imposed  \n
!              time-lag.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CovarianceW(col1, col2, nrow, lag, cov)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: lag
    real(kind = dbl), intent(in) :: col1(nrow)
    real(kind = dbl), intent(in) :: col2(nrow)
    real(kind = dbl), intent(out) :: cov
    !> local variables
    integer :: i
    integer :: N2
    real(kind = dbl) ::sum1
    real(kind = dbl) ::sum2

    sum1 = 0d0
    sum2 = 0d0
    Cov = 0d0
    N2 = 0
    do i = 1, nrow - lag
        if (col1(i) /= error .and. col2(i+lag) /= error) then
            N2 = N2 + 1
            Cov = Cov + col1(i) * col2(i+lag)
            sum1 = sum1 + col1(i)
            sum2 = sum2 + col2(i+lag)
        end if
    end do

    if (N2 /= 0) then
        sum1 = sum1 / dble(N2)
        sum2 = sum2 / dble(N2)
        cov = cov / dble(N2)
        cov = cov - sum1 * sum2
    else
        cov = error
    end if
end subroutine CovarianceW

!***************************************************************************
!
! \brief       Stochastic Detrending
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine VariableStochasticDetrending(Var, Primes, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Primes(N)
    !> local variables
    integer :: i

    Primes(1) = error
    do i = 2, N
        if (Var(i) /= error .and. Var(i-1) /= error) then
            Primes(i) = Var(i) - Var(i-1)
        else 
            Primes(i) = error
        end if
    end do
end subroutine VariableStochasticDetrending

!***************************************************************************
!
! \brief       Linear detrending of one time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine VariableLinearDetrending(Var, Primes, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Primes(N)
    !> Local variables
    real(kind = dbl) :: Trend(N)

    call CalculateTrend(Var, Trend, N)
    call Detrend(Var, Trend, Primes, N)

end subroutine VariableLinearDetrending

!***************************************************************************
!
! \brief       Remove trend from time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Detrend(Var, Trend, Primes, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(in) :: Trend(N)
    real(kind = dbl), intent(out) :: Primes(N)


    Primes = error
    where (Var /= error .and. Trend /= error)
        Primes = Var - Trend
    end where
end subroutine Detrend

!***************************************************************************
!
! \brief       Calculate linear trend in time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CalculateTrend(Var, Trend, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Trend(N)
    !> local variables
    integer :: i
    integer :: nn
    integer :: mm
    real(kind = dbl) :: sumx1
    real(kind = dbl) :: sumx2
    real(kind = dbl) :: mean
    real(kind = dbl) :: sumtime
    real(kind = dbl) :: sumtime2
    real(kind = dbl) :: b


    !> Linear regression
    sumx1 = 0d0
    sumx2 = 0d0
    sumtime = 0d0
    sumtime2 = 0d0
    nn = 0
    do i = 1, N
        if (Var(i) /= error) then
            nn = nn + 1
            sumx1 = sumx1 + (Var(i) * (dble(nn - 1)))
            sumx2 = sumx2 + Var(i)
            sumtime = sumtime + (dble(nn - 1))
            sumtime2 = sumtime2 + (dble(nn - 1))**2
        end if
    end do
    if (nn /= 0) then
        mean = sumx2 / dble(nn)
    end if

    !> Trend
    mm = 0
    b = (sumx1 - (sumx2 * sumtime) / dble(nn)) / (sumtime2 - (sumtime * sumtime) / dble(nn))
    do i = 1, N
        mm = mm + 1
        if (Var(i) /= error) then
            Trend(i) = mean + b * (dble(mm - 1) - sumtime / dble(nn))
        else
            Trend(i) = error
        end if
    end do        
end subroutine CalculateTrend
