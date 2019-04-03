!***************************************************************************
! random_error_handle.f90
! -----------------------
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
! \brief       Estimate flux random uncertainty according to the selected method
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RandomUncertaintyHandle(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)


    write(*, '(a)') '  Estimating random uncertainty..'

    !> Calculate random uncertainty
    select case (RUsetup%meth)
        case('finkelstein_sims_01')
            call IntegralTurbulenceScale(Set, size(Set, 1), size(Set, 2))
            call RU_Finkelstein_Sims_01(Set, nrow, ncol)
        case('mann_lenschow_94')
            call IntegralTurbulenceScale(Set, size(Set, 1), size(Set, 2))
            call RU_Mann_Lenschow_04(nrow)
        case('none')
            Essentials%rand_uncer(u:gas4) = error
            Essentials%rand_uncer_LE = error
            Essentials%rand_uncer_ET = error
        case('mahrt_98')
            !> Mahrt has been calculated already, so don't need to do anything
            continue
        case default
            call ExceptionHandler(42)
            Essentials%rand_uncer(u:gas4) = error
            Essentials%rand_uncer_LE = error
            Essentials%rand_uncer_ET = error
            return
    end select
    write(*, '(a)') '  Done.'
end subroutine RandomUncertaintyHandle

!***************************************************************************
!
! \brief       Estimate random error according to \n
!              Finkelstein and Sims (2001), Eq. 8- 10
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RU_Finkelstein_Sims_01(Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    !> local variables
    integer :: var
    integer :: lag
    integer :: LagMax(M)
    integer :: errcnt
    real(kind = dbl), allocatable :: gam(:, :, :)
    real(kind = dbl) :: varcov
    real(kind = dbl), external :: LaggedCovarianceNoError

    !> Define max lag based on ITS
    LagMax(u:gas4) = nint(ITS(u:gas4) * Metadata%ac_freq)
    where (LagMax < 0) LagMax = nint(error)
    do var = u, gas4
        if (var == v .or. var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error .and. LagMax(var) /= nint(error)) then
            allocate (gam(0:LagMax(var), 2, 2))
            gam = 0d0
            do lag = 0, LagMax(var)
                gam(lag, 1, 1) = &
                    LaggedCovarianceNoError(Set(:, w), Set(:, w), &
                                            size(Set, 1), lag, error)
                gam(lag, 2, 2) = &
                    LaggedCovarianceNoError(Set(:, var), Set(:, var), &
                                            size(Set, 1), lag, error)
                gam(lag, 1, 2) = &
                    LaggedCovarianceNoError(Set(:, w), Set(:, var), &
                                            size(Set, 1), lag, error)
                gam(lag, 2, 1) = &
                    LaggedCovarianceNoError(Set(:, w), Set(:, var), &
                                            size(Set, 1), -lag, error)
            end do

            !> variance of covariances, Eq. 8  in Finkelstein & Sims (2001, JGR)
            !> Initialize the value for lag = 0
            varcov = 0d0
            if (gam(0, 1, 1) /= error .and. gam(0, 2, 2) /= error) &
                varcov = gam(0, 1, 1) * gam(0, 2, 2) + gam(0, 1, 2) * gam(0, 2, 1)

            !> Now cycle on lag. Do it one sided and multiply by 2 (Eq. 9 and 10)
            errcnt = 0
            do lag = 1, LagMax(var)
                if (gam(lag, 1, 1) /= error .and. gam(0, 2, 2) /= error &
                    .and. gam(0, 1, 2) /= error .and. gam(0, 2, 1) /= error) then
                    varcov = varcov + 2d0 * gam(lag, 1, 1) * gam(lag, 2, 2) &
                        + 2d0 * gam(lag, 1, 2) * gam(lag, 2, 1)
                else
                    errcnt = errcnt + 1
                end if
            enddo
            deallocate (gam)
            !> Normalization (see Eq. 8 for why using N here)
            varcov = varcov / dfloat(N - errcnt)

            !> Random error is the square root of this variance
            if (varcov /= 0) then
                Essentials%rand_uncer(var) = dsqrt(abs(varcov))
            else
                Essentials%rand_uncer(var) = error
            end if
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RU_Finkelstein_Sims_01

!***************************************************************************
!
! \brief       Estimate random error according to Mann and Lenschow (1994)
!              See e.g. Eq. 5 in Finkelstein and Sims (2001)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RU_Mann_Lenschow_04(N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: var
    real(kind = dbl) :: corr_coeff(E2NumVar)

    do var = u, gas4
        if (var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error) then

            !> Correlation coefficient
            corr_coeff(var) = dabs(Stats%cov(w, var)) &
                / (dsqrt(Stats%cov(w, w)) * dsqrt(Stats%cov(var, var)))

            !> Random uncertainty
            Essentials%rand_uncer(var) = abs(Stats%cov(w, var)) &
                * dsqrt((1d0 + corr_coeff(var)**2) / corr_coeff(var)**2) &
                * dsqrt (2d0 * ITS(var) / (N / Metadata%ac_freq))
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RU_Mann_Lenschow_04


!***************************************************************************
!
! \brief       Estimate random error according to \n
!              Mahrt (1998), Eqs. 8 - 9
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RU_Mahrt_98(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    integer :: j
    integer :: Ni
    integer :: Nj
    integer :: var
    integer, parameter :: nrec = 6
    integer, parameter :: nsubrec = 6
    real(kind = dbl) :: covmat(GHGNumVar, GHGNumVar)
    real(kind = dbl)  :: Fij(nsubrec, GHGNumVar)
    real(kind = dbl)  :: Fijs(nrec * nsubrec, GHGNumVar)
    real(kind = dbl)  :: Fi_bar(GHGNumVar)
    real(kind = dbl)  :: Fi_bars(nrec, GHGNumVar)
    real(kind = dbl)  :: F_bar(GHGNumVar)
    real(kind = dbl) :: SumSquares(GHGNumVar)
    real(kind = dbl) :: sigma_wis(nrec, GHGNumVar)
    real(kind = dbl) :: sigma_btw(GHGNumVar)
    real(kind = dbl), allocatable :: sSet(:, :)
    real(kind = dbl), allocatable :: ssSet(:, :)

    Ni = nrow / nrec
    Nj = Ni / nsubrec
    if (.not. allocated(sSet)) allocate(sSet(Ni, GHGNumVar))
    do i = 1, nrec
        sSet(:, :) = Set(Ni * (i-1) + 1: Ni * i, 1:GHGNumVar)
        !> Compute covariance matrices on sub-sub-periods
        do j = 1, nsubrec
            if (.not. allocated(ssSet)) allocate(ssSet(Nj, GHGNumVar))
            ssSet(:, :) = sSet(Nj * (j-1) + 1: Nj * j, :)
            call CovarianceMatrixNoError(ssSet, size(ssSet, 1), size(ssSet, 2), covmat, error)
            Fij(j, :) = covmat(w, :)
            if (allocated(ssSet)) deallocate(ssSet)
        end do

        !> Mean of covariances on sub-sub-periods, per variable, per sub-period, F(i)bar in Eq. 8
        call AverageNoError(Fij, nsubrec, GHGNumVar, Fi_bar, error)

        !> Accumulate covariances on sub-sub-records, needed to compute Fbar in Eq. 10
        Fijs(nsubrec * (i-1) + 1: nsubrec * i, :) = Fij(:, :)

        !> Accumulate mean covariances, per superiod, needed in Eq. 10
        Fi_bars(i, :) = Fi_bar(:)

        !> Sum of squares of residuals, per variable, per sub-period, Eq. 8
        SumSquares = 0d0
        do j = 1, nsubrec
            SumSquares(:) = SumSquares(:) + (Fij(j, :) - Fi_bar(:))**2
        end do

        !> Standard deviation within, per variable, per sub-period, sig_wi(i) in Eq. 8
        sigma_wis(i, :) = dsqrt(SumSquares(:) / (nsubrec - 1))
    end do
    if (allocated(sSet)) deallocate(sSet)
    
    !> Mean of covariances on sub-sub-periods and sub-periods, per variable, Fbar in Eq. 10
    call AverageNoError(Fijs, nrec * nsubrec, GHGNumVar, F_bar, error)

    !> Random errore RE, Eq. 9
    do var = u, gas4
        if (E2Col(var)%present) then
            Essentials%rand_uncer(var) = sum(sigma_wis(:, var)) / nrec / sqrt(float(nsubrec))
        else
            Essentials%rand_uncer(var) = error
        end if
    end do

    !> Sum of squares of residuals, per variable, Eq. 10
    SumSquares = 0d0
    do i = 1, nrec
        SumSquares(:) = SumSquares(:) + (Fi_bars(i, :) - F_bar(:))**2
    end do

    !> Between-records standard deviation, sig_btw in Eq. 10
    sigma_btw = dsqrt(SumSquares(:) / (nrec - 1))

    !> Non stationarity ratio
    where (E2Col(u:GHGNumVar)%present)
        Essentials%mahrt98_NR(:) = sigma_btw(:) / Essentials%rand_uncer(u:gas4)
    else where
        Essentials%mahrt98_NR(:) = error
    end where
end subroutine RU_Mahrt_98    

!***************************************************************************
!
! \brief       Estimate random instrument noise (RIN) according to
!              Billesbach (2011), Eq. 3
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Under development
!***************************************************************************
subroutine RIN_Billesbach_11(Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(inout) :: Set(N, M)
    !> local variables
    integer :: i
    integer :: ntimes = 1
    real(kind = dbl) :: tmpW(N)

    !> Calculate ntimes (given) times the relevant covariances
    !> with 1 time series shuffled (w, so that shuffling is done only once)
    tmpW = Set(w, 1:N)
    do i = 1, ntimes
        call RandomShuffle(tmpW, Set(w, 1:N), N)
        call CovarianceMatrixNoError(Set, size(Set, 1), size(Set, 2), &
            Stats%Cov, error)
    end do
    !> Reset w to its original shape
    Set(w, 1:N) = tmpW
end subroutine RIN_Billesbach_11
!***************************************************************************
!
! \brief       shuffle array elements randomly
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Under development
!***************************************************************************
subroutine RandomShuffle(arr, arrout, N)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: arr(N)
    real(kind = dbl), intent(out) :: arrout(N)
    !> Local variables
    integer :: work
    integer :: ix(N)
    integer :: i
    integer :: j
    integer, external :: RandomBetween

    !> Create array of indexes from 1 to size(arr)
    do i = 1, size(arr)
        ix(i) = i
    end do

    !> Shuffle indexes
    do i = N, 2, -1
        j = RandomBetween(1, i)
        !> swap
        work = ix(j)
        ix(j) = ix(i)
        ix(i) = work
    end do
    !> Assign to shuffled array
    arrout = arr(ix)
end subroutine RandomShuffle

!***************************************************************************
!
! \brief       Generate a random number between min and max
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Under development
!***************************************************************************
integer function RandomBetween(min, max)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: min
    integer, intent(in) :: max
    real(kind = dbl) :: x

    call random_number(x)
    RandomBetween =int((max - min) * x + min)
end function RandomBetween

!***************************************************************************
!
! \brief       Initialize random number generation
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Under development
!***************************************************************************
subroutine InitRandomSeed()
    use m_rp_global_var
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, dt(8), pid, t(2), s
    integer(8) :: count, tms


    call random_seed(size = n)
    allocate(seed(n))

    !> XOR:ing the current time and pid. The PID is
    !> useful in case one launches multiple instances of the same
    !> program in parallel.
    call system_clock(count)
    if (count /= 0) then
        t = transfer(count, t)
    else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
        t = transfer(tms, t)
    end if
    s = ieor(t(1), t(2))
    pid = getpid() + 1099279 ! Add a prime
    s = ieor(s, pid)
    if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
    else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
    end if
    call random_seed(put=seed)
end subroutine InitRandomSeed
