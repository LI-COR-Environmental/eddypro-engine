!***************************************************************************
! random_uncertainty_handle.f90
! -----------------------------
! Copyright (C) 2012-2014, LI-COR Biosciences
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

    if (RUsetup%meth == 'none') then
        Essentials%rand_uncer = aflx_error
        Essentials%rand_uncer_LE = aflx_error
        return
    end if

    write(*, '(a)') '  Estimating random uncertainty..'

    !> Calculate Integral turbulence scale
    call IntegralTurbulenceScale(Set, size(Set, 1), size(Set, 2))

    !> Calculate random uncertainty
    Essentials%rand_uncer(u:gas4) = error
    Essentials%rand_uncer_LE = error
    select case (RUsetup%meth)
        case('finkelstein_sims_01')
            call RU_Finkelstein_Sims_01(Set, nrow, ncol)
        case('mann_lenschow_94')
            call RU_Mann_Lenschow_04(nrow)
        case('tbd')
            !call RE_Lenschow(Set, nrow, ncol)
        case default
            call ExceptionHandler(42)
            Essentials%rand_uncer(u:gas4) = aflx_error
            Essentials%rand_uncer_LE = aflx_error
            return
    end select
    write(*, '(a)') '  done.'
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


    !> Define mm based on ITS
    LagMax(u:gas4) = nint(ITS(u:gas4) * Metadata%ac_freq)
    where (LagMax < 0) LagMax = nint(error)
    do var = u, gas4
        if (var == v .or. var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error .and. LagMax(var) /= nint(error)) then
            allocate (gam(0:LagMax(var), 2, 2))
            gam = 0d0
            do lag = 0, LagMax(var)
                gam(lag, 1, 1) = LaggedCovarianceNoError(Set(:, w), Set(:, w), size(Set, 1), lag, error)
                gam(lag, 2, 2) = LaggedCovarianceNoError(Set(:, var), Set(:, var), size(Set, 1), lag, error)
                gam(lag, 1, 2) = LaggedCovarianceNoError(Set(:, w), Set(:, var), size(Set, 1), lag, error)
                gam(lag, 2, 1) = LaggedCovarianceNoError(Set(:, var), Set(:, w), size(Set, 1), lag, error)
            end do

            !> variance of covariances, Eq. 8  in Finkelstein & Sims (2001, JGR)
            !> Initialize the value for h = 0
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
    real(kind=dbl) :: corr_coeff(E2NumVar)

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
    real(kind=dbl), intent(in) :: arr(N)
    real(kind=dbl), intent(out) :: arrout(N)
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
    real(kind=dbl) :: x

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
