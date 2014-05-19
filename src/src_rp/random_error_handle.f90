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
            call RE_Finkelstein_Sims_01(Set, nrow, ncol)
        case('mann_lenschow_94')
            call RE_Mann_Lenschow_04(nrow)
        case('tbd')
            !call RE_Lenschow(Set, nrow, ncol)
        case default
            call ErrorHandle(0, 0, 42)
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
subroutine RE_Finkelstein_Sims_01(Set, N, M)
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
    integer :: nn
    integer :: cnt
    real(kind = dbl), allocatable :: gam(:, :, :)
    real(kind = dbl) :: varcov(M)
    integer :: i


    !> Define mm based on ITS
    LagMax(u:gas4) = nint(ITS(u:gas4) * Metadata%ac_freq)
    where (LagMax < 0) LagMax = nint(error)
    do var = u, gas4
        if (var == v .or. var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error .and. LagMax(var) /= nint(error)) then
            allocate (gam(0:LagMax(var), M, M))
            gam = 0d0
            do lag = 0, LagMax(var)
                nn = N - lag
                cnt = 0
                do i = 1, nn
                    if (Set(i, w) /= error .and. Set(i, var) /= error &
                        .and. Set(i + lag, w) /= error .and. Set(i + lag, var) /= error) then
                        cnt = cnt + 1
                        gam(lag,   w,   w) = gam(lag,   w,   w) + Set(i, w)   * Set(i + lag, w)
                        gam(lag, var, var) = gam(lag, var, var) + Set(i, var) * Set(i + lag, var)
                        gam(lag,   w, var) = gam(lag,   w, var) + Set(i, w)   * Set(i + lag, var)
                        gam(lag, var,   w) = gam(lag, var,   w) + Set(i, var) * Set(i + lag, w)
                    end if
                end do
                if (cnt /= 0) then
                    gam(lag,   w,   w) = gam(lag,   w,   w) / dfloat(cnt)
                    gam(lag, var, var) = gam(lag, var, var) / dfloat(cnt)
                    gam(lag,   w, var) = gam(lag,   w, var) / dfloat(cnt)
                    gam(lag, var,   w) = gam(lag, var,   w) / dfloat(cnt)
                else
                    gam(lag,   w,   w) = error
                    gam(lag, var, var) = error
                    gam(lag,   w, var) = error
                    gam(lag, var,   w) = error
                end if
            end do

            !> variance of covariances, Eq. 8  in Finkelstein & Sims (2001, JGR)
            !> Initialize the value for h = 0
            varcov(var) = 0d0
            if (gam(0, w, w) /= error) &
                varcov(var) = gam(0, w, w) * gam(0, var, var) + gam(0, w, var) * gam(0, var, w)

            !> Now cycle on lag. Do it one sided and multiply by 2 (Eq. 9 and 10)
            cnt = 0
            do lag = 0, LagMax(var)
                if (gam(lag, w, w) /= error) then
                    varcov(var) = varcov(var) + 2d0 * gam(lag, w, w) * gam(lag, var, var) &
                        + 2d0 * gam(lag, w, var) * gam(lag, var, w)
                else
                    cnt = cnt + 1
                end if
            enddo
            deallocate (gam)
            !> Normalization (see Eq. 8)
            varcov(var) = varcov(var) / dfloat(N - cnt)

            !> Random error is the square root of this variance
            if (varcov(var) /= 0) then
                Essentials%rand_uncer(var) = dsqrt(abs(varcov(var)))
            else
                Essentials%rand_uncer(var) = error
            end if
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RE_Finkelstein_Sims_01

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
subroutine RE_Mann_Lenschow_04(N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: var
    real(kind = dbl) :: r2(E2NumVar)

    do var = u, gas4
        if (var == w) cycle
        if (E2Col(var)%present .and. ITS(var) /= error) then
            r2(var) = dabs(Stats%cov(w, var)) &
                / (dsqrt(Stats%cov(w, w)) * dsqrt(Stats%cov(var, var)))
            Essentials%rand_uncer(var) = abs(Stats%cov(w, var)) &
                * dsqrt((1d0 + r2(var)**2) / r2(var)**2) * dsqrt (2d0 * ITS(var) / (N / Metadata%ac_freq))
        else
            Essentials%rand_uncer(var) = error
        end if
    end do
end subroutine RE_Mann_Lenschow_04
