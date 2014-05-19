!***************************************************************************
! integral_turbulence_scale.f90
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
! \brief       Calculate integral turbulence time scale
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine IntegralTurbulenceScale(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Local variables
    integer :: lag
    integer :: var
    integer :: nn
    integer :: LagMax
    real(kind = dbl) :: dt
    real(kind = dbl) :: ITS_bill
    real(kind = dbl) :: cov(ncol, ncol)
    real(kind = dbl), allocatable :: xs(:,:)
    real(kind = dbl), allocatable :: w_cross_corr(:,:)
    logical :: w_cross_corr_failed(ncol)


    write(*, '(a)', advance = 'no') '   Estimating integral turbulence scale..'

    !> Initializations
    LagMax = nint(RUsetup%tlag_max * Metadata%ac_freq)
    allocate(w_cross_corr(0:LagMax, ncol))
    dt = 1d0 / Metadata%ac_freq

    !> Cross-correlation (w and all variables) functions
    do lag = 0, LagMax
        nn = nrow - lag
        allocate(xs(nn, ncol))
        !> Copy dataset, translated for all variables but w
        xs(1:nn, w) = set(1:nn, w)
        do var = u, gas4
            if (var /= w .and. E2Col(var)%present) xs(1:nn, var) = set(1+lag: nn+lag, var)
        end do

        !> Cross-correlation function
        call CovarianceMatrixNoError(xs, size(xs, 1), size(xs, 2), cov, error)
        where (E2Col(u:gas4)%present)
            w_cross_corr(lag, u:gas4) = cov(w, u:gas4)
        end where
        deallocate(xs)
    end do

    !> Normalize cross-correlation function
    w_cross_corr_failed = .false.
    do var = u, gas4
        if (var /= w .and. E2Col(var)%present) then
            if (w_cross_corr(0, var) /= 0d0 .and. w_cross_corr(0, var) /= error) then
                w_cross_corr(0:lagMax, var) = w_cross_corr(0:lagMax, var) / w_cross_corr(0, var)
            else
                w_cross_corr_failed(var) = .true.
            end if
        end if
    end do

    !> Integral turbulence scale function
    ITS = 0d0
    select case(RUsetup%its_meth)
        case('cross_0')
            !> First crossing of y = 0
            do var = u, gas4
                if (var /= w .and. E2Col(var)%present) then
                    if (.not. w_cross_corr_failed(var)) then
                        do lag = 0, LagMax
                            if (w_cross_corr(lag, var) < 0d0) exit
                            ITS(var) = ITS(var) + dt * dabs(w_cross_corr(lag, var))
                        end do
                    else
                        ITS(var) = error
                    end if

                else
                    ITS(var) = error
                end if
            enddo
        case('cross_e')
            !> First crossing of y = 1/e
            do var = u, gas4
                if (var /= w .and. E2Col(var)%present) then
                    if (.not. w_cross_corr_failed(var)) then
                        do lag = 0, LagMax
                            if (w_cross_corr(lag, var) <  1d0/exp(1d0)) exit
                            ITS(var) = ITS(var) + dt * dabs(w_cross_corr(lag, var))
                        end do
                    else
                        ITS(var) = error
                    end if
                else
                    ITS(var) = error
                end if
            end do
        case('full_integral')
            !> Intgrate over the full range of variation of time lag
            do var = u, gas4
                if (var /= w .and. E2Col(var)%present) then
                    if (.not. w_cross_corr_failed(var)) then
                        do lag = 0, LagMax
                            ITS(var) = ITS(var) + dt * dabs(w_cross_corr(lag, var))
                        end do
                    else
                        ITS(var) = error
                    end if
                else
                    ITS(var) = error
                end if
            end do
    end select
    deallocate(w_cross_corr)

    !> Filter reasonable values of ITS, in case anything went wrong.
    !> ITS shouldn't be higher than the integral of "1" over the whole time lag period.
    !> Use a factor of 2 to account for anomalies.

    !> For badly estimated ITS, uses simple formula from Wyngaard (1973), as cited in Billesbach (2012):
    !> "The integraltimescale can be estimated (under neutral stability) as the
    !> instrument height divided by the mean wind speed".
    if (Stats%Mean(u) /= error) then
        ITS_bill = (E2Col(u)%Instr%height - Metadata%d) / Stats%Mean(u)
    else
        ITS_bill = error
    end if

    where (ITS(u:gas4) > 2. * RUsetup%tlag_max .or. ITS(u:gas4) == error)
       ITS(u:gas4) = ITS_bill
    end where
    write(*, '(a)') ' done.'
end subroutine IntegralTurbulenceScale
