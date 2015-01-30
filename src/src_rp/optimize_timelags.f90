!***************************************************************************
! optimize_timelags.f90
! ---------------------
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
! \brief       Calculate most likely time lag and range of variation for
!              all gases
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OptimizeTimelags(toSet, nrow, actn, M, h2o_n, MM, cls_size)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: M
    integer, intent(in) :: MM
    real(kind = dbl), intent(in) :: cls_size
    type (TimeLagDatasetType), intent(in) :: toSet(nrow)
    integer, intent(in) :: actn(M)
    integer, intent(out) :: h2o_n(MM)
    !> local variables
    integer :: gas
    integer :: cls
    integer :: i
    integer :: N
    integer :: nn
    integer :: first
    integer :: last
    integer :: nup, ndw
    integer, parameter :: min_numerosity = 15
    real(kind = dbl) :: medx
    real(kind = dbl) :: medup, meddw
    real(kind = dbl), allocatable :: tmpx(:)
    real(kind = dbl), allocatable :: devx(:)
    real(kind = dbl), allocatable :: tmpup(:),tmpdw(:)
    real(kind = dbl), allocatable :: devup(:),devdw(:)
    real(kind = dbl) :: MAD
    real(kind = dbl) :: MADup, MADdw
    real(kind = dbl) :: mvec
    real(kind = dbl) :: sdvec
    real(kind = dbl) :: tmpvec(nrow)
    real(kind = dbl) ,parameter :: min_range = 0.3d0

!TO REFINE integer :: read_status
!TO REFINE integer :: h2on

    E2Col(h2o)%present = .true.
    do gas = co2, gas4
        if (E2Col(gas)%present) then
            !> All gases, including H2O, are treated here
            toPasGas(gas)%def = error
            toPasGas(gas)%min = error
            toPasGas(gas)%max = error
            N = actn(gas)
            allocate (tmpx(N), devx(N))
            tmpx(1:N) = toSet(1:N)%tlag(gas)
            call median(tmpx, N, medx)
            devx(1:N) = dabs(toSet(1:N)%tlag(gas) - medx)
            call median(devx, N, MAD)
            if (MAD < 0.1 ) MAD = 0.1 !< Set a minimum value for MAD
            toPasGas(gas)%def = medx
            toPasGas(gas)%max = medx + (TOSetup%pg_range * MAD / 0.6745d0)
            toPasGas(gas)%min = medx - (TOSetup%pg_range * MAD / 0.6745d0)
            deallocate (tmpx, devx)

            !> If H2O was split in classes, now make H2O calculations
            if (gas == h2o .and. MM > 1) then
                !> Water vapour, the same as above, but for RH classes
                toH2O%def=error
                toH2O%min=error
                toH2O%max=error
                do cls = 1, MM
                    h2o_n(cls) = 0
                    tmpvec = 0d0
                    do i = 1, actn(gas)
                        if(toSet(i)%RH >= dfloat(cls - 1) * cls_size &
                            .and. toSet(i)%RH <= dfloat(cls) * cls_size) then
                            h2o_n(cls) = h2o_n(cls) + 1
                            tmpvec(h2o_n(cls)) = toSet(i)%tlag(h2o)
                        end if
                    end do
                    N = h2o_n(cls)
                    if (N < min_numerosity) cycle
                    !> Eliminate outliers in each class
                    !> and redefine "short" time lag set (without outliers)
                    mvec = sum(tmpvec(:)) / N
                    sdvec = dsqrt(sum((tmpvec(:) - mvec)**2))
                    nn = 0
                    do i = 1, N
                        if (tmpvec(i) > mvec + 5.0d0 * sdvec &
                            .or. tmpvec(i) < mvec - 5.0d0 * sdvec) cycle
                        nn = nn + 1
                        tmpvec(nn) = tmpvec(i)
                    end do

                    !> Calculate plausibility range
                    N = nn
                    allocate(tmpx(N), devx(N))
                    allocate(tmpup(N), tmpdw(N))
                    allocate(devup(N), devdw(N))
                    tmpx(1:N) = tmpvec(1:N)
                    call median(tmpx, N, medx)
                    devx(1:N) = tmpvec(1:N) - medx
                    nup = 0
                    ndw = 0
                    do i = 1, N
                        if (devx(i) > 0d0) then
                            nup = nup + 1
                            tmpup(nup) = devx(i)
                        elseif (devx(i) < 0d0) then
                            ndw = ndw + 1
                            tmpdw(ndw) = dabs(devx(i))
                        end if
                    end do
                    call median(tmpup(1:nup), nup, medup)
                    call median(tmpdw(1:ndw), ndw, meddw)
                    devup(1:nup) = dabs(tmpup(1:nup) - medup)
                    devdw(1:ndw) = dabs(tmpdw(1:ndw) - meddw)
                    call median(devup(1:nup), nup, MADup)
                    call median(devdw(1:ndw), ndw, MADdw)
                    toH2O(cls)%def = medx
                    toH2O(cls)%max = medx + (TOSetup%pg_range * MADup / 0.6745d0)
                    toH2O(cls)%min = medx - (TOSetup%pg_range * MADdw / 0.6745d0)
                    if (TOSetup%pg_range * MADup / 0.6745d0 < min_range) &
                        toH2O(cls)%max = medx + min_range
                    if (TOSetup%pg_range * MADdw / 0.6745d0 < min_range) &
                        toH2O(cls)%min = medx - min_range
                    deallocate (tmpx, devx)
                    deallocate (tmpup, tmpdw)
                    deallocate (devup, devdw)
                end do

                !> Adjust time lags for classes not filled
                !> Detects first good class
                first = 1
                do cls = 2, MM
                    if (h2o_n(cls) > min_numerosity) then
                        first = cls
                        exit
                    end if
                end do
                !> Detects last good class
                last = MM
                do cls = MM - 1, 1, -1
                    if (h2o_n(cls) > min_numerosity) then
                        last = cls
                        exit
                    end if
                end do
                !> Set initial classes equal to first good class
                if (first > 1) toH2O(1:first - 1) = toH2O(first)
                !> For intermediate classes, averages before and after
                if (last > first + 1) then
                    do cls = first + 1, last - 1
                        if (h2o_n(cls) == min_numerosity) then
                            toH2O(cls)%def = (toH2O(cls-1)%def + toH2O(cls+1)%def) * 0.5d0
                            toH2O(cls)%min = (toH2O(cls-1)%min + toH2O(cls+1)%min) * 0.5d0
                            toH2O(cls)%max = (toH2O(cls-1)%max + toH2O(cls+1)%max) * 0.5d0
                        end if
                    end do
                end if
                !> Set traling classes to linear extrapolation of last 2 good ones
                if (last < MM .and. last > 2) then
                    do cls = last + 1, MM
                        toH2O(cls)%def = 2d0 * toH2O(cls - 1)%def - toH2O(cls - 2)%def
                        toH2O(cls)%min = 2d0 * toH2O(cls - 1)%min - toH2O(cls - 2)%min
                        toH2O(cls)%max = 2d0 * toH2O(cls - 1)%max - toH2O(cls - 2)%max
                    end do
                end if
            end if
        end if
    end do

    !> If time lag optimization failed, switch to covariance maximization
    if (toH2O(1)%def == error .and. toH2O(MM)%def == error) then
        call ExceptionHandler(43)
        Meth%tlag = 'maxcov'
    end if
end subroutine OptimizeTimelags
