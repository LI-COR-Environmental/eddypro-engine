!***************************************************************************
! test_spike_detection.f90
! ------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
! \brief       Detects  and count spikes, and replace by \n
!              linear interpolation if requested. \n
!              Hard-flags file for too many spikes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestSpikeDetection(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: win_len
    integer, parameter :: step = 100 !window advancement through the file, in samples
    integer :: max_pass = 10
    real(kind = dbl), parameter :: lim_step = 0.1d0 !increase of inliers range
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: nn = 0
    integer :: imin = 0
    integer :: imax = 0
    integer :: cnt = 0
    integer :: passes = 0
    integer :: wdw_num = 0
    integer :: wdw = 0
    integer :: nspikes(E2NumVar)
    integer :: nspikes_sng(E2NumVar)
    integer :: tot_spikes(E2NumVar)
    integer :: tot_spikes_sng(E2NumVar)
    integer :: hflags(u:gas4)
    real(kind = dbl) :: Mean(E2NumVar) = 0.d0
    real(kind = dbl) :: StDev(E2NumVar) = 0.d0
    real(kind = dbl) :: LocMean(N, E2NumVar)
    real(kind = dbl) :: LocStDev(N, E2NumVar)
    logical :: IsSpike(N, E2NumVar)
    real(kind = dbl) :: m = 0.d0
    real(kind = dbl) :: q = 0.d0
    real(kind = dbl) :: adv_lim(E2NumVar)
    real(kind = dbl), allocatable :: XX(:, :)
    logical :: again = .false.
    logical :: new_spike


    if (printout) write(*, '(a)', advance = 'no') '   Spike detection/removal test..'

    if (.not. RPsetup%filter_sr) max_pass = 1

    win_len = RPsetup%avrg_len / 6
    if (win_len == 0) win_len = 1
    nn = idint(dble(win_len) * Metadata%ac_freq * 60.d0) !> win length in samples
    wdw_num = idint(dble(N - nn) / 1d2) + 1              !> number of wins for the current file

    !> initialisations
    allocate(XX(nn, E2NumVar))
    LocMean = 0d0
    LocStDev = 0d0
    passes = 0
    nspikes = 0
    nspikes_sng = 0
    tot_spikes = 0
    tot_spikes_sng = 0

    !> Set different threshold for different variables. Specifically, w, ch4 and gas4 have their own tresholds
    adv_lim(u:pe) = sr%lim_u
    adv_lim(w)    = sr%lim_w
    adv_lim(co2)  = sr%lim_co2
    adv_lim(h2o)  = sr%lim_h2o
    adv_lim(ch4)  = sr%lim_ch4
    adv_lim(gas4)  = sr%lim_gas4

    IsSpike = .false.
    !> main cycle, looping over the moving window
100 continue
    passes = passes + 1
    do wdw = 1, wdw_num
        !> pick up the dataset from Set for the current win
        do i = 1, nn
            where(E2Col(:)%present)
                XX(i, :) = Set(i + step * (wdw - 1), :)
            end where
        end do
        !> define min and max of central points
        imin = nn/2 - step/2 + step * (wdw - 1)
        imax = nn/2 + step/2 -1 + step * (wdw - 1)
        !> mean window values
        Mean = sum(XX, dim = 1)
        Mean = Mean / dble(nn)
        !> window standard deviations
        StDev = 0.d0
        do i = 1, nn
            where(E2Col(:)%present)
                StDev(:) = StDev(:) + (XX(i, :) - Mean(:)) **2
            end where
        end do
        where(E2Col(:)%present)
            StDev(:) = StDev(:) / dble(nn - 1)
            StDev(:) = dsqrt(StDev(:))
        end where
        !> stick window values only to central elements
        do i = imin, imax
            where(E2Col(:)%present)
                LocMean(i, :)  = Mean(:)
                LocStDev(i, :) = StDev(:)
            end where
        end do
        !> special case first window
        if (wdw == 1) then
            do i = 1, imin
                where(E2Col(:)%present)
                    LocMean(i, :)  = Mean(:)
                    LocStDev(i, :) = StDev(:)
                end where
            end do
        end if
        !> special case last window
        if (wdw == wdw_num) then
            do i = imax, N
                where(E2Col(:)%present)
                    LocMean(i, :)  = Mean(:)
                    LocStDev(i, :) = StDev(:)
                end where
            end do
        end if
    end do
    !> spikes detection and removal (if requested) in the whole file
    do j = u, pe
        if (E2Col(j)%present) then
            cnt = 0
            !> special case first record in the file
            if (Set(1, j) > LocMean(1, j) + adv_lim(j) * LocStDev(1, j)) then
                Set(1, j) = LocMean(1, j) + adv_lim(j) * LocStDev(1, j)
                    if (.not. IsSpike (1, j)) then
                        nspikes(j) = nspikes(j) + 1
                        nspikes_sng(j) = nspikes_sng(j) + 1
                        IsSpike (1, j) = .true.
                    end if
            end if
            if (Set(1, j) < LocMean(1, j) - adv_lim(j) * LocStDev(1, j)) then
                Set(1, j) = LocMean(1, j) - adv_lim(j) * LocStDev(1, j)
                    if (.not. IsSpike (1, j)) then
                        nspikes(j) = nspikes(j) + 1
                        nspikes_sng(j) = nspikes_sng(j) + 1
                        IsSpike (1, j) = .true.
                    end if
            end if
            !> Following lines
            if (RPsetup%filter_sr) then
                do i = 2, N
                    if((Set(i, j) > LocMean(i, j) + adv_lim(j) * LocStDev(i, j)) .or. &
                     (Set(i, j) < LocMean(i, j) - adv_lim(j) * LocStDev(i, j))) then
                        cnt = cnt + 1
                    else
                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
                            !> check whether it was a spike already, if not increment the
                            !> number of spikes found
                            new_spike = .true.
                            do k = 1, cnt
                                if (IsSpike(i-k, j)) new_spike = .false.
                                exit
                            enddo
                            if (new_spike) then
                                nspikes(j) = nspikes(j) + 1
                                nspikes_sng(j) = nspikes_sng(j) + cnt
                            end if
                            !> update spike flags
                            do k = 1, cnt
                            IsSpike(i-k, j) = .true.
                            enddo
                            !> replace with linear interpolation if requested
                            m = (Set(i, j) - Set(i - (cnt + 1), j)) / (dble(cnt + 1))
                            q = Set(i - (cnt + 1), j)
                            do k = i - cnt, i - 1
                                Set(k, j) = (m * (dble(k - (i - cnt - 1))) + q)
                            end do
                            cnt = 0

                        else if (cnt > sr%num_spk) then
                            cnt = 0
                        end if
                    end if
                end do
            else
                do i = 2, N
                    if((Set(i, j) > LocMean(i, j) + adv_lim(j) * LocStDev(i, j)) .or. &
                     (Set(i, j) < LocMean(i, j) - adv_lim(j) * LocStDev(i, j))) then
                        cnt = cnt + 1
                    else
                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
                            !> check whether it was a spike already, if not increment the
                            !> number of spikes found
                            new_spike = .true.
                            do k = 1, cnt
                                if (IsSpike(i-k, j)) new_spike = .false.
                                exit
                            enddo
                            if (new_spike) then
                                nspikes(j) = nspikes(j) + 1
                                nspikes_sng(j) = nspikes_sng(j) + cnt
                            end if
                            !> update spike flags
                            do k = 1, cnt
                            IsSpike(i-k, j) = .true.
                            enddo
                            cnt = 0
                        else if (cnt > sr%num_spk) then
                            cnt = 0
                        end if
                    end if
                end do
            end if
        end if
    end do

    !> accumulates the number of spikes
    tot_spikes = tot_spikes + nspikes
    tot_spikes_sng = tot_spikes_sng + nspikes_sng

    !> if spikes have been detected (and removed) does the whole process again
    again = .false.
    do j = u, Pe
        if (E2Col(j)%present) then
            if (nspikes(j) /= 0) again = .true.
            nspikes(j) = 0
            nspikes_sng(j) = 0
        end if
    end do
    if (again .and. (passes < max_pass)) then
        adv_lim(:) = adv_lim(:) + lim_step
        goto 100
    end if
    deallocate(XX)

    !> hflags the variable if nspikes is larger than a prescribed threshold
    !> For flagging, limits attention to variables u to gas4
    hflags = 9
    do j = u, gas4
        if (E2Col(j)%present) then
            if(100.d0 * (dble(tot_spikes(j)) / dble(N)) >= sr%hf_lim) then
                hflags(j) = 1
            else
                hflags(j) = 0
            end if
        end if
    end do

    !> Create an 8-digits number containing the values of the hflags
    IntHF%sr = 900000000
    do j = u, gas4
        IntHF%sr = IntHF%sr + hflags(j) * 10 ** (gas4 - j)
    end do

    !> Write on output variable
    Essentials%e2spikes(u:pe) = tot_spikes(u:pe)
    if (printout) write(*,'(a)') ' done.'
end subroutine TestSpikeDetection
