!***************************************************************************
! despike_user_set.f90
! --------------------
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
! \brief       Detects spikes in insensitive (user) variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DespikeUserSet(UserSet, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: UserSet(nrow, ncol)
    !> local variables
    integer, parameter :: win_len = 5 !window (later win) length in minutes
    integer, parameter :: step = 100 !window advancement through the file, in samples
    integer, parameter :: max_pass = 10 !max number of passes through the file
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
    integer :: num_spk = 4
    integer :: nspikes(ncol)
    integer :: nspikes_sng(ncol)
    real(kind = dbl) :: Mean(ncol)
    real(kind = dbl) :: StDev(ncol)
    real(kind = dbl) :: LocMean(nrow, ncol)
    real(kind = dbl) :: LocStDev(nrow, ncol)
    logical :: IsSpike(nrow, ncol)
    real(kind = dbl) :: mm = 0.d0
    real(kind = dbl) :: q = 0.d0
    real(kind = dbl) :: adv_lim
    real(kind = dbl), allocatable :: XX(:, :)
    logical :: again = .false.
    logical :: new_spike


    Mean(:) = 0.d0
    StDev(:) = 0.d0
    write(*, '(a)', advance = 'no') '  Despiking user set..'
    nn = idint(dble(win_len) * Metadata%ac_freq * 60.d0) !win length in samples
    wdw_num = idint(dble(nrow - nn) / 1d2) + 1 !number of wins for the current file

    !> initialisations
    allocate(XX(nn, ncol))
    LocMean = 0d0
    LocStDev = 0d0
    passes = 0
    nspikes = 0
    nspikes_sng = 0
    adv_lim = 3.5d0
    IsSpike = .false.

    !> main cycle, looping over the moving window
100 continue
    passes = passes + 1
    do wdw = 1, wdw_num
        !> pick up the dataset from Set for the current win
        do i = 1, nn
            XX(i, :) = UserSet(i + step * (wdw - 1), :)
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
            StDev(:) = StDev(:) + (XX(i, :) - Mean(:)) **2
        end do
        StDev(:) = StDev(:) / dble(nn - 1)
        StDev = dsqrt(StDev)
        !> stick window values only to central elements
        do i = imin, imax
            LocMean(i, :)  = Mean(:)
            LocStDev(i, :) = StDev(:)
        end do
        !> special case first window
        if (wdw == 1) then
            do i = 1, imin
                LocMean(i, :)  = Mean(:)
                LocStDev(i, :) = StDev(:)
            end do
        end if
        !> special case last window
        if (wdw == wdw_num) then
            do i = imax, nrow
                LocMean(i, :)  = Mean(:)
                LocStDev(i, :) = StDev(:)
            end do
        end if
    end do
    !> spikes detection and removal in the whole file
    do j = 1, ncol
        cnt = 0
        !> special case first record in the file
        if (UserSet(1, j) > LocMean(1, j) + adv_lim * LocStDev(1, j)) then
            UserSet(1, j) = LocMean(1, j) + adv_lim * LocStDev(1, j)
                if (.not. IsSpike (1, j)) then
                    nspikes(j) = nspikes(j) + 1
                    nspikes_sng(j) = nspikes_sng(j) + 1
                    IsSpike (1, j) = .true.
                end if
        end if
        if (UserSet(1, j) < LocMean(1, j) - adv_lim * LocStDev(1, j)) then
            UserSet(1, j) = LocMean(1, j) - adv_lim * LocStDev(1, j)
                if (.not. IsSpike (1, j)) then
                    nspikes(j) = nspikes(j) + 1
                    nspikes_sng(j) = nspikes_sng(j) + 1
                    IsSpike (1, j) = .true.
                end if
        end if

        do i = 2, nrow
            if((UserSet(i, j) > LocMean(i, j) + adv_lim * LocStDev(i, j)) .or. &
             (UserSet(i, j) < LocMean(i, j) - adv_lim * LocStDev(i, j))) then
                cnt = cnt + 1
            else
                if ((cnt /= 0) .and. (cnt <= num_spk)) then
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
                    !> replace with linear interpolation
                    mm = (UserSet(i, j) - UserSet(i - (cnt + 1), j)) / (dble(cnt + 1))
                    q = UserSet(i - (cnt + 1), j)
                    do k = i - cnt, i - 1
                        UserSet(k, j) = (mm * (dble(k - (i - cnt - 1))) + q)
                    end do
                    cnt = 0
                else if (cnt > num_spk) then
                    cnt = 0
                end if
            end if
        end do
    end do

    !> if spikes have been detected (and removed) does the whole process again
    again = .false.
    do j = 1, ncol
        if (nspikes(j) /= 0) again = .true.
        nspikes(j) = 0
        nspikes_sng(j) = 0
    end do
    if (again .and. (passes < max_pass)) then
        adv_lim = adv_lim + lim_step
        goto 100
    end if
    deallocate(XX)
    write(*,'(a)') ' Done.'
end subroutine DespikeUserSet
