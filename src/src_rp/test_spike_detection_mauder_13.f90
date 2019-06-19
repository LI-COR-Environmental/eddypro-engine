!***************************************************************************
! test_spike_detection_mauder_13.f90
! ----------------------------------
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
! \brief       Detects and replace spikes by linear interpolation \n
!              Hard-flags file for too many spikes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestSpikeDetectionMauder13(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: max_pass = 10
    integer :: passes = 0
    integer :: i, ii, lN
    integer :: j
    integer :: k
    integer :: cnt = 0
    integer :: nspikes(E2NumVar)
    integer :: nspikes_sng(E2NumVar)
    integer :: tot_spikes(E2NumVar)
    integer :: tot_spikes_sng(E2NumVar)
    integer :: hflags(u:gas4)
    real(kind = dbl) :: MAD
    real(kind = dbl) :: medx
    real(kind = dbl) :: zlim
    real(kind = dbl) :: m
    real(kind = dbl) :: q
    logical :: again
    logical :: IsSpike(N, E2NumVar)
    logical :: new_spike
    real(kind = dbl), allocatable :: tmpx(:)


    if (printout) write(*, '(a)', advance = 'no') '   Spike detection/removal test..'
    zlim = 7d0
    passes = 0
    nspikes = 0
    nspikes_sng = 0
    tot_spikes = 0
    tot_spikes_sng = 0
    IsSpike = .false.
100 continue
    passes = passes + 1
    do j = u, pe
        if (E2Col(j)%present) then
            cnt = 0
            !> Define array tmpx made of Set(:, j) stripped of potential error values.
            lN = N - count(Set(:, j) == error)
            allocate(tmpx(lN))
            ii = 0
            do i = 1, N
                if (Set(i, j) == error) cycle
                ii = ii + 1
                tmpx(ii) = Set(i, j)
            end do

            !> Calculate median (medx) and MAD
            call median(tmpx, lN, medx)
            tmpx = dabs(tmpx - medx)
            call median(tmpx, lN, MAD)
            deallocate(tmpx)

            !> Spike detection (and replacement if needed)
            if (RPsetup%filter_sr) then
                do i = 2, N-1
                    if (Set(i, j) == error) cycle

                    if(Set(i, j) > medx + (zlim * MAD / 0.6745d0) .or. &
                        Set(i, j) < medx - (zlim * MAD / 0.6745d0)) then
                        cnt = cnt + 1
                    else
!                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
                        if (cnt /= 0) then
                            !> check whether it was a spike already,
                            !> if not increment the number of spikes found
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
                            if (i - (cnt + 1) > 0 .and. i - (cnt + 1) <= N) then
                                m = (Set(i, j) - Set(i - (cnt + 1), j)) / (dble(cnt + 1))
                                q = Set(i - (cnt + 1), j)
                                do k = i - cnt, i - 1
                                    Set(k, j) = (m * (dble(k - (i - cnt - 1))) + q)
                                end do
                            end if
                            cnt = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In Mauder's method there is no concept of max num of outliers to define a spike
! so these two lines need eliminating.
!                        else if (cnt > sr%num_spk) then
!                            cnt = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        end if
                    end if
                end do
            else
                do i = 2, N-1
                    if (Set(i, j) == error) cycle

                    if(Set(i, j) > medx + (zlim * MAD / 0.6745d0) .or. &
                        Set(i, j) < medx - (zlim * MAD / 0.6745d0)) then
                        cnt = cnt + 1
                    else
                        if (cnt /= 0) then
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In Mauder's method there is no concept of max num of outliers to define a spike
! so these two lines need eliminating.
!                        else if (cnt > sr%num_spk) then
!                            cnt = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        goto 100
    end if

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
    if (.not. RPsetup%filter_sr) tot_spikes_sng(u:pe) = 0
    where (E2Col(u:pe)%present) 
        Essentials%e2spikes(u:pe) = tot_spikes(u:pe)
        Essentials%m_despiking(u:pe) = tot_spikes_sng(u:pe)
    elsewhere
        Essentials%e2spikes(u:pe) = ierror
        Essentials%m_despiking(u:pe) = ierror
    endwhere
    
    if (printout) write(*,'(a)') ' Done.'
end subroutine TestSpikeDetectionMauder13
