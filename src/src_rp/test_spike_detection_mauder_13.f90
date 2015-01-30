!***************************************************************************
! test_spike_detection_mauder_13.f90
! ----------------------------------
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
subroutine TestSpikeDetectionMauder13(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: max_pass = 10
    integer :: passes = 0
    integer :: i
    integer :: j
    integer :: k
    integer :: cnt = 0
    integer :: nspikes(E2NumVar)
    integer :: nspikes_sng(E2NumVar)
    integer :: tot_spikes(E2NumVar)
    integer :: tot_spikes_sng(E2NumVar)
    integer :: hflags(u:gas4)
    real(kind = dbl) :: tmpx(N)
    real(kind = dbl) :: devx(N)
    real(kind = dbl) :: MAD
    real(kind = dbl) :: medx
    real(kind = dbl) :: zlim
    real(kind = dbl) :: m
    real(kind = dbl) :: q
    logical :: again
    logical :: IsSpike(N, E2NumVar)
    logical :: new_spike

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
            tmpx(1:N) = Set(1:N, j)
            call median(tmpx, N, medx)
            devx(1:N) = dabs(Set(1:N, j) - medx)
            call median(devx, N, MAD)

            if (RPsetup%filter_sr) then
                do i = 1, N
                    if(Set(i, j) > medx + (zlim * MAD / 0.6745d0) .or. &
                        Set(i, j) < medx - (zlim * MAD / 0.6745d0)) then
                        cnt = cnt + 1
                    else
                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
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
                        else if (cnt > sr%num_spk) then
                            cnt = 0
                        end if
                    end if
                end do
            else
                do i = 2, N
                    if(Set(i, j) > medx + (zlim * MAD / 0.6745d0) .or. &
                        Set(i, j) < medx - (zlim * MAD / 0.6745d0)) then
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
    Essentials%e2spikes(u:pe) = tot_spikes(u:pe)
end subroutine TestSpikeDetectionMauder13
