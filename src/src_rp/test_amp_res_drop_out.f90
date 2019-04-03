!***************************************************************************
! test_amp_res_drop_out.f90
! -------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Checks for poor amplitude resolution and for dropouts \n
!              Flags file accordingly
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestAmpResDropOut(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: win_len = 0
    integer :: j = 0
    integer :: nn = 0
    integer :: wdw = 0
    integer :: bin = 0
    integer :: drop = 0
    integer :: wdw_num = 0
    integer :: binmin(GHGNumVar)
    integer :: binmax(GHGNumVar)
    integer :: nmptbins(GHGNumVar)
    integer :: tot_mptbins(GHGNumVar)
    integer :: amp_res(N, GHGNumVar)
    integer :: SumX(GHGNumVar)
    integer :: ndrops(N, GHGNumVar)
    integer :: ndrops_ext(GHGNumVar)
    integer :: ndrops_ctr(GHGNumVar)
    integer :: tot_drops_ext(GHGNumVar)
    integer :: tot_drops_ctr(GHGNumVar)
    integer :: ar_hflags(GHGNumVar)
    integer :: do_hflags(GHGNumVar)
    real(kind = dbl), allocatable :: XX(:, :)
    real(kind = dbl) :: MaxX(GHGNumVar) = 0.d0
    real(kind = dbl) :: MinX(GHGNumVar) = 0.d0
    real(kind = dbl) :: RangeX(GHGNumVar) = 0.d0
    real(kind = dbl) :: Mean(GHGNumVar) = 0.d0
    real(kind = dbl) :: Var(GHGNumVar) = 0.d0
    real(kind = dbl) :: StDev(GHGNumVar) = 0.d0
    real(kind = dbl) :: IntX(GHGNumVar)
    real(kind = dbl) :: MaxIntX(GHGNumVar)
    real(kind = dbl) :: MinIntX(GHGNumVar)


    write(*, '(a)', advance = 'no') '   Amplitude resolution and/or dropouts test..'

    !> Initializations
    win_len = RPsetup%avrg_len / 15
    if (win_len == 0) win_len = 1
    nn = idint(dble(win_len) * Metadata%ac_freq * 60.d0) !win length in samples

    do%extlim_up = 100.d0 - do%extlim_dw
    wdw_num = idint(dble((N - nn) / (nn / 2))) + 1
    tot_drops_ext = 0
    tot_drops_ctr = 0
    tot_mptbins = 0
    allocate(XX(nn, GHGNumVar))
    !> Main cycle, looping over the moving window
    do wdw = 1, wdw_num
        amp_res = 0
        nmptbins = 0
        ndrops = 0
        do i = 1, nn
            XX(i, u:GHGNumVar) = Set(i + idint(dble(nn / 2)) * (wdw - 1), u:GHGNumVar)
        end do
        !> Ranges
        MaxX = maxval(XX, dim = 1)
        MinX = minval(XX, dim = 1)
        RangeX(:) = MaxX(:) - MinX(:)
        !> Mean values
        Mean = sum(XX, dim = 1)
        Mean = Mean / dble(nn)
        !> Standard deviations
        Var = 0.d0
        do i = 1, nn
            Var(:) = Var(:) + (XX(i, :) - Mean(:))**2
        end do
        Var(:) = Var(:) / dble(nn - 1)
        StDev = dsqrt(Var)
        !> Interval amplitudes and limits
        do j = u, GHGNumVar
            if(RangeX(j) < (dble(ar%lim) * StDev(j))) then
                MinIntX(j) = MinX(j)
                MaxIntX(j) = MaxX(j)
                IntX(j) = RangeX(j) / (dble(ar%bins))
            else
                MinIntX(j) = Mean(j) - (dble(ar%lim) / 2.d0) * StDev(j)
                MaxIntX(j) = Mean(j) + (dble(ar%lim) / 2.d0) * StDev(j)
                IntX(j) = dble(ar%lim)*StDev(j) / dble(ar%bins)
            end if
        end do
        !> Count data and drop-outs in the bins
        do j = u, GHGNumVar
            drop = 0
            do i = 1, nn
                put_in_bin: do bin = 1, ar%bins
                    if ((XX(i, j) >= (MinIntX(j) + (dble(bin - 1))*IntX(j))).and. &
                        (XX(i, j) < (MinIntX(j) +  (dble(bin)) * intX(j)))) then
                        if (bin == drop) then
                            ndrops(bin, j) = ndrops(bin, j) + 1
                        end if
                        drop = bin
                        amp_res(bin, j) = amp_res(bin, j) + 1
                        exit put_in_bin
                    end if
                end do put_in_bin
            end do
        end do

        !> Dropouts
        !> Setect extreme bins ( x<10% and x>90% ) and
        !> central bins (10%<x<90%); here, 10 = do_extlim_dw, 90 = do_extlim_up
        SumX = 0
        binmax(:) = ar%bins - (ar%bins / 10)
        do j = u, GHGNumVar
            ext_bin_detect: do bin = 1, ar%bins
                SumX(j) = SumX(j) + amp_res(bin, j)
                if (dble(SumX(j)) <= (do%extlim_dw*(dble(nn / 100)))) then
                    binmin(j) = bin
                end if
                if (dble(SumX(j)) >= (do%extlim_up*(dble(nn / 100)))) then
                    binmax(j) = bin - 1
                    exit ext_bin_detect
                end if
            end do ext_bin_detect
        end do

        !> Summation of dropout-flagged values, in the extreme and central intervals
        ndrops_ext = 0
        ndrops_ctr = 0
        do j = u, GHGNumVar
            if(binmin(j) > 0 .and. binmax(j) >= binmin(j)) then
                do bin = 1, binmin(j)
                    ndrops_ext(j) = ndrops_ext(j) + ndrops(bin, j)
                end do
                do bin = binmax(j), ar%bins
                    ndrops_ext(j) = ndrops_ext(j) + ndrops(bin, j)
                end do
                do bin = binmin(j), binmax(j)
                    ndrops_ctr(j) = ndrops_ctr(j) + ndrops(bin, j)
                end do
            end if
        end do
        !> Accumulate number dropouts
        tot_drops_ext = tot_drops_ext + ndrops_ext
        tot_drops_ctr = tot_drops_ctr + ndrops_ctr

        !> Amplitude resolution
        !> Store information about empty bins
        do j = u, GHGNumVar
            do bin = 1, ar%bins
                if (amp_res(bin, j) == 0) nmptbins(j) = nmptbins(j) + 1
            end do
        end do
        !> Accumulate number of empty bins
        tot_mptbins = tot_mptbins + nmptbins
    end do
    deallocate(XX)

    !> Flag the variable if the num of empty bins is larger than threshold
    ar_hflags = 9
    do j = u, GHGNumVar
        if (E2Col(j)%present) then
            if(1d2 * (dble(tot_mptbins(j))) / dble(ar%bins * wdw_num) >= ar%hf_lim) then
                ar_hflags(j) = 1
            else
                ar_hflags(j) = 0
            end if
            Essentials%ar_s(j) = 1d2 * (dble(tot_mptbins(j))) / dble(ar%bins * wdw_num)
        else
            Essentials%ar_s(j) = error
        end if
    end do
    !> Flags the variable if the num of dropouts is larger than threshold
    do_hflags = 9
    do j = 1, GHGNumVar
        if (E2Col(j)%present) then
            if (100.d0 * dble(idint(dble(tot_drops_ctr(j) / (nn * wdw_num)))) >= do%hf1_lim) then
                do_hflags(j) = 1
            else
                do_hflags(j) = 0
            end if
            if (100.d0 * dble(idint(dble(tot_drops_ext(j) / (nn * wdw_num)))) >= do%hf2_lim) then
                do_hflags(j) = 1
            else
                do_hflags(j) = 0
            end if
            Essentials%do_s_ctr(j) = 100.d0 * dble(idint(dble(tot_drops_ctr(j) / (nn * wdw_num))))
            Essentials%do_s_ext(j) = 100.d0 * dble(idint(dble(tot_drops_ext(j) / (nn * wdw_num))))
        else
            Essentials%do_s_ctr(j) = error
            Essentials%do_s_ext(j) = error
        end if
    end do

    !>  Create 8-digits numbers containing hflag/sflag values
    IntHF%ar = 900000000
    IntHF%do = 900000000
    do j = 1, GHGNumVar
        IntHF%ar = IntHF%ar + ar_hflags(j)*10**(GHGNumVar - j)
        IntHF%do = IntHF%do + do_hflags(j)*10**(GHGNumVar - j)
    end do
    write(*,'(a)') ' Done.'
end subroutine TestAmpResDropOut
