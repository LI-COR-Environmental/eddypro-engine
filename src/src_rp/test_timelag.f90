!***************************************************************************
! test_timelag.f90
! ----------------
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
! \brief       Check for scalar time-lags out of suggested ranges    \n
!              hard-flags and soft-flags file accordingly.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestTimeLag(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: rlag(GHGNumVar)
    integer :: hflags(4)
    integer :: sflags(4)
    integer :: min_rl(GHGNumVar)
    integer :: max_rl(GHGNumVar)
    integer :: def_rl(GHGNumVar)
    real(kind = dbl) :: FirstCol(N)
    real(kind = dbl) :: SecondCol(N)
    real(kind = dbl) :: tlag(GHGNumVar) = 0.d0
    real(kind = dbl) :: DefCov(GHGNumVar)
    real(kind = dbl) :: MaxCov(GHGNumVar)

    write(*, '(a)', advance = 'no') '   Time lag test..'

    !> Initializations
    hflags = 0
    sflags = 0

    !> Define min e max "row-lags" for scalars, using timelags retrieved from metadata file
    do j = co2, GHGNumVar
        min_rl(j) = nint(E2Col(j)%min_tl * Metadata%ac_freq)
        max_rl(j) = nint(E2Col(j)%max_tl * Metadata%ac_freq)
    end do
    !> Default values are taken from EddyPro settings
    def_rl(co2)  = nint(tl%def_co2 * Metadata%ac_freq)
    def_rl(h2o)  = nint(tl%def_h2o * Metadata%ac_freq)
    def_rl(ch4)  = nint(tl%def_ch4 * Metadata%ac_freq)
    def_rl(gas4) = nint(tl%def_n2o * Metadata%ac_freq)

    !> Actual time-lags (tlag), maximum of the cov. (Rmax) \n
    !>  and cov. for default timelag (R0)
    !> Flags if the difference is too high
    do i = co2, GHGNumVar
        if (E2Col(i)%present) then
            FirstCol(:)  = Set(:, w)
            SecondCol(:) = Set(:, i)
            call CovMaxRS(def_rl(i), min_rl(i), max_rl(i), FirstCol, SecondCol &
                , MaxCov(i), DefCov(i), tlag(i), rlag(i), N)
            if((MaxCov(i) - DefCov(i)) * 1d2 / DefCov(i) >= tl%hf_lim) hflags(i-ts) = 1
            if((MaxCov(i) - DefCov(i)) * 1d2 / DefCov(i) >= tl%sf_lim) sflags(i-ts) = 1
        end if
    end do

    !> Creates 4-digits numbers containing the flags; 4 is the number of EddyPro gases
    IntHF%tl = 90000
    IntSF%tl = 90000
    do j = 1, 4
        IntHF%tl = IntHF%tl + hflags(j)*10**(4 - j)
        IntSF%tl = IntSF%tl + sflags(j)*10**(4 - j)
    end do
    write(*,'(a)') ' Done.'
end subroutine TestTimeLag

!***************************************************************************
!
! \brief       Performs covariance analysis for determinig the "optimal" time-lag
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CovMaxRS(lagctr, lagmin, lagmax, Col1, Col2, MaxCov, DefCov, TLag, RLag, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N !< Number of raw data
    integer, intent(in) :: lagmin
    integer, intent(in) :: lagmax
    integer, intent(in) :: lagctr
    real(kind = dbl), intent(in) :: Col1(N) !< Fluctuations
    real(kind = dbl), intent(in) :: Col2(N) !< Fluctuations
    integer, intent(out) :: RLag
    real(kind = dbl), intent(out) :: TLag
    real(kind = dbl), intent(out) :: DefCov
    real(kind = dbl), intent(out) :: MaxCov
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: ii = 0
    integer :: N2 !< Number of rows of the raw dataset after time-lag compensation
    real(kind = dbl), allocatable :: ShLocSet(:, :)
    real(kind = dbl) :: Mean(2)
    real(kind = dbl) :: Cov

    Cov = 0.d0
    MaxCov = 0.d0
    DefCov = 0.d0
    TLag = 0.d0
    do i = lagmin, lagmax
        N2 = N - abs(i)
        allocate(ShLocSet(N2, 2))

        do ii = 1, N2
            if (i < 0) then
                ShLocSet(ii, 1) = Col1(ii - i)
                ShLocSet(ii, 2) = Col2(ii)
            else
                ShLocSet(ii, 1) = Col1(ii)
                ShLocSet(ii, 2) = Col2(ii + i)
            end if
        end do
        Mean = sum(ShLocSet, dim = 1)
        Mean = Mean / dble(N2)
        ShLocSet(:, 1) = ShLocSet(:, 1) - Mean(1)
        ShLocSet(:, 2) = ShLocSet(:, 2) - Mean(2)
        !> covariance
        do j = 1, N2
            Cov = Cov + ShLocSet(j, 1) * ShLocSet(j, 2)
        end do
        Cov = Cov / dble(N2 - 1)
        if (i == lagctr) DefCov = Cov
        !> max cov and actual time-lag
        if (abs(Cov) > abs(MaxCov)) then
            MaxCov = Cov
            TLag = dble(i) / Metadata%ac_freq
            RLag = i
        end if
        deallocate(ShLocSet)
    end do
end subroutine CovMaxRS
