!***************************************************************************
! test_discontinuities.f90
! ------------------------
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
! \brief       Checks for unphysical discontinuities in the time series \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestDiscontinuities(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: win_len
    integer :: nn = 0
    integer :: wdw_num = 0
    integer :: wdw = 0
    integer :: npoints_par
    integer :: hflags(GHGNumVar)
    integer :: sflags(GHGNumVar)
    real(kind = dbl) :: Mean(GHGNumVar)
    real(kind = dbl) :: Mean_up(GHGNumVar)
    real(kind = dbl) :: Mean_dw(GHGNumVar)
    real(kind = dbl) :: Var(GHGNumVar)
    real(kind = dbl) :: Var_up(GHGNumVar)
    real(kind = dbl) :: Var_dw(GHGNumVar)
    real(kind = dbl) :: HaarAvr(GHGNumVar)
    real(kind = dbl) :: HaarVar(GHGNumVar)
    real(kind = dbl), allocatable :: XX(:, :)
    real(kind = dbl), allocatable :: XX_dw(:, :)
    real(kind = dbl), allocatable :: XX_up(:, :)


    write(*, '(a)', advance = 'no') '   Discontinuities test..'

    !> Additional control parameters
    win_len = RPsetup%avrg_len / 6
    if (win_len == 0) win_len = 1

    !> Initializations
    nn = idint((dble(win_len)) * Metadata%ac_freq * 6d1)
    wdw_num = idint(dble(N - nn) / 1d2) + 1

    hflags = 0
    sflags = 0
    allocate(XX(nn, GHGNumVar))
    allocate(XX_dw(nn/2, GHGNumVar))
    allocate(XX_up(nn/2, GHGNumVar))
    
    do wdw = 1, wdw_num
        npoints_par = 0
        !> Full window
        do i = 1, nn
            XX(i, u:GHGNumVar) = Set(i + 100 * (wdw - 1), u:GHGNumVar)
        end do
        !> Half windows
        do i = 1, nn / 2
            XX_dw(i, :) = XX(i, :)
            XX_up(i, :) = XX(nn/2 + i, :)
        end do
        !> Convert instantaneous molar densities into mole fractions using standard air molar volume
        do i = 1, nn
            if(E2Col(co2)%measure_type == 'molar_density') XX(i, co2) = XX(i, co2) * StdVair * 1d3
            if(E2Col(h2o)%measure_type == 'molar_density') XX(i, h2o) = XX(i, h2o) * StdVair
            if(E2Col(ch4)%measure_type == 'molar_density') XX(i, ch4) = XX(i, ch4) * StdVair * 1d3
            if(E2Col(gas4)%measure_type == 'molar_density') XX(i, gas4) = XX(i, gas4) * StdVair * 1d3
        end do

        !> Whole window mean values
        call AverageNoError(XX, size(XX, 1), size(XX, 2), Mean, error)
        !> Half windows mean values
        call AverageNoError(XX_dw, size(XX_dw, 1), size(XX_dw, 2), Mean_dw, error)
        call AverageNoError(XX_up, size(XX_up, 1), size(XX_up, 2), Mean_up, error)

        !> Whole window variance
        call StDevNoError(XX, size(XX, 1), size(XX, 2), Var, error)
        Var = Var**2
        !> Half windows variances
        call StDevNoError(XX_dw, size(XX_dw, 1), size(XX_dw, 2), Var_dw, error)
        call StDevNoError(XX_up, size(XX_up, 1), size(XX_up, 2), Var_up, error)

        !> Haar functions
        HaarAvr(:) = Mean_dw(:) - Mean_up(:)
        HaarVar(:) = (Var_dw(:) - Var_up(:)) / Var(:)

        !> Hard/soft flags for discontinuities beyond prescribed thresholds
        do j = u, v
            if (HaarAvr(j) > ds%hf_uv)  hflags(j) = 1
            if (HaarAvr(j) > ds%sf_uv)  sflags(j) = 1
            if (HaarVar(j) > ds%hf_var) hflags(j) = 1
            if (HaarVar(j) > ds%sf_var) sflags(j) = 1
        end do
        if (HaarAvr(w) > ds%hf_w)       hflags(w) = 1
        if (HaarAvr(w) > ds%sf_w)       sflags(w) = 1
        if (HaarVar(w) > ds%hf_var)     hflags(w) = 1
        if (HaarVar(w) > ds%sf_var)     sflags(w) = 1
        if (HaarAvr(ts) > ds%hf_t)      hflags(ts) = 1
        if (HaarAvr(ts) > ds%sf_t)      sflags(ts) = 1
        if (HaarVar(ts) > ds%hf_var)    hflags(ts) = 1
        if (HaarVar(ts) > ds%sf_var)    sflags(ts) = 1
        if (HaarAvr(co2) > ds%hf_co2)   hflags(co2) = 1
        if (HaarAvr(co2) > ds%sf_co2)   sflags(co2) = 1
        if (HaarVar(co2) > ds%hf_var)   hflags(co2) = 1
        if (HaarVar(co2) > ds%sf_var)   sflags(co2) = 1
        if (HaarAvr(h2o) > ds%hf_h2o)   hflags(h2o) = 1
        if (HaarAvr(h2o) > ds%sf_h2o)   sflags(h2o) = 1
        if (HaarVar(h2o) > ds%hf_var)   hflags(h2o) = 1
        if (HaarVar(h2o) > ds%sf_var)   sflags(h2o) = 1
        if (HaarAvr(ch4) > ds%hf_ch4)   hflags(ch4) = 1
        if (HaarAvr(ch4) > ds%sf_ch4)   sflags(ch4) = 1
        if (HaarVar(ch4) > ds%hf_var)   hflags(ch4) = 1
        if (HaarVar(ch4) > ds%sf_var)   sflags(ch4) = 1
        if (HaarAvr(gas4) > ds%hf_gas4)   hflags(gas4) = 1
        if (HaarAvr(gas4) > ds%sf_gas4)   sflags(gas4) = 1
        if (HaarVar(gas4) > ds%hf_var)   hflags(gas4) = 1
        if (HaarVar(gas4) > ds%sf_var)   sflags(gas4) = 1

        if((sum(hflags) == GHGNumVar) .and. (sum(sflags) == GHGNumVar)) exit
    end do
    if(allocated(XX)) deallocate(XX)
    if(allocated(XX_dw)) deallocate(XX_dw)
    if(allocated(XX_up)) deallocate(XX_up)

    ! creates a 8-digits number containing - in each digit -
    ! the values of the h/s flags:
    IntHF%ds = 900000000
    IntSF%ds = 900000000
    do j = 1, GHGNumVar
        if (E2Col(j)%present) then
            IntHF%ds = IntHF%ds + hflags(j) * 10 **(GHGNumVar - j)
            IntSF%ds = IntSF%ds + sflags(j) * 10 **(GHGNumVar - j)
        else
            IntHF%ds = IntHF%ds + 9 * 10 **(GHGNumVar - j)
            IntSF%ds = IntSF%ds + 9 * 10 **(GHGNumVar - j)
        end if
    end do
    write(*,'(a)') ' Done.'
end subroutine TestDiscontinuities
