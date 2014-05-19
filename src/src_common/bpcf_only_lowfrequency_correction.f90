!***************************************************************************
! bpcf_only_lowfrequency_correction.f90
! -------------------------------------
! Copyright (C) 2013-2014, LI-COR Biosciences
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
! \brief       Calculate correction factors only for high-pass filtering.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_OnlyLowFrequencyCorrection(measuring_height, displ_height, loc_var_present, wind_speed, zL, &
    ac_frequency, avrg_length, detrending_time_constant, detrending_method)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(8), intent(in) :: detrending_method
    !> local variables
    integer :: i
    integer :: gas
    integer :: nseconds
    integer :: nfreq
    type(BPTFType), allocatable   :: BPTF(:)
    real(kind = dbl), allocatable :: kf(:)
    real(kind = dbl), allocatable :: nf(:)
    type(SpectralType), allocatable :: Cospectrum(:)

    !> Number of natural frequencies in an artificial freq range
    !> f_min = 1 / 4h --> f_max = 5 Hz
    nseconds = 4 * 60 * 60
    nfreq = (5 * nseconds) / 2
    allocate (BPTF(nfreq), Cospectrum(nfreq), nf(nfreq), kf(nfreq))

    do i = 1, nfreq
        nf(i) = dble(i) / dble(nseconds)
    end do

    !> Initialize all transfer functions to 1
    call SetTransferFunctionsToValue(BPTF, nfreq, 1d0)

    !> normalized frequency vector, bkf
    kf(:) = nf(:) * dabs((measuring_height - displ_height) / wind_speed)

    !> analytical co-spectra after Moncrieff et al. 1997
    call CospectraMoncrieff97(nf, kf, Cospectrum, zL, nfreq)

    !> add analytic high-pass transfer functions
    call AnalyticHighPassTransferFunction(nf, size(nf), w, ac_frequency, avrg_length, &
        detrending_method, detrending_time_constant, BPTF)

    do gas = co2, gas4
        if (loc_var_present(gas)) call AnalyticHighPassTransferFunction(nf, size(nf), gas, ac_frequency, avrg_length, &
            detrending_method, detrending_time_constant, BPTF)
    end do

    !> combined TF (only high-pass analytic here)
    if (loc_var_present(co2))  call BandPassTransferFunction(BPTF, w, co2,  w_co2,  nfreq)
    if (loc_var_present(h2o))  call BandPassTransferFunction(BPTF, w, h2o,  w_h2o,  nfreq)
    if (loc_var_present(ch4))  call BandPassTransferFunction(BPTF, w, ch4,  w_ch4,  nfreq)
    if (loc_var_present(gas4)) call BandPassTransferFunction(BPTF, w, gas4, w_gas4, nfreq)

    !> calculate correction factors
    if(loc_var_present(co2)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_co2),  co2,  nf, nfreq, BPTF)
    if(loc_var_present(h2o)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_h2o),  h2o,  nf, nfreq, BPTF)
    if(loc_var_present(ch4)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_ch4),  ch4,  nf, nfreq, BPTF)
    if(loc_var_present(gas4)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_gas4), gas4, nf, nfreq, BPTF)
    deallocate(nf, kf, BPTF, Cospectrum)
end subroutine BPCF_OnlyLowFrequencyCorrection
