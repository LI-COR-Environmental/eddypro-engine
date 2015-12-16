
!***************************************************************************
! BPCF_LI7550_analog_filters.f90
! ------------------------------
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
! \brief       Calculate spectral corrections factors related to spectral
!              attenuations due to analog signal filters in the LI-7550.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_LI7550AnalogFilters(measuring_height, displ_height, loc_var_present, &
    wind_speed, zL, ac_frequency, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: ac_frequency
    logical, intent(in) :: printout
    !> local variables
    integer :: i
    integer, parameter :: nseconds = 7200
    integer, parameter :: nfreq = 500
    real(kind = dbl) :: kf(nfreq)
    real(kind = dbl) :: nf(nfreq)
    type(SpectralType) :: Cospectrum(nfreq)
    type(BPTFType) :: BPTF(nfreq)


    if (printout) write(*,'(a)', advance='no') &
        '   Calculating correction for LI-7550 analog signals filtering..'

    !> Log natural frequencies in an artificial freq range
    !> f_min = 1 / 2h --> f_max = 10 Hz
    !> Natural frequency array
    nf(1) = 1d0 / nseconds
    do i = 2, nfreq
        nf(i) = dexp(dlog(nf(1)) + (dlog(10d0)-dlog(nf(1)))/dfloat(nfreq) * dfloat(i))
    end do

    !> Normalized frequency vector, kf
    kf(:) = nf(:) * dabs((measuring_height - displ_height) / wind_speed)

    !> Initialize all transfer functions to 1
    call SetTransferFunctionsToValue(BPTF, nfreq, 1d0)

    !> Analytic co-spectra after Moncrieff et al. (1997, JH)
    call CospectraMoncrieff97(nf, kf, Cospectrum, zL, nfreq)

    !> Analytic low-pass transfer function
    call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), u, ac_frequency, &
        loc_var_present, BPTF)
    call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), w, ac_frequency, &
        loc_var_present, BPTF)
    call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), ts, ac_frequency, &
        loc_var_present, BPTF)
    if (loc_var_present(co2)) &
        call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), co2, ac_frequency, &
            loc_var_present, BPTF)
    if (loc_var_present(h2o)) &
        call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), h2o, ac_frequency, &
            loc_var_present, BPTF)

    !> reset to 1 BA and ZOH low-pass transfer functions if the case
    if (.not. EddyProProj%hf_correct_ghg_ba) then
        do i = 1, nfreq
            BPTF(i)%LP%ba_sonic = 1d0
            BPTF(i)%LP%ba_irga = 1d0
        end do
    end if
    if (.not. EddyProProj%hf_correct_ghg_zoh) then
        do i = 1, nfreq
            BPTF(i)%LP%zoh_sonic = 1d0
        end do
    end if

    !> combined transfer functions (low-pass analytic + high-pass analytic)
    !> combined tf (low-pass analytic + high-pass analytic)
    !> calculate correction factors
    call BandPassTransferFunction(BPTF, w,  u,  w_u, nfreq)
    call BandPassTransferFunction(BPTF, w, ts, w_ts, nfreq)
    if (loc_var_present(co2)) &
        call BandPassTransferFunction(BPTF, w, co2,  w_co2,  nfreq)
    if (loc_var_present(h2o))  &
        call BandPassTransferFunction(BPTF, w, h2o,  w_h2o,  nfreq)
    if (loc_var_present(ch4))  &
        call BandPassTransferFunction(BPTF, w, ch4,  w_ch4,  nfreq)
    if (loc_var_present(gas4)) &
        call BandPassTransferFunction(BPTF, w, gas4, w_gas4, nfreq)

    !> calculate correction factors
    call SpectralCorrectionFactors(Cospectrum%of(w_u),  u,  nf, nfreq, BPTF)
    call SpectralCorrectionFactors(Cospectrum%of(w_ts), ts, nf, nfreq, BPTF)
    if(loc_var_present(co2)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_co2), co2, nf, nfreq, BPTF)
    if(loc_var_present(h2o)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_h2o), h2o, nf, nfreq, BPTF)
    if(loc_var_present(ch4)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_ch4), ch4, nf, nfreq, BPTF)
    if(loc_var_present(gas4)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_gas4), gas4, nf, nfreq, BPTF)

    if (printout) write(*,'(a)') ' Done.'
end subroutine BPCF_LI7550AnalogFilters

