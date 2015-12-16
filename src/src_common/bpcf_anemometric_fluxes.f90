!***************************************************************************
! bpcf_anemometric_fluxes.f90
! ---------------------------
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
! \brief      Calculates fully analytical correction factor for anemometric fluxes \n
!             Based on Moncrieff et al. 1997, Journal of Hydrology 188-189 pp. 589-611.
! \author     Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_AnemometricFluxes(measuring_height, displ_height, loc_var_present, &
        LocInstr, wind_speed, t_air, zL, ac_frequency, avrg_length, &
        detrending_time_constant, detrending_method, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: t_air
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(2), intent(in) :: detrending_method
    logical, intent(in) :: printout
    !> local variables
    integer :: i
    integer, parameter :: nseconds = 7200
    integer, parameter :: nfreq = 500
    real(kind = dbl) :: kf(nfreq)
    real(kind = dbl) :: nf(nfreq)
    type(SpectralType) :: Cospectrum(nfreq)
    type(BPTFType) :: BPTF(nfreq)


    !> Log natural frequencies in an artificial freq range
    !> f_min = 1 / 2h --> f_max = 10 Hz
    nf(1) = 1d0 / nseconds
    do i = 2, nfreq
        nf(i) = dexp(dlog(nf(1)) + (dlog(10d0)-dlog(nf(1)))/dfloat(nfreq) * dfloat(i))
    end do

    !> normalized frequency vector, bkf
    kf(:) = nf(:) * dabs((measuring_height - displ_height) / wind_speed)

    !> Initialize all transfer functions to 1
    call SetTransferFunctionsToValue(BPTF, nfreq, 1d0)

    if (EddyProProj%lf_meth == 'analytic') then
        !> Add analytic high-pass transfer functions
        if (printout) write(*,'(a)') '   High-pass correction for anemometric &
            &fluxes. Method: Moncrieff et al. (2004)..'
        call AnalyticHighPassTransferFunction(nf, size(nf), u, ac_frequency, &
            avrg_length, detrending_method, detrending_time_constant, BPTF)
        call AnalyticHighPassTransferFunction(nf, size(nf), w, ac_frequency, &
            avrg_length, detrending_method, detrending_time_constant, BPTF)
        call AnalyticHighPassTransferFunction(nf, size(nf), ts, ac_frequency, &
            avrg_length, detrending_method, detrending_time_constant, BPTF)
        if (printout) write(*,'(a)') '   Done.'
    end if

    !> analytical cospectra after Moncrieff et al. (1997, JH)
    call CospectraMoncrieff97(nf, kf, Cospectrum, zL, nfreq)

    if (EddyProProj%hf_meth /= 'none') then
        !> Analytical low-pass transfer function
        if (printout) write(*,'(a)') '   Low-pass correction for anemometric &
            &fluxes. Method: Moncrieff et al. (1997)..'
        call AnalyticLowPassTransferFunction(nf, size(nf),  u, LocInstr, &
            loc_var_present, wind_speed, t_air, BPTF)
        call AnalyticLowPassTransferFunction(nf, size(nf),  w, LocInstr, &
            loc_var_present, wind_speed, t_air, BPTF)
        call AnalyticLowPassTransferFunction(nf, size(nf), ts, LocInstr, &
            loc_var_present, wind_speed, t_air, BPTF)
        if (printout) write(*,'(a)') '   Done.'
    end if

    !> combined tf (low-pass analytic + high-pass analytic)
    call BandPassTransferFunction(BPTF, w,  u,  w_u, nfreq)
    call BandPassTransferFunction(BPTF, w, ts, w_ts, nfreq)

    !> calculate correction factors
    call SpectralCorrectionFactors(Cospectrum%of(w_u),  u,  nf, nfreq, BPTF)
    call SpectralCorrectionFactors(Cospectrum%of(w_ts), ts, nf, nfreq, BPTF)
end subroutine BPCF_AnemometricFluxes


