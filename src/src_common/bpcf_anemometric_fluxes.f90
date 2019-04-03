!***************************************************************************
! bpcf_anemometric_fluxes.f90
! ---------------------------
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
