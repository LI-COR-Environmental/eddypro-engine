!***************************************************************************
! bpcf_moncrieff_97.f90
! ---------------------
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
! \brief       Calculate fully analytical correction factor, based on theoretical \n
!              shapes of co-spectra and of system transfer function. See \n
!              Moncrieff et al. 1997, J. of Hydrology.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_Moncrieff97(measuring_height, displ_height, loc_var_present, &
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
    integer :: gas
    integer, parameter :: nseconds = 7200
    integer, parameter :: nfreq = 500
    real(kind = dbl) :: kf(nfreq)
    real(kind = dbl) :: nf(nfreq)
    type(SpectralType) :: Cospectrum(nfreq)
    type(BPTFType) :: BPTF(nfreq)


    !> Log natural frequencies in an artificial freq range
    !> f_min = 1 / 2h --> f_max = 10 Hz
    !> Natural frequency array
    nf(1) = 1d0 / nseconds
    do i = 2, nfreq
        nf(i) = dexp(dlog(nf(1)) + (dlog(10d0)-dlog(nf(1)))/dfloat(nfreq) * dfloat(i))
    end do

    !> Initialize all transfer functions to 1
    call SetTransferFunctionsToValue(BPTF, nfreq, 1d0)

    if (EddyProProj%lf_meth == 'analytic') then
        !> Add analytic high-pass transfer functions
        if (printout) write(*,'(a)') '   High-pass correction for gas fluxes. &
            &Method: Moncrieff et al. (2004)'
        call AnalyticHighPassTransferFunction(nf, size(nf), w, ac_frequency, &
            avrg_length, detrending_method, detrending_time_constant, BPTF)

        do gas = co2, gas4
            if (loc_var_present(gas)) &
                call AnalyticHighPassTransferFunction(nf, size(nf), gas, &
                ac_frequency, avrg_length, detrending_method, &
                detrending_time_constant, BPTF)
        end do
        if (printout) write(*,'(a)') '   Done.'
    end if

    if (printout) write(*,'(a)') '   Low-pass correction for gas fluxes. &
        &Method: Moncrieff et al. (1997)'

    !> normalized frequency vector, kf
    kf(:) = nf(:) * dabs((measuring_height - displ_height) / wind_speed)

    !> analytical co-spectra after Moncrieff et al. (1997, JH)
    call CospectraMoncrieff97(nf, kf, Cospectrum, zL, nfreq)

    !> analytic low-pass transfer function
    call AnalyticLowPassTransferFunction(nf, size(nf), w, LocInstr, &
        loc_var_present, wind_speed, t_air, BPTF)
    do gas = co2, gas4
        if(loc_var_present(gas)) call AnalyticLowPassTransferFunction(nf, size(nf),  &
            gas, LocInstr, loc_var_present, wind_speed, t_air, BPTF)
    end do

    !> combined transfer functions (low-pass analytic + high-pass analytic)
    if (loc_var_present(co2))  call BandPassTransferFunction(BPTF, w, co2,  w_co2,  nfreq)
    if (loc_var_present(h2o))  call BandPassTransferFunction(BPTF, w, h2o,  w_h2o,  nfreq)
    if (loc_var_present(ch4))  call BandPassTransferFunction(BPTF, w, ch4,  w_ch4,  nfreq)
    if (loc_var_present(gas4)) call BandPassTransferFunction(BPTF, w, gas4, w_gas4, nfreq)

    !> calculate correction factors
    if(loc_var_present(co2)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_co2), co2, nf, nfreq, BPTF)
    if(loc_var_present(h2o)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_h2o), h2o, nf, nfreq, BPTF)
    if(loc_var_present(ch4)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_ch4), ch4, nf, nfreq, BPTF)
    if(loc_var_present(gas4)) &
        call SpectralCorrectionFactors(Cospectrum%of(w_gas4), gas4, nf, nfreq, BPTF)
    if (printout) write(*,'(a)') '   Done.'
end subroutine BPCF_Moncrieff97
