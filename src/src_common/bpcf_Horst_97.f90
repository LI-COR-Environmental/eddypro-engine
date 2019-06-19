!***************************************************************************
! bpcf_Horst_97.f90
! -----------------
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
! \brief       Calculates spectral correction factors based on the procedure \n
!              described in Horst, 1997, Boundary-Layer Meteorology 82: 219-233, 1997. \n
!              Based on  analytical cospectra and \n
!              analytical transfer functions, parameterized in-situ.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_Horst97(measuring_height, displ_height, loc_var_present, wind_speed, zL, &
    ac_frequency, avrg_length, detrending_time_constant, detrending_method, lEx, LocSetup)
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
    character(2), intent(in) :: detrending_method
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    integer :: gas
    integer, parameter :: nseconds = 7200
    integer, parameter :: nfreq = 500
    real(kind = dbl) :: kf(nfreq)
    real(kind = dbl) :: nf(nfreq)
    type(SpectralType) :: Cospectrum(nfreq)
    type(BPTFType) :: BPTF(nfreq)
    type(SpectralType) :: BPCF_Horst
    !> local variables
    integer :: i

    !> add analytic high-pass transfer functions, if requested
    if (EddyProProj%lf_meth == 'analytic') then
        !> Log natural frequencies in an artificial freq range
        !> f_min = 1 / 2h --> f_max = 10 Hz
        !> Natural frequency array
        nf(1) = 1d0 / nseconds
        do i = 2, nfreq
            nf(i) = dexp(dlog(nf(1)) + (dlog(10d0)-dlog(nf(1)))/dfloat(nfreq) * dfloat(i))
        end do

        !> Initialize all transfer functions to 1
        call SetTransferFunctionsToValue(BPTF, nfreq, 1d0)

        !> Calculate analytic high-pass transfer functions
        call AnalyticHighPassTransferFunction(nf, size(nf), w, ac_frequency, avrg_length, &
            detrending_method, detrending_time_constant, BPTF)

        do gas = co2, gas4
            if (loc_var_present(gas)) call AnalyticHighPassTransferFunction(nf, size(nf), gas, ac_frequency, avrg_length, &
                detrending_method, detrending_time_constant, BPTF)
        end do

        !> normalized frequency vector, kf
        kf(:) = nf(:) * dabs((measuring_height - displ_height) / wind_speed)

        !> analytical cospectra after Moncrieff et al. (1997, JH)
        call CospectraMoncrieff97(nf, kf, Cospectrum, zL, nfreq)

        !> combined tf (only high-pass analytic)
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
    end if

    !> file-specific cut-off frequencies
    call RetrieveLPTFpars(lEx, 'iir', LocSetup)

    !> calculate correction factors
    BPCF_Horst%of(w_co2:w_gas4) = 1d0
    call CorrectionFactorsHorst97(BPCF_Horst, lEx)

    !> Combine low freq (analytic) + high freq (in situ) correction factors
    where (BPCF%of(w_co2:w_gas4) /= error .and. BPCF_Horst%of(w_co2:w_gas4) /= error)
        BPCF%of(w_co2:w_gas4) = BPCF%of(w_co2:w_gas4) * BPCF_Horst%of(w_co2:w_gas4)
    elsewhere (BPCF_Horst%of(w_co2:w_gas4) /= error)
        BPCF%of(w_co2:w_gas4) = BPCF_Horst%of(w_co2:w_gas4)
    end where
end subroutine BPCF_Horst97
