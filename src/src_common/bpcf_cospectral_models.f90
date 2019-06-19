!***************************************************************************
! bpcf_cospectral_models.f90
! --------------------------
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
! \brief       Calculates theoretical co-spectra, after \n
!              Moncrieff et al. (1997), Eq 13-16
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CospectraMoncrieff97(nf, kf, Cospectrum, zL, N)
    use m_common_global_var
    implicit none
    ! in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: kf(N)
    real(kind = dbl), intent(in) :: nf(N)
    real(kind = dbl), intent(in) :: zL
    type(SpectralType), intent(out)  :: Cospectrum(N)
    ! local variables
    integer :: i = 0
    real(kind = dbl) :: Ac
    real(kind = dbl) :: Au
    real(kind = dbl) :: Bc
    real(kind = dbl) :: Bu


    if(zL > 0.d0) then
        do i = 1, N
            !> 1.1 stable conditions
            !> Scalar fluxes (T, CO2, H2O)
            Ac = 0.284d0 * ((1.d0 + 6.4d0 * zL)**0.75d0)
            Bc = 2.34d0 * (Ac**(-1.1d0))
            Cospectrum(i)%of(w_co2) = kf(i) / (nf(i) * (Ac + Bc * kf(i)**2.1d0))
            ! Reynolds stress, however not corrected so far
            Au = 0.124d0 * ((1.d0 + 7.9d0 * zL)**0.75d0)
            Bu = 2.34d0 * (Au**(-1.1d0))
            Cospectrum(i)%of(w_u) = kf(i) / (nf(i) * (Au + Bu * kf(i)**2.1d0))
        end do
    else
        do i = 1, N
            !> 1.2 unstable conditions
            !> Scalar fluxes (T, CO2, H2O)
            if(kf(i) < 0.54d0) then
                Cospectrum(i)%of(w_co2) = 12.92d0 * kf(i) / (nf(i) * (1.d0 + 26.7d0 * kf(i))**1.375d0)
            else
                Cospectrum(i)%of(w_co2) = 4.378d0 * kf(i) / (nf(i) * (1.d0 + 3.8d0 * kf(i))**2.4d0)
            end if
            ! reynolds stress, however not corrected so far
            if(kf(i) < 0.24d0) then
                Cospectrum(i)%of(w_u) = 20.78d0 * kf(i) / (nf(i) * (1.d0 + 31.0d0 * kf(i))**1.575d0)
            else
                Cospectrum(i)%of(w_u) = 12.66d0 * kf(i) / (nf(i) * (1.d0 + 9.6d0 * kf(i))**2.4d0)
            end if
        end do
    end if
    !> Set cospectral model for temperature, h2o, ch4 and gas4 equal to co2
    Cospectrum(:)%of(w_ts)  = Cospectrum(:)%of(w_co2)
    Cospectrum(:)%of(w_h2o) = Cospectrum(:)%of(w_co2)
    Cospectrum(:)%of(w_ch4) = Cospectrum(:)%of(w_co2)
    Cospectrum(:)%of(w_gas4) = Cospectrum(:)%of(w_co2)
end subroutine CospectraMoncrieff97
