!***************************************************************************
! spectral_density.f90
! --------------------
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
! \brief       Calculate power (co)spectra from normalised Fourier coefficients \n
!              (Co)spectra are intended as one-sided, i.e. the negative frequency range \n
!              is folded into the positive thus spectral densities for 1 < f < +inf \n
!              are multiplied by 2 to account for the respective negative frequencies. \n
!              Normalization is done first with the "window squared and summed", \n
!              as defined in Numerical Recipes in C/Fortran \n
!              eq. 13.4.11. This equals N**2 if no tapering is applied. \n
!              Further normalization is done with "df", such that: \n
!              int_fmin^fmax{S(f)df} = 1.
!              "df" is given by df = (fmin-fmax) / (N/2) ca. = AcFreq/Ns
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OneSidedPowerSpectrum(xx, yy, ac_freq, sumw, co, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: ac_freq
    real(kind = dbl), intent(in) :: xx(N)
    real(kind = dbl), intent(in) :: yy(N)
    real(kind = dbl), intent(in) :: sumw
    real(kind = dbl), intent(out) :: co(N/2+1)
    !> local variables
    integer :: i = 0

    !> Cospectral densities
	!> In general: Co(a,b)=(Re(fft(a))*Re(fft(b))+Im(fft(a))*Im(fft(b))
    co(1) = xx(1) * yy(1)
    do i = 2, N - 2, 2
        co(i/2 + 1) = 2d0 * (xx(i) * yy(i) + xx(i + 1) * yy(i + 1))
    end do
    co(N/2 + 1) = xx(N) * yy(N)
    co(:) = co(:) * (dble(N)/ac_freq) / sumw
end subroutine OneSidedPowerSpectrum


!***************************************************************************
!
! \brief       Raw power spectral density
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine rawPSD(xx, yy, psd, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: xx(N)
    real(kind = dbl), intent(in) :: yy(N)
    real(kind = dbl), intent(out) :: psd(N/2+1)
    !> local variables
    integer :: i

    !> (co)spectral density
	!> Co(a,b)=(Re(fft(a))*Re(fft(b))+Im(fft(a))*Im(fft(b))
    psd(1) = xx(1) * yy(1)
    do i = 2, N - 2, 2
        psd(i/2 + 1) = xx(i) * yy(i) + xx(i + 1) * yy(i + 1)
    end do
    psd(N/2 + 1) = xx(N) * yy(N)
end subroutine rawPSD
