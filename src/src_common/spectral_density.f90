!***************************************************************************
! spectral_density.f90
! --------------------
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
