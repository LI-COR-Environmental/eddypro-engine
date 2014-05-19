!***************************************************************************
! cospectral_models_rp.f90
! ------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
    !> As for now, set cospectral model for temperature, h2o, ch4 and gas4 equal to co2
    Cospectrum(:)%of(w_ts)  = Cospectrum(:)%of(w_co2)
    Cospectrum(:)%of(w_h2o) = Cospectrum(:)%of(w_co2)
    Cospectrum(:)%of(w_ch4) = Cospectrum(:)%of(w_co2)
    Cospectrum(:)%of(w_n2o) = Cospectrum(:)%of(w_co2)
end subroutine CospectraMoncrieff97
