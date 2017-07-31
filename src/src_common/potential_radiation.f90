!***************************************************************************
! potential_radiation.f90
! -----------------------
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
! \brief       Computes potential radiation, needed for identification \n
!              of night and daytime
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
function PotentialRadiation(latit) result(RP)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: latit
    !> local variables
    integer :: i
    integer :: lDOY
    real(kind = dbl) :: lat_rad
    real(kind = dbl) :: rpot
    real(kind = dbl) :: theta       !< Some sort of solar angular position
    real(kind = dbl) :: theta_rad       !< Some sort of solar angular position
    real(kind = dbl) :: delta       !< solar declination (in radiants)
    real(kind = dbl) :: ET          !< Equation of Time
    real(kind = dbl) :: omega       !< solar hour angle
    real(kind = dbl) :: LAS         !< boh
    real(kind = dbl), parameter :: solar_constant = 1376d0
    real(kind = dbl) :: RP(17568)

    !> Compute potential radiation, based on lat/long, on a 30 min basis
    do i = 1, 17568
        !> Compute day
        lDOY = (i / 48) + 1
        !> Compute theta
        theta = 2d0 * p * (lDOY - 1d0) / 365d0;
        !> Compute equation of time
        ET = 0.000075d0 + 0.001868d0 * dcos(theta) - 0.032077d0 * dsin(theta) &
          - 0.014615d0 * dcos(2d0 * theta) - 0.040849d0 * dsin(2d0 * theta)
        !> Compute solar declination in radiants
        delta = 0.006918d0 - 0.399912d0 * dcos(theta) + 0.070257d0 * dsin(theta) &
             - 0.006758d0 * dcos(2d0 * theta) + 0.000907d0 * dsin(2d0 * theta) &
             -0.002697d0 * dcos(3d0 * theta) + 0.00148d0 * dsin(3d0 * theta)
        !> Compute LAS
        LAS = 12.d0 - (mod(i, 48) / 2.d0 + 0.25d0)
        LAS = dabs(LAS)
        !> Compute omega
        omega = -15d0 * LAS
        !> Latitude and longitude in radiants
        lat_rad = latit * p / 180d0;
        !> Compute theta_rad
        theta_rad = dacos(dsin(delta) * dsin(lat_rad) &
                 + dcos(delta) * dcos(lat_rad) * dcos(omega * p / 180d0))
        !> Compute extraterrestrial solar radiation
        rpot = solar_constant * (1.00011d0 + 0.034221d0 * dcos(theta) + 0.00128d0 * dsin(theta) &
             + 0.000719d0 * dcos(2d0 * theta) + 0.000077d0 * dsin(2d0 * theta))
        RP(i) = rpot*dcos(theta_rad)
        if (RP(i) < 0d0) RP(i) = 0d0
    end do
end function PotentialRadiation
