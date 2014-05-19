!***************************************************************************
! calibrate_gas4.f90
! ------------------
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
! \brief       Calibrate 4th gas if a cal-ref column is available.
!              Note that so far the calibration procedure is fully customized
!              on the needs of a specific O3 analyzer
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CalibrateGas4(Set, nrow, ncol)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: j

    do j = 1, NumUserVar
        if (UserCol(j)%var == 'cal-ref') then
            !> Converts mV in ppb and then to ppm (with "/ 1d3")
            Set(:, gas4) = Set(:, gas4) * UserStats%Mean(j) / Stats%Mean(gas4) / 1d3
        end if
    end do
end subroutine CalibrateGas4
