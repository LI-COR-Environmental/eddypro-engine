!***************************************************************************
! cross_wind_corr.f90
! -------------------
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
! \brief       Calculates cross-wind correction on sonic temperature \n
!              according to Liu et al. (2001, BLM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CrossWindCorr(LocCol, Set, nrow, ncol, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    logical, intent(in) :: printout
    type(ColType), intent(in) :: LocCol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i = 0
    real(kind = dbl) :: A(3)
    real(kind = dbl) :: B(3)
    real(kind = dbl) :: C(3)
    real(kind = dbl) :: Un
    real(kind = dbl) :: Vn
    real(kind = dbl) :: Wn


    if (printout) write(*, '(a)', advance = 'no') '  Cross-wind correction..'

    !> Define correction coefficients, depending on anemometer model
    call CrossWindCoeff(LocCol, A, B, C)
    do i = 1, nrow
        !> Define 3D cross-wind velocities for sonic temperature correction
        if (Set(i, u) /= error .and. Set(i, v) /= error .and. Set(i, ts) /= error) then
            Un = A(1)*Set(i, u)**2 + B(1)*Set(i, v)**2 + C(1)*Set(i, u)*Set(i, v)
            Vn = A(2)*Set(i, u)**2 + B(2)*Set(i, v)**2 + C(2)*Set(i, u)*Set(i, v)
            Wn = A(3)*Set(i, u)**2 + B(3)*Set(i, v)**2 + C(3)*Set(i, u)*Set(i, v)
            !> Temperature correction for cross-wind components
            Set(i, ts) = Set(i, ts) + (Un + Vn + Wn) / 1209.d0
        end if
    end do
    if (printout) write(*,'(a)') ' Done.'
end subroutine CrossWindCorr

!***************************************************************************
!
! \brief       Calculates cross-wind correction on sonic temperature \n
!              according to Liu et al. (2001, BLM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CrossWindCoeff(LocCol, A, B, C)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(in) :: LocCol
    real(kind = dbl), intent(out) :: A(3)
    real(kind = dbl), intent(out) :: B(3)
    real(kind = dbl), intent(out) :: C(3)
    !> local variables
    real(kind = dbl) :: phi


    A = 0.d0
    B = 0.d0
    C = 0.d0
    select case (LocCol%instr%firm(1:len_trim(LocCol%instr%firm)))
    case('gill')
        select case (LocCol%instr%model(1:len_trim(LocCol%instr%model) - 2))
        case('hs_50','hs_100')
            phi  = 48.75d0 * p / 180.d0
            A = (/ (1.d0 - dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2) /)
            B = (/1.d0, (1.d0 - 0.75d0 * dcos(phi)**2), (1.d0 - 0.75d0 * dcos(phi)**2) /)
            C = (/0.d0, (0.5d0 * dsqrt(3.d0) * dcos(phi)**2), (-0.5d0 * dsqrt(3.d0) * dcos(phi)**2) /)
        case('r2')
            A = (/ 0.5d0, 0.d0, 0.d0 /)
            B = (/  1.d0, 0.d0, 0.d0 /)
            C = (/  0.d0, 0.d0, 0.d0 /)
            !> Multiply by 3, because the R2 only uses 1 axis, while the equation in CrossWindCorr assumes average of 3 axis
            A = A * 3d0
            B = B * 3d0
        case('r3_50', 'r3_100', 'wm', 'wmpro')
            phi  = 45.d0 * p / 180.d0
            A = (/ (1.d0 - dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2) /)
            B = (/1.d0, (1.d0 - 0.75d0 * dcos(phi)**2), (1.d0 - 0.75d0 * dcos(phi)**2) /)
            C = (/0.d0, (0.5d0 * dsqrt(3.d0) * dcos(phi)**2), (-0.5d0 * dsqrt(3.d0) * dcos(phi)**2) /)
        case('r3a_100')
            !> TO BE UPDATED, THIS SETTINGS ARE THOSE OF SYMMETRIC GILL'S *****
            phi  = 45.d0 * p / 180.d0
            A = (/ (1.d0 - dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2) /)
            B = (/1.d0, (1.d0 - 0.75d0 * dcos(phi)**2), (1.d0 - 0.75d0 * dcos(phi)**2) /)
            C = (/0.d0, (0.5d0 * dsqrt(3.d0) * dcos(phi)**2), (-0.5d0 * dsqrt(3.d0) * dcos(phi)**2) /)
        end select
    case('metek')
        A = (/  1.d0, (5.d0 / 8.d0), (5.d0 / 8.d0) /)
        B = (/ 0.5d0, (7.d0 / 8.d0), (7.d0 / 8.d0) /)
        C = (/  0.d0, (0.25d0 * dsqrt(3.d0)), (-0.25d0 * dsqrt(3.d0)) /)
    case('csi')
        A = (/ 0.75d0, (15.d0 / 16.d0), (15.d0 / 16.d0) /)
        B = (/   1.d0, (13.d0 / 16.d0), (13.d0 / 16.d0) /)
        C = (/   0.d0, (dsqrt(3.d0) / 8.d0), (-dsqrt(3.d0) / 8.d0) /)
    case('young')
        !> TO BE UPDATED, THIS SETTINGS ARE THOSE OF R3  ****
        phi  = 45.d0 * p / 180.d0
        A = (/ (1.d0 - dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2), (1.d0 - 0.25d0 * dcos(phi)**2) /)
        B = (/1.d0, (1.d0 - 0.75d0 * dcos(phi)**2), (1.d0 - 0.75d0 * dcos(phi)**2) /)
        C = (/0.d0, (0.5d0 * dsqrt(3.d0) * dcos(phi)**2), (-0.5d0 * dsqrt(3.d0) * dcos(phi)**2) /)
    end select
end subroutine CrossWindCoeff
