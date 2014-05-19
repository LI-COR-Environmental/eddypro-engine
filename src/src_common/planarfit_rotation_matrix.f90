!***************************************************************************
! planar_fit_rotation_matrix.f90
! ------------------------------
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
! \brief       Calculate planar fit rotation matrix, to rotate winds into \n
!              the mean-streamline coordinate system
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PlanarFitRotationMatrix(sec, PP)
    use m_common_global_var
    implicit none
    ! in/out variables
    integer, intent(in) :: sec
    real(kind = dbl), intent(out) :: PP(3, 3)
    ! local variables
    real(kind = dbl) :: aux
    real(kind = dbl) :: cos_a
    real(kind = dbl) :: cos_b
    real(kind = dbl) :: sin_a
    real(kind = dbl) :: sin_b


    !---------------------------------------------------------------------------------
    !> \warning: element PFb(1) is coefficient b0, PFb(2) is b1 and PFb(3) is b2
    !>           of the equation w = b0 + u * b1 + v * b2 (Wilczak et al. 2001, eq. 39)
    !---------------------------------------------------------------------------------

    !> Define the last row of PP (See Wilczak et al. 2001, eqs. 42)
    aux = dsqrt(PFb(2, sec)**2 + PFb(3, sec)**2 + 1.d0)
    PP(3, 1) = - PFb(2, sec) / aux
    PP(3, 2) = - PFb(3, sec) / aux
    PP(3, 3) = 1.d0/ aux
    !> Define pitch (a) and roll (b) angles (See Wilczak et al. 2001, eqs. 44)
    aux = dsqrt(PP(3, 2)**2 + PP(3, 3)**2)
    sin_a = PP(3, 1)
    cos_a = aux
    sin_b = - PP(3, 2) / aux
    !sin_b = PFb(3, sec) / dsqrt(PFb(3, sec)**2 + 1d0)  !< equivalent alternative (van Dijk et al. 2004, eq. 3.41)
    cos_b =   PP(3, 3) / aux
    !cos_b = 1d0 / dsqrt(PFb(3, sec)**2 + 1d0)          !< equivalent alternative (van Dijk et al. 2004, eq. 3.41)

    !> Define remaining elements of PP (See Wilczak et al. 2001, eq. 36)
    PP(1, 1) =   cos_a
    PP(1, 2) =   sin_a * sin_b
    PP(1, 3) = - sin_a * cos_b
    PP(2, 1) =   0.d0
    PP(2, 2) =   cos_b
    PP(2, 3) =   sin_b
end subroutine PlanarFitRotationMatrix
