!***************************************************************************
! gill_wm_w_boost.f90
! -------------------
! Copyright (C) 2016, LI-COR Biosciences
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
! \brief       Fix firmware bug for WindMaster and WindMaster Pro with firmware
!              versions of the form 2329.x.y with x < 700.
! \author      Gerardo Fratini
! \notes
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ApplyGillWmWBoost(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    include '../src_common/interfaces.inc'

    Set(:, w) = asymmetric_linear_transformation(Set(:, w), size(Set(:, w), 1), &
                                                1.166d0, 0d0, 1.289d0, 0d0)
end subroutine ApplyGillWmWBoost


!***************************************************************************
!
! \brief       "Apply" firmware bug for WindMaster and WindMaster Pro with firmware
!              versions other than [2329.x.y with x < 700]. Needed if one
!              wants to apply the AoA correction of Nakai and Shimoyama 2012.
! \author      Gerardo Fratini
! \notes
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RemoveGillWmWBoost(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    include '../src_common/interfaces.inc'

    Set(:, w) = asymmetric_linear_transformation(Set(:, w), size(Set(:, w), 1), &
                                                1d0/1.166d0, 0d0, 1d0/1.289d0, 0d0)
end subroutine RemoveGillWmWBoost
