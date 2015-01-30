!***************************************************************************
! fix_dataset_for_spectra.f90
! ---------------------------
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
! \brief       replace gaps with linear interpolation of neighboring data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixDatasetForSpectra(Set, nrow, ncol, nrow2)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(out) :: nrow2
    real(kind = dbl) :: Set(nrow, ncol)
    !> local variables
    integer :: j
    integer :: tnrow


    !> Check that each variable doesn't have too many error codes
    do j = 1, ncol
        if (count(Set(1:nrow, j) == error) > nrow / 100) SpecCol(j)%present = .false.
    end do

    !> nrow2 is the smallest nrow of all columns
    nrow2 = nrow
    do j = 1, GHGNumVar
        if (SpecCol(j)%present) then
            call ReplaceGapWithLinearInterpolation(Set(1:nrow, j), size(Set, 1), tnrow, error)
            if (tnrow < nrow2) nrow2 = tnrow
        end if
    end do
end subroutine FixDatasetForSpectra
