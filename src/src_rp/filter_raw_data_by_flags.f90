!***************************************************************************
! filter_raw_data_by_flags.f90
! ----------------------------
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
! \brief       Eliminate individual raw records corresponding to out-ranged flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilterRawDataByFlags(LocCol, Raw, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: Raw(nrow, ncol)
    !> local variables
    integer :: i
    integer :: j
    logical :: filtered(nrow)

    filtered = .false.
    !> External cycle on all columns
    do j = 1, ncol
        !> Detect if column is a flag
        if (LocCol(j)%flag%col > 0) then
            !> If column is a flag, filters accordingly
            if (LocCol(j)%flag%upper) then
                do i = 1, nrow
                    if ((.not. filtered(i)) .and. Raw(i, j) > LocCol(j)%flag%threshold) then
                        Raw(i, 1:ncol) = error
                        filtered(i) = .true.
                    end if
                end do
            else
                do i = 1, nrow
                    if ((.not. filtered(i)) .and. Raw(i, j) < LocCol(j)%flag%threshold) then
                        Raw(i, 1:ncol) = error
                        filtered(i) = .true.
                    end if
                end do
            end if
        end if
    end do
end subroutine FilterRawDataByFlags
