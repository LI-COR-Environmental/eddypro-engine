!***************************************************************************
! count_records_and_values.f90
! ----------------------------
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
! \brief       Count either:
!              Available records (any record with at least one non-error value)
!              Available values for a passed variable
!              Available values for a pair of passed variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function CountRecordsAndValues(Set, nrow, ncol, var1, var2)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow, ncol
    integer, optional, intent(in) :: var1, var2
    real(kind = dbl), intent(in) :: Set(nrow, ncol)


    if (.not. present(var1) .and. .not. present(var2)) then
        !> No var passed, count whole records
        CountRecordsAndValues = count(any(Set(:, 1:ncol) /= error, dim = 2))
        return
    else if (present(var1) .and. .not. present(var2)) then
        CountRecordsAndValues = count(Set(:, var1) /= error)
        return
    else
        CountRecordsAndValues = count(Set(:, var1) /= error .and. Set(:, var2) /= error)
        return
    end if

end function CountRecordsAndValues