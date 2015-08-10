!***************************************************************************
! eliminate_corrupted_variables.f90
! ---------------------------------
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
! \brief       If a variable as more than 30% values set to error, set it
!              as not present.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine EliminateCorruptedVariables(LocSet, nrow, ncol, skip_period, logout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: LocSet(nrow, ncol)
    logical, intent(in)  :: logout
    logical, intent(out) :: skip_period
    !> local variables
    integer :: i
    logical :: mask(nrow)


    if (logout) write(*,'(a)', advance = 'no') '  Verifying time series integrity..'

    do i = 1, ncol
        mask(:) = Locset(:, i) == error
        if (count(mask) > MaxPeriodNumRecords * RPsetup%max_lack/1d2) E2Col(i) = NullCol
    end do

    skip_period = .false.
    if ((.not. E2Col(u)%present) .or. &
        (.not. E2Col(v)%present) .or. &
        (.not. E2Col(w)%present)) skip_period = .true.

    if (logout) write(*,'(a)') ' Done.'
end subroutine EliminateCorruptedVariables
