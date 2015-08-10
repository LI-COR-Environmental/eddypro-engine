!***************************************************************************
! generate_t_cell.f90
! -------------------
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
! \brief       Generate LI-7200 cell temperature timeseries based on availability
!              Tcell = Tcell if Tcell is available
!              Tcell = 0.2 * Tin + 0.8 * Tout if Tin and Tout available
!              Tcell = Tin or Tout, if either is missing
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine GenerateTcell(Set, nrow, ncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables


    if(.not. E2Col(tc)%present) then
        if (E2Col(ti1)%present .and. E2Col(ti2)%present) then
            E2Col(tc) = E2Col(ti1)
            E2Col(tc)%var = 'tc'
            where (Set(1:nrow, ti1) /= error .and. Set(1:nrow, ti2) /= error)
                Set(1:nrow, tc) = Set(1:nrow, ti1) * 0.2d0 + Set(1:nrow, ti2) * 0.8d0
            else where (Set(1:nrow, ti1) /= error)
                Set(1:nrow, tc) = Set(1:nrow, ti1)
            else where (Set(1:nrow, ti2) /= error)
                Set(1:nrow, tc) = Set(1:nrow, ti2)
            else where
                Set(1:nrow, tc) = error
            end where

        elseif(E2Col(ti1)%present) then
            E2Col(tc) = E2Col(ti1)
            E2Col(tc)%var = 'tc'
            Set(1:nrow, tc) = Set(1:nrow, ti1)

        elseif(E2Col(ti2)%present) then
            E2Col(tc) = E2Col(ti2)
            E2Col(tc)%var = 'tc'
            Set(1:nrow, tc) = Set(1:nrow, ti2)
        end if
    end if
    E2Col(ti1)%present = .false.
    E2Col(ti2)%present = .false.
end subroutine GenerateTcell
