!***************************************************************************
! def_cross_indx_rng.f90
! ----------------------
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
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefCrossIndxRng(N, RowLags1, RowLags2, Nmin, Nmax, Nred)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: RowLags1
    integer, intent(in) :: RowLags2
    integer, intent(out) :: Nmin
    integer, intent(out) :: Nmax
    integer, intent(out) :: Nred

    if (RowLags1 >= 0 .and. RowLags2 >= 0) then
        Nred = N - max(RowLags1, RowLags2)
        Nmin = 1
        Nmax = Nred
    elseif (RowLags1 < 0 .and. RowLags2 >= 0) then
        Nred = N - (RowLags2 - RowLags1)
        Nmin = 1 - RowLags1
        Nmax = N - RowLags2
    elseif (RowLags1 >= 0 .and. RowLags2 < 0) then
        Nred = N - (RowLags1 - RowLags2)
        Nmin = 1 - RowLags2
        Nmax = N - RowLags1
    elseif (RowLags1 < 0 .and. RowLags2 < 0) then
        Nred = N - abs(min(RowLags1, RowLags2))
        Nmin = 1 + abs(min(RowLags1, RowLags2))
        Nmax = N
    end if
end subroutine DefCrossIndxRng
