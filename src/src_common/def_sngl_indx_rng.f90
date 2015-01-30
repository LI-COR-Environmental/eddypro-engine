!***************************************************************************
! def_sngl_indx_rng.f90
! ---------------------
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
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefSnglIndxRng(N, RowLags1, Nmin, Nmax, Nred)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: RowLags1
    integer, intent(out) :: Nmin
    integer, intent(out) :: Nmax
    integer, intent(out) :: Nred

    if (RowLags1 >= 0) then
        Nred = N - RowLags1
        Nmin = 1
        Nmax = Nred
    else
        Nred = N - abs(RowLags1)
        Nmin = 1 - RowLags1
        Nmax = N
    end if
end subroutine DefSnglIndxRng
