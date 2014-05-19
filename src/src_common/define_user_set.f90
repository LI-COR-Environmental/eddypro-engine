!***************************************************************************
! define_user_set.f90
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
! \brief       Define "UserSet", the pre-defined set of variables  \n
!              needed for any following processing. \n
!              Variables are: u, v, w, ts, co2, h2o, ch4, gas4, tc, tc, \n
!              ti1, ti2, pi, te, pe
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineUserSet(LocCol, Raw, nrow, ncol, UserSet, unrow, uncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: unrow, uncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in) :: Raw(nrow, ncol)
    real(kind = dbl), intent(out) :: UserSet(unrow, uncol)
    !> local variables
    integer :: j
    integer :: jj


    UserCol%var = 'none'
    UserCol%measure_type = 'none'
    UserSet = error
    jj = 0
    do j = 1, ncol
        select case (LocCol(j)%var(1:len_trim(LocCol(j)%var)))
            !> Sonic and irga variables without property "useit"
            case('co2','h2o','ch4','n2o', 'cell_t', 'int_t_1', 'int_t_2', 'int_p', &
                'air_t', 'air_p', 'u','v','w','ts','sos', &
                'flag_1', 'flag_2')
                if (.not. LocCol(j)%useit) then
                    jj = jj + 1
                    if (jj > uncol) exit
                    UserCol(jj) = LocCol(j)
                    UserCol(jj)%present = .true.
                    UserSet(1:unrow, jj) = Raw(1:unrow, j)
                end if
            case default
                !> Variables with a custom label
                if (.not. LocCol(j)%useit) then
                    jj = jj + 1
                    if (jj > uncol) exit
                    UserCol(jj) = LocCol(j)
                    UserCol(jj)%present = .true.
                    UserSet(1:unrow, jj) = Raw(1:unrow, j)
                end if
                !> Special case of 4th gas calibration reference
                if (j == Gas4CalRefCol) UserCol(jj)%var = 'cal-ref'
        end select
    end do
end subroutine DefineUserSet
