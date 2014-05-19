!***************************************************************************
! filter_data_for_diagnostics.f90
! -------------------------------
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
! \brief       Filters dataset according to diagnostic flag values
!              on bad flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Allow user to select a different threshold (now 10%) to decide
!              whether to calculate flux or not or to switch to WPL
!***************************************************************************
subroutine FilterDatasetForDiagnostics(Set, nrow, ncol, DiagSet, dnrow, dncol)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: dnrow, dncol
    real(kind = dbl), intent(in) :: DiagSet(dnrow, dncol)
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> lcoal variables
    integer :: var
!    logical :: mask(nrow)


    do var = co2, pe
        select case(var)
            case (co2:gas4, pi:pe)
                !> For LI-7200 flag is OK if = 1, so checks that all bits are
                !> set to 1, which means 7 integer for 3 bits words and 15 for four bits words
                if (index(E2Col(var)%instr%model, 'li7200') /= 0) then
                    do i = 1, nrow
                        if (DiagSet(i, diag72) /= error .and. &
                            (ibits(nint(DiagSet(i, diag72)),  5, 4) < 15 .or. &
                            ibits(nint(DiagSet(i, diag72)), 12, 1)  == 0)) &
                            Set(i, var) = error
                    end do
                !> For LI-7500 flag is OK if = 1, so checks that all bits are
                !> set to 1, which means 7 integer for three bits words and 15 for four bits words
                elseif (index(E2Col(var)%instr%model, 'li7500') /= 0) then
                    do i = 1, nrow
                        if (DiagSet(i, diag75) /= error .and. &
                            ibits(nint(DiagSet(i, diag75)),  5, 3) < 7) &
                            Set(i, var) = error
                    end do
                !> For LI-7700 flag is OK if = 0, so checks that all bits are
                !> set to 0, which means 0 integer for both 3 and 4 bits words
                elseif (index(E2Col(var)%instr%model, 'li7700') /= 0) then
                    do i = 1, nrow
                        if (DiagSet(i, diag77) /= error .and. &
                            (ibits(nint(DiagSet(i, diag77)), 5, 1) /= 0 .or. &
                            ibits(nint(DiagSet(i, diag77)),  8, 4) /= 0 .or. &
                            ibits(nint(DiagSet(i, diag77)), 13, 3) /= 0)) &
                            Set(i, var) = error
                    end do
                end if
        end select
    end do

!!    > Special case of Tin/Tout for LI-7200
!!    > Count occurrences of bad flags for either temperature reading
!    mask = .false.
!    mask(:) = (DiagSet(:, diag72) /= error .and. ibits(nint(DiagSet(:, diag72)),  10, 1) == 0) .or. &
!              (DiagSet(:, diag72) /= error .and. ibits(nint(DiagSet(:, diag72)),  11, 1) == 0)
!!    > If there are too many bad flags (the 0.1 here should become user selectable
!    if (count(mask) >= nrow * 1d-1) then
!!        > If CO2 and H2O from LI-7200 are mixing ratio, avoid calculating temperature term for WPL
!        if ((index(E2Col(co2)%instr%model, 'li7200') /= 0 .and. E2Col(co2)%measure_type /= 'molar_density') &
!            .or. (index(E2Col(h2o)%instr%model, 'li7200') /= 0 .and.  E2Col(h2o)%measure_type /= 'molar_density')) then
!            E2Col(co2)%present = .false.
!            E2Col(h2o)%present = .false.
!        else
!            E2Col(ti1)%present = .false.
!            E2Col(ti2)%present = .false.
!            E2Col(tc)%present  = .false.
!        end if
!    else
!        where (mask(:))
!            Set(:, tc) = error
!            Set(:, ti1) = error
!            Set(:, ti2) = error
!        end where
!    end if
end subroutine FilterDatasetForDiagnostics
