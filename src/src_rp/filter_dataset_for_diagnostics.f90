!***************************************************************************
! filter_dataset_for_diagnostics.f90
! ----------------------------------
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
subroutine FilterDatasetForDiagnostics(Set, nrow, ncol, DiagSet, dnrow, dncol, &
    lDiagAnemometer, filter_for_diag_irga)
    use m_rp_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: dnrow, dncol
    real(kind = dbl), intent(in) :: DiagSet(dnrow, dncol)
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    logical, intent(in) :: filter_for_diag_irga
    type(DiagAnemType), intent(in) :: lDiagAnemometer
    !> lcoal variables
    integer :: var
    logical :: mask(nrow)


    Essentials%m_diag_anem = 0
    Essentials%m_diag_irga = 0

    !> Anemometer diagnostics (binary 0/1 flag)
    if (lDiagAnemometer%binary_flag_present) then
        if (index(MasterSonic%model, 'wm') /= 0 &
            .or. index(MasterSonic%model, 'hs') /= 0) then
            !> Special case for some Gill sonics, for which values of 10 (0A) and
            !> 11 (0B) indicate good data
            mask = DiagSet(1:nrow, diagAnem) /= error .and. DiagSet(1:nrow, diagAnem) /= 0 .and. &
                           DiagSet(1:nrow, diagAnem) /= 10 .and. DiagSet(1:nrow, diagAnem) /= 11
            Essentials%m_diag_anem = count(mask)
            where (DiagSet(1:nrow, diagAnem) /= error .and. DiagSet(1:nrow, diagAnem) /= 0 .and. &
                   DiagSet(1:nrow, diagAnem) /= 10 .and. DiagSet(1:nrow, diagAnem) /= 11)
                Set(1:nrow, u) = error
                Set(1:nrow, v) = error
                Set(1:nrow, w) = error
                Set(1:nrow, ts) = error
            endwhere
        else
            mask = DiagSet(1:nrow, diagAnem) /= error .and. DiagSet(1:nrow, diagAnem) /= 0
            Essentials%m_diag_anem = count(mask)
            !> Otherwise filter all data when flag is different from zero.
            where (DiagSet(1:nrow, diagAnem) /= error .and. DiagSet(1:nrow, diagAnem) /= 0)
                Set(1:nrow, u) = error
                Set(1:nrow, v) = error
                Set(1:nrow, w) = error
                Set(1:nrow, ts) = error
            endwhere
        end if
    end if
    !> Anemometer diagnostics (StaA for some Gill models)
    if (lDiagAnemometer%staa_present) then
        mask = DiagSet(1:nrow, diagStaA) == 0
        Essentials%m_diag_anem = Essentials%m_diag_anem + count(mask)
        where (DiagSet(1:nrow, diagStaA) == 0)
            Set(1:nrow, u) = error
            Set(1:nrow, v) = error
            Set(1:nrow, w) = error
            Set(1:nrow, ts) = error
        endwhere
    end if

    !> IRGA diagnostics
    if (filter_for_diag_irga) then
        do var = co2, pe
            select case(var)
                case (co2:gas4, pi:pe)
                    !> For LI-7200 flag is OK if = 1, so checks that all bits are
                    !> set to 1, which means 7 integer for 3 bits words and 15 for four bits words
                    if (index(E2Col(var)%instr%model, 'li7200') /= 0) then
                        do i = 1, nrow
                            if (DiagSet(i, diag72) /= error .and. &
                                (ibits(nint(DiagSet(i, diag72)),  5, 4) < 15 .or. &
                                ibits(nint(DiagSet(i, diag72)), 12, 1)  == 0)) then
                                Set(i, var) = error
                                Essentials%m_diag_irga(var) = Essentials%m_diag_irga(var) + 1
                            end if
                        end do
                    !> For LI-7500 flag is OK if = 1, so checks that all bits are
                    !> set to 1, which means 7 integer for three bits words and 15 for four bits words
                    elseif (index(E2Col(var)%instr%model, 'li7500') /= 0) then
                        do i = 1, nrow
                            if (DiagSet(i, diag75) /= error .and. &
                                ibits(nint(DiagSet(i, diag75)),  5, 3) < 7) then
                                Set(i, var) = error
                                Essentials%m_diag_irga(var) = Essentials%m_diag_irga(var) + 1
                            end if
                        end do
                    !> For LI-7700 flag is OK if = 0, so checks that all bits are
                    !> set to 0, which means 0 integer for both 3 and 4 bits words
                    elseif (index(E2Col(var)%instr%model, 'li7700') /= 0) then
                        do i = 1, nrow
                            if (DiagSet(i, diag77) /= error .and. &
                                (ibits(nint(DiagSet(i, diag77)), 5, 1) /= 0 .or. &
                                ibits(nint(DiagSet(i, diag77)),  8, 4) /= 0 .or. &
                                ibits(nint(DiagSet(i, diag77)), 13, 3) /= 0)) then
                                Set(i, var) = error
                                Essentials%m_diag_irga(var) = Essentials%m_diag_irga(var) + 1
                            end if
                        end do
                    end if
            end select
        end do
    end if
    where (.not. E2Col(co2:gas4)%present)
        Essentials%m_diag_irga(co2:gas4) = ierror
    endwhere

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
