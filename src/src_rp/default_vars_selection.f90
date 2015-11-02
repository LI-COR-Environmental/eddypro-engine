!***************************************************************************
! default_vars_selection.f90
! --------------------------
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
! \brief       Automatic selection of variables for fluxes, among those available
!              for use in embeeded mode
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefaultVarsSelection(LocCol)
    use m_rp_global_var
    !> in/out variables
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    !> local variables
    integer :: i
    integer :: nsonics
    integer :: nirgas_co2
    integer :: nirgas_ch4
    integer :: co2_instr_indx
    integer :: ch4_instr_indx
    integer :: co2r_col, co2f_col, co2d_col
    integer :: h2or_col, h2of_col, h2od_col
    integer :: tmp


    !> Define master anemometer
    nsonics   = 0
    nirgas_co2 = 0
    nirgas_ch4 = 0
    do i = 1, NumInstruments
        if (Instr(i)%category == 'sonic') nsonics = nsonics + 1
        if (Instr(i)%category == 'irga') then
            select case(Instr(i)%model(1:len_trim(Instr(i)%model) - 2))
                case('li7500', 'li7500a', 'li7500rs', 'li7200', 'li7200rs')
                nirgas_co2 = nirgas_co2 + 1
                case('li7700')
                nirgas_ch4 = nirgas_ch4 + 1
            end select
        end if
    end do


    co2_instr_indx = 0
    ch4_instr_indx = 0
    do i = 1, NumInstruments

        !> Attach master sonic property to the first anemometer, regardless of which one it is,
        !> Implement this: ..under the only condition that all anemometric variables are available for it.
        if (Instr(i)%category == 'sonic') &
            EddyProProj%master_sonic = Instr(i)%model(1:len_trim(Instr(i)%model))

        !> Define analyser for co2/h2o to be used for fluxes, if any
        select case (nirgas_co2)
            case(1)
                !> if there's only one co2/h2o irga, pick that one for fluxes
                if (index(Instr(i)%model, 'li72') /= 0 .or. index(Instr(i)%model, 'li75') /= 0) co2_instr_indx = i
            case(2:)
                !> if there's more than one co2/h2o irga, selects LI-7200 for fluxes
                if (index(Instr(i)%model, 'li72') /= 0) co2_instr_indx = i
        end select

        !> Define analyser for ch4 to be used for fluxes, if any
        if (nirgas_ch4 /= 0 .and. index(Instr(i)%model, 'li77') /= 0) ch4_instr_indx = i
    end do

    !> If EddyProProj%col(i) = 0, that's an user decision not to use
    !> the variable, so does not change it
    !> If EddyProProj%col(i) > 0, that's an user decision to use
    !> that column, so does not change it
    !> In embedded mode, EddyProProj%col(i) < 0 means that the user is leaving
    !> the decision on which variable to use to EddyPro, so select most
    !> appropriate variable now


    !> co2 and h2o variables
    co2r_col = 0
    co2f_col = 0
    co2d_col = 0
    h2or_col = 0
    h2of_col = 0
    h2od_col = 0
    do i = 1, NumCol
        !> co2
        if (EddyProProj%col(co2) < 0 .and. LocCol(i)%var == 'co2') then
            !> A co2 reading was detected, see if it comes from the right instrument
            if (LocCol(i)%Instr%model == Instr(co2_instr_indx)%model) then
                !> If it's from a LI-7500(A), must be a molar density, otherwise nothing
                if (index(Instr(co2_instr_indx)%model, 'li75') /= 0 .and. LocCol(i)%measure_type == 'molar_density') &
                    EddyProProj%col(co2) = i
                if (index(Instr(co2_instr_indx)%model, 'li72') /= 0 .and. LocCol(i)%measure_type == 'mixing_ratio') &
                    co2r_col = i
                if (index(Instr(co2_instr_indx)%model, 'li72') /= 0 .and. LocCol(i)%measure_type == 'mole_fraction') &
                    co2f_col = i
                if (index(Instr(co2_instr_indx)%model, 'li72') /= 0 .and. LocCol(i)%measure_type == 'molar_density') &
                    co2d_col = i
            end if
        end if

        !> h2o (note: co2_instr_indx is valid also for h2o!)
        if (EddyProProj%col(h2o) < 0 .and. LocCol(i)%var == 'h2o') then
            !> An h2o reading was detected, see if it comes from the right instrument
            if (LocCol(i)%Instr%model == Instr(co2_instr_indx)%model) then
                !> If it's from a LI-7500(A), must be a molar density, otherwise nothing
                if (index(Instr(co2_instr_indx)%model, 'li75') /= 0 .and. LocCol(i)%measure_type == 'molar_density') &
                    EddyProProj%col(h2o) = i
                !> If it's from a LI-7200, order is (1) mixing ratio (2) mole fraction and (3) molar density
                if (index(Instr(co2_instr_indx)%model, 'li72') /= 0 .and. LocCol(i)%measure_type == 'mixing_ratio') &
                    h2or_col = i
                if (index(Instr(co2_instr_indx)%model, 'li72') /= 0 .and. LocCol(i)%measure_type == 'mole_fraction') &
                    h2of_col = i
                if (index(Instr(co2_instr_indx)%model, 'li72') /= 0 .and. LocCol(i)%measure_type == 'molar_density') &
                    h2od_col = i
            end if
        end if

        !> ch4
        if (EddyProProj%col(ch4) < 0 .and. LocCol(i)%var == 'ch4') then
            !> An ch4 reading was detected, see if it comes from the right instrument
            if (LocCol(i)%Instr%model == Instr(ch4_instr_indx)%model) then
                !> If it's from a LI-7500(A), must be a molar density, otherwise nothing
                if (index(Instr(ch4_instr_indx)%model, 'li77') /= 0 .and. LocCol(i)%measure_type == 'molar_density') &
                    EddyProProj%col(ch4) = i
            end if
        end if
    end do

    !> In case of LI-7200, need to select the most appropriate column according to the priority list
    !> Order is (1) mixing ratio (2) mole fraction and (3) molar density
    if (co2d_col /= 0) EddyProProj%col(co2) = co2d_col
    if (co2f_col /= 0) EddyProProj%col(co2) = co2f_col
    if (co2r_col /= 0) EddyProProj%col(co2) = co2r_col
    if (h2od_col /= 0) EddyProProj%col(h2o) = h2od_col
    if (h2of_col /= 0) EddyProProj%col(h2o) = h2of_col
    if (h2or_col /= 0) EddyProProj%col(h2o) = h2or_col

    !> Check if internal temperatures and pressure are available
    do i = 1, NumCol
        !> In/out and average cell temperatures
        if (EddyProProj%col(ti1) < 0 .and. LocCol(i)%var == 'int_t_1') EddyProProj%col(ti1) = i
        if (EddyProProj%col(ti2) < 0 .and. LocCol(i)%var == 'int_t_2') EddyProProj%col(ti2) = i
        if (EddyProProj%col(tc)  < 0 .and. LocCol(i)%var == 'cell_t')  EddyProProj%col(tc)  = i
        !> Cell pressure
        if (EddyProProj%col(pc) < 0  .and. LocCol(i)%var == 'int_p')   EddyProProj%col(pc)  = i
    end do

    !> Ambient temperature. If available, better using it, prioritarily from the LI-7700
    tmp = -1
    do i = 1, NumCol
        if (EddyProProj%col(te) < 0 .and. LocCol(i)%var == 'air_t' .and. index(LocCol(i)%Instr%model, 'li77') /= 0) then
            EddyProProj%col(te) = i
            exit
        end if
        if (EddyProProj%col(te) < 0 .and. LocCol(i)%var == 'air_t') tmp = i
    end do
    if (EddyProProj%col(te) < 0 .and. tmp > 0) EddyProProj%col(te) = tmp

    !> Ambient pressure. If available, better using it, prioritarily from the LI-7700
    tmp = -1
    do i = 1, NumCol
        if (EddyProProj%col(pe) < 0 .and. LocCol(i)%var == 'air_p' .and. index(LocCol(i)%Instr%model, 'li77') /= 0) then
            EddyProProj%col(pe) = i
            exit
        end if
        if (EddyProProj%col(pe) < 0 .and. LocCol(i)%var == 'air_p') tmp = i
    end do
    if (EddyProProj%col(pe) < 0 .and. tmp > 0) EddyProProj%col(pe) = tmp

    !> Diagnostic flags
    do i = 1, NumCol
        if (co2_instr_indx /= 0) then
            if (LocCol(i)%var == 'diag_75' .and. LocCol(i)%Instr%model == Instr(co2_instr_indx)%model) &
                EddyProProj%col(E2NumVar + diag75) = i
            if (LocCol(i)%var == 'diag_72' .and. LocCol(i)%Instr%model == Instr(co2_instr_indx)%model) &
                EddyProProj%col(E2NumVar + diag72) = i
        end if
        if (ch4_instr_indx /= 0) then
            if (LocCol(i)%var == 'diag_77' .and. LocCol(i)%Instr%model == Instr(ch4_instr_indx)%model) &
                EddyProProj%col(E2NumVar + diag77) = i
            end if
    end do
end subroutine DefaultVarsSelection
