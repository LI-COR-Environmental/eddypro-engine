!***************************************************************************
! define_e2_set.f90
! -----------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
subroutine DefineVars(LocCol, ncol, uncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: ncol
    integer, intent(in) :: uncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    !> local variables
    integer :: j, jj
    character(len(LocCol%label)), external :: replace


    E2Col = NullCol
    !> First, master sonic wind components are handled
    do j = 1, ncol
        !> u-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'u' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(u) = LocCol(j)
            E2Col(u)%present = .true.
            cycle
        end if
        !> v-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'v' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(v) = LocCol(j)
            E2Col(v)%present = .true.
            cycle
        end if
        !> w-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'w' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(w) = LocCol(j)
            E2Col(w)%present = .true.
            cycle
        end if
        !> sonic/fast temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'ts' .and. LocCol(j)%useit) then
            E2Col(ts) = LocCol(j)
            E2Col(ts)%present = .true.
            cycle
        end if
        !> speed of sound (already converted into sonic temperature)
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'sos' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(ts) = LocCol(j)
            E2Col(ts)%present = .true.
            cycle
        end if
    end do

    !> Now, all remaining EddyPro standard variables (concentrations, temperatures and pressure)
    do j = 1, ncol
        !> master co2 concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'co2' .and. LocCol(j)%useit) then
            E2Col(co2) = LocCol(j)
            E2Col(co2)%present = .true.
            cycle
        end if
        !> master h2o concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'h2o' .and. LocCol(j)%useit) then
            E2Col(h2o) = LocCol(j)
            E2Col(h2o)%present = .true.
            cycle
        end if
        !> master ch4 concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'ch4' .and. LocCol(j)%useit) then
            E2Col(ch4) = LocCol(j)
            E2Col(ch4)%present = .true.
            cycle
        end if
        !> master 4th gas concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'n2o' .and. LocCol(j)%useit) then
            E2Col(gas4) = LocCol(j)
            E2Col(gas4)%present = .true.
            cycle
        end if
        !> master cell temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'cell_t' .and. LocCol(j)%useit) then
            E2Col(tc) = LocCol(j)
            E2Col(tc)%present = .true.
            cycle
        end if
        !> master internal temperature 1
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_t_1' .and. LocCol(j)%useit) then
            E2Col(ti1) = LocCol(j)
            E2Col(ti1)%present = .true.
            cycle
        end if
        !> master internal temperature 2
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_t_2' .and. LocCol(j)%useit) then
            E2Col(ti2) = LocCol(j)
            E2Col(ti2)%present = .true.
            cycle
        end if
        !> master internal pressure
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_p' .and. LocCol(j)%useit) then
            E2Col(pi) = LocCol(j)
            E2Col(pi)%present = .true.
            cycle
        end if
        !> master air temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'air_t' .and. LocCol(j)%useit) then
            E2Col(te) = LocCol(j)
            E2Col(te)%present = .true.
            cycle
        end if
        !> master air pressure
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'air_p' .and. LocCol(j)%useit) then
            E2Col(pe) = LocCol(j)
            E2Col(pe)%present = .true.
            cycle
        end if
    end do

    UserCol%var = 'none'
    UserCol%measure_type = 'none'
    jj = 0
    do j = 1, ncol
        select case (LocCol(j)%var(1:len_trim(LocCol(j)%var)))
            !> Sonic and irga variables without property "useit"
            case('co2','h2o','ch4','n2o', 'cell_t', 'int_t_1', 'int_t_2', &
                'int_p', 'air_t', 'air_p', 'u','v','w','ts','sos', &
                'flag_1', 'flag_2')
                if (.not. LocCol(j)%useit) then
                    jj = jj + 1
                    if (jj > uncol) exit
                    UserCol(jj) = LocCol(j)
                    UserCol(jj)%present = .true.
                end if
            case default
                !> Variables with a custom label
                if (.not. LocCol(j)%useit) then
                    jj = jj + 1
                    if (jj > uncol) exit
                    UserCol(jj) = LocCol(j)
                    UserCol(jj)%present = .true.
                    !> Replace spaces with underscores
                    UserCol(jj)%label = replace(UserCol(jj)%label, &
                        ' ', '_', len(UserCol(jj)%label))
                end if
                !> Special case of 4th gas calibration reference
                if (j == Gas4CalRefCol) UserCol(jj)%var = 'cal-ref'
        end select
    end do
end subroutine DefineVars