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
! \brief       Define "E2Set", the pre-defined set of variables needed for any following \n
!              processing. Variables are: u, v, w, ts, co2, h2o, ch4, gas4, tc, ti1, ti2, pi, te, pe
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineE2Set(LocCol, Raw, nrow, ncol, E2Set, e2nrow, e2ncol, DiagSet, dnrow, dncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: e2nrow, e2ncol
    integer, intent(in) :: dnrow, dncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in)  :: Raw(nrow, ncol)
    real(kind = dbl), intent(out) :: E2Set(e2nrow, e2ncol)
    real(kind = dbl), intent(out) :: DiagSet(dnrow, dncol)
    !> local variables
    integer :: j


    E2Col = NullCol
    E2Set = error
    DiagSet = error
    !> First, master sonic wind components are handled
    do j = 1, ncol
        !> u-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'u' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(u) = LocCol(j)
            E2Col(u)%present = .true.
            E2Set(1:e2nrow, u) = Raw(1:e2nrow, j)
            cycle
        end if
        !> v-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'v' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(v) = LocCol(j)
            E2Col(v)%present = .true.
            E2Set(1:e2nrow, v) = Raw(1:e2nrow, j)
            cycle
        end if
        !> w-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'w' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(w) = LocCol(j)
            E2Col(w)%present = .true.
            E2Set(1:e2nrow, w) = Raw(1:e2nrow, j)
            cycle
        end if
        !> sonic/fast temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'ts' .and. LocCol(j)%useit) then
            E2Col(ts) = LocCol(j)
            E2Col(ts)%present = .true.
            E2Set(1:e2nrow, ts) = Raw(1:e2nrow, j)
            cycle
        end if
        !> speed of sound (already converted into sonic temperature)
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'sos' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(ts) = LocCol(j)
            E2Col(ts)%present = .true.
            E2Set(1:e2nrow, ts) = Raw(1:e2nrow, j)
            cycle
        end if
    end do

    !> Now, all remaining EddyPro standard variables (concentrations, temperatures and pressure)
    do j = 1, ncol
        !> master co2 concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'co2' .and. LocCol(j)%useit) then
            E2Col(co2) = LocCol(j)
            E2Col(co2)%present = .true.
            E2Set(1:e2nrow, co2) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master h2o concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'h2o' .and. LocCol(j)%useit) then
            E2Col(h2o) = LocCol(j)
            E2Col(h2o)%present = .true.
            E2Set(1:e2nrow, h2o) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master ch4 concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'ch4' .and. LocCol(j)%useit) then
            E2Col(ch4) = LocCol(j)
            E2Col(ch4)%present = .true.
            E2Set(1:e2nrow, ch4) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master 4th gas concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'n2o' .and. LocCol(j)%useit) then
            E2Col(gas4) = LocCol(j)
            E2Col(gas4)%present = .true.
            E2Set(1:e2nrow, gas4) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master cell temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'cell_t' .and. LocCol(j)%useit) then
            E2Col(tc) = LocCol(j)
            E2Col(tc)%present = .true.
            E2Set(1:e2nrow, tc) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master internal temperature 1
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_t_1' .and. LocCol(j)%useit) then
            E2Col(ti1) = LocCol(j)
            E2Col(ti1)%present = .true.
            E2Set(1:e2nrow, ti1) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master internal temperature 2
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_t_2' .and. LocCol(j)%useit) then
            E2Col(ti2) = LocCol(j)
            E2Col(ti2)%present = .true.
            E2Set(1:e2nrow, ti2) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master internal pressure
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_p' .and. LocCol(j)%useit) then
            E2Col(pi) = LocCol(j)
            E2Col(pi)%present = .true.
            E2Set(1:e2nrow, pi) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master air temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'air_t' .and. LocCol(j)%useit) then
            E2Col(te) = LocCol(j)
            E2Col(te)%present = .true.
            E2Set(1:e2nrow, te) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master air pressure
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'air_p' .and. LocCol(j)%useit) then
            E2Col(pe) = LocCol(j)
            E2Col(pe)%present = .true.
            E2Set(1:e2nrow, pe) = Raw(1:e2nrow, j)
            cycle
        end if
    end do

    do j = 1, ncol
        !> Diagnostic flags set
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'diag_72' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diag72) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'diag_75' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diag75) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'diag_77' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diag77) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'anemometer_diagnostic' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diagAnem) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'Gill_StaA' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diagStaA) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'Gill_StaD' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diagStaD) = Raw(1:dnrow, j)
            cycle
        end if
    end do
end subroutine DefineE2Set
