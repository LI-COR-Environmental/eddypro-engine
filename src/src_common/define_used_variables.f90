!***************************************************************************
! define_used_variables.f90
! -------------------------
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
! \brief       Determine whether to use individual data column or not, based either on the
!              fact that it's coming from a master sonic, or on the value assigned
!              to %useit (for variables other than sonics).
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineUsedVariables(LocCol)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    !> local variables
    integer :: i
    logical :: ts_found


    NumUserVar = 0

    !> Associate the "flag" properties to columns selected as such by user
    do i = 1, NumRawFlags
        LocCol(RawFlag(i)%col)%flag = RawFlag(i)
    end do

    !> Associate the "master_sonic" property to the relevant instrument
    !> This automatically associate the master_sonic property to the
    !> relevant sonic variables
    LocCol%Instr%master_sonic = .false.
    do i = 1, NumCol
        if (index(EddyProProj%master_sonic, &
            trim(adjustl(LocCol(i)%Instr%model))) /= 0) then
            LocCol(i)%Instr%master_sonic = .true.
        end if
    end do

    LocCol%useit = .false.
    !> Information in EddyPro project file (user explicitly selects which
    !> variables are to be used)
    where (EddyProProj%Col(co2:E2NumVar) > 0)
        LocCol(EddyProProj%Col(co2:E2NumVar))%useit = .true.
    endwhere

    where (EddyProProj%Col(E2NumVar + diag72 :E2NumVar + diag77) > 0)
        LocCol(EddyProProj%Col(E2NumVar + diag72 :E2NumVar + diag77))%useit = .true.
    endwhere

    !> If gas4 column was selected, change its name to 'n2o', to be treated
    !> as such. The column label still holds the actual variable name
    !> as selected/entered in the Metadat File Editor
    if (EddyProProj%Col(gas4) > 0) LocCol(EddyProProj%Col(gas4))%var = 'n2o'

    !> Diagnostic flags
    NumDiag = 0
    Diag7200%present = .false.
    Diag7500%present = .false.
    Diag7700%present = .false.
    if (EddyProProj%Col(E2NumVar + diag72) > 0) then
        LocCol(EddyProProj%Col(E2NumVar + diag72))%useit = .true.
        NumDiag = NumDiag + 1
        Diag7200%present = .true.
    end if
    if (EddyProProj%Col(E2NumVar + diag75) > 0) then
        LocCol(EddyProProj%Col(E2NumVar + diag75))%useit = .true.
        NumDiag = NumDiag + 1
        Diag7500%present = .true.
    end if
    if (EddyProProj%Col(E2NumVar + diag77) > 0) then
        LocCol(EddyProProj%Col(E2NumVar + diag77))%useit = .true.
        NumDiag = NumDiag + 1
        Diag7700%present = .true.
    end if

    !> Loop on the actual number of columns and determine
    !> whether to use them or not
    Gas4CalRefCol = 0
    do i = 1, NumCol
        !> Variables from the master_sonic are to be used
        if (LocCol(i)%instr%master_sonic) then
            LocCol(i)%useit = .true.
            cycle
        end if
        !> Count users variables, made up of: sonic variables from a non-master
        !> sonic; irga variables without the property "use_it", and
        !> those with a custom label
        select case (LocCol(i)%var(1:len_trim(LocCol(i)%var)))
            case('ignore', 'not_numeric')
                !> Skip flags and columns to be ignored
                continue
            case('co2','h2o','ch4','n2o', 'cell_t','int_t_1', 'int_t_2', &
                'int_p', 'air_t', 'air_p', 'u','v','w','ts','sos', &
                'flag_1', 'flag_2') !< this two are for back-compatibility
                !> Sonic and irga variables without property "use_it"
                if (.not. LocCol(i)%useit .and. NumUserVar < MaxUserVar - 1) &
                    NumUserVar = NumUserVar + 1
            case default
                !> Variables with a custom label and without property "use_it"
                if (.not. LocCol(i)%useit .and. NumUserVar < MaxUserVar - 1) &
                    NumUserVar = NumUserVar + 1
        end select

        !> Detect whether an 4th gas calibration data column is available
        if (index(LocCol(i)%var, 'cal-ref') /= 0) Gas4CalRefCol = i
    end do

    !> If user selects a different temperature reading
    !> (instead of sonic temperature) for sensible heat flux, redefine Ts
    !> as that column, but changes instrument category to "fast_t_sensor"
    !> to remember that it does not need water vapor correction and
    !> (sonic-specific) spectral corrections.
    if (EddyProProj%Col(ts) > 0) then
        LocCol(EddyProProj%Col(ts))%useit = .true.
        LocCol(EddyProProj%Col(ts))%var = 'ts'
        LocCol(EddyProProj%Col(ts))%instr%category = 'fast_t_sensor'
        !> Search Ts or SoS from master sonic and change property in
        !> "don't use it", so now it will fall into the "non sensitive"
        !> variables group. Note that the total number of User Variables
        !> did not change
        ts_found = .false.
        do i = 1, NumCol
            if (LocCol(i)%instr%master_sonic &
                .and. (LocCol(i)%var == 'ts' .or. LocCol(i)%var == 'sos' ) &
                .and. LocCol(i)%useit) then
                LocCol(i)%useit = .false.
                ts_found = .true.
                !exit
            end if
        end do
        !> If Ts was not there instead, the number of user variables must
        !> be reduced by one, because one was used
        !> as a fast temperature
        if (.not. ts_found) NumUserVar = NumUserVar - 1
    end if
end subroutine DefineUsedVariables
