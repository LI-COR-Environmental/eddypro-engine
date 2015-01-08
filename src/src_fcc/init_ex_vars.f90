!***************************************************************************
! init_ex_vars.f90
! ----------------
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
! \brief       Reads essentials file, retrieving all information that might \n
!              be useful to other programs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitExVars(StartTimestamp, EndTimestamp, NumRecords, NumValidRecords)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(out) :: NumRecords
    integer, intent(out) :: NumValidRecords
    type(DateType), intent(out) :: StartTimestamp
    type(DateType), intent(out) :: EndTimestamp
    !> local variables
    integer :: open_status
    integer :: j
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: InitializationPerformed
    type (EXType) :: lEX
    character(4000) :: string
    character(100) :: substr

    write(*,'(a)') &
        ' Initializing retrieval of EddyPro-RP results from file: '
    write(*,'(a)') '  "' // trim(adjustl(AuxFile%ex)) // '"..'

    !> Open EX file
    open(udf, file = AuxFile%ex, status = 'old', iostat = open_status)

    !> Exit with error in case of problems opening the file
    if (open_status /= 0) call ExceptionHandler(60)

    write(*, '(a)') '  File found, importing content..'
    !> Retrieve label of forth gas from header
    read(udf, '(a)') string
    substr = string(index(string, 'ru_ch4'):index(string, 'ru_ch4') + 30)
    g4lab = substr(8: index(substr, '_flux') - 1)
    g4l = len_trim(g4lab)
    !> Retrieve names of user variables from header
    if (len_trim(string) >= index(string, 'num_user_var') + 13) then
        UserVarHeader = &
            string(index(string, 'num_user_var') + 13: len_trim(string))
    else
        UserVarHeader = ''
    end if

    !> Initialize variables that are determined for the whole
    !> dataset (presence of certain variables)
    Diag7200%present = .false.
    Diag7500%present = .false.
    Diag7700%present = .false.

    !> Cycle on all records
    NumRecords = 0
    NumValidRecords = 0
    InitializationPerformed = .false.
    FCCMetadata%ru = .true.
    do
        !> Read essentials record
        call ReadExRecord('', udf, -1, lEx, ValidRecord, EndOfFileReached)
        if (EndOfFileReached) exit

        !> Counts
        NumRecords = NumRecords + 1
        if (ValidRecord) NumValidRecords = NumValidRecords + 1

        !> Handles dates
        if (ValidRecord .and. NumValidRecords == 1) &
            call DateTimeToDateType(lEX%date, lEX%time, StartTimestamp)
        if (ValidRecord) &
            call DateTimeToDateType(lEX%date, lEX%time, EndTimestamp)

        !> Some initializations
        if (ValidRecord .and. .not. InitializationPerformed) then
            !> Determine whether LI-COR's flags are available
            do j = 1, 9
                if (lEx%licor_flags(j) /= error) then
                    Diag7200%present = .true.
                    exit
                end if
            end do
            do j = 10, 13
                if (lEx%licor_flags(j) /= error) then
                    Diag7500%present = .true.
                    exit
                end if
            end do
            do j = 14, 29
                if (lEx%licor_flags(j) /= error) then
                    Diag7700%present = .true.
                    exit
                end if
            end do

            !> Reads DateStep
            DateStep = DateType(0, 0, 0, 0, nint(lEx%avrg_length))

            !> Define whether random uncertainty was calculated by
            !> looking at only 1 value (if one value is -6999d0, all
            !> of them are the same)
            if (lEx%rand_uncer(u) == aflx_error) FCCMetadata%ru = .false.

            !> Acquisition frequency and gas analyser path type for H2O
            FCCMetadata%ac_freq = lEx%ac_freq
            FCCMetadata%H2oPathType = lEx%instr(ih2o)%path_type

            InitializationPerformed = .true.
        end if
    end do
    close(udf)

    !> Adjust Start/End timestamps to define the
    !> boundaries of the MasterTimeseries
    StartTimestamp = StartTimestamp - DateStep
    EndTimestamp = EndTimestamp - DateStep
    write(*,'(a)') ' Done.'
end subroutine InitExVars
