!***************************************************************************
! init_ex_vars.f90
! ----------------
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
    integer :: gas
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: InitializationPerformed
    type (ExType) :: lEX

    write(*,'(a)') &
        ' Initializing retrieval of EddyPro-RP results from file: '
    write(*,'(a)') '  "' // trim(adjustl(AuxFile%ex)) // '"..'

    !> Open EX file
    open(udf, file = AuxFile%ex, status = 'old', iostat = open_status)

    !> Exit with error in case of problems opening the file
    if (open_status /= 0) call ExceptionHandler(60)

    write(*, '(a)') '  File found, importing content..'

    !> Store header to string, for writing it on output
    read(udf, '(a)') fluxnet_header

    !> Retrieve label of forth gas from header
    ! read(udf, '(a)') dataline
    ! substr = dataline(index(dataline, 'ru_ch4'):index(dataline, 'ru_ch4') + 30)
    ! g4lab = substr(8: index(substr, '_flux') - 1)
    ! g4l = len_trim(g4lab)
    ! !> Retrieve names of user variables from header
    ! if (len_trim(dataline) >= index(dataline, 'num_user_var') + 13) then
    !     UserVarHeader = &
    !         dataline(index(dataline, 'num_user_var') + 13: len_trim(dataline))
    ! else
    !     UserVarHeader = ''
    ! end if


                                !*********************************************** Extract UserVarHeader
                                !*********************************************** Determine how to handle g4lab
    g4lab = 'GS4'               !*********************************************** Temporary
    g4l = len_trim(g4lab)
    UserVarHeader = 'XXXXXX'

    !> Initialize variables that are determined for the whole
    !> dataset (presence of certain variables)
    Diag7200%present = .false.
    Diag7500%present = .false.
    Diag7700%present = .false.
    fcc_var_present = .false.
    FCCMetadata%ru = .false.
    FCCMetadata%ac_freq = -1
    DateStep = DateType(0, 0, 0, 0, ierror)

    !> Cycle on all records
    NumRecords = 0
    NumValidRecords = 0
    InitializationPerformed = .false.

    do
        !> Read essentials record
        call ReadExRecord('', udf, -1, lEx, ValidRecord, EndOfFileReached)
        if (EndOfFileReached) exit

        !> Counts
        NumRecords = NumRecords + 1
        if (ValidRecord) NumValidRecords = NumValidRecords + 1

        !> Handles dates
        if (ValidRecord .and. NumValidRecords == 1) &
            call DateTimeToDateType(lEx%end_date, lEX%end_time, StartTimestamp)
        if (ValidRecord) &
            call DateTimeToDateType(lEx%end_date, lEX%end_time, EndTimestamp)

        !> Initializations
        if (ValidRecord .and. .not. InitializationPerformed) then

            !> Look for variable presence (u thru GS4)
            if (lEx%WS /= error) fcc_var_present(u:w) = .true.
            if (lEx%Ts /= error) fcc_var_present(ts)  = .true.
            do gas = co2, gas4
                fcc_var_present(gas) = lEx%measure_type_int(gas) /= ierror .or. fcc_var_present(gas)  
            end do
                
            !> Determine whether LI-COR's flags are available
            if (.not. Diag7200%present) then
                do j = 1, 9
                    if (lEx%licor_flags(j) /= error) then
                        Diag7200%present = .true.
                        exit
                    end if
                end do
            end if

            if (.not. Diag7500%present) then
                do j = 10, 13
                    if (lEx%licor_flags(j) /= error) then
                        Diag7500%present = .true.
                        exit
                    end if
                end do
            end if

            if (.not. Diag7700%present) then
                do j = 14, 29
                    if (lEx%licor_flags(j) /= error) then
                        Diag7700%present = .true.
                        exit
                    end if
                end do
            end if

            !> Reads DateStep
            if (DateStep == DateType(0, 0, 0, 0, ierror)) DateStep = DateType(0, 0, 0, 0, nint(lEx%avrg_length))

            !> Define whether random uncertainty was calculated by
            !> looking at only 1 value (if one value is -6999d0, all
            !> of them are the same)
            if (lEx%rand_uncer(u) == aflx_error) FCCMetadata%ru = .false.

            !> Acquisition frequency and gas analyser path type for H2O
            if (FCCMetadata%ac_freq <= 0) FCCMetadata%ac_freq = lEx%ac_freq
            FCCMetadata%H2oPathType = lEx%instr(ih2o)%path_type
        end if

        if (all(fcc_var_present) .and. Diag7200%present .and. Diag7500%present .and. Diag7700%present .and. &
           FCCMetadata%ac_freq > 0 .and. FCCMetadata%ru .and. DateStep /= DateType(0, 0, 0, 0, ierror)) then
            InitializationPerformed = .true.
        end if
    end do
    close(udf)

    !> Adjust start timestamp so that Start/End define the whole period
    !> From beginning of first period to end of last period
    StartTimestamp = StartTimestamp - DateStep
    write(*,'(a)') ' Done.'
end subroutine InitExVars
