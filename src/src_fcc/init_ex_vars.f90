!***************************************************************************
! init_ex_vars.f90
! ----------------
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
    integer :: st
    integer :: en
    integer :: j
    integer :: gas
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: InitializationPerformed
    type (ExType) :: lEX
    include '../src_common/interfaces_1.inc'

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

    st = index(fluxnet_header, ',FCH4,') + 6
    en = st + index(fluxnet_header(st:), ',') - 2
    g4lab = fluxnet_header(st+1:en)
    call lowercase(g4lab)
    g4l = len_trim(g4lab)
    

    st = index(fluxnet_header, 'NUM_CUSTOM_VARS') + 16
    en = index(fluxnet_header, 'NUM_BIOMET_VARS') - 2
    UserVarHeader = fluxnet_header(st:en)
    UserVarHeader = replace2(UserVarHeader, 'CUSTOM_', '')
    call lowercase(UserVarHeader)

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
