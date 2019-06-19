!***************************************************************************
! make_dataset.f90
! ----------------
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
! \brief       Create dataset filling timestamp gaps with error codes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine MakeDataset(PathIn, MasterTimeSeries, nrow, StartIndx, &
    EndIndx, AddNoFile, hnrow)
    use m_common_global_Var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: StartIndx
    integer, intent(in) :: EndIndx
    integer, intent(in) :: hnrow
    type (DateType), intent(in) :: MasterTimeSeries(nrow)
    character(*), intent(in) :: PathIn
    logical, intent(in) :: AddNoFile
    !> local variables
    integer :: i
    integer :: tmp_indx
    integer :: open_status
    integer :: read_status
    character(PathLen) :: PathOut
    character(LongOutstringLen) :: dataline
    character(LongOutstringLen) :: dataline_utf8
    character(LongOutstringLen) :: Tmpdataline
    character(LongOutstringLen) :: ErrString
    character(10) :: fdate
    character(5)  :: ftime
    type (DateType) :: fTimestamp
    logical :: add_error


    tmp_indx = index(PathIn, TmpExt)
    PathOut  = PathIn(1:tmp_indx - 1)

    !> Initialize ErrString by reading one valid line from the file
    ErrString = ''

    call InitContinuousDataset(PathIn, ErrString, hnrow)

    !> Open input file
    open(udf, file = PathIn(1:len_trim(PathIn)), status = 'old', &
        iostat = open_status, encoding = 'utf-8')
    if (open_status /= 0 ) then
        write(*,'(a)')
        write(*,'(a)') '  A problem occurred while opening file: ', &
            trim(adjustl(PathIn))
        write(*,'(a)') '   File not imported.'
        return
    end if

    !> Open output file
    open(udf2, file = PathOut, iostat = open_status, encoding = 'utf-8')
    if (open_status /= 0) then
        call ExceptionHandler(67)
        return
    end if

    !> Copy header from input file and paste it to output file
    do i = 1, hnrow
        read(udf, '(a)', iostat = read_status) dataline
        call latin1_to_utf8(dataline, dataline_utf8)
        write(udf2, '(a)') trim(adjustl(dataline_utf8))
    end do

    !> Now creates actual dataset, inserting either actual
    !> results or ErrString depending on whether results are available
    !> for each time step of the MasterTimeSeries
    add_error = .false.
    periods_loop: do i = StartIndx, EndIndx
        !> Search for current time step in the file
        do
            !> Read the data string
            read(udf, '(a)', iostat = read_status) dataline

            if(read_status /= 0) then
                add_error = .true.
                !> Insert here a control
                exit
            end if

            !> Retrieve date and time from dataline
            if (AddNoFile) then
                Tmpdataline = &
                dataline(index(dataline, separator) + 1: len_trim(dataline))
            else
                Tmpdataline = dataline
            end if
            fdate = Tmpdataline(1:index(Tmpdataline, separator) - 1)
            Tmpdataline = &
            Tmpdataline(index(Tmpdataline, separator) + 1: len_trim(Tmpdataline))
            ftime = Tmpdataline(1:index(Tmpdataline, separator) - 1)

            !> Convert into timestamp and take it back to the
            !> beginning of the averaging period
            call DateTimeToDateType(fdate, ftime, fTimestamp)
            fTimestamp = fTimestamp - DateStep
            if (MasterTimeSeries(i) > fTimestamp + DateStep) then
                call AddErrorString(udf2, MasterTimeSeries(i), &
                    ErrString, len(ErrString), AddNoFile)
                cycle periods_loop
            end if

            !> Check if fTimestamp is equal to current time step timestamp
            if(fTimestamp == MasterTimeSeries(i)) then
                !> If yes, time to write on output file
                write(udf2,'(a)') trim(adjustl(dataline))
                cycle periods_loop
            elseif (fTimestamp > MasterTimeSeries(i)) then
                !> If not, checks if timestamp in file is later than current
                !> in MasterTimeSeries. If so, writes error string on output,
                !> backspaces and cycle periods_loops
                if (i < EndIndx) then
                    call AddErrorString(udf2, MasterTimeSeries(i+1), &
                        ErrString, len(ErrString), AddNoFile)
                end if
                backspace(udf)
                cycle periods_loop
            end if
        end do

        !> This takes care of time periods following the end of the raw data period
        if (add_error .and. i < EndIndx) then
            call AddErrorString(udf2, MasterTimeSeries(i+1), &
                ErrString, len(ErrString), AddNoFile)
            add_error = .false.
        end if
    end do periods_loop

    close(udf)
    close(udf2)
end subroutine MakeDataset


!***************************************************************************
! \file        src/make_dataset.f90
! \brief       Create dataset filling timestamp gaps with error code
! \version     4.5.0
! \date        2012-10-22
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddErrorString(unt, Timestamp, ErrString, LenErrStr, AddNoFile)
    use m_common_global_Var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    integer, intent(in) :: LenErrStr
    character(LenErrStr) :: ErrString
    logical, intent(in) :: AddNoFile
    type (DateType), intent(in) :: Timestamp
    !> local variables
    character(10) :: date
    character(5) :: time
    integer :: int_doy
    real(kind = dbl) :: float_doy
    character(32) :: char_doy
    character(LenErrStr) :: dataline

    !> Convert Timestamp to date and time and calculate doy
    call DateTypeToDateTime(Timestamp, date, time)
    call DateTimeToDOY(date, time, int_doy, float_doy)
    write(char_doy, '(G20.7)') float_doy
    call ShrinkString(char_doy)

    !> Create output string
    if (AddNoFile) then
        dataline = 'not_enough_data,' // date(1:10) // ',' // time // ',' &
                    // char_doy(1: index(char_doy, '.')+ 3) &
                    // ErrString(index(ErrString, ',,') + 1: len_trim(ErrString))
    else
        dataline = date(1:10) // ',' // time // ',' &
                // char_doy(1: index(char_doy, '.')+ 3) &
                // ErrString (index(ErrString, ',,') + 1: len_trim(ErrString))
    end if

    !> write on file
    write(unt, '(a)') trim(adjustl(dataline))

end subroutine AddErrorString
