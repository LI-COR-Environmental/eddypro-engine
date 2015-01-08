!***************************************************************************
! make_dataset.f90
! ----------------
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
    character(256) :: PathOut
    character(MaxOutstringLen) :: DataString
    character(MaxOutstringLen) :: DataString_utf8
    character(MaxOutstringLen) :: TmpDataString
    character(MaxOutstringLen) :: ErrString
    character(10) :: fdate
    character(5)  :: ftime
    type (DateType) :: fTimestamp

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

    !> Copy header from input file and copy it to output file
    do i = 1, hnrow
        read(udf, '(a)', iostat = read_status) DataString
        call latin1_to_utf8(DataString, DataString_utf8)
        write(udf2, '(a)') trim(adjustl(DataString_utf8))
    end do

    !> Now creates actual dataset, inserting either actual
    !> results or ErrString depending on whether results are available
    !> for each time step of the MasterTimeSeries
    periods_loop: do i = StartIndx, EndIndx
        !> Search for current time step in the file
        do

            if (MasterTimeSeries(i) > fTimestamp + DateStep) then
                call AddErrorString(udf2, MasterTimeSeries(i), &
                    ErrString, len(ErrString), &
                    PathIn == GHGEUROPE_Path(1:len_trim(GHGEUROPE_Path)), &
                    AddNoFile)
                cycle periods_loop
            end if

            !> Read the data string
            read(udf, '(a)', iostat = read_status) DataString

            if(read_status /= 0) then
                !> Insert here a control
                exit
            end if

            !> Retrieve date and time from Datastring
            if (AddNoFile) then
                TmpDataString = &
                DataString(index(DataString, separator) + 1: len_trim(DataString))
            else
                TmpDataString = DataString
            end if
            fdate = TmpDataString(1:index(TmpDataString, separator) - 1)
            TmpDataString = &
            TmpDataString(index(TmpDataString, separator) + 1: len_trim(TmpDataString))
            ftime = TmpDataString(1:index(TmpDataString, separator) - 1)

            !> Convert into timestamp and take it back to the
            !> beginning of the averaging period
            call DateTimeToDateType(fdate, ftime, fTimestamp)
            fTimestamp = fTimestamp - DateStep

            !> Check if fTimestamp is equal to current time step timestamp
            if(fTimestamp == MasterTimeSeries(i)) then
                !> If yes, time to write on output file
                write(udf2,'(a)') trim(adjustl(DataString))
                cycle periods_loop
            elseif (fTimestamp > MasterTimeSeries(i)) then
                !> If not, checks if timestamp in file is later than current
                !> in MasterTimeSeries. If so, writes error string on output,
                !> backspaces and cycle periods_loops
                if (i < EndIndx) then
                    call AddErrorString(udf2, MasterTimeSeries(i+1), &
                        ErrString, len(ErrString), &
                        PathIn == GHGEUROPE_Path(1:len_trim(GHGEUROPE_Path)), &
                        AddNoFile)
                end if
                backspace(udf)
                cycle periods_loop
            end if
        end do
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
subroutine AddErrorString(unt, Timestamp, ErrString, LenErrStr, &
    IsGhgEuropeFile, AddNoFile)
    use m_common_global_Var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    integer, intent(in) :: LenErrStr
    character(LenErrStr) :: ErrString
    logical, intent(in) :: IsGhgEuropeFile
    logical, intent(in) :: AddNoFile
    type (DateType), intent(in) :: Timestamp
    !> local variables
    character(10) :: date
    character(5) :: time
    integer :: int_doy
    real(kind = dbl) :: float_doy
    character(32) :: char_doy
    character(LenErrStr) :: String

    !> Convert Timestamp to date and time and calculate doy
    call DateTypeToDateTime(Timestamp, date, time)
    call DateTimeToDOY(date, time, int_doy, float_doy)
    call WriteDatumFloat(float_doy, char_doy, EddyProProj%err_label)
    call ShrinkString(char_doy)

    !> Create output string
    if (AddNoFile) then
        if (IsGhgEuropeFile) then
            String = 'not_enough_data,' // date(1:10) // ',' // time // ',' &
                     // EddyProProj%err_label(1:len_trim(EddyProProj%err_label)) &
                     // ErrString(index(ErrString, ',,') + 1: len_trim(ErrString))
        else
            String = 'not_enough_data,' // date(1:10) // ',' // time // ',' &
                     // char_doy(1: index(char_doy, '.')+ 3) &
                     // ErrString(index(ErrString, ',,') + 1: len_trim(ErrString))
        end if
    else
        String = date(1:10) // ',' // time // ',' &
                // char_doy(1: index(char_doy, '.')+ 3) &
                // ErrString (index(ErrString, ',,') + 1: len_trim(ErrString))
    end if

    !> write on file
    write(unt, '(a)') trim(adjustl(String))


end subroutine AddErrorString
