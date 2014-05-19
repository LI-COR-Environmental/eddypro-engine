!***************************************************************************
! import_native_data.f90
! ----------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Open data files depending on file type
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ImportNativeData(Filepath, FirstRecord, LastRecord, LocCol, fRaw, nrow, ncol, skip_file, N, FileEndReached)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow, ncol
    character(*), intent(in) :: Filepath
    integer, intent(out) :: N
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: skip_file
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    logical, intent(out) :: FileEndReached
    !> local variables
    integer :: i
    integer :: io_status
    integer :: read_status
    integer(kind = 1) :: rec_len
    character(256) :: datastring


    skip_file = .false.
    !> Open native data file
    select case (EddyProProj%ftype(1:len_trim(EddyProProj%ftype)))

        case ('alteddy_bin')
            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', recl = 12)
            write(LogLogical, '(L1)') io_status
            if (io_status /= 0) then
                call log_msg(' Error while opening native binary file. file skipped.')
                call ErrorHandle(1, 0, 4)
                skip_file = .true.
                return
            end if

        case ('eddymeas_bin')
            !> Open raw file in binary mode
            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', recl = 8 + (NumCol - 4) * 2)
            write(LogLogical, '(L1)') io_status
            if (io_status /= 0) then
                call log_msg(' Error while opening native binary file. file skipped.')
                call ErrorHandle(1, 0, 4)
                skip_file = .true.
                return
            end if

        case ('edisol_bin')
            !> Reads record length from first byte of the file

            open(udf, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', recl = 1)
            write(LogLogical, '(L1)') io_status
            if (io_status /= 0) then
                call log_msg(' Error while opening native binary file. file skipped.')
                call ErrorHandle(1, 0, 4)
                skip_file = .true.
                return
            end if
            read(udf, rec=1, iostat = read_status) rec_len

            !> Check that rec_len is consistent with number of variables in files, if not skip
            if (read_status /= 0 .or. rec_len < 0 .or. rec_len > 24 .or. rec_len /= size(fRaw, 2) * 2) then
                call log_msg(' Error while reading native binary file. file skipped.')
                call ErrorHandle(1, 0, 5)
                skip_file = .true.
                return
            end if
            close(udf)

            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', recl = rec_len)

        case ('generic_bin')
            open(unat, file = trim(adjustl(Filepath)), access='direct', &
                form = 'unformatted', recl = 1, iostat = io_status)
            write(LogLogical, '(L1)') io_status
            if (io_status /= 0) then
                call log_msg(' Error while opening native binary file. file skipped.')
                call ErrorHandle(1, 0, 4)
                skip_file = .true.
                return
            end if

        case ('tob1')
            !> If number of header rows is /= 0, open file in TEXT mode to read data format (IEEE4 or FP2)
            if (FileInterpreter%tob1_format == 'none' .and. FileInterpreter%header_rows > 0) then
                open(udf, file = trim(adjustl(Filepath)), status = 'old', iostat = io_status)
                write(LogLogical, '(L1)') io_status
                if (io_status /= 0) then
                    call log_msg(' Error while opening native binary file. file skipped.')
                    call ErrorHandle(1, 0, 4)
                    skip_file = .true.
                    return
                end if
                do i = 1, FileInterpreter%header_rows
                    !> Read file as text and looks for IEEE4 or FP2 in any row of the file
                    read(udf, '(a)') datastring
                    if (index(datastring, '"IEEE4"') /= 0 .or. index(datastring, '"ieee4"') /= 0) then
                        FileInterpreter%tob1_format = 'IEEE4'
                        exit
                    elseif (index(datastring, '"FP2"') /= 0 .or. index(datastring, '"fp2"') /= 0) then
                        FileInterpreter%tob1_format = 'FP2'
                        exit
                    end if
                end do
                !> In the same datastring, counts the number of ULONG variables
                FileInterpreter%ulongs = 0
                do
                    if(index(datastring, '"ULONG"') /= 0) then
                        FileInterpreter%ulongs = FileInterpreter%ulongs + 1
                        datastring = datastring(index(datastring, '"ULONG"') + 1: len_trim(datastring))
                    else
                        exit
                    end if
                end do
                close(udf)
            end if
            !> Open raw file in binary mode
            open(udf, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', recl = 8192)
            write(LogLogical, '(L1)') io_status
            if (io_status /= 0) then
                call log_msg(' Error while opening native binary file. file skipped.')
                call ErrorHandle(1, 0, 4)
                skip_file = .true.
                return
            end if

        case default
            open(unat, file = trim(adjustl(Filepath)), status = 'old', iostat = io_status)
            write(LogLogical, '(L1)') io_status
            if (io_status /= 0) then
                call log_msg(' Error while opening generic ASCII file. file skipped.')
                call ErrorHandle(1, 0, 6)
                skip_file = .true.
                return
            end if
    end select

    !> Read native file
    call ReadNativeFile(Filepath, FirstRecord, LastRecord, rec_len, &
        LocCol, fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)

    !> Define fRaw, containing both EddyPro and user-defined variables,
    !> but no 'ignore' or 'not_numeric' columns.
    !> Values are now converted into standard physical units, if needed.
    call DefineAllVarSet(LocCol, fRaw, size(fRaw, 1), size(fRaw, 2), N)
    close(unat)

end subroutine ImportNativeData

!***************************************************************************
!
! \brief       Interface to actual file reading routines
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadNativeFile(Filepath, FirstRecord, LastRecord, rec_len, &
    LocCol, fRaw, nrow, ncol, N, FileEndReached)

    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow, ncol
    integer(kind = 1), intent(in) :: rec_len
    character(*), intent(in) :: Filepath
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    integer, intent(out) :: N
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: FileEndReached


    !> select routine to read file, depending on the native file format
    select case (trim(adjustl(EddyProProj%ftype)))

        !> generic ascii files, with separated values
        case('generic_ascii', 'licor_ghg', 'toa')
            if (FileInterpreter%file_with_text) then
                call ImportAsciiWithText(FirstRecord, LastRecord, LocCol, &
                    fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)
            else
                call ImportAscii(FirstRecord, LastRecord, LocCol, &
                    fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)
            end if

        !> TOB1 files (Campbell(R) Scientific proprietary format)
        case('tob1')
            call ImportTOB1(Filepath, FirstRecord, LastRecord, LocCol, &
                fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)

        !> SLT files created with by EddySoft
        case('eddymeas_bin')
            call ImportSLTEddySoft(FirstRecord, LastRecord, LocCol, &
                fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)

        !> SLT files created with by EdiSol
        case('edisol_bin')
            call ImportSLTEdiSol(rec_len, FirstRecord, LastRecord, LocCol, &
                fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)

        case('generic_bin')
            call ImportBinary(FirstRecord, LastRecord, LocCol, &
                fRaw, size(fRaw, 1), size(fRaw, 2), N, FileEndReached)
        end select
end subroutine ReadNativeFile
