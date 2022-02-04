!***************************************************************************
! import_native_data.f90
! ----------------------
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
! \brief       Open data files depending on file type
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ImportNativeData(Filepath, FirstRecord, LastRecord, LocCol, &
    fRaw, nrow, ncol, skip_file, N, FileEndReached)
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
    !MC integer(kind = 1) :: rec_len
    integer(i1) :: rec_len
    character(ShortInstringLen) :: dataline


    skip_file = .false.
    !> Open native data file
    select case (trim(adjustl(EddyProProj%ftype)))

        case ('eddymeas_bin')
            !> Open raw file in binary mode
            !MC record length should be ncol*2 bytes,
            !   which should be the same as 8 + (NumCol - 4) * 2
            !   because NumCol is ncol without the 4 anemometer variables
            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', &
                recl = 8 + (NumCol - 4) * 2)
            if (io_status /= 0) then
                call ExceptionHandler(54)
                skip_file = .true.
                return
            end if

        case ('edisol_bin')
            !> Reads record length from first byte of the file

            open(udf, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', &
                form = 'unformatted', recl = 1)
            if (io_status /= 0) then
                call ExceptionHandler(54)
                skip_file = .true.
                return
            end if
            read(udf, rec=1, iostat = read_status) rec_len

            !> Check that rec_len is consistent with number of
            !> variables in files, if not skip
            if (read_status /= 0 .or. rec_len < 0 .or. rec_len > 24 &
                .or. rec_len /= size(fRaw, 2) * 2) then
                call ExceptionHandler(54)
                skip_file = .true.
                return
            end if
            close(udf)

            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', form = 'unformatted', &
                recl = rec_len)

        case ('generic_bin')
            open(unat, file = trim(adjustl(Filepath)), access='direct', &
                form = 'unformatted', recl = 1, iostat = io_status)
            if (io_status /= 0) then
                call ExceptionHandler(55)
                skip_file = .true.
                return
            end if

        case ('tob1')
            !> If number of header rows is /= 0, open file in TEXT mode to
            !> read data format (IEEE4 or FP2)
            if (FileInterpreter%tob1_format == 'none' &
                .and. FileInterpreter%header_rows > 0) then
                open(udf, file = trim(adjustl(Filepath)), &
                    status = 'old', iostat = io_status)
                if (io_status /= 0) then
                    call ExceptionHandler(56)
                    skip_file = .true.
                    return
                end if
                do i = 1, FileInterpreter%header_rows
                    !> Read file as text and looks for IEEE4 or FP2 in
                    !> any row of the file
                    read(udf, '(a)') dataline
                    if (index(dataline, '"IEEE4"') /= 0 &
                        .or. index(dataline, '"ieee4"') /= 0) then
                        FileInterpreter%tob1_format = 'IEEE4'
                        exit
                    elseif (index(dataline, '"FP2"') /= 0 &
                        .or. index(dataline, '"fp2"') /= 0) then
                        FileInterpreter%tob1_format = 'FP2'
                        exit
                    end if
                end do
                !> In the same dataline, counts the number of ULONG variables
                FileInterpreter%ulongs = 0
                do
                    if(index(dataline, '"ULONG"') /= 0) then
                        FileInterpreter%ulongs = FileInterpreter%ulongs + 1
                        dataline = &
                        dataline(index(dataline, '"ULONG"') + 1: len_trim(dataline))
                    else
                        exit
                    end if
                end do
                close(udf)
            end if
            !> Open raw file in binary mode
            open(udf, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', &
                form = 'unformatted', recl = 8192)
            if (io_status /= 0) then
                call ExceptionHandler(56)
                skip_file = .true.
                return
            end if

        case ('alteddy_bin')
            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status, access='direct', &
                form = 'unformatted', recl = 12)
            if (io_status /= 0) then
                call ExceptionHandler(54)
                skip_file = .true.
                return
            end if

        case default
            open(unat, file = trim(adjustl(Filepath)), status = 'old', &
                iostat = io_status)
            if (io_status /= 0) then
                call ExceptionHandler(57)
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
    !MC integer(kind = 1), intent(in) :: rec_len
    integer(i1), intent(in) :: rec_len
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
