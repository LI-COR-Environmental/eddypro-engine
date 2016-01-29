!***************************************************************************
! import_tob1.f90
! ---------------
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
! \brief       Reads in TOB1 binary file \n
!              Selects only hot columns, discarding "ignore" \n
!              and "not_numeric" columns
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ImportTOB1(Filepath, FirstRecord, LastRecord, LocCol, fRaw, nrow, ncol, N, FileEndReached)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: FirstRecord
    integer, intent(in) :: LastRecord
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(out) :: N
    character(*), intent(in) :: Filepath
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(out) :: fRaw(nrow, ncol)
    logical, intent(out) :: FileEndReached
    !> local variables
    integer :: i
    integer :: j
    integer :: jj
    integer :: read_status
    integer :: rec_num
    integer :: nhlines
    integer :: first_data_byte
    integer(kind = 2) :: int2_fp2
    integer :: int4_fp2
    character(8192) :: chunk_head
    character(8195) :: chunk_head2
    character(8192) :: chunk
    character(8195) :: chunk2
    character(8192) :: newchunk
    character(8195) :: newchunk2
    character(PathLen)  :: OnlyDataPath
    character(128)  :: TmpFilepath
    real(kind = sgl) :: TmpfRaw(nrow, NumCol)
    real(kind = sgl) :: Dataline(NumCol)
    logical :: repeat
    logical :: TmpFileExists
    logical :: DeleteFile
    type(ColType) :: TmpCol(MaxNumCol)


    FileEndReached = .false.

    if (FileInterpreter%header_rows > 0) then
        DeleteFile = .True.
        !> Create path of tmp file, that either already exists or must be created
        !> by copying only data from original TOB1 file
        TmpFilepath = trim(adjustl(TmpDir)) // &
                      Filepath(index(Filepath, slash, .true.) + 1: len_trim(Filepath))

        !> Check if Tmp file exists
        inquire(file = TmpFilepath, exist = TmpFileExists)
        if (TmpFileExists) then
            close(udf)
            OnlyDataPath = TmpFilepath
        else
            open(udf2, file = TmpFilepath)
            !> first chunk to eliminate header
            read(udf, rec = 1, iostat = read_status) chunk_head
            nhlines = 0
            first_data_byte = 1
            do i = 1, 8191
                if(chunk_head(i:i+1) == char(13) // char(10)) then
                    nhlines = nhlines + 1
                    if (nhlines == FileInterpreter%header_rows) then
                        first_data_byte = i + 2
                        exit
                    end if
                end if
            end do
            write(udf2, '(a)', advance = 'no') chunk_head(first_data_byte:len_trim(chunk_head))

            !> Header is out, copies the rest in scratch file
            rec_num = 2
            read(udf, rec = rec_num, iostat = read_status) chunk
            do
                rec_num = rec_num + 1
                newchunk = ''
                read(udf, rec = rec_num, iostat = read_status) newchunk
                if(read_status /= 0) then
                    do i = len(chunk), 1, -1
                        if(chunk(i:i) /= char(0)) exit
                        if(chunk(i:i) == char(0)) chunk(i:i) = ''
                    end do
                    write(udf2,'(a)') chunk(1:len_trim(chunk))
                    exit
                else
                    write(udf2,'(a)', advance = 'no') chunk(1:len(chunk))
                    chunk = newchunk
                end if
            end do
            close(udf)
            close(udf2)
            OnlyDataPath = TmpFilepath
        end if
    else
        close(udf)
        OnlyDataPath = Filepath
        DeleteFile = .False.
    end if

    !> Read scratch file without header
    select case (FileInterpreter%tob1_format)
        case('IEEE4')
            open(udf, file = OnlyDataPath(1:len_trim(OnlyDataPath)), status = 'old', &
                access='direct', form = 'unformatted', recl = 4)
            i = 0
            N = 0
            rec_num = 0
            record_loop: do
                i = i + 1
                !> Normal exit instruction
                if (N > LastRecord - FirstRecord) exit record_loop
                !> Read one line of data
                do j = 1, NumCol
                    read(udf, rec = rec_num + j, iostat = read_status) Dataline(j)
                    if(read_status /= 0) then
                        FileEndReached = .true.
                        if (DeleteFile) then
                            close(udf,  status = 'delete')
                        else
                            close(udf)
                        end if
                        exit record_loop
                    end if
                end do
                rec_num = rec_num + NumCol
                !> Cycle until FirstRecord
                if (i < FirstRecord) cycle record_loop
                !> If data line is good, copy into TmpfRaw
                N = N + 1
                TmpfRaw(N, :) = Dataline(:)

                !> If only weird values were read, exit loop
                if(maxval(abs(TmpfRaw(N, 1:NumCol))) < 1e-15) exit record_loop
            end do record_loop
            close(udf)

        case('FP2')
            open(udf, file = OnlyDataPath(1:len_trim(OnlyDataPath)), status = 'old', &
                access='direct', form = 'unformatted', recl = 2)
            i = 0
            N = 0
            rec_num = 0
            record_loop2: do
                i = i + 1
               !> Normal exit instruction
                if (N > LastRecord - FirstRecord) exit record_loop2
                 !> Read one line of data (skipping ulong values if present)
                do j = 1, FileInterpreter%ulongs
                    read(udf, rec = rec_num + j, iostat = read_status) int2_fp2
                end do
                rec_num = rec_num + FileInterpreter%ulongs
                do j = 1, NumCol
                    read(udf, rec = rec_num + j, iostat = read_status) int2_fp2
                    if(read_status /= 0) then
                        FileEndReached = .true.
                        if (DeleteFile) then
                            close(udf,  status = 'delete')
                        else
                            close(udf)
                        end if

                        exit record_loop2
                    end if
                    if (int2_fp2 >= 0) then
                        int4_fp2 = int2_fp2
                    else
                        int4_fp2 = int2_fp2 + 65536_4
                    end if
                    Dataline(j) = FP2(int4_fp2)
                end do
                rec_num = rec_num + NumCol
                !> Cycle until FirstRecord
                if (i < FirstRecord) cycle record_loop2
                N = N + 1
                TmpfRaw(N, :) = Dataline(:)
                !> If only weird values were read, exit loop
                if(maxval(abs(TmpfRaw(N, 1:NumCol))) < 1e-15) exit record_loop2
            end do record_loop2
            close(udf)

        case default
            call ExceptionHandler(31)
    end select

    !> Check if most data were imported correctly
    repeat = .false.
!    if (N /= 0) then
!        if (.not. allocated(mask)) allocate(mask(N))
!        do j = 1, NumCol
!            if (LocCol(j)%var /= 'ignore' .and. LocCol(j)%var /= 'not_numeric') then
!                mask(1:N) = abs(TmpfRaw(1:N, j)) < 1d6
!                if (count(mask) < N * 7d-1) then
!                    repeat = .true.
!                    if (allocated(mask)) deallocate(mask)
!                    exit
!                end if
!            end if
!        end do
!        if (allocated(mask)) deallocate(mask)
!    end if


    !> If a column was found with weird values, reads the file again with another chunk size
    if(repeat) then
        open(udf, file = Filepath(1:len_trim(Filepath)), status = 'old', &
            access='direct', form = 'unformatted', recl = 8193)

        if (FileInterpreter%header_rows > 0) then

            !> Create path of tmp file, that either already exists or must be created
            !> by copying only data from original TOB1 file
            TmpFilepath = trim(adjustl(TmpDir)) // &
                          Filepath(index(Filepath, slash, .true.) + 1: len_trim(Filepath))

            !> Check if Tmp file is open
            inquire(file = TmpFilepath, exist = TmpFileExists)
            if (TmpFileExists) then
                close(udf)
                OnlyDataPath = TmpFilepath
            else
                open(udf2, file = TmpFilepath)
                !> first chunk to eliminate header
                read(udf, rec = 1, iostat = read_status) chunk_head2
                nhlines = 0
                first_data_byte = 1
                do i = 1, 8194
                    if(chunk_head2(i:i+1) == char(13) // char(10)) then
                        nhlines = nhlines + 1
                        if (nhlines == FileInterpreter%header_rows) then
                            first_data_byte = i + 2
                            exit
                        end if
                    end if
                end do
                write(udf2, '(a)', advance = 'no') chunk_head2(first_data_byte:len_trim(chunk_head2))

                !> Header is out, copies the rest in scratch file
                rec_num = 2
                read(udf, rec = rec_num, iostat = read_status) chunk2
                do
                    rec_num = rec_num + 1
                    newchunk2 = ''
                    read(udf, rec = rec_num, iostat = read_status) newchunk2
                    if(read_status /= 0) then
                        do i = len(chunk2), 1, -1
                            if(chunk2(i:i) /= char(0)) exit
                            if(chunk2(i:i) == char(0)) chunk2(i:i) = ''
                        end do
                        write(udf2,'(a)') chunk2(1:len_trim(chunk2))
                        exit
                    else
                        write(udf2,'(a)', advance = 'no') chunk2(1:len(chunk2))
                        chunk2 = newchunk2
                    end if
                end do
                close(udf)
                close(udf2)
                OnlyDataPath = TmpFilepath
            end if
        else
            close(udf)
            OnlyDataPath = Filepath
        end if

        !> Read scratch file without header
        select case (FileInterpreter%tob1_format)

            case('IEEE4')
                open(udf, file = OnlyDataPath(1:len_trim(OnlyDataPath)), status = 'old', &
                    access='direct', form = 'unformatted', recl = 4)
                i = 0
                N = 0
                rec_num = 0
                record_loop3: do
                    i = i + 1
                    !> Normal exit instruction
                    if (N > LastRecord - FirstRecord) exit record_loop3
                    !> Read one line of data
                    do j = 1, NumCol
                        read(udf, rec = rec_num + j, iostat = read_status) Dataline(j)
                        if(read_status /= 0) then
                            FileEndReached = .true.
                            if (DeleteFile) then
                                close(udf,  status = 'delete')
                            else
                                close(udf)
                            end if
                            exit record_loop3
                        end if
                    end do
                    rec_num = rec_num + NumCol
                    !> Cycle until FirstRecord
                    if (i < FirstRecord) cycle record_loop3
                    !> If data line is good, copy into TmpfRaw
                    N = N + 1
                    TmpfRaw(N, :) = Dataline(:)
                    !> If only weird values were read, exit loop
                    if(maxval(abs(TmpfRaw(N, 1:NumCol))) < 1e-15) exit record_loop3
                end do record_loop3
                close(udf)

            case('FP2')
                open(udf, file = OnlyDataPath(1:len_trim(OnlyDataPath)), status = 'old', &
                    access='direct', form = 'unformatted', recl = 2)
                i = 0
                N = 0
                rec_num = 0
                record_loop4: do
                    i = i + 1
                   !> Normal exit instruction
                    if (N > LastRecord - FirstRecord) exit record_loop4
                     !> Read one line of data (skipping ulong values if present)
                    do j = 1, FileInterpreter%ulongs
                        read(udf, rec = rec_num + j, iostat = read_status) int2_fp2
                    end do
                    rec_num = rec_num + FileInterpreter%ulongs
                    do j = 1, NumCol
                        read(udf, rec = rec_num + j, iostat = read_status) int2_fp2
                        if(read_status /= 0) then
                            FileEndReached = .true.
                            if (DeleteFile) then
                                close(udf,  status = 'delete')
                            else
                                close(udf)
                            end if
                            exit record_loop4
                        end if
                        if (int2_fp2 >= 0) then
                            int4_fp2 = int2_fp2
                        else
                            int4_fp2 = int2_fp2 + 65536_4
                        end if
                        Dataline(j) = FP2(int4_fp2)
                    end do
                    rec_num = rec_num + NumCol
                    !> Cycle until FirstRecord
                    if (i < FirstRecord) cycle record_loop4
                    N = N + 1
                    TmpfRaw(N, :) = Dataline(:)
                    !> If only weird values were read, exit loop
                    if(maxval(abs(TmpfRaw(N, 1:NumCol))) < 1e-15) exit record_loop4
                end do record_loop4
                close(udf)

            case default
                call ExceptionHandler(31)
        end select
    end if

    !> Store only hot columns
    jj = 0
    do j = 1, NumCol
        if (LocCol(j)%var /= 'ignore' .and. LocCol(j)%var /= 'not_numeric') then
            jj = jj + 1
            TmpCol(jj) = LocCol(j)
            fRaw(1:N, jj) = TmpfRaw(1:N, j)
            if (Gas4CalRefCol == j) Gas4CalRefCol = jj
        end if
    end do
    LocCol = TmpCol
end subroutine ImportTOB1
