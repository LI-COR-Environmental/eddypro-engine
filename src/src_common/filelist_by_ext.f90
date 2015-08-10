!***************************************************************************
! filelist_by_ext.f90
! -------------------
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
! \brief       Detects and stores names of all files with extension Ext in DirIn
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FileListByExt(DirIn, Ext, MatchTemplate, HardMatch, Template, &
        doy_format, GetTimestamp, Recurse, FileList, nrow, printout, indent)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    character(*), intent(in) :: indent
    character(*), intent(in) :: Ext
    character(*), intent(in) :: DirIn
    character(*), intent(in) :: Template
    logical, intent(in) :: MatchTemplate
    logical, intent(in) :: HardMatch
    logical, intent(in) :: doy_format
    logical, intent(in) :: GetTimestamp
    logical, intent(in) :: Recurse
    logical, intent(in) :: printout
    type(FileListType), intent(out) :: FileList(nrow)
    !> local variables
    integer :: i
    integer :: cnt
    integer :: open_status
    integer :: read_status
    integer :: dir_status
    character(CommLen) :: comm
    character(PathLen) :: dataline
    character(FilenameLen) :: fname
    logical, external :: NameMatchesTemplate


    if (printout) then
        write(*,'(a)') indent // ' Retrieving file names from directory:'
        write(*,'(a)') indent // '  "' // trim(adjustl(DirIn)) // '"'
        if (Recurse) &
            write(*,'(a)') indent // '   and its sub-directories..'
    end if

    !> List files, recursively in all cases
    select case (OS(1:len_trim(OS)))
        case('win')
            comm = 'dir "' // DirIn(1:len_trim(DirIn)) // '*' &
                // Ext(1:len_trim(Ext)) // '" ' // ' /O:N /S /B > ' // '"' &
                // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
        case('linux')
            comm = 'find "' // DirIn(1:len_trim(DirIn)) // '" -iname *' &
                // Ext(1:len_trim(Ext)) // ' > ' // '"' // trim(adjustl(TmpDir)) &
                // 'flist.tmp" ' // comm_err_redirect
        case('mac')
            comm = 'find "' // DirIn(1:len_trim(DirIn)-1) // '" -iname *' &
                // Ext(1:len_trim(Ext)) // ' > ' // '"' // trim(adjustl(TmpDir)) &
                // 'flist.tmp" ' // comm_err_redirect
    end select
    dir_status = system(comm)

    call system(comm_copy // '"' // trim(adjustl(TmpDir)) // 'flist.tmp" ' &
        // '"' // trim(adjustl(TmpDir)) // 'flist2.tmp" ' &
        // comm_out_redirect // comm_err_redirect)

    open(udf, file = trim(adjustl(TmpDir)) // 'flist2.tmp', &
        iostat = open_status)
    !> control on temporary file
    if (open_status /= 0) then
        close(udf)
        call system(comm_del // '"' // trim(adjustl(TmpDir)) // 'flist*.tmp"')
        call ExceptionHandler(1)
        cnt = 0
        return
    end if

    !> Read file paths
    cnt = 0
    if (Recurse) then
        do
            read(udf, '(a)', iostat = read_status) dataline
            if (read_status /= 0) exit
            dataline = trim(adjustl(dataline))
            fname =  &
                dataline(index(dataline, slash, .true.) + 1: len_trim(dataline))
            if (MatchTemplate) then
                if (NameMatchesTemplate(fname, Template, HardMatch)) then
                    cnt = cnt + 1
                    FileList(cnt)%path = trim(adjustl(dataline))
                end if
            else
                cnt = cnt + 1
                FileList(cnt)%path = trim(adjustl(dataline))
            end if
        end do
    else
        do
            read(udf, '(a)', iostat = read_status) dataline
            if (read_status /= 0) exit
            if(dataline(1:index(dataline, slash, .true.)) &
                /= DirIn(1:len_trim(DirIn))) cycle
            dataline = trim(adjustl(dataline))
            fname =  &
                dataline(index(dataline, slash, .true.) + 1: len_trim(dataline))
            if (MatchTemplate) then
                if (NameMatchesTemplate(fname, Template, HardMatch)) then
                    cnt = cnt + 1
                    FileList(cnt)%path = trim(adjustl(dataline))
                end if
            else
                cnt = cnt + 1
                FileList(cnt)%path = trim(adjustl(dataline))
            end if
        end do
    end if
    close(udf)

    !> Control on number of files found
    if (cnt == 0) then
        if (MatchTemplate) then
            call ExceptionHandler(78)
        else
            call ExceptionHandler(8)
        end if
    end if

    !> Define file names and timestamp if requested
    do i = 1, cnt
        call basename(FileList(i)%path, FileList(i)%name, slash)
    end do

    !> Some logging
    if (printout) then
        write(LogInteger, '(i8)') cnt
        write(*,'(a)') indent // '  ' // trim(adjustl(LogInteger)) &
            // ' files found.'
    end if

    !> Retrieve timestamps from file names if requested
    if(GetTimestamp) then
        if (printout) write(*, '(a)', advance = 'no') &
            indent // '  Retrieving timestamps from &
            &file names..'
        do i = 1, cnt
            call FilenameToTimestamp(FileList(i)%name, Template, doy_format, &
                FileList(i)%Timestamp)
        end do
        if (printout) write(*, '(a)') ' Done.'
    end if

    !> Delete temporary file
    call system(comm_del // '"' // trim(adjustl(TmpDir)) // 'flist*.tmp"' &
        // comm_err_redirect)

    if (printout) write(*,'(a)')   indent // ' Done.'
end subroutine FileListByExt
