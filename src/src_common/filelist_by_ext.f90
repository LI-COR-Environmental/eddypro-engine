!***************************************************************************
! filelist_by_ext.f90
! -------------------
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
! \brief       Detects and stores names of all files with extension Ext in DirIn
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FileListByExt(DirIn, Ext, MatchTemplate, Template, doy_format, GetTimestamp, Recurse, FileList, nrow, logout, indent)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    character(*), intent(in) :: indent
    character(*), intent(in) :: Ext
    character(*), intent(in) :: DirIn
    character(*), intent(in) :: Template
    logical, intent(in) :: MatchTemplate
    logical, intent(in) :: doy_format
    logical, intent(in) :: GetTimestamp
    logical, intent(in) :: Recurse
    logical, intent(in) :: logout
    type(FileListType), intent(out) :: FileList(nrow)
    !> local variables
    integer :: i
    integer :: cnt
    integer :: open_status
    integer :: read_status
    integer :: dir_status
    character(1024) :: comm
    character(256) :: String
    character(64) :: TmpFileName
    logical, external :: NameMatchesTemplate


    if (logout) then
        if (Recurse) then
            write(*,'(a)') indent // ' Retrieving file names from directory: "' &
                // DirIn(1:len_trim(DirIn)) // '" and its sub-directories..'
        else
            write(*,'(a)') indent // ' Retrieving file names from directory: "' &
                // DirIn(1:len_trim(DirIn)) // '"..'
        end if
    end if

    !> List files, recursively in all cases
    select case (OS(1:len_trim(OS)))
        case('win')
            comm = 'dir "' // DirIn(1:len_trim(DirIn)) // '*'// Ext(1:len_trim(Ext)) &
                // '" ' // ' /O:N /S /B > ' // '"' // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
        case('linux')
            comm = 'find "' // DirIn(1:len_trim(DirIn)) // '" -iname *' &
                // Ext(1:len_trim(Ext)) // ' > ' // '"' // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
        case('mac')
            comm = 'find "' // DirIn(1:len_trim(DirIn)) // '" -iname *' &
                // Ext(1:len_trim(Ext)) // ' > ' // '"' // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
    end select
    dir_status = system(comm)

    call system(comm_copy // '"' // trim(adjustl(TmpDir)) // 'flist.tmp" ' // '"' // trim(adjustl(TmpDir)) // 'flist2.tmp" ' &
        // comm_out_redirect // comm_err_redirect)

    open(udf, file = trim(adjustl(TmpDir)) // 'flist2.tmp', iostat = open_status)
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
            read(udf, '(a)', iostat = read_status) String
            if (read_status /= 0) exit
            String = trim(adjustl(String))
            TmpFileName =  String(index(String, slash, .true.) + 1: len_trim(String))
            if (MatchTemplate) then
                if (NameMatchesTemplate(TmpFileName, Template)) then
                    cnt = cnt + 1
                    FileList(cnt)%path = trim(adjustl(String))
                end if
            else
                cnt = cnt + 1
                FileList(cnt)%path = trim(adjustl(String))
            end if
        end do
    else
        do
            read(udf, '(a)', iostat = read_status) String
            if (read_status /= 0) exit
            if(String(1:index(String, slash, .true.)) /= DirIn(1:len_trim(DirIn))) cycle
            String = trim(adjustl(String))
            TmpFileName =  String(index(String, slash, .true.) + 1: len_trim(String))
            if (MatchTemplate) then
                if (NameMatchesTemplate(TmpFileName, Template)) then
                    cnt = cnt + 1
                    FileList(cnt)%path = trim(adjustl(String))
                end if
            else
                cnt = cnt + 1
                FileList(cnt)%path = trim(adjustl(String))
            end if
        end do
    end if
    close(udf)

    !> Control on number of files found
    if (cnt == 0) call ExceptionHandler(8)

    !> Define file names and timestamp if requested
    do i = 1, cnt
        FileList(i)%name = FileList(i)%path &
            (index(FileList(i)%path, slash, .true.) + 1: len_trim(FileList(i)%path))
    end do

    !> Some logging
    if (logout) then
        write(LogInteger, '(i8)') cnt
        write(*,'(a)') indent // '  ' // trim(adjustl(LogInteger)) // ' files found.'
    end if

    !> Retrieve timestamps from file names if requested
    if(GetTimestamp) then
        if (logout) write(*, '(a)', advance = 'no') '  Retrieving timestamps from &
            &file names..'
        do i = 1, cnt
            call FilenameToTimestamp(FileList(i)%name, Template, doy_format, FileList(i)%Timestamp)
        end do
        if (logout) write(*, '(a)') ' Done.'
    end if

    !> Delete temporary file
    call system(comm_del // '"' // trim(adjustl(TmpDir)) // 'flist*.tmp"' // comm_err_redirect)

    if (logout) write(*,'(a)')   indent // ' Done.'
end subroutine FileListByExt
