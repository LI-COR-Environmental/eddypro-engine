!***************************************************************************
! dir_sub.f90
! -----------
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
! \brief       Collection of subs for handling directories and files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************

!***************************************************************************
!
! \brief       Create a new directory
! \author      Antonio Forgione
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function CreateDir(directory)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: directory
    !> local variables
    character(1024) :: command


    !> create dir and if one already exists, skip obvious system message
    !> redirecting through windows NUL (equivalent to linux /dev/null)
    command = 'mkdir ' // directory(1: len_trim(directory)) // comm_err_redirect
    CreateDir = system(command)
end function CreateDir

!***************************************************************************
! \file        src/dir_subs.f90
! \brief       Test if file name matches Template
! \version     5.2.0
! \date        2013-03-12
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function NameMatchesTemplate(FileName, Pattern)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: Filename
    character(*), intent(in) :: Pattern
    !> Local variables
    integer :: start


    !> Initialization
    NameMatchesTemplate = .false.

    !> Check on the length of the file name
    if(EddyProProj%ftype /= 'licor_ghg' .and. len_trim(FileName) /= len_trim(Pattern)) return

    !> Year must be done by numbers
    if (index(Pattern, 'yyyy') /= 0) then
        start = index(Pattern, 'yyyy')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
        if(FileName(start + 2 :start + 2) > '9' .or. FileName(start + 2 :start + 2) < '0') return
        if(FileName(start + 3 :start + 3) > '9' .or. FileName(start + 3 :start + 3) < '0') return
    else if (index(Pattern, 'yy') /= 0) then
        start = index(Pattern, 'yy')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
    end if
    !> Month must be done by numbers
    if (index(Pattern, 'mm') /= 0) then
        start = index(Pattern, 'mm')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
    end if
    !> Day or DOY must be done by numbers
    if (index(Pattern, 'ddd') /= 0) then
        start = index(Pattern, 'ddd')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
        if(FileName(start + 2 :start + 2) > '9' .or. FileName(start + 2 :start + 2) < '0') return
    else if (index(Pattern, 'dd') /= 0) then
        start = index(Pattern, 'dd')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
    end if
    !> Hour must be done by numbers
    if (index(Pattern, 'HH') /= 0) then
        start = index(Pattern, 'HH')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
    end if
    !> Minute must be done by numbers
    if (index(Pattern, 'MM') /= 0) then
        start = index(Pattern, 'MM')
        if(FileName(start:start) > '9' .or. FileName(start:start) < '0') return
        if(FileName(start + 1 :start + 1) > '9' .or. FileName(start + 1 :start + 1) < '0') return
    end if
    NameMatchesTemplate = .true.

end function NameMatchesTemplate

!***************************************************************************
!
! \brief       Counts the number of files in a dir \n
!               (and its subfolders if requested)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine NumberOfFilesInDir(DirIn, ext, MatchTemplate, Template, N, rN)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: DirIn
    character(*), intent(in) :: Ext
    character(*), intent(in) :: Template
    logical, intent(in) :: MatchTemplate
    integer, intent(out) :: N
    integer, intent(out) :: rN
    !> local variables
    integer :: dir_status
    integer :: open_status
    character(1024) :: comm
    character(256) :: string
    character(64) :: TmpFileName
    logical, external :: NameMatchesTemplate

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
        call system(comm_del // '"' // trim(adjustl(TmpDir)) // '"flist*.tmp')
        call ExceptionHandler(1)
        N = 0
        rN = 0
        return
    end if

    !> Open temporary file and counts how many files are present.
    N = 0
    rN = 0
    do
        read(udf, '(a)', iostat = open_status) string
        if (open_status /= 0) exit
            TmpFileName =  string(index(string, slash, .true.) + 1: len_trim(string))
            if (MatchTemplate) then
                if (NameMatchesTemplate(TmpFileName, Template)) then
                    N = N + 1
                    if (string(1:index(string, slash, .true.)) == trim(adjustl(DirIn))) rN = rN + 1
                end if
            else
                N = N + 1
                if (string(1:index(string, slash, .true.)) == trim(adjustl(DirIn))) rN = rN + 1
            end if
    end do
    close(udf)
    call system(comm_del // '"' // trim(adjustl(TmpDir)) // '"*.tmp' // comm_err_redirect)
end subroutine NumberOfFilesInDir

!***************************************************************************
!
! \brief       Counts the number of files in a selected subperiod
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function NumberOfFilesInSubperiod(FileList, nrow, StartTimestamp, EndTimestamp)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    type(FileListType), intent(in) :: FileList(nrow)
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(in) :: EndTimestamp
    !> Local variables
    integer :: i
    integer :: cnt

    !> Write in output filelist only files within the selected interval
    cnt = 0
    do i = 1, nrow
        if (FileList(i)%timestamp >= StartTimestamp .and. FileList(i)%timestamp <= EndTimestamp) &
            cnt = cnt + 1
    end do
    NumberOfFilesInSubperiod = cnt

end function NumberOfFilesInSubperiod

!***************************************************************************
!
! \brief       Adjust path spelling, replacing slashes with backslashes \n
!              or vice-versa, based on OS-specific settings.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjFilePath(FilePath, slash)
    implicit none
    !> in/out variables
    character(*), intent(in) :: slash
    character(*), intent(inout) :: FilePath
    !> local variables
    integer :: i = 0


    do i = 1, len_trim(FilePath)
        if(FilePath(i:i) == '/' .or. FilePath(i:i) == '\') FilePath(i:i) = slash
    end do
end subroutine AdjFilePath

!***************************************************************************
!
! \brief       Adjust dir spelling, replacing slashes with backslashes \n
!              or vice-versa, based on OS-specific settings, and adding
!              ending (back)slash if not present
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjDir(LocDir, slash)
    implicit none
    !> in/out variables
    character(*), intent(in) :: slash
    character(*), intent(inout) :: LocDir
    integer, external :: dirlen
    !> local variables
    integer :: i = 0


    do i = 1, dirlen(LocDir)
        if(LocDir(i:i) == '/' .or. LocDir(i:i) == '\' ) LocDir(i:i) = slash
    end do
    if(LocDir(dirlen(LocDir):dirlen(LocDir)) /= slash) &
        LocDir(dirlen(LocDir) + 1:dirlen(LocDir) + 1) = slash
end subroutine AdjDir

!***************************************************************************
!
! \brief       Forces backslashes to slashes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ForceForwardSlash(string)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: string
    !> local variables
    integer :: i = 0


    do i = 1, len_trim(string)
        if(string(i:i) == '\') string(i:i) = '/'
    end do
end subroutine ForceForwardSlash

!***************************************************************************
!
! \brief       Clean file name by eliminating possible OS-dependent extra characters
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StripFilename(filename)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: filename


    if (filename(1:2) == './') filename = filename(3:len_trim(filename))
end subroutine StripFilename

!***************************************************************************
!
! \brief       Quick look at CSV file to provide (max) number of rows and cols
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine scanCsvFile(fpath, separator, cols_from_header, nrow, ncol, failed)
    implicit none
    !> in/out variables
    integer, intent(in) :: cols_from_header
    character(*), intent(in) :: fpath
    character(*), intent(in) :: separator
    integer, intent(out) :: ncol
    integer, intent(out) :: nrow
    logical, intent(out) :: failed
    !> Local variables
    integer :: i
    integer :: tncol
    integer :: io_status
    logical :: count_cols
    character(2048) :: row
    integer, external :: countsubstring
    integer, external :: SplitCount


    failed = .false.

    !> Open file
    open(10, file=fpath, iostat=io_status)
    if (io_status /=0) then
        failed = .true.
        return
    end if

   !> Read number of columns, from provided header row if any
    nrow = 0
    ncol = 0
    i = 0
    count_cols = .true.
    do
        read(10, '(a)', iostat = io_status) row
        i = i + 1
        if (io_status > 0) cycle
        if (io_status < 0) exit
        nrow = nrow + 1
        if (count_cols) then
            tncol = SplitCount(trim(row), separator, '', .false.)
            ncol = max(tncol, ncol)
            if (i == cols_from_header) then
                ncol = tncol
                count_cols = .false.
            end if
        end if
    end do
    close(10)

    if (ncol * nrow == 0) failed = .true.
end subroutine scanCsvFile



