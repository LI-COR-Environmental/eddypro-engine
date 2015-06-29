!***************************************************************************
! dir_sub.f90
! -----------
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
    character(CommLen) :: comm


    !> create dir and if one already exists, skip obvious system message
    !> redirecting through windows NUL (equivalent to linux /dev/null)
    comm = 'mkdir ' // directory(1: len_trim(directory)) // comm_err_redirect
    CreateDir = system(comm)
end function CreateDir

!***************************************************************************
! \file
! \brief       Returns whether FileName matches Template
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function NameMatchesTemplate(FileName, Template, HardMatch)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: Filename
    character(*), intent(in) :: Template
    logical, intent(in) :: HardMatch
    !> Local variables
    integer :: s
    integer :: s_year, s_month, s_day, s_hour, s_minute
    integer :: e_year, e_month, e_day, e_hour, e_minute
    integer :: s_ts, e_ts
    logical, external :: is_not_numeric
    logical, external :: strings_match
    character(1), parameter :: wild_card = '?'

    !> Initialization
    NameMatchesTemplate = .false.

    !> Check on the length of filename
    if(EddyProProj%run_env /= 'embedded' &
        .and. len_trim(FileName) /= len_trim(Template)) return

    !> Check timestamp
    !> Year must be done by numbers
    if (index(Template, 'yyyy') /= 0) then
        s = index(Template, 'yyyy')
        s_year = index(Template, 'yyyy')
        e_year = index(Template, 'yyyy') + 3
        if (is_not_numeric(Filename(s:s+3))) return
    else if (index(Template, 'yy') /= 0) then
        s = index(Template, 'yy')
        s_year = index(Template, 'yy')
        e_year = index(Template, 'yy') + 1
        if (is_not_numeric(Filename(s:s+1))) return
    end if

    !> Month must be done by numbers
    if (index(Template, 'mm') /= 0) then
        s = index(Template, 'mm')
        s_month = index(Template, 'mm')
        e_month = index(Template, 'mm') + 1
        if (is_not_numeric(Filename(s:s+1))) return
    else
        s_month = nint(error)
        e_month = - nint(error)
    end if

    !> Day or DOY must be done by numbers
    if (index(Template, 'ddd') /= 0) then
        s = index(Template, 'ddd')
        s_day = index(Template, 'ddd')
        e_day = index(Template, 'ddd') + 2
        if (is_not_numeric(Filename(s:s+2))) return
    else if (index(Template, 'dd') /= 0) then
        s = index(Template, 'dd')
        s_day = index(Template, 'dd')
        e_day = index(Template, 'dd') + 1
        if (is_not_numeric(Filename(s:s+1))) return
    end if

    !> Hour must be done by numbers
    if (index(Template, 'HH') /= 0) then
        s = index(Template, 'HH')
        s_hour = index(Template, 'HH')
        e_hour = index(Template, 'HH') + 1
        if (is_not_numeric(Filename(s:s+1))) return
    end if

    !> Minute must be done by numbers
    if (index(Template, 'MM') /= 0) then
        s = index(Template, 'MM')
        s_minute = index(Template, 'MM')
        e_minute = index(Template, 'MM') + 1
        if (is_not_numeric(Filename(s:s+1))) return
    end if

    !> If running in embedded mode, template only contains timestamp
    !> so if it got to here, test is passed
    if (EddyProProj%run_env == 'embedded') then
        NameMatchesTemplate = .true.
        return
    endif

    !> Check rest of the template if hard-match was requested
    if(HardMatch) then
        s_ts = min(s_year, s_month, s_day, s_hour, s_minute)

        !> Check prefix
        if (s_ts > 1 &
            .and. .not. strings_match(Template(1:s_ts-1), &
            Filename(1:s_ts-1), wild_card)) return

        !> Check suffix
        e_ts = max(e_year, e_month, e_day, e_hour, e_minute)
        if (e_ts < len_trim(Template) - 1 &
            .and. .not. strings_match(Template(e_ts+1:len_trim(Template)), &
            Filename(e_ts+1:len_trim(Filename)), wild_card)) return
    end if

    !> If all tests were passed, filename matches template
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
    character(CommLen) :: comm
    character(PathLen) :: dataline
    character(FilenameLen) :: TmpFileName
    logical, external :: NameMatchesTemplate

    !> List files, recursively in all cases
    select case (OS(1:len_trim(OS)))
        case('win')
            comm = 'dir "' // DirIn(1:len_trim(DirIn)) // '*' &
                // Ext(1:len_trim(Ext)) &
                // '" ' // ' /O:N /S /B > ' // '"' &
                // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
        case('linux')
            comm = 'find "' // DirIn(1:len_trim(DirIn)) &
                // '" -iname *' &
                // Ext(1:len_trim(Ext)) // ' > ' // '"' &
                // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
        case('mac')
            comm = 'find "' // DirIn(1:len_trim(DirIn)-1) &
                // '" -iname *' &
                // Ext(1:len_trim(Ext)) // ' > ' // '"' &
                // trim(adjustl(TmpDir)) // 'flist.tmp" ' // comm_err_redirect
    end select
    dir_status = system(comm)

    !> Exit with error if dir command failed
    if (dir_status /= 0) call ExceptionHandler(86)

    call system(comm_copy // '"' // trim(adjustl(TmpDir)) &
        // 'flist.tmp" ' // '"' // trim(adjustl(TmpDir)) // 'flist2.tmp" ' &
        // comm_out_redirect // comm_err_redirect)

    open(udf, file = trim(adjustl(TmpDir)) // 'flist2.tmp', iostat = open_status)

    !> Initialize
    N = 0
    rN = 0

    !> control on temporary file
    if (open_status /= 0) then
        close(udf)
        call system(comm_del // '"' // trim(adjustl(TmpDir)) // '"flist*.tmp')
        call ExceptionHandler(1)
        return
    end if

    !> Open temporary file and counts files
    do
        read(udf, '(a)', iostat = open_status) dataline
        if (open_status /= 0) exit
            TmpFileName = &
                dataline(index(dataline, slash, .true.) + 1: len_trim(dataline))
            if (MatchTemplate) then
                if (NameMatchesTemplate(TmpFileName, &
                    Template, MatchTemplate)) then
                    N = N + 1
                    if (dataline(1:index(dataline, slash, .true.)) &
                        == trim(adjustl(DirIn))) rN = rN + 1
                end if
            else
                N = N + 1
                if (dataline(1:index(dataline, slash, .true.)) &
                    == trim(adjustl(DirIn))) rN = rN + 1
            end if
    end do
    close(udf)
    call system(comm_del // '"' // trim(adjustl(TmpDir)) &
        // '"*.tmp' // comm_err_redirect)
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
integer function NumberOfFilesInSubperiod(FileList, nrow, &
    StartTimestamp, EndTimestamp)
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
        if (FileList(i)%timestamp >= StartTimestamp &
            .and. FileList(i)%timestamp <= EndTimestamp) &
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
! \brief       Forces slashes to either forward (default) or backward.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ForceSlash(string, backslash)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: string
    logical, intent(in) :: backslash
    !> local variables
    integer :: i = 0
    character(1) :: slash

    slash = '/'
    if (backslash) slash = '\'
    do i = 1, len_trim(string)
        if(string(i:i) == '\' .or. string(i:i) == '/' ) string(i:i) = slash
    end do
end subroutine ForceSlash

!***************************************************************************
!
! \brief       Clean file name by eliminating possible
!              OS-dependent extra characters
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
    use m_typedef
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
    character(LongInstringLen) :: row
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
