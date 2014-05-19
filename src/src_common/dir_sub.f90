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
subroutine NumberOfFilesInDir(DirIn, ext, N, rN)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: DirIn
    character(*), intent(in) :: Ext
    integer, intent(out) :: N
    integer, intent(out) :: rN
    !> local variables
    integer :: dir_status
    integer :: open_status
    character(1024) :: comm
    character(256) :: string

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
    write(LogLogical, '(L1)') dir_status
    LogString = ' dir_error=' // LogLogical
    call log_msg(LogString)

    call system(comm_copy // '"' // trim(adjustl(TmpDir)) // 'flist.tmp" ' // '"' // trim(adjustl(TmpDir)) // 'flist2.tmp" ' &
        // comm_out_redirect // comm_err_redirect)

    open(udf, file = trim(adjustl(TmpDir)) // 'flist2.tmp', iostat = open_status)

    !> control on temporary file
    write(LogLogical, '(L1)') open_status
    LogString = ' open_flist_error=' // LogLogical
    call log_msg(LogString)
    if (open_status /= 0) then
        close(udf)
        call system(comm_del // '"' // trim(adjustl(TmpDir)) // '"flist*.tmp')
        call ErrorHandle(0, 0, 1)
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
        N = N + 1
        if (string(1:index(string, slash, .true.)) == trim(adjustl(DirIn))) rN = rN + 1
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



