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
subroutine FileListByExt(DirIn, Ext, MatchPrototype, Prototype, doy_format, &
    GetTimestamp, Recurse, FileList, nrow, logout, indent)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    character(*), intent(in) :: indent
    character(*), intent(in) :: Ext
    character(*), intent(in) :: DirIn
    character(*), intent(in) :: Prototype
    logical, intent(in) :: MatchPrototype
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
    logical, external :: NameMatchesPrototype

    call log_msg(' inf=creating list of files to be processed')
    LogString = ' data_path=' // DirIn(1:len_trim(DirIn))
    call DoubleCharInString(LogString, slash)
    call log_msg(LogString)

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
        call system(comm_del // '"' // trim(adjustl(TmpDir)) // 'flist*.tmp"')
        call ErrorHandle(0, 0, 1)
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
            if (MatchPrototype) then
                if (NameMatchesPrototype(TmpFileName, Prototype)) then
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
            if (MatchPrototype) then
                if (NameMatchesPrototype(TmpFileName, Prototype)) then
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

    !> Define file names and timestamp if requested
    do i = 1, cnt
        FileList(i)%name = FileList(i)%path &
            (index(FileList(i)%path, slash, .true.) + 1: len_trim(FileList(i)%path))
    end do

    !> Some logging
    write(LogInteger, '(i6)') cnt
    call SchrinkString(LogInteger)
    LogString = ' num_file=' // LogInteger(1:len_trim(LogInteger))
    call log_msg(LogString)
    if (logout) write(*,'(a)') indent // '  ' // LogInteger(1:len_trim(LogInteger)) // ' files found.'

    !> Retrieve timestamps from file names if requested
    if(GetTimestamp) then
        if (logout) write(*, '(a)', advance = 'no') '  Retrieving timestamps from &
            &file names..'
        do i = 1, cnt
            call FilenameToTimestamp(FileList(i)%name, Prototype, doy_format, FileList(i)%Timestamp)
        end do
        if (logout) write(*, '(a)') ' Done.'
    end if

    !> Control on number of files found
    if (cnt == 0) then
        call log_msg(' stop=no files to process. terminating execution.')
        call log_shutdown
        call ErrorHandle(0, 0, 8)
    end if

    !> Delete temporary file
    call system(comm_del // '"' // trim(adjustl(TmpDir)) // 'flist*.tmp"' // comm_err_redirect)

    if (logout) write(*,'(a)')   indent // ' Done.'
end subroutine FileListByExt

!***************************************************************************
!
! \brief       Test if file name matches prototype
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function NameMatchesPrototype(FileName, Pattern)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*), intent(in) :: Filename
    character(*), intent(in) :: Pattern
    !> Local variables
    integer :: start


    !> Initialization
    NameMatchesPrototype = .false.

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
    NameMatchesPrototype = .true.

end function NameMatchesPrototype
