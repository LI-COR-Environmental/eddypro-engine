!***************************************************************************
! parse_ini_file.f90
! ------------------
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
! \brief       Parse any files in EddyPro INI-format and stores detected tags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ParseIniFile(IniFile, key, NumTags, CharTags, nnum, nchar,&
           NumTagFound, CharTagFound, IniFileNotFound)
    use m_common_global_var
    implicit none
    !in/out variables
    !character(*), intent(in) :: IniFile
    character(*) :: IniFile
    character(*), intent(in) :: key
    integer, intent(in) :: nnum
    integer, intent(in) :: nchar
    logical, intent(out) :: NumTagFound(nnum)
    logical, intent(out) :: CharTagFound(nchar)
    logical, intent(out) :: IniFileNotFound
    type(Numerical), intent(inout) :: NumTags(nnum)
    type(Text), intent(inout) :: CharTags(nchar)
    !> local variables
    integer :: nlines
    integer :: io_status
    type(text) :: Tags(MaxNLinesIni)


    IniFileNotFound = .false.
    open(udf, file = IniFile, status = 'old', iostat = io_status)
    if (io_status == 0) then
        !> parse the ini file and store all tags found in it
        call StoreIniTags(udf, key, Tags, nlines)
        close(udf)
        !> search relevant tags among those found in the file
        call SearchLocalTags(Tags, nlines, NumTags, CharTags, nnum, nchar,&
             NumTagFound, CharTagFound)
             return
    else
        call ExceptionHandler(7)
        IniFileNotFound = .true.
    end if
end subroutine ParseIniFile

!***************************************************************************
!
! \brief       Store all tags and values detected in the INI file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StoreIniTags(uini, key, Tags, nlines)
    use m_common_global_var
    implicit none
    !in/out variables
    integer, intent(in) :: uini
    character(*), intent(in) :: key
    type(Text), intent(out) :: Tags(MaxNLinesIni)
    integer, intent(out) :: nlines
    !> local variables
    integer :: io_status
    integer :: separ
    integer :: com
    character(512) :: String


    nlines = 0
    !> if key is not given
    if (len_trim(key) == 0) then
        do
            read(uini, '(a)', iostat = io_status) String
            if (io_status /= 0) exit
            call stripstr(String)
            if ((String(1:1) /= "[").and.&
                (String(1:1) /= ";").and.&
                (len_trim(String) /= 0)) then
                com = index(String, ";")
                if (com /= 0) then
                    String(com:len_trim(String)) = ""
                    call stripstr(String)
                end if
                nlines = nlines + 1
                separ = index(String, "=")
                Tags(nlines)%Label = String(1:separ - 1)
                Tags(nlines)%Value = String(separ + 1:len_trim(String))
             end if
            call clearstr(String)
        end do
        return
    end if

    !> if key is given
    outloop: do
        read(uini, '(a)', iostat = io_status) String
        if (io_status /= 0) exit outloop
        call stripstr(String)
        if (String(1:1) == "[" .and. &
            index(String, key(1:len_trim(key))) == 2) then
            ! if key was found, start to store Tags
            inloop: do
                read(uini, '(a)', iostat = io_status) String
                if (io_status /= 0) exit outloop
                call stripstr(String)
                if (len_trim(String) == 0) cycle inloop
                ! if another group is found, check if it belongs to the key
                if (String(1:1) == "[") then
                    if (index(String, key(1:len_trim(key))) /= 2) then
                        exit outloop
                    else
                        cycle inloop
                    end if
                end if
                ! if string is meaningful, store Tag
                if (String(1:1) /= ";" .and. len_trim(String) /= 0) then
                    com = index(String, ";")
                    if (com /= 0) then
                        String(com:len_trim(String)) = ""
                        call stripstr(String)
                    end if
                    nlines = nlines + 1
                    separ = index(String, "=")
                    Tags(nlines)%Label = String(1:separ - 1)
                    Tags(nlines)%Value = String(separ + 1:len_trim(String))
                end if
                call clearstr(String)
            end do inloop
        else
            cycle outloop
        end if
    end do outloop
end subroutine StoreIniTags

!***************************************************************************
!
! \brief       Search relevant tags, among those found in the INI file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SearchLocalTags(Tags, nlines, NumTags, CharTags, nnum, nchar,&
           NumTagFound, CharTagFound)
    use m_common_global_var
    implicit none
    !in/out variables
    integer, intent(in) :: nlines
    integer, intent(in) :: nnum
    integer, intent(in) :: nchar
    type(Text), intent(in) :: Tags(MaxNLinesIni)
    type(Numerical), intent(out) :: NumTags(nnum)
    type(Text), intent(out) :: CharTags(nchar)
    logical, intent(out) :: NumTagFound(nnum)
    logical, intent(out) :: CharTagFound(nchar)
    !> local variables
    integer :: read_stat
    integer :: i
    integer :: j


    !> search numeric variables
    NumTagFound = .false.
    do i = 1, nnum
        loop1: do j = 1, nlines
            if(index(Tags(j)%Label(1:len_trim(Tags(j)%Label)), &
                NumTags(i)%label(1:len_trim(NumTags(i)%label))) > 0 ) then
                NumTagFound(i) = .true.
                read(Tags(j)%value, *, iostat=read_stat) NumTags(i)%value
                if (read_stat /= 0) NumTags(i)%value = error
                exit loop1
            end if
        end do loop1
    end do

    !> search text variables
    CharTagFound = .false.
    CharTags%value = ''
    do i = 1, nchar
        loop2: do j = 1, nlines
            if(index(Tags(j)%Label(1:len_trim(Tags(j)%Label)), &
                CharTags(i)%label(1:len_trim(CharTags(i)%label))) > 0 ) then
                CharTagFound(i) = .true.
                CharTags(i)%value = Tags(j)%value
                exit loop2
           end if
       end do loop2
    end do
end subroutine SearchLocalTags
