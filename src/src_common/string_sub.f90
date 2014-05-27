!***************************************************************************
! string_sub.f90
! --------------
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
! \brief       Collection of subroutines for handling strings of characters
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
! \brief       Clear a string, substituting characters with blank spaces
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine clearstr(string)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: string
    !> local variables
    integer :: i

    do i = 1, len(string)
         string(i:i) = ' '
    end do
end subroutine clearstr

!***************************************************************************
!
! \brief       Extract file name from a path, if a slash is identified
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine basename(path, filename, slash)
    implicit none
    !> in/out variables
    character(*), intent(in) :: path
    character(*), intent(in) :: slash
    character(*), intent(out) :: filename


    if (index(path, slash) == 0) then
        filename = path
    else
        filename = path(index(path, slash, .true.) + 1: len_trim(path))
    end if
end subroutine basename

!***************************************************************************
!
! \brief       Determine length of a directory path
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated  Soon to be deprecated
! \test
! \todo
!***************************************************************************
integer function dirlen(string)
    implicit none
    !> in/out variables
    character(*), intent(in) :: string
    !> local variables
    integer :: i

    i = 0
    do while ((string(i + 1:i + 2) /= '  ').and.(string(i + 1:i + 1) /= char(0)))
        i = i + 1
    end do
    dirlen = i
end

!***************************************************************************
!
! \brief       Compare two strings for identity
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function eqstring(string1, string2)
    implicit none
    !> in/out variables
    character(*), intent(in) :: string1
    character(*), intent(in) :: string2


    eqstring = ((index(string1(:len_trim(string1)), string2(:len_trim(string2))) > 0) .and. &
                (len_trim(string1) == len_trim(string2)))
end function eqstring

!***************************************************************************
!
! \brief       Remove trailing and leading spaces and \n
!              quotation marks from string
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine stripstr(string)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: string
    !> local variables
    character :: quote(2)
    integer :: i
    integer :: j
    integer :: k
    integer :: changed
    integer :: oldlen

    !> Remove trailing spaces
    do i = len(string), 1, -1
        if(string(i:i) /= ' ') exit
        if(string(i:i) == ' ') string(i:i) = ''
    end do

    !> Remove leading spaces
    changed = 1
    i = 1
    k = 0
    oldlen = len(string)
    do while((changed == 1).and.(k <= oldlen))
        changed = 0
        if(string(i:i) == ' ' .or. string(i:i) == char(0)) then
            do j = i + 1, len(string) - 1
                string(j - 1:j - 1) = string(j:j)
            end do
            changed = 1
            k = k + 1
        else
            i = i + 1
        end if
    end do

    !> Remove outer quotes
    quote(1) = "'"
    quote(2) = '"'
    do k = 1, 2
        if(string(1:1) == quote(k)) then
            do i = len(string), 1, - 1
                if(string(i:i) == quote(k)) then
                    string(i:i) = char(0)
                    do j = 1, len(string)
                        string(j:j) = string(j + 1:j + 1)
                    end do
                 end if
              end do
        end if
    end do
end subroutine stripstr

!***************************************************************************
!
! \brief       Eliminate blank characters from a string
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ShrinkString(string)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: string
    !> local variables
    integer :: i

    i = 0
    do
        i = i + 1
        if (string(i:i) == ' ') then
            string = string(1:i-1) // string(i+1:len_trim(string))
            i = i - 1
        end if
        if (i >= len_trim(string)) exit
    end do
end subroutine ShrinkString

!***************************************************************************
!
! \brief       Append string2 to string1
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine strcat(string1, string2)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: string1
    character(*), intent(in) :: string2
    !> local variables
    integer :: i
    integer :: len1

    len1 = len_trim(string1)
    do i = 1, len_trim(string2)
        string1(len1 + i:len1 + i) = string2(i:i)
    end do
end subroutine strcat

!***************************************************************************
!
! \brief       Convert a character of a specified length \n
!              into an integer number
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine char2int(string, num, len)
    implicit none
    !> in/out variables
    integer, intent(in) :: len
    character(len), intent(in) :: string
    integer, intent(out) :: num
    !> local variables
    integer :: i
    integer :: j

    num = 0
    do j = 1, len
        do i = 48, 57
            if(string(j:j) == char(i)) num = num + (i - 48)*10**(len - j)
        end do
    end do
end subroutine char2int

!***************************************************************************
!
! \brief       Convert an integer into a character of specified length
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine int2char(num, string, len)
    implicit none
    !> in/out variables
    integer, intent(in) :: len
    integer, intent(in) :: num
    character(len), intent(inout) :: string
    !> local variables
    integer :: i
    integer :: j
    integer :: aux
    logical :: check

    check = .false.
    aux = num
    do j = 1, len
        do i = 1, 9
            if (((aux / (i*10**(len - j))) >= 1).and.((aux / ((i + 1)*10**(len - j))) < 1)) then
            check = .true.
            string(j:j) = char(i + 48)
            aux = aux - i*10**(len - j)
            end if
        end do
        if(.not. check) string(j:j) = char(48)
        check = .false.
    end do
end subroutine int2char

!***************************************************************************
!
! \brief       Add a datum to a "separator"-separated string, and \n
!              adds the separator
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddDatum(dataline, datum, separator)
    implicit none
    !> in/out variables
    character(*), intent(inout) :: datum
    character(*), intent(in) :: separator
    character(*), intent(inout) :: dataline

    if(len_trim(datum) == 0) then
        dataline = dataline(1:len_trim(dataline)) // ','
    else
        call stripstr(datum)
        call strcat(dataline, datum(1:len_trim(datum)))
        call strcat(dataline, separator(1:len_trim(separator)))
    end if
end subroutine AddDatum

!***************************************************************************
!
! \brief       Write an integer datum as a string, \n
!              or replace it with the user-defined error string if the case
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteDatumInt(int_datum, char_datum, err_label)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: int_datum
    character(*), intent(out) :: char_datum
    character(*), intent(in) :: err_label

    if (int_datum /= nint(error)) then
        write(char_datum, *) int_datum
        call ShrinkString(char_datum)
    else
        char_datum = err_label(1:len_trim(err_label))
    end if
end subroutine WriteDatumInt

!***************************************************************************
!
! \brief       Write a real datum as a string, \n
!              or replace it with the user-defined error string if the case
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteDatumFloat(float_datum, char_datum, err_label)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: float_datum
    character(*), intent(out) :: char_datum
    character(*), intent(in) :: err_label

    if (float_datum /= error) then
        write(char_datum, *) float_datum
    else
        char_datum = err_label(1:len_trim(err_label))
    end if
end subroutine WriteDatumFloat

!***************************************************************************
!
! \brief       Tranforms string into upper case
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine uppercase(str)
    implicit none
    character(*), intent(inout) :: str
    integer :: i, del

    del = iachar('a') - iachar('A')

    do i = 1, len_trim(str)
        if (lge(str(i:i), 'a') .and. lle(str(i:i), 'z')) then
            str(i:i) = achar(iachar(str(i:i)) - del)
        end if
    end do
end subroutine uppercase

!***************************************************************************
!
! \brief       Tranforms string into lower case
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine lowercase(str)
    implicit none
    character(*), intent(inout) :: str
    integer :: i, del

    del = iachar('a') - iachar('A')

    do i = 1, len_trim(str)
        if (lge(str(i:i), 'A') .and. lle(str(i:i), 'Z')) then
            str(i:i) = achar(iachar(str(i:i)) + del)
        end if
    end do
end subroutine lowercase

!***************************************************************************
!
! \brief       Double 'char' anytime it is encountered in 'string'
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DoubleCharInString(string, char)
    implicit none
    ! in/out variables
    character(*), intent(in) :: char
    character(*), intent(inout) :: string
    ! local variables
    integer :: i

    i = 1
    do while (i <= len_trim(string))
        if (string(i:i) == char) then
            string = string(1:i) // char // string(i+1:len_trim(string))
            i = i + 1
        end if
        i = i + 1
    end do
end subroutine DoubleCharInString

!***************************************************************************
!
! \brief       Eliminate repeated character in string, most notably useful
!              to eliminate multiple contiguous delimiters
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine StripConsecutiveChar(string, char)
    implicit none
    !> in/out variables
    character(*), intent(in) :: char
    character(*), intent(inout) :: string
    !> local variables
    integer :: i

    i = 2
    do while (i <= len_trim(string) - 1)
        if ((string(i:i) == char) .and. string(i:i) == string(i-1:i-1)) then
            string = string(1:i-1) // string(i+1:len_trim(string))
            i = i - 1
        end if
        i = i + 1
    end do
    !> Special case of last character
    if (string(len_trim(string):len_trim(string)) == &
        string(len_trim(string)-1:len_trim(string)-1)) then
        string = string(1:len_trim(string)-1)
    end if
end subroutine StripConsecutiveChar

!***************************************************************************
! \file        src/string_sub.f90
! \brief       Returns True is software version "ver_new" is
!              strictly more recent than "ver_old"
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function NewerSwVer(sw_ver_new, sw_ver_old)
    implicit none
    !> in/out variables
    character(*), intent(in) :: sw_ver_new
    character(*), intent(in) :: sw_ver_old
    !> local variables
    integer :: new
    integer :: old
    integer :: dot
    character(8) :: str
    character(16) :: ver_new
    character(16) :: ver_old


    ver_new = sw_ver_new
    ver_old = sw_ver_old
    !> Compare Major version
    !> New
    dot = index(ver_new, '.')
    str = ver_new(1: dot - 1)
    call char2int(str, new, len_trim(str))
    !> Old
    dot = index(ver_old, '.')
    str = ver_old(1: dot - 1)
    call char2int(str, old, len_trim(str))
    !> comparison
    if (new > old) then
        NewerSwVer = .true.
        return
    elseif (new < old) then
        NewerSwVer = .False.
        return
    end if

    !> Compare Minor version
    !> New
    ver_new = ver_new(dot + 1: len_trim(ver_new))
    dot = index(ver_new, '.')
    str = ver_new(1: dot - 1)
    call char2int(str, new, len_trim(str))
    !> Old
    ver_old = ver_old(dot + 1: len_trim(ver_old))
    dot = index(ver_old, '.')
    str = ver_old(1: dot - 1)
    call char2int(str, old, len_trim(str))
    !> comparison
    if (new > old) then
        NewerSwVer = .true.
        return
    elseif (new < old) then
        NewerSwVer = .False.
        return
    end if

    !> Compare build version
    !> New
    ver_new = ver_new(dot + 1: len_trim(ver_new))
    call char2int(ver_new, new, len_trim(str))
    !> Old
    ver_old = ver_old(dot + 1: len_trim(ver_old))
    call char2int(ver_old, old, len_trim(str))
    !> comparison
    if (new > old) then
        NewerSwVer = .true.
        return
    else
        NewerSwVer = .false.
        return
    end if
end function NewerSwVer

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine latin1_to_utf8(latin,utf8)
    character(*), intent(in)  :: latin
    character(*), intent(out) :: utf8
    integer :: unicode,k,n,i
    character(4) :: cutf8
    i=0
    do k=1,LEN_TRIM(latin)
      call latin1_to_unicode(latin(k:k),unicode)
      call unicode_to_utf8(unicode,cutf8,n)
      if(i+n > len(utf8)) exit
      utf8(i+1:i+n)=cutf8(1:n)
      i=i+n
    end do
    utf8(i+1:)=' '
end subroutine

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine latin1_to_unicode(latin,unicode)
    character, intent(in)  :: latin
    integer, intent(out) :: unicode
    unicode=IACHAR(latin)
end subroutine

!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine unicode_to_utf8(unicode,utf8,n)
    integer, parameter :: i4=selected_int_kind(9)
!    integer, parameter :: h80=Z'80',hC0=Z'C0',hE0=Z'E0',hF0=Z'F0',hF4=Z'F4',h800=Z'800',h8000=Z'8000'
    integer, parameter :: h80   = 128
    integer, parameter :: hC0   = 192
    integer, parameter :: hE0   = 224
    integer, parameter :: hF0   = 240
    !integer, parameter :: hF4   = 244
    integer, parameter :: h800  = 2048
    integer, parameter :: h8000 = 32768

    integer(i4), intent(in)  :: unicode
    character, intent(out)   :: utf8(*)
    integer(i4), intent(out) :: n
    integer :: ia,ib,ic,id
    if(unicode < h80) then
      utf8(1)=ACHAR(unicode)
      n=1
    else if(unicode < h800) then
      ia=RSHIFT(unicode,6)
      ib=unicode-LSHIFT(ia,6)
      utf8(1)=ACHAR(hC0+ia)
      utf8(2)=ACHAR(h80+ib)
      n=2
    else if(unicode < h8000) then
      ia=RSHIFT(unicode,12)
      ib=RSHIFT(unicode-LSHIFT(ia,12),6)
      ic=unicode-LSHIFT(ia,12)-LSHIFT(ib,6)
      utf8(1)=ACHAR(hE0+ia)
      utf8(2)=ACHAR(h80+ib)
      utf8(3)=ACHAR(h80+ic)
      n=3
    else
      ia=RSHIFT(unicode,18)
      ib=RSHIFT(unicode-LSHIFT(ia,18),12)
      ic=RSHIFT(unicode-LSHIFT(ia,18)-LSHIFT(ib,12),6)
      id=unicode-LSHIFT(ia,18)-LSHIFT(ib,12)-LSHIFT(ic,6)
      utf8(1)=ACHAR(hF0+ia)
      utf8(2)=ACHAR(h80+ib)
      utf8(3)=ACHAR(h80+ic)
      utf8(4)=ACHAR(h80+id)
      n=4
    end if
end subroutine
