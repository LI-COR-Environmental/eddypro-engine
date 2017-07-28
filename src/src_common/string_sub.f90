!***************************************************************************
! string_sub.f90
! --------------
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
! \brief       Convert an integer into a string of max 32 chars
!              If a 'pad' value is specified, pads with leading zeros
!              to reach a string length equal to 'pad'.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine int2char(num, string, pad)
    implicit none
    !> in/out variables
    integer, intent(in) :: pad
    integer, intent(in) :: num
    character(*) :: string
    !> local variables
    integer :: i
    character(32) :: str


    !> Convert integer to string and adjust left
    call clearstr(str)
    write(str, '(i32)') num
    str = adjustl(str)

    !> Pad with zeros on the left if
    !> passed length is /= zero and > length of str
    if (pad > len_trim(str)) then
        i = pad - len_trim(str)
        str(i + 1:) = str
        str(1:i) = repeat('0', i)
    end if

    string = trim(adjustl(str))
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
!
! \brief       Retrieve software version from string
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
function SwVerFromString(string)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: string
    type(SwVerType) :: SwVerFromString
    !> Local variables
    character(32) :: swVerString
    character(6) :: chunk
    logical, external :: is_not_numeric
    integer, external :: CountCharInString


    swVerString = trim(adjustl(string))

    !> If swVerString is empty, return error software version
    if (len_trim(swVerString) == 0) then
        SwVerFromString = errSwVer
        return
    end if

    !> If swVerString has no '.', return error software version
    if (CountCharInString(swVerString, '.') == 0) then
        SwVerFromString = errSwVer
        return
    end if

    !> If swVerString does not contain revision, add a dummy one
    if (CountCharInString(swVerString, '.') == 1) &
        swVerString = trim(swVerString) // '.01'

    !> If swVerString has 4 fields, crop to 3 fields
    if (CountCharInString(swVerString, '.') == 3) &
        swVerString = swVerString(1:index(swVerString, '.', .true.) - 1)

    !> Major version
    chunk = swVerString(1:index(swVerString, '.')-1)
    if (is_not_numeric(chunk)) then
        SwVerFromString = errSwVer
        return
    end if
    call char2int(trim(chunk), SwVerFromString%major, len_trim(chunk))

    !> Minor version
    chunk = swVerString(index(swVerString, '.')+1: index(swVerString, '.', .true.)-1)
    if (is_not_numeric(chunk)) then
        SwVerFromString = errSwVer
        return
    end if
    call char2int(trim(chunk), SwVerFromString%minor, len_trim(chunk))

    !> Revision
    chunk = swVerString(index(swVerString, '.', .true.)+1:len(swVerString))
    if (is_not_numeric(chunk)) then
        SwVerFromString = errSwVer
        return
    end if
    call char2int(trim(chunk), SwVerFromString%revision, len_trim(chunk))
end function

!***************************************************************************
!
! \brief       Compare software versions and returns .true. if v1 is
!              strictly newer than v2
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function CompareSwVer(v1, v2)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(SwVerType), intent(in) :: v1
    type(SwVerType), intent(in) :: v2


    !>Compare Major
    if (v1%major > v2%major) then
        CompareSwVer = .true.
        return
    elseif (v1%major < v2%major) then
        CompareSwVer = .false.
        return
    end if
    !>If Major is equal, compare Minor
    if (v1%minor > v2%minor) then
        CompareSwVer = .true.
        return
    elseif (v1%minor < v2%minor) then
        CompareSwVer = .false.
        return
    end if
    !>If Major and Minor are equal, compare Revision
    if (v1%revision > v2%revision) then
        CompareSwVer = .true.
    else
        CompareSwVer = .false.
    end if
end function CompareSwVer

!***************************************************************************
!
! \brief       Compare software versions and returns .true. if v1 is
!              strictly newer than v2
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
logical function EqualSwVer(v1, v2)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(SwVerType), intent(in) :: v1
    type(SwVerType), intent(in) :: v2


    !>Compare Major
    if (v1%major == v2%major .and. &
        v1%minor == v2%minor .and. &
        v1%revision == v2%revision) then
        EqualSwVer = .true.
    else
        EqualSwVer = .false.
    end if
end function EqualSwVer

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

!***************************************************************************
!
! \brief       Counts how many times string s2 occurs in string s1
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function countsubstring(s1, s2) result(c)
    character(*), intent(in) :: s1, s2
    integer :: p, posn

    c = 0
    if(len(s2) == 0) return
    p = 1
    do
        posn = index(s1(p:), s2)
        if(posn == 0) return
        c = c + 1
        p = p + posn + len(s2)
    end do
end function


!***************************************************************************
!
! \brief       Reads row and returns number of items excluding those
!              matching substring (case insensitive)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
integer function SplitCount(string, delimiter, exclude, caseSensitive) result(cnt)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: string
    character(*), intent(in) :: delimiter
    character(*), intent(in) :: exclude
    logical, intent(in) :: caseSensitive
    !> Local variables
    integer :: del
    character(len(string)) :: item
    character(len(string)) :: lstring
    character(len(exclude)) :: lexclude


    lstring = trim(string)
    lexclude = trim(exclude)

    !> If exclusion is not case sensitive, convert everything to upper-case
    if (.not. caseSensitive) then
        call uppercase(lstring)
        call uppercase(lexclude)
    end if

    cnt = 0
    do
        del = index(lstring, delimiter)
        if (del == 1) then
            lstring = lstring(2: len_trim(lstring))
            cycle
        end if
        if (del > 1 .and. (len_trim(lstring) >= del)) then
            item = lstring(1:del)
            cnt = cnt + 1
            if (len(lexclude) > 0 .and. &
                index(item, lexclude) /= 0) cnt = cnt - 1
            lstring = lstring(del+1: len_trim(lstring))
        else
            if (len_trim(lstring) > 0) then
                item = trim(lstring)
                cnt = cnt + 1
                if (len(lexclude) > 0 .and. &
                    index(item, lexclude) /= 0) cnt = cnt - 1
            end if
            exit
        end if
    end do
end function SplitCount

!***************************************************************************
!
! \brief       In 'string', replace "replace_what" with "replace_with"
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
! function replace(string, what, with, outlen) result(nstring)
!     use m_common_global_var
!     implicit none
!     !> in/out variables
!     integer, intent(in) :: outlen
!     character(*), intent(in) :: string
!     character(*), intent(in) :: what
!     character(*), intent(in) :: with
!     character(:), allocatable :: nstring
!     !> Local variables
!     character(len(string) * 10) :: tstring
!     integer :: i

!     i = outlen
!     tstring = string
!     do
!         i = index(tstring, what)
!         if (i == 0) exit
!         if (i == 1) then
!             tstring = with // trim(tstring(1+len(what):))
!         else if (i == len_trim(tstring)) then
!             tstring = trim(tstring(:i-1)) // with
!         else
!             tstring = tstring(:i-1) // with // trim(tstring(i+len(what):))
!         end if
!     end do
!     nstring = trim(adjustl(tstring))
! end function replace

function replace2(string, what, with) result(nstring)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: string
    character(*), intent(in) :: what
    character(*), intent(in) :: with
    character(:), allocatable :: nstring
    !> Local variables
    character(len(string) * 10) :: tstring
    integer :: i
    integer :: j
    integer :: cnt
    integer :: start

    tstring = string
    start = 1
    cnt = 0
    do
        i = index(string(start:), what)
        if (i <= 0) exit
        i = i + start - 1
        j = i + (len(with) - len(what)) * cnt 
        if (j == 1) then
            tstring = with // tstring(len(what)+1:)
        elseif (j == len(string) - len(what) + 1) then
            tstring = tstring(:j-1) // with
        else 
            tstring = tstring(:j-1) // with // tstring(j+len(what):)
        end if
        cnt = cnt + 1
        start = i + 1
    end do
    nstring = trim(tstring)
end function replace2

function replace(string, what, with, outlen) result(nstring)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: outlen
    character(*), intent(in) :: string
    character(*), intent(in) :: what
    character(*), intent(in) :: with
    character(outlen) :: nstring
    !> Local variables
    integer :: i, nnr

    nstring = string
    nnr = len_trim(with)
    do
        i = index(trim(nstring), what)
        if (i == 0) return
        nstring = nstring(:i-1) // with(:nnr) // nstring(i+len(what):)
    end do
end function replace
!***************************************************************************
!
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
integer function CountCharInString(s, c) result(n)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*) :: s
    character(*) :: c
    !> Local variables
    integer :: i


    n = 0
    do i = 1, len_trim(s)
        if (s(i:i) == c) n = n + 1
    end do
end function

!***************************************************************************
!
! \brief       Return whether passed string s is (not) numeric
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
logical function is_numeric(s)
    use m_numeric_kinds
    implicit none
    !> In/out variables
    character(*),intent(in) :: s
    !> Local variables
    real(kind = dbl) :: x
    integer :: e

    read(s, *, iostat=e) x
    is_numeric = e == 0
end function is_numeric

logical function is_not_numeric(s)
    use m_numeric_kinds
    implicit none
    !> In/out variables
    character(*),intent(in) :: s
    !> Local variables
    logical, external :: is_numeric


    is_not_numeric = .not. is_numeric(s)
end function is_not_numeric

!***************************************************************************
!
! \brief       Returns whether a string has format dddd-dd-dd with d a digit
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
logical function is_date(s)
    use m_numeric_kinds
    implicit none
    !> In/out variables
    character(*),intent(in) :: s
    !> Local variables
    logical, external :: is_numeric

    is_date = len_trim(s) == 10 &
        .and. is_numeric(s(1:4)) &
        .and. is_numeric(s(6:7)) &
        .and. is_numeric(s(9:10)) &
        .and. s(5:5) == '-' .and. s(8:8) == '-'
end function is_date

!***************************************************************************
!
! \brief       Returns whether a string has format dd:dd with d a digit
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
logical function is_time(s)
    use m_numeric_kinds
    implicit none
    !> In/out variables
    character(*),intent(in) :: s
    !> Local variables
    logical, external :: is_numeric

    is_time = len_trim(s) == 5 &
        .and. is_numeric(s(1:2)) &
        .and. is_numeric(s(4:5)) &
        .and. s(3:3) == ':'
end function is_time

!***************************************************************************
!
! \brief       Returns whether two strings match except for passed wildcard
!              which is intended to be in s1 (the template)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
logical function strings_match(s1, s2, wcard)
    use m_numeric_kinds
    implicit none
    !> In/out variables
    character(*),intent(in) :: s1
    character(*),intent(in) :: s2
    character(*),intent(in) :: wcard
    !> Local variables
    integer :: i


    !> Check on strings length
    if (len_trim(s1) /= len_trim(s2)) then
        strings_match = .false.
        return
    end if

    !> Check on strings content
    do i = 1, len(s1)
        if (s1(i:i) == s2(i:i) .or. s1(i:i) == wcard(1:1)) then
            cycle
        else
            strings_match = .false.
            return
        end if
    end do
    strings_match = .true.
end function strings_match
!***************************************************************************
!
! \brief       Returns index of n-th occurrence of c in s. If s doesn't have
!              n times c, returns ierror
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
integer function strCharIndex(s, c, n) result(ix)
    use m_common_global_var
    implicit none
    !> In/out variables
    character(*),intent(in) :: s
    character(*),intent(in) :: c
    integer, intent(in) :: n
    !> Local variables
    integer :: i
    integer :: cnt

    ix = ierror
    cnt = 0
    do i = 1, len(s)
        if (s(i:i) == c) then 
            cnt = cnt + 1
            if (cnt == n) then 
                ix = i
                exit
            end if
        end if
    end do
end function strCharIndex



function strByInt(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
end function strByInt
