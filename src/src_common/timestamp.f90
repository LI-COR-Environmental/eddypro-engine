!>*****************************************************************************80
!> HMS_CURRENT_HMS returns the current HMS time as integers.
!>  Example:
!>    If the current time is 9:45:54.872 AM, then
!>    H = 9, M = 45, S = 54, MM = 872
!>
!>  Licensing:
!>    This code is distributed under the GNU LGPL license.
!>
!>  Modified: 26 February 2005
!>  Author: John Burkardt
!>  Patched: Gerardo Fratini, 21 July 2011
!>
!>  Parameters:
!>    Output, integer (kind = 4) H, M, S, MM, the current hour, minute,
!>    second, and thousandths of a second.
!***************************************************************************
subroutine hms_current_hms(h, m, s, mm)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    integer, intent(out) :: h
    integer, intent(out) :: mm
    integer, intent(out) :: m
    integer, intent(out) :: s
    !> local variables
    integer :: values(8)

    call date_and_time(values=values)
    h = values(5)
    m = values(6)
    s = values(7)
    mm = values(8)
    return
end subroutine hms_current_hms


!>***************************************************************************
!> HMS_CURRENT_PRINT prints the current HMS time, and a user specified string.
!>  Example:
!>     Wallclock:  9:45:54.872 AM  Started determinant calculation.
!>     Wallclock:  9:47:32.738 AM  Finished determinant calculation.
!>
!>  Licensing:
!>     This code is distributed under the GNU LGPL license.
!>
!>  Modified: 05 May 2003
!>  Author: John Burkardt
!>  Patched: Gerardo Fratini, 21 July 2011
!>
!>  Parameters:
!>     Input, character (len = *) STRING, the string to be printed.
!>***************************************************************************
subroutine hms_current_print(string1, string3, adv)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    character(*), intent(in) ::  string1
    character(*), intent(in) ::  string3
    logical, intent(in) ::  adv
    !> local variables
    character(23) ::string2

    call hms_current_string(string2)
    if (adv) then
        write (*, '(a,a,a)') string1, string2, trim(string3)
    else
        write (*, '(a,a,a)', advance = 'no') string1, string2, trim(string3)
    end if
    return
end

!>***************************************************************************
!>  HMS_CURRENT_STRING writes the current HMS data into a string.
!>  Example: STRING = ' 9:45:54.872 AM'
!>
!>  Licensing:
!>    This code is distributed under the GNU LGPL license.
!>
!>  Modified: 26 February 2005
!>  Author: John Burkardt
!>  Patched: Gerardo Fratini, 21 July 2011
!>
!>  Parameters:
!>    Output, character (len = *) STRING, contains the HMS information.
!>    A character length of 15 should always be sufficient.
!>***************************************************************************
subroutine hms_current_string(string)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    character(*), intent(out) ::  string
    !> local variables
    !character(2) :: ampm
    integer :: h
    integer :: mm
    integer :: n
    integer :: s
    integer :: year
    integer :: day
    integer :: month
    integer :: values(8)

    call date_and_time (values = values)
    year = values(1)
    month = values(2)
    day = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

!    if (h < 12) then
!        ampm = 'AM'
!    else if (h == 12) then
!        if (n == 0 .and. s == 0) then
!            ampm = 'Nn'
!        else
!            ampm = 'PM'
!        end if
!    else
!        h = h - 12
!        if (h < 12) then
!            ampm = 'PM'
!        else if (h == 12) then
!            if (n == 0 .and. s == 0) then
!                ampm = 'Md'
!            else
!                ampm = 'AM'
!            end if
!        end if
!    end if
    write (string, '(i4,a1,i2.2,a1,i2.2,a1,i2,a1,i2.2,a1,i2.2,a1,i3.3)') &
        year, '-', month, '-', day, ' ',h, ':', n, ':', s, '.', mm
    return
end subroutine hms_current_string

!>***************************************************************************
!> HMS_DELTA_PRINT prints the change in HMS time, and a user specified string.
!>  Example:
!>    Delta Wallclock:  0:01:37.966 AM  Determinant calculation.
!>
!>  Licensing:
!>    This code is distributed under the GNU LGPL license.
!>
!>  Modified: 06 May 2003
!>  Author: John Burkardt
!>  Patched: Gerardo Fratini, 22 July 2011
!>
!>  Parameters:
!>    Input, character (len = *) STRING, the string to be printed.
!>***************************************************************************
subroutine hms_delta_print(string1, string2)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    character(*), intent(in) ::  string1
    character(*), intent(in) ::  string2
    !> local variables
    integer :: h_del
    integer :: h_old
    integer :: m_del
    integer :: m_old
    integer :: mm_del
    integer :: mm_old
    integer ::s_del
    integer ::s_old
    integer, save :: s = 0
    integer, save :: mm = 0
    integer, save :: m = 0
    integer, save :: h = -1

    !> Back up the previous time.
    if (h == -1) then
        call hms_current_hms (h, m, s, mm)
        h_old = h
        m_old = m
        s_old = s
        mm_old = mm
    else
        h_old = h
        m_old = m
        s_old = s
        mm_old = mm
        call hms_current_hms (h, m, s, mm)
    end if

    h_del = h - h_old
    m_del = m - m_old
    s_del = s - s_old
    mm_del = mm - mm_old
    if (mm_del < 0) then
        s_del = s_del - 1
        mm_del = mm_del + 1000
    end if
    if (s_del < 0) then
        m_del = m_del - 1
        s_del = s_del + 60
    end if
    if (m_del < 0) then
        m_del = m_del + 60
        h_del = h_del - 1
    end if
    if (h_del < 0) then
        h_del = h_del + 24
    end if

    write (*, '(a,i2,a1,i2.2,a1,i2.2,a1,i3.3,2x,a)') &
        trim(string1), h_del, ':', m_del, ':', s_del, '.', mm_del, &
        trim (string2)
!    write (123, '(a,i2,a1,i2.2,a1,i2.2,a1,i3.3,2x,a)') &
!        trim(String1), h_del, ':', m_del, ':', s_del, '.', mm_del, &
!        trim (string2)
    return
end subroutine hms_delta_print

!>***************************************************************************
!> HMS_DELTA_PRINT returns the change in HMS time.
!>  Example:
!>    Delta Wallclock:  0:01:37.966 AM  Determinant calculation.
!>
!>  Licensing:
!>    This code is distributed under the GNU LGPL license.
!>
!>  Modified: 05 August 2011
!>  Author: Gerardo Fratini
!>
!>***************************************************************************
subroutine hms_delta(h_del, m_del, s_del, mm_del)
    use m_numeric_kinds
    implicit none
    !> i/o variables
    integer, intent(out) :: h_del
    integer, intent(out) :: m_del
    integer, intent(out) ::s_del
    integer, intent(out) :: mm_del
    !> local variables
    integer :: h_old
    integer :: m_old
    integer :: mm_old
    integer ::s_old
    integer, save :: s = 0
    integer, save :: mm = 0
    integer, save :: m = 0
    integer, save :: h = -1

    !> Back up the previous time.
    if (h == -1) then
        call hms_current_hms (h, m, s, mm)
        h_old = h
        m_old = m
        s_old = s
        mm_old = mm
    else
        h_old = h
        m_old = m
        s_old = s
        mm_old = mm
        call hms_current_hms (h, m, s, mm)
    end if

    h_del = h - h_old
    m_del = m - m_old
    s_del = s - s_old
    mm_del = mm - mm_old


    if (mm_del < 0) then
        s_del = s_del - 1
        mm_del = mm_del + 1000
    end if
    if (s_del < 0) then
        m_del = m_del - 1
        s_del = s_del + 60
    end if
    if (m_del < 0) then
        m_del = m_del + 60
        h_del = h_del - 1
    end if
    if (h_del < 0) then
        h_del = h_del + 24
    end if

end subroutine hms_delta


!>***************************************************************************
!>  TIMESTAMP prints the current YMDHMS date as a time stamp.
!>  Example:
!>    31 May 2001   9:45:54.872 AM
!>
!>  Licensing:
!>    This code is distributed under the GNU LGPL license.
!>
!>  Modified: 06 August 2005
!>  Author: John Burkardt
!>  Patched: Gerardo Fratini, 21 July 2011
!>
!>  Parameters:
!>    none.
!>***************************************************************************
subroutine timestamp()
    use m_numeric_kinds
    implicit none
    !> local variables
    integer :: d
    integer :: h
    integer :: m
    integer :: mm
    integer :: n
    integer :: s
    integer :: values(8)
    integer :: y
    character(8) :: ampm
    character(9), parameter, dimension(12) :: month = (/  &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)

    call date_and_time (values = values)
    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)
    if (h < 12) then
        ampm = 'AM'
    else if (h == 12) then
        if (n == 0 .and. s == 0) then
            ampm = 'Noon'
        else
            ampm = 'PM'
        end if
    else
        h = h - 12
        if (h < 12) then
            ampm = 'PM'
        else if (h == 12) then
            if (n == 0 .and. s == 0) then
                ampm = 'Midnight'
            else
                ampm = 'AM'
            end if
        end if
    end if

    write (*, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
        d, trim (month(m)), y, h, ':', n, ':', s, '.', mm, trim (ampm)
  return
end subroutine timestamp

!>***************************************************************************
!> TIMESTRING writes the current YMDHMS timestamp into a string.
!>  Example:
!>    STRING = '31 May 2001   9:45:54.872 AM'
!>
!>  Licensing:
!>    This code is distributed under the GNU LGPL license.
!>
!>  Modified: 06 August 2005
!>  Author: John Burkardt
!>  Patched: Gerardo Fratini, 21 July 2011
!>
!> Parameters:
!>    Output, character (len = *) STRING, contains the date information.
!>    A character length of 40 should always be sufficient.
!>***************************************************************************
subroutine timestring(string)
    use m_numeric_kinds
    implicit none
    !> in/out variables
    character(*), intent(inout) ::  string
    !> local variables
    integer :: d
    integer :: h
    integer :: m
    integer :: mm
    integer :: n
    integer :: s
    integer :: values(8)
    integer :: y
    character(8) :: ampm
    character(9), parameter, dimension(12) :: month = (/  &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)

    call date_and_time (values = values)

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if (h < 12) then
        ampm = 'AM'
    else if (h == 12) then
        if (n == 0 .and. s == 0) then
            ampm = 'Noon'
        else
            ampm = 'PM'
        end if
    else
        h = h - 12
        if (h < 12) then
            ampm = 'PM'
        else if (h == 12) then
            if (n == 0 .and. s == 0) then
                ampm = 'Midnight'
            else
                ampm = 'AM'
            end if
        end if
    end if

    write (string, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
        d, trim (month(m)), y, h, ':', n, ':', s, '.', mm, trim (ampm)
  return
end subroutine timestring
