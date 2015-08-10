!***************************************************************************
! write_out_user_stats.f90
! ------------------------
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
! \brief       Writes insensitive variables statistics on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutUserStats(unt, string, N, AddHeader)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    integer, intent(in) :: N
    character(*), intent(in) :: string
    logical, intent(inout) :: AddHeader
    !> local variables
    integer :: j = 0
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum = ''


    if (AddHeader) then
        call AddUserStatsHeader()
        AddHeader = .false.
    end if

    call clearstr(dataline)
    !> add file info
    call AddDatum(dataline, string(1:len_trim(string)), separator)
    call WriteDatumInt(N, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    do j = 1, NumUserVar
        call WriteDatumFloat(UserStats%Mean(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(UserStats%StDev(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(UserStats%Skw(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(UserStats%Kur(j), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    write(unt, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutUserStats

!***************************************************************************
!
! \brief       Adds header on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddUserStatsHeader()
    use m_rp_global_var
    implicit none
    !> local variables
    character(LongOutstringLen) :: headerline
    integer :: j

    call clearstr(headerline)
    headerline = 'filename,date,time,DOY,n_samples'
    do j = 1, NumUserVar
        call strcat(headerline, ',mean(')
        call strcat(headerline, trim(adjustl(UserCol(j)%label)))
        call strcat(headerline, '),')
        call strcat(headerline, 'st_dev(')
        call strcat(headerline, trim(adjustl(UserCol(j)%label)))
        call strcat(headerline, '),')
        call strcat(headerline, 'skw(')
        call strcat(headerline, trim(adjustl(UserCol(j)%label)))
        call strcat(headerline, '),')
        call strcat(headerline, 'kur(')
        call strcat(headerline, trim(adjustl(UserCol(j)%label)))
        call strcat(headerline, ')')
    end do
    if(RPsetup%out_st(1)) write(u_user_st1, '(a)') headerline(1:len_trim(headerline))
    if(RPsetup%out_st(2)) write(u_user_st2, '(a)') headerline(1:len_trim(headerline))
    if(RPsetup%out_st(3)) write(u_user_st3, '(a)') headerline(1:len_trim(headerline))
    if(RPsetup%out_st(4)) write(u_user_st4, '(a)') headerline(1:len_trim(headerline))
    if(RPsetup%out_st(5)) write(u_user_st5, '(a)') headerline(1:len_trim(headerline))
    if(RPsetup%out_st(6)) write(u_user_st6, '(a)') headerline(1:len_trim(headerline))
    if(RPsetup%out_st(7)) write(u_user_st7, '(a)') headerline(1:len_trim(headerline))
end subroutine AddUserStatsHeader
