!***************************************************************************
! read_biomet_file.f90
! --------------------
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
! \brief       Reads biomet data file and retrieve only necessary data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadBiometFile(BiometFile, bN, M, skip_file)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: M
    character(*), intent(in) :: BiometFile
    integer, intent(out) :: bN
    logical, intent(out) :: skip_file
    !> local variables
    integer :: i
    integer :: open_status
    integer :: io_status
    character(64) :: timestring
    character(32) :: caux

    skip_file = .false.
    open(udf, file = BiometFile, iostat = open_status)
    if (open_status /= 0) then
        skip_file = .true.
        return
    end if

    !> Skip header
    if (BiometMeta%nhead > 0) then
        do i = 1, BiometMeta%nhead
            read(udf, *)
        end do
    end if

    !> Read data file
    bN = 0
    BiometSet = error
    if (BiometMeta%BiometCol(1)%var == 'DATE' .and. BiometMeta%BiometCol(2)%var == 'TIME') then
        !> Read file with both date and time in first columns
        do
            bN = bN + 1
            if (len_trim(BiometMeta%data_label) == 0) then
                read(udf, *, iostat = io_status) dvec(bN), timestring, BiometSet(bN, 1:M)
            else
                read(udf, *, iostat = io_status) caux, dvec(bN), timestring, BiometSet(bN, 1:M)
            end if
            if (io_status < 0 .or. io_status == 5001) exit
            if (io_status == 5010) then
                bN = bN - 1
                cycle
            end if
            if (io_status > 0) then
                bN = bN - 1
                cycle
            end if
            tvec(bN) = timestring(1:5)
        end do
        bN = bN - 1
        dvec(1:bN) = dvec(1:bN)(1:10)
    elseif (BiometMeta%BiometCol(1)%var == 'DATE') then
        !> Read file with only date in first column
        do
            bN = bN + 1
            if (len_trim(BiometMeta%data_label) == 0) then
                read(udf, *, iostat = io_status) dvec(bN), BiometSet(bN, 1:M)
            else
                read(udf, *, iostat = io_status) caux, dvec(bN), BiometSet(bN, 1:M)
            end if
            if (io_status < 0 .or. io_status == 5001) exit
            if (io_status == 5010) then
                bN = bN - 1
                cycle
            end if
            if (io_status > 0) then
                bN = bN - 1
                cycle
            end if
        end do
        bN = bN - 1
        dvec(1:bN) = dvec(1:bN)(1:10)
        call BuildDvecTvecFromBiometFileName(BiometFile, .false., .true., bN)
    elseif (BiometMeta%BiometCol(1)%var == 'TIME') then
        !> Read file with only time in first column
        do
            bN = bN + 1
            if (len_trim(BiometMeta%data_label) == 0) then
                read(udf, *, iostat = io_status) tvec(bN), BiometSet(bN, 1:M)
            else
                read(udf, *, iostat = io_status) caux, tvec(bN), BiometSet(bN, 1:M)
            end if
            if (io_status < 0 .or. io_status == 5001) exit
            if (io_status == 5010) then
                bN = bN - 1
                cycle
            end if
            if (io_status > 0) then
                bN = bN - 1
                cycle
            end if
        end do
        bN = bN - 1
        tvec(1:bN) = tvec(1:bN)(1:5)
        call BuildDvecTvecFromBiometFileName(BiometFile, .true., .false., bN)
    else
        !> Read file without timestamp
        do
            bN = bN + 1
            if (len_trim(BiometMeta%data_label) == 0) then
                read(udf, *, iostat = io_status) BiometSet(bN, 1:M)
            else
                read(udf, *, iostat = io_status) caux, BiometSet(bN, 1:M)
            end if
            if (io_status < 0 .or. io_status == 5001) exit
            if (io_status == 5010) then
                bN = bN - 1
                cycle
            end if
            if (io_status > 0) then
                bN = bN - 1
                cycle
            end if
        end do
        bN = bN - 1
        call BuildDvecTvecFromBiometFileName(BiometFile, .true., .true., bN)
    end if
    close(udf)
    if (bN <= 0) skip_file = .true.
end subroutine ReadBiometFile

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
subroutine BuildDvecTvecFromBiometFileName(fname, dvecIsNeeded, tvecIsNeeded, bN)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bN
    character(*), intent(in) :: fname
    logical, intent(in) :: dvecIsNeeded
    logical, intent(in) :: tvecIsNeeded
    !> local variables
    integer :: i
    type(datetype) :: ld_date
    type(datetype) :: init_ld_date
    character(10) :: date
    character(5) :: time
    character(13) :: datestring

    datestring = fname(1:4) // fname(6:7) // fname(9:10) // '-' // fname(12:13) // fname(14:15)
    call DateStringToTimestamp(datestring, .true., init_ld_date)
    BiometMeta%step = 1d0 / BiometMeta%rate / 60d0 !< time step in minutes
    do i = 1, bN
        ld_date = init_ld_date + datetype(0, 0, 0, 0, floor(BiometMeta%step * i))
        call DateTypeToDateTime(ld_date, date, time)
        if (dvecIsNeeded) dvec(i) = date
        if (tvecIsNeeded) tvec(i) = time
    end do
end subroutine BuildDvecTvecFromBiometFileName

