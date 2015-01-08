!***************************************************************************
! init_dynamic_metadata.f90
! -------------------------
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
!*******************************************************************************
!
! \brief       Read dynamic metadata file and figure out available parameters
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!*******************************************************************************
subroutine InitDynamicMetadata(N)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(out) :: N
    !> Local variables
    integer :: open_status
    integer :: io_status


    write(*, '(a)', advance = 'no') ' Initializing dynamic metadata usage..'

    !> Open file
    open(udf, file = AuxFile%DynMD, status = 'old', iostat = open_status)

    !> Interpret dynamic metadata file header and control in case of error
    if (open_status == 0) then
        call ReadDynamicMetadataHeader(udf)
    else
        call ExceptionHandler(68)
        EddyProProj%use_dynmd_file = .false.
    end if

    !> Count number of rows in the file (all of them, no matter if well formed),
    !> to give a maximum number to calibration data arrays
    N = 0
    countloop: do
        read(udf, *, iostat = io_status)
        if (io_status < 0 .or. io_status == 5001 .or. io_status == 5008) exit
        N = N + 1
    end do countloop
    close(udf)

    write(*, '(a)') ' Done.'
end subroutine InitDynamicMetadata

!***************************************************************************
!
! \brief       Reads and interprets file header, searching for knwon variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine ReadDynamicMetadataHeader(unt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    !> local variables
    character(1024) :: datastring
    character(64) :: Headerlabels(NumStdDynMDVars)
    integer :: read_status
    integer :: sepa
    integer :: cnt
    integer :: i
    integer :: j


    read(unt, '(a)', iostat = read_status) datastring
    cnt = 0
    do
        sepa = index(datastring, ',')
        if (sepa == 0) sepa = len_trim(datastring) + 1
        if (len_trim(datastring) == 0) exit
        cnt = cnt + 1
        Headerlabels(cnt) = datastring(1:sepa - 1)
        datastring = datastring(sepa + 1: len_trim(datastring))
    end do

    DynamicMetadataOrder = nint(error)
    do i = 1, cnt
        do j = 1, NumStdDynMDVars
        if(trim(adjustl(StdDynMDVars(j))) &
            == trim(adjustl(Headerlabels(i)))) then
            DynamicMetadataOrder(j) = i
            exit
        end if
        end do
    end do
end subroutine ReadDynamicMetadataHeader
