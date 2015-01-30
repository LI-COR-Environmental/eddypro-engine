!***************************************************************************
! read_planar_fit_file.f90
! ------------------------
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
! \brief       Read planar fit file and import rotation matrices
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadPlanarFitFile()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: io_status
    integer :: open_status
    integer :: sec
    integer :: i
    integer :: zero
    integer :: j
    character(ShortInstringLen) :: dataline
    character(64) :: strg


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading planar-fit file: ' // AuxFile%pf(1:len_trim(AuxFile%pf))
    open(udf, file = AuxFile%pf, status = 'old', iostat = open_status)

    if (open_status == 0) then
        write(*, '(a)') ' planar fit file found, reading rotation matrices..'

        !> Read number of sectors from relevant line
        do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            if (index(dataline, 'Number_of_selected_wind_sectors') == 0) cycle
            read(dataline(index(dataline, ':') + 1: len_trim(dataline)), *) PFSetup%num_sec
            exit
        end do

        !> Skip remaining lines until beginning of rotation matrices
        do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            if (index(dataline, 'Rotation matrices') == 0) cycle
            exit
        end do

        !> Read rotation matrices
        do sec = 1, PFSetup%num_sec
            call clearstr(strg)
            read(udf, '(a)', iostat = io_status) strg
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if

            if(sec == 1) then
                !> From first sector retrieve north offset
                read(strg(26:30), '(i5)', iostat = io_status) zero
                if (zero >= 0) then
                    PFSetup%north_offset = zero
                else
                    PFSetup%north_offset = zero - 360d0
                end if
            end if

            !> Read end-of-sector angle
            read(strg(32:35), '(i4)', iostat = io_status) PFSetup%wsect_end(sec)
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            !> Determine whether sector is to be excluded
            PFSetup%wsect_exclude(sec) = index(strg,'excluded') /= 0

            inloop: do i = 1, 3
                read(udf, *, iostat = io_status) (PFMat(i, j, sec), j = 1, 3)
                if (io_status /= 0) then
                    Meth%rot = 'double_rotation'
                    call ExceptionHandler(29)
                    return
                end if
            enddo inloop
        end do
        write(LogInteger, '(i6)') PFSetup%num_sec
        write(*,'(a)') '  ' // trim(adjustl(LogInteger)) // ' sector(s) found.'
    else
       !> If the specified planar-fit file is not found or is empty, switches to double rotations
        Meth%rot = 'double_rotation'
        call ExceptionHandler(30)
    end if

    write(*,'(a)')   ' Done.'
end subroutine ReadPlanarFitFile
