!***************************************************************************
! read_planar_fit_file.f90
! ------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
    integer :: start
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

        !> Skip remaining lines until beginning of fitting plane coefficients
        do
            read(udf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            if (index(dataline, 'WindSector') == 0) cycle
            exit
        end do

        !> Read fitting plane coefficients
        do sec = 1, PFSetup%num_sec
            call clearstr(strg)
            read(udf, '(a)', iostat = io_status) strg
            if (io_status /= 0) then
                Meth%rot = 'double_rotation'
                call ExceptionHandler(29)
                return
            end if
            start = index(strg, '-') + 6
            read(strg(start:), *) (PFb(i, sec), i = 1, 3)
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
