!***************************************************************************
! write_out_planar_fit.f90
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
! \brief       Write planar fit results on output file \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutPlanarFit(NumElem, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: NumElem(N)
    integer, external :: CreateDir
    !> local variables
    integer :: open_status = 1
    integer :: i
    integer :: j
    integer :: zero
    integer :: sec


    !> Create output file
    PlanarFit_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
              // EddyProProj%id(1:len_trim(EddyProProj%id)) &
              // PlanarFit_FilePadding // Timestamp_FilePadding // TxtExt
    open(upf, file = PlanarFit_Path, iostat = open_status, encoding = 'utf-8')

    !> Write on output file sector-wise plane coefficients
    if (PFSetup%north_offset >= 0) then
        zero = nint(PFSetup%north_offset)
    else
        zero = nint(360 + PFSetup%north_offset)
    end if
    write(upf, '(a)') 'Planar_fit_results'
    write(upf, '(a, i3)') 'Number_of_selected_wind_sectors: ', PFSetup%num_sec
    write(upf, '(a, i4)') 'Minimum_number_of_data_per_wind_sector: ', PFSetup%min_per_sec
    write(upf, '(a, f5.2)') 'Maximum_average_vertical_wind_component_(m/s): ', PFSetup%w_max
    write(upf, '(a, f5.2)') 'Minimum_average_horizontal_wind_component_(m/s): ', PFSetup%u_min
    write(upf, '(a, a)') 'Beginning_of_planar_fit_determination_period: ', PFSetup%start_date
    write(upf, '(a, a)') 'End_of_planar_fit_determination_period: ', PFSetup%end_date
    write(upf, '(a)')

    write(upf, '(a)') 'Fitting planes coefficients'
    write(upf, '(a11, 3x, a11, 3x, 5(a14))') 'WindSector', 'Degrees', 'B0', 'B1', 'B2'
    write(upf, '(i11, 5x, i4, a1, i4, 3x, 3(f14.7))') 1, zero, '-', &
        PFSetup%wsect_end(1) + nint(PFSetup%north_offset), (PFb(i, 1), i = 1, 3)
    if (PFSetup%num_sec > 1) then
        do sec = 2, PFSetup%num_sec
            write(upf, '(i11, 5x, i4, a1, i4, 3x, 3(f14.7))') &
            sec, PFSetup%wsect_end(sec - 1) + nint(PFSetup%north_offset), '-', &
                PFSetup%wsect_end(sec) + nint(PFSetup%north_offset), (PFb(i, sec), i = 1, 3)
        end do
    end if
    write(upf, *)

    !> Write on output file sector-wise rotation matrix
    write(upf, '(a)') 'Rotation matrices'
    if (PFSetup%wsect_exclude(1)) then
        write(upf, '(a, i2, a, i4, a2, i4, a)') 'Sector number ', 1,'; angle: ', zero, '-', &
            PFSetup%wsect_end(1) + nint(PFSetup%north_offset), '  (excluded)'
    else
        write(upf, '(a, i2, a, i4, a2, i4, a, i4)') 'Sector number ', 1,'; angle: ', zero, '-', &
            PFSetup%wsect_end(1) + nint(PFSetup%north_offset), '  sector numerosity: ', NumElem(1)
    end if
    do j = 1, 3
        write(upf, '(3(f12.6, 1x))') PFMat(j, 1:3, 1)
    end do

    if (PFSetup%num_sec > 1) then
        do sec = 2, PFSetup%num_sec
            if (PFSetup%wsect_exclude(sec)) then
                write(upf, '(a, i2, a, i4, a2, i4, a)') 'Sector number ', sec, &
                '; angle: ', PFSetup%wsect_end(sec - 1), '-', &
                    PFSetup%wsect_end(sec) + nint(PFSetup%north_offset), '  (excluded)'
            else
                write(upf, '(a, i2, a, i4, a2, i4, a, i4)') 'Sector number ', sec, &
                '; angle: ', PFSetup%wsect_end(sec - 1) + nint(PFSetup%north_offset), '-', &
                    PFSetup%wsect_end(sec) + nint(PFSetup%north_offset), '  sector numerosity: ', NumElem(sec)
            end if
            do j = 1, 3
                write(upf, '(3(f12.6, 1x))') PFMat(j, 1:3, sec)
            end do
        end do
    end if
    close(udf)
    write(*,'(a)') '  Results written on file: ' &
        // PlanarFit_Path(1:len_trim(PlanarFit_Path))
end subroutine WriteOutPlanarFit
