!***************************************************************************
! write_out_metadata.f90
! ----------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutMetadata(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    character(16000) :: dataline

    !> local variables
    integer :: var
    integer :: i
    integer :: gas
    integer :: igas
    character(DatumLen) :: datum
    include '../src_common/interfaces_1.inc'


    call clearstr(dataline)
    !> Preliminary timestmap information
    ! write(datum, *) lEx%fname(1:len_trim(lEx%fname))
    ! call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%date(1:10)
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%time(1:5)
    call AddDatum(dataline, datum, separator)

    !> Site location and characteristics
    write(datum, *) lEx%lat
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%lon
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%alt
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%canopy_height
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%disp_height
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%rough_length
    call AddDatum(dataline, datum, separator)

    !> Acquisition setup
    write(datum, *) lEx%file_length
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%ac_freq
    call AddDatum(dataline, datum, separator)
    !> Master sonic height and north offset
    write(datum, *) lEx%instr(sonic)%firm(1:len_trim(lEx%instr(sonic)%firm))
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%model(1:len_trim(lEx%instr(sonic)%model))
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%height
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%wformat
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%wref
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%north_offset
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%hpath_length
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%vpath_length
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%instr(sonic)%tau
    call AddDatum(dataline, datum, separator)
    !> irgas
    do igas = ico2, igas4
        if (fcc_var_present(3 + igas)) then
            write(datum, *) lEx%instr(igas)%firm(1:len_trim(lEx%instr(igas)%firm))
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%model(1:len_trim(lEx%instr(igas)%model))
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%measure_type(3 + igas)
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%nsep
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%esep
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%vsep
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%tube_l
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%tube_d
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%tube_f
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%kw
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%ko
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%hpath_length
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%vpath_length
            call AddDatum(dataline, datum, separator)
            write(datum, *) lEx%instr(igas)%tau
            call AddDatum(dataline, datum, separator)
        end if
    end do

    write(umd,*) dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutMetadata
