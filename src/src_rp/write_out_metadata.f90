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
! \brief       Write all results on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutMetadata(init_string)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    !> local variables
    integer :: gas
!    integer :: prof
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    include '../src_common/interfaces.inc'


    call clearstr(dataline)
    call AddDatum(dataline, init_string, separator)

    !> Site location and characteristics
    write(datum, *) Metadata%lat
    call AddDatum(dataline, datum, separator)
    write(datum, *) Metadata%lon
    call AddDatum(dataline, datum, separator)
    write(datum, *) Metadata%alt
    call AddDatum(dataline, datum, separator)
    write(datum, *) Metadata%canopy_height
    call AddDatum(dataline, datum, separator)
    write(datum, *) Metadata%d
    call AddDatum(dataline, datum, separator)
    write(datum, *) Metadata%z0
    call AddDatum(dataline, datum, separator)
    !> Acquisition setup
    write(datum, *) Metadata%file_length
    call AddDatum(dataline, datum, separator)
    write(datum, *) Metadata%ac_freq
    call AddDatum(dataline, datum, separator)
    !> Master sonic height and north offset
    write(datum, *) E2Col(u)%instr%firm(1:len_trim(E2Col(u)%Instr%firm))
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%model(1:len_trim(E2Col(u)%Instr%model))
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%height
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%wformat
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%wref
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%north_offset
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%hpath_length * 1d2
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%vpath_length  * 1d2
    call AddDatum(dataline, datum, separator)
    write(datum, *) E2Col(u)%instr%tau
    call AddDatum(dataline, datum, separator)
    !> irgas
    do gas = co2, gas4
        if (OutVarPresent(gas)) then
            write(datum, *) E2Col(gas)%instr%firm(1:len_trim(E2Col(gas)%Instr%firm))
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%model(1:len_trim(E2Col(gas)%Instr%model))
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%measure_type
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%nsep * 1d2
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%esep * 1d2
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%vsep * 1d2
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%tube_l * 1d2
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%tube_d * 1d3
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%tube_f * 6d4
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%kw
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%ko
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%hpath_length * 1d2
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%vpath_length * 1d2
            call AddDatum(dataline, datum, separator)
            write(datum, *) E2Col(gas)%instr%tau
            call AddDatum(dataline, datum, separator)
        end if
    end do

    write(umd, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutMetadata
