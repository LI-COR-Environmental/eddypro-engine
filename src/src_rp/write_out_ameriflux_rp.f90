!***************************************************************************
! write_out_ameriflux_rp.f90
! --------------------------
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
! \brief       Writes outputs on AmeriFlux style
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutAmeriFlux_rp(date, time)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    !> local variables
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    integer :: int_doy
    real(kind = dbl) :: float_doy


    call DateTimeToDOY(date, time, int_doy, float_doy)

    call clearstr(dataline)
    call AddDatum(dataline, date(1:4) // ',0', separator)
    write(datum, *) float_doy
    call AddDatum(dataline, datum, separator)
    write(datum, *) int_doy
    call AddDatum(dataline, datum, separator)
    if(time(1:1) == '0') call AddDatum(dataline, time(2:2) // time(4:5), separator)
    if(time(1:1) /= '0') call AddDatum(dataline, time(1:2) // time(4:5), separator)
    write(datum, *) Ambient%us
    call AddDatum(dataline, datum, separator)
    write(datum, *) Ambient%Ta - 273.15d0
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stats4%wind_dir
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stats5%Mean(u)
    call AddDatum(dataline, datum, separator)
    if (Stor%of(co2) /= error .and. Flux3%co2 /= error) then
        write(datum, *) Flux3%co2 + Stor%of(co2)
        call AddDatum(dataline, datum, separator)
    else
        write(datum, *) aflx_error
        call AddDatum(dataline, datum, separator)
    end if
    write(datum, *) Flux3%co2
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stor%of(co2)
    call AddDatum(dataline, datum, separator)
    write(datum, *) Flux3%H
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stor%H
    call AddDatum(dataline, datum, separator)
    write(datum, *) Flux3%LE
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stor%LE
    call AddDatum(dataline, datum, separator)
    write(datum, *) aflx_error
    call AddDatum(dataline, datum, separator)
    write(datum, *) aflx_error
    call AddDatum(dataline, datum, separator)
    write(datum, *) aflx_error
    call AddDatum(dataline, datum, separator)
    write(datum, *) aflx_error
    call AddDatum(dataline, datum, separator)
    write(datum, *) aflx_error
    call AddDatum(dataline, datum, separator)
    write(datum, *) aflx_error
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stats%RH
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stats%Pr * 1d-3
    call AddDatum(dataline, datum, separator)
    write(datum, *) Stats%chi(co2)
    call AddDatum(dataline, datum, separator)
    write(datum, *) Ambient%VPD  * 1d-3
    call AddDatum(dataline, datum, separator)
    call AddDatum(dataline, '-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.', separator)
    write(datum, *) Stats%chi(h2o)
    call AddDatum(dataline, datum, separator)
    call AddDatum(dataline, '-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.', separator)
    write(datum, *) Ambient%zL
    call AddDatum(dataline, datum, separator)
    write(uaflx,*) dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutAmeriFlux_rp
