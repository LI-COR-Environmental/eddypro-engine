!***************************************************************************
! write_ameriflux_output.f90
! --------------------------
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
! \brief       Writes outputs on AmeriFlux style
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteAmeriFluxOutput(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    type(ExType), intent(in) :: lEx
    !> local variables
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    type(datetype) :: now
    integer :: doy_now
    real(kind=dbl) :: dtime

    call DateTimeToDateType(lEx%date, lEx%time, now)
    doy_now = doy(now)
    dtime = dfloat(doy_now) + 0.02083d0 * (2d0 * now%hour)
    select case (now%minute)
        case(0:15)
            continue
        case(16:45)
            dtime = dtime +  0.02083d0
        case(46:59)
            dtime = dtime +  2d0 * 0.02083d0
    end select

    call clearstr(dataline)
    call AddDatum(dataline, lEx%date(1:4) // ',0', separator)
    write(datum, *) dtime
    call AddDatum(dataline, datum, separator)
    write(datum, *) doy_now
    call AddDatum(dataline, datum, separator)
    if(lEx%time(1:1) == '0') call AddDatum(dataline, lEx%time(2:2) // lEx%time(4:5), separator)
    if(lEx%time(1:1) /= '0') call AddDatum(dataline, lEx%time(1:2) // lEx%time(4:5), separator)
    write(datum, *) lEx%ustar
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%Ta - 273.16d0
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%WD
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%WS
    call AddDatum(dataline, datum, separator)
    if (lEx%Stor%of(co2) /= error .and. Flux3%co2 /= error) then
        write(datum, *) Flux3%co2 + lEx%Stor%of(co2)
        call AddDatum(dataline, datum, separator)
    else
        write(datum, *) aflx_error
        call AddDatum(dataline, datum, separator)
    end if
    write(datum, *) Flux3%co2
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%Stor%of(co2)
    call AddDatum(dataline, datum, separator)
    write(datum, *) Flux3%H
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%Stor%H
    call AddDatum(dataline, datum, separator)
    write(datum, *) Flux3%LE
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%Stor%LE
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
    write(datum, *) lEx%RH
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%Pa * 1d-3
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%chi(co2)
    call AddDatum(dataline, datum, separator)
    write(datum, *) lEx%VPD  * 1d-3
    call AddDatum(dataline, datum, separator)
    call AddDatum(dataline, '-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.', separator)
    write(datum, *) lEx%chi(h2o)
    call AddDatum(dataline, datum, separator)
    call AddDatum(dataline, '-6999.,-6999.,-6999.,-6999.,-6999.,-6999.,-6999.', separator)
    write(datum, *) lEx%zL
    call AddDatum(dataline, datum, separator)

    write(uaflx, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteAmeriFluxOutput

