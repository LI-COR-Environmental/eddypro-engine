!***************************************************************************
! air_and_cell_parameters.f90
! ---------------------------
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
! \brief       Calculate ambient and cell average T, P and molar volume
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AirAndCellParameters()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas

    !> Air temperature/pressure estimates
    !> Last true condition determines which temperature is used
    Stats%T = Stats%Mean(ts)
    if(Stats%Mean(te) > 220d0 .and. Stats%Mean(te) < 340d0) Stats%T = Stats%Mean(te)
    if(Stats%mT       > 220d0 .and. Stats%mT       < 340d0) Stats%T = Stats%mT

    !> Last true condition determines which pressure is used
    Stats%Pr = Metadata%bar_press
    if(Stats%Mean(pe) > 40000 .and. Stats%Mean(pe) < 110000) Stats%Pr = Stats%Mean(pe)
    if(Stats%mPr      > 40000 .and. Stats%mPr      < 110000) Stats%Pr = Stats%mPr

    !> Ambient air molar volume [m+3 mol-1] and air mass density [kg m-3]
    if (Stats%Pr > 0d0 .and. Stats%T /= error) then
        LitePar%Va = Ru * Stats%T / Stats%Pr
    else
        LitePar%Va = error
    end if

    !> Cell Temperature, if applicable \n
    if (Stats%Mean(tc) /= error) then
        LitePar%Tcell = Stats%Mean(tc)
    else
        LitePar%Tcell = Stats%T
    end if

    !> Cell pressure, if applicable \n
    if(Stats%Mean(pi) /= error) then
        LitePar%Pcell = Stats%Mean(pi)
    elseif (Stats%Pr /= error) then
        LitePar%Pcell = Stats%Pr
    else
        LitePar%Pcell = Metadata%bar_press
    end if

    !> Using cell temperature, for each gas column related to a closed-path analyser,
    !> determine cell air molar volume. If Tcell == Stats%T (that is, if no internal temperature
    !> is provided) and Pcell == Stats%Pr, then cell air molar volume is equal to ambient air molar volume
    E2Col%Va = error
    do gas = co2, gas4
        if (E2Col(gas)%instr%path_type == 'closed') then
            if (LitePar%Pcell > 0d0 .and. LitePar%Tcell /= error) then
                E2Col(gas)%Va = Ru * LitePar%Tcell / LitePar%Pcell
            else
                E2Col(gas)%Va = error
            end if
        end if
    end do
end subroutine AirAndCellParameters
