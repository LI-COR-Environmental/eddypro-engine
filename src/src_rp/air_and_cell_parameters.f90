!***************************************************************************
! air_and_cell_parameters.f90
! ---------------------------
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
    if(Stats%Mean(te)  > 220d0 .and. Stats%Mean(te) < 340d0) Stats%T = Stats%Mean(te)
    if(biomet%val(bTa) > 220d0 .and. biomet%val(bTa) < 340d0) Stats%T = biomet%val(bTa)

    !> Last true condition determines which pressure is used
    Stats%Pr = Metadata%bar_press
    if(Stats%Mean(pe)  > 40000 .and. Stats%Mean(pe)  < 110000) Stats%Pr = Stats%Mean(pe)
    if(biomet%val(bPa) > 40000 .and. biomet%val(bPa) < 110000) Stats%Pr = biomet%val(bPa)

    !> Ambient air molar volume [m+3 mol-1] and air mass density [kg m-3]
    if (Stats%Pr > 0d0 .and. Stats%T /= error) then
        Ambient%Va = Ru * Stats%T / Stats%Pr
    else
        Ambient%Va = error
    end if

    !> Cell Temperature, if applicable \n
    if (Stats%Mean(tc) /= error) then
        Ambient%Tcell = Stats%Mean(tc)
    else
        Ambient%Tcell = Stats%T
    end if

    !> Cell pressure, if applicable \n
    if(Stats%Mean(pi) /= error) then
        Ambient%Pcell = Stats%Mean(pi)
    elseif (Stats%Pr /= error) then
        Ambient%Pcell = Stats%Pr
    else
        Ambient%Pcell = Metadata%bar_press
    end if

    !> Using cell temperature, for each gas column related to a closed-path analyser,
    !> determine cell air molar volume. If Tcell == Stats%T (that is, if no internal temperature
    !> is provided) and Pcell == Stats%Pr, then cell air molar volume is equal to ambient air molar volume
    E2Col%Va = error
    do gas = co2, gas4
        if (E2Col(gas)%instr%path_type == 'closed') then
            if (Ambient%Pcell > 0d0 .and. Ambient%Tcell /= error) then
                E2Col(gas)%Va = Ru * Ambient%Tcell / Ambient%Pcell
            else
                E2Col(gas)%Va = error
            end if
        end if
    end do
end subroutine AirAndCellParameters
