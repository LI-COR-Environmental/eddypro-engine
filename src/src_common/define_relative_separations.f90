!***************************************************************************
! define_relative_separations.f90
! -------------------------------
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
! \brief       Define relative sensor separations starting from absolute ones
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineRelativeSeparations()
    use m_common_global_var
    !> in/out variables
    !> local variables
    integer :: gas


    do gas = co2, gas4
        if (E2Col(gas)%present) then
            !> Separations relative to selected anemometer
            E2Col(gas)%instr%nsep = E2Col(gas)%instr%nsep - E2Col(u)%instr%nsep
            E2Col(gas)%instr%esep = E2Col(gas)%instr%esep - E2Col(u)%instr%esep
            E2Col(gas)%instr%vsep = E2Col(gas)%instr%vsep - E2Col(u)%instr%vsep
            E2Col(gas)%instr%hsep = dsqrt(E2Col(gas)%instr%nsep**2 + E2Col(gas)%instr%esep**2)
            !> Absolute height of instruments
            E2Col(gas)%instr%height = E2Col(u)%instr%height + E2Col(gas)%instr%vsep
        end if
    end do
end subroutine
