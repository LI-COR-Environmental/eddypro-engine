!***************************************************************************
! define_relative_separations.f90
! -------------------------------
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
