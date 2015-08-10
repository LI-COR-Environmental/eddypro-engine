!***************************************************************************
! retrieve_vars_selection.f90
! ---------------------------
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
! \brief       Based on selected variables stored in BypassCol, figure out \n
!              which columns should be used in the current GHG file \n
!              column selection for flux computation \n
!              join data creating the AllRaw dataset \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveVarsSelection(LocBypassCol, LocCol)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    type(ColType), intent(in) :: LocBypassCol(MaxNumCol)
    !> local variables
    integer ::i
    integer ::j

    LocCol%useit = .false.

    !> Associate the "flag" properties to columns selected as such by user
    do i = 1, NumRawFlags
        LocCol(RawFlag(i)%col)%flag = RawFlag(i)
    end do

    !> Outer cycle on reference column information
    ol: do i = 1, MaxNumCol
        !> Inner cycle on newly imported column information
        il: do j = 1, MaxNumCol
            !> Work only on columns with attribute "use it"
            if (LocBypassCol(i)%useit) then
                !> Look in the new columns info for a variable with same
                !> properties as the current reference one
                Select case (LocBypassCol(i)%measure_type)
                    case('molar_density', 'mole_fraction', 'mixing_ratio')
                        if (LocBypassCol(i)%var == LocCol(j)%var .and. &
                            LocBypassCol(i)%instr%firm == LocCol(j)%instr%firm .and. &
                            LocBypassCol(i)%instr%model == LocCol(j)%instr%model .and. &
                            LocBypassCol(i)%measure_type == LocCol(j)%measure_type) &
                            LocCol(j)%useit = .true.
                    case default
                        if (LocBypassCol(i)%var == LocCol(j)%var .and. &
                            LocBypassCol(i)%instr%firm == LocCol(j)%instr%firm .and. &
                            LocBypassCol(i)%instr%model == LocCol(j)%instr%model) then
                            LocCol(j)%useit = .true.
                            LocCol(j)%instr%model = LocBypassCol(i)%instr%model
                            LocCol(j)%instr%master_sonic = LocBypassCol(i)%instr%master_sonic
                        end if
                end select
            end if
        end do il
    end do ol
end subroutine RetrieveVarsSelection
