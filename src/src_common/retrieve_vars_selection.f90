!***************************************************************************
! retrieve_vars_selection.f90
! ---------------------------
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
