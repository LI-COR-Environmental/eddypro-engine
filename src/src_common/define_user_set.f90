!***************************************************************************
! define_user_set.f90
! -------------------
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
! \brief       Define "UserSet", the pre-defined set of variables  \n
!              needed for any following processing. \n
!              Variables are: u, v, w, ts, co2, h2o, ch4, gas4, tc, tc, \n
!              ti1, ti2, pi, te, pe
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineUserSet(LocCol, Raw, nrow, ncol, UserSet, unrow, uncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: unrow, uncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in) :: Raw(nrow, ncol)
    real(kind = dbl), intent(out) :: UserSet(unrow, uncol)
    !> local variables
    integer :: j
    integer :: jj
    character(len(LocCol%label)), external :: replace


    UserCol%var = 'none'
    UserCol%measure_type = 'none'
    UserSet = error
    jj = 0
    do j = 1, ncol
        select case (LocCol(j)%var(1:len_trim(LocCol(j)%var)))
            !> Sonic and irga variables without property "useit"
            case('co2','h2o','ch4','n2o', 'cell_t', 'int_t_1', 'int_t_2', &
                'int_p', 'air_t', 'air_p', 'u','v','w','ts','sos', &
                'flag_1', 'flag_2')
                if (.not. LocCol(j)%useit) then
                    jj = jj + 1
                    if (jj > uncol) exit
                    UserCol(jj) = LocCol(j)
                    UserCol(jj)%present = .true.
                    UserSet(1:unrow, jj) = Raw(1:unrow, j)
                end if
            case default
                !> Variables with a custom label
                if (.not. LocCol(j)%useit) then
                    jj = jj + 1
                    if (jj > uncol) exit
                    UserCol(jj) = LocCol(j)
                    UserCol(jj)%present = .true.
                    UserSet(1:unrow, jj) = Raw(1:unrow, j)
                    !> Replace spaces with underscores
                    UserCol(jj)%label = replace(UserCol(jj)%label, &
                        ' ', '_', len(UserCol(jj)%label))
                end if
                !> Special case of 4th gas calibration reference
                if (j == Gas4CalRefCol) UserCol(jj)%var = 'cal-ref'
        end select
    end do
end subroutine DefineUserSet
