!***************************************************************************
! filter_dataset_for_flags.f90
! ----------------------------
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
! \brief       Eliminate individual raw records corresponding to out-ranged flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilterDatasetForFlags(LocCol, Raw, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: Raw(nrow, ncol)
    !> local variables
    integer :: i
    integer :: j
    logical :: filtered(nrow)


    ! write(*, '(a)', advance='no') '  Filtering raw data for custom flags..'
    filtered = .false.
    !> External cycle on all columns
    do j = 1, ncol
        !> Detect if column is a flag
        if (LocCol(j)%flag%col > 0) then
            !> If column is a flag, filters accordingly
            if (LocCol(j)%flag%upper) then
                do i = 1, nrow
                    if ((.not. filtered(i)) .and. Raw(i, j) > LocCol(j)%flag%threshold) then
                        Raw(i, 1:ncol) = error
                        filtered(i) = .true.
                    end if
                end do
            else
                do i = 1, nrow
                    if ((.not. filtered(i)) .and. Raw(i, j) < LocCol(j)%flag%threshold) then
                        Raw(i, 1:ncol) = error
                        filtered(i) = .true.
                    end if
                end do
            end if
        end if
    end do
    Essentials%m_custom_flags = count(filtered)
    ! write(*, '(a)') '  Done.'
    if (trim(EddyProProj%ftype) == 'licor_ghg') then
        ! Native formats read quickly so progress report after
        ! one day is fine. ghg files read very slow so that
        ! a sign of life after each file is useful
        write(*, '(a, i6)') '   Number of records eliminated for custom flags: ',  &
            Essentials%m_custom_flags
    end if

end subroutine FilterDatasetForFlags
