!***************************************************************************
! fix_timelag_opt_dataset.f90
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
! \brief       Eliminate error codes for easier following processing
!              Needs to create a new RH column for each gas
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixTimelagOptDataset(TimelagOpt, nrow, toSet, ton, actn, tlncol)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ton
    integer, intent(in) :: tlncol
    type(TimeLagOptType), intent(in):: TimelagOpt(nrow)
    type(TimeLagDatasetType), intent(out):: toSet(ton)
    integer, intent(out) :: actn(tlncol)
    !> Local variables
    integer :: i
    integer :: gas

    toSet = TimelagDatasetType(0d0, 0d0)
    actn = 0
    do i = 1, ton
        do gas = co2, gas4
            if(gas == h2o .and. TimelagOpt(i)%RH == error) cycle
            if(TimelagOpt(i)%tlag(gas) /= error) then
                actn(gas) = actn(gas) + 1
                toSet(actn(gas))%tlag(gas) = TimelagOpt(i)%tlag(gas)
                if (gas == h2o) toSet(actn(gas))%RH = TimelagOpt(i)%RH
            end if
        end do
    end do
end subroutine FixTimelagOptDataset

