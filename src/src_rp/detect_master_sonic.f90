!***************************************************************************
! detect_master_sonic.f90
! -----------------------
! Copyright (C) 2015-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
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
! \brief
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DetectMasterSonic(LocCol, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    integer, intent(in) :: ncol
    !> local variables
    integer :: i

    MasterSonic = NullInstrument
    do i = 1, ncol
        !> u-component of wind vector
        if (trim(adjustl(LocCol(i)%var)) == 'u' &
            .and. LocCol(i)%Instr%master_sonic) then
            MasterSonic = LocCol(i)%Instr
            exit
        end if
    end do

    !> Convenient variable that tells if sonic data is biased by the w-boost bug
    !> Note: If SwVer is not available for the sonic, and the sonic is a WM/WMP, then
    !> assume the bug is not present
    if ((MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'wm' .or. &
         MasterSonic%model(1:len_trim(MasterSonic%model) - 2) == 'wmpro') .and. &
         (MasterSonic%sw_ver%major == 2329 .and. MasterSonic%sw_ver%minor < 700)) then
         SonicDataHasWBug = .true.
    else
         SonicDataHasWBug = .false.
    end if
end subroutine DetectMasterSonic
