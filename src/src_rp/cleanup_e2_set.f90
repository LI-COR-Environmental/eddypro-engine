!***************************************************************************
! cleanup_e2_set.f90
! ------------------
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
! \brief       Set implausible values to error code in EddyPro main dataset \n
!              variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CleanUpE2Set(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: ncol
    integer, intent(in) :: nrow
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)


    !> General
    where(Set(:, :) < -300d0)
        Set(:, :) = error
    end where
    !> Specific for wind components
    where(Set(:, u:w) < -50d0 .or. Set(:, u:w) > 50d0)
        Set(:, u:w) = error
    end where
    !> Specific for temperatures
    where(Set(:, ts) < 200d0 .or. Set(:, ts) > 350d0)
        Set(:, ts) = error
    end where
    where(Set(:, ti1) < 200d0 .or. Set(:, ti1) > 350d0)
        Set(:, ti1) = error
    end where
    where(Set(:, ti2) < 200d0 .or. Set(:, ti2) > 350d0)
        Set(:, ti2) = error
    end where
    where(Set(:, te) < 200d0 .or. Set(:, te) > 350d0)
        Set(:, te) = error
    end where
    !> Specific for concentrations
    !> TBD
end subroutine CleanUpE2Set
