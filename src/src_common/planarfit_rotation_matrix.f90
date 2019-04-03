!***************************************************************************
! planarfit_rotation_matrix.f90
! -----------------------------
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
! \brief       Calculate planar fit rotation matrix, to rotate winds into \n
!              the mean-streamline coordinate system
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PlanarFitRotationMatrix(sec, PP)
    use m_common_global_var
    implicit none
    ! in/out variables
    integer, intent(in) :: sec
    real(kind = dbl), intent(out) :: PP(3, 3)
    ! local variables
    real(kind = dbl) :: aux
    real(kind = dbl) :: cos_a
    real(kind = dbl) :: cos_b
    real(kind = dbl) :: sin_a
    real(kind = dbl) :: sin_b


    !---------------------------------------------------------------------------------
    !> \warning: element PFb(1) is coefficient b0, PFb(2) is b1 and PFb(3) is b2
    !>           of the equation w = b0 + u * b1 + v * b2 (Wilczak et al. 2001, eq. 39)
    !---------------------------------------------------------------------------------

    !> Define the last row of PP (See Wilczak et al. 2001, eqs. 42)
    aux = dsqrt(PFb(2, sec)**2 + PFb(3, sec)**2 + 1.d0)
    PP(3, 1) = - PFb(2, sec) / aux
    PP(3, 2) = - PFb(3, sec) / aux
    PP(3, 3) = 1.d0/ aux
    !> Define pitch (a) and roll (b) angles (See Wilczak et al. 2001, eqs. 44)
    aux = dsqrt(PP(3, 2)**2 + PP(3, 3)**2)
    sin_a = PP(3, 1)
    cos_a = aux
    sin_b = - PP(3, 2) / aux
    !sin_b = PFb(3, sec) / dsqrt(PFb(3, sec)**2 + 1d0)  !< equivalent alternative (van Dijk et al. 2004, eq. 3.41)
    cos_b =   PP(3, 3) / aux
    !cos_b = 1d0 / dsqrt(PFb(3, sec)**2 + 1d0)          !< equivalent alternative (van Dijk et al. 2004, eq. 3.41)

    !> Define remaining elements of PP (See Wilczak et al. 2001, eq. 36)
    PP(1, 1) =   cos_a
    PP(1, 2) =   sin_a * sin_b
    PP(1, 3) = - sin_a * cos_b
    PP(2, 1) =   0.d0
    PP(2, 2) =   cos_b
    PP(2, 3) =   sin_b
end subroutine PlanarFitRotationMatrix
