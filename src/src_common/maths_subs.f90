!***************************************************************************
! maths_subs.f90
! --------------
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
! \brief       Collection of math functions
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!
! \brief       Calculate polynomial value, given coefficients and scalar x
!              for a polynom of 6th degree (or less, if used appropriately)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PolyVal(coeff, deg, x, N, polval)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: N
    integer, intent(in) :: deg
    real(kind = dbl), intent(in) :: coeff(0:deg)
    real(kind = dbl), intent(in) :: x(N)
    real(kind = dbl) :: polval(N)
    !> Local variables
    integer :: i
    integer :: j


    do i = 1, N
        if (x(i) /= error) then
            polval(i) = 0d0
            do j = 0, deg
                polval(i) = polval(i) + coeff(j) * x(i)**j
            end do
        else
            polval(i) = error
        end if
    end do
end subroutine PolyVal

!***************************************************************************
!
! \brief       Calculate sinc(x) = sin(pi*x)/(pi*x).
!              Handle singularity at x = 0
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
function sinc(x, N)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: N
    real(kind = dbl) , intent(in) :: x(N)
    real(kind = dbl)  :: sinc(N)


    where(x(:) /= 0d0)
        sinc(:) = sin(p*x(:))/ (p*x(:))
    elsewhere
        sinc(:) = error
    end where
end function

!***************************************************************************
!
! \brief       Apply fixed gains/offsets to an array. Gains/offsets can be
!              specified differently for positive and negative array
!              values
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
function asymmetric_linear_transformation(x, N, pgain, poffset, ngain, noffset)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: N
    real(kind = dbl) , intent(in) :: x(N)
    real(kind = dbl) , intent(in) :: pgain
    real(kind = dbl) , intent(in) :: poffset
    real(kind = dbl) , intent(in) :: ngain
    real(kind = dbl) , intent(in) :: noffset
    real(kind = dbl)  :: asymmetric_linear_transformation(N)

    where (x(:) /= error)
        where (x(:) >= 0d0)
            asymmetric_linear_transformation(:) = x(:) * pgain + poffset
        elsewhere
            asymmetric_linear_transformation(:) = x(:) * ngain + noffset
        end where
    elsewhere
        asymmetric_linear_transformation(:) = error
    end where
end function
