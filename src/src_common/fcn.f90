!***************************************************************************
! fcn.f90
! -------
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
! \brief       Calculate either residuals or the Jacobian matrix in a \n
!              non linear regression procedure
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine fcn(m, npar, FcnPar, fvec, fjac, iflag)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: m
    integer, intent(in) :: npar
    real(kind = dbl) , intent(in)    :: FcnPar(:)
    real(kind = dbl) , intent(out)   :: fjac(:,:)
    real(kind = dbl) , intent(inout) :: fvec(:)
    integer, intent(inout)    :: iflag
    !> Local variables
    real(kind = dbl) , external :: func
    real(kind = dbl) , external :: ridders_diff
    integer :: i
    real(kind = dbl)  :: err


	!> Fn = FcnPar(1), f_co = FcnPar(2)\n
	!> m = no. of cases, npar = no. of parameters (2)\n
    select case(TFShape(1:len_trim(TFShape)))
    case ('iir')
        if (iflag == 1) then
          fvec = yFit(1:m) - (FcnPar(1) * zFit(1:m) / (1d0 + (xFit(1:m)/FcnPar(2))**2))
        elseif (iflag == 2) then
          do i = 1, m
            fjac(i, 1) = - zFit(i) / (1d0 + (xFit(i)/FcnPar(2))**2)
            fjac(i, 2) = - (2d0 * FcnPar(1) * zFit(i) * FcnPar(2) * xFit(i)**2) / (FcnPar(2)**2 + xFit(i)**2)**2
          end do
        end if
    case('sigma')
        if (iflag == 1) then
          fvec = yFit(1:m) - (zFit(1:m) * dexp(-0.346574d0 * (xFit(1:m)/FcnPar(1))**2))
        elseif (iflag == 2) then
          do i = 1, m
            fjac(i, 1) = - (2d0 * 0.346574d0 * zFit(i) * xFit(i)**2 &
                       * zFit(i) * dexp(-0.346574d0 * (xFit(i)/FcnPar(1))**2) / FcnPar(1)**3)
          end do
        end if
    case ('exponential')
        ddum(1:m) = dexp(FcnPar(1) * xFit(1:m)**2 + FcnPar(2) * xFit(1:m) + FcnPar(3))
        if (iflag == 1) then
          fvec = yFit(1:m) - ddum(1:m)
        elseif (iflag == 2) then
          do i = 1, m
            fjac(i, 1) = - ddum(i) * xFit(i)**2
            fjac(i, 2) = - ddum(i) * xFit(i)
            fjac(i, 3) = - ddum(i)
          end do
        end if
    case ('hyperbole')
        ddum(1:m) = zFit(1:m) / (FcnPar(2) + xFit(1:m))
        if (iflag == 1) then
          fvec = yFit(1:m) - (FcnPar(1) * ddum(1:m) + 1d0) * zzFit(1:m)
        elseif (iflag == 2) then
          do i = 1, m
            fjac(i, 1) = - ddum(i) * zzFit(i)
            fjac(i, 2) = FcnPar(1) * ddum(i) * zzFit(i) / (FcnPar(2) + xFit(i))
          end do
        end if
    case ('power')
        if (iflag == 1) then
          fvec = yFit(1:m) - (FcnPar(1) * xFit(1:m) ** FcnPar(2))
        elseif (iflag == 2) then
          do i = 1, m
            fjac(i, 1) = - (xFit(i) ** FcnPar(2))
            fjac(i, 2) = - (FcnPar(1) * xFit(i) ** FcnPar(2)) * dlog(xFit(i))
          end do
        end if
    case ('offset_power')
        if (iflag == 1) then
          fvec = yFit(1:m) - (FcnPar(1) + FcnPar(2) * xFit(1:m) ** FcnPar(3))
        elseif (iflag == 2) then
          do i = 1, m
            fjac(i, 1) = - 1d0
            fjac(i, 2) = - (xFit(i) ** FcnPar(3))
            fjac(i, 3) = - (FcnPar(2) * xFit(i) ** FcnPar(3)) * dlog(xFit(i))
          end do
        end if
    case ('cospectra_massman')
        if (iflag == 1) then
            do i = 1, m
                fvec(i) = yFit(i) - func(xFit(i), FcnPar(1), FcnPar(2), FcnPar(3))
            end do
        else
            do i = 1, m
                fjac(i, 1) = -1d0 * ridders_diff(func, xFit(i), FcnPar(1), FcnPar(2), FcnPar(3), 2, 0.1d0, err)
                fjac(i, 2) = -1d0 * ridders_diff(func, xFit(i), FcnPar(1), FcnPar(2), FcnPar(3), 3, 0.1d0, err)
                fjac(i, 3) = -1d0 * ridders_diff(func, xFit(i), FcnPar(1), FcnPar(2), FcnPar(3), 4, 0.1d0, err)
            end do
        end if
    end select
    return
end subroutine fcn
