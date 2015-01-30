!***************************************************************************
! aoa_cal_nakai_2012.f90
! ----------------------
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
! \brief       Working subroutine in Nakai's procedure, Nakai and Shimoyama (2012, AFM)
! \author      Gerardo Fratini
! \notes       This subroutine is a port from the original C code developed by Taro Nakai
!              distributed under a Creative Commons license. As of August 2012, the original
!              code is available at:
!              https://sites.google.com/site/micrometeorologist/aoa2012.h?attredirects=0
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AoaSteffensen2012(VelU, VelV, VelW)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(inout) :: VelU
    real(kind = dbl), intent(inout) :: VelV
    real(kind = dbl), intent(inout) :: VelW
    !> local variables
    real(kind = dbl) :: ws, wd, aoa
    real(kind = dbl), external :: sinerr12
    real(kind = dbl), external :: coserr12


    !> Steffensen's linear iteration
    ws = dsqrt(VelU**2 + VelV**2)
    if (ws == 0d0) then
        if (VelW >= 0d0) aoa =  90d0
        if (VelW < 0d0)  aoa = -90d0
    else
        aoa = datan(VelW/ws) * 180d0 / p
    end if

    !> Wind direction
    if (ws == 0d0) then
        wd = 0d0
    else
        if (velV >= 0d0) then
            wd = 180d0 - acos(velU / ws) * 180d0 / p;
        else
            wd = 180d0 + acos(velU / ws) * 180d0 / p;
        end if
    end if

    !> Steffensen's method to find true AoA
    if (ws /= 0d0) then
        call Steffensen2012(aoa, wd, VelW/ws)
    else
        if (velW > 0d0) then
            aoa = -90d0
        else
            aoa = 90d0
        end if
    end if
    !> calibration
    if (coserr12(aoa, wd) /= 0d0) VelU = VelU / coserr12(aoa, wd)
    if (coserr12(aoa, wd) /= 0d0) VelV = VelV / coserr12(aoa, wd)
    if (sinerr12(aoa, wd) /= 0d0) VelW = VelW / sinerr12(aoa, wd)
end subroutine AoaSteffensen2012

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai and Shimoyama (2012, AFM)
! \author      Gerardo Fratini
! \notes       This subroutine is a port from the original C code developed by Taro Nakai
!              distributed under a Creative Commons license. As of August 2012, the original
!              code is available at:
!              https://sites.google.com/site/micrometeorologist/aoa2012.h?attredirects=0
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function sinerr12(aoa, wd)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(inout) :: wd
    real(kind = dbl), intent(inout) :: aoa
    real(kind = dbl), parameter :: a1 = -3.19818998552857d-10
    real(kind = dbl), parameter :: a2 = -2.69824417931343d-8
    real(kind = dbl), parameter :: a3 =  4.16728613218081d-6
    real(kind = dbl), parameter :: a4 =  4.85252964763967d-4
    real(kind = dbl), parameter :: a5 =  1.67354200080193d-2
    real(kind = dbl), parameter :: b1 =  5.92731123831391d-10
    real(kind = dbl), parameter :: b2 =  1.44129103378194d-7
    real(kind = dbl), parameter :: b3 =  1.20670183305798d-5
    real(kind = dbl), parameter :: b4 =  3.92584527104954d-4
    real(kind = dbl), parameter :: b5 =  3.82901759130896d-3
    real(kind = dbl) :: A_aoa
    real(kind = dbl) :: B_aoa
    real(kind = dbl) :: aoa_orig, wd_orig

    aoa_orig = aoa
    wd_orig = wd
    if (aoa > 0d0) then
        aoa = - aoa
        wd = wd + 180
    end if

    !> Calculate sine correction factor as from Eq. 12
    A_aoa = a1 * aoa**5 + a2 * aoa**4 + a3 * aoa**3 + a4 * aoa**2 + a5 * aoa + 1d0 !< Eq. 15
    B_aoa = b1 * aoa**5 + b2 * aoa**4 + b3 * aoa**3 + b4 * aoa**2 + b5 * aoa       !< Eq. 16
    sinerr12 = A_aoa - B_aoa * dsin(3d0 * wd * p / 180d0) !< Eq. 12

    aoa = aoa_orig
    wd = wd_orig
end function sinerr12

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai and Shimoyama (2012, AFM)
! \author      Gerardo Fratini
! \notes       This subroutine is a port from the original C code developed by Taro Nakai
!              distributed under a Creative Commons license. As of August 2012, the original
!              code is available at:
!              https://sites.google.com/site/micrometeorologist/aoa2012.h?attredirects=0
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function coserr12(aoa, wd)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(inout) :: wd
    real(kind = dbl), intent(inout) :: aoa
    !> local variables
    real(kind = dbl), parameter :: c1 = -1.20804470033571d-9
    real(kind = dbl), parameter :: c2 = -1.58051314507891d-7
    real(kind = dbl), parameter :: c3 = -4.95504975706944d-6
    real(kind = dbl), parameter :: c4 =  1.60799801968464d-5
    real(kind = dbl), parameter :: c5 =  1.28143810766839d-3
    real(kind = dbl), parameter :: d1 =  2.27154016448720d-9
    real(kind = dbl), parameter :: d2 =  3.85646200219364d-7
    real(kind = dbl), parameter :: d3 =  2.03402753902096d-5
    real(kind = dbl), parameter :: d4 =  3.94248403622007d-4
    real(kind = dbl), parameter :: d5 =  9.18428193641156d-4
    real(kind = dbl) :: C_aoa
    real(kind = dbl) :: D_aoa
    real(kind = dbl) :: aoa_orig, wd_orig

    aoa_orig = aoa
    wd_orig = wd
    if (aoa > 0d0) then
        aoa = - aoa
        wd = wd + 180
    end if
    if (aoa < -70d0) aoa = -70d0  !< see Eq. 14

    !> Calculate sine correction factor as from Eq. 12
    C_aoa = c1 * aoa**5 + c2 * aoa**4 + c3 * aoa**3 + c4 * aoa**2 + c5 * aoa + 1d0  !> Eq. 17
    D_aoa = d1 * aoa**5 + d2 * aoa**4 + d3 * aoa**3 + d4 * aoa**2 + d5 * aoa        !> Eq. 18
    coserr12 = C_aoa + D_aoa * dsin(3d0 * wd * p / 180d0) !< Eq. 13

    aoa = aoa_orig
    wd = wd_orig
end function coserr12

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai and Shimoyama (2012, AFM)
! \author      Gerardo Fratini
! \notes       This subroutine is a port from the original C code developed by Taro Nakai
!              distributed under a Creative Commons license. As of August 2012, the original
!              code is available at:
!              https://sites.google.com/site/micrometeorologist/aoa2012.h?attredirects=0
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Steffensen2012(aoa, wd, a)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: wd, a
    real(kind = dbl), intent(out) :: aoa
    real(kind = dbl) :: x0, x1, x2, x3
    real(kind = dbl) :: key
    real(kind = dbl), external :: gx12

    x0 = aoa
    do
        x1 = gx12(x0, wd, a)
        x2 = gx12(x1, wd, a)
        x3 = x2
        key = x2 - 2d0 * x1 + x0
        if (dabs(key) < 0.01d0) exit
        x3 = x0 - (x1 - x0)**2d0 / key
        x0 = x3
    end do
    aoa = x3
end subroutine Steffensen2012

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai and Shimoyama (2012, AFM)
! \author      Gerardo Fratini
! \notes       This subroutine is a port from the original C code developed by Taro Nakai
!              distributed under a Creative Commons license. As of August 2012, the original
!              code is available at:
!              https://sites.google.com/site/micrometeorologist/aoa2012.h?attredirects=0
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function gx12(xx, wd, a)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(inout) :: xx
    real(kind = dbl), intent(in) :: a
    real(kind = dbl) :: wd
    !> local variables
    real(kind = dbl), external :: sinerr12
    real(kind = dbl), external :: coserr12


    gx12 = datan(a * coserr12(xx, wd) / sinerr12(xx, wd)) * 180d0 / p  !< Eq. 20
end function gx12
