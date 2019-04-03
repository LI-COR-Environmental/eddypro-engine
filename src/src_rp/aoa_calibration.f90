!***************************************************************************
! aoa_calibration_nakai.f90
! -------------------------
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
! \brief       Applies angle-of-attack correction (or flow distorsion correction, \n
!              or head correction) according to Nakai et al. (2006, AFM)
! \author      Taro Nakai, Gerardo Fratini
! \notes       This subroutine is taken from Taro Nakai's web page at:
!              http://todomatsu.lowtem.hokudai.ac.jp/~taro/download/dlcount.php?fname=aoa.f
!              Patched 04.2009 by G. Fratini to be included in ECO2S/EddyPro.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AoaCalibration(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    real(kind = dbl) :: wind(3)
    include '../src_common/interfaces.inc'


    !> No AoA to apply
    if (RPsetup%calib_aoa == 'none') return

    !> Nakai et al. 2006
    if (RPsetup%calib_aoa == 'nakai_06') then
        !> Apply Nakai et al. 2006
        do i = 1, nrow
            wind(u:w) = Set(i, u:w)
            if (all(wind /= error)) then
                call AoaSteffensen(wind(u), wind(v), wind(w))
            else
                wind = error
            end if
            Set(i, u:w) = wind(u:w)
        end do
        return
    end if

    !> Nakai et al. 2012
    if (RPsetup%calib_aoa == 'nakai_12') then
        !> If WindMaster has w-boost fixed, need to remove the w-boost
        !> before AoA after Nakai and Shimoyama 2012 can be applied
        if (.not. SonicDataHasWBug) &
            call RemoveGillWmWBoost(Set, nrow, ncol)

        !>Now apply correction
        do i = 1, nrow
            wind(u:w) = Set(i, u:w)
            if (all(wind /= error)) then
                call AoaSteffensen2012(wind(u), wind(v), wind(w))
            else
                wind = error
            end if
            Set(i, u:w) = wind(u:w)
        end do
        return
    end if
end subroutine AoaCalibration

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai et al. (2006)
! \author      Taro Nakai, Gerardo Fratini
! \notes       This subroutine is taken from Taro Nakai's web page at:
!              http://todomatsu.lowtem.hokudai.ac.jp/~taro/download/dlcount.php?fname=aoa.f
!              Patched 04.2009 by G. Fratini to be included in ECO2S/EddyPro.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AoaSteffensen(VelU, VelV, VelW)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(inout) :: VelU    !< Original and corrected VelU
    real(kind = dbl), intent(inout) :: VelV    !< Original and corrected VelV
    real(kind = dbl), intent(inout) :: VelW    !< Original and corrected VelW
    !> local variables
    real(kind = dbl) :: ws, wd, aoa
    real(kind = dbl) :: sinerr, coserr

    !> Horizontal wind speed
    ws = dsqrt(VelU**2 + VelV**2)
    if (ws == 0d0) then
        if (VelW >= 0d0) aoa =  90d0
        if (VelW < 0d0)  aoa = -90d0
    else
        aoa = atan(VelW/ws) * 180d0 / p
    end if

    !> Wind direction
    if (VelU > 0d0) then
        wd = 180d0
    else if (VelV < 0d0) then
         wd = 360d0
    else
         wd = 0d0
    end if
    if (VelU == 0d0) then
        if (VelV > 0d0) wd =  90d0
        if (VelV < 0d0) wd = 270d0
    else
        wd = wd - (atan(VelV/VelU) * 180d0 / p)
    end if

    !> Steffensen's method to find true AoA
    if (ws /= 0d0) call Steffensen(aoa, wd, VelW/ws)

    !> sine error calculation
    call RetrieveSinErr(aoa, sinerr)
    sinerr = dsin(aoa * p / 180d0) * sinerr + 0.0195d0

    !> cosine error calculation
    call RetrieveCosErr(aoa, wd, coserr)

    !> result output
    VelU = VelU * (dcos(aoa * p / 180d0) / coserr)
    VelV = VelV * (dcos(aoa * p / 180d0) / coserr)
    if ( ws /= 0d0 ) then
        VelW = dtan(aoa * p / 180d0) * dsqrt(VelU**2 + VelV**2)
    else
        VelW = VelW * (dsin(aoa * p / 180d0) / sinerr)
    end if
end subroutine AoaSteffensen

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai et al. (2006)
! \author      Taro Nakai, Gerardo Fratini
! \notes       This subroutine is taken from Taro Nakai's web page at:
!              http://todomatsu.lowtem.hokudai.ac.jp/~taro/download/dlcount.php?fname=aoa.f
!              Patched 04.2009 by G. Fratini to be included in ECO2S/EddyPro.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Steffensen(aoa, wd, a)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: wd, a
    real(kind = dbl), intent(out) :: aoa
    real(kind = dbl) :: x0, x1, x2, x3
    real(kind = dbl) :: key

    x0 = aoa
    do
        call RetrieveGx(x0, wd, a, x1)
        call RetrieveGx(x1, wd, a, x2)
        x3 = x2
        key = x2 - 2d0 * x1 + x0
        if (dabs(key) < 0.01d0) exit
        x3 = x0 - (x1 - x0)**2d0 / key
        x0 = x3
    end do
    aoa = x3
end subroutine Steffensen

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai et al. (2006)
! \author      Taro Nakai, Gerardo Fratini
! \notes       This subroutine is taken from Taro Nakai's web page at:
!              http://todomatsu.lowtem.hokudai.ac.jp/~taro/download/dlcount.php?fname=aoa.f
!              Patched 04.2009 by G. Fratini to be included in ECO2S/EddyPro.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveSinErr(xx, sinerr)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: xx
    real(kind = dbl), intent(out) :: sinerr
    real(kind = dbl) :: pn(4), pp(4)

    data pn /0.428727148d0, 55.59348879d0, 0.222867784d0, 0.4882d0/
    data pp /0.570590482d0, 1610.881585d0, 0.111150653d0, 0.972080458d0/

    if (xx < 0d0) then
        sinerr =  pn(1) / (1d0 + pn(2) * exp(-pn(3) * (xx + 90d0))) + pn(4)
    else
        sinerr = -pp(1) / (1d0 + pp(2) * exp(-pp(3) * xx)) + pp(4)
    end if
end subroutine RetrieveSinErr

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai et al. (2006)
! \author      Taro Nakai, Gerardo Fratini
! \notes       This subroutine is taken from Taro Nakai's web page at:
!              http://todomatsu.lowtem.hokudai.ac.jp/~taro/download/dlcount.php?fname=aoa.f
!              Patched 04.2009 by G. Fratini to be included in ECO2S/EddyPro.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveCosErr(xx, wd, coserr)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: wd
    real(kind = dbl), intent(out) :: coserr
    real(kind = dbl), intent(inout) :: xx
    !> local variables
    real(kind = dbl) :: q(3)
    real(kind = dbl) :: r
    real(kind = dbl) :: x_orig, f_aoa
    data q  /  1.41546d-6, 8.51092d-4, 1.00672d0/
    data r  / 6.28032d0 /

    x_orig = xx
    if (xx < -70d0) xx = -70d0
    if (xx >  70d0) xx =  70d0
    f_aoa = q(1)*xx**3 + q(2)*xx**2 + q(3)*xx + r * sin(3d0 * wd * p/180d0)
    xx = x_orig
    if (xx < -70d0) f_aoa = -90d0 * (1d0 - (90d0 + xx) / 20d0) + (90d0 + xx) / 20d0 * f_aoa;
    if (xx >  70d0) f_aoa =  90d0 * (1d0 - (90d0 - xx) / 20d0) + (90d0 - xx) / 20d0 * f_aoa;
    coserr   = cos(f_aoa * p/180d0)
end subroutine RetrieveCosErr

!***************************************************************************
!
! \brief       Working subroutine in Nakai's procedure, Nakai et al. (2006)
! \author      Taro Nakai, Gerardo Fratini
! \notes       This subroutine is taken from Taro Nakai's web page at:
!              http://todomatsu.lowtem.hokudai.ac.jp/~taro/download/dlcount.php?fname=aoa.f
!              Patched 04.2009 by G. Fratini to be included in ECO2S/EddyPro.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveGx(xx, wd, a, gx)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(inout) :: xx
    real(kind = dbl), intent(in) :: wd
    real(kind = dbl), intent(in) :: a
    real(kind = dbl), intent(out) :: gx
    !> local variables
    real(kind = dbl) :: sinerr
    real(kind = dbl) :: coserr

    call RetrieveSinErr(xx, sinerr)
    call RetrieveCosErr(xx, wd, coserr)
    gx = datan((a * coserr - 0.0195d0) / sinerr / cos(xx * p / 180d0)) * 180d0 / p
end subroutine RetrieveGx
