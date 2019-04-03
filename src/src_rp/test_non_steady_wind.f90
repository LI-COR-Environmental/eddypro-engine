!***************************************************************************
! test_non_steady_wind.f90
! ------------------------
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
! \brief       Checks for non steadiness of horizontal wind \n
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestNonSteadyWind(Set, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    real(kind = dbl) :: sv1(3) = 0.d0
    real(kind = dbl) :: sv2(3) = 0.d0
    real(kind = dbl) :: bv(3) = 0.d0
    real(kind = dbl) :: stime = 0.d0
    real(kind = dbl) :: stime2 = 0.d0
    real(kind = dbl) :: Vmin(3) = 0.d0
    real(kind = dbl) :: Vmax(3) = 0.d0
    real(kind = dbl) :: RNv(2) = 0.d0
    real(kind = dbl) :: RNS = 0.d0
    real(kind = dbl) :: CosPhi = 0.d0
    real(kind = dbl) :: SinPhi = 0.d0
    real(kind = dbl) :: DirYaw = 0.d0
    real(kind = dbl) :: Yaw(3, 3) = 0.d0
    real(kind = dbl) :: Vel(N, 3), Vm(3) = 0.d0
    real(kind = dbl) :: Vmi(N, 3), TmpV(N, 3)


    write(*, '(a)', advance = 'no') '   Non-steady horizontal wind test..'

    !> Initializations
    IntHF%ns = 0
    do i = 1, N
        Vel(i, u) = Set(i, u)
        Vel(i, v) = Set(i, v)
        Vel(i, w) = Set(i, w)
    end do

    !> Mean values
    Vm = sum(Vel, dim = 1)
    Vm = Vm / dble(N)

    !> Rotations, to define alongwind and crosswind components
    !> Yaw angle
    SinPhi = Vm(v) / ((Vm(u) **2 + Vm(v) **2) **0.5d0)
    CosPhi = Vm(u) / ((Vm(u) **2 + Vm(v) **2) **0.5d0)
    DirYaw = 180.d0 * acos(CosPhi) / p
    if (SinPhi < 0.d0) DirYaw = 360.d0 - DirYaw
    call YawMtx(DirYaw, Yaw)
    TmpV = 0.d0
    do i = 1, N
        do j = U, V
            do k = U, W
                TmpV(i, j) = TmpV(i, j) + Yaw(j, k) * Vel(i, k)
            end do
        end do
        Vel(i, :) = TmpV(i, :)
    end do
    Vm = sum(Vel, dim = 1)
    Vm = Vm / dble(N)

    !> Linear regressions
    sv1 = 0.d0
    sv2 = 0.d0
    stime = 0.d0
    stime2 = 0.d0
    do i = 1, N
        sv1(:) = sv1(:) + (Vel(i, :) * (dble(i - 1)))
        sv2(:) = sv2(:) + Vel(i, :)
        stime = stime + (dble(i - 1))
        stime2 = stime2 + (dble(i - 1)) **2
    end do
    bv(:) = (sv1(:) - (sv2(:) * stime) / dble(N)) / &
            (stime2 - (stime * stime) / dble(N))

    !> Extreme values of the regression
    Vmin(:) = Vm(:) + bv(:) * ((dble(1 - 1)) - stime / dble(N))
    Vmax(:) = Vm(:) + bv(:) * ((dble(N - 1)) - stime / dble(N))

    !> Local mean
    do i = 1, N
        Vmi(i, :) = Vm(:) + bv(:) * ((dble(i - 1)) - stime / dble(N))
    end do
    !> Relative variations
    do j = 1, 2
        RNv(j) = abs((Vmax(j) - Vmin(j)) / Vm(U))
    end do
    RNS = abs(sqrt((Vmax(1) - Vmin(1)) **2 + (Vmax(2) - Vmin(2)) **2) / Vm(U))
    !> Hard-flagging if one of RNV or RNS is beyond specified thresholds
    Essentials%ns_s_rnv(1) = RNv(1)
    Essentials%ns_s_rnv(2) = RNv(2)
    Essentials%ns_s_rns = RNS
    if ((RNv(1) >= ns%hf_lim) .or. (RNv(2) >= ns%hf_lim) .or. (RNS >= ns%hf_lim)) then
        IntHF%ns = 1
    end if
    write(*,'(a)') ' Done.'
end subroutine TestNonSteadyWind
