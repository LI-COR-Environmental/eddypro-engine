!***************************************************************************
! tilt_correction.f90
! -------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2014, LI-COR Biosciences
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
! \brief       Apply the selected method of rotation for tilt correction
! \author      Gerardo Fratini
! \note
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TiltCorrection(RotMeth, GoPlanarFit, Set, nrow, ncol, nsec, DirYaw, DirPitch, DirRoll, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: nsec
    logical, intent(in) :: GoPlanarFit(nsec)
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    character(*), intent(in) :: RotMeth
    real(kind = dbl), intent(out) :: DirYaw
    real(kind = dbl), intent(out) :: DirPitch
    real(kind = dbl), intent(out) :: DirRoll
    !> Local variables
    integer :: WindSec


    if (printout) write(*, '(a)') '  Performing tilt correction by: ' &
        // RotMeth(1:len_trim(RotMeth))

    !> If any wind component is error code, set the others to error too, for the sake of a proper
    !> application of the rotation methods
    where (Set(1:nrow, u) == error .or. Set(1:nrow, v) == error .or. Set(1:nrow, w) == error)
        Set(1:nrow, u) = error
        Set(1:nrow, v) = error
        Set(1:nrow, w) = error
    end where

    !> Initializations
    Essentials%yaw = error
    Essentials%pitch = error
    Essentials%roll = error
    select case(RotMeth(1:len_trim(RotMeth)))
        case ('double_rotation')
            call DoubleRotation(Set, size(Set, 1), size(Set, 2), DirYaw, DirPitch)
        case ('triple_rotation')
            call DoubleRotation(Set, size(Set, 1), size(Set, 2), DirYaw, DirPitch)
            call ThirdRotation(Set, size(Set, 1), size(Set, 2), DirRoll)
        case ('planar_fit', 'planar_fit_no_bias')
            WindSec = nint(error)
            call WindSector(Stats%wind_dir, WindSec)
            if (WindSec /= error) then
                if (GoPlanarFit(WindSec)) then
                    call PlanarFitByWindSector(WindSec, Set, size(Set, 1), size(Set, 2), DirYaw)
                else
                    call DoubleRotation(Set, size(Set, 1), size(Set, 2), DirYaw, DirPitch)
                end if
            else
                call DoubleRotation(Set, size(Set, 1), size(Set, 2), DirYaw, DirPitch)
            end if
        case ('not')
            continue
    end select
    if (printout) write(*, '(a)') '  Done.'
end subroutine TiltCorrection

!***************************************************************************
!
! \brief       Double rotation, based on Wiczak et al. (2001, BLM, eq. 22-29)
! \author      Gerardo Fratini
! \note
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DoubleRotation(Set, nrow, ncol, DirYaw, DirPitch)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: DirYaw
    real(kind = dbl), intent(out) :: DirPitch
    !> local variables
    integer :: i = 0
    real(kind = dbl) :: Mean(ncol)
    real(kind = dbl) :: SinTheta
    real(kind = dbl) :: CosTheta
    real(kind = dbl) :: SinPhi
    real(kind = dbl) :: Yaw(3, 3)
    real(kind = dbl) :: Pitch(3, 3)

    call AverageNoError(Set, size(Set, 1), size(Set, 2), Mean, error)

    !> yaw angle (note that cos(theta) = u / sqrt(u^2+v^2)
    !> gives the same theta as tan(theta) = v / u)
    SinTheta = Mean(v) / dsqrt(Mean(u)**2 + Mean(v)**2)
    CosTheta = Mean(u) / dsqrt(Mean(u)**2 + Mean(v)**2)
    DirYaw = 180d0 * dacos(CosTheta) / p
    if (SinTheta < 0d0) then
        DirYaw = 360d0 - DirYaw
    end if
    !> yaw matrix
    call YawMtx(DirYaw, Yaw)
    !> apply yaw matrix
    do i = 1, nrow
        if(Set(i, u) /= error .and. Set(i, v) /= error .and. Set(i, w) /= error) then
            Set(i, u:w) = matmul(Yaw, Set(i, u:w))
        end if
    end do
    call AverageNoError(Set, size(Set, 1), size(Set, 2), Mean, error)

    !> pitch angle (note that sin(phi) = w / sqrt(u^2+w^2)
    !> gives the same phi as tan(phi) = w / u)
    SinPhi = Mean(w) / dsqrt(Mean(u)**2 + Mean(w)**2)
    DirPitch = 180d0 * dasin(SinPhi) / p
    !> pitch matrix
    call PitchMtx(DirPitch, Pitch)
    !> apply pitch matrix
    do i = 1, nrow
        if(Set(i, u) /= error .and. Set(i, v) /= error .and. Set(i, w) /= error) then
            Set(i, u:w) = matmul(Pitch,Set(i, u:w))
        end if
    end do
    call AverageNoError(Set, size(Set, 1), size(Set, 2), Mean, error)
end subroutine DoubleRotation

!***************************************************************************
!
! \brief       Applies third rotation (nullifying cross-stream stess)
!              Based on Wiczak et al. (2001, BLM, eq. 30-33)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ThirdRotation(Set, nrow, ncol, DirRoll)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: DirRoll
    !> local variables
    integer :: i = 0
    real(kind = dbl) :: Mean(ncol)
    real(kind = dbl) :: Cov(ncol, ncol)
    real(kind = dbl) :: Dumm
    real(kind = dbl) :: Dum1
    real(kind = dbl) :: Dum2
    real(kind = dbl) :: Roll(3, 3)

    call AverageNoError(Set, size(Set, 1), size(Set, 2), Mean, error)
    call CovarianceMatrixNoError(Set, size(Set, 1), size(Set, 2), Cov, error)

    !> Roll angle
    DirRoll = error / 180d0 * p
    Dum1 = 2.d0 * Cov(v, w)
    Dum2 = Cov(v, v) - Cov(w, w)
    if (Dum2 /= 0d0) then
        Dumm = Dum1 / Dum2
    else
        return
    end if
    DirRoll = (5d-1 * datan(Dumm)) * 180d0 / p

    !> Roll correction is applied only if Roll < 10 degrees
    if (abs(DirRoll) < 10d0) then
        !> Roll matrix
        call RollMtx(DirRoll, Roll)
        !> apply pitch matrix
        do i = 1, nrow
            if(Set(i, u) /= error .and. Set(i, v) /= error .and. Set(i, w) /= error) then
                Set(i, u:w) = matmul(Roll,Set(i, u:w))
            end if
        end do
    end if
end subroutine ThirdRotation

!***************************************************************************
!
! \brief       Select the wind sector, applies planar fit rotation matricx \n
!              PFMat appplies Yaw rotation.
!              Based on Wiczak et al. (2001, BLM, eq. 45-46)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PlanarFitByWindSector(WindSec, Set, nrow, ncol, DirYaw)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: WindSec
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    real(kind = dbl) :: Mean(ncol)
    real(kind = dbl) :: Cov(ncol, ncol)
    real(kind = dbl):: DirYaw
    real(kind = dbl) :: Mat(3, 3)
    real(kind = dbl):: Yaw(3, 3)
    real(kind = dbl) :: SinTheta
    real(kind = dbl) :: CosTheta

    call AverageNoError(Set, size(Set, 1), size(Set, 2), Mean, error)
    call CovarianceMatrixNoError(Set, size(Set, 1), size(Set, 2), Cov, error)

    !> define the sector-specific rotation matrix for the current wind sector
    if (WindSec /= nint(error)) then
        Mat(:,:) = PFMat(:, :, WindSec)
        !> apply planar fit matrix
        do i = 1, nrow
            if(Set(i, u) /= error .and. Set(i, v) /= error .and. Set(i, w) /= error) then
                Set(i, u:w) = matmul(Mat, Set(i, u:w))
            end if
        end do
        call AverageNoError(Set, size(Set, 1), size(Set, 2), Mean, error)
        !> new yaw matrix
        if (dabs(Mean(u)) > 0d0 .and. dabs(Mean(v)) > 0d0) then
            !> The 2 following commented lines use the defintion of gamma as from
            !> Wilczak et al. 2001, BLM, eq. 36
            !DirYaw =  datan2(LocStats%Mean(V), LocStats%Mean(U))
            !DirYaw = 180.d0 * DirYaw / p
            !> Here, instead, gamma (DirYaw) is defined as the "real" wind direction
            SinTheta = Mean(v) / dsqrt(Mean(u)**2 + Mean(v)**2)
            CosTheta = Mean(u) / dsqrt(Mean(u)**2 + Mean(v)**2)
            DirYaw = 180d0 * dacos(CosTheta) / p
            if (SinTheta < 0d0) then
                DirYaw = 360d0 - DirYaw
            end if
            call YawMtx(DirYaw, Yaw)
            !> apply yaw matrix
            do i = 1, nrow
                if(Set(i, u) /= error .and. Set(i, v) /= error .and. Set(i, w) /= error) then
                    Set(i, u:w) = matmul(Yaw, Set(i, u:w))
                end if
            end do
            !> Subtract b0 from individual values of vertical wind
            if (RPsetup%pf_subtract_b0) then
                where(Set(:, w) /= error)
                    Set(:, w) = Set(:, w) - PFb(1, WindSec)
                end where
            end if
        endif
    end if
end subroutine PlanarFitByWindSector

!***************************************************************************
!
! \brief       Build rotation matrix for Yaw angle
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine YawMtx(DirYaw, Yaw)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: DirYaw
    real(kind = dbl), intent(out) :: Yaw(3, 3)
    !> local variables
    real(kind = dbl) :: CosTheta
    real(kind = dbl) :: SinTheta

    CosTheta = dcos(p * DirYaw / 180d0)
    SinTheta = dsin(p * DirYaw / 180d0)

    Yaw(1, 1:3) = (/CosTheta , SinTheta, 0d0/)
    Yaw(2, 1:3) = (/-SinTheta, CosTheta, 0d0/)
    Yaw(3, 1:3) = (/0d0      ,      0d0, 1d0/)
end subroutine YawMtx

!***************************************************************************
!
! \brief       Build rotation matrix for Pitch angle
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PitchMtx(DirPitch, Pitch)
    use m_rp_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: DirPitch
    real(kind = dbl), intent(out) :: Pitch(3, 3)
    !> local variables
    real(kind = dbl) :: CosPhi
    real(kind = dbl) :: SinPhi

    CosPhi = dcos(p * DirPitch / 180d0)
    SinPhi = dsin(p * DirPitch / 180d0)
    Pitch(1, 1:3) = (/ CosPhi, 0d0, SinPhi/)
    Pitch(2, 1:3) = (/    0d0, 1d0,    0d0/)
    Pitch(3, 1:3) = (/-SinPhi, 0d0, CosPhi/)
end subroutine PitchMtx

!***************************************************************************
!
! \brief       Build rotation matrix for Roll angle
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RollMtx(DirRoll, Roll)
    use m_rp_global_var
    implicit none
    ! in/out variables
    real(kind = dbl), intent(in) :: DirRoll
    real(kind = dbl), intent(out) :: Roll(3, 3)
    ! local variables
    real(kind = dbl) :: CosPsi
    real(kind = dbl) :: SinPsi

    CosPsi = dcos(p * DirRoll / 180d0)
    SinPsi = dsin(p * DirRoll / 180d0)
    Roll(1, 1:3) = (/1d0,     0d0,    0d0/)
    Roll(2, 1:3) = (/0d0,  CosPsi, SinPsi/)
    Roll(3, 1:3) = (/0d0, -SinPsi, CosPsi/)
end subroutine RollMtx
