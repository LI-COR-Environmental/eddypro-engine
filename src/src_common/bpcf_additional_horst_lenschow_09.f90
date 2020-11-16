!***************************************************************************
! bpcf_additional_horst_lenschow_09.f90
! -------------------------------------
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
! \brief       Calculate spectral correction factors for instruments
!              separation according to Horst and Lenschow (2009, BLM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CF_HorstLenschow09(lEx, LocSetup)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ExType), optional, intent(in) :: lEx
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    real(kind = dbl) :: n_mx, k_mx, r_x(GHGNumVar)
    real(kind = dbl) :: n_my, k_my, r_y(GHGNumVar)
    real(kind = dbl) :: n_mz(GHGNumVar), k_mz(GHGNumVar), r_z(GHGNumVar)
    real(kind = dbl) :: zL
    real(kind = dbl) :: Ax
    real(kind = dbl) :: Ay
    real(kind = dbl) :: alpha
    real(kind = dbl) :: direc
    real(kind = dbl) :: r
 
    integer :: igas
    integer :: gas

    !> Initialization
    ADDCF%of(co2:gas4) = 1d0

    if (lEx%Flux0%zL /= error) then
        !> normalized wavenumber corresponding to cospectrum peak
        !> in streamwise direction as a function of stability (Eq. 15)
        if (lEx%Flux0%zL <= -0.1d0) then
            n_mx = 0.07d0
        else
            n_mx = 2.31d0 - 2.24d0 / (1.015d0 + 0.15d0 * lEx%Flux0%zL)**2
        end if

        !> normalized wavenumber corresponding to cospectrum peak
        !> in cross-stream direction as a function of stability (Eq. 18)
        if (lEx%Flux0%zL <= -0.05d0) then
            n_my = 0.15d0
        else
            n_my = 2.43d0 - 2.28d0 / (1.01d0 + 0.2d0 * lEx%Flux0%zL)**2
        end if

        !> normalized wavenumber corresponding to cospectrum peak
        !> in vertical direction as a function of stability (Eqs. 29 and 30)
        do igas = ico2, igas4
            gas = igas + 3
            if (lEx%var_present(gas)) then
                if (lEx%instr(igas)%vsep >= 0) then
                    zL = lEx%Flux0%zL
                    if (zL <= 0.03d0) then
                        n_mz(gas) = 0.1d0
                    else
                        n_mz(gas) = 0.43d0 - 0.33d0 / (0.964d0 + 1.2d0 * zL)**2
                    end if
                else
                    if (lEx%Flux0%L /= 0 .and. lEx%Flux0%L /= error) then
                        zL = (lEx%instr(sonic)%height + lEx%instr(igas)%vsep - lEx%disp_height) / lEx%Flux0%L
                        if (zL <= -0.03d0) then
                            n_mz(gas) = 0.013d0
                        else
                            n_mz(gas) = 0.3d0 - 0.287d0 / (1.051d0 + 1.7d0 * zL)**2
                        end if
                    else
                        n_mz(gas) = error
                    end if
                end if
            else
                n_mz(gas) = error
            end if
        end do

        !> peak wavenumbers in the three directions
        !> Note the factor Ua/U = 1.1 for the streamwise wavenumber
        k_mx = 2d0 * p * n_mx / 1.1d0 / (lEx%instr(sonic)%height - lEx%disp_height)
        k_my = 2d0 * p * n_my / (lEx%instr(sonic)%height - lEx%disp_height)
        do igas = ico2, igas4
            gas = igas + 3
            if (n_mz(gas) /= error) then
                if (lEx%instr(igas)%vsep >= 0) then
                    k_mz(gas) = 2d0 * p * n_mz(gas) / (lEx%instr(sonic)%height - lEx%disp_height)
                else
                    k_mz(gas) = 2d0 * p * n_mz(gas) / (lEx%instr(sonic)%height + lEx%instr(igas)%vsep - lEx%disp_height)
                end if
            end if
        end do

        !> Distances in stream, cross and vertical directions
        do gas = co2, gas4
            if (lEx%var_present(gas)) then
                igas = gas - 3
                call direction(lEx%instr(igas)%nsep, - lEx%instr(igas)%esep, direc)
                r = dsqrt(lEx%instr(igas)%nsep**2 + lEx%instr(igas)%esep**2)
                alpha = (lEx%WD - direc - 180d0) * p / 180d0
                r_x(gas) = dabs(r * dcos(alpha))
                r_y(gas) = dabs(r * dsin(alpha))
                r_z(gas) = dabs(lEx%instr(igas)%vsep)
            else
                r_x(gas) = error
                r_y(gas) = error
                r_z(gas) = error
            end if
        end do

        !> correction factors
        if(LocSetup%SA%horst_lens09 == 'full') then
            do gas = co2, gas4
                if (k_mx /= error .and. r_x(gas) /= error) then
                    Ax = dexp(-k_mx * r_x(gas))
                else
                    Ax = error
                end if
                if (k_my /= error .and. r_y(gas) /= error) then
                    Ay = dexp(-(k_my * r_y(gas))**1.2)       ! Eq. 16 instead of 13 - After notification by M. Aubinet, Nov. 2020
                else
                    Ay = error
                end if
                !> Horizontal correction factor
                if (Ax /= error .and. Ay /= error) then
                    ADDCF%of(gas) = dexp(dsqrt(dlog(Ax)**2 + dlog(Ay)**2))
                elseif(Ay /= error) then
                    ADDCF%of(gas) = Ay
                elseif(Ax /= error) then
                    ADDCF%of(gas) = Ax
                else
                    ADDCF%of(gas) = 1d0
                end if
                !> Add vertical correction factor
                if (k_mz(gas) /= error .and. r_z(gas) /= error) &
                    ADDCF%of(gas) = ADDCF%of(gas) * dexp(k_mz(gas) * r_z(gas))
            end do
        elseif(LocSetup%SA%horst_lens09 == 'cross_and_vertical') then
            do gas = co2, gas4
                if (k_my /= error .and. r_y(gas) /= error .and. &
                    k_mz(gas) /= error .and. r_z(gas) /= error) then
                    ADDCF%of(gas) = dexp(k_my * r_y(gas)) * dexp(k_mz(gas) * r_z(gas))
                end if
            end do
        end if
    end if
end subroutine CF_HorstLenschow09

subroutine direction(rx, ry, direc)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: rx
    real(kind = dbl), intent(in) :: ry
    real(kind = dbl), intent(out) :: direc

    if (rx > 0.d0) then
        direc = 180.d0
    else if (ry < 0.d0) then
        direc = 360.d0
    else
        direc = 0.d0
    end if
    if (rx == 0.d0) then
        if (ry > 0.d0) direc = 90.d0
        if (ry < 0.d0) direc = 270.d0
    else
        direc = direc - (datan(ry / rx) * 180.d0 / p)
    end if
end subroutine direction
