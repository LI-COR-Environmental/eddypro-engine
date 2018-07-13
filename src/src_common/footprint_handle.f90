!***************************************************************************
! footprint_handle.f90
! --------------------
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
! \brief       Hub to and implementation of to several cross-wind \n
!              integrated footprint models
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FootprintHandle(var_w, ustar, zL, wind_speed, MO_length, sonic_height, &
        disp_height, rough_length)
    use m_common_global_var
    implicit none
    !> In/out variables
    real(kind = dbl), intent(in) :: var_w
    real(kind = dbl), intent(in) :: ustar
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: MO_length
    real(kind = dbl), intent(in) :: sonic_height
    real(kind = dbl), intent(in) :: disp_height
    real(kind = dbl), intent(in) :: rough_length
    !> Local variables
    real(kind = dbl) :: std_w


    if (foot_model_used == 'none') then
        Foot = errFootprint
        return
    end if

    !> If Kljun model was chosen, but conditions are outside those stated at Pag. 512 of the paper
    !> shift to Kormann and Meixner model.
    if (foot_model_used == 'kljun_04' .and. &
        (var_w <= 0d0 .or. ustar < kj_us_min .or. &
        zL < kj_zL_min .or. zL > kj_zL_max .or. sonic_height < 1d0)) then
        if (EddyProProj%fluxnet_mode) then
            Foot = errFootprint
            return
        else
            foot_model_used = 'kormann_meixner_01'
        end if
    end if

    !> DSZKP  mmzdxgc c: <-- code contribution by Luna Marie Fratini, Jul 2018
    !> Calculate std_w
    if (var_w >= 0d0) then
        std_w = dsqrt(var_w)
    else
        std_w = error
    end if

    select case(foot_model_used)
        case('kljun_04')
            call Kljun04(std_w, ustar, zL, sonic_height, disp_height, rough_length)
        case('kormann_meixner_01')
            call KormannMeixner01(ustar, zL, wind_speed, sonic_height, disp_height)
        case('hsieh_00')
            call Hsieh00(MO_length, sonic_height, disp_height, rough_length)
    end select
end subroutine FootprintHandle

!***************************************************************************
!
! \brief       Footprint esitmations based on Kljun et al. (2004, BLM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Kljun04(std_w, ustar, zL, sonic_height, disp_height, rough_length)
    use m_common_global_var
    implicit none
    !> In/out variables
    real(kind = dbl), intent(in) :: std_w
    real(kind = dbl), intent(in) :: ustar
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: sonic_height
    real(kind = dbl), intent(in) :: disp_height
    real(kind = dbl), intent(in) :: rough_length
    !> local variables
    real(kind = dbl) :: xstarmax
    real(kind = dbl) :: xstar
    real(kind = dbl) :: af, bb, ac, ad
    real(kind = dbl) :: a, b, c, d
    real(kind = dbl) :: zm
    real(kind = dbl) , external :: gammln

    real(kind = dbl) :: L(96)

    !> Values of L' for varying R (percentage of footprint) for b = 3.70 (see Fig. A1 in Kljun et al. 2004)
    data L(1:96) / 0.000000d0, 0.302000d0, 0.368000d0, 0.414000d0, 0.450000d0, 0.482000d0, 0.510000d0, &
        0.536000d0, 0.560000d0, 0.579999d0, 0.601999d0, 0.621999d0, 0.639999d0, 0.657998d0, 0.675998d0, &
        0.691998d0, 0.709998d0, 0.725998d0, 0.741997d0, 0.755997d0, 0.771997d0, 0.785997d0, 0.801997d0, &
        0.815996d0, 0.829996d0, 0.843996d0, 0.857996d0, 0.871996d0, 0.885995d0, 0.899995d0, 0.911995d0, &
        0.925995d0, 0.939995d0, 0.953995d0, 0.965994d0, 0.979994d0, 0.993994d0, 1.005994d0, 1.019994d0, &
        1.033994d0, 1.045993d0, 1.059993d0, 1.073993d0, 1.085993d0, 1.099993d0, 1.113993d0, 1.127992d0, &
        1.141992d0, 1.155992d0, 1.169992d0, 1.183992d0, 1.197991d0, 1.211991d0, 1.225991d0, 1.239991d0, &
        1.253991d0, 1.269991d0, 1.283990d0, 1.299990d0, 1.315990d0, 1.329990d0, 1.345990d0, 1.361989d0, &
        1.379989d0, 1.395989d0, 1.411989d0, 1.429989d0, 1.447988d0, 1.465988d0, 1.483988d0, 1.501988d0, &
        1.521987d0, 1.539987d0, 1.559987d0, 1.581987d0, 1.601986d0, 1.623986d0, 1.647986d0, 1.669985d0, &
        1.693985d0, 1.719985d0, 1.745984d0, 1.773984d0, 1.801984d0, 1.831983d0, 1.863983d0, 1.895983d0, &
        1.931982d0, 1.969982d0, 2.009982d0, 2.053984d0, 2.101986d0, 2.153988d0, 2.211991d0, 2.279994d0, 2.355998d0 /

    !> Initialization to error
    Foot = ErrFootprint

    !> Height above displacement height
    zm = sonic_height - disp_height

    !> Check on retrieved parameters
    if (zm < 1d0 .or. rough_length <= 0d0) return

    !> Calculate a, b, c, d, depending only on z0 (Eq. 13-16 in Kljun et al. 2004)
    af = 0.175d0
    bb = 3.418d0
    ac = 4.277d0
    ad = 1.685d0
    b  = 3.69895d0
    a = af / (bb - dlog(rough_length))
    c = ac * (bb - dlog(rough_length))
    d = ad * (bb - dlog(rough_length))

    if (std_w == error .or. ustar < kj_us_min .or. zL < kj_zL_min .or. zL > kj_zL_max) then
        return
    else
        !> Calculate location of peak influence
        xstarmax = c - d
        Foot%peak = xstarmax * zm *(std_w / ustar)**(-0.8d0)

        !> Calculate offset from tower: location of 1% contribution
        xstar = L(2) * c - d
        Foot%offset = xstar * zm *(std_w / ustar)**(-0.8d0)

        !> Calculate distances including increasing percentages of the footprint
        xstar = L(11) * c - d
        Foot%x10 = xstar * zm *(std_w / ustar)**(-0.8d0)
        xstar = L(31) * c - d
        Foot%x30 = xstar * zm *(std_w / ustar)**(-0.8d0)
        xstar = L(51) * c - d
        Foot%x50 = xstar * zm *(std_w / ustar)**(-0.8d0)
        xstar = L(71) * c - d
        Foot%x70 = xstar * zm *(std_w / ustar)**(-0.8d0)
        xstar = L(81) * c - d
        Foot%x80 = xstar * zm *(std_w / ustar)**(-0.8d0)
        xstar = L(91) * c - d
        Foot%x90 = xstar * zm *(std_w / ustar)**(-0.8d0)
    end if
end subroutine Kljun04

!***************************************************************************
!
! \brief       Footprint esitmations based on Kormann and Meixner, 2001
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine KormannMeixner01(ustar, zL, wind_speed, sonic_height, disp_height)
    use m_common_global_var
    implicit none
    !> In/out variables
    real(kind = dbl), intent(in) :: ustar
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: sonic_height
    real(kind = dbl), intent(in) :: disp_height
    !> local variables
    integer :: i
    real(kind = dbl) :: n
    real(kind = dbl) :: phi_m
    real(kind = dbl) :: phi_c
    real(kind = dbl) :: psi_m
    real(kind = dbl) :: eta
    real(kind = dbl) :: key
    real(kind = dbl) :: m
    real(kind = dbl) :: UU
    real(kind = dbl) :: r
    real(kind = dbl) :: mmu
    real(kind = dbl) :: zeta
    real(kind = dbl) :: zm
    real(kind = dbl) :: int_foot
    real(kind = dbl), parameter :: di = 1d0
    logical :: do_offset
    logical :: do10
    logical :: do30
    logical :: do50
    logical :: do70
    logical :: do80


    !> Initialization to error
    Foot = ErrFootprint

    zm = sonic_height - disp_height

    !> ALTERNATIVE u(z) for z=zm (Monin-Obukhov similarity profile)
    !wind_speed = ustar *(dlog(zm / rough_length) + psi_m) / vk

    !> Similarity relations (Paulson, 1970)
    if (zL > 0) then
        phi_m = 1d0 + 5d0 * zL
        phi_c = phi_m
        psi_m = - 5d0 * zL
    else
        phi_m = (1d0 - 16d0 * zL)**(-1d0 / 4d0)
        phi_c = (1d0 - 16d0 * zL)**(-1d0 / 2d0)
        eta = (1d0 - 16d0 * zL)**(1d0 / 4d0)
        psi_m = 2d0 * dlog((1d0 + eta) / 2d0) + dlog((1d0 + eta**2) / 2d0) - 2d0 * datan(eta) + p / 2d0  !< K&M2001
        !psi_m = 0.0954d0 - 1.86d0 * (zm/MO_length) - 1.07d0 * (zm/MO_length)**2 - 0.249 * (zm/MO_length)**3  !< Zhang & Anthes 1983, polynomial interpolation
        !psi_m = dlog(((1 + eta) / 2d0)**2 * (1d0 + eta**2) / 2d0) - 2d0 * datan(eta) + p / 2d0 !< alternative
    end if
    psi_m = -1d0 * psi_m  !< change sign to conform with K&M usage

    !> Intermediate parameters for K&M2001
    !> exponent of the diffusivity power law
    if (zL > 0) then
        n  = 1d0 / phi_m
    else
        n = (1d0 - 24d0 * zL) / (1d0 - 16d0 * zL)
    end if

    !> proportionality constant of the diffusivity power law (Eqs. 11 and 32)
    key  = vk * ustar * zm / (phi_c * zm**n)

    !> exponent of the wind speed power law
    m = ustar * phi_m / (vk * wind_speed)

    !> proportionality constant of the wind speed power law (Eqs. 11 and 31)
    !UU = ustar * (dlog(zm / rough_length) + psi_m) / (vk * zm**m)
    UU = wind_speed / zm**m

    !> Intermediate parameters
    r = 2d0 + m - n
    mmu = (1 + m) / r
    zeta = UU * zm**r / (r**2 * key)

    !> Footprint according to Kormann and Meixner, 2001
    do_offset = .true.
    do10 = .true.
    do30 = .true.
    do50 = .true.
    do70 = .true.
    do80 = .true.
    int_foot = 0d0

    do i = 1, 10000
        !> Cross-wind integrated 1D function
        int_foot = int_foot + di * (zeta**mmu * dexp(-zeta / (i * di)) / ((i * di)**(1 + mmu) * gamma(mmu)))
        if (do_offset .and. int_foot > 0.01d0) then
            Foot%offset = i * di
            do_offset = .false.
        end if
        if (do10 .and. int_foot > 0.1d0) then
            Foot%x10 = i * di
            do10 = .false.
        end if
        if (do30 .and. int_foot > 0.3d0) then
            Foot%x30 = i * di
            do30 = .false.
        end if
        if (do50 .and. int_foot > 0.5d0) then
            Foot%x50 = i * di
            do50 = .false.
        end if
        if (do70 .and. int_foot > 0.7d0) then
            Foot%x70 = i * di
            do70 = .false.
        end if
        if (do80 .and. int_foot > 0.8d0) then
            Foot%x80 = i * di
            do80 = .false.
        end if
        if (int_foot > 0.9d0) then
            Foot%x90 = i * di
            exit
        end if
    end do
    !> Peak value
    Foot%peak = zeta / (1d0 + mmu)
end subroutine KormannMeixner01

!***************************************************************************
!
! \brief       Footprint esitmations based on Hsieh et al. 2000
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Hsieh00(MO_length, sonic_height, disp_height, rough_length)
    use m_common_global_var
    implicit none
    !> In/out variables
    real(kind = dbl), intent(in) :: MO_length
    real(kind = dbl), intent(in) :: sonic_height
    real(kind = dbl), intent(in) :: disp_height
    real(kind = dbl), intent(in) :: rough_length
    !> local variables
    integer :: i
    real(kind = dbl) :: a1
    real(kind = dbl) :: DD
    real(kind = dbl) :: PP
    real(kind = dbl) :: fact
    real(kind = dbl) :: p1
    real(kind = dbl) :: zu
    real(kind = dbl) :: z0m
    real(kind = dbl) :: zm
    real(kind = dbl) :: zL
    real(kind = dbl) :: int_foot
    real(kind = dbl), parameter :: di = 5d0
    logical :: do_offset
    logical :: do10
    logical :: do30
    logical :: do50
    logical :: do70
    logical :: do80


    !> Initialization to error
    Foot = ErrFootprint

    zm = sonic_height - disp_height

    !> Model by Hsieh et al. (2000)
    !> Intermediate parameters
    a1 = 0.3d0
    p1 = 0.86d0
    z0m = rough_length
    zu = zm * (dlog(zm/rough_length) - 1d0 + rough_length/zm)
    zL = zu / MO_length

    !> Parameers D and P in Eq. 17
    DD = 0.97d0
    PP = 1d0
    if (dabs(zL) < 0.04d0) then
        !> Neutral and near neutral conditions
        DD = 0.97d0
        PP = 1d0
    elseif(zL < 0d0) then
        !> Unstable conditions
        DD = 0.28d0
        PP = 0.59d0
    elseif(zL > 0d0) then
        !> Stable conditions
        DD = 2.44d0
        PP = 1.33d0
    end if

    !> Footprint according to Hsieh et al., 2000
    do_offset = .true.
    do10 = .true.
    do30 = .true.
    do50 = .true.
    do70 = .true.
    do80 = .true.
    do i = 1, 10000
        !> Cross-wind integrated 1D function
        fact = DD * zu**PP * dabs(MO_length)**(1d0 - PP) / (vk**2 * (i * di))
        int_foot = dexp(-fact)
        if (do_offset .and. int_foot > 0.01d0) then
            Foot%offset = i * di
            do_offset = .false.
        end if
        if (do10 .and. int_foot > 0.1d0) then
            Foot%x10 = i * di
            do10 = .false.
        end if
        if (do30 .and. int_foot > 0.3d0) then
            Foot%x30 = i * di
            do30 = .false.
        end if
        if (do50 .and. int_foot > 0.5d0) then
            Foot%x50 = i * di
            do50 = .false.
        end if
        if (do70 .and. int_foot > 0.7d0) then
            Foot%x70 = i * di
            do70 = .false.
        end if
        if (do80 .and. int_foot > 0.8d0) then
            Foot%x80 = i * di
            do80 = .false.
        end if
        if (int_foot > 0.9d0) then
            Foot%x90 = i * di
            exit
        end if
    end do

    !> Peak distance
    Foot%peak = DD * zu**PP * dabs(MO_length)**(1d0 - PP) / (2d0 * vk**2)
end subroutine Hsieh00
