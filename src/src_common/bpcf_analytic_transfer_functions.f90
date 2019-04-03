!***************************************************************************
! bpcf_analytic_transfer_functions.f90
! ------------------------------------
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
! \brief       Calculates analytic form of the low-pass transfer function \n
!              depending on the instrumental setup and average wind-speed
!              after, Moncrieff et al. (1997), J. of Hydrology
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AnalyticLowPassTransferFunction(nf, N, var, LocInstr, loc_var_present, &
    wind_speed, t_air, BPTF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: var
    real(kind = dbl), intent(in) :: nf(N)
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: t_air
    type(BPTFType), intent(out) :: BPTF(N)
    !> local variables
!    integer :: Nba
    real(kind = dbl) :: fp_irga(N)
    real(kind = dbl) :: fp_sonic(N)
    real(kind = dbl) :: fs_ver(N)
    real(kind = dbl) :: fs_lat(N)
    real(kind = dbl) :: tube_velocity
    real(kind = dbl) :: tube_time
    real(kind = dbl) :: Re
    real(kind = dbl) :: sonic_freq
    real(kind = dbl) :: irga_freq
    real(kind = dbl) :: air_viscosity


    if (loc_var_present(var)) then
        select case(var)
            case (u, v, w, ts)
                sonic_freq = 1d0 / LocInstr(sonic)%tau
                fp_sonic(1:N) = &
                    nf(1:N) * dabs(LocInstr(sonic)%vpath_length / wind_speed)
                !> sonic dynamic frequency response
                BPTF(1:N)%LP(var)%dsonic = &
                    1.d0 / dsqrt(1.d0 + (2.d0 * p * nf(1:N) / sonic_freq)**2 )
                !> sonic path-averaging
                BPTF(1:N)%LP(var)%wsonic = &
                    (2.d0 / (p * fp_sonic(1:N))) * (1.d0 + dexp(-2.d0 * p * fp_sonic(1:N)) / 2.d0 &
                    - 3.d0 * (1.d0 - dexp(-2.d0 * p * fp_sonic(1:N))) / (4.d0 * p * fp_sonic(1:N)))

                if (LocInstr(var)%category == 'fast_t_sensor') then
                    !> For a fast temperature sensor (typically a thermocouple)
                    !> spectral losses can be safely neglected
                    BPTF(1:N)%LP(var)%dsonic = 1d0
                    BPTF(1:N)%LP(var)%wsonic = 1d0
                end if

            case (co2, h2o, ch4, gas4)
                irga_freq = 1d0 / LocInstr(var)%tau
                fp_irga(1:N) = nf(1:N) * dabs(LocInstr(var)%vpath_length / wind_speed)
                fs_ver(1:N)  = nf(1:N) * dabs(LocInstr(var)%vsep / wind_speed)
                fs_lat(1:N)  = nf(1:N) * dabs(LocInstr(var)%hsep / wind_speed)
                if (LocInstr(var)%path_type == 'closed') then
                    !> average speed [m/s], time [s] and Reynolds number in the tube
                    !> (air viscosity as a function of temperature)
                    air_viscosity  = (-1.1555d-14 * t_air**3) + &
                        (9.5728d-11 * t_air**2) + (3.7604d-8 * t_air) - 3.4484d-6
                    tube_velocity  = LocInstr(var)%tube_f / (p * (LocInstr(var)%tube_d / 2d0)**2)
                    ! if (var == h2o) &
                        !tube_velocity  = tube_velocity * Essentials%used_timelag(co2) / Essentials%used_timelag(h2o)
                    tube_time = LocInstr(var)%tube_l / tube_velocity
                    Re       = tube_velocity * LocInstr(var)%tube_d / air_viscosity
                    !> attenuantions in the intake tube
                    if(Re < 2300.d0) then
                        !> for laminar regimes
                        BPTF(1:N)%LP(var)%t = &
                            dexp(- tube_time * (p * LocInstr(var)%tube_d / 2d0 &
                            * nf(1:N))**2 / (6.d0 * Dc(var)))
                    else
                        !> for turbulent regimes (Lenshow and Raupach (1991), JGR, eq. 7)
                        BPTF(1:N)%LP(var)%t = &
                            dexp((-80.d0 * LocInstr(var)%tube_d / 2.d0 * (nf(1:N)**2) &
                            * (Re**(-0.125d0)) * LocInstr(var)%tube_l) / (tube_velocity**2))
                    end if
                else
                    BPTF(1:N)%LP(var)%t = 1d0
                end if
                !> transfer functions common to open and closed path irgas
                !> irga dynamic frequency response
                BPTF(1:N)%LP(var)%dirga  = 1.d0 / dsqrt(1.d0 + (2.d0 * p * nf(1:N) /  irga_freq)**2 )
                !> irga path-averaging
                BPTF(1:N)%LP(var)%wirga = &
                    dsqrt((3d0 + dexp(-2d0 * p * fp_irga(1:N)) - 4d0 / (2d0 * p * fp_irga(1:N)) * &
                    (1d0 - dexp(-2d0 * p * fp_irga(1:N)))) / (2d0 * p * fp_irga(1:N)))
                !> In case something went wrong with the sqrt above,
                !> set NaN to error. Problem seen in support case
                !> with exceedingly low (implausible) instrument path length.
                where (isnan(BPTF(1:N)%LP(var)%wirga))
                    BPTF(1:N)%LP(var)%wirga = error
                end where

                !> irga separations (vertical and horizontal)
                BPTF(1:N)%LP(var)%sver = dexp(-9.9d0 * fs_ver(1:N)**1.5d0)
                BPTF(1:N)%LP(var)%shor = dexp(-9.9d0 * fs_lat(1:N)**1.5d0)
        end select
    end if

!!> irga/sonic response mismatch
!BPTF(i)%LP(var)%m = (1.d0 + (irga_freq**(-1)) * (sonic_freq**(-1)) &
!                  * (2.d0 * p * nf(i))**2) &
!                  / dsqrt( (1.d0 + (2.d0 * p * nf(i) / irga_freq)**2) &
!                         * (1.d0 + (2.d0 * p * nf(i) / sonic_freq)**2))

end subroutine AnalyticLowPassTransferFunction

!***************************************************************************
! \brief       Calculates theoretical high-pass transfer function.
!              See Moncrieff et al. (2004, HoM),
!              Kristensen (1998), Aubinet et al. (2000, EuroFlux meth.) and
!              Geissler & Ibrom (2008)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AnalyticHighPassTransferFunction(nf, N, var, ac_frequency, avrg_length, &
        detrending_method, detrending_time_constant, BPTF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: var
    real(kind = dbl), intent(in) :: nf(N)
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(2), intent(in) :: detrending_method
    type(BPTFType), intent(out) :: BPTF(N)
    !> local variables
    real(kind = dbl) :: dt
    real(kind = dbl) :: total_time
    real(kind = dbl) :: eta
    real(kind = dbl) :: tau
    real(kind = dbl) :: phase(N)
    real(kind = dbl) :: phase_ld(N)


    !> Low-frequency transfer function depends on de-trending method and related time constant
    dt = 1d0 / ac_frequency              !< in seconds
    total_time = dble(avrg_length) * 60d0  !< in seconds
    if (total_time /= 0) then
        select case(trim(adjustl(detrending_method)))
            case('ba')
                phase(1:N) = p * nf(1:N) * total_time
                BPTF(1:N)%HP(var) = 1d0 - dsin(phase(1:N))**2 / phase(1:N)**2
            case('ld')
                !> Note that the time constant for TF part that refers to the linear detrending
                !> in general is not the averaging period (total_time),
                !> but the actual user-defined time constant
                !> For the definition see Geissler & Ibrom 2008, Kristensen 1998 and Aubinet 2000.
                if (detrending_time_constant /= 0 .and. detrending_time_constant /= nint(error)) then
                    phase(1:N)    = p * nf(1:N) * total_time
                    phase_ld(1:N) = p * nf(1:N) * dble(detrending_time_constant)
                    BPTF(1:N)%HP(var) = 1d0 &
                        - dsin(phase(1:N))**2 / phase(1:N)**2 &
                        - 3d0 * (dsin(phase_ld(1:N)) / phase_ld(1:N) - dcos(phase_ld(1:N)))**2 &
                        / phase_ld(1:N)**2
                else
                    BPTF(1:N)%HP(var) = 1d0
                end if
            case('rm', 'ew')
                !> For the definition of eta and tau (a and t_d) see Geissler & Ibrom 2008
                if (detrending_time_constant /= 0 .and. detrending_time_constant /= nint(error)) then
                    eta = dexp(-dt / dble(detrending_time_constant))
                    tau = eta / (ac_frequency * (1d0 - eta))
                    phase(1:N) = 2d0 * p * nf(1:N) * tau
                    BPTF(1:N)%HP(var) = phase(1:N)**2 / (1d0 + phase(1:N)**2 / eta)
                else
                    BPTF(1:N)%HP(var) = 1d0
                end if
        end select
    end if
end subroutine AnalyticHighPassTransferFunction

!***************************************************************************
! \brief       Calculates low-pass transfer functions for block-averaging
!              and DAC zero-order hold filtering in the LI-7550 data logger
!              occurring prior to software version 7.7.0
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LI7550_AnalogSignalsTransferFunctions(nf, N, var, ac_frequency, &
        loc_var_present, BPTF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: var
    real(kind = dbl), intent(in) :: nf(N)
    real(kind = dbl), intent(in) :: ac_frequency
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(BPTFType), intent(out) :: BPTF(N)
    !> local variables
    integer :: Nba
    real(kind = dbl) :: Tba
    real(kind = dbl), parameter :: li7550_sonic_sampling_frequency = 20d0
    include 'interfaces.inc'


    !> TF related to analog signals filtering in the LI-7550
    Nba = nint(20 / ac_frequency)
    Tba = dfloat(Nba) / li7550_sonic_sampling_frequency

    if (loc_var_present(var)) then
        select case(var)
            case(u, v, w, ts)
                !> Block averaging
                BPTF(1:N)%LP(var)%ba_sonic = dsqrt(dabs(sinc(nf(1:N)*Tba, N)))
                !> ZOH
                BPTF(1:N)%LP(var)%zoh_sonic = &
                    dsqrt(dabs(sinc(nf(1:N)/EddyProProj%sonic_output_rate/2d0, N)))
            case(co2, h2o)
                BPTF(1:N)%LP(var)%ba_irga = dsqrt(dabs(sinc(nf(1:N)*Tba, N)))
        end select
    end if
end subroutine LI7550_AnalogSignalsTransferFunctions
