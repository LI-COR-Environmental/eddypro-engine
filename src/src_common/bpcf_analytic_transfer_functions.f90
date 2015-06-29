!***************************************************************************
! bpcf_analytic_transfer_functions.f90
! ------------------------------------
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
subroutine AnalyticLowPassTransferFunction(nf, N, var, LocInstr, loc_var_present, wind_speed, t_air, BPTF)
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
                fp_sonic(1:N) = nf(1:N) * dabs(LocInstr(sonic)%vpath_length / wind_speed)
                !> sonic dynamic frequency response
                BPTF(1:N)%LP(var)%dsonic = 1.d0 / dsqrt(1.d0 + (2.d0 * p * nf(1:N) / sonic_freq)**2 )
                !> sonic path-averaging
                BPTF(1:N)%LP(var)%wsonic = (2.d0 / (p * fp_sonic(1:N))) * (1.d0 + dexp(-2.d0 * p * fp_sonic(1:N)) / 2.d0 &
                                        - 3.d0 * (1.d0 - dexp(-2.d0 * p * fp_sonic(1:N))) / (4.d0 * p * fp_sonic(1:N)))

                if (LocInstr(var)%category == 'fast_t_sensor') then
                    !> For a fast temperature sensor (typically a thermocouple) spectral losses can be safely neglected
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
                        !tube_velocity  = tube_velocity * Essentials%timelag(co2) / Essentials%timelag(h2o)
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
                !> In case something went wrong with the sqrt above, set NaN to error. Problem seen in support case
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
