!***************************************************************************
! bpcf_massman_00.f90
! -------------------
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
! \brief       Calculates spectral correction factors for relevant fluxes
!              according to Massman (2000, 2001, AFM)
! \author      Gerardo Fratini
! \note        Corrected and refined by Stephen Chan, Lawrence Berkeley National Lab
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine bpcf_Massman00(measuring_height, displ_height, loc_var_present, LocInstr, wind_speed, t_air, zL, &
        avrg_length, detrending_time_constant, detrending_method, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: t_air
    real(kind = dbl), intent(in) :: zL
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(2), intent(in) :: detrending_method
    logical, intent(in) :: printout
    !> local variables
    integer :: var
    real(kind = dbl) :: t_sonic_hla_4tau
    real(kind = dbl) :: t_sonic_vla_4tau
    real(kind = dbl) :: t_sonic_la_4scalar
    real(kind = dbl) :: t_lat(GHGNumVar)
    real(kind = dbl) :: t_irga_la(GHGNumVar)
    real(kind = dbl) :: t_irga_va(GHGNumVar)
    real(kind = dbl) :: t_bw(GHGNumVar)
    real(kind = dbl) :: t_e(GHGNumVar)
    real(kind = dbl) :: t_tube(GHGNumVar)
    real(kind = dbl) :: t_ba
    real(kind = dbl) :: t_det
    real(kind = dbl) :: fx
    real(kind = dbl) :: aaa
    real(kind = dbl) :: bbb
    real(kind = dbl) :: pp(GHGNumVar)
    real(kind = dbl) :: alpha
    real(kind = dbl) :: unstable_corr_fact(GHGNumVar)
    real(kind = dbl) :: AirVisc
    real(kind = dbl) :: lambda(GHGNumVar)
    real(kind = dbl) :: Re
    real(kind = dbl) :: TubeVel
    real(kind = dbl), external :: LUT_delta


    if (printout) write(*,'(a)') '   Band-pass correction for all fluxes. Method: Massman (2000, 2001)..'

    !> Initialization
    t_e = 1d-10

    !> Equivalent time constants, Massman (2000, Table 1).
    t_sonic_hla_4tau    = LocInstr(sonic)%hpath_length / (2.8d0 * wind_speed) !< sonic line averaging (horizontal, for momentum)
    t_sonic_vla_4tau    = LocInstr(sonic)%vpath_length / (5.7d0 * wind_speed) !< sonic line averaging (vertical, for momentum)
    t_sonic_la_4scalar  = LocInstr(sonic)%vpath_length / (8.4d0 * wind_speed) !< sonic line averaging (for scalar flux)
    t_lat(co2:gas4)     = LocInstr(co2:gas4)%hsep  / (1.1d0 * wind_speed) !< lateral separation for all gases
    !t_lon(co2:gas4)     = LocInstr(co2:gas4)%lsep / (1.05d0 * wind_speed)  !< longitudinal separation (not used for now)
    t_irga_la(co2:gas4) = LocInstr(co2:gas4)%vpath_length / (4d0 * wind_speed) !< Line averaging scalar sensor
    t_irga_va(co2:gas4) = (0.2d0 + 0.4d0 * LocInstr(co2:gas4)%hpath_length / LocInstr(co2:gas4)%vpath_length) &
        *(LocInstr(co2:gas4)%vpath_length / wind_speed)   !< Volume averaging scalar sensor

    !> Bandwidth limitation             !<<<<<< add bandwidth limitations for analysers other than 7500x
    where (index(LocInstr(co2:gas4)%model, '7500') /= 0)
        t_bw(co2:gas4) = 0.016d0
    elsewhere
        t_bw(co2:gas4) = 1d-10  !< a very small value.
    endwhere

    !> Tube attenuation
    !delta_d = 2d0  !< to be expressed as a function of Re (should come from Massman 1991)
    !alpha_1 = 1d0  !< to be defined (should come from Massman 1991)
    AirVisc  = (-1.1555d-14 * t_air**3) &
        + (9.5728d-11 * t_air**2) + (3.7604d-8 * t_air) - 3.4484d-6

    lambda = error
    do var = co2, gas4
        if (loc_var_present(var)) then
            if (LocInstr(var)%path_type == 'closed') then
                TubeVel     = LocInstr(var)%tube_f / (p * (LocInstr(var)%tube_d / 2d0)**2)
                Re          = TubeVel * LocInstr(var)%tube_d / AirVisc
                lambda(var) = LUT_delta(Re, var)
                !lambda(var) = 0.5d0 * dabs(alpha_1) * delta_d**(-1) * Re**1.8  !< ****** Massman and Ibrom (2008, Eq. 8 and subsequent text), not used though
                if (lambda(var) /= error) then
                    t_tube(var) = dsqrt(lambda(var) * LocInstr(var)%tube_d / 2d0 * LocInstr(var)%tube_l) / (0.83d0 * TubeVel)
                else
                    call ExceptionHandler(51)
                    t_tube(var) = 1d-10 !< a very small value.
                end if
            else
                t_tube(var) = 1d-10 !< a very small value, open path case
            end if
        end if
    end do

    !> High pass filtering
    select case (trim(adjustl(detrending_method)))
        case ('ld')
            t_det = dble(avrg_length) * 60d0 / 5.3d0
        case ('ba')
            t_det = 1d10    !< a very large value.
        case ('rm', 'ew')
            t_det = detrending_time_constant
        case default
            t_det = 1d10    !< a very large value.
    end select

    !> Block averaging, which is done always
    t_ba = dble(avrg_length) * 60d0 / 2.8d0

    !> Combination of all time constants (Massman 2000, Eq. 9)
    !> For scalar fluxes
    where (loc_var_present(co2:gas4))
        t_e(w_co2:w_gas4) = dsqrt(t_sonic_la_4scalar**2 + t_lat(co2:gas4)**2 &
            + t_irga_la(co2:gas4)**2 + t_irga_va(co2:gas4)**2 + t_tube(co2:gas4)**2 + t_bw(co2:gas4)**2)
    elsewhere
        t_e(w_co2:w_gas4) = error
    endwhere

    !> For sensible heat
    t_e(w_ts) = dsqrt(t_sonic_hla_4tau**2 + t_sonic_vla_4tau**2)
    !> For momentum flux
    t_e(w_u)  = dsqrt(t_sonic_hla_4tau**2 + t_sonic_vla_4tau**2)

    !> Peak frequency (nx) and normalized frequency (fx), Massman (2000, Eq. 7 and following text)
    if(zL < 0d0) then
        fx = 0.085d0 * wind_speed / (measuring_height - displ_height)
        alpha = 0.925d0
    else
        fx = (2d0 - 1.915d0 /(1d0 + 0.5 * zL)) * wind_speed / (measuring_height - displ_height)
        alpha = 1d0
    end if

    !> Parameters a, b and p, Massman(2000, Eq. 11)
    aaa = (2d0 * p * fx * t_det)**alpha
    bbb = (2d0 * p * fx * t_ba )**alpha
    where (t_e(w_u:w_gas4) /= error)
        pp(w_u:w_gas4) = (2d0 * p * fx * t_e(w_u:w_gas4))**alpha
    elsewhere
        pp(w_u:w_gas4) = error
    endwhere

    !> Correction factors, Massman(2001, Table 1; refining Table 2 in Massman, 2000)
    do var = w_u, w_gas4
        if (pp(var) /= error .and. unstable_corr_fact(var) /= error) then
            BPCF%of(var) = (aaa*bbb/((aaa+1d0)*(bbb+1d0))) * (aaa*bbb/((aaa+pp(var))*(bbb+pp(var)))) &
                         * (1d0/(pp(var)+1d0)) * (1d0+(pp(var)+1d0)/(aaa+bbb))
        else
            BPCF%of(var) = error
        end if
    end do

    where (BPCF%of(w_u:w_gas4) /= 0d0 .and. BPCF%of(w_u:w_gas4) /= error)
        BPCF%of(w_u:w_gas4) = 1d0 / BPCF%of(w_u:w_gas4)
    elsewhere
        BPCF%of(w_u:w_gas4) = error
    end where

    if (printout) write(*,'(a)') '   Done.'
end subroutine BPCF_Massman00

!***************************************************************************
!
! \brief       Retrieve value of parameter Lambda needed in \n
!              Massman (2000, 2001, AFM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
double precision function LUT_delta(Re, var)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: Re
    integer, intent(in) :: var
    !> Local variables
    integer, parameter :: naux = 28
    integer :: i
    integer :: indx
    real(kind = dbl) :: Lambda_co2(naux)
    real(kind = dbl) :: Lambda_h2o(naux)
    real(kind = dbl) :: Lambda_ch4(naux)
    real(kind = dbl) :: Reynolds(naux)
    real(kind = dbl) :: Lambda(naux)

    data (Lambda_co2(mmm), mmm = 1, naux) &
        / 24.39, 12.21, 10.33, 8.13, 6.83, 5.3, 4.05, 2.97, 2.37, 1.99, &
         1.71, 1.51, 1.36, 1.23, 1.13, 1.04, 0.97, 0.91, 0.86, 0.81,   &
         0.77, 0.73, 0.7, 0.49, 0.39, 0.3, 0.25, 0.22/
    data (Lambda_h2o(mmm), mmm = 1, naux) &
        / 14.03, 8.22, 7.17, 5.89, 5.09, 4.1, 3.25, 2.47, 2.02, 1.72, &
         1.5, 1.34, 1.21, 1.11, 1.02, 0.95, 0.89, 0.83, 0.79, 0.75,  &
         0.71, 0.68, 0.65, 0.47, 0.38, 0.29, 0.24, 0.21/
    data (Lambda_ch4(mmm), mmm = 1, naux) &
        / 14.73, 8.52, 7.42, 6.07, 5.23, 4.2, 3.32, 2.52, 2.05, 1.74, 1.52, &
         1.35, 1.22, 1.12, 1.03, 0.96, 0.9, 0.84, 0.79, 0.75, 0.72, 0.68,  &
         0.65, 0.47, 0.38, 0.29, 0.24, 0.22/
    data (Reynolds(mmm), mmm = 1, naux) &
        / 2300, 2500, 2600, 2800, 3000, 3400, 4000, 5000, 6000, 7000, &
         8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, &
         17000, 18000, 19000, 20000, 30000, 40000, 60000, 80000, 100000 /


    select case(var)
        case (co2, gas4)
            Lambda = Lambda_co2
        case (h2o)
            Lambda = Lambda_h2o
        case (ch4)
            Lambda = Lambda_ch4
    end select

    !> Detect closest (smaller) Reynolds
    if (Re < 2300) then
        LUT_delta = error
        return
    end if

    indx = 1
    do i = 2, naux
        if (Reynolds(i) > Re) then
            indx = i - 1
            exit
        end if
    end do
    LUT_delta = (Lambda(indx + 1) - Lambda(indx)) &
        / (Reynolds(indx + 1) - Reynolds(indx)) * (Re - Reynolds(indx)) + Lambda(indx)
end function LUT_delta



