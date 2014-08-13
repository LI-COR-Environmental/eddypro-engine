!***************************************************************************
! fluxes1.f90
! -----------
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
! \brief       Calculates fluxes at Level 1. Mainly spectral corrections
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Fluxes1(lEx)
    use m_fx_global_var
    implicit none
    !> In/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    real (kind = dbl) :: Cox

    Flux1 = fluxtype('', '', error, error, error, error, error, error, error, error, error, &
        error, error, error, error, error, error, error, error)

    !> First, apply oxygen correction to Krypton and Lyman-alpha hygrometers,
    !> according to van Dijk et al. (2003, JAOT, eq. 13b)
    select case (lEx%instr(ih2o)%model(1:len_trim(lEx%instr(ih2o)%model) - 2))
        case('open_path_krypton','closed_path_krypton', 'open_path_lyman','closed_path_lyman')
            if (lEx%instr(ih2o)%kw /= 0d0 .and. lEx%Ta > 0d0 &
                .and. lEx%Bowen /= error .and. lEx%lambda > 0d0) then
                Cox = 1d0 + 0.23d0 * lEx%instr(ih2o)%ko / lEx%instr(ih2o)%kw &
                    * lEx%Bowen * lEx%lambda / lEx%Ta
                lEx%cov_w(h2o) = Cox * lEx%cov_w(h2o)
                lEx%var(h2o) = Cox**2 * lEx%var(h2o)
                !> Alternative formulation by T.W. Horst
                !> http://www.eol.ucar.edu/instrumentation/sounding/isfs/isff-support-center&
                !> &/how-tos/corrections-to-sensible-and-latent-heat-flux-measurements
                !lEx%cov_w(h2o) = lEx%cov_w(h2o) / (1 - 8d0 * 0.23d0 * lEx%instr(ih2o)%ko &
                ! / lEx%instr(ih2o)%kw * lEx%bowen)
            endif
    end select

    !> Sensible heat flux, H in [W m-2]
    Flux1%H = lEx%Flux0%H

    !> Internal sensible heat flux, Hint in [W m-2]
    Flux1%Hi_co2 = lEx%Flux0%Hi_co2
    Flux1%Hi_h2o = lEx%Flux0%Hi_h2o
    Flux1%Hi_ch4 = lEx%Flux0%Hi_ch4
    Flux1%Hi_gas4 = lEx%Flux0%Hi_gas4

    !> Level 1 all gases
    !> co2
    if (lEx%instr(ico2)%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%co2 = lEx%Flux0%co2
    else
        !> For all open-path gases, apply BPCF to L0: get L1
        if (BPCF%of(w_co2) /= error) then
            Flux1%co2 = lEx%Flux0%co2 * BPCF%of(w_co2)
        else
            Flux1%co2 = lEx%Flux0%co2
        end if
    end if
    if (lEx%Flux0%co2 == error) Flux1%co2 = error

    !> h2o
    lEx%Flux0%E = lEx%Flux0%LE / lEx%lambda
    if (lEx%instr(ih2o)%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%LE  = lEx%Flux0%LE
        Flux1%h2o = lEx%Flux0%h2o
        Flux1%E   = lEx%Flux0%E
    else
        !> For all open-path gases, applied BPCF to LO: get L1
        if (BPCF%of(w_h2o) /= error) then
            Flux1%h2o = lEx%Flux0%h2o * BPCF%of(w_h2o)
            Flux1%E   = lEx%Flux0%E   * BPCF%of(w_h2o)
            Flux1%LE  = lEx%Flux0%LE  * BPCF%of(w_h2o)
        else
            Flux1%h2o = lEx%Flux0%h2o
            Flux1%E   = lEx%Flux0%E
            Flux1%LE  = lEx%Flux0%LE
        end if
    end if
    if (lEx%Flux0%h2o == error) Flux1%h2o = error
    if (lEx%Flux0%h2o == error) Flux1%E   = error
    if (lEx%Flux0%h2o == error) Flux1%LE  = error

    !> ch4
    if (lEx%instr(ich4)%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%ch4 = lEx%Flux0%ch4
    else
        !> For all open-path gases, applied BPCF to LO: get L1
        if (BPCF%of(w_ch4) /= error) then
            Flux1%ch4 = lEx%Flux0%ch4 * BPCF%of(w_ch4)
        else
            Flux1%ch4 = lEx%Flux0%ch4
        end if
    end if
    if (lEx%Flux0%ch4 == error) Flux1%ch4 = error

    !> n2o
    if (lEx%instr(igas4)%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%gas4 = lEx%Flux0%gas4
    else
        !> For all open-path gases, applied BPCF to L0: get L1
        if (BPCF%of(w_gas4) /= error) then
            Flux1%gas4 = lEx%Flux0%gas4 * BPCF%of(w_gas4)
        else
            Flux1%gas4 = lEx%Flux0%gas4
        end if
    end if
    if (lEx%Flux0%gas4 == error) Flux1%gas4 = error

    !> Level 1 evapotranspiration fluxes with H2O covariances at timelags of other scalars
    !> Do nothing, no spectral correction needed
    Flux1%E_co2 = lEx%Flux0%E_co2
    Flux1%E_ch4 = lEx%Flux0%E_ch4
    Flux1%E_gas4 = lEx%Flux0%E_gas4

    !> Momentum flux [kg m-1 s-2] and friction velocity [m s-1]
    if (BPCF%of(w_u) /= error) then
        Flux1%tau = lEx%Flux0%tau * BPCF%of(w_u)
        Ambient%us    = Ambient%us * dsqrt(BPCF%of(w_u))
    else
        Flux1%tau = lEx%Flux0%tau
    end if
    if (lEx%Flux0%tau == error) Flux1%tau = error
end subroutine Fluxes1
