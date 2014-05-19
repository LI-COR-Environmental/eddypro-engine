!***************************************************************************
! fluxes1_rp.f90
! --------------
!Copyright (C) 2011, LI-COR Biosciences
!
!This file is part of EddyPro (TM).
!
!EddyPro (TM) is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!EddyPro (TM) is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!***************************************************************************

!***************************************************************************
! \file        src/fluxes1_rp.f90
! \brief       Calculates fluxes at Level 1. Mainly spectral corrections
! \version     3.0.0
! \date        2011-03-24
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Merge with corresponding FCC sub into a common one
!***************************************************************************
subroutine Fluxes1_rp()
    use m_rp_global_var
    implicit none
    real (kind = dbl) :: Cox

    write(*,'(a)', advance = 'no') '  Calculating fluxes Level 1..'

    Flux1 = fluxtype('', '', error, error, error, error, error, error, error, error, error, &
        error, error, error, error, error, error, error, error)

    !> First, apply oxygen correction to Krypton and Lyman-alpha hygrometers,
    !> according to van Dijk et al. (2003, JAOT, eq. 13b)
    select case (E2Col(h2o)%Instr%model(1:len_trim(E2Col(h2o)%Instr%model) - 2))
        case('open_path_krypton','closed_path_krypton', 'open_path_lyman','closed_path_lyman')
            if (E2Col(h2o)%Instr%kw /= 0d0 .and. LitePar%Ta > 0d0 &
                .and. LitePar%Bowen /= error .and. LitePar%lambda > 0) then
                Cox = 1d0 + 0.23d0 * E2Col(h2o)%Instr%ko / E2Col(h2o)%Instr%kw &
                    * LitePar%Bowen * LitePar%lambda / LitePar%Ta
                Stats%Cov(w, h2o) = Cox * Stats%Cov(w, h2o)
                Stats%Cov(h2o, h2o) = Cox**2 * Stats%Cov(h2o, h2o)
                !> Alternative formulation by T.W. Horst
                !> http://www.eol.ucar.edu/instrumentation/sounding/isfs/isff-support-center&
                !> &/how-tos/corrections-to-sensible-and-latent-heat-flux-measurements
                !Stats%Cov(w, h2o) = Stats%Cov(w, h2o) / (1 - 8d0 * 0.23d0 * E2Col(h2o)%Instr%ko &
                ! / E2Col(h2o)%Instr%kw * LitePar%bowen)
            endif
    end select

    !> Sensible heat flux, H in [W m-2]
    Flux1%H = Flux0%H

    !> Internal sensible heat flux, Hint in [W m-2]
    Flux1%Hi_co2 = Flux0%Hi_co2
    Flux1%Hi_h2o = Flux0%Hi_h2o
    Flux1%Hi_ch4 = Flux0%Hi_ch4
    Flux1%Hi_gas4 = Flux0%Hi_gas4

    !> Level 1 all gases
    !> co2
    if (E2Col(co2)%Instr%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%co2 = Flux0%co2
    else
        !> For all open-path gases, apply BPCF to L0: get L1
        if (BPCF%of(w_co2) /= error) then
            Flux1%co2 = Flux0%co2 * BPCF%of(w_co2)
        else
            Flux1%co2 = Flux0%co2
        end if
    end if
    if (Flux0%co2 == error) Flux1%co2 = error

    !> h2o
    if (E2Col(h2o)%Instr%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%h2o = Flux0%h2o
        Flux1%E   = Flux0%E
        Flux1%LE  = Flux0%LE
    else
        !> For all open-path gases, applied BPCF to LO: get L1
        if (BPCF%of(w_h2o) /= error) then
            Flux1%h2o = Flux0%h2o * BPCF%of(w_h2o)
            Flux1%E   = Flux0%E   * BPCF%of(w_h2o)
            Flux1%LE  = Flux0%LE  * BPCF%of(w_h2o)
        else
            Flux1%h2o = Flux0%h2o
            Flux1%E   = Flux0%E
            Flux1%LE  = Flux0%LE
        end if
    end if
    if (Flux0%h2o == error) Flux1%h2o = error
    if (Flux0%h2o == error) Flux1%E   = error
    if (Flux0%h2o == error) Flux1%LE  = error

    !> ch4
    if (E2Col(ch4)%Instr%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%ch4 = Flux0%ch4
    else
        !> For all open-path gases, applied BPCF to LO: get L1
        if (BPCF%of(w_ch4) /= error) then
            Flux1%ch4 = Flux0%ch4 * BPCF%of(w_ch4)
        else
            Flux1%ch4 = Flux0%ch4
        end if
    end if
    if (Flux0%ch4 == error) Flux1%ch4 = error

    !> gas4
    if (E2Col(gas4)%Instr%path_type == 'closed') then
        !> For all closed-path gases, Level 1 is same as Level 0
        Flux1%gas4 = Flux0%gas4
    else
        !> For all open-path gases, applied BPCF to L0: get L1
        if (BPCF%of(w_n2o) /= error) then
            Flux1%gas4 = Flux0%gas4 * BPCF%of(w_n2o)
        else
            Flux1%gas4 = Flux0%gas4
        end if
    end if
    if (Flux0%gas4 == error) Flux1%gas4 = error

    !> Level 1 evapotranspiration fluxes with H2O covariances at timelags of other scalars
    !> Do nothing, no spectral correction needed
    Flux1%E_co2 = Flux0%E_co2
    Flux1%E_ch4 = Flux0%E_ch4
    Flux1%E_gas4 = Flux0%E_gas4

    !> Momentum flux [kg m-1 s-2] and friction velocity [m s-1]
    if (BPCF%of(w_u) /= error) then
        Flux1%tau = Flux0%tau * BPCF%of(w_u)
        LitePar%us    = LitePar%us * dsqrt(BPCF%of(w_u))
    else
        Flux1%tau = Flux0%tau
    end if
    if (Flux0%tau == error) Flux1%tau = error

    write(*,'(a)')   ' done.'
end subroutine Fluxes1_rp
