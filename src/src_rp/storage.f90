!***************************************************************************
! storage.f90
! -----------
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
! \brief       Estimates storage fluxes based on single-point measurements or
!              on profiles, as available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Storage(PrevStats, prevAmbient)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(StatsType), intent(in) :: PrevStats
    type(AmbientStateType), intent(in) :: prevAmbient
    !> local variables
    integer, parameter :: stH = 1
    integer, parameter :: stCO2 = 2
    integer, parameter :: stH2O = 3
    integer, parameter :: stCH4 = 4
    integer, parameter :: stGAS4 = 5
    integer :: gas
    integer :: i
    integer :: j
    integer :: var
    real (kind = dbl) :: seconds
    real (kind = dbl) :: dcdt(5, MaxProfNodes)
    real (kind = dbl) :: dStor(5, MaxProfNodes - 1)
    real (kind = dbl) :: dz(5, MaxProfNodes - 1)
    character(10) tmp_date
    character(5) tmp_time


    write(*, '(a)', advance = 'no') '  Calculating storage terms..'

    !> Check that time periods are consecutive. If not, set storage to error and exit
    call AddDateStep(PrevStats%date, PrevStats%time, tmp_date, tmp_time, DateStep)

    if (tmp_date /= Stats%date .or. tmp_time /= Stats%time) then
        Stor%H  = error
        Stor%LE = error
        Stor%of(co2:gas4)  = error
        return
    end if

!    seconds = RPsetup%avrg_len * 6d1
!    !> define concentration differences in time, at each height
!    where (PrevSlowVar%prof_t(:) /= error .and. BiometVar%prof_t(:) /= error)
!        dcdt(stH, :) = (BiometVar%prof_t(:) - PrevSlowVar%prof_t(:)) / seconds
!    elsewhere
!        dcdt(stH, :) = error
!    end where
!    where (PrevSlowVar%prof_co2(:) /= error .and. BiometVar%prof_co2(:) /= error)
!        dcdt(stCO2, :) = (BiometVar%prof_co2(:) - PrevSlowVar%prof_co2(:)) / seconds
!    elsewhere
!        dcdt(stCO2, :) = error
!    end where
!    where (PrevSlowVar%prof_h2o(:) /= error .and. BiometVar%prof_h2o(:) /= error)
!        dcdt(stH2O, :) = (BiometVar%prof_h2o(:) - PrevSlowVar%prof_h2o(:)) / seconds
!    elsewhere
!        dcdt(stH2O, :) = error
!    end where
!    where (PrevSlowVar%prof_ch4(:) /= error .and. BiometVar%prof_ch4(:) /= error)
!        dcdt(stCH4, :) = (BiometVar%prof_ch4(:) - PrevSlowVar%prof_ch4(:)) / seconds
!    elsewhere
!        dcdt(stCH4, :) = error
!    end where
!    where (PrevSlowVar%prof_gas4(:) /= error .and. BiometVar%prof_gas4(:) /= error)
!        dcdt(stGAS4, :) = (BiometVar%prof_gas4(:) - PrevSlowVar%prof_gas4(:)) / seconds
!    elsewhere
!        dcdt(stGAS4, :) = error
!    end where

    !dcdt(stCO2, 2) = error
!
!    !> Initialize dz at their natural values
!    dz = bSetup%dz
!
!    !> Storage fluxs as integration of profile differences (using the trapezi formula)
!    Stor%H = 0d0
!    Stor%LE = 0d0
!    Stor%of(:) = 0d0
!    dStor = 0d0
!    do var = stH, stGAS4
!        ol: do i = 1, MaxProfNodes - 1
!            if (dcdt(var, i) /= error) then
!                il: do j = i + 1, MaxProfNodes
!                    if (dcdt(var, j) /= error) then
!                        dStor(var, i) = (dcdt(var, j) + dcdt(var, i)) * dz(var, i) * 0.5d0
!                        exit il
!                    else
!                        dz(var, i) = dz(var, i) + bSetup%dz(var, j)
!                    end if
!                end do il
!                Stor%H = Stor%H + dStor(stH, i)
!                Stor%of(co2) = Stor%of(co2) + dStor(stCO2, i)
!                Stor%of(h2o) = Stor%of(h2o) + dStor(stH2O, i)
!                Stor%of(ch4) = Stor%of(ch4) + dStor(stCH4, i)
!                Stor%of(gas4) = Stor%of(gas4) + dStor(stGAS4, i)
!            end if
!        end do ol
!    end do
!
!    !> Units adjustments
!    if (Stor%H /= 0) Stor%H = Ambient%RhoCp * Stor%H
!    !> Gas from mixing ratio to molar density (should use Vd, not Va. But Vd not necessarily available at this stage)
!    do gas = co2, gas4
!        select case(gas)
!            case (co2, ch4, gas4)
!                if (Stor%of(gas) /= 0 .and. Ambient%Va /= 0) &
!                    Stor%of(gas) = Stor%of(gas) * 1d3 / Ambient%Va
!            case(h2o)
!                if (Stor%of(gas) /= 0 .and. Ambient%Va /= 0) &
!                    Stor%of(gas) = Stor%of(gas) / Ambient%Va
!        end select
!    end do
!
!    if (Stor%of(h2o) /= 0) then
!        Stor%LE = Stor%of(h2o) * MW(h2o) * Ambient%lambda * 1d-3
!    else
!        Stor%LE = 0d0
!    end if
!
    !> If Stor = 0, it means no profile was available or selected, so calculate it with 1-point formula
    !> Storage for sensible heat
    if (Stor%H == 0) then
        if(Ambient%RhoCp /= error .and. Ambient%Ta /= error .and. prevAmbient%Ta /= error) then
            Stor%H = Ambient%RhoCp * (Ambient%Ta - prevAmbient%Ta) * E2Col(u)%Instr%height / seconds
        else
            Stor%H = error
        end if
    end if

    !> for all gases
    do gas = co2, gas4
        if (Stor%of(gas) == 0) then
            select case(gas)
                case (co2, ch4, gas4)
                    if (Stats%chi(gas) /= error .and. PrevStats%chi(gas) /= error) then
                        Stor%of(gas) = (Stats%chi(gas) - PrevStats%chi(gas)) / Ambient%Va * 1d-3 &
                                     * E2Col(gas)%Instr%height * 1d3 / seconds
                    else
                        Stor%of(gas) = error
                    end if
                case(h2o)
                    if (Stats%chi(gas) /= error .and. PrevStats%chi(gas) /= error) then
                        Stor%of(gas) = (Stats%chi(gas) - PrevStats%chi(gas)) / Ambient%Va &
                                     * E2Col(gas)%Instr%height / seconds
                        Stor%LE = Stor%of(h2o) * MW(h2o) * Ambient%lambda * 1d-3
                    else
                        Stor%of(gas) = error
                        Stor%LE      = error
                    end if
            end select
        end if
    end do

    write(*, '(a)') ' Done.'
end subroutine Storage
