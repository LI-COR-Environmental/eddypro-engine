!***************************************************************************
! bpcf_bandpass_spectral_corrections.f90
! --------------------------------------
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
! \brief       Hub to several spectral correction methods
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BandPassSpectralCorrections(measuring_height, displ_height, loc_var_present, wind_speed, t_air, zL, &
    ac_frequency, avrg_length, detrending_method, detrending_time_constant, nfull, printout, LocInstr, &
    LocFileList, nrow_full, lEx, LocSetup)
    use m_common_global_var
    implicit none
    !> In/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: t_air
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: nrow_full
    character(8), intent(in) :: detrending_method
    integer, intent(in) :: detrending_time_constant
    logical, intent(in) :: printout
    integer, intent(in) :: nfull
    !> Optional variables
    type(FileListType), optional, intent(in) :: LocFileList(nfull)
    type(ExType), optional, intent(in) :: lEx
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    integer :: month
    integer :: gas
    character(32) :: actual_hf_method


    !> Checks that parameters are passed correctly
    if (app == 'EddyPro-FCC') then
        if (.not. present(lEx) .or. .not. present(LocSetup) .or. .not. present(LocFileList)) &
        call ExceptionHandler(52)
    end if

    !> Spectral correction factors for anemometric fluxes (tau, H)
    BPCF%of(w_u) = 1d0
    BPCF%of(w_ts) = 1d0
    call BPCF_AnemometricFluxes(measuring_height, displ_height, loc_var_present, LocInstr, wind_speed, t_air, zL, &
        ac_frequency, avrg_length, detrending_time_constant, detrending_method, printout)

    !> Spectral correction factors for all gases
    BPCF%of(w_co2: w_gas4)  = 1d0

    !> Relevant only to FCC - Before entering correction method, check if the selected method
    !> can be implemented in the current situation
    !> (defined by RH for H2O and by the month for other gases). If not, sets the method to Moncrieff et al. 1997
    !> Note that even if only one condition fails, the method is set to Moncrieff for all gases
    actual_hf_method = trim(adjustl(EddyProProj%hf_meth))
    if (app == 'EddyPro-FCC') then
        select case (trim(adjustl(EddyProProj%hf_meth)))
            case('horst_97', 'ibrom_07', 'fratini_12')
                call char2int(lEx%date(6:7), month, 2)
                if(lEx%var_present(h2o) .and. (RegPar(dum, dum)%e1 == error &
                    .or. RegPar(dum, dum)%e2 == error .or. RegPar(dum, dum)%e3 == error)) then
                    actual_hf_method = 'moncrieff_97'
                    write(*,*)
                    call ExceptionHandler(69)
                end if
                if (actual_hf_method /= 'moncrieff_97') then
                    do gas = co2, gas4
                        if(gas /= h2o .and. lEx%var_present(gas)) then
                            if(RegPar(gas,  LocSetup%SA%class(gas, month))%fc == error) then
                                actual_hf_method = 'moncrieff_97'
                                call ExceptionHandler(69)
                                exit
                            end if
                        end if
                    end do
                end if
        end select
    end if
    EddyProProj%hf_meth = actual_hf_method

    select case(trim(adjustl(actual_hf_method)))
        case('none', 'not')
            if (EddyProProj%lf_meth == 'analytic') &
                call BPCF_OnlyLowFrequencyCorrection(measuring_height, displ_height, loc_var_present, wind_speed, zL, &
                    ac_frequency, avrg_length, detrending_time_constant, detrending_method)
        case('moncrieff_97')
            !> Correction after Moncrieff et al (1997, JH) fully analytical
            call BPCF_Moncrieff97(measuring_height, displ_height, loc_var_present, LocInstr, wind_speed, t_air, zL, &
                ac_frequency, avrg_length, detrending_time_constant, detrending_method, printout)
        case('massman_00')
            !> Correction after Massman (2000, 2001) fully analytical
            call BPCF_Massman00(measuring_height, displ_height, loc_var_present, LocInstr, wind_speed, t_air, zL, &
                avrg_length, detrending_time_constant, detrending_method, printout)
        case('horst_97')
            if (app == 'EddyPro-FCC') then
                !> Correction after Horst (1997, BLM), in-situ/analytical
                call BPCF_Horst97(measuring_height, displ_height, loc_var_present, wind_speed, zL, &
                    ac_frequency, avrg_length, detrending_time_constant, detrending_method, lEx, LocSetup)

                if (LocSetup%SA%horst_lens09 /= 'none') then
                    call CF_HorstLenschow09(lEx, LocSetup)
                    where (ADDCF%of(co2:gas4) < dabs(error))
                        BPCF%of(co2:gas4) = BPCF%of(co2:gas4) * ADDCF%of(co2:gas4)
                    end where
                end if
            end if
        case('ibrom_07')
            if (app == 'EddyPro-FCC') then
                !> Correction after Ibrom et al (2007, AFM), fully in-situ
                call BPCF_Ibrom07(measuring_height, displ_height, loc_var_present, wind_speed, zL, &
                    ac_frequency, avrg_length, detrending_time_constant, detrending_method, lEx, LocSetup)
                if (LocSetup%SA%horst_lens09 /= 'none') then
                    call CF_HorstLenschow09(lEx, LocSetup)
                    where (ADDCF%of(co2:gas4) < dabs(error))
                        BPCF%of(co2:gas4) = BPCF%of(co2:gas4) * ADDCF%of(co2:gas4)
                    end where
                end if
            end if
        case('fratini_12')
            if (app == 'EddyPro-FCC') then
                !> Correction after Fratini et al. 2012, AFM
                call BPCF_Fratini12(loc_var_present, LocInstr, wind_speed, t_air, ac_frequency, avrg_length, &
                    detrending_time_constant, detrending_method, nfull, nrow_full, LocFileList, lEx, LocSetup)

                if (LocSetup%SA%horst_lens09 /= 'none') then
                    call CF_HorstLenschow09(lEx, LocSetup)
                    where (ADDCF%of(co2:gas4) < dabs(error))
                        BPCF%of(co2:gas4) = BPCF%of(co2:gas4) * ADDCF%of(co2:gas4)
                    end where
                end if
            end if
    end select
end subroutine BandPassSpectralCorrections
