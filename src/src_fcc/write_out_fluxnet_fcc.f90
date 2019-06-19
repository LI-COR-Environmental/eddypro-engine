!***************************************************************************
! write_out_fluxnet_fcc.f90
! -------------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFluxnetFcc(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    character(16000) :: dataline

    !> local variables
    integer :: var
    integer :: i
    integer :: gas
    integer :: igas
    character(9) :: vm97flags(GHGNumVar)
    include '../src_common/interfaces_1.inc'


    call clearstr(dataline)
    !> Timestamp
    !> Start/end imestamps
    call AddDatum(dataline, trim(adjustl(lEx%start_timestamp)), separator)
    call AddDatum(dataline, trim(adjustl(lEx%end_timestamp)), separator)
    call AddFloatDatumToDataline(lEx%DOY_start, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%DOY_end, dataline, EddyProProj%err_label)

    !> Filename
    call AddCharDatumToDataline(lEx%fname, dataline, EddyProProj%err_label)

    !> Potential radiation and daytime
    call AddFloatDatumToDataline(lEx%RP, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(lEx%nighttime_int, dataline, EddyProProj%err_label)

    !> Number of records
    call AddIntDatumToDataline(lEx%nr_theor, dataline, EddyProProj%err_label)        
    call AddIntDatumToDataline(lEx%nr_files, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(lEx%nr_after_custom_flags, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(lEx%nr_after_wdf, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(lEx%nr(u), dataline, EddyProProj%err_label)
    do var = ts, gas4
        call AddIntDatumToDataline(lEx%nr(var), dataline, EddyProProj%err_label)
    end do
    call AddIntDatumToDataline(lEx%nr_w(u), dataline, EddyProProj%err_label)
    do var = ts, gas4
        call AddIntDatumToDataline(lEx%nr_w(var), dataline, EddyProProj%err_label)
    end do

    !> Final fluxes
    call AddFloatDatumToDataline(Flux3%Tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%ET, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux3%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux3%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)

    !> Random uncertainties
    call AddFloatDatumToDataline(lEx%rand_uncer(u), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(ts), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer_LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer_ET, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(co2), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(h2o), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(ch4), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(lEx%rand_uncer(gas4), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)

    !> Storage fluxes
    call AddFloatDatumToDataline(lEx%Stor%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Stor%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Stor%ET, dataline, EddyProProj%err_label)
    do gas = co2, h2o
        call AddFloatDatumToDataline(lEx%Stor%of(gas), dataline, EddyProProj%err_label)
        end do
    do gas = ch4, gas4
        call AddFloatDatumToDataline(lEx%Stor%of(gas), dataline, EddyProProj%err_label)
    end do

    !> Advection fluxes
    do gas = co2, gas4
        if (lEx%rot_w /= error .and. lEx%d(gas) >= 0d0) then
            if (lEx%rot_w /= error .and. lEx%d(gas) /= error) then
                if (gas == co2) then
                    call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), &
                        dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
                else if (gas == h2o) then
                    call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), dataline, EddyProProj%err_label)
                else if (gas == ch4 .or. gas == gas4) then
                    call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), &
                        dataline, EddyProProj%err_label, gain=1d6, offset=0d0)
                end if
            else
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

    !> Turbulence and micromet
    !> Unrotated and rotated wind components
    call AddFloatDatumToDataline(lEx%unrot_u, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%unrot_v, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%unrot_w, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_u, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_v, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_w, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%WS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%MWS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%WD, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%WD_SIGMA, dataline, EddyProProj%err_label)

    !> Turbulence
    call AddFloatDatumToDataline(Flux3%ustar, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%TKE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%L, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%zL, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%bowen, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Tstar, dataline, EddyProProj%err_label)

    !> Thermodynamics
    !> Temperature, pressure, RH, VPD, e, es, etc.
    call AddFloatDatumToDataline(lEx%Ts, dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(lEx%Ta, dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(lEx%Pa, dataline, EddyProProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(lEx%RH, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Va, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%RHO%a, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%RhoCp, dataline, EddyProProj%err_label)
    !> Water
    call AddFloatDatumToDataline(lEx%RHO%w, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%e, dataline, EddyProProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%es, dataline, EddyProProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%Q, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%VPD, dataline, EddyProProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%Tdew, dataline, EddyProProj%err_label, gain=1d0, offset=-273.15d0)
    !> Dry air
    call AddFloatDatumToDataline(lEx%Pd, dataline, EddyProProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(lEx%RHO%d, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Vd, dataline, EddyProProj%err_label)
    !> Specific heat of evaporation
    call AddFloatDatumToDataline(lEx%lambda, dataline, EddyProProj%err_label)
    !> Wet to dry air density ratio
    call AddFloatDatumToDataline(lEx%sigma, dataline, EddyProProj%err_label)

    !> Gas concentrations/densities
    do gas = co2, gas4
        call AddIntDatumToDataline(lEx%measure_type_int(gas), dataline, EddyProProj%err_label)
        call AddFloatDatumToDataline(lEx%d(gas), dataline, EddyProProj%err_label)
        if (gas == ch4 .or. gas == gas4) then
            call AddFloatDatumToDataline(lEx%r(gas), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
            call AddFloatDatumToDataline(lEx%chi(gas), dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
        else
            call AddFloatDatumToDataline(lEx%r(gas), dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(lEx%chi(gas), dataline, EddyProProj%err_label)
        end if
    end do

    !> Time lags
    do gas = co2, gas4
        call AddFloatDatumToDataline(lEx%act_tlag(gas), dataline, EddyProProj%err_label)
        call AddFloatDatumToDataline(lEx%used_tlag(gas), dataline, EddyProProj%err_label)
        call AddFloatDatumToDataline(lEx%nom_tlag(gas), dataline, EddyProProj%err_label)
        call AddFloatDatumToDataline(lEx%min_tlag(gas), dataline, EddyProProj%err_label)
        call AddFloatDatumToDataline(lEx%max_tlag(gas), dataline, EddyProProj%err_label)
    end do

    !> Stats
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%median(var), dataline, &
                EddyProProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%median(var), dataline, EddyProProj%err_label)
        end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%Q1(var), dataline, &
                EddyProProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%Q1(var), dataline, EddyProProj%err_label)
        end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%Q3(var), dataline, &
                EddyProProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%Q3(var), dataline, EddyProProj%err_label)
        end if
    end do
    do var = u, gas4
        call AddFloatDatumToDataline(sqrt(lEx%stats%Cov(var, var)), dataline, EddyProProj%err_label)
    end do
    do var = u, gas4
        call AddFloatDatumToDataline(lEx%stats%Skw(var), dataline, EddyProProj%err_label)
    end do
    do var = u, gas4
        call AddFloatDatumToDataline(lEx%stats%Kur(var), dataline, EddyProProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%stats%Cov(w, u), dataline, EddyProProj%err_label)
    do var = ts, gas4
        call AddFloatDatumToDataline(lEx%stats%Cov(w, var), dataline, EddyProProj%err_label)
    end do
    do var = h2o, gas4
        call AddFloatDatumToDataline(lEx%stats%Cov(co2, var), dataline, EddyProProj%err_label)
    end do
    do var = ch4, gas4
        call AddFloatDatumToDataline(lEx%stats%Cov(h2o, var), dataline, EddyProProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%stats%Cov(ch4, gas4), dataline, EddyProProj%err_label)

    !> Footprint
    call AddFloatDatumToDataline(Foot%peak, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%offset, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x10, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x30, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x50, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x70, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x80, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Foot%x90, dataline, EddyProProj%err_label)

    !> Fluxes Level 0 (uncorrected)
    call AddFloatDatumToDataline(lEx%Flux0%ustar, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%L, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%zL, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%ET, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(lEx%Flux0%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    !> Fluxes Level 1 
    call AddFloatDatumToDataline(Flux1%Tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%ET, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux1%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux1%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    !> Fluxes Level 2
    call AddFloatDatumToDataline(Flux2%Tau, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%H, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%LE, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%ET, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(Flux2%ch4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux2%gas4, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)

    !> Cell values
    call AddFloatDatumToDataline(lEx%Tcell, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Pcell, dataline, EddyProProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(lEx%Vcell(gas), dataline, EddyProProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%Flux0%E_co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%E_ch4, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%E_gas4, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_co2, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_h2o, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_ch4, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_gas4, dataline, EddyProProj%err_label)

    !> Burba terms
    call AddFloatDatumToDataline(lEx%Burba%h_bot, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Burba%h_top, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Burba%h_spar, dataline, EddyProProj%err_label)

    !> LI-7700 multipliers
    call AddFloatDatumToDataline(lEx%Mul7700%A, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Mul7700%B, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%Mul7700%C, dataline, EddyProProj%err_label)

    !> Spectral correction factors
    call AddFloatDatumToDataline(BPCF%of(w_u), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_ts), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_h2o), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_h2o), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_co2), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_h2o), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_ch4), dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_gas4), dataline, EddyProProj%err_label)

    !> Degraded covariances
    call AddFloatDatumToDataline(lEx%degT%cov, dataline, EddyProProj%err_label)
    do i = 1, 9
        call AddFloatDatumToDataline(lEx%degT%dcov(i), dataline, EddyProProj%err_label)
    end do
    do var = u, gas4
        call AddIntDatumToDataline(lEx%spikes(var), dataline, EddyProProj%err_label)
    end do

    !> Write first string from Chunks
    !> M_CUSTOM_FLAGS thru VM97_NSW_RNS
    call AddDatum(dataline, trim(fluxnetChunks%s(1)), separator)

    !> VM97 flags, here organized per variable instead of per test
    if (lEx%vm_flags(1) == '-9999') then
        do var = u, gas4
            call AddCharDatumToDataline(EddyProProj%err_label, dataline, EddyProProj%err_label)
        end do
    else
        do var = u, gas4
            vm97flags(var)(1 : 1) = '8' 
            vm97flags(var)(2 : 2) = lEx%vm_flags(1)(var + 1 : var + 1)
            vm97flags(var)(3 : 3) = lEx%vm_flags(2)(var + 1 : var + 1)
            vm97flags(var)(4 : 4) = lEx%vm_flags(3)(var + 1 : var + 1)
            vm97flags(var)(5 : 5) = lEx%vm_flags(4)(var + 1 : var + 1)
            vm97flags(var)(6 : 6) = lEx%vm_flags(5)(var + 1 : var + 1)
            vm97flags(var)(7 : 7) = lEx%vm_flags(6)(var + 1 : var + 1)
            vm97flags(var)(8 : 8) = lEx%vm_flags(7)(var + 1 : var + 1)
            vm97flags(var)(9 : 9) = lEx%vm_flags(8)(var + 1 : var + 1)
            call AddCharDatumToDataline(trim(vm97flags(var)), dataline, EddyProProj%err_label)
        end do
    end if

    !> Uncomment to reintroduce flags for last 3 tests
    call AddCharDatumToDataline(lEx%vm_tlag_hf, dataline, separator)
    call AddCharDatumToDataline(lEx%vm_tlag_sf, dataline, separator)
    call AddCharDatumToDataline(lEx%vm_aoa_hf, dataline, separator)
    call AddCharDatumToDataline(lEx%vm_nshw_hf, dataline, separator)

    !> Write second string from Chunks
    call AddDatum(dataline, fluxnetChunks%s(2), separator)

    !> Foken's QC details
    call AddFloatDatumToDataline(lEx%TAU_SS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%H_SS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%FC_SS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%FH2O_SS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%FCH4_SS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%FGS4_SS, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%U_ITC, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%W_ITC, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%TS_ITC, dataline, EddyProProj%err_label)

    !> Write second string from Chunks
    !> FK04_ST_FLAG_W_U thru ...
    call AddDatum(dataline, fluxnetChunks%s(3), separator)

    !> LI-COR's IRGAs diagnostics breakdown
    do i = 1, 29
        call AddFloatDatumToDataline(lEx%licor_flags(i), dataline, EddyProProj%err_label)
    end do

    !> AGC/RSSI
    call AddFloatDatumToDataline(lEx%agc72, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%agc75, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rssi77, dataline, EddyProProj%err_label)

    !> Write third string from Chunks
    !> WBOOST_APPLIED thru AXES_ROTATION_METHOD
    call AddDatum(dataline, fluxnetChunks%s(4), separator)

    !> Rotation angles
    call AddFloatDatumToDataline(lEx%yaw, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%pitch, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%roll, dataline, EddyProProj%err_label)

    !> Detrending method and time constant
    call AddIntDatumToDataline(lEx%det_meth_int, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%det_timec, dataline, EddyProProj%err_label)

    !> Write forth string from Chunks
    !> TIMELAG_DETECTION_METHOD thru FOOTPRINT_MODEL
    call AddDatum(dataline, fluxnetChunks%s(5), separator)

    !> Metadata
    call AddIntDatumToDataline(lEx%logger_swver%major, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(lEx%logger_swver%minor, dataline, EddyProProj%err_label)
    call AddIntDatumToDataline(lEx%logger_swver%revision, dataline, EddyProProj%err_label)
    !>> Site info
    call AddFloatDatumToDataline(lEx%lat, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%lon, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%alt, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%canopy_height, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%disp_height, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%rough_length, dataline, EddyProProj%err_label)
    !>> Acquisition setup
    call AddFloatDatumToDataline(lEx%file_length, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%ac_freq, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%avrg_length, dataline, EddyProProj%err_label)
    !>> Master sonic height and north offset
    call AddDatum(dataline, trim(lEx%instr(sonic)%firm), separator)
    call AddDatum(dataline, trim(lEx%instr(sonic)%model), separator)
    call AddFloatDatumToDataline(lEx%instr(sonic)%height, dataline, EddyProProj%err_label)
    call AddDatum(dataline, lEx%instr(sonic)%wformat, separator)
    call AddDatum(dataline, lEx%instr(sonic)%wref, separator)
    call AddFloatDatumToDataline(lEx%instr(sonic)%north_offset, dataline, EddyProProj%err_label)
    call AddFloatDatumToDataline(lEx%instr(sonic)%hpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(lEx%instr(sonic)%vpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(lEx%instr(sonic)%tau, dataline, EddyProProj%err_label)

    !>> irgas
    do igas = ico2, igas4
        call AddDatum(dataline, trim(lEx%instr(igas)%firm), separator)
        call AddDatum(dataline, trim(lEx%instr(igas)%model), separator)
        call AddFloatDatumToDataline(lEx%instr(igas)%nsep, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%esep, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%vsep, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tube_l, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tube_d, dataline, EddyProProj%err_label, gain=1d3, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tube_f, dataline, EddyProProj%err_label, gain=6d4, offset=0d0)
        if (igas == ih2o) then
            call AddFloatDatumToDataline(lEx%instr(igas)%kw, dataline, EddyProProj%err_label)
            call AddFloatDatumToDataline(lEx%instr(igas)%ko, dataline, EddyProProj%err_label)
        end if
        call AddFloatDatumToDataline(lEx%instr(igas)%hpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%vpath_length, dataline, EddyProProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tau, dataline, EddyProProj%err_label)
    end do

    !> Custom variables

    call AddIntDatumToDataline(lEx%ncustom, dataline, EddyProProj%err_label)
    if (lEx%ncustom > 0) then
        do i = 1, lEx%ncustom
            call AddFloatDatumToDataline(lEx%user_var(i), dataline, EddyProProj%err_label)
        end do
    end if

    !> Write sisxth string from Chunks
    !> Biomet data
    call AddDatum(dataline, fluxnetChunks%s(6), separator)

    !> Replace error codes with user-defined error code
    dataline = replace2(dataline, ',-9999,', ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',NaN,',   ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',+Inf,', ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',-Inf,', ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',Inf,', ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',+Infinity,', ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',-Infinity,', ',' // trim(EddyProProj%err_label) // ',')
    dataline = replace2(dataline, ',Infinity,', ',' // trim(EddyProProj%err_label) // ',')

    write(uflxnt, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutFluxnetFcc
