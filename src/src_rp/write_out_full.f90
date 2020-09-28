!***************************************************************************
! write_out_full.f90
! ------------------
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
! \brief       Write all results on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFull(init_string, PeriodRecords, PeriodActualRecords)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    integer, intent(in) :: PeriodRecords
    integer, intent(in) :: PeriodActualRecords
    !> local variables
    integer :: var
    integer :: gas
!    integer :: prof
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum
    include '../src_common/interfaces.inc'

    !> Preliminary file and timestamp information
    call clearstr(dataline)
    call AddDatum(dataline, trim(adjustl(init_string)), separator)
    if (Stats%daytime) then
        call AddDatum(dataline, '1', separator)
    else
        call AddDatum(dataline, '0', separator)
    endif
    call WriteDatumInt(PeriodRecords, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(PeriodActualRecords, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> Corrected fluxes (Level 3)
    !> Tau
    call WriteDatumFloat(Flux3%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (RUsetup%meth /= 'none') then
        call WriteDatumFloat(Essentials%rand_uncer(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> H
    call WriteDatumFloat(Flux3%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(QCFlag%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (RUsetup%meth /= 'none') then
        call WriteDatumFloat(Essentials%rand_uncer(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> LE
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux3%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer_LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> Gases
    if(OutVarPresent(co2)) then
        call WriteDatumFloat(Flux3%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(co2), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux3%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    if(OutVarPresent(ch4)) then
        call WriteDatumFloat(Flux3%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(ch4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    if(OutVarPresent(gas4)) then
        call WriteDatumFloat(Flux3%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (RUsetup%meth /= 'none') then
            call WriteDatumFloat(Essentials%rand_uncer(gas4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> storage fluxes
    call WriteDatumFloat(Stor%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Stor%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            call WriteDatumFloat(Stor%of(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

    !> vertical advection fluxes
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            if (Stats5%Mean(w) /= error .and. Stats%d(gas) >= 0d0) then
                if (gas /= h2o) then
                    call WriteDatumFloat(Stats5%Mean(w) * Stats%d(gas) * 1d3, datum, EddyProProj%err_label)
                    call AddDatum(dataline, datum, separator)
                else
                    call WriteDatumFloat(Stats5%Mean(w) * Stats%d(gas), datum, EddyProProj%err_label)
                    call AddDatum(dataline, datum, separator)
                end if
            else
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

    !> Gas concentrations, densities and timelags
    do gas = co2, gas4
        if (OutVarPresent(gas)) then
            call WriteDatumFloat(Stats%d(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(Stats%chi(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(Stats%r(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(Essentials%used_timelag(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (Essentials%def_tlag(gas)) then
                call AddDatum(dataline, '1', separator)
            else
                call AddDatum(dataline, '0', separator)
            endif
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, '9', separator)
        end if
    end do

    !> Air properties
    call WriteDatumFloat(Stats7%Mean(ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Ta, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats%Pr, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(RHO%a, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (RHO%a /= 0d0 .and. RHO%a /= error) then
        call WriteDatumFloat(Ambient%RhoCp / RHO%a, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    call WriteDatumFloat(Ambient%Va, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (Flux3%h2o /= error) then
        call WriteDatumFloat(Flux3%h2o * h2o_to_ET, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    call WriteDatumFloat(RHO%w, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%e, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%es, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Q, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats%RH, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%VPD, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Td, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> Unrotated and rotated wind components
    call WriteDatumFloat(Stats4%Mean(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats4%Mean(v), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats4%Mean(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats5%Mean(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats5%Mean(v), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats5%Mean(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%WS, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%MWS, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats4%wind_dir, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> rotation angles
    call WriteDatumFloat(Essentials%yaw, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%pitch, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Essentials%roll, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> turbulence
    call WriteDatumFloat(Ambient%us, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Stats7%TKE, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%L, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%zL, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%bowen, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(Ambient%Ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> footprint
    if (Meth%foot /= 'none') then
        select case(foot_model_used(1:len_trim(foot_model_used)))
            case('kljun_04')
            call AddDatum(dataline, '0', separator)
            case('kormann_meixner_01')
            call AddDatum(dataline, '1', separator)
            case('hsieh_00')
            call AddDatum(dataline, '2', separator)
        end select
        call WriteDatumFloat(Foot%peak, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%offset, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x10, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x30, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x50, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x70, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x90, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, EddyProProj%err_label, separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> Uncorrected fluxes (Level 0)
    !> Tau
    call WriteDatumFloat(Flux0%tau, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> H
    call WriteDatumFloat(Flux0%H, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumFloat(BPCF%of(w_ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> LE
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux0%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    !> Gases
    if(OutVarPresent(co2)) then
        call WriteDatumFloat(Flux0%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_co2), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if(OutVarPresent(h2o)) then
        call WriteDatumFloat(Flux0%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if(OutVarPresent(ch4)) then
        call WriteDatumFloat(Flux0%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_ch4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if(OutVarPresent(gas4)) then
        call WriteDatumFloat(Flux0%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> Vickers and Mahrt 97 hard flags
    call AddDatum(dataline, '8'//CharHF%sr(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%ar(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%do(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%al(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%sk(2:9), separator)
    call AddDatum(dataline, '8'//CharSF%sk(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%ds(2:9), separator)
    call AddDatum(dataline, '8'//CharSF%ds(2:9), separator)
    call AddDatum(dataline, '8'//CharHF%tl(6:9), separator)
    call AddDatum(dataline, '8'//CharSF%tl(6:9), separator)
    call AddDatum(dataline, '8'//CharHF%aa(9:9), separator)
    call AddDatum(dataline, '8'//CharHF%ns(9:9), separator)

    !> Spikes for EddyPro variables
    call WriteDatumInt(Essentials%e2spikes(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Essentials%e2spikes(v), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Essentials%e2spikes(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(Essentials%e2spikes(ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do var = co2, gas4
        if(OutVarPresent(var)) then
            call WriteDatumInt(Essentials%e2spikes(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do

    !> LI-COR's diagnostic flags
    if (Diag7200%present) then
        call WriteDatumInt(Diag7200%head_detect, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%t_out, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%t_in, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%aux_in, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%delta_p, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%chopper, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%detector, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%pll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7200%sync, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (Diag7500%present) then
        call WriteDatumInt(Diag7500%chopper, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7500%detector, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7500%pll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7500%sync, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (Diag7700%present) then
        call WriteDatumInt(Diag7700%not_ready, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%no_signal, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%re_unlocked, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_temp, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%laser_temp_unregulated, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%block_temp_unregulated, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%motor_spinning, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%pump_on, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%top_heater_on, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bottom_heater_on, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%calibrating, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%motor_failure, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc1, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%bad_aux_tc3, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(Diag7700%box_connected, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> AGCs and RSSIs for LI-7200 and LI-7500
    if (Diag7200%present) then
        call WriteDatumInt(nint(Essentials%AGC72), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if
    if (Diag7500%present) then
        call WriteDatumInt(nint(Essentials%AGC75), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    elseif(EddyProProj%fix_out_format) then
        call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
    end if

    !> Variances
    do var = u, ts
        call WriteDatumFloat(Stats%Cov(var, var), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end do
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            call WriteDatumFloat(Stats%Cov(gas, gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
    end do
    !> w-covariances
    call WriteDatumFloat(Stats%Cov(w, ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    do gas = co2, gas4
        if(OutVarPresent(gas)) then
            call WriteDatumFloat(Stats%Cov(w, gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, &
            trim(adjustl(EddyProProj%err_label)), separator)
        end if
    enddo

    !> Mean values of user variables
    if (NumUserVar > 0) then
        do var = 1, NumUserVar
            call WriteDatumFloat(UserStats%Mean(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
    end if

    write(uflx, '(a)') dataline(1:len_trim(dataline) - 1)

end subroutine WriteOutFull
