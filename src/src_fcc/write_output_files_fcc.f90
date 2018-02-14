!***************************************************************************
! write_output_files.f90
! ----------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutputFiles(lEx)
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
    character(DatumLen) :: datum
    character(14) :: iso_basic
    include '../src_common/interfaces_1.inc'

    !>***************************************************************
    !>***************************************************************

    !>Write out full output file (main output)
    if (EddyProProj%out_full) then
        call clearstr(dataline)
        !> Preliminary file and timestamp information
        ! call AddDatum(dataline, lEx%fname(1:len_trim(lEx%fname)), separator)
        call AddDatum(dataline, lEx%date(1:10), separator)
        call AddDatum(dataline, lEx%time(1:5), separator)
        call WriteDatumFloat(float_doy, datum, EddyProProj%err_label)
        call stripstr(datum)  !< Added to fix a strange behaviour
        call AddDatum(dataline, datum(1:index(datum, '.') + 3), separator)
        if (lEx%daytime) then
            call AddDatum(dataline, '1', separator)
        else
            call AddDatum(dataline, '0', separator)
        endif
        call WriteDatumInt(lEx%file_records, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%used_records, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Corrected fluxes (Level 3)
        !> Tau
        call WriteDatumFloat(Flux3%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (FCCMetadata%ru) then
            call WriteDatumFloat(lEx%rand_uncer(u), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> H
        call WriteDatumFloat(Flux3%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (FCCMetadata%ru) then
            call WriteDatumFloat(lEx%rand_uncer(ts), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> LE
        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(Flux3%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (FCCMetadata%ru) then
                call WriteDatumFloat(lEx%rand_uncer_LE, datum, EddyProProj%err_label)
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
        if(fcc_var_present(co2)) then
            call WriteDatumFloat(Flux3%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (FCCMetadata%ru) then
                call WriteDatumFloat(lEx%rand_uncer(co2), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(Flux3%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (FCCMetadata%ru) then
                call WriteDatumFloat(lEx%rand_uncer(h2o), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(ch4)) then
            call WriteDatumFloat(Flux3%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (FCCMetadata%ru) then
                call WriteDatumFloat(lEx%rand_uncer(ch4), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(gas4)) then
            call WriteDatumFloat(Flux3%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            if (FCCMetadata%ru) then
                call WriteDatumFloat(lEx%rand_uncer(gas4), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> storage
        call WriteDatumFloat(lEx%Stor%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(lEx%Stor%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        do gas = co2, n2o
            if(fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%Stor%of(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do

        !> vertical advection fluxes
        do gas = co2, n2o
            if(fcc_var_present(gas)) then
                if (lEx%rot_w /= error .and. lEx%d(gas) >= 0d0) then
                    if (gas /= h2o) then
                        call WriteDatumFloat(lEx%rot_w * lEx%d(gas) * 1d3, datum, EddyProProj%err_label)
                        call AddDatum(dataline, datum, separator)
                    else
                        call WriteDatumFloat(lEx%rot_w * lEx%d(gas), datum, EddyProProj%err_label)
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
        do gas = co2, n2o
            if (fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%d(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(lEx%chi(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(lEx%r(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                call WriteDatumFloat(lEx%tlag(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
                if (lEx%def_tlag(gas)) then
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
        call WriteDatumFloat(lEx%Ts, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Ta, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Pa, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%RHO%a, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (lEx%RHO%a /= 0d0 .and. lEx%RHO%a /= error) then
            call WriteDatumFloat(lEx%RhoCp /lEx%RHO%a, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        call WriteDatumFloat(lEx%Va, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (Flux3%h2o /= error) then
            call WriteDatumFloat(Flux3%h2o * 0.0648d0, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        call WriteDatumFloat(lEx%RHO%w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%e, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%es, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Q, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%RH, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%VPD, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Tdew, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Unrotated and rotated wind components
        call WriteDatumFloat(lEx%unrot_u, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%unrot_v, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%unrot_w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rot_u, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rot_v, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rot_w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%WS, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%MWS, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%WD, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> rotation angles
        call WriteDatumFloat(lEx%yaw, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%pitch, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%roll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> turbulence
        call WriteDatumFloat(lEx%ustar, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%TKE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%L, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%zL, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%bowen, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Tstar, datum, EddyProProj%err_label)
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
            call AddDatum(dataline, '9', separator)
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
        call WriteDatumFloat(lEx%Flux0%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> H
        call WriteDatumFloat(lEx%Flux0%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(BPCF%of(w_ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> LE
        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(lEx%Flux0%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        !> Gases
        if(fcc_var_present(co2)) then
            call WriteDatumFloat(lEx%Flux0%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_co2), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(lEx%Flux0%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(fcc_var_present(ch4)) then
            call WriteDatumFloat(lEx%Flux0%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_ch4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(fcc_var_present(gas4)) then
            call WriteDatumFloat(lEx%Flux0%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_gas4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Vickers and Mahrt 97 flags
        do i = 1, 12
            write(datum, *) lEx%vm_flags(i)
            call AddDatum(dataline, datum, separator)
        end do

        !> Spikes for EddyPro variables
        call WriteDatumInt(lEx%spikes(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%spikes(v), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%spikes(w), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%spikes(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do var = co2, gas4
            if(fcc_var_present(var)) then
                call WriteDatumInt(lEx%spikes(var), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do

        !> LI-COR's diagnostic flags
        if (Diag7200%present) then
            do i = 1, 9
                call WriteDatumInt(nint(lEx%licor_flags(i)), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        elseif(EddyProProj%fix_out_format) then
            do i = 1, 9
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end do
        end if
        if (Diag7500%present) then
            do i = 10, 13
                call WriteDatumInt(nint(lEx%licor_flags(i)), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        elseif(EddyProProj%fix_out_format) then
            do i = 1, 4
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end do
        end if
        if (Diag7700%present) then
            do i = 14, 29
                call WriteDatumInt(nint(lEx%licor_flags(i)), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        elseif(EddyProProj%fix_out_format) then
            do i = 1, 16
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end do
        end if

        !> AGCs
        if (Diag7200%present) then
            call WriteDatumInt(nint(lEx%agc72), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if (Diag7500%present) then
            call WriteDatumInt(nint(lEx%agc75), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
!        if (Diag7700%present) then
!            call WriteDatumInt(nint(lEx%rssi77), datum, EddyProProj%err_label)
!            call AddDatum(dataline, datum, separator)
!        elseif(EddyProProj%fix_out_format) then
!            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
!        end if

        !> Variances
        do var = u, ts
            call WriteDatumFloat(lEx%var(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do gas = co2, gas4
            if(fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%var(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        end do
        !> w-covariances
        call WriteDatumFloat(lEx%cov_w(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do gas = co2, gas4
            if(fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%cov_w(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
            end if
        enddo

        !> Mean values of user variables
        if (NumUserVar > 0) then
            do var = 1, NumUserVar
                call WriteDatumFloat(lEx%user_var(var), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        end if

        write(uflx, '(a)')   dataline(1:len_trim(dataline) - 1)
    end if

    !>****************************************************************
    !>****************************************************************
    !> FLUXNET output
    if (EddyProProj%out_fluxnet) then
        call clearstr(dataline)

        !> derive ISO basic format timestamp
        iso_basic = lEx%date(1:4) // lEx%date(6:7) &
            // lEx%date(9:10) // lEx%time(1:2) // lEx%time(4:5) // '00'

        call clearstr(dataline)
        call AddDatum(dataline, trim(adjustl(iso_basic)), separator)

        !> Gas concentrations
        do gas = co2, h2o
            if (fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%chi(gas), datum, '-9999.')
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do gas = ch4, gas4
            if (fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%chi(gas) * 1d3, datum, '-9999.')
                call AddDatum(dataline, datum, separator)
            end if
        end do

        !> Corrected fluxes (Level 3)
        !> Tau
        call WriteDatumFloat(Flux3%tau, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%tau, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        !> H
        call WriteDatumFloat(Flux3%H, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(QCFlag%H, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        !> LE
        if(fcc_var_present(h2o)) then
            write(datum, *) Flux3%LE
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%h2o, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        end if
        !> Gases
        if(fcc_var_present(co2)) then
            call WriteDatumFloat(Flux3%co2, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%co2, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        end if
        if(fcc_var_present(ch4)) then
            call WriteDatumFloat(Flux3%ch4 * 1d3, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%ch4, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        end if
        if(fcc_var_present(gas4)) then
            call WriteDatumFloat(Flux3%gas4 * 1d3, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
            call WriteDatumInt(QCFlag%gas4, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        end if

        !> Turbulence
        call WriteDatumFloat(lEx%WD, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%WS, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%MWS, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        if(lEx%var(u) > 0d0) then
            call WriteDatumFloat(dsqrt(lEx%var(u)), datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(lEx%var(v) > 0d0) then
            call WriteDatumFloat(dsqrt(lEx%var(v)), datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        if(lEx%var(w) > 0d0) then
            call WriteDatumFloat(dsqrt(lEx%var(w)), datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if
        call WriteDatumFloat(lEx%ustar, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%L, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%zL, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%peak, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x70, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x80, datum, '-9999.')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Foot%x90, datum, '-9999.')
        call AddDatum(dataline, datum, separator)

        !> Ambient pressure in kPa
        if (lEx%Pa /= error) then
            call WriteDatumFloat(lEx%Pa * 1d-3, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> RH
        call WriteDatumFloat(lEx%RH, datum, '-9999.')
        call AddDatum(dataline, datum, separator)

        !> Ambient temperature in degC
        if (lEx%Ta /= error) then
            call WriteDatumFloat(lEx%Ta - 273.15d0, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> VPD in hPa
        if (lEx%VPD /= error) then
            call WriteDatumFloat(lEx%VPD * 1d-2, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Sonic temperature in degC
        if (lEx%Ts /= error) then
            call WriteDatumFloat(lEx%Ts - 273.15d0, datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Sonic temperature standard deviation
        if(lEx%var(ts) > 0d0) then
            call WriteDatumFloat(dsqrt(lEx%var(ts)), datum, '-9999.')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        write(ufnet_e, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>****************************************************************
    !>****************************************************************
    !> write to metadata output file
    if (EddyProProj%out_md) then
        call clearstr(dataline)
        !> Preliminary timestmap information
        ! write(datum, *) lEx%fname(1:len_trim(lEx%fname))
        ! call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%date(1:10)
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%time(1:5)
        call AddDatum(dataline, datum, separator)

        !> Site location and characteristics
        write(datum, *) lEx%lat
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%lon
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%alt
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%canopy_height
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%disp_height
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%rough_length
        call AddDatum(dataline, datum, separator)

        !> Acquisition setup
        write(datum, *) lEx%file_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%ac_freq
        call AddDatum(dataline, datum, separator)
        !> Master sonic height and north offset
        write(datum, *) lEx%instr(sonic)%firm(1:len_trim(lEx%instr(sonic)%firm))
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%model(1:len_trim(lEx%instr(sonic)%model))
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%height
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%wformat
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%wref
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%north_offset
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%hpath_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%vpath_length
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%instr(sonic)%tau
        call AddDatum(dataline, datum, separator)
        !> irgas
        do igas = ico2, igas4
            if (fcc_var_present(3 + igas)) then
                write(datum, *) lEx%instr(igas)%firm(1:len_trim(lEx%instr(igas)%firm))
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%model(1:len_trim(lEx%instr(igas)%model))
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%measure_type(3 + igas)
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%nsep
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%esep
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%vsep
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%tube_l
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%tube_d
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%tube_f
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%kw
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%ko
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%hpath_length
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%vpath_length
                call AddDatum(dataline, datum, separator)
                write(datum, *) lEx%instr(igas)%tau
                call AddDatum(dataline, datum, separator)
            end if
        end do

        write(umd,*) dataline(1:len_trim(dataline) - 1)
    end if

    !>****************************************************************
    !>****************************************************************
    !>Write out ICOS output file
    if (EddyProProj%out_icos) then
        call clearstr(dataline)
        !> Timestamp
        call AddDatum(dataline, trim(lEx%timestamp), separator)

        !> Potential radiation and daytime
        call WriteDatumFloat(lEx%RP, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%daytime_int, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Number of records
        call WriteDatumInt(lEx%nr_theor, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)        
        call WriteDatumInt(lEx%nr_files, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%nr_after_custom_flags, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%nr_after_wdf, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%nr(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do var = ts, gas4
            call WriteDatumInt(lEx%nr(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        call WriteDatumInt(lEx%nr_w(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do var = ts, gas4
            call WriteDatumInt(lEx%nr_w(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> Final fluxes
        call WriteDatumFloat(Flux3%Tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux3%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux3%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux3%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux3%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux3%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux3%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Random uncertainties
        call WriteDatumFloat(lEx%rand_uncer(u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rand_uncer(ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rand_uncer_LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do gas = co2, gas4
            call WriteDatumFloat(lEx%rand_uncer(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> Storage fluxes
        call WriteDatumFloat(lEx%Stor%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Stor%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do gas = co2, gas4
            call WriteDatumFloat(lEx%Stor%of(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> vertical advection fluxes
        do gas = co2, n2o
            if (lEx%rot_w /= error .and. lEx%d(gas) >= 0d0) then
                if (gas /= h2o) then
                    call WriteDatumFloat(lEx%rot_w * lEx%d(gas) * 1d3, datum, EddyProProj%err_label)
                    call AddDatum(dataline, datum, separator)
                else
                    call WriteDatumFloat(lEx%rot_w * lEx%d(gas), datum, EddyProProj%err_label)
                    call AddDatum(dataline, datum, separator)
                end if
            else
                call AddDatum(dataline, trim(EddyProProj%err_label), separator)
            end if
        end do
        
        !> Unrotated and rotated wind components
        call WriteDatumFloat(lEx%unrot_u, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%unrot_v, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%unrot_w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rot_u, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rot_v, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rot_w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%WS, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%MWS, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%WD, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Turbulence
        call WriteDatumFloat(lEx%ustar, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%TKE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%L, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%zL, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%bowen, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Tstar, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Thermodynamics
        call WriteDatumFloat(lEx%Ts, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Ta, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Pa, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%RH, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Va, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%RHO%a, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%RhoCp, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Water
        call WriteDatumFloat(lEx%RHO%w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%e, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%es, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Q, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%VPD, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Tdew, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Dry air
        call WriteDatumFloat(lEx%Pd, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%RHO%d, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Vd, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%lambda, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%sigma, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Gas concentrations/densities
        do gas = co2, gas4
            call WriteDatumInt(lEx%measure_type_int(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%d(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%r(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%chi(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> Time lags
        do gas = co2, gas4
            call WriteDatumFloat(lEx%act_tlag(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%used_tlag(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%nom_tlag(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%min_tlag(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%max_tlag(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> Stats
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%mean(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%median(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%Q1(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%Q3(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%Cov(var, var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%Skw(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = u, gas4
            call WriteDatumFloat(lEx%stats%Kur(var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        call WriteDatumFloat(lEx%stats%Cov(w, u), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do var = ts, gas4
            call WriteDatumFloat(lEx%stats%Cov(w, var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = h2o, gas4
            call WriteDatumFloat(lEx%stats%Cov(co2, var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        do var = ch4, gas4
            call WriteDatumFloat(lEx%stats%Cov(h2o, var), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        call WriteDatumFloat(lEx%stats%Cov(ch4, gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Footprint
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

        !> Fluxes Level 0 (uncorrected)
        call WriteDatumFloat(lEx%Flux0%L, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%zL, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%Tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> Fluxes Level 1 
        call WriteDatumFloat(Flux1%Tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux1%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux1%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux1%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux1%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux1%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux1%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !> Fluxes Level 2
        call WriteDatumFloat(Flux2%Tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux2%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux2%LE, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux2%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux2%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux2%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(Flux2%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Cell values
        call WriteDatumFloat(lEx%Tcell, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Pcell, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do gas = co2, gas4
            call WriteDatumFloat(lEx%Vcell(gas), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do
        call WriteDatumFloat(lEx%Flux0%E_co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%E_ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%E_gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%Hi_co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%Hi_h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%Hi_ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Flux0%Hi_gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Burba terms
        call WriteDatumFloat(lEx%Burba%h_bot, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Burba%h_top, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Burba%h_spar, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> LI-7700 multipliers
        call WriteDatumFloat(lEx%Mul7700%A, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Mul7700%B, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%Mul7700%C, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Spectral correction factors
        if(fcc_var_present(u)) then
            call WriteDatumFloat(BPCF%of(w_u), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(ts)) then
            call WriteDatumFloat(BPCF%of(w_ts), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(co2)) then
            call WriteDatumFloat(BPCF%of(w_co2), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(ch4)) then
            call WriteDatumFloat(BPCF%of(w_ch4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        if(fcc_var_present(gas4)) then
            call WriteDatumFloat(BPCF%of(w_gas4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, trim(adjustl(EddyProProj%err_label)), separator)
        end if

        !> Degraded covariances
        call WriteDatumFloat(lEx%degT%cov, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        do i = 1, 9
            call WriteDatumFloat(lEx%degT%dcov(i), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> Write first string from Chunks
        !> M_CUSTOM_FLAGS thru VM97_NSW_RNS
        call AddDatum(dataline, trim(icosChunks%s(1)), separator)
        !> VM97 flags and Foken's QC details
        do i = 1, 12
            call AddDatum(dataline, trim(lEx%vm_flags(i)), separator)
        end do
        call WriteDatumFloat(lEx%st_w_u, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%st_w_ts, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%st_w_co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%st_w_h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%st_w_ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%st_w_gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%dt_u, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%dt_w, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%dt_ts, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Write second string from Chunks
        !> FK04_ST_FLAG_W_U thru ...
        call AddDatum(dataline, icosChunks%s(2), separator)

        !> LI-COR's IRGAs diagnostics breakdown
        do i = 1, 29
            call WriteDatumFloat(lEx%licor_flags(i), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do

        !> AGC/RSSI
        call WriteDatumFloat(lEx%agc72, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%agc75, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rssi77, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Write third string from Chunks
        !> WBOOST_APPLIED thru AXES_ROTATION_METHOD
        call AddDatum(dataline, icosChunks%s(3), separator)

        !> Rotation angles
        call WriteDatumFloat(lEx%yaw, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%pitch, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%roll, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Detrending method and time constant
        call WriteDatumInt(lEx%det_meth_int, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%det_timec, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !> Write forth string from Chunks
        !> TIMELAG_DETECTION_METHOD thru FOOTPRINT_MODEL
        call AddDatum(dataline, icosChunks%s(4), separator)

        !> Metadata
        call WriteDatumInt(lEx%logger_swver%major, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%logger_swver%minor, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumInt(lEx%logger_swver%revision, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !>> Site info
        call WriteDatumFloat(lEx%lat, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%lon, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%alt, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%canopy_height, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%disp_height, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%rough_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !>> Acquisition setup
        call WriteDatumFloat(lEx%file_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%ac_freq, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%avrg_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        !>> Master sonic height and north offset
        call AddDatum(dataline, trim(lEx%instr(sonic)%firm), separator)
        call AddDatum(dataline, trim(lEx%instr(sonic)%model), separator)
        call WriteDatumFloat(lEx%instr(sonic)%height, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call AddDatum(dataline, lEx%instr(sonic)%wformat, separator)
        call AddDatum(dataline, lEx%instr(sonic)%wref, separator)
        call WriteDatumFloat(lEx%instr(sonic)%north_offset, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%instr(sonic)%hpath_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%instr(sonic)%vpath_length, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(lEx%instr(sonic)%tau, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)

        !>> irgas
        do igas = ico2, igas4
            call AddDatum(dataline, trim(lEx%instr(igas)%firm), separator)
            call AddDatum(dataline, trim(lEx%instr(igas)%model), separator)
            call WriteDatumFloat(lEx%instr(igas)%nsep, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%esep, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%vsep, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%tube_l, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%tube_d, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%tube_f, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%kw, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%ko, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%hpath_length, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%vpath_length, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(lEx%instr(igas)%tau, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        end do  

        !> Write fifth string from Chunks
        !> Custom variables and biomet data
        call AddDatum(dataline, icosChunks%s(5), separator)

        !> Replace NaN or -9999 with user-defined error code
        dataline = replace2(dataline, ',-9999,', ',' // trim(EddyProProj%err_label) // ',')
        dataline = replace2(dataline, ',NaN,',   ',' // trim(EddyProProj%err_label) // ',')
    end if
    write(uicos, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutputFiles
