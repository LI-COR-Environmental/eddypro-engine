!***************************************************************************
! write_out_files.f90
! -------------------
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
    character(10000) :: dataline

    !> local variables
    integer :: var
    integer :: i
    integer :: gas
    integer :: igas
    character(64) :: datum

    !>***************************************************************
    !>***************************************************************

    !>Write out full output file (main output)
    if (EddyProProj%out_full) then
        call clearstr(dataline)
        !> Preliminary file and timestamp information
        call AddDatum(dataline, lEx%fname(1:len_trim(lEx%fname)), separator)
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
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end if
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if

        !> storage
        call WriteDatumFloat(lEx%Stor%H, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(lEx%Stor%LE, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        do gas = co2, n2o
            if(fcc_var_present(gas)) then
                call WriteDatumFloat(lEx%Stor%of(gas), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                    call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                end if
            elseif(EddyProProj%fix_out_format) then
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        call WriteDatumFloat(lEx%Va, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
        if (Flux3%h2o /= error) then
            call WriteDatumFloat(Flux3%h2o * 0.0648d0, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        !> Gases
        if(fcc_var_present(co2)) then
            call WriteDatumFloat(lEx%Flux0%co2, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_co2), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        if(fcc_var_present(h2o)) then
            call WriteDatumFloat(lEx%Flux0%h2o, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_h2o), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        if(fcc_var_present(ch4)) then
            call WriteDatumFloat(lEx%Flux0%ch4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_ch4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        if(fcc_var_present(gas4)) then
            call WriteDatumFloat(lEx%Flux0%gas4, datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
            call WriteDatumFloat(BPCF%of(w_gas4), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end do
        end if
        if (Diag7500%present) then
            do i = 10, 13
                call WriteDatumInt(nint(lEx%licor_flags(i)), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        elseif(EddyProProj%fix_out_format) then
            do i = 1, 4
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end do
        end if
        if (Diag7700%present) then
            do i = 14, 29
                call WriteDatumInt(nint(lEx%licor_flags(i)), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        elseif(EddyProProj%fix_out_format) then
            do i = 1, 16
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
            end do
        end if

        !> AGCs
        if (Diag7200%present) then
            call WriteDatumInt(nint(lEx%agc72), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
        if (Diag7500%present) then
            call WriteDatumInt(nint(lEx%agc75), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        elseif(EddyProProj%fix_out_format) then
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
!        if (Diag7700%present) then
!            call WriteDatumInt(nint(lEx%rssi77), datum, EddyProProj%err_label)
!            call AddDatum(dataline, datum, separator)
!        elseif(EddyProProj%fix_out_format) then
!            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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
                call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
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

    !> write to GHG-Europe style output file
    if (EddyProProj%out_ghg_eu) then
        call clearstr(dataline)
        !> Preliminary timestmap information
        write(datum, *) lEx%fname(1:len_trim(lEx%fname))
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%date(1:10)
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%time(1:5)
        call AddDatum(dataline, datum, separator)

        !> Air properties
        write(datum, *) lEx%Ta - 273.16d0     !< celsius
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%Pa * 1d-3           !< kPa
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%RH                  !< %
        call AddDatum(dataline, datum, separator)

        !> Gas concentrations
        do gas = co2, gas4
            if (fcc_var_present(gas)) then
                write(datum, *) lEx%chi(gas)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do gas = ch4, gas4
            if (fcc_var_present(gas)) then
                write(datum, *) lEx%chi(gas) * 1d3
                call AddDatum(dataline, datum, separator)
            end if
        end do

        !> Corrected fluxes (Level 3)
        !> Tau
        write(datum, *) Flux3%tau
        call AddDatum(dataline, datum, separator)
        write(datum, *) QCFlag%tau
        call AddDatum(dataline, datum, separator)
        !> H
        write(datum, *) Flux3%H
        call AddDatum(dataline, datum, separator)
        write(datum, *) QCFlag%H
        call AddDatum(dataline, datum, separator)
        !> LE
        if(fcc_var_present(h2o)) then
            write(datum, *) Flux3%LE
            call AddDatum(dataline, datum, separator)
            write(datum, *) QCFlag%h2o
            call AddDatum(dataline, datum, separator)
        end if
        !> Gases
        if(fcc_var_present(co2)) then
            write(datum, *) Flux3%co2
            call AddDatum(dataline, datum, separator)
            write(datum, *) QCFlag%co2
            call AddDatum(dataline, datum, separator)
        end if
        if(fcc_var_present(h2o)) then
            write(datum, *) Flux3%h2o
            call AddDatum(dataline, datum, separator)
            write(datum, *) QCFlag%h2o
            call AddDatum(dataline, datum, separator)
        end if
        if(fcc_var_present(ch4)) then
            write(datum, *) Flux3%ch4 * 1d3  !< expressed here in nmol+1m-2s-1
            call AddDatum(dataline, datum, separator)
            write(datum, *) QCFlag%ch4
            call AddDatum(dataline, datum, separator)
        end if
        if(fcc_var_present(gas4)) then
            write(datum, *) Flux3%gas4 * 1d3  !< expressed here in nmol+1m-2s-1
            call AddDatum(dataline, datum, separator)
            write(datum, *) QCFlag%gas4
            call AddDatum(dataline, datum, separator)
        end if

        !> Storage
        write(datum, *) lEx%Stor%H
        call AddDatum(dataline, datum, separator)
        if(fcc_var_present(h2o)) then
            write(datum, *) lEx%Stor%LE
            call AddDatum(dataline, datum, separator)
        end if
        if(fcc_var_present(co2)) then
            write(datum, *) lEx%Stor%of(co2)
            call AddDatum(dataline, datum, separator)
        end if

        !> Turbulence
        write(datum, *) lEx%WS
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%WD
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%ustar
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%L
        call AddDatum(dataline, datum, separator)
        write(datum, *) lEx%zL
        call AddDatum(dataline, datum, separator)

        !> footprint
        write(datum, *) Foot%peak
        call AddDatum(dataline, datum, separator)
        write(datum, *) Foot%x70
        call AddDatum(dataline, datum, separator)
        write(datum, *) Foot%x90
        call AddDatum(dataline, datum, separator)

        write(ughgeu, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>****************************************************************
    !>****************************************************************

    !> write to metadata output file
    if (EddyProProj%out_md) then
        call clearstr(dataline)
        !> Preliminary timestmap information
        write(datum, *) lEx%fname(1:len_trim(lEx%fname))
        call AddDatum(dataline, datum, separator)
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
end subroutine WriteOutputFiles
