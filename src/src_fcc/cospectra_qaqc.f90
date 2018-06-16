!*******************************************************************************
! cospectra_qaqc.f90
! ------------------
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
!*******************************************************************************
!
! \brief       Set (co)spectra to error if user-provided quality criteria \n
!              are not met
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine CospectraQAQC(BinSpec, BinCosp, nrow, lEx, &
    BinCospForStable, BinCospForUnstable, skip_spectra, skip_cospectra)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(ExType), intent(in) :: lEx
    type(SpectraSetType), intent(inout) :: BinSpec(nrow)
    type(SpectraSetType), intent(inout) :: BinCosp(nrow)
    type(SpectraSetType) :: BinCospForStable(nrow)
    type(SpectraSetType) :: BinCospForUnstable(nrow)
    logical, intent(out) :: skip_spectra
    logical, intent(out) :: skip_cospectra
    !> Local variables
    integer :: i
    character(9) :: hf_sr
    character(9) :: hf_do
    character(9) :: hf_sk, sf_sk
    character(9) :: hf_ds, sf_ds
    integer :: STFlg(GHGNumVar)
    integer :: DTFlg(GHGNumVar)
    integer :: qc_tau, qc_H, qc_co2, qc_h2o, qc_ch4, qc_gas4


    !> Initialization
    skip_spectra   = .false.
    skip_cospectra = .false.
    BinCospForStable = BinCosp
    BinCospForUnstable = BinCosp


    !> Razor blade on spectra and co-spectra for unstable case, \n
    !> based on corresponding fluxes.
    !> Unstable case
    if (dabs(lEx%Flux0%H) > FCCsetup%SA%min_un_H &
        .and. dabs(lEx%Flux0%H) < FCCsetup%SA%max_H) then

        if (dabs(lEx%Flux0%LE) < FCCsetup%SA%min_un_LE &
            .or. dabs(lEx%Flux0%LE) > FCCsetup%SA%max_LE) then
            BinSpec%of(h2o) = error
            BinCospForUnstable%of(h2o) = error
        end if
        if (dabs(lEx%Flux0%co2) < FCCsetup%SA%min_un_co2 &
            .or. dabs(lEx%Flux0%co2) > FCCsetup%SA%max_co2)  then
            BinSpec%of(co2) = error
            BinCospForUnstable%of(co2) = error
        end if
        if (dabs(lEx%Flux0%ch4) < FCCsetup%SA%min_un_ch4 &
            .or. dabs(lEx%Flux0%ch4) > FCCsetup%SA%max_ch4)  then
            BinSpec%of(ch4) = error
            BinCospForUnstable%of(ch4) = error
        end if
        if (dabs(lEx%Flux0%gas4) < FCCsetup%SA%min_un_gas4 &
            .or. dabs(lEx%Flux0%gas4) > FCCsetup%SA%max_gas4) then
            BinSpec%of(gas4) = error
            BinCospForUnstable%of(gas4) = error
        end if
        if (dabs(lEx%Flux0%LE) < FCCsetup%SA%min_un_LE .and. &
            dabs(lEx%Flux0%co2) < FCCsetup%SA%min_un_co2 .and. &
            dabs(lEx%Flux0%ch4) < FCCsetup%SA%min_un_ch4 .and. &
            dabs(lEx%Flux0%gas4) < FCCsetup%SA%min_un_gas4) then
            BinSpec = ErrSpec
            BinCospForUnstable = ErrSpec
            skip_spectra = .true.
        end if
    else
        BinSpec = ErrSpec
        BinCospForUnstable = ErrSpec
        skip_spectra = .true.
    end if

    !> Filter co-spectra for u*
    if (lEx%ustar < FCCsetup%SA%min_un_ustar &
        .or. lEx%ustar > FCCsetup%SA%max_ustar) then
        BinCospForUnstable = ErrSpec
        BinSpec = ErrSpec
    end if

    !> Stable case
    if (dabs(lEx%Flux0%H) > FCCsetup%SA%min_st_H &
        .and. dabs(lEx%Flux0%H) < FCCsetup%SA%max_H) then

        if (dabs(lEx%Flux0%LE) < FCCsetup%SA%min_st_LE &
            .or. dabs(lEx%Flux0%LE) > FCCsetup%SA%max_LE) &
            BinCospForStable%of(h2o) = error

        if (dabs(lEx%Flux0%co2) < FCCsetup%SA%min_st_co2 &
            .or. dabs(lEx%Flux0%co2) > FCCsetup%SA%max_co2)  &
            BinCospForStable%of(co2) = error

        if (dabs(lEx%Flux0%ch4) < FCCsetup%SA%min_st_ch4 &
            .or. dabs(lEx%Flux0%ch4) > FCCsetup%SA%max_ch4)  &
            BinCospForStable%of(ch4) = error

        if (dabs(lEx%Flux0%gas4) < FCCsetup%SA%min_un_gas4 &
            .or. dabs(lEx%Flux0%gas4) > FCCsetup%SA%max_gas4) &
            BinCospForStable%of(gas4) = error

        if (dabs(lEx%Flux0%LE) < FCCsetup%SA%min_st_LE .and. &
            dabs(lEx%Flux0%co2) < FCCsetup%SA%min_un_co2 .and. &
            dabs(lEx%Flux0%ch4) < FCCsetup%SA%min_un_ch4 .and. &
            dabs(lEx%Flux0%gas4) < FCCsetup%SA%min_un_gas4) then
            BinCospForStable = ErrSpec
            skip_cospectra = .true.
        end if
    else
        BinCospForStable = ErrSpec
        skip_cospectra = .true.
    end if

    !> Filter co-spectra for u*
    if (lEx%ustar < FCCsetup%SA%min_st_ustar &
        .or. lEx%ustar > FCCsetup%SA%max_ustar) then
        BinCospForStable = ErrSpec
        skip_cospectra = .true.
    end if

    !> Filter based on results of Vickers and Mahrt (1997) quality tests
    !> if requested
    if (FCCsetup%SA%filter_cosp_by_vm_flags) then
        hf_sr(1:8) = lEx%vm_flags(1)(2:9)
        hf_do(1:8) = lEx%vm_flags(3)(2:9)

        hf_sk(1:8) = lEx%vm_flags(5)(2:9)
        sf_sk(1:8) = lEx%vm_flags(6)(2:9)

        hf_ds(1:8) = lEx%vm_flags(7)(2:9)
        sf_ds(1:8) = lEx%vm_flags(8)(2:9)

        !> If vertical wind speed is flagged, all cospectra are eliminated
        if (hf_sr(w:w) == '1' .or. hf_do(w:w) == '1' &
            .or. hf_sk(w:w) == '1' .or. hf_ds(w:w) == '1') then
            BinCospForUnstable = ErrSpec
        end if

        !> Elimination of individual (co)spectra based on the flags on
        !> the relevant variable
        do i = u, gas4
            if (hf_sr(i:i) == '1' .or. hf_do(i:i) == '1' &
                .or. hf_sk(i:i) == '1' .or. hf_ds(i:i) == '1') then
                BinSpec%of(i) = error
                BinCospForUnstable%of(i) = error
            end if
        end do
    end if

    !> Filter based on results of Foken quality tests if requested.
    !> Regardless of user's choice on how to flag fluxes, here the 0/1/2 scheme
    !> of Mauder and Foken 2004 is used
    if (FCCsetup%SA%foken_lim >= 0) then
        !> Partial flags
        !> Stationarity flags
        call PartialFlagLF(nint(lEx%FC_SS), STFlg(w_co2))
        call PartialFlagLF(nint(lEx%FH2O_SS), STFlg(w_h2o))
        call PartialFlagLF(nint(lEx%FCH4_SS), STFlg(w_ch4))
        call PartialFlagLF(nint(lEx%FGS4_SS), STFlg(w_gas4))
        call PartialFlagLF(nint(lEx%H_SS),  STFlg(w_ts))
        call PartialFlagLF(nint(lEx%TAU_SS),   STFlg(w_u))
        !> Developed turbulence flags
        call PartialFlagLF(nint(lEx%U_ITC), DTFlg(u))
        call PartialFlagLF(nint(lEx%W_ITC), DTFlg(w))
        call PartialFlagLF(nint(lEx%TS_ITC), DTFlg(ts))
        DTFlg(u)  = max(DTFlg(u),  DTFlg(w))

        !> Composite flags
        call GTK2Flag(STFlg(w_u),   DTFlg(u), qc_tau)
        call GTK2Flag(STFlg(w_ts),  DTFlg(w), qc_H)
        call GTK2Flag(STFlg(w_co2), DTFlg(w), qc_co2)
        call GTK2Flag(STFlg(w_h2o), DTFlg(w), qc_h2o)
        call GTK2Flag(STFlg(w_ch4), DTFlg(w), qc_ch4)
        call GTK2Flag(STFlg(w_gas4), DTFlg(w), qc_gas4)

        !> Actual (co)spectra elimination
        if (qc_H < FCCsetup%SA%foken_lim &
            .and. qc_tau < FCCsetup%SA%foken_lim) then
            if (qc_h2o >= FCCsetup%SA%foken_lim) then
                BinSpec%of(h2o) = error
                BinCospForUnstable%of(h2o) = error
            end if
            if (qc_co2 >= FCCsetup%SA%foken_lim)  then
                BinSpec%of(co2) = error
                BinCospForUnstable%of(co2) = error
            end if
            if (qc_ch4 >= FCCsetup%SA%foken_lim)  then
                BinSpec%of(ch4) = error
                BinCospForUnstable%of(ch4) = error
            end if
            if (qc_gas4 >= FCCsetup%SA%foken_lim) then
                BinSpec%of(gas4) = error
                BinCospForUnstable%of(gas4) = error
            end if
            if (qc_h2o >= FCCsetup%SA%foken_lim &
                .and. qc_co2 >= FCCsetup%SA%foken_lim &
                .and. qc_ch4 >= FCCsetup%SA%foken_lim &
                .and. qc_gas4 >= FCCsetup%SA%foken_lim) then
                BinSpec = ErrSpec
                BinCospForUnstable = ErrSpec
                skip_spectra = .true.
            end if
        else
            BinSpec = ErrSpec
            BinCospForUnstable = ErrSpec
            skip_spectra = .true.
        end if
    end if

    !> For time sorted cospectra use milder filtering
    BinCosp = BinCospForStable

end subroutine CospectraQAQC
