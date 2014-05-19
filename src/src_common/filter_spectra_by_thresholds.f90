!***************************************************************************
! filter_spectra_by_thresholds.f90
! --------------------------------
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
! \brief       Retrieve "essentials" information from file, based on
!              timestamp provided on input
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FilterSpectraByThresholds(BinSpec, BinCosp, nrow, lEx, &
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


    !> Initialization
    skip_spectra   = .false.
    skip_cospectra = .false.
    BinCospForStable = BinCosp
    BinCospForUnstable = BinCosp

    !> Razor blade on spectra and cospectra for unstable case,based on corresponding flux values
    if (lEx%Flux0%H > FCCsetup%SA%trshld_H .and. dabs(lEx%Flux0%H) < MaxHFlux) then
        if (lEx%Flux0%LE  <= FCCsetup%SA%trshld_LE .or. dabs(lEx%Flux0%LE) > MaxLEFlux) then
            BinSpec%of(h2o) = error
            BinCospForUnstable%of(h2o) = error
        end if
        if (dabs(lEx%Flux0%co2)  <= FCCsetup%SA%trshld_co2 .or. dabs(lEx%Flux0%co2)  > MaxCO2Flux)  then
            BinSpec%of(co2)  = error
            BinCospForUnstable%of(co2)  = error
        end if
        if (dabs(lEx%Flux0%ch4)  <= FCCsetup%SA%trshld_ch4 .or. dabs(lEx%Flux0%ch4)  > MaxCH4Flux)  then
            BinSpec%of(ch4)  = error
            BinCospForUnstable%of(ch4)  = error
        end if
        if (dabs(lEx%Flux0%gas4) <= FCCsetup%SA%trshld_gas4 .or. dabs(lEx%Flux0%gas4)  > MaxGAS4Flux) then
            BinSpec%of(gas4) = error
            BinCospForUnstable%of(gas4) = error
        end if

        if (lEx%Flux0%LE  <= FCCsetup%SA%trshld_LE .and. &
            dabs(lEx%Flux0%co2)  <= FCCsetup%SA%trshld_co2 .and. &
            dabs(lEx%Flux0%ch4)  <= FCCsetup%SA%trshld_ch4 .and. &
            dabs(lEx%Flux0%gas4) <= FCCsetup%SA%trshld_gas4) then
            BinSpec = ErrSpec
            BinCospForUnstable = ErrSpec
            skip_spectra = .true.
        end if
    else
        BinSpec = ErrSpec
        BinCospForUnstable = ErrSpec
        skip_spectra = .true.
    end if
        !> Further filter on u* for cospectra
    if (lEx%ustar < minUstar) then
        BinCospForUnstable = ErrSpec
    end if

    !> Milder razor blade for stable cospectra
    if (dabs(lEx%Flux0%H) > 5d0 .and. dabs(lEx%Flux0%H) < MaxHFlux) then
        if (lEx%Flux0%LE  <= 3d0 .or. dabs(lEx%Flux0%LE) > MaxLEFlux) &
            BinCospForStable%of(h2o) = error

        if (dabs(lEx%Flux0%co2)  <= 2d0 .or. dabs(lEx%Flux0%co2)  > MaxCO2Flux)  &
            BinCospForStable%of(co2)  = error

        if (dabs(lEx%Flux0%ch4)  <= 1d-3 .or. dabs(lEx%Flux0%ch4)  > MaxCH4Flux)  &
            BinCospForStable%of(ch4)  = error

        if (dabs(lEx%Flux0%gas4)  > MaxGAS4Flux) &
            BinCospForStable%of(gas4) = error

        if (lEx%Flux0%LE  <= 3d0 .and. &
            dabs(lEx%Flux0%co2)  <= 2d0) then
            BinCospForStable = ErrSpec
            skip_cospectra = .true.
        end if
    else
        BinCospForStable = ErrSpec
        skip_cospectra = .true.
    end if

    !> Further filter on u* for cospectra
    if (lEx%ustar < 0.05d0) then
        BinCospForStable = ErrSpec
        skip_cospectra = .true.
    end if

    !> For time sorted cospectra use milder filtering
    BinCosp = BinCospForStable
end subroutine FilterSpectraByThresholds
