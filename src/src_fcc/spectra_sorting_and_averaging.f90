!***************************************************************************
! spectra_sorting_and_averaging.f90
! ---------------------------------
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
! \brief       Sort current spectra in the respective classes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SpectraSortingAndAveraging(lEx, BinSpec, nrow, nbins)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: nbins
    type(ExType), intent(in) :: lEx
    type(SpectraSetType), intent(in) :: BinSpec(nrow)
    !> local variables
    integer :: sort
    integer :: i
    integer :: gas
    integer :: month

    !> H2O spectrum
    sort = 0
    if (lEx%RH > 5d0 .and. lEx%RH < 15d0) then
        sort = RH10
    elseif (lEx%RH > 15d0 .and. lEx%RH < 25d0) then
        sort = RH20
    elseif (lEx%RH > 25d0 .and. lEx%RH < 35d0) then
        sort = RH30
    elseif (lEx%RH > 35d0 .and. lEx%RH < 45d0) then
        sort = RH40
    elseif (lEx%RH > 45d0 .and. lEx%RH < 55d0) then
        sort = RH50
    elseif (lEx%RH > 55d0 .and. lEx%RH < 65d0) then
        sort = RH60
    elseif (lEx%RH > 65d0 .and. lEx%RH < 75d0) then
        sort = RH70
    elseif (lEx%RH > 75d0 .and. lEx%RH < 85d0) then
        sort = RH80
    elseif (lEx%RH > 85d0 .and. lEx%RH < 95d0) then
        sort = RH90
    end if
    if (sort /= 0) then
        do i = 1, nbins
            if (BinSpec(i)%fnum /= 0 .and. BinSpec(i)%of(h2o) /= error .and. BinSpec(i)%of(ts) /= error) then
                MeanBinSpec(i, sort)%cnt(h2o)  = MeanBinSpec(i, sort)%cnt(h2o)  + 1
                MeanBinSpec(i, sort)%fnum(h2o) = MeanBinSpec(i, sort)%fnum(h2o) + BinSpec(i)%fnum
                MeanBinSpec(i, sort)%fn(h2o)   = MeanBinSpec(i, sort)%fn(h2o)   + BinSpec(i)%fn
                MeanBinSpec(i, sort)%of(h2o)   = MeanBinSpec(i, sort)%of(h2o)   + BinSpec(i)%of(h2o)
                MeanBinSpec(i, sort)%ts(h2o)   = MeanBinSpec(i, sort)%ts(h2o)   + BinSpec(i)%of(ts)
            end if
        end do
    end if

    !> Passive gases spectra
    !> select month class
    sort = 0
    call char2int(lEx%date(6:7), month, 2)

    do gas = co2, gas4
        if (gas == h2o) cycle
        sort = FCCsetup%SA%class(gas, month)
        if (sort /= 0) then
            do i = 1, nbins
                if (BinSpec(i)%fnum /= 0 .and. BinSpec(i)%of(gas) /= error .and. BinSpec(i)%of(ts) /= error) then
                    MeanBinSpec(i, sort)%cnt(gas)  = MeanBinSpec(i, sort)%cnt(gas)  + 1
                    MeanBinSpec(i, sort)%fnum(gas) = MeanBinSpec(i, sort)%fnum(gas) + BinSpec(i)%fnum
                    MeanBinSpec(i, sort)%fn(gas)   = MeanBinSpec(i, sort)%fn(gas)   + BinSpec(i)%fn
                    MeanBinSpec(i, sort)%of(gas)   = MeanBinSpec(i, sort)%of(gas)   + BinSpec(i)%of(gas)
                    MeanBinSpec(i, sort)%ts(gas)   = MeanBinSpec(i, sort)%ts(gas)   + BinSpec(i)%of(ts)
                end if
            end do
        end if
    end do
end subroutine SpectraSortingAndAveraging
