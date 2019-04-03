!***************************************************************************
! spectra_sorting_and_averaging.f90
! ---------------------------------
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


    call char2int(lEx%end_date(6:7), month, 2)
    do gas = co2, gas4
        sort = 0
        if (gas /= h2o) then
            sort = FCCsetup%SA%class(gas, month)
        else
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
        end if

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
