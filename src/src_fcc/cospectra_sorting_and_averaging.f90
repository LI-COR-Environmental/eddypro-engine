!***************************************************************************
! cospectra_sorting_and_averaging.f90
! -----------------------------------
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
subroutine CospectraSortingAndAveraging(BinCosp, nrow, time, nbins)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: nbins
    character(*), intent(in) :: time
    type(SpectraSetType), intent(in) :: BinCosp(nrow)
    !> Local variables
    integer :: sort
    integer :: i
    integer :: gas

    !> Sort and average all cospectra in 8 classes of 3 hours
    !> for a sort of daily cospectral course
    sort = 0
    select case (time(1:2))
        case('00','01','02')
        sort = 1
        case('03','04','05')
        sort = 2
        case('06','07','08')
        sort = 3
        case('09','10','11')
        sort = 4
        case('12','13','14')
        sort = 5
        case('15','16','17')
        sort = 6
        case('18','19','20')
        sort = 7
        case('21','22','23')
        sort = 8
    end select

    if (sort == 0) return
    do gas = w_ts, w_gas4
        do i = 1, nbins
            if (BinCosp(i)%fnum /= 0 .and. BinCosp(i)%of(gas) /= error) then
                MeanBinCosp(i, sort)%cnt(gas)  = MeanBinCosp(i, sort)%cnt(gas)  + 1
                MeanBinCosp(i, sort)%fnum(gas) = MeanBinCosp(i, sort)%fnum(gas) + BinCosp(i)%fnum
                MeanBinCosp(i, sort)%fn(gas)   = MeanBinCosp(i, sort)%fn(gas)   + BinCosp(i)%fn
                MeanBinCosp(i, sort)%of(gas)   = MeanBinCosp(i, sort)%of(gas)   + BinCosp(i)%of(gas)
            end if
        end do
    end do
end subroutine CospectraSortingAndAveraging

