!***************************************************************************
! biomet_retrieve_embedded_data.f90
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
! \brief       Finalize retrieval of biomet data from already created bSet \n
!              in case of embedded biomet files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometRetrieveEmbeddedData(proceed, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: proceed
    logical, intent(in) :: printout
    !> local variables
    integer :: i

    !> Retrieve embedded biomet data if they exist (the option was
    !> selected and data was successfully read with at least one
    !> valid biomet record)
    if (printout) write(*,'(a)') '  Retrieving biomet data..'

    !> Initialize biomet data to error
    if (allocated(bAggr)) bAggr = error
    if (allocated(bAggrFluxnet)) bAggrFluxnet = error
    if (allocated(bAggrEddyPro)) bAggrEddyPro = error

    if (proceed) then
        if (printout) write(LogInteger, '(i3)') nbRecs
        if (printout) write(*, '(a)') '   ' // trim(adjustl(LogInteger)) &
            // ' biomet records imported.'

        !> Aggregate biomet variables over the averaging interval
        call BiometAggregate(bSet, size(bSet, 1), size(bSet, 2), bAggr)

        !> Convert data to standard units
        call BiometStandardEddyProUnits()

        !> Aggregate biomet variables over the averaging interval
        call BiometAggregate(bSet, size(bSet, 1), size(bSet, 2), bAggrEddyPro)

        !> Convert aggregated values to FLUXNET units
        call BiometStandardFluxnetUnits()
    else
        if (printout) call ExceptionHandler(72)
    end if

    !> Associate values to variables, as selected by user
    !> The - 2 is to account for the DATE and TIME columns in the file, which
    !> are not included in bAggr. The 2 shall eventually be replaced by nbTimestamp
    !> as per read_biomet_meta_file.f90
    do i = bTa, bRg
        if (bSetup%sel(i) > 0) biomet%val(i) = bAggrEddyPro(bSetup%sel(i) - 2)
    end do

    !> Deallocate variables no longer used
    if (allocated(bSet)) deallocate(bSet)
    if (allocated(bTs)) deallocate(bTs)
    if (printout) write(*,'(a)') '  Done.'

end subroutine BiometRetrieveEmbeddedData



