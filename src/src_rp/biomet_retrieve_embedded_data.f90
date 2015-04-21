!***************************************************************************
! biomet_retrieve_embedded_data.f90
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

    if (proceed) then
        if (printout) write(LogInteger, '(i3)') nbRecs
        if (printout) write(*, '(a)') '   ' // trim(adjustl(LogInteger)) &
            // ' biomet records imported.'

        !> Convert data to standard units
        call BiometStandardEddyProUnits()

        !> Aggregate biomet variables over the averaging interval
        call BiometAggregate(bSet, size(bSet, 1), size(bSet, 2), bAggr)

        !> Convert aggregated values to FLUXNET units
        call BiometStandardFluxnetUnits()
    else
        if (printout) call ExceptionHandler(72)
    end if

    !> Associate values to variables, as selected by user
    do i = bTa, bRg
        if (bSetup%sel(i) > 0) biomet%val(i) = bAggr(bSetup%sel(i))
    end do

    !> Deallocate variables no longer used
    if (allocated(bSet)) deallocate(bSet)
    if (allocated(bTs)) deallocate(bTs)
    if (printout) write(*,'(a)') '  Done.'

end subroutine BiometRetrieveEmbeddedData



