!***************************************************************************
! retrieve_ex_var_by_timestamp.f90
! --------------------------------
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
subroutine RetrieveExVarsByTimestamp(unt, Timestamp, lEx, endReached, skip)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    type(DateType), intent(in) :: Timestamp
    logical, intent(out) :: endReached
    logical, intent(out) :: skip
    type(ExType), intent(out) :: lEx
    !> Local variables
    type(DateType) :: ExTimestamp
    logical :: EndOfFileReached
    logical :: ValidRecord

    skip = .false.
    endReached = .false.
    do
        call ReadExRecord('', unt, -1, lEx, ValidRecord, EndOfFileReached)

        !> If end of files was reached, exit routine with error flag
        if (EndOfFileReached) then
            endReached = .true.
            skip = .true.
            return
        end if

        !> If timestamp matches, exit routine (with error flag if the case)
        call DateTimeToDateType(lEx%end_date, lEx%end_time, ExTimestamp)
        if (ExTimestamp == Timestamp) then
            if (.not. ValidRecord) skip = .true.
            return
        end if

        !> If timestamp exceeds the one looked for, backwards essentials unit
        !> and exit with error code
        if (ExTimestamp >= Timestamp) then
            backspace(unt)
            skip = .true.
            return
        end if
    end do
end subroutine RetrieveExVarsByTimestamp
