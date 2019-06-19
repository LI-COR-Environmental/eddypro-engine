!***************************************************************************
! extract_column_from_essentials.f90
! ----------------------------------
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
! \brief       Extract data column(s) from essentials files, based on what is
!              requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExtractColumnFromEssentials(ExFilename, NumExRecords, column, array, nrow, ncol, NumActRecords)
    use m_common_global_var
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: NumExRecords
    integer, intent(out) :: NumActRecords
    character(*), intent(in) :: ExFilename
    character(*), intent(in) :: column
    real(kind = dbl), intent(out) :: array(nrow, ncol)
    !> Local variables
    integer :: open_status
    integer :: cnt
    logical :: ValidRecord
    logical :: EndOfFileReached
    type(ExType) :: lEx

    !> Open file, stop execution in case of problems and skip header
    open(uex, file = ExFilename, status = 'old', iostat = open_status)
    if (open_status /= 0) call ExceptionHandler(60)
    !> Skip header
    read(uex, *)

    !> Based on what was selected, perform appropriate part of the routine
    array = error
    cnt = 0
    select case (trim(adjustl(column)))
        case ('degraded_wT_covariances')
            do i = 1, NumExRecords
                call ReadExRecord('', uex, -1, lEx, ValidRecord, EndOfFileReached)
                if (EndOfFileReached) exit
                !> Skip implausible values from the dataset
                if (dabs(lEx%WS) > MaxWindIntensity .or. dabs(lEx%degT%cov) > MaxWTCov &
                    .or. dabs(lEx%degT%dcov(1)) > MaxWTCov .or. dabs(lEx%degT%dcov(Nt)) > MaxWTCov) cycle
                if (ValidRecord) then
                    cnt = cnt + 1
                    array(cnt, 1:Nt) = lEx%degT%dcov(1:Nt)
                    array(cnt, Nt + 1) = lEx%degT%cov
                    array(cnt, Nt + 2) = lEx%WS
                    array(cnt, Nt + 3) = lEx%zL
                end if
            end do
            NumActRecords = cnt
    end select
    close(uex)
end subroutine ExtractColumnFromEssentials
