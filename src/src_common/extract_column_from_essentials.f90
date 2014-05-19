!***************************************************************************
! extract_column_from_essentials.f90
! ----------------------------------
! Copyright (C) 2013-2014, LI-COR Biosciences
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
    if (open_status /= 0) call ErrorHandle(2, 0, 1)
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
