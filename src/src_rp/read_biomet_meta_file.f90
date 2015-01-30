!***************************************************************************
! read_biomet_meta_file.f90
! -------------------------
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
! \brief       Reads biomet metadata file for interpreting paired data file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadBiometMetaFile(MetaFile, skip_file)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: MetaFile
    logical, intent(out) :: skip_file


    !> Initialize aux variables needed to read biomet file
    call biometSniffMetaFile(MetaFile, skip_file)
    if (skip_file .or. nbVars <= 0) return
    call biometInitEmbedded()

    !> parse ini file and store all numeric and character tags
    call ParseIniFile(MetaFile, '', BiometNTags, BiometCTags, &
        size(BiometNTags), size(BiometCTags), &
         BiometNTagFound, BiometCTagFound, skip_file)

    !> selects only tags needed in this software, and store
    !> them in relevant variables
    call WriteBiometMetaVariables(skip_file)

end subroutine ReadBiometMetaFile

!***************************************************************************
!
! \brief       Retrieves relevant variables from the tags found in the \n
!              metadata file expected in EddyPro .metadata format
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteBiometMetaVariables(skip_file)
    use m_rp_global_var
    implicit none
    !> In/out variables
    logical, intent(out) :: skip_file
    !> Local variables
    integer :: initn_col
    integer :: initc_col
    integer :: leapn_col
    integer :: leapc_col
    integer :: i
    integer :: cnt
    integer :: tsCnt
    integer :: ix
    integer :: nbTimestamp
    character(32) :: label


    !> File general features
    skip_file = .false.
    bFileMetadata%nhead     = nint(BiometNTags(1)%value)
    bFileMetadata%duration  = BiometNTags(2)%value
    bFileMetadata%time_step = nint(BiometNTags(3)%value)

    select case (trim(adjustl(BiometCTags(1)%value)))
        case ('tab')
            bFileMetadata%separator = char(9)
        case ('comma')
            bFileMetadata%separator = ','
        case ('semicolon')
            bFileMetadata%separator = ';'
        case ('space')
            bFileMetadata%separator = ' '
        case default
            bFileMetadata%separator = BiometCTags(1)%value(1:1)
    end select
    bFileMetadata%data_label = trim(adjustl(BiometCTags(2)%value))

    !> Determine number of timestamp variables
    leapc_col = 7
    initc_col = 3 - leapc_col
    nbTimestamp = 0
    do i = 1, nbVars
        if(BiometCTagFound(initc_col + i*leapc_col)) then
            label = trim(adjustl(BiometCTags(initc_col + i*leapc_col)%value))
            call uppercase(label)
            if (label == 'DATE' .or. label == 'TIME') then
                nbTimestamp = nbTimestamp + 1
                cycle
            end if
        else
            exit
        end if
    end do

    !> Reduce nbVars to the number of variables excluding timestamp
    nbVars = nbVars - nbTimestamp

    !> Control on number of biomet variables found
    if (nbVars <= 0) then
        skip_file = .true.
        return
    end if

    !> Allocate and initialize bVars
    if (allocated(bVars)) deallocate(bVars)
    allocate(bVars(nbVars))
    if (allocated(bAggr)) deallocate(bAggr)
    allocate(bAggr(nbVars))
    bVars = nullbVar

    !> Variables description
    leapn_col = 6
    leapc_col = 7
    initn_col = 4 - leapn_col
    initc_col = 3 - leapc_col
    cnt = 0
    tsCnt = 0
    bFileMetadata%tsPattern = ''
    do i = 1, nbVars + nbTimestamp
        ix = initc_col + i*leapc_col
        if(BiometCTagFound(ix)) then
            !> if
            label = trim(adjustl(BiometCTags(ix)%value))
            call uppercase(label)
            if (label == 'DATE') then
                tsCnt = tsCnt + 1
                bFileMetadata%tsCols(tsCnt) = i
                bFileMetadata%tsPattern = &
                    trim(adjustl(bFileMetadata%tsPattern)) &
                    // 'yyyy-mm-dd'
            else if (label == 'TIME') then
                tsCnt = tsCnt + 1
                bFileMetadata%tsCols(tsCnt) = i
                bFileMetadata%tsPattern = &
                    trim(adjustl(bFileMetadata%tsPattern)) &
                    // 'HH:MM'
            else
                cnt = cnt + 1
                bVars(cnt)%label    = trim(adjustl(BiometCTags(ix)%value))
                bVars(cnt)%id       = trim(adjustl(BiometCTags(ix + 1)%value))
                bVars(cnt)%instr    = trim(adjustl(BiometCTags(ix + 2)%value))
                bVars(cnt)%unit_in  = trim(adjustl(BiometCTags(ix + 3)%value))
            end if
        end if
    end do
    bFileMetadata%numTsCol = tsCnt

    !> If variable has no label, stick ID to label
    do i = 1, nbVars
        if (len_trim(bVars(i)%label) == 0) bVars(i)%label = bVars(i)%id
    end do

    !> Append suffix if variables have not
    call BiometAppendLocationSuffix()

    !> Fill variables information based on label and other available fields
    call BiometEnrichVarsDescription()
end subroutine WriteBiometMetaVariables
