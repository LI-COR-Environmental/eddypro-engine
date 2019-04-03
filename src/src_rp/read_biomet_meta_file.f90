!***************************************************************************
! read_biomet_meta_file.f90
! -------------------------
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
    logical, external :: BiometValidateVar


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
    if (allocated(bAggrFluxnet)) deallocate(bAggrFluxnet)
    allocate(bAggrFluxnet(nbVars))
    if (allocated(bAggrEddyPro)) deallocate(bAggrEddyPro)
    allocate(bAggrEddyPro(nbVars))
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
                call uppercase(bVars(cnt)%unit_in)

                !> Check validity of biomet variable label
                !if (.not. BiometValidateVar(bVars(cnt))) then
                !    call ExceptionHandler(73)
                !    skip_file = .true.
                !    return
                !end if

                if (len_trim(bVars(cnt)%label) == 0) bVars(cnt)%label = bVars(cnt)%id
                if (len_trim(bVars(cnt)%label) == 0) bVars(cnt)%label = 'UNNAMED'

                !> Retrieve variable base name
                call biometBaseNameAndPositionalQualifierFromLabel(bVars(cnt)%label, &
                    bVars(cnt)%base_name, bVars(cnt)%pq_string)

           end if
        end if
    end do

    !> Validate timestamp pattern
    call tsValidateTemplate(bFileMetadata%tsPattern)

    !> Set ISO timestamp or not
    if (index(bFileMetadata%tsPattern, 'ddd') /= 0) then
        bFileMetadata%tsIso = .false.
    else
        bFileMetadata%tsIso = .true.
    end if
    bFileMetadata%numTsCol = tsCnt

    !> Fill variables information based on label and other available fields
    call BiometEnrichVarsDescription()

    ! !> Append suffix if variables have not
    ! if (EddyProProj%fluxnet_standardize_biomet) &
    !     call BiometAppendDefaultPositionalQualifier()

end subroutine WriteBiometMetaVariables
