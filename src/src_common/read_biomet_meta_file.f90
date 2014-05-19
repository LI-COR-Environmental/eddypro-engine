!***************************************************************************
! read_biomet_meta_file.f90
! -------------------------
! Copyright (C) 2011-2014, LI-COR Biosciences
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
subroutine ReadBiometMetaFile(MetaFile, IniFileNotFound)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: MetaFile
    logical, intent(out) :: IniFileNotFound


    !> parse ini file and store all numeric and character tags
    call ParseIniFile(MetaFile, '', BiometNTags, BiometCTags, size(BiometNTags), size(BiometCTags), &
         BiometNTagFound, BiometCTagFound, IniFileNotFound)

    !> selects only tags needed in this software, and store them in relevant variables
    call WriteBiometMetaVariables()

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
subroutine WriteBiometMetaVariables()
    use m_common_global_var
    implicit none
    !> in/out variables
    !> local variables
    integer :: initn_col
    integer :: initc_col
    integer :: leapn_col
    integer :: leapc_col
    integer :: i = 0


    !> File general features
    BiometMeta%nhead      = nint(BiometNTags(1001)%value)
    BiometMeta%duration   = BiometNTags(1002)%value
    BiometMeta%rate       = BiometNTags(1003)%value

    select case (BiometCTags(1001)%value(1:len_trim(BiometCTags(1001)%value)))
        case ('tab')
        BiometMeta%separator = char(9)
        case ('comma')
        BiometMeta%separator = ','
        case ('semicolon')
        BiometMeta%separator = ';'
        case ('space')
        BiometMeta%separator = ' '
        case default
        BiometMeta%separator = BiometCTags(1001)%value(1:1)
    end select
    BiometMeta%data_label = BiometCTags(1002)%value(1:len_trim(BiometCTags(1002)%value))

    !> Initialization of column information
    BiometMeta%BiometCol = NullBiometCol
    BiometMeta%BiometCol(:)%var = 'IGNORE'

    !> Variables description
    leapn_col = 6
    leapc_col = 7
    initn_col = 1 - leapn_col
    initc_col = 1 - leapc_col
    NumBiometCol = 0
    do i = 1, MaxNumBiometCol
        if(BiometCTagFound(initc_col + i*leapc_col)) then
            NumBiometCol = NumBiometCol + 1
            !> Textual fields
            !> Variable name
            BiometMeta%BiometCol(i)%var = BiometCTags(initc_col + i*leapc_col)%value &
                (1:len_trim(BiometCTags(initc_col + i*leapc_col)%value))
            !> Variable ID
            BiometMeta%BiometCol(i)%id = BiometCTags(initc_col + i*leapc_col + 1)%value &
                (1:len_trim(BiometCTags(initc_col + i*leapc_col + 1)%value))
            !> Variable Instrument
            BiometMeta%BiometCol(i)%instr = BiometCTags(initc_col + i*leapc_col + 2)%value &
                (1:len_trim(BiometCTags(initc_col + i*leapc_col + 2)%value))
            !> Variable unit in
            BiometMeta%BiometCol(i)%unit_in = BiometCTags(initc_col + i*leapc_col + 3)%value &
                (1:len_trim(BiometCTags(initc_col + i*leapc_col + 3)%value))
            !> Variable unit out
            BiometMeta%BiometCol(i)%unit_out = BiometCTags(initc_col + i*leapc_col + 4)%value &
                (1:len_trim(BiometCTags(initc_col + i*leapc_col + 4)%value))
            !> Numeric fields
            !> Variable name
            BiometMeta%BiometCol(i)%gain   = BiometNTags(initn_col + i*leapn_col)%value
            BiometMeta%BiometCol(i)%offset = BiometNTags(initn_col + i*leapn_col + 1)%value
        end if
    end do

    !> Upper-case everything
    do i = 1, NumBiometCol
        call uppercase(BiometMeta%BiometCol(i)%var(1:len_trim(BiometMeta%BiometCol(i)%var)))
        call uppercase(BiometMeta%BiometCol(i)%id(1:len_trim(BiometMeta%BiometCol(i)%id)))
        call uppercase(BiometMeta%BiometCol(i)%unit_in(1:len_trim(BiometMeta%BiometCol(i)%unit_in)))
        call uppercase(BiometMeta%BiometCol(i)%unit_out(1:len_trim(BiometMeta%BiometCol(i)%unit_out)))
        call uppercase(BiometMeta%BiometCol(i)%instr(1:len_trim(BiometMeta%BiometCol(i)%instr)))
    end do

    !> Caclulate how many date-related columns are present in the files
    numBiometDateCol = 0
    if (BiometMeta%BiometCol(1)%var == 'DATE' .or. BiometMeta%BiometCol(1)%var == 'TIME') numBiometDateCol = 1
    if (BiometMeta%BiometCol(1)%var == 'DATE' .and. BiometMeta%BiometCol(2)%var == 'TIME') numBiometDateCol = 2
    NumBiometVar = NumBiometCol - numBiometDateCol

    !> Calculate time step
    if (BiometMeta%rate > 0) then
        BiometMeta%step = 1d0 / BiometMeta%rate / 60d0
    else
        BiometMeta%step = error
    end if
end subroutine WriteBiometMetaVariables
