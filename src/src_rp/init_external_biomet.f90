!***************************************************************************
! init_external_biomet.f90
! ------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Read 1 biomet file and figure out time step with \n
!              a call to BiometTimeStep()
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine InitExternalBiomet()
    use m_rp_global_var
    implicit none
    integer :: i
    integer :: ii
    integer :: j
    integer :: jj
    integer :: nbiomet
    integer :: nfile
    integer :: nfl
    integer :: open_status
    integer :: read_status
    integer :: sepa
    integer :: var_num
    character(32) :: text_vars(1000)
    character(1024) :: datastring
    character(64) :: tstamp_string
    type (FilelistType) :: BiometFileList(1000)


    call log_msg(' inf=initializing external biomet usage.')
    write(*, '(a)', advance = 'no') ' Initializing external biomet usage..'

    if (EddyProProj%biomet_data == 'ext_file') then
        LogString = ' biomet_file_in=' // AuxFile%biomet(1:len_trim(AuxFile%biomet))
        call DoubleCharInString(LogString, slash)
        call log_msg(LogString)
        nfile = 1
        BiometFileList(1)%path = AuxFile%biomet
    elseif (EddyProProj%biomet_data == 'ext_dir') then
        call FileListByExt(Dir%biomet, trim(adjustl(EddyProProj%biomet_tail)), .false., 'none', .false., &
            .false., EddyProProj%biomet_recurse, BiometFileList, size(BiometFileList), .false., ' ')
    end if

    BiometUnits = NullBiometUnits
    ProfileUnits = NullProfileUnits
    i = 0
    file_loop: do nfl = 1, nfile
        !> Open biomet measurement file(s) and read data, selecting those in plausible ranges
        open(udf, file = BiometFileList(nfl)%path, status = 'old', iostat = open_status)
        write(LogLogical, '(L1)') open_status
        LogString = ' open_error=' //Loglogical
        call log_msg(LogString)

        Biomet(1:MaxNumBiometRow) = ErrBiomet
        if (open_status == 0) then
            !> Interpret file header or "header string"
            call ReadBiometHeader(udf)
            NumBiometVar = NumSlowVar
            !> Preliminarly define date prototype
            BiometSetup%tstamp_prototype = ''
            do ii = 1, NumSlowVar
                if (index(BiometOrd(ii)%var, 'TIMESTAMP') /= 0) then
                    BiometSetup%tstamp_prototype = BiometSetup%tstamp_prototype(1:len_trim(BiometSetup%tstamp_prototype)) &
                        // BiometOrd(ii)%units(1:len_trim(BiometOrd(ii)%units))
                    NumBiometVar = NumBiometVar - 1
                end if
            end do

            !> Start loop on datalines
            rec_loop: do

                !> Exit instruction
                if (i > MaxNumBiometRow - 1) exit file_loop

                tstamp_string = ''
                datastring = ''
                var_num = 0
                read(udf, '(a)', iostat = read_status) datastring
                if (read_status /= 0) then
                    close(udf)
                    exit rec_loop
                end if

                i = i + 1
                text_vars = 'none'
                do
                    sepa = index(datastring, BiometSetup%separator)
                    if (sepa == 0) sepa = len_trim(datastring) + 1
                    if (len_trim(datastring) == 0) exit
                    var_num = var_num + 1
                    text_vars(var_num) = datastring(1:sepa - 1)
                    datastring = datastring(sepa + 1: len_trim(datastring))
                end do
                if (var_num /= NumSlowVar)then
                    i = i - 1
                    cycle rec_loop
                end if

                var_loop: do j = 1, var_num
                    if (len_trim(text_vars(j)) == 0) cycle
                    !> Import timestamp
                    do jj = TIMESTAMP_1, TIMESTAMP_7
                        if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                            if (BiometOrd(j)%units == 'ddd' .and. len_trim(text_vars(j)) == 2) &
                                text_vars(j) = '0' // text_vars(j)(1:len_trim(text_vars(j)))
                            if (BiometOrd(j)%units == 'ddd' .and. len_trim(text_vars(j)) == 1) &
                                text_vars(j) = '00' // text_vars(j)(1:len_trim(text_vars(j)))
                            if (BiometOrd(j)%units == 'dd' .and. len_trim(text_vars(j)) == 1) &
                                text_vars(j) = '0' // text_vars(j)(1:len_trim(text_vars(j)))
                            if (BiometOrd(j)%units == 'HH' .and. len_trim(text_vars(j)) == 1) &
                                text_vars(j) = '0' // text_vars(j)(1:len_trim(text_vars(j)))
                            if (BiometOrd(j)%units == 'MM' .and. len_trim(text_vars(j)) == 1) &
                                text_vars(j) = '0' // text_vars(j)(1:len_trim(text_vars(j)))
                            !> Special case of Excel messing up timestamp in american style mm/dd/yyyy
                            if (BiometOrd(j)%units == 'mm/dd/yyyy') then
                                if (text_vars(j)(2:2) == '/') &
                                    text_vars(j) = '0' // text_vars(j)(1:len_trim(text_vars(j)))
                                if (text_vars(j)(5:5) == '/') &
                                    text_vars(j) = text_vars(j)(1:3) // '0' // text_vars(j)(4:len_trim(text_vars(j)))
                            end if
                            !> Special case of Excel messing up timestamp in american style HH:MM:SS
                            if (BiometOrd(j)%units == 'HH:MM:SS') then
                                if (text_vars(j)(2:2) == ':') &
                                    text_vars(j) = '0' // text_vars(j)(1:len_trim(text_vars(j)))
                            end if
                            tstamp_string = tstamp_string(1:len_trim(tstamp_string)) &
                                // text_vars(j)(1:len_trim(text_vars(j)))
                            cycle var_loop
                        end if
                    end do

                    !> Read Biomet units
                    if (i == 1) then
                        do jj = TA_1_1_1, TA_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Ta(jj - TA_1_1_1 + 1) = BiometOrd(j)%units
                                if (BiometSetup%Ta == j) BiometSetup%Ta = jj - TA_1_1_1 + 1
                                cycle var_loop
                            end if
                        end do
                        do jj = PA_1_1_1, PA_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Pa(jj - PA_1_1_1 + 1) = BiometOrd(j)%units
                                if (BiometSetup%Pa == j) BiometSetup%Pa = jj - PA_1_1_1 + 1
                                cycle var_loop
                            end if
                        end do
                        do jj = RH_1_1_1, RH_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%RH(jj - RH_1_1_1 + 1) = BiometOrd(j)%units
                                if (BiometSetup%RH == j) BiometSetup%RH = jj - RH_1_1_1 + 1
                                cycle var_loop
                            end if
                        end do
                        do jj = RG_1_1_1, RG_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Rg(jj - RG_1_1_1 + 1) = BiometOrd(j)%units
                                if (BiometSetup%Rg == j) BiometSetup%Rg = jj - RG_1_1_1 + 1
                                cycle var_loop
                            end if
                        end do
                        do jj = RN_1_1_1, RN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Rn(jj - RN_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = RD_1_1_1, RD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Rd(jj - RD_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = RR_1_1_1, RR_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Rr(jj - RR_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = R_UVA_1_1_1, R_UVA_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Ruva(jj - R_UVA_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = R_UVB_1_1_1, R_UVB_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Ruvb(jj - R_UVB_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = LWIN_1_1_1, LWIN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%LWin(jj - LWIN_1_1_1 + 1) = BiometOrd(j)%units
                                if (BiometSetup%LWin == j) BiometSetup%LWin = jj - LWIN_1_1_1 + 1
                                cycle var_loop
                            end if
                        end do
                        do jj = LWOUT_1_1_1, LWOUT_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%LWout(jj - LWOUT_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = TC_1_1_1, TC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Tc(jj - TC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = TBOLE_1_1_1, TBOLE_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Tbole(jj - TBOLE_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = P_1_1_1, P_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%P(jj - P_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = P_RAIN_1_1_1, P_RAIN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Prain(jj - P_RAIN_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = P_SNOW_1_1_1, P_SNOW_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Psnow(jj - P_SNOW_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SNOWD_1_1_1, SNOWD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SNOWD(jj - SNOWD_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = PPFD_1_1_1, PPFD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%PPFD(jj - PPFD_1_1_1 + 1) = BiometOrd(j)%units
                                if (BiometSetup%PPFD == j) BiometSetup%PPFD = jj - PPFD_1_1_1 + 1
                                cycle var_loop
                            end if
                        end do
                        do jj = PPFDd_1_1_1, PPFDd_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%PPFDd(jj - PPFDd_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = PPFDr_1_1_1, PPFDr_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%PPFDr(jj - PPFDr_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = PPFDbc_1_1_1, PPFDbc_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                 BiometUnits%PPFDbc(jj - PPFDbc_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = APAR_1_1_1, APAR_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%APAR(jj - APAR_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = ALB_1_1_1, ALB_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%ALB(jj - ALB_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = PRI_1_1_1, PRI_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%PRI(jj - PRI_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = LAI_1_1_1, LAI_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%LAI(jj - LAI_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SWIN_1_1_1, SWIN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SWin(jj - SWIN_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SWOUT_1_1_1, SWOUT_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SWout(jj - SWOUT_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SWBC_1_1_1, SWBC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SWbc(jj - SWBC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SWDIF_1_1_1, SWDIF_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SWdif(jj - SWDIF_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = TBC_1_1_1, TBC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%Tbc(jj - TBC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = WS_1_1_1, WS_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%WS(jj - WS_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = MWS_1_1_1, MWS_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%MWS(jj - MWS_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = WD_1_1_1, WD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%WD(jj - WD_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SAPFLOW_1_1_1, SAPFLOW_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SAPFLOW(jj - SAPFLOW_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = STEMFLOW_1_1_1, STEMFLOW_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = TR_1_1_1, TR_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%TR(jj - TR_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SWC_1_1_1, SWC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SWC(jj - SWC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = SHF_1_1_1, SHF_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%SHF(jj - SHF_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                        do jj = TS_1_1_1, TS_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                BiometUnits%TS(jj - TS_1_1_1 + 1) = BiometOrd(j)%units
                                cycle var_loop
                            end if
                        end do
                    end if

                end do var_loop

                call BiometDateTime(trim(adjustl(BiometSetup%tstamp_prototype)), &
                    tstamp_string, Biomet(i)%date, Biomet(i)%time)

            end do rec_loop
            !exit file_loop
        else
            call log_msg( ' err=error while reading biomet file. file skipped.')
            call ErrorHandle(0, 0, 2)
        end if
    end do file_loop
    nbiomet = i
    close(udf)

    !> Assess Timestep
    call BiometTimeStep(nbiomet)
    write(*, '(a)') ' done'
end subroutine InitExternalBiomet

!***************************************************************************
!
! \brief       Reads and interprets file header, searching for known variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine ReadBiometHeader(unt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    !> local variables
    character(2048) :: varstring
    character(2048) :: unitstring
    character(32) :: HeaderVars(NumStdSlow + NumStdProfile + NumStdCustom)
    character(32) :: HeaderUnits(NumStdSlow + NumStdProfile + NumStdCustom)
    integer :: read_status
    integer :: var_sep
    integer :: unit_sep
    integer :: cnt
    integer :: i
    integer :: j

    !> Read strings containing variables and units
    if (BiometSetup%use_header) then
        !> Import varstring from Biomet measurement file
        do
            read(unt, '(a)', iostat = read_status) varstring
            if (read_status /= 0) then
                varstring  = 'none'
                exit
            end if
            if (index(varstring, '_1_1_1') /= 0) exit
        end do
        !> Import unitstring from Biomet measurement file
        do
            read(unt, '(a)', iostat = read_status) unitstring
            if (read_status /= 0) then
                unitstring = 'none'
                exit
            end if
            if (index(unitstring, 'yy') /= 0 .and. index(unitstring, 'dd') /= 0 .and. index(unitstring, 'HH') /= 0) exit
        end do
    else
        !> Take string from ".eddypro" file
!        varstring = BiometSetup%var_string
!        unitstring = BiometSetup%unit_string
    end if

    !> Retrieve variables and units from read string
    cnt = 0
    do
        var_sep = index(varstring, BiometSetup%separator)
        unit_sep = index(unitstring, BiometSetup%separator)
        if (var_sep == 0) var_sep = len_trim(varstring) + 1
        if (unit_sep == 0) unit_sep = len_trim(unitstring) + 1
        if (len_trim(varstring) == 0) exit
        cnt = cnt + 1
        HeaderVars(cnt) = varstring(1:var_sep - 1)
        HeaderUnits(cnt) = unitstring(1:unit_sep - 1)
        varstring = varstring(var_sep + 1: len_trim(varstring))
        unitstring = unitstring(unit_sep + 1: len_trim(unitstring))
        !> Support for Sutron files
        if (HeaderVars(cnt) == 'Date') HeaderVars(cnt) = 'TIMESTAMP_1'
        if (HeaderVars(cnt) == 'Time') HeaderVars(cnt) = 'TIMESTAMP_2'
    end do
    !> Adjust all Header labels: capitalize and schrinks
    do i = 1, cnt
        call uppercase(HeaderVars(i))
        call SchrinkString(HeaderVars(i))
        !> But not for timestamp labels (capitalization meaningful)
        if(index(HeaderVars(i), 'TIMESTAMP') == 0) then
            call uppercase(HeaderUnits(i))
            call SchrinkString(HeaderUnits(i))
        end if
    end do

    !> Retrieve available variables
    n_cstm_biomet = 0
    CstmOrd(1:300)    = NullOrd
    BiometOrd(1:300)    = NullOrd
    ProfileOrd(1:300) = NullOrd
    ol: do i = 1, cnt
        !> Biomet variables
        do j = 1, NumStdSlow
            if(HeaderVars(i)(1:len_trim(HeaderVars(i))) == &
                StdSlow(j)(1:len_trim(StdSlow(j)))) then
                BiometOrd(i)%var = StdSlow(j)(1:len_trim(StdSlow(j)))
                !> Retrieve units
                BiometOrd(i)%units =  HeaderUnits(i)(1:len_trim(HeaderUnits(i)))
                cycle ol
            end if
        end do
        !> Profile variables
        do j = 1, NumStdProfile
            if(StdProfile(j)(1:len_trim(StdProfile(j))) == &
                HeaderVars(i)(1:len_trim(HeaderVars(i)))) then
                ProfileOrd(i)%var = StdProfile(j)(1:len_trim(StdProfile(j)))
                !> Retrieve units
                ProfileOrd(i)%units =  HeaderUnits(i)(1:len_trim(HeaderUnits(i)))
                cycle ol
            end if
        end do
        !> Custom variables
        n_cstm_biomet = n_cstm_biomet + 1
        CstmOrd(n_cstm_biomet)%var = HeaderVars(i)(1:len_trim(HeaderVars(i)))
        CstmOrd(n_cstm_biomet)%units =  HeaderUnits(i)(1:len_trim(HeaderUnits(i)))
    end do ol
    NumSlowVar = cnt

    !> Now skip remaining header lines, if any
    if (BiometSetup%use_header) then
        if (BiometSetup%head_lines > 2) then
            do i = 3, BiometSetup%head_lines
                read(unt, *)
            end do
        end if
    else
        do i = 1, BiometSetup%head_lines
            read(unt, *)
        end do
    end if
end subroutine ReadBiometHeader
!***************************************************************************
!
! \brief       Retrieve time step intrinsic in Biomet file (in minutes)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometTimeStep(N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: i
    integer :: j
    type(DateType) :: curr_date
    type(DateType) :: prev_date
    integer :: step
    integer :: VarStep(1000)
    integer :: maxi
    integer :: max

    VarStep = 0
    do i = 2, N
        call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, curr_date)
        call DateTimeToDateType(Biomet(i-1)%date, Biomet(i-1)%time, prev_date)
        step = nint(timelag(curr_date, prev_date) * 24d0 * 60d0) !< in minutes
        do j = 1, 1000
            if (step == j) then
                VarStep(j) = VarStep(j) + 1
                exit
            end if
        end do
    end do
    maxi = 0
    max = 0
    do i = 1, 1000
        if (VarStep(i) > max) then
            max = VarStep(i)
            maxi = i
        end if
    end do
    BiometSetup%tstep = maxi !< in minutes
end subroutine BiometTimeStep
