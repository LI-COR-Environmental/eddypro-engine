!***************************************************************************
! read_ext_biomet_files.f90
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
! \brief       Interpret and read a generic "Biomet variables" data file, \n
!              meant to containg meteo, ecological and profile measurements. \n
!              Import all recognized variables contained in it.
! \author      Gerardo Fratini
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine ReadExtBiometFiles(BiometDataExist, last_nfl, last_nrec, CurrentTimestamp, nbiomet)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type (DateType), intent(in) :: CurrentTimestamp
    integer, intent(out) :: nbiomet
    logical, intent(out) :: BiometDataExist
    integer, intent(inout) :: last_nfl
    integer, intent(inout) :: last_nrec
    !> local variables
    integer :: i
    integer :: j
    integer :: ii
    integer :: jj
    integer :: iii
!    integer :: kkk
!    integer :: lll
    integer :: nfl
    integer :: nrec
    integer :: open_status
    integer :: read_status
    integer :: var_num
    integer :: sepa
    integer :: nfile
    integer :: ncstm
    character(32) :: text_vars(1000)
    character(1024) :: datastring
    character(64) :: tstamp_string
    type (FilelistType) :: BiometFileList(1000)
    type (DateType) :: tol
    type (DateType) :: win
    type (DateType) :: biomet_date
    logical :: BiometPeriodHooked
    logical :: skip_init


    call log_msg(' inf=importing biomet data from external file(s).')

    !> Initializations
    skip_init = .true.

    !> Define time periods for which biomet data are needed
    win = datetype(0, 0, 0, 0, RPsetup%avrg_len)
    tol = datetype(0, 0, 0, 0, max(BiometSetup%tstep, nint(BiometMeta%step)) / 2)

    if (EddyProProj%biomet_data == 'ext_file') then
        nfile = 1
        BiometFileList(1)%path = AuxFile%biomet
    elseif (EddyProProj%biomet_data == 'ext_dir') then
        call FileListByExt(Dir%biomet, trim(adjustl(EddyProProj%biomet_tail)), .false., 'none', .false., &
            .false., EddyProProj%biomet_recurse, BiometFileList, size(BiometFileList), .false., ' ')
    end if

    !> Initialization
    Biomet(1:MaxNumBiometRow) = ErrBiomet
    Profile(1:MaxNumBiometRow) = ErrProfile

    i = 0
    nbiomet = 0
    BiometPeriodHooked = .false.
    file_loop: do nfl = last_nfl, nfile
        nrec = 0
        !> Open biomet measurement file(s) and read data, selecting those in plausible ranges
        open(udf, file = BiometFileList(nfl)%path, status = 'old', iostat = open_status)
        write(LogLogical, '(L1)') open_status
        LogString = ' open_error=' //Loglogical
        call log_msg(LogString)

        if (open_status == 0) then
            BiometDataExist = .true.

            if(nfl == last_nfl) write(*, '(a)') '  Searching biomet data in file: '
            write(*, '(a)') '   ' // BiometFileList(nfl)%path(1:len_trim(BiometFileList(nfl)%path))

            !> Skip header
            if (BiometSetup%head_lines > 0) then
                do ii = 1, BiometSetup%head_lines
                    read(udf, *)
                end do
            end if

            !> Start loop on data lines
            rec_loop: do
                if (skip_init .and. last_nrec >= 1) then
                    skip_init = .false.
                    nrec = last_nrec
                    do iii = 1, last_nrec
                        read(udf,*)
                    end do
                end if
                i = i + 1
                tstamp_string = ''
                datastring = ''
                var_num = 0
                read(udf, '(a)', iostat = read_status) datastring
                if (read_status < 0) then
                    close(udf)
                    cycle file_loop
                end if
                if (read_status > 0) then
                    i = i - 1
                    cycle rec_loop
                end if

                nrec = nrec + 1
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

                ncstm = 0
                ol: do j = 1, var_num
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
                            cycle ol
                        end if
                    end do

                    !> Retrieve timestamp from tstamp_string
                    call BiometDateTime(BiometSetup%tstamp_prototype(1:len_trim(BiometSetup%tstamp_prototype)), &
                        tstamp_string, Biomet(i)%date, Biomet(i)%time)
                    call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, biomet_date)

                    !> Check if biomet_date is strictly within the averaging period, otherwise cycles
                    if (biomet_date <= CurrentTimestamp + tol .and. biomet_date >= CurrentTimestamp - win + tol) then
                        BiometPeriodHooked = .true.
                        last_nfl = nfl

                        !> Import Biomet measurements
                        do jj = TA_1_1_1, TA_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Ta(jj - TA_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Ta(jj - TA_1_1_1 + 1) = BiometOrd(j)%units
!                                if (i == 1 .and. BiometSetup%Ta == j) BiometSetup%Ta = jj - TA_1_1_1 + 1
                                cycle ol
                            end if
                        end do
                        do jj = PA_1_1_1, PA_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Pa(jj - PA_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Pa(jj - PA_1_1_1 + 1) = BiometOrd(j)%units
!                                if (i == 1 .and. BiometSetup%Pa == j) BiometSetup%Pa = jj - PA_1_1_1 + 1
                                cycle ol
                            end if
                        end do
                        do jj = RH_1_1_1, RH_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%RH(jj - RH_1_1_1 + 1)
!                                if (i == 1) BiometUnits%RH(jj - RH_1_1_1 + 1) = BiometOrd(j)%units
!                                if (i == 1 .and. BiometSetup%RH == j) BiometSetup%RH = jj - RH_1_1_1 + 1
                                cycle ol
                            end if
                        end do
                        do jj = RG_1_1_1, RG_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Rg(jj - RG_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Rg(jj - RG_1_1_1 + 1) = BiometOrd(j)%units
!                                if (i == 1 .and. BiometSetup%Rg == j) BiometSetup%Rg = jj - RG_1_1_1 + 1
                                cycle ol
                            end if
                        end do
                        do jj = RN_1_1_1, RN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Rn(jj - RN_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Rn(jj - RN_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = RD_1_1_1, RD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Rd(jj - RD_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Rd(jj - RD_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = RR_1_1_1, RR_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%RR(jj - RR_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Rr(jj - RR_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = R_UVA_1_1_1, R_UVA_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Ruva(jj - R_UVA_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Ruva(jj - R_UVA_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = R_UVB_1_1_1, R_UVB_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Ruvb(jj - R_UVB_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Ruvb(jj - R_UVB_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = LWIN_1_1_1, LWIN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%LWin(jj - LWIN_1_1_1 + 1)
                                BiometUnits%LWin(jj - LWIN_1_1_1 + 1) = BiometOrd(j)%units
 !                               if (i == 1 .and. BiometSetup%LWin == j) BiometSetup%LWin = jj - LWIN_1_1_1 + 1
                                cycle ol
                            end if
                        end do
                        do jj = LWOUT_1_1_1, LWOUT_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%LWout(jj - LWOUT_1_1_1 + 1)
!                                if (i == 1) BiometUnits%LWout(jj - LWOUT_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = TC_1_1_1, TC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Tc(jj - TC_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Tc(jj - TC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = TBOLE_1_1_1, TBOLE_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Tbole(jj - TBOLE_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Tbole(jj - TBOLE_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = P_1_1_1, P_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%P(jj - P_1_1_1 + 1)
!                                if (i == 1) BiometUnits%P(jj - P_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = P_RAIN_1_1_1, P_RAIN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Prain(jj - P_RAIN_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Prain(jj - P_RAIN_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = P_SNOW_1_1_1, P_SNOW_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Psnow(jj - P_SNOW_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Psnow(jj - P_SNOW_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SNOWD_1_1_1, SNOWD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SNOWD(jj - SNOWD_1_1_1 + 1)
!                                if (i == 1) BiometUnits%SNOWD(jj - SNOWD_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = PPFD_1_1_1, PPFD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%PPFD(jj - PPFD_1_1_1 + 1)
!                                if (i = 1) BiometUnits%PPFD(jj - PPFD_1_1_1 + 1) = BiometOrd(j)%units
!                                if (i == 1 .and. BiometSetup%PPFD == j) BiometSetup%PPFD = jj - PPFD_1_1_1 + 1
                                cycle ol
                            end if
                        end do
                        do jj = PPFDd_1_1_1, PPFDd_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%PPFDd(jj - PPFDd_1_1_1 + 1)
!                                if (i == 1) BiometUnits%PPFDd(jj - PPFDd_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = PPFDr_1_1_1, PPFDr_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%PPFDr(jj - PPFDr_1_1_1 + 1)
!                                if (i == 1) BiometUnits%PPFDr(jj - PPFDr_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = PPFDbc_1_1_1, PPFDbc_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%PPFDbc(jj - PPFDbc_1_1_1 + 1)
!                                if (i == 1) BiometUnits%PPFDbc(jj - PPFDbc_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = APAR_1_1_1, APAR_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%APAR(jj - APAR_1_1_1 + 1)
!                                if (i == 1) BiometUnits%APAR(jj - APAR_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = ALB_1_1_1, ALB_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Alb(jj - ALB_1_1_1 + 1)
!                                if (i == 1) BiometUnits%ALB(jj - ALB_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = PRI_1_1_1, PRI_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%PRI(jj - PRI_1_1_1 + 1)
!                                if (i == 1) BiometUnits%PRI(jj - PRI_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = LAI_1_1_1, LAI_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%LAI(jj - LAI_1_1_1 + 1)
!                                if (i == 1) BiometUnits%LAI(jj - LAI_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SWIN_1_1_1, SWIN_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SWin(jj - SWIN_1_1_1 + 1)
                                if (i == 1) BiometUnits%SWin(jj - SWIN_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SWOUT_1_1_1, SWOUT_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SWout(jj - SWOUT_1_1_1 + 1)
!                                if (i == 1) BiometUnits%SWout(jj - SWOUT_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SWBC_1_1_1, SWBC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SWbc(jj - SWBC_1_1_1 + 1)
!                                if (i == 1) BiometUnits%SWbc(jj - SWBC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SWDIF_1_1_1, SWDIF_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SWdif(jj - SWDIF_1_1_1 + 1)
!                                if (i == 1) BiometUnits%SWdif(jj - SWDIF_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = TBC_1_1_1, TBC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%Tbc(jj - TBC_1_1_1 + 1)
!                                if (i == 1) BiometUnits%Tbc(jj - TBC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = WS_1_1_1, WS_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%WS(jj - WS_1_1_1 + 1)
!                                if (i == 1) BiometUnits%WS(jj - WS_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = MWS_1_1_1, MWS_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%MWS(jj - MWS_1_1_1 + 1)
!                                if (i == 1) BiometUnits%MWS(jj - MWS_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = WD_1_1_1, WD_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%WD(jj - WD_1_1_1 + 1)
!                                if (i == 1) BiometUnits%WD(jj - WD_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SAPFLOW_1_1_1, SAPFLOW_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SAPFLOW(jj - SAPFLOW_1_1_1 + 1)
 !                               if (i == 1) BiometUnits%SAPFLOW(jj - SAPFLOW_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = STEMFLOW_1_1_1, STEMFLOW_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%STEMFLOW(jj - STEMFLOW_1_1_1 + 1)
 !                               if (i == 1) BiometUnits%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = TR_1_1_1, TR_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%TR(jj - TR_1_1_1 + 1)
 !                               if (i == 1) BiometUnits%TR(jj - TR_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = SWC_1_1_1, SWC_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SWC(jj - SWC_1_1_1 + 1)
 !                               if (i == 1) BiometUnits%SWC(jj - SWC_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                       do jj = SHF_1_1_1, SHF_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%SHF(jj - SHF_1_1_1 + 1)
!                                if (i == 1) BiometUnits%SHF(jj - SHF_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do
                        do jj = TS_1_1_1, TS_10_1_1
                            if (index(BiometOrd(j)%var, StdSlow(jj)) == 1) then
                                read(text_vars(j), *) Biomet(i)%TS(jj - TS_1_1_1 + 1)
!                                 if (i == 1) BiometUnits%TS(jj - TS_1_1_1 + 1) = BiometOrd(j)%units
                                cycle ol
                            end if
                        end do

                        !> If current var was not found to be among recognized one, then it's a custom variable
                        ncstm = ncstm + 1
                        read(text_vars(j), *) CstmBiometSet(i, ncstm)

!                    TAKE THE DETECTION OF NUMBER AND UNITS OF PROFILE TO INIT_EXT_BIOMET!!!!!!!!!!!!!!!!!!!!!!!!
!                    !> Import profile measurements
!                    do jj = SWC_1_1_1, SWC_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - SWC_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - SWC_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - SWC_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%SWC(kkk, lll)
!                            if (i == 1) ProfileUnits%SWC(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = SHF_1_1_1, SHF_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - SHF_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - SHF_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - SHF_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%SHF(kkk, lll)
!                            if (i == 1) ProfileUnits%SHF(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = TS_1_1_1, TS_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - TS_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - TS_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - TS_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%ST(kkk, lll)
!                            if (i == 1) ProfileUnits%ST(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = CO2_1_1_1, CO2_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - CO2_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - CO2_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - CO2_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%CO2(kkk, lll)
!                            if (i == 1) ProfileUnits%CO2(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = H2O_1_1_1, H2O_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - H2O_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - H2O_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - H2O_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%H2O(kkk, lll)
!                            if (i == 1) ProfileUnits%H2O(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = CH4_1_1_1, CH4_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - CH4_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - CH4_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - CH4_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%CH4(kkk, lll)
!                            if (i == 1) ProfileUnits%CH4(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = GAS4_1_1_1, GAS4_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - GAS4_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - GAS4_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - GAS4_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%GAS4(kkk, lll)
!                            if (i == 1) ProfileUnits%GAS4(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
!                    do jj = T_1_1_1, T_5_7_1
!                        if (index(ProfileOrd(j)%var, StdProfile(jj)) == 1) then
!                            if (MOD(jj,MaxProfRep) /= 0) then
!                                kkk = MOD((jj - T_1_1_1 + 1), MaxProfRep)
!                                lll = (jj - T_1_1_1 + 1) / MaxProfRep + 1
!                            else
!                                kkk = MaxProfRep
!                                lll = (jj - T_1_1_1 + 1) / MaxProfRep
!                            end if
!                            read(text_vars(j), *) Profile(i)%T(kkk, lll)
!                            if (i == 1) ProfileUnits%T(kkk) = ProfileOrd(j)%units
!                            cycle ol
!                        end if
!                    end do
                    else
                        if (BiometPeriodHooked) then
                            last_nrec = nrec - 1
                            close(udf)
                            exit file_loop
                        else
                            i = i - 1
                            cycle rec_loop
                        end if
                    end if
                end do ol
            end do rec_loop
        else
            call log_msg( ' err=error while reading biomet file. file skipped.')
            call ErrorHandle(0, 0, 2)
        end if
    end do file_loop
    nbiomet = i - 1

    if (nbiomet < 1) then
        BiometDataExist = .false.
        write(*,'(a)') '  No valid biomet records found for this averaging period. Continuing without biomet data.'
    else
        write(LogInteger, '(i6)') nbiomet
        call SchrinkString(LogInteger)
        write(*,'(a)') '  ' // LogInteger(1:len_trim(LogInteger)) &
            // ' biomet record(s) imported correctly for this averaging period.'
    end if

    !> Adjust units as needed
    call BiometStandardUnits()

    !> Adjust timesteps of Biomet data if needed
    call AdjustBiometTimestamps(nbiomet)
end subroutine ReadExtBiometFiles

!***************************************************************************
!
! \brief       Convert input units into standard units
! \author      Gerardo Fratini
! \note        not part of EddyPro Express
!              Radiations (Rg, Rn, Rd, Rr, LWin, LWout, Ruva, Ruvb) are not expected to need unit conversion
!              Photons flux densities (PPFD, PPFDd, PPFDr, PPFDbc, APAR) are not expected to need unit conversion
!              Albedo (Alb) is not expected to need unit conversion
!              Photochemical reflectance index is not expected to need unit conversion
!              Soil water content (SWC) is not expected to need unit conversion
!              Soil heat flux (SHF) is not expected to need unit conversion
! \sa
! \bug
! \deprecated
! \test
!***************************************************************************
subroutine BiometStandardUnits()
    use m_rp_global_var
    implicit none
    !> local variables
!    integer :: kkk
!    integer :: lll
    integer :: jj

    !> Temperatures
    do jj = TA_1_1_1, TA_10_1_1
        select case(BiometUnits%Ta(jj - TA_1_1_1 + 1))
            case('C','°C')
                Biomet%Ta(jj - TA_1_1_1 + 1) = Biomet%Ta(jj - TA_1_1_1 + 1) + 273.16d0
            case('F','°F')
                Biomet%Ta(jj - TA_1_1_1 + 1) = (Biomet%Ta(jj - TA_1_1_1 + 1) - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case('CK')
                Biomet%Ta(jj - TA_1_1_1 + 1) = Biomet%Ta(jj - TA_1_1_1 + 1) * 1d-2
            case('CC','C°C')
                Biomet%Ta(jj - TA_1_1_1 + 1) = Biomet%Ta(jj - TA_1_1_1 + 1) * 1d-2 + 273.16d0
            case('CF','C°F')
                Biomet%Ta(jj - TA_1_1_1 + 1) = (Biomet%Ta(jj - TA_1_1_1 + 1) * 1d-2 - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case default
        end select
    end do
    do jj = TC_1_1_1, TC_10_1_1
        select case(BiometUnits%Tc(jj - TC_1_1_1 + 1))
            case('C','°C')
                Biomet%Tc(jj - TC_1_1_1 + 1) = Biomet%Tc(jj - TC_1_1_1 + 1) + 273.16d0
            case('F','°F')
                Biomet%Tc(jj - TC_1_1_1 + 1) = (Biomet%Tc(jj - TC_1_1_1 + 1) - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case('CK')
                Biomet%Tc(jj - TC_1_1_1 + 1) = Biomet%Tc(jj - TC_1_1_1 + 1) * 1d-2
            case('CC','C°C')
                Biomet%Tc(jj - TC_1_1_1 + 1) = Biomet%Tc(jj - TC_1_1_1 + 1) * 1d-2 + 273.16d0
            case('CF','C°F')
                Biomet%Tc(jj - TC_1_1_1 + 1) = (Biomet%Tc(jj - TC_1_1_1 + 1) * 1d-2 - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case default
        end select
    end do
    do jj = TS_1_1_1, TS_10_1_1
        select case(BiometUnits%Ts(jj - TS_1_1_1 + 1))
            case('C','°C')
                Biomet%Ts(jj - TS_1_1_1 + 1) = Biomet%Ts(jj - TS_1_1_1 + 1) + 273.16d0
            case('F','°F')
                Biomet%Ts(jj - TS_1_1_1 + 1) = (Biomet%Ts(jj - TS_1_1_1 + 1) - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case('CK')
                Biomet%Ts(jj - TS_1_1_1 + 1) = Biomet%Ts(jj - TS_1_1_1 + 1) * 1d-2
            case('CC','C°C')
                Biomet%Ts(jj - TS_1_1_1 + 1) = Biomet%Ts(jj - TS_1_1_1 + 1) * 1d-2 + 273.16d0
            case('CF','C°F')
                Biomet%Ts(jj - TS_1_1_1 + 1) = (Biomet%Ts(jj - TS_1_1_1 + 1) * 1d-2 - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case default
        end select
    end do
    do jj = TBC_1_1_1, TBC_10_1_1
        select case(BiometUnits%Tbc(jj - TBC_1_1_1 + 1))
            case('C','°C')
                Biomet%Tbc(jj - TBC_1_1_1 + 1) = Biomet%Tbc(jj - TBC_1_1_1 + 1) + 273.16d0
            case('F','°F')
                Biomet%Tbc(jj - TBC_1_1_1 + 1) = (Biomet%Tbc(jj - TBC_1_1_1 + 1) - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case('CK')
                Biomet%Tbc(jj - TBC_1_1_1 + 1) = Biomet%Tbc(jj - TBC_1_1_1 + 1) * 1d-2
            case('CC','C°C')
                Biomet%Tbc(jj - TBC_1_1_1 + 1) = Biomet%Tbc(jj - TBC_1_1_1 + 1) * 1d-2 + 273.16d0
            case('CF','C°F')
                Biomet%Tbc(jj - TBC_1_1_1 + 1) = (Biomet%Tbc(jj - TBC_1_1_1 + 1) * 1d-2 - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case default
        end select
    end do
    do jj = TBOLE_1_1_1, TBOLE_10_1_1
        select case(BiometUnits%Tc(jj - TBOLE_1_1_1 + 1))
            case('C','°C')
                Biomet%Tbole(jj - TBOLE_1_1_1 + 1) = Biomet%Tbole(jj - TBOLE_1_1_1 + 1) + 273.16d0
            case('F','°F')
                Biomet%Tbole(jj - TBOLE_1_1_1 + 1) = (Biomet%Tbole(jj - TBOLE_1_1_1 + 1) - 32d0) &
                    * 5d0 / 9d0 + 273.16d0
            case('CK')
                Biomet%Tbole(jj - TBOLE_1_1_1 + 1) = Biomet%Tbole(jj - TBOLE_1_1_1 + 1) * 1d-2
            case('CC','C°C')
                Biomet%Tbole(jj - TBOLE_1_1_1 + 1) = Biomet%Tbole(jj - TBOLE_1_1_1 + 1) * 1d-2 + 273.16d0
            case('CF','C°F')
                Biomet%Tbole(jj - TBOLE_1_1_1 + 1) = (Biomet%Tbole(jj - TBOLE_1_1_1 + 1) &
                    * 1d-2 - 32d0) * 5d0 / 9d0 + 273.16d0
            case default
        end select
    end do
    !> Pressures
    do jj = PA_1_1_1, PA_10_1_1
        select case(BiometUnits%Pa(jj - PA_1_1_1 + 1))
            case('HPA')
                Biomet%Pa(jj - PA_1_1_1 + 1) = Biomet%Pa(jj - PA_1_1_1 + 1) * 1d2
            case('KPA')
                Biomet%Pa(jj - PA_1_1_1 + 1) = Biomet%Pa(jj - PA_1_1_1 + 1) * 1d3
            case('MMHG', 'TORR')
                Biomet%Pa(jj - PA_1_1_1 + 1) = Biomet%Pa(jj - PA_1_1_1 + 1) * 133.32d0
            case('PSI')
                Biomet%Pa(jj - PA_1_1_1 + 1) = Biomet%Pa(jj - PA_1_1_1 + 1) * 6894.6d0
            case('BAR')
                Biomet%Pa(jj - PA_1_1_1 + 1) = Biomet%Pa(jj - PA_1_1_1 + 1) * 1d5
            case('ATM')
                Biomet%Pa(jj - PA_1_1_1 + 1) = Biomet%Pa(jj - PA_1_1_1 + 1) * 0.980665d5
            case default
                continue
        end select
    end do
    !> Relative humidities
    do jj = RH_1_1_1, RH_10_1_1
        select case(BiometUnits%RH(jj - RH_1_1_1 + 1))
            case('NUMBER','#','DIMENSIONLESS')
                Biomet%RH(jj - RH_1_1_1 + 1) = Biomet%RH(jj - RH_1_1_1 + 1) * 1d2
            case default
                continue
        end select
    end do
    !> Lengths
    do jj = P_1_1_1, P_10_1_1
        select case(BiometUnits%P(jj - P_1_1_1 + 1))
            case('NM')
                Biomet%P(jj - P_1_1_1 + 1) = Biomet%P(jj - P_1_1_1 + 1) * 1d-9
            case('UM')
                Biomet%P(jj - P_1_1_1 + 1) = Biomet%P(jj - P_1_1_1 + 1) * 1d-6
            case('MM')
                Biomet%P(jj - P_1_1_1 + 1) = Biomet%P(jj - P_1_1_1 + 1) * 1d-3
            case('CM')
                Biomet%P(jj - P_1_1_1 + 1) = Biomet%P(jj - P_1_1_1 + 1) * 1d-2
            case('KM')
                Biomet%P(jj - P_1_1_1 + 1) = Biomet%P(jj - P_1_1_1 + 1) * 1d3
            case default
                continue
        end select
    end do
    do jj = P_RAIN_1_1_1, P_RAIN_10_1_1
        select case(BiometUnits%Prain(jj - P_RAIN_1_1_1 + 1))
            case('NM')
                Biomet%Prain(jj - P_RAIN_1_1_1 + 1) = Biomet%Prain(jj - P_RAIN_1_1_1 + 1) * 1d-9
            case('UM')
                Biomet%Prain(jj - P_RAIN_1_1_1 + 1) = Biomet%Prain(jj - P_RAIN_1_1_1 + 1) * 1d-6
            case('MM')
                Biomet%Prain(jj - P_RAIN_1_1_1 + 1) = Biomet%Prain(jj - P_RAIN_1_1_1 + 1) * 1d-3
            case('CM')
                Biomet%Prain(jj - P_RAIN_1_1_1 + 1) = Biomet%Prain(jj - P_RAIN_1_1_1 + 1) * 1d-2
            case('KM')
                Biomet%Prain(jj - P_RAIN_1_1_1 + 1) = Biomet%Prain(jj - P_RAIN_1_1_1 + 1) * 1d3
            case default
                continue
        end select
    end do
    do jj = P_SNOW_1_1_1, P_SNOW_10_1_1
        select case(BiometUnits%Psnow(jj - P_SNOW_1_1_1 + 1))
            case('NM')
                Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) = Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) * 1d-9
            case('UM')
                Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) = Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) * 1d-6
            case('MM')
                Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) = Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) * 1d-3
            case('CM')
                Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) = Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) * 1d-2
            case('KM')
                Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) = Biomet%Psnow(jj - P_SNOW_1_1_1 + 1) * 1d3
            case default
                continue
        end select
    end do
    do jj = SNOWD_1_1_1, SNOWD_10_1_1
        select case(BiometUnits%SNOWD(jj - SNOWD_1_1_1 + 1))
            case('NM')
                Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) = Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) * 1d-9
            case('UM')
                Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) = Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) * 1d-6
            case('MM')
                Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) = Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) * 1d-3
            case('CM')
                Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) = Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) * 1d-2
            case('KM')
                Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) = Biomet%SNOWD(jj - SNOWD_1_1_1 + 1) * 1d3
            case default
                continue
        end select
    end do
    do jj = STEMFLOW_1_1_1, STEMFLOW_10_1_1
        select case(BiometUnits%STEMFLOW(jj - STEMFLOW_1_1_1 + 1))
            case('NM')
                Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) * 1d-9
            case('UM')
                Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) * 1d-6
            case('MM')
                Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) * 1d-3
            case('CM')
                Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) * 1d-2
            case('KM')
                Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) = Biomet%STEMFLOW(jj - STEMFLOW_1_1_1 + 1) * 1d3
            case default
                continue
        end select
    end do

    !> Velocities
    do jj = WS_1_1_1, WS_10_1_1
        select case(BiometUnits%WS(jj - WS_1_1_1 + 1))
            case('CM+1S-1','CM/S','CMS^-1','CMS-1')
                Biomet%WS(jj - WS_1_1_1 + 1) = Biomet%WS(jj - WS_1_1_1 + 1) * 1d-2
            case('MM+1S-1','MM/S','MMS^-1','MMS-1')
                Biomet%WS(jj - WS_1_1_1 + 1) = Biomet%WS(jj - WS_1_1_1 + 1) * 1d-3
            case default
                continue
        end select
    end do
    do jj = MWS_1_1_1, MWS_10_1_1
        select case(BiometUnits%MWS(jj - MWS_1_1_1 + 1))
            case('CM+1S-1','CM/S','CMS^-1','CMS-1')
                Biomet%MWS(jj - MWS_1_1_1 + 1) = Biomet%MWS(jj - MWS_1_1_1 + 1) * 1d-2
            case('MM+1S-1','MM/S','MMS^-1','MMS-1')
                Biomet%MWS(jj - MWS_1_1_1 + 1) = Biomet%MWS(jj - MWS_1_1_1 + 1) * 1d-3
            case default
                continue
        end select
    end do

    !> Profile temperatures
!    do jj = TS_1_1_1, TS_5_7_1
!        if (MOD(jj,MaxProfRep) /= 0) then
!            kkk = MOD((jj - TS_1_1_1 + 1), MaxProfRep)
!            lll = (jj - TS_1_1_1 + 1) / MaxProfRep + 1
!        else
!            kkk = MaxProfRep
!            lll = (jj - TS_1_1_1 + 1) / MaxProfRep
!        end if
!        select case(ProfileUnits%ST(kkk))
!            case('C','°C')
!                Profile%ST(kkk, lll) = Profile%ST(kkk, lll) + 273.16d0
!            case('F','°F')
!                Profile%ST(kkk, lll) = (Profile%ST(kkk, lll) &
!                    - 32d0) * 5d0 / 9d0 + 273.16d0
!            case('CK')
!                Profile%ST(kkk, lll) = Profile%ST(kkk, lll) * 1d-2
!            case('CC','C°C')
!                Profile%ST(kkk, lll) = Profile%ST(kkk, lll) * 1d-2 + 273.16d0
!            case('CF','C°F')
!                Profile%ST(kkk, lll) = (Profile%ST(kkk, lll) &
!                    * 1d-2 - 32d0) * 5d0 / 9d0 + 273.16d0
!            case default
!        end select
!    end do
!    do jj = T_1_1_1, T_5_7_1
!        if (MOD(jj,MaxProfRep) /= 0) then
!            kkk = MOD((jj - T_1_1_1 + 1), MaxProfRep)
!            lll = (jj - T_1_1_1 + 1) / MaxProfRep + 1
!        else
!            kkk = MaxProfRep
!            lll = (jj - T_1_1_1 + 1) / MaxProfRep
!        end if
!        select case(ProfileUnits%T(kkk))
!            case('C','°C')
!                Profile%T(kkk, lll) = Profile%T(kkk, lll) + 273.16d0
!            case('F','°F')
!                Profile%T(kkk, lll) = (Profile%T(kkk, lll) &
!                    - 32d0) * 5d0 / 9d0 + 273.16d0
!            case('CK')
!                Profile%T(kkk, lll) = Profile%T(kkk, lll) * 1d-2
!            case('CC','C°C')
!                Profile%T(kkk, lll) = Profile%T(kkk, lll) * 1d-2 + 273.16d0
!            case('CF','C°F')
!                Profile%T(kkk, lll) = (Profile%T(kkk, lll) * 1d-2 - 32d0) &
!                    * 5d0 / 9d0 + 273.16d0
!            case default
!        end select
!    end do
!
!    !> Profile concentrations
!    !> CO2 is taken to ppm
!    do jj = CO2_1_1_1, CO2_5_7_1
!        if (MOD(jj,MaxProfRep) /= 0) then
!            kkk = MOD((jj - CO2_1_1_1 + 1), MaxProfRep)
!            lll = (jj - CO2_1_1_1 + 1) / MaxProfRep + 1
!        else
!            kkk = MaxProfRep
!            lll = (jj - CO2_1_1_1 + 1) / MaxProfRep
!        end if
!        select case(ProfileUnits%CO2(kkk))
!            case('PPT','PPTD','MMOLMOL-1','MMOL/MOL','MMOLMOL^-1')
!                Profile%CO2(kkk, lll) = Profile%CO2(kkk, lll) * 1d3
!            case('PPBD','NMOLMOL-1','NMOL/MOL','NMOLMOL^-1')
!                Profile%CO2(kkk, lll) = Profile%CO2(kkk, lll) * 1d-3
!            case default
!        end select
!    end do
!    !> H2O is taken to ppt
!    do jj = H2O_1_1_1, H2O_5_7_1
!        if (MOD(jj,MaxProfRep) /= 0) then
!            kkk = MOD((jj - H2O_1_1_1 + 1), MaxProfRep)
!            lll = (jj - H2O_1_1_1 + 1) / MaxProfRep + 1
!        else
!            kkk = MaxProfRep
!            lll = (jj - H2O_1_1_1 + 1) / MaxProfRep
!        end if
!        select case(ProfileUnits%H2O(kkk))
!            case('PPM','PPMD','UMOLMOL-1','UMOL/MOL','UMOLMOL^-1')
!                Profile%H2O(kkk, lll) = Profile%H2O(kkk, lll) * 1d-3
!            case('PPBD','NMOLMOL-1','NMOL/MOL','NMOLMOL^-1')
!                Profile%H2O(kkk, lll) = Profile%H2O(kkk, lll) * 1d-6
!            case default
!        end select
!    end do
!    !> CH4 is taken to ppm
!    do jj = CH4_1_1_1, CH4_5_7_1
!        if (MOD(jj,MaxProfRep) /= 0) then
!            kkk = MOD((jj - CH4_1_1_1 + 1), MaxProfRep)
!            lll = (jj - CH4_1_1_1 + 1) / MaxProfRep + 1
!        else
!            kkk = MaxProfRep
!            lll = (jj - CH4_1_1_1 + 1) / MaxProfRep
!        end if
!        select case(ProfileUnits%CH4(kkk))
!            case('PPT','PPTD','MMOLMOL-1','MMOL/MOL','MMOLMOL^-1')
!                Profile%CH4(kkk, lll) = Profile%CH4(kkk, lll) * 1d3
!            case('PPBD','NMOLMOL-1','NMOL/MOL','NMOLMOL^-1')
!                Profile%CH4(kkk, lll) = Profile%CH4(kkk, lll) * 1d-3
!            case default
!        end select
!    end do
!    !> GAS4 is taken to ppm
!    do jj = GAS4_1_1_1, GAS4_5_7_1
!        if (MOD(jj,MaxProfRep) /= 0) then
!            kkk = MOD((jj - GAS4_1_1_1 + 1), MaxProfRep)
!            lll = (jj - GAS4_1_1_1 + 1) / MaxProfRep + 1
!        else
!            kkk = MaxProfRep
!            lll = (jj - GAS4_1_1_1 + 1) / MaxProfRep
!        end if
!        select case(ProfileUnits%GAS4(kkk))
!            case('PPT','PPTD','MMOLMOL-1','MMOL/MOL','MMOLMOL^-1')
!                Profile%GAS4(kkk, lll) = Profile%GAS4(kkk, lll) * 1d3
!            case('PPBD','NMOLMOL-1','NMOL/MOL','NMOLMOL^-1')
!                Profile%GAS4(kkk, lll) = Profile%GAS4(kkk, lll) * 1d-3
!            case default
!        end select
!    end do
end subroutine BiometStandardUnits

!***************************************************************************
!
! \brief       Retrieve timestamp info from timestamp strings in the file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BiometDateTime(pattern, string, date, time)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: pattern
    character(*), intent(in) :: string
    character(*), intent(out) :: date
    character(*), intent(out) :: time
    !> local variables
    integer :: loc_doy
    integer :: start
    integer :: int_year
    logical :: isleap


    date(5:5) = '-'
    date(8:8) = '-'
    time(3:3) = ':'

    !> year
    if (index(Pattern, 'yyyy') /= 0) then
        start = index(Pattern, 'yyyy')
        date(1:4) = string(start: start + 3)
    else
        if (index(Pattern, 'yy') /= 0) then
            start = index(Pattern, 'yy')
            if (string(start: start + 1) > '70') then
                date(1:4) = '19' // string(start: start + 1)
            else
                date(1:4) = '20' // string(start: start + 1)
            end if
        else
            date(1:4) = 'xxxx'
        end if
    end if

    if (date(1:4) /= 'xxxx') then
        read(date(1:4), '(i4)') int_year
        isleap = leapyear(int_year)
    else
        isleap = .false.
    end if

    !> month
    if (index(Pattern, 'mm') /= 0) then
        start = index(Pattern, 'mm')
        date(6:7) = string(start : start + 3)
    else
        date(6:7) = 'xx'
    end if

    !> day or DOY
    if (index(Pattern, 'ddd') /= 0) then
        start = index(Pattern, 'ddd')
        call Char2Int(string(start:start + 2), loc_doy, 3)
        if (isleap) then
            date(6:10) = DayOfLeapYear(loc_doy)
        else
            date(6:10) = DayOfYear(loc_doy)
        end if
    elseif (index(Pattern, 'dd') /= 0) then
        start = index(Pattern, 'dd')
        date(9:10) = string(start : start + 1)
    end if

    !> hour
    if (index(Pattern, 'HH') /= 0) then
        start = index(Pattern, 'HH')
        time(1:2) = string(start : start + 1)
    else
        time(1:2) = 'xx'
    end if

    !> minute
    if (index(Pattern, 'MM') /= 0) then
        start = index(Pattern, 'MM')
        time(4:5) = string(start : start + 1)
    else
        time(4:5) = 'xx'
    end if
end subroutine BiometDateTime


!***************************************************************************
!
! \brief       Retrieve time step intrinsic in Biomet file (in minutes)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjustBiometTimestamps(N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: i
    type(datetype) :: BioometTimestamp

    select case (BiometSetup%tstamp_ref)
        case ('begin')
            do i = 1, N
                call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, BioometTimestamp)
                BioometTimestamp = BioometTimestamp + datetype(0, 0, 0, 0, BiometSetup%tstep)
                call DateTypeToDateTime(BioometTimestamp, Biomet(i)%date, Biomet(i)%time)
            end do

        case ('middle')
            do i = 1, N
                call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, BioometTimestamp)
                BioometTimestamp = BioometTimestamp + datetype(0, 0, 0, 0, BiometSetup%tstep/2)
                call DateTypeToDateTime(BioometTimestamp, Biomet(i)%date, Biomet(i)%time)
            end do
        case ('end')
        return
    end select
end subroutine AdjustBiometTimestamps
