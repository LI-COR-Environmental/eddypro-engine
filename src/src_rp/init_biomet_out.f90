!***************************************************************************
! init_biomet_out.f90
! -------------------
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
! \brief       Initializes biomet output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitBiometOut()
    use m_rp_global_var
    use iso_fortran_env
    implicit none

    !> Local variables
    integer :: open_status
    integer :: i
    integer :: dot
    integer :: rep
    character(2) :: repc
    character(256) :: Test_Path
    character(10000) :: header1 = ''
    character(10000) :: header2 = ''
    character(10000) :: head1_utf8 = ''
    character(10000) :: head2_utf8 = ''
    character(10000) :: head3_utf8 = ''
!    character(1) :: profc
!    integer :: prof

    !> Biomet measurements
    if (EddyProProj%out_biomet .and. NumBiometVar > 0) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyProProj%id(1:len_trim(EddyProProj%id)) &
                  // Biomet_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Slow_Path = Test_Path(1:dot) // CsvTmpExt
        open(uslow, file = Slow_Path, iostat = open_status, encoding = 'utf-8')

        !> Initialize string to void
        call Clearstr(header1)
        call Clearstr(header2)
        call Clearstr(head1_utf8)
        call Clearstr(head2_utf8)
        call Clearstr(head3_utf8)

        call AddDatum(header1,'date,time,DOY', separator)
        call AddDatum(header2,'[yyyy-mm-dd],[HH:MM],[ddd.ddd]', separator)

        do rep = 1, MaxBiometRep
            if (BiometUnits%Ta(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Ta_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[K]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Pa(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Pa_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[Pa]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%RH(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'RH_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[%]', separator)
            end if
        end do
        !> Radiations
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rg(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Rg_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rn(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Rn_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rd(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Rd_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rr(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Rr_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Ruva(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Ruva_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Ruvb(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Ruvb_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%LWin(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'LWin_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%LWout(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'LWout_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWin(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SWin_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWout(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SWout_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWbc(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SWbc_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWdif(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SWdif_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do

        !> Related to radiations
        do rep = 1, MaxBiometRep
            if (BiometUnits%ALB(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'ALB_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[#]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PRI(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'PRI_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[#]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%LAI(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'LAI_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m+2m-2]', separator)
            end if
        end do

        !> PPFD's
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFD(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'PPFD_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[' // char(181) // 'mol+1m-2s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFDd(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'PPFDd_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[' // char(181) // 'mol+1m-2s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFDr(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'PPFDr_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[' // char(181) // 'mol+1m-2s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFDbc(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'PPFDbc_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[' // char(181) // 'mol+1m-2s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%APAR(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'APAR_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[' // char(181) // 'mol+1m-2s-1]', separator)
            end if
        end do

        !> Temperatures
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tc(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Tc_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[K]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tbole(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Tbole_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[K]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tbc(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Tbc_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[K]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tr(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Tr_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[K]', separator)
            end if
        end do

        !> Precipitations
        do rep = 1, MaxBiometRep
            if (BiometUnits%P(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'P_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Prain(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Prain_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Psnow(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Psnow_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SNOWD(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SNOWD_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m]', separator)
            end if
        end do

        !> Others
        do rep = 1, MaxBiometRep
            if (BiometUnits%WS(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'WS_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m+1s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%MWS(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'MWS_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m+1s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%WD(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'WD_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[deg]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SAPFLOW(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SAPFLOW_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[mmol+1m-2s-1]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%STEMFLOW(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'STEMFLOW_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWC(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SWC_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[m+3m-3]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SHF(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'SHF_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[W+1m-2]', separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%TS(rep) /= 'none') then
                call int2char(rep, repc, 2)
                if (rep < 10) repc = repc(2:2)
                call AddDatum(header1,'Ts_' // repc(1:len_trim(repc)) // '_1_1', separator)
                call AddDatum(header2,'[K]', separator)
            end if
        end do

        !> Profiles
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%SWC(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'SWC_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[%]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%SHF(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'SHF_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[W+1m-2]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%ST(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'Tsoil_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[K]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%T(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'Tprof_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[K]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%CO2(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'CO2_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[ppm]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%H2O(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'H2O_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[ppt]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%CH4(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,'CH4_' // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[ppm]', separator)
!                end if
!            end do
!        end do
!        do rep = 1, MaxProfRep
!            do prof = 1, MaxProfNodes
!                if (ProfileUnits%GAS4(rep) /= 'none') then
!                    call int2char(rep, repc, 2)
!                    call int2char(prof, profc, 1)
!                    if (rep < 10) repc = repc(2:2)
!                    call AddDatum(header1,e2sg(gas4)(1:len_trim(e2sg(gas4))) &
!                        // repc(1:len_trim(repc)) // '_' // profc(1:len_trim(profc)) // '_1', separator)
!                    call AddDatum(header2,'[ppm]', separator)
!                end if
!            end do
!        end do

        !> Custom variables
        if (n_cstm_biomet > 0) then
            do i = 1, n_cstm_biomet
                call AddDatum(header1,CstmOrd(i)%var(1:len_trim(CstmOrd(i)%var)), separator)
                call AddDatum(header2,'[' // CstmOrd(i)%units(1:len_trim(CstmOrd(i)%units)) // ']', separator)
            end do
        end if
        call latin1_to_utf8(header1, head1_utf8)
        call latin1_to_utf8(header2, head2_utf8)

        !> Write on output file
        write(uslow, '(a)') head1_utf8(1:len_trim(head1_utf8) - 1)
        write(uslow, '(a)') head2_utf8(1:len_trim(head2_utf8) - 1)
    end if

end subroutine InitBiometOut
