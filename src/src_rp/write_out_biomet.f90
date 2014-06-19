!***************************************************************************
! write_out_biomet.f90
! --------------------
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
! \brief       Write biomet data on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBiomet(init_string)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    !> local variables
    integer :: rep
    character(10000) :: dataline
    character(64) :: datum

    if (EddyProProj%out_biomet .and. NumBiometVar > 0) then
        call clearstr(dataline)
        call AddDatum(dataline, init_string(index(init_string, ',') + 1:len_trim(init_string)), separator)
        do rep = 1, MaxBiometRep
            if (BiometUnits%Ta(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Ta(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Pa(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Pa(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%RH(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%RH(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        !> Radiations
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rg(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Rg(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rn(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Rn(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rd(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Rd(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Rr(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Rr(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Ruva(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Ruva(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Ruvb(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Ruvb(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%LWin(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%LWin(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%LWout(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%LWout(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWin(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SWin(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWout(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SWout(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWbc(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SWbc(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWdif(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SWdif(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        !> Related to radiations
        do rep = 1, MaxBiometRep
            if (BiometUnits%ALB(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%ALB(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PRI(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%PRI(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%LAI(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%LAI(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        !> PPFD's
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFD(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%PPFD(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFDd(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%PPFDd(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFDr(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%PPFDr(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%PPFDbc(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%PPFDbc(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%APAR(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%APAR(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        !> Temperatures
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tc(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Tc(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tbole(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Tbole(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tbc(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Tbc(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Tr(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Tr(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        !> Precipitations
        do rep = 1, MaxBiometRep
            if (BiometUnits%P(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%P(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Prain(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Prain(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%Psnow(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%Psnow(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SNOWD(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SNOWD(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        !> Others
        do rep = 1, MaxBiometRep
            if (BiometUnits%WS(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%WS(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%MWS(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%MWS(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%WD(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%WD(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SAPFLOW(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SAPFLOW(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%STEMFLOW(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%STEMFLOW(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SWC(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SWC(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%SHF(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%SHF(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do
        do rep = 1, MaxBiometRep
            if (BiometUnits%TS(rep) /= 'none') then
                call WriteDatumFloat(E2Biomet%TS(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end if
        end do

        !> Profiles
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%SWC(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%SWC(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%SHF(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%SHF(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%ST(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%ST(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%T(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%T(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%CO2(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%CO2(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%H2O(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%H2O(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%CH4(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%CH4(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
!        do rep = 1, MaxProfRep
!            if (ProfileUnits%GAS4(rep) /= 'none') then
!                do prof = 1, MaxProfNodes
!                    call WriteDatumFloat(E2Profile%GAS4(rep, prof), datum, EddyProProj%err_label)
!                    call AddDatum(dataline, datum, separator)
!                end do
!            end if
!        end do
        !> Custom variables
        if (n_cstm_biomet > 0) then
            do rep = 1, n_cstm_biomet
                call WriteDatumFloat(CstmBiomet(rep), datum, EddyProProj%err_label)
                call AddDatum(dataline, datum, separator)
            end do
        end if

        write(uslow, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

end subroutine WriteOutBiomet
