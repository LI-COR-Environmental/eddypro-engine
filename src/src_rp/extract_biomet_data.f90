!***************************************************************************
! extract_biomet_data.f90
! -----------------------
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
! \brief       Interpret content of embedded biomet file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExtractEmbeddedBiometData(Set, N, M, bN)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    integer, intent(in) :: bN
    real(kind=dbl), intent(in) :: Set(N, M)
    !> in/out variables
    integer :: i
    integer :: cnt(38)
    integer, parameter :: lTA       = 1
    integer, parameter :: lTC       = 2
    integer, parameter :: lTBC      = 3
    integer, parameter :: lTS       = 4
    integer, parameter :: lTBOLE    = 5
    integer, parameter :: lRH       = 6
    integer, parameter :: lPA       = 7
    integer, parameter :: lRG       = 8
    integer, parameter :: lRN       = 9
    integer, parameter :: lRD       = 10
    integer, parameter :: lRR       = 11
    integer, parameter :: lLWIN     = 12
    integer, parameter :: lLWOUT    = 13
    integer, parameter :: lSWIN     = 14
    integer, parameter :: lSWout    = 15
    integer, parameter :: lSWBC     = 16
    integer, parameter :: lSWDIF    = 17
    integer, parameter :: lRUVA     = 18
    integer, parameter :: lRUVB     = 19
    integer, parameter :: lP        = 20
    integer, parameter :: lPRAIN    = 21
    integer, parameter :: lPSNOW    = 22
    integer, parameter :: lPPFD     = 23
    integer, parameter :: lPPFDD    = 24
    integer, parameter :: lPPFDR    = 25
    integer, parameter :: lPPFDBC   = 26
    integer, parameter :: lAPAR     = 27
    integer, parameter :: lALB      = 28
    integer, parameter :: lPRI      = 29
    integer, parameter :: lLAI      = 30
    integer, parameter :: lSAPFLOW  = 31
    integer, parameter :: lSTEMFLOW = 32
    integer, parameter :: lTR       = 33
    integer, parameter :: lWS       = 34
    integer, parameter :: lMWS      = 35
    integer, parameter :: lWD       = 36
    integer, parameter :: lSWC      = 37
    integer, parameter :: lSHF      = 38

    !> Initialization
    biomet(1:MaxNumBiometRow) = ErrBiomet
    BiometUnits = NullBiometUnits
    ProfileUnits = NullProfileUnits

    !> Extract time stamps
    Biomet(1:bN)%date = ddvec(1:bN)
    Biomet(1:bN)%time = ttvec(1:bN)

    !> Extract variables
    cnt = 0
    NumSlowVar = 0
    n_cstm_biomet = 0
    do i = NumBiometDateCol + 1, NumBiometCol
        !> Temperatures
        if (BiometMeta%BiometCol(i)%var == 'TA') then
            cnt(lTA) = cnt(lTA) + 1
            if (cnt(lTA) <= 10) then
                Biomet(1:bN)%Ta(cnt(lTA)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Ta(cnt(lTA)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
                !> Check if it was selected for use in computations
                if (BiometSetup%Ta == i) BiometSetup%Ta = cnt(lTA)
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'TC') then
            cnt(lTC) = cnt(lTC) + 1
            if (cnt(lTC) <= 10) then
                Biomet(1:bN)%Tc(cnt(lTC)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Tc(cnt(lTC)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'TBC') then
            cnt(lTBC) = cnt(lTBC) + 1
            if (cnt(lTBC) <= 10) then
                Biomet(1:bN)%Tbc(cnt(lTBC)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Tbc(cnt(lTBC)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'TS') then
            cnt(lTS) = cnt(lTS) + 1
            if (cnt(lTS) <= 10) then
                Biomet(1:bN)%Ts(cnt(lTS)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Ts(cnt(lTS)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'TBOLE') then
            cnt(lTBOLE) = cnt(lTBOLE) + 1
            if (cnt(lTBOLE) <= 10) then
                Biomet(1:bN)%Tbole(cnt(lTBOLE)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Tbole(cnt(lTBOLE)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> RH and Pa
        if (BiometMeta%BiometCol(i)%var == 'RH') then
            cnt(lRH) = cnt(lRH) + 1
            if (cnt(lRH) <= 10) then
                Biomet(1:bN)%RH(cnt(lRH)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%RH(cnt(lRH)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
                !> Check if it was selected for use in computations
                if (BiometSetup%RH == i) BiometSetup%RH = cnt(lRH)
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'PA') then
            cnt(lPA) = cnt(lPA) + 1
            if (cnt(lPA) <= 10) then
                Biomet(1:bN)%Pa(cnt(lPA)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Pa(cnt(lPA)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
                !> Check if it was selected for use in computations
                if (BiometSetup%Pa == i) BiometSetup%Pa = cnt(lPA)
            end if
            cycle
        end if

        !> Radiations
        if (BiometMeta%BiometCol(i)%var == 'RG') then
            cnt(lRG) = cnt(lRG) + 1
            if (cnt(lRG) <= 10) then
                Biomet(1:bN)%Rg(cnt(lRG)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Rg(cnt(lRG)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
                !> Check if it was selected for use in computations
                if (BiometSetup%Rg == i) BiometSetup%Rg = cnt(lRG)
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'RN') then
            cnt(lRN) = cnt(lRN) + 1
            if (cnt(lRN) <= 10) then
                Biomet(1:bN)%Rn(cnt(lRN)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Rn(cnt(lRN)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'RD') then
            cnt(lRD) = cnt(lRD) + 1
            if (cnt(lRD) <= 10) then
                Biomet(1:bN)%Rd(cnt(lRD)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Rd(cnt(lRD)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'RR') then
            cnt(lRR) = cnt(lRR) + 1
            if (cnt(lRR) <= 10) then
                Biomet(1:bN)%Rr(cnt(lRR)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Rr(cnt(lRR)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'LWIN') then
            cnt(lLWIN) = cnt(lLWIN) + 1
            if (cnt(lLWIN) <= 10) then
                Biomet(1:bN)%LWin(cnt(lLWIN)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%LWin(cnt(lLWIN)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
                !> Check if it was selected for use in computations
                if (BiometSetup%LWin == i) BiometSetup%LWin = cnt(lLWIN)
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'LWOUT') then
            cnt(lLWOUT) = cnt(lLWOUT) + 1
            if (cnt(lLWOUT) <= 10) then
                Biomet(1:bN)%LWout(cnt(lLWOUT)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%LWout(cnt(lLWOUT)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'SWIN') then
            cnt(lSWIN) = cnt(lSWIN) + 1
            if (cnt(lSWIN) <= 10) then
                Biomet(1:bN)%SWin(cnt(lSWIN)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SWin(cnt(lSWIN)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'SWOUT') then
            cnt(lSWOUT) = cnt(lSWOUT) + 1
            if (cnt(lSWOUT) <= 10) then
                Biomet(1:bN)%SWout(cnt(lSWOUT)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SWout(cnt(lSWOUT)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'SWBC') then
            cnt(lSWBC) = cnt(lSWBC) + 1
            if (cnt(lSWBC) <= 10) then
                Biomet(1:bN)%SWbc(cnt(lSWBC)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SWbc(cnt(lSWBC)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'SWDIF') then
            cnt(lSWDIF) = cnt(lSWDIF) + 1
            if (cnt(lSWDIF) <= 10) then
                Biomet(1:bN)%SWdif(cnt(lSWDIF)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SWdif(cnt(lSWDIF)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'RUVA') then
            cnt(lRUVA) = cnt(lRUVA) + 1
            if (cnt(lRUVA) <= 10) then
                Biomet(1:bN)%Ruva(cnt(lRUVA)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Ruva(cnt(lRUVA)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'RUVB') then
            cnt(lRUVB) = cnt(lRUVB) + 1
            if (cnt(lRUVB) <= 10) then
                Biomet(1:bN)%Ruvb(cnt(lRUVB)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Ruvb(cnt(lRUVB)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> Precipitations
        if (BiometMeta%BiometCol(i)%var == 'P') then
            cnt(lP) = cnt(lP) + 1
            if (cnt(lP) <= 10) then
                Biomet(1:bN)%P(cnt(lP)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%P(cnt(lP)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'P_RAIN') then
            cnt(lPRAIN) = cnt(lPRAIN) + 1
            if (cnt(lPRAIN) <= 10) then
                Biomet(1:bN)%Prain(cnt(lPRAIN)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Prain(cnt(lPRAIN)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'P_SNOW') then
            cnt(lPSNOW) = cnt(lPSNOW) + 1
            if (cnt(lPSNOW) <= 10) then
                Biomet(1:bN)%Psnow(cnt(lPSNOW)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Psnow(cnt(lPSNOW)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> PPFDs
        if (BiometMeta%BiometCol(i)%var == 'PPFD') then
            cnt(lPPFD) = cnt(lPPFD) + 1
            if (cnt(lPPFD) <= 10) then
                Biomet(1:bN)%PPFD(cnt(lPPFD)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%PPFD(cnt(lPPFD)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
                !> Check if it was selected for use in computations
                if (BiometSetup%PPFD == i) BiometSetup%PPFD = cnt(lPPFD)
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'PPFDD') then
            cnt(lPPFDD) = cnt(lPPFDD) + 1
            if (cnt(lPPFDD) <= 10) then
                Biomet(1:bN)%PPFDd(cnt(lPPFDD)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%PPFDd(cnt(lPPFDD)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'PPFDR') then
            cnt(lPPFDR) = cnt(lPPFDR) + 1
            if (cnt(lPPFDR) <= 10) then
                Biomet(1:bN)%PPFDr(cnt(lPPFDR)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%PPFDr(cnt(lPPFDR)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'PPFDBC') then
            cnt(lPPFDBC) = cnt(lPPFDBC) + 1
            if (cnt(lPPFDBC) <= 10) then
                Biomet(1:bN)%PPFDbc(cnt(lPPFDBC)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%PPFDbc(cnt(lPPFDBC)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'APAR') then
            cnt(lAPAR) = cnt(lAPAR) + 1
            if (cnt(lAPAR) <= 10) then
                Biomet(1:bN)%APAR(cnt(lAPAR)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%APAR(cnt(lAPAR)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> Mixed adimensional params
        if (BiometMeta%BiometCol(i)%var == 'ALB') then
            cnt(lALB) = cnt(lALB) + 1
            if (cnt(lALB) <= 10) then
                Biomet(1:bN)%Alb(cnt(lALB)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%Alb(cnt(lALB)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'PRI') then
            cnt(lPRI) = cnt(lPRI) + 1
            if (cnt(lPRI) <= 10) then
                Biomet(1:bN)%PRI(cnt(lPRI)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%PRI(cnt(lPRI)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'LAI') then
            cnt(lLAI) = cnt(lLAI) + 1
            if (cnt(lLAI) <= 10) then
                Biomet(1:bN)%LAI(cnt(lLAI)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%LAI(cnt(lLAI)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> Trees
        if (BiometMeta%BiometCol(i)%var == 'SAPFLOW') then
            cnt(lSAPFLOW) = cnt(lSAPFLOW) + 1
            if (cnt(lSAPFLOW) <= 10) then
                Biomet(1:bN)%SapFlow(cnt(lSAPFLOW)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SapFlow(cnt(lSAPFLOW)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'STEMFLOW') then
            cnt(lSTEMFLOW) = cnt(lSTEMFLOW) + 1
            if (cnt(lSTEMFLOW) <= 10) then
                Biomet(1:bN)%StemFlow(cnt(lSTEMFLOW)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%StemFlow(cnt(lSTEMFLOW)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'TR') then
            cnt(lTR) = cnt(lTR) + 1
            if (cnt(lTR) <= 10) then
                Biomet(1:bN)%TR(cnt(lTR)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%TR(cnt(lTR)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> Wind
        if (BiometMeta%BiometCol(i)%var == 'WS') then
            cnt(lWS) = cnt(lWS) + 1
            if (cnt(lWS) <= 10) then
                Biomet(1:bN)%WS(cnt(lWS)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%WS(cnt(lWS)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'MWS') then
            cnt(lMWS) = cnt(lMWS) + 1
            if (cnt(lMWS) <= 10) then
                Biomet(1:bN)%MWS(cnt(lMWS)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%MWS(cnt(lMWS)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'WD') then
            cnt(lWD) = cnt(lWD) + 1
            if (cnt(lWD) <= 10) then
                Biomet(1:bN)%WD(cnt(lWD)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%WD(cnt(lWD)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        !> Soil
        if (BiometMeta%BiometCol(i)%var == 'SWC') then
            cnt(lSWC) = cnt(lSWC) + 1
            if (cnt(lSWC) <= 10) then
                Biomet(1:bN)%SWC(cnt(lSWC)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SWC(cnt(lSWC)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if
        if (BiometMeta%BiometCol(i)%var == 'SHF') then
            cnt(lSHF) = cnt(lSHF) + 1
            if (cnt(lSHF) <= 10) then
                Biomet(1:bN)%SHF(cnt(lSHF)) = Set(1:bN, i - NumBiometDateCol)
                BiometUnits%SHF(cnt(lSHF)) = BiometMeta%BiometCol(i)%unit_in
                NumSlowVar = NumSlowVar + 1
            end if
            cycle
        end if

        !> If variable was not found among known ones, then it's a custom variable
        !> so use the ID (not the var) to identify it
        n_cstm_biomet = n_cstm_biomet + 1
        if (len_trim(BiometMeta%BiometCol(i)%id) /= 0 .and. BiometMeta%BiometCol(i)%id /= 'IGNORE') then
            CstmOrd(n_cstm_biomet)%var = BiometMeta%BiometCol(i)%id(1:len_trim(BiometMeta%BiometCol(i)%id))
            CstmOrd(n_cstm_biomet)%units =  BiometMeta%BiometCol(i)%unit_in(1:len_trim(BiometMeta%BiometCol(i)%unit_in))
            CstmBiometSet(1:bN, n_cstm_biomet) = Set(1:bN, i - NumBiometDateCol)
        end if
    end do

    call BiometStandardUnits()
end subroutine ExtractEmbeddedBiometData
