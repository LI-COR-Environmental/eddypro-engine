!***************************************************************************
! retrieve_biomet_var.f90
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
! \brief       Retrieve Biomet variables for current averaging interval, \n
!              based on timestamp
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveBiometData(EmbBiometDataExist, BiometFileList, NumBiometFiles, LastBiometFile, &
        LastBiometRecord, InitialTimestamp, FinalTimestamp, bN, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bN
    integer, intent(in) :: NumBiometFiles
    type (FileListType), intent(in) :: BiometFileList(NumBiometFiles)
    type (DateType), intent(in) :: InitialTimestamp
    type (DateType), intent(in) :: FinalTimestamp
    integer, intent(inout) :: LastBiometFile
    integer, intent(inout) :: LastBiometRecord
    logical, intent(in) :: EmbBiometDataExist
    logical, intent(in) :: printout
    !> local variables
    integer :: i
    integer :: rep
    integer :: ncstm
!    integer :: co2_ii, co2_jj, h2o_ii, h2o_jj
!    integer :: ch4_ii, ch4_jj, gas4_ii, gas4_jj
    integer :: cslow(4)
!    integer :: cnt(MaxProfNodes)
    integer :: N
    real(kind = dbl) :: tslow(4)
    type (DateType) :: BiometTimestamp
    type (DateType) :: tol
    type (DateType) :: win
    !> Parameters defining plausible ranges for biomet variables
    real(kind = dbl), parameter :: Tmin  = 220d0  !< in K
    real(kind = dbl), parameter :: Tmax  = 330d0  !< in K
    real(kind = dbl), parameter :: PAmin =  40000d0  !< in Pa
    real(kind = dbl), parameter :: PAmax = 120000d0  !< in Pa
    real(kind = dbl), parameter :: RHmin =    0d0  !< in %
    real(kind = dbl), parameter :: lRHmax = 100d0  !< in %
    real(kind = dbl), parameter :: Rmin  =   -10d0  !< in W m-2
    real(kind = dbl), parameter :: Rmax  = 2000d0  !< in W m-2
    real(kind = dbl), parameter :: RNmin = -2000d0  !< in W m-2
    real(kind = dbl), parameter :: RNmax =  2000d0  !< in W m-2
    real(kind = dbl), parameter :: Pmin  =  0d0  !< in m
    real(kind = dbl), parameter :: LAImin =   0d0  !< in #
    real(kind = dbl), parameter :: LAImax =  10d0  !< in #
    real(kind = dbl), parameter :: WSmin  =   0d0  !< in m/s
    real(kind = dbl), parameter :: WSmax  = 100d0  !< in m/s
    real(kind = dbl), parameter :: PRImin  =  -1d0  !< in #
    real(kind = dbl), parameter :: PRImax  =   1d0  !< in #
    real(kind = dbl), parameter :: ALBmin  =   0d0  !< in #
    real(kind = dbl), parameter :: ALBmax  =   1d0  !< in #
    real(kind = dbl), parameter :: PPFDmin = -1000d0  !< in umol/m2/sec to be verified
    real(kind = dbl), parameter :: PPFDmax = 10000d0  !< in umol/m2/sec to be verified
    real(kind = dbl), parameter :: SWCmin  =   0d0  !< in %
    real(kind = dbl), parameter :: SWCmax  = 100d0  !< in %
    logical :: BiometDataExist


    !> If requested, read external file containing biomet measurements
    BiometDataExist = .false.
    if (index(EddyProProj%biomet_data, 'ext_') /= 0) then
        call ReadExtBiometFiles(BiometDataExist, BiometFileList, NumBiometFiles, &
            LastBiometFile, LastBiometRecord, FinalTimestamp, N, printout)
    else
        N = bN
    end if

    if (EmbBiometDataExist .or. BiometDataExist) then
        win = datetype(0, 0, 0, 0, RPsetup%avrg_len)
        tol = datetype(0, 0, 0, 0, max(BiometSetup%tstep, nint(BiometMeta%step)) / 2)

!        !> Retrieve indexes for identifying CO2..GAS4 data from profiles
!        co2_ii = BiometSetup%co2 / 10
!        co2_jj = mod(BiometSetup%co2, 10)
!        h2o_ii = BiometSetup%h2o / 10
!        h2o_jj = mod(BiometSetup%h2o, 10)
!        ch4_ii = BiometSetup%ch4 / 10
!        ch4_jj = mod(BiometSetup%ch4, 10)
!        gas4_ii = BiometSetup%gas4 / 10
!        gas4_jj = mod(BiometSetup%gas4, 10)

        !> Initializations
        cslow = 0
        tslow = 0d0
        CstmBiomet = 0d0
        ncstm = 0
        E2Biomet = NullBiomet
        CountBiomet = NullCountBiomet
        !E2Profile = NullProfile
        !CountProfile = NullCountProfile
        BiometVar = BiometVarType(error, error, error, error, error, error, error, error, &
            error, error, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0)

        do i = 1, N
            !> Convert current Biomet timestamp into datetype
            call DateTimeToDateType(Biomet(i)%date, Biomet(i)%time, BiometTimestamp)

            !> Check if BiometTimestamp is strictly within the averaging period
            if (BiometTimestamp <= FinalTimestamp + tol .and. &
                BiometTimestamp >= FinalTimestamp - win + tol) then
                do rep = 1, MaxBiometRep
                    if (BiometUnits%Ta(rep) /= 'none') then
                        if(Biomet(i)%Ta(rep) >= Tmin .and. Biomet(i)%Ta(rep) <= Tmax) then
                            E2Biomet%Ta(rep) = E2Biomet%Ta(rep) + Biomet(i)%Ta(rep)
                            CountBiomet%Ta(rep) = CountBiomet%Ta(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Pa(rep) /= 'none') then
                        if(Biomet(i)%Pa(rep) >= PAmin .and. Biomet(i)%Pa(rep) <= PAmax) then
                            E2Biomet%Pa(rep) = E2Biomet%Pa(rep) + Biomet(i)%Pa(rep)
                            CountBiomet%Pa(rep) = CountBiomet%Pa(rep) + 1
                        end if
                    end if
                    if (BiometUnits%RH(rep) /= 'none') then
                        if(Biomet(i)%RH(rep) >= RHmin .and. Biomet(i)%RH(rep) <= lRHmax) then
                            E2Biomet%RH(rep) = E2Biomet%RH(rep) + Biomet(i)%RH(rep)
                            CountBiomet%RH(rep) = CountBiomet%RH(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Rg(rep) /= 'none') then
                        if(Biomet(i)%Rg(rep) >= Rmin .and. Biomet(i)%Rg(rep) <= Rmax) then
                            E2Biomet%Rg(rep) = E2Biomet%Rg(rep) + Biomet(i)%Rg(rep)
                            CountBiomet%Rg(rep) = CountBiomet%Rg(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Rn(rep) /= 'none') then
                        if(Biomet(i)%Rn(rep) >= RNmin .and. Biomet(i)%Rn(rep) <= RNmax) then
                            E2Biomet%Rn(rep) = E2Biomet%Rn(rep) + Biomet(i)%Rn(rep)
                            CountBiomet%Rn(rep) = CountBiomet%Rn(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Rd(rep) /= 'none') then
                        if(Biomet(i)%Rd(rep) >= Rmin .and. Biomet(i)%Rd(rep) <= Rmax) then
                            E2Biomet%Rd(rep) = E2Biomet%Rd(rep) + Biomet(i)%Rd(rep)
                            CountBiomet%Rd(rep) = CountBiomet%Rd(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Rr(rep) /= 'none') then
                        if(Biomet(i)%Rr(rep) >= Rmin .and. Biomet(i)%Rr(rep) <= Rmax) then
                            E2Biomet%Rr(rep) = E2Biomet%Rr(rep) + Biomet(i)%Rr(rep)
                            CountBiomet%Rr(rep) = CountBiomet%Rr(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Ruva(rep) /= 'none') then
                        if(Biomet(i)%Ruva(rep) >= Rmin .and. Biomet(i)%Ruva(rep) <= Rmax) then
                            E2Biomet%Ruva(rep) = E2Biomet%Ruva(rep) + Biomet(i)%Ruva(rep)
                            CountBiomet%Ruva(rep) = CountBiomet%Ruva(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Ruvb(rep) /= 'none') then
                        if(Biomet(i)%Ruvb(rep) >= Rmin .and. Biomet(i)%Ruvb(rep) <= Rmax) then
                            E2Biomet%Ruvb(rep) = E2Biomet%Ruvb(rep) + Biomet(i)%Ruvb(rep)
                            CountBiomet%Ruvb(rep) = CountBiomet%Ruvb(rep) + 1
                        end if
                    end if
                    if (BiometUnits%LWin(rep) /= 'none') then
                        if(Biomet(i)%LWin(rep) > Rmin .and. Biomet(i)%LWin(rep) < Rmax) then
                            E2Biomet%LWin(rep) = E2Biomet%LWin(rep) + Biomet(i)%LWin(rep)
                            CountBiomet%LWin(rep) = CountBiomet%LWin(rep) + 1
                        end if
                    end if
                    if (BiometUnits%LWout(rep) /= 'none') then
                        if(Biomet(i)%LWout(rep) >= Rmin .and. Biomet(i)%LWout(rep) <= Rmax) then
                            E2Biomet%LWout(rep) = E2Biomet%LWout(rep) + Biomet(i)%LWout(rep)
                            CountBiomet%LWout(rep) = CountBiomet%LWout(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Tc(rep) /= 'none') then
                        if(Biomet(i)%Tc(rep) >= Tmin .and. Biomet(i)%Tc(rep) <= Tmax) then
                            E2Biomet%Tc(rep) = E2Biomet%Tc(rep) + Biomet(i)%Tc(rep)
                            CountBiomet%Tc(rep) = CountBiomet%Tc(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Tbole(rep) /= 'none') then
                        if(Biomet(i)%Tbole(rep) >= Tmin .and. Biomet(i)%Tbole(rep) <= Tmax) then
                            E2Biomet%Tbole(rep) = E2Biomet%Tbole(rep) + Biomet(i)%Tbole(rep)
                            CountBiomet%Tbole(rep) = CountBiomet%Tbole(rep) + 1
                        end if
                    end if
                    if (BiometUnits%P(rep) /= 'none') then
                        if(Biomet(i)%P(rep) >= Pmin) then
                            E2Biomet%P(rep) = E2Biomet%P(rep) + Biomet(i)%P(rep)
                            CountBiomet%P(rep) = CountBiomet%P(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Prain(rep) /= 'none') then
                        if(Biomet(i)%Prain(rep) >= Pmin) then
                            E2Biomet%Prain(rep) = E2Biomet%Prain(rep) + Biomet(i)%Prain(rep)
                            CountBiomet%Prain(rep) = CountBiomet%Prain(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Psnow(rep) /= 'none') then
                        if(Biomet(i)%Psnow(rep) >= Pmin) then
                            E2Biomet%Psnow(rep) = E2Biomet%Psnow(rep) + Biomet(i)%Psnow(rep)
                            CountBiomet%Psnow(rep) = CountBiomet%Psnow(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SNOWD(rep) /= 'none') then
                        if(Biomet(i)%SNOWD(rep) /= error) then
                            E2Biomet%SNOWD(rep) = E2Biomet%SNOWD(rep) + Biomet(i)%SNOWD(rep)
                            CountBiomet%SNOWD(rep) = CountBiomet%SNOWD(rep) + 1
                        end if
                    end if
                    if (BiometUnits%PPFD(rep) /= 'none') then
                        if(Biomet(i)%PPFD(rep) >= PPFDmin .and. Biomet(i)%PPFD(rep) <= PPFDmax) then
                            E2Biomet%PPFD(rep) = E2Biomet%PPFD(rep) + Biomet(i)%PPFD(rep)
                            CountBiomet%PPFD(rep) = CountBiomet%PPFD(rep) + 1
                        end if
                    end if
                    if (BiometUnits%PPFDd(rep) /= 'none') then
                        if(Biomet(i)%PPFDd(rep) >= PPFDmin .and. Biomet(i)%PPFDd(rep) <= PPFDmax) then
                            E2Biomet%PPFDd(rep) = E2Biomet%PPFDd(rep) + Biomet(i)%PPFDd(rep)
                            CountBiomet%PPFDd(rep) = CountBiomet%PPFDd(rep) + 1
                        end if
                    end if
                    if (BiometUnits%PPFDr(rep) /= 'none') then
                        if(Biomet(i)%PPFDr(rep) >= PPFDmin .and. Biomet(i)%PPFDr(rep) <= PPFDmax) then
                            E2Biomet%PPFDr(rep) = E2Biomet%PPFDr(rep) + Biomet(i)%PPFDr(rep)
                            CountBiomet%PPFDr(rep) = CountBiomet%PPFDr(rep) + 1
                        end if
                    end if
                    if (BiometUnits%PPFDbc(rep) /= 'none') then
                        if(Biomet(i)%PPFDbc(rep) /= error) then
                            E2Biomet%PPFDbc(rep) = E2Biomet%PPFDbc(rep) + Biomet(i)%PPFDbc(rep)
                            CountBiomet%PPFDbc(rep) = CountBiomet%PPFDbc(rep) + 1
                        end if
                    end if
                    if (BiometUnits%APAR(rep) /= 'none') then
                        if(Biomet(i)%APAR(rep) /= error) then
                            E2Biomet%APAR(rep) = E2Biomet%APAR(rep) + Biomet(i)%APAR(rep)
                            CountBiomet%APAR(rep) = CountBiomet%APAR(rep) + 1
                        end if
                    end if
                    if (BiometUnits%ALB(rep) /= 'none') then
                        if(Biomet(i)%ALB(rep) > ALBmin .and. Biomet(i)%ALB(rep) < ALBmax) then
                            E2Biomet%ALB(rep) = E2Biomet%ALB(rep) + Biomet(i)%ALB(rep)
                            CountBiomet%ALB(rep) = CountBiomet%ALB(rep) + 1
                        end if
                    end if
                    if (BiometUnits%PRI(rep) /= 'none') then
                        if(Biomet(i)%PRI(rep) >= PRImin .and. Biomet(i)%PRI(rep) <= PRImax) then
                            E2Biomet%PRI(rep) = E2Biomet%PRI(rep) + Biomet(i)%PRI(rep)
                            CountBiomet%PRI(rep) = CountBiomet%PRI(rep) + 1
                        end if
                    end if
                    if (BiometUnits%LAI(rep) /= 'none') then
                        if(Biomet(i)%LAI(rep) >= LAImin .and. Biomet(i)%LAI(rep) <= LAImax) then
                            E2Biomet%LAI(rep) = E2Biomet%LAI(rep) + Biomet(i)%LAI(rep)
                            CountBiomet%LAI(rep) = CountBiomet%LAI(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SWin(rep) /= 'none') then
                        if(Biomet(i)%SWin(rep) >= Rmin .and. Biomet(i)%SWin(rep) <= Rmax) then
                            E2Biomet%SWin(rep) = E2Biomet%SWin(rep) + Biomet(i)%SWin(rep)
                            CountBiomet%SWin(rep) = CountBiomet%SWin(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SWout(rep) /= 'none') then
                        if(Biomet(i)%SWout(rep) >= Rmin .and. Biomet(i)%SWout(rep) <= Rmax) then
                            E2Biomet%SWout(rep) = E2Biomet%SWout(rep) + Biomet(i)%SWout(rep)
                            CountBiomet%SWout(rep) = CountBiomet%SWout(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SWbc(rep) /= 'none') then
                        if(Biomet(i)%SWbc(rep) >= Rmin .and. Biomet(i)%SWbc(rep) <= Rmax) then
                            E2Biomet%SWbc(rep) = E2Biomet%SWbc(rep) + Biomet(i)%SWbc(rep)
                            CountBiomet%SWbc(rep) = CountBiomet%SWbc(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SWdif(rep) /= 'none') then
                        if(Biomet(i)%SWdif(rep) >= Rmin .and. Biomet(i)%SWdif(rep) <= Rmax) then
                            E2Biomet%SWdif(rep) = E2Biomet%SWdif(rep) + Biomet(i)%SWdif(rep)
                            CountBiomet%SWdif(rep) = CountBiomet%SWdif(rep) + 1
                        end if
                    end if
                    if (BiometUnits%TBC(rep) /= 'none') then
                        if(Biomet(i)%TBC(rep) >= Tmin .and. Biomet(i)%Tbc(rep) <= Tmax) then
                            E2Biomet%TBC(rep) = E2Biomet%TBC(rep) + Biomet(i)%TBC(rep)
                            CountBiomet%TBC(rep) = CountBiomet%TBC(rep) + 1
                        end if
                    end if
                    if (BiometUnits%WS(rep) /= 'none') then
                        if(Biomet(i)%WS(rep) >= WSmin .and. Biomet(i)%WS(rep) <= WSmax) then
                            E2Biomet%WS(rep) = E2Biomet%WS(rep) + Biomet(i)%WS(rep)
                            CountBiomet%WS(rep) = CountBiomet%WS(rep) + 1
                        end if
                    end if
                    if (BiometUnits%MWS(rep) /= 'none') then
                        if(Biomet(i)%MWS(rep) >= WSmin .and. Biomet(i)%MWS(rep) <= WSmax) then
                            E2Biomet%MWS(rep) = E2Biomet%MWS(rep) + Biomet(i)%MWS(rep)
                            CountBiomet%MWS(rep) = CountBiomet%MWS(rep) + 1
                        end if
                    end if
                    if (BiometUnits%WD(rep) /= 'none') then
                        if(Biomet(i)%WD(rep) >= 0d0 .and. Biomet(i)%WD(rep) < 360d0 ) then
                            E2Biomet%WD(rep) = E2Biomet%WD(rep) + Biomet(i)%WD(rep)
                            CountBiomet%WD(rep) = CountBiomet%WD(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SAPFLOW(rep) /= 'none') then
                        if(Biomet(i)%SAPFLOW(rep) /= error) then
                            E2Biomet%SAPFLOW(rep) = E2Biomet%SAPFLOW(rep) + Biomet(i)%SAPFLOW(rep)
                            CountBiomet%SAPFLOW(rep) = CountBiomet%SAPFLOW(rep) + 1
                        end if
                    end if
                    if (BiometUnits%STEMFLOW(rep) /= 'none') then
                        if(Biomet(i)%STEMFLOW(rep) /= error) then
                            E2Biomet%STEMFLOW(rep) = E2Biomet%STEMFLOW(rep) + Biomet(i)%STEMFLOW(rep)
                            CountBiomet%STEMFLOW(rep) = CountBiomet%STEMFLOW(rep) + 1
                        end if
                    end if
                    if (BiometUnits%TR(rep) /= 'none') then
                        if(Biomet(i)%TR(rep) /= error) then
                            E2Biomet%TR(rep) = E2Biomet%TR(rep) + Biomet(i)%TR(rep)
                            CountBiomet%TR(rep) = CountBiomet%TR(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SWC(rep) /= 'none') then
                        if(Biomet(i)%SWC(rep) > SWCmin .and. Biomet(i)%SWC(rep) < SWCmax) then
                            E2Biomet%SWC(rep) = E2Biomet%SWC(rep) + Biomet(i)%SWC(rep)
                            CountBiomet%SWC(rep) = CountBiomet%SWC(rep) + 1
                        end if
                    end if
                    if (BiometUnits%SHF(rep) /= 'none') then
                        if(Biomet(i)%SHF(rep) /= error) then
                            E2Biomet%SHF(rep) = E2Biomet%SHF(rep) + Biomet(i)%SHF(rep)
                            CountBiomet%SHF(rep) = CountBiomet%SHF(rep) + 1
                        end if
                    end if
                    if (BiometUnits%Ts(rep) /= 'none') then
                        if(Biomet(i)%Ts(rep) >= Tmin .and. Biomet(i)%Ts(rep) <= Tmax) then
                            E2Biomet%Ts(rep) = E2Biomet%Ts(rep) + Biomet(i)%Ts(rep)
                            CountBiomet%Ts(rep) = CountBiomet%Ts(rep) + 1
                        end if
                    end if
                end do

!                !> Profiles
!                do rep = 1, MaxProfRep
!                    if (ProfileUnits%SWC(rep) /= 'none') then
!                        E2Profile%SWC(rep, :) = E2Profile%SWC(rep, :) + Profile(i)%SWC(rep, :)
!                        CountProfile%SWC(rep, :) = CountProfile%SWC(rep, :) + 1
!                    end if
!                    if (ProfileUnits%SHF(rep) /= 'none') then
!                        E2Profile%SHF(rep, :) = E2Profile%SHF(rep, :) + Profile(i)%SHF(rep, :)
!                        CountProfile%SHF(rep, :) = CountProfile%SHF(rep, :) + 1
!                    end if
!                    if (ProfileUnits%ST(rep) /= 'none') then
!                        E2Profile%ST(rep, :) = E2Profile%ST(rep, :) + Profile(i)%ST(rep, :)
!                        CountProfile%ST(rep, :) = CountProfile%ST(rep, :) + 1
!                    end if
!                    if (ProfileUnits%T(rep) /= 'none') then
!                        E2Profile%T(rep, :) = E2Profile%T(rep, :) + Profile(i)%T(rep, :)
!                        CountProfile%T(rep, :) = CountProfile%T(rep, :) + 1
!                    end if
!                    if (ProfileUnits%CO2(rep) /= 'none') then
!                        E2Profile%CO2(rep, :) = E2Profile%CO2(rep, :) + Profile(i)%CO2(rep, :)
!                        CountProfile%CO2(rep, :) = CountProfile%CO2(rep, :) + 1
!                    end if
!                    if (ProfileUnits%H2O(rep) /= 'none') then
!                        E2Profile%H2O(rep, :) = E2Profile%H2O(rep, :) + Profile(i)%H2O(rep, :)
!                        CountProfile%H2O(rep, :) = CountProfile%H2O(rep, :) + 1
!                    end if
!                    if (ProfileUnits%CH4(rep) /= 'none') then
!                        E2Profile%CH4(rep, :) = E2Profile%CH4(rep, :) + Profile(i)%CH4(rep, :)
!                        CountProfile%CH4(rep, :) = CountProfile%CH4(rep, :) + 1
!                    end if
!                    if (ProfileUnits%GAS4(rep) /= 'none') then
!                        E2Profile%GAS4(rep, :) = E2Profile%GAS4(rep, :) + Profile(i)%GAS4(rep, :)
!                        CountProfile%GAS4(rep, :) = CountProfile%GAS4(rep, :) + 1
!                    end if
!                end do

!                    !> Concentrations
!                    if (BiometSetup%CO2 > 0) then
!                        if (Profile(i)%CO2(co2_ii, co2_jj) /= error) then
!                        tslow(1) = tslow(1)  + Profile(i)%CO2(co2_ii, co2_jj)
!                        cslow(1) = cslow(1) + 1
!                        end if
!                    end if
!                    if (BiometSetup%H2O > 0) then
!                        if (Profile(i)%H2O(h2o_ii, h2o_jj) /= error) then
!                        tslow(2) = tslow(2)  + Profile(i)%H2O(h2o_ii, h2o_jj)
!                        cslow(2) = cslow(2) + 1
!                        end if
!                    end if
!                    if (BiometSetup%CH4 > 0) then
!                        if (Profile(i)%CH4(ch4_ii, ch4_jj) /= error) then
!                        tslow(3) = tslow(3)  + Profile(i)%CH4(ch4_ii, ch4_jj)
!                        cslow(3) = cslow(3) + 1
!                        end if
!                    end if
!                    if (BiometSetup%GAS4 > 0) then
!                        if (Profile(i)%GAS4(gas4_ii, gas4_jj) /= error) then
!                        tslow(4) = tslow(4)  + Profile(i)%GAS4(gas4_ii, gas4_jj)
!                        cslow(4) = cslow(4) + 1
!                        end if
!                    end if

                    !> Custom variables
                    if (n_cstm_biomet > 0) then
                        ncstm = ncstm + 1
                        CstmBiomet(1:n_cstm_biomet) = &
                            CstmBiomet(1:n_cstm_biomet) + CstmBiometSet(i, 1:n_cstm_biomet)
                    end if
            end if
        end do

        !> Normalize to get the averages
        where (CountBiomet%Ta(:) /= 0)
            E2Biomet%Ta(:) = E2Biomet%Ta(:) / CountBiomet%Ta(:)
        elsewhere
            E2Biomet%Ta(:) = error
        end where
        where (CountBiomet%Pa(:) /= 0)
            E2Biomet%Pa(:) = E2Biomet%Pa(:) / CountBiomet%Pa(:)
        elsewhere
            E2Biomet%Pa(:) = error
        end where
        where (CountBiomet%RH(:) /= 0)
            E2Biomet%RH(:) = E2Biomet%RH(:) / CountBiomet%RH(:)
        elsewhere
            E2Biomet%RH(:) = error
        end where
        where (CountBiomet%Rg(:) /= 0)
            E2Biomet%Rg(:) = E2Biomet%Rg(:) / CountBiomet%Rg(:)
        elsewhere
            E2Biomet%Rg(:) = error
        end where
        where (CountBiomet%Rn(:) /= 0)
            E2Biomet%Rn(:) = E2Biomet%Rn(:) / CountBiomet%Rn(:)
        elsewhere
            E2Biomet%Rn(:) = error
        end where
        where (CountBiomet%Rd(:) /= 0)
            E2Biomet%Rd(:) = E2Biomet%Rd(:) / CountBiomet%Rd(:)
        elsewhere
            E2Biomet%Rd(:) = error
        end where
        where (CountBiomet%Rr(:) /= 0)
            E2Biomet%Rr(:) = E2Biomet%Rr(:) / CountBiomet%Rr(:)
        elsewhere
            E2Biomet%Rr(:) = error
        end where
        where (CountBiomet%Ruva(:) /= 0)
            E2Biomet%Ruva(:) = E2Biomet%Ruva(:) / CountBiomet%Ruva(:)
        elsewhere
            E2Biomet%Ruva(:) = error
        end where
        where (CountBiomet%Ruvb(:) /= 0)
            E2Biomet%Ruvb(:) = E2Biomet%Ruvb(:) / CountBiomet%Ruvb(:)
        elsewhere
            E2Biomet%Ruvb(:) = error
        end where
        where (CountBiomet%LWin(:) /= 0)
            E2Biomet%LWin(:) = E2Biomet%LWin(:) / CountBiomet%LWin(:)
        elsewhere
            E2Biomet%LWin(:) = error
        end where
        where (CountBiomet%LWout(:) /= 0)
            E2Biomet%LWout(:) = E2Biomet%LWout(:) / CountBiomet%LWout(:)
        elsewhere
            E2Biomet%LWout(:) = error
        end where
        where (CountBiomet%Tc(:) /= 0)
            E2Biomet%Tc(:) = E2Biomet%Tc(:) / CountBiomet%Tc(:)
        elsewhere
            E2Biomet%Tc(:) = error
        end where
        where (CountBiomet%Tbole(:) /= 0)
            E2Biomet%Tbole(:) = E2Biomet%Tbole(:) / CountBiomet%Tbole(:)
        elsewhere
            E2Biomet%Tbole(:) = error
        end where
        where (CountBiomet%P(:) /= 0)
!            E2Biomet%P(:) = E2Biomet%P(:) / CountBiomet%P(:)
        elsewhere
            E2Biomet%P(:) = error
        end where
        where (CountBiomet%Prain(:) /= 0)
!            E2Biomet%Prain(:) = E2Biomet%Prain(:) / CountBiomet%Prain(:)
        elsewhere
            E2Biomet%Prain(:) = error
        end where
        where (CountBiomet%Psnow(:) /= 0)
!            E2Biomet%Psnow(:) = E2Biomet%Psnow(:) / CountBiomet%Psnow(:)
        elsewhere
            E2Biomet%Psnow(:) = error
        end where
        where (CountBiomet%SNOWD(:) /= 0)
            E2Biomet%SNOWD(:) = E2Biomet%SNOWD(:) / CountBiomet%SNOWD(:)
        elsewhere
            E2Biomet%SNOWD(:) = error
        end where
        where (CountBiomet%PPFD(:) /= 0)
            E2Biomet%PPFD(:) = E2Biomet%PPFD(:) / CountBiomet%PPFD(:)
        elsewhere
            E2Biomet%PPFD(:) = error
        end where
        where (CountBiomet%PPFDd(:) /= 0)
            E2Biomet%PPFDd(:) = E2Biomet%PPFDd(:) / CountBiomet%PPFDd(:)
        elsewhere
            E2Biomet%PPFDd(:) = error
        end where
        where (CountBiomet%PPFDr(:) /= 0)
            E2Biomet%PPFDr(:) = E2Biomet%PPFDr(:) / CountBiomet%PPFDr(:)
        elsewhere
            E2Biomet%PPFDr(:) = error
        end where
        where (CountBiomet%PPFDbc(:) /= 0)
            E2Biomet%PPFDbc(:) = E2Biomet%PPFDbc(:) / CountBiomet%PPFDbc(:)
        elsewhere
            E2Biomet%PPFDbc(:) = error
        end where
        where (CountBiomet%APAR(:) /= 0)
            E2Biomet%APAR(:) = E2Biomet%APAR(:) / CountBiomet%APAR(:)
        elsewhere
            E2Biomet%APAR(:) = error
        end where
        where (CountBiomet%ALB(:) /= 0)
            E2Biomet%ALB(:) = E2Biomet%ALB(:) / CountBiomet%ALB(:)
        elsewhere
            E2Biomet%ALB(:) = error
        end where
        where (CountBiomet%PRI(:) /= 0)
            E2Biomet%PRI(:) = E2Biomet%PRI(:) / CountBiomet%PRI(:)
        elsewhere
            E2Biomet%PRI(:) = error
        end where
        where (CountBiomet%LAI(:) /= 0)
            E2Biomet%LAI(:) = E2Biomet%LAI(:) / CountBiomet%LAI(:)
        elsewhere
            E2Biomet%LAI(:) = error
        end where
        where (CountBiomet%SWin(:) /= 0)
            E2Biomet%SWin(:) = E2Biomet%SWin(:) / CountBiomet%SWin(:)
        elsewhere
            E2Biomet%SWin(:) = error
        end where
        where (CountBiomet%SWout(:) /= 0)
            E2Biomet%SWout(:) = E2Biomet%SWout(:) / CountBiomet%SWout(:)
        elsewhere
            E2Biomet%SWout(:) = error
        end where
        where (CountBiomet%SWbc(:) /= 0)
            E2Biomet%SWbc(:) = E2Biomet%SWbc(:) / CountBiomet%SWbc(:)
        elsewhere
            E2Biomet%SWbc(:) = error
        end where
        where (CountBiomet%SWdif(:) /= 0)
            E2Biomet%SWdif(:) = E2Biomet%SWdif(:) / CountBiomet%SWdif(:)
        elsewhere
            E2Biomet%SWdif(:) = error
        end where
        where (CountBiomet%TBC(:) /= 0)
            E2Biomet%TBC(:) = E2Biomet%TBC(:) / CountBiomet%TBC(:)
        elsewhere
            E2Biomet%TBC(:) = error
        end where
        where (CountBiomet%WS(:) /= 0)
            E2Biomet%WS(:) = E2Biomet%WS(:) / CountBiomet%WS(:)
        elsewhere
            E2Biomet%WS(:) = error
        end where
        where (CountBiomet%MWS(:) /= 0)
            E2Biomet%MWS(:) = E2Biomet%MWS(:) / CountBiomet%MWS(:)
        elsewhere
            E2Biomet%MWS(:) = error
        end where
        where (CountBiomet%WD(:) /= 0)
            E2Biomet%WD(:) = E2Biomet%WD(:) / CountBiomet%WD(:)
        elsewhere
            E2Biomet%WD(:) = error
        end where
        where (CountBiomet%SAPFLOW(:) /= 0)
            E2Biomet%SAPFLOW(:) = E2Biomet%SAPFLOW(:) / CountBiomet%SAPFLOW(:)
        elsewhere
            E2Biomet%SAPFLOW(:) = error
        end where
        where (CountBiomet%STEMFLOW(:) /= 0)
            E2Biomet%STEMFLOW(:) = E2Biomet%STEMFLOW(:) / CountBiomet%STEMFLOW(:)
        elsewhere
            E2Biomet%STEMFLOW(:) = error
        end where
        where (CountBiomet%TR(:) /= 0)
            E2Biomet%TR(:) = E2Biomet%TR(:) / CountBiomet%TR(:)
        elsewhere
            E2Biomet%TR(:) = error
        end where
        where (CountBiomet%SWC(:) /= 0)
            E2Biomet%SWC(:) = E2Biomet%SWC(:) / CountBiomet%SWC(:)
        elsewhere
            E2Biomet%SWC(:) = error
        end where
        where (CountBiomet%SHF(:) /= 0)
            E2Biomet%SHF(:) = E2Biomet%SHF(:) / CountBiomet%SHF(:)
        elsewhere
            E2Biomet%SHF(:) = error
        end where
        where (CountBiomet%TS(:) /= 0)
            E2Biomet%TS(:) = E2Biomet%TS(:) / CountBiomet%TS(:)
        elsewhere
            E2Biomet%TS(:) = error
        end where
        !> Special case of concentrations
!        where (cslow(:) > 0)
!            tslow(:) = tslow(:) / cslow(:)
!        elsewhere
!            tslow(:) = error
!        end where
!        !> Profiles
!        where (CountProfile%SWC /= 0)
!            E2Profile%SWC = E2Profile%SWC / CountProfile%SWC
!        elsewhere
!            E2Profile%SWC = error
!        end where
!        where (CountProfile%SHF /= 0)
!            E2Profile%SHF = E2Profile%SHF / CountProfile%SHF
!        elsewhere
!            E2Profile%SHF = error
!        end where
!        where (CountProfile%ST /= 0)
!            E2Profile%ST = E2Profile%ST / CountProfile%ST
!        elsewhere
!            E2Profile%ST = error
!        end where
!        where (CountProfile%T /= 0)
!            E2Profile%T = E2Profile%T / CountProfile%T
!        elsewhere
!            E2Profile%T = error
!        end where
!        where (CountProfile%CO2 /= 0)
!            E2Profile%CO2 = E2Profile%CO2 / CountProfile%CO2
!        elsewhere
!            E2Profile%CO2 = error
!        end where
!        where (CountProfile%H2O /= 0)
!            E2Profile%H2O = E2Profile%H2O / CountProfile%H2O
!        elsewhere
!            E2Profile%H2O = error
!        end where
!        where (CountProfile%CH4 /= 0)
!            E2Profile%CH4 = E2Profile%CH4 / CountProfile%CH4
!        elsewhere
!            E2Profile%CH4 = error
!        end where
!        where (CountProfile%GAS4 /= 0)
!            E2Profile%GAS4 = E2Profile%GAS4 / CountProfile%GAS4
!        elsewhere
!            E2Profile%GAS4 = error
!        end where
        !> Custom variables
        if (n_cstm_biomet > 0 .and. ncstm > 0) then
            CstmBiomet = CstmBiomet / ncstm
        end if

        !> Associate values to variables, as selected by user
        if (BiometSetup%Ta > 0) BiometVar%Ta = E2Biomet%Ta(BiometSetup%Ta)
        if (BiometSetup%Pa > 0) BiometVar%Pa = E2Biomet%Pa(BiometSetup%Pa)
        if (BiometSetup%RH > 0) BiometVar%RH = E2Biomet%RH(BiometSetup%RH)
        if (BiometSetup%Rg > 0) BiometVar%Rg = E2Biomet%Rg(BiometSetup%Rg)
        if (BiometSetup%LWin > 0) BiometVar%LWin = E2Biomet%LWin(BiometSetup%LWin)
        if (BiometSetup%PPFD > 0) BiometVar%PPFD = E2Biomet%PPFD(BiometSetup%PPFD)
        BiometVar%CO2  = tslow(1)
        BiometVar%H2O  = tslow(2)
        BiometVar%CH4  = tslow(3)
        BiometVar%GAS4 = tslow(4)

        !> Profiles
!        if (BiometSetup%prof_swc > 0) then
!            BiometVar%prof_swc(:) = E2Profile%SWC(BiometSetup%prof_swc, :)
!        elseif (BiometSetup%prof_swc == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%SWC(j, k) /= error) then
!                        BiometVar%prof_swc(k) = BiometVar%prof_swc(k) + E2Profile%SWC(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_swc(:) = BiometVar%prof_swc(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_swc(:) = error
!        end if
!        if (BiometSetup%prof_shf > 0) then
!            BiometVar%prof_shf(:) = E2Profile%shf(BiometSetup%prof_shf, :)
!        elseif (BiometSetup%prof_shf == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%shf(j, k) /= error) then
!                        BiometVar%prof_shf(k) = BiometVar%prof_shf(k) + E2Profile%shf(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_shf(:) = BiometVar%prof_shf(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_shf(:) = error
!        end if
!        if (BiometSetup%prof_ts > 0) then
!            BiometVar%prof_ts(:) = E2Profile%st(BiometSetup%prof_ts, :)
!        elseif (BiometSetup%prof_ts == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%st(j, k) /= error) then
!                        BiometVar%prof_ts(k) = BiometVar%prof_ts(k) + E2Profile%st(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_ts(:) = BiometVar%prof_ts(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_ts(:) = error
!        end if
!        if (BiometSetup%prof_t > 0) then
!            BiometVar%prof_t(:) = E2Profile%t(BiometSetup%prof_t, :)
!        elseif (BiometSetup%prof_t == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%t(j, k) /= error) then
!                        BiometVar%prof_t(k) = BiometVar%prof_t(k) + E2Profile%t(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_t(:) = BiometVar%prof_t(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_t(:) = error
!        end if
!        if (BiometSetup%prof_co2 > 0) then
!            BiometVar%prof_co2(:) = E2Profile%co2(BiometSetup%prof_co2, :)
!        elseif (BiometSetup%prof_co2 == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%co2(j, k) /= error) then
!                        BiometVar%prof_co2(k) = BiometVar%prof_co2(k) + E2Profile%co2(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_co2(:) = BiometVar%prof_co2(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_co2(:) = error
!        end if
!        if (BiometSetup%prof_h2o > 0) then
!            BiometVar%prof_h2o(:) = E2Profile%h2o(BiometSetup%prof_h2o, :)
!        elseif (BiometSetup%prof_h2o == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%h2o(j, k) /= error) then
!                        BiometVar%prof_h2o(k) = BiometVar%prof_h2o(k) + E2Profile%h2o(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_h2o(:) = BiometVar%prof_h2o(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_h2o(:) = error
!        end if
!        if (BiometSetup%prof_ch4 > 0) then
!            BiometVar%prof_ch4(:) = E2Profile%ch4(BiometSetup%prof_ch4, :)
!        elseif (BiometSetup%prof_ch4 == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%ch4(j, k) /= error) then
!                        BiometVar%prof_ch4(k) = BiometVar%prof_ch4(k) + E2Profile%ch4(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_ch4(:) = BiometVar%prof_ch4(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_ch4(:) = error
!        end if
!        if (BiometSetup%prof_gas4 > 0) then
!            BiometVar%prof_gas4(:) = E2Profile%gas4(BiometSetup%prof_gas4, :)
!        elseif (BiometSetup%prof_gas4 == 0) then
!            cnt = 0
!            do j = 1, MaxProfRep
!                do k = 1, MaxProfNodes
!                    if (E2Profile%gas4(j, k) /= error) then
!                        BiometVar%prof_gas4(k) = BiometVar%prof_gas4(k) + E2Profile%gas4(j, k)
!                        cnt(k) = cnt(k) + 1
!                    end if
!                end do
!            end do
!            where (cnt(:) /= 0)
!                BiometVar%prof_gas4(:) = BiometVar%prof_gas4(:) / cnt(:)
!            endwhere
!        else
!            BiometVar%prof_gas4(:) = error
!        end if
    else
        E2Biomet = ErrBiomet
        BiometVar = ErrBiometVar
    end if
end subroutine RetrieveBiometData
