!***************************************************************************
! bpcf_aux_subs.f90
! -----------------
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
! \brief      Contains auxiliary subroutines for BPCF calculations.
! \author     Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************

!***************************************************************************
! \brief       Set all band-pass TF components to given value.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SetTransferFunctionsToValue(BPTF, nfreq, val)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(inout) :: nfreq
    type(BPTFType), intent(out) :: BPTF(nfreq)
    real(kind = dbl), intent(in) :: val
    !> local variables
    integer :: var

    do var = u, gas4
        BPTF(1:nfreq)%HP(var)  = val
        BPTF(1:nfreq)%EXP(var) = val
        BPTF(1:nfreq)%LP(var)  = LPTFType(val, val, val, val, val, val, val, val)
    end do
end subroutine SetTransferFunctionsToValue

!***************************************************************************
! \brief       Calculates total transfer function (band-pass transfer function)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BandPassTransferFunction(BPTF, var1, var2, varout, nfreq)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: var1
    integer, intent(in) :: var2
    integer, intent(in) :: varout
    integer, intent(in) :: nfreq
    type(BPTFType), intent(inout) :: BPTF(nfreq)

    BPTF%BP(varout) = BPTF%LP(var1)%dirga  * BPTF%LP(var2)%dirga  &
                    * BPTF%LP(var1)%dsonic * BPTF%LP(var2)%dsonic &
                    * BPTF%LP(var1)%wirga  * BPTF%LP(var2)%wirga  &
                    * BPTF%LP(var1)%wsonic * BPTF%LP(var2)%wsonic &
                    * BPTF%LP(var1)%sver   * BPTF%LP(var2)%sver   &
                    * BPTF%LP(var1)%shor   * BPTF%LP(var2)%shor   &
                    * BPTF%LP(var1)%t      * BPTF%LP(var2)%t      &
                    * BPTF%HP(var1)        * BPTF%HP(var2) &
                    * BPTF%EXP(var1)       * BPTF%EXP(var2)
end subroutine BandPassTransferFunction

!***************************************************************************
! \brief       Calculate spectral correction factors for all fluxes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SpectralCorrectionFactors(Cosp, var, nf, nfreq, BPTF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nfreq
    integer, intent(in) :: var
    real(kind = dbl), intent(in) :: nf(nfreq)
    real(kind = dbl), intent(in) :: Cosp(nfreq)
    type(BPTFType), intent(in) :: BPTF(nfreq)
    !> local variables
    integer :: k = 0! \file        src/bpcf_aux_subs.f90

    integer :: err_cnt = 0
    real(kind = dbl) :: IntCO
    real(kind = dbl) :: IntTFCO
    real(kind = dbl) :: nf_min, nf_max
    real(kind = dbl) :: df


    !> If cospectrum is made up of only error codes \n
    !> set correction factors to error as well
    err_cnt = 0
    do k = 1, nfreq
        if (Cosp(k) == error) err_cnt = err_cnt + 1
    end do
    if (err_cnt == nfreq) then
        BPCF%of(var) = error
        return
    end if

    !> Artificial frequency range, large enough to accomodate all cases
    nf_min = 1d0/5000d0
    nf_max = 100d0

    !> Integrals of cospectrum and filtered cospectrum
    IntCO = 0d0
    IntTFCO = 0d0
    do k = 1, nfreq - 1
        if (nf(k) > nf_min .and. nf(k + 1) < nf_max .and. &
            Cosp(k) /= error .and. BPTF(k)%BP(var) /= error) then
            df = nf(k + 1) - nf(k)
            IntCO = IntCO + Cosp(k) * df
            IntTFCO = IntTFCO + BPTF(k)%BP(var) * Cosp(k) * df
        end if
    end do
    if (IntTFCO /= 0d0) then
        BPCF%of(var) = IntCO / IntTFCO
    else
        BPCF%of(var) = error
    end if
end subroutine SpectralCorrectionFactors

!***************************************************************************
! \brief       Retrieve transfer function parameters, based on current month or RH
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveLPTFpars(lEx, tf_shape, LocSetup)
    use m_common_global_var
    implicit none
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    character(*), optional, intent(in) :: tf_shape
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    integer :: RH
    integer :: month
    real (kind = dbl) :: A
    real (kind = dbl) :: B
    real (kind = dbl) :: C
    real (kind = dbl) :: lRH

    f_c(co2:gas4) = error
    f_2(co2:gas4) = error

    select case (tf_shape)
        case('iir')
            !> calculate H2O cut-off frequency from current RH, using
            !> exponential fit parameters
            if (lEx%var_present(h2o)) then
                A = RegPar(dum, dum)%e1
                B = RegPar(dum, dum)%e2
                C = RegPar(dum, dum)%e3
                lRH = lEx%RH * 1d-2
                f_c(h2o) = dexp(A * lRH**2 + B * lRH + C)
            end if
            !> select relevant tranfer function parameters
            !> according to the month, for CO2, CH4, GAS4
            call char2int(lEx%date(6:7), month, 2)
            if(lEx%var_present(co2))  f_c(co2)  = RegPar(co2,  LocSetup%SA%class(co2,  month))%fc
            if(lEx%var_present(ch4))  f_c(ch4)  = RegPar(ch4,  LocSetup%SA%class(ch4,  month))%fc
            if(lEx%var_present(gas4)) f_c(gas4) = RegPar(gas4, LocSetup%SA%class(gas4, month))%fc

        case('sigma')
            !> select relevant tranfer function parameters
            !> according to the RH-class, for H2O
            if (lEx%var_present(h2o)) then
                do RH = RH10, RH90
                    if(lEx%RH > dfloat(RH)*10d0 - 5d0 &
                       .and. lEx%RH < dfloat(RH)*10d0 + 5d0) then
                        f_2(h2o) = RegPar(h2o, RH)%f2
                        exit
                    end if
                end do
            end if
            call char2int(lEx%date(6:7), month, 2)
            if(lEx%var_present(co2))  f_2(co2)  = RegPar(co2,  LocSetup%SA%class(co2,  month))%f2
            if(lEx%var_present(ch4))  f_2(ch4)  = RegPar(ch4,  LocSetup%SA%class(ch4,  month))%f2
            if(lEx%var_present(gas4)) f_2(gas4) = RegPar(gas4, LocSetup%SA%class(gas4,  month))%f2
    end select
end subroutine RetrieveLPTFpars

!***************************************************************************
! \brief       Calculates spectral correction factors based on the procedure \n
!              described in Horst, 1997, Boundary-Layer Meteorology 82: 219–233, 1997. \n
!              Based on  analytical cospectra and \n
!              analytical transfer functions, parameterized in-situ.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CorrectionFactorsHorst97(lBPCF, lEx)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(SpectralType), intent(inout) :: lBPCF
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    !> local variables
    real(kind = dbl) :: Nm
    real(kind = dbl) :: alpha
    real(kind = dbl) :: zeta
    real(kind = dbl) :: t_c

    !> Define peak frequency and alpha exponent, basing on stability, \n
    !> see Horst (1997, BLM), eq. 11. Note that his eq. 5 is equal to eq. 3 \n
    !> in Ibrom et al (2007, AFM), with 2*pi*t_c = 1/fc. Thus, the fc derived with \n
    !> the procedure described in Ibrom et al. 2007 can be used to calculate \n
    !> tau_c, and hence to derive the BPCF with Horst's method.
    if (lEx%zL <= 0d0) then
        Nm = 0.085d0
        alpha = 7d0 / 8d0
    else
        Nm = 2d0 - 1.915d0 / (1d0 + 0.5d0 * lEx%zL)
        alpha = 1d0
    end if
    zeta = lEx%instr(sonic)%height - lEx%disp_height

    !> Correction factor for co2
    if (lEx%var_present(co2)) then
        t_c = 1d0 / (2d0 * p * f_c(co2))
        lBPCF%of(w_co2) = (1d0 + 2d0 * p * lEx%WS * t_c * Nm / zeta )**alpha
    end if
    !> Correction factor for h2o
    if (lEx%var_present(h2o)) then
        t_c = 1d0 / (2d0 * p * f_c(h2o))
        lBPCF%of(w_h2o) = (1d0 + 2d0 * p * lEx%WS * t_c * Nm / zeta )**alpha
    end if
    !> Correction factor for ch4
    if (lEx%var_present(ch4)) then
        t_c = 1d0 / (2d0 * p * f_c(ch4))
        lBPCF%of(w_ch4) = (1d0 + 2d0 * p * lEx%WS * t_c * Nm / zeta )**alpha
    end if
    !> Correction factor for gas4
    if (lEx%var_present(gas4)) then
        t_c = 1d0 / (2d0 * p * f_c(gas4))
        lBPCF%of(w_gas4) = (1d0 + 2d0 * p * lEx%WS * t_c * Nm / zeta )**alpha
    end if
end subroutine CorrectionFactorsHorst97

!***************************************************************************
! \brief       Calculate correction factors, according to Ibrom et al (2007) \n
!              Agr. For. Meteorol., 147, 149-156. \n
!              Integral method. Derive cut-off frequencies from ensemble averages \n
!              of all available spectra, sorted for RH (H2) and months (other gases), eq. 6 \n
!              Derive low-pass correction factors from degraded temperature time-series, eq. 9. \n
!              analytical transfer functions, parameterized in-situ.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CorrectionFactorsIbrom07(do_co2, do_h2o, do_ch4, do_gas4, lBPCF, lEx)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(SpectralType), intent(inout) :: lBPCF
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    logical, intent(in) :: do_co2
    logical, intent(in) :: do_h2o
    logical, intent(in) :: do_ch4
    logical, intent(in) :: do_gas4


    if (lEx%zL >= 0d0) then
        if (lEx%var_present(co2) .and. do_co2) &
            lBPCF%of(w_co2) = StPar(1) * lEx%WS  / (StPar(2) + f_c(co2))  + 1d0
        if (lEx%var_present(h2o) .and. do_h2o) &
            lBPCF%of(w_h2o) = StPar(1) * lEx%WS  / (StPar(2) + f_c(h2o))  + 1d0
        if (lEx%var_present(ch4) .and. do_ch4) &
            lBPCF%of(w_ch4) = StPar(1) * lEx%WS  / (StPar(2) + f_c(ch4))  + 1d0
        if (lEx%var_present(gas4) .and. do_gas4) &
            lBPCF%of(w_gas4) = StPar(1) * lEx%WS / (StPar(2) + f_c(gas4)) + 1d0
    else
        if (lEx%var_present(co2) .and. do_co2) &
            lBPCF%of(w_co2) = UnPar(1) * lEx%WS  / (UnPar(2) + f_c(co2))  + 1d0
        if (lEx%var_present(h2o) .and. do_h2o) &
            lBPCF%of(w_h2o) = UnPar(1) * lEx%WS  / (UnPar(2) + f_c(h2o))  + 1d0
        if (lEx%var_present(ch4) .and. do_ch4) &
            lBPCF%of(w_ch4) = UnPar(1) * lEx%WS  / (UnPar(2) + f_c(ch4))  + 1d0
        if (lEx%var_present(gas4) .and. do_gas4) &
            lBPCF%of(w_gas4) = UnPar(1) * lEx%WS / (UnPar(2) + f_c(gas4)) + 1d0
    end if
end subroutine CorrectionFactorsIbrom07

!***************************************************************************
! \brief       Build up IIR and SIGMA transfer functions, \n
!              with file specific cut-offs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExperimentalLPTF(shape, nf, N, BPTF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: nf(N)
    character(*), intent(in) :: shape
    type(BPTFType), intent(out) :: BPTF(N)

    !> experimental transfer function, Fratini et al. 2012, Eq. 1 and 3
    if (shape == 'iir') then
        if (f_c(co2) /= error) then
            BPTF(1:N)%EXP(w_co2) = 1d0 / (1d0 + (nf(1:N) / f_c(co2))**2)
        else
            BPTF(:)%EXP(w_co2) = 1d0
        end if
        if (f_c(h2o) /= error) then
            BPTF(1:N)%EXP(w_h2o) = 1d0 / (1d0 + (nf(1:N) / f_c(h2o))**2)
        else
            BPTF(:)%EXP(w_h2o) = 1d0
        end if
        if (f_c(ch4) /= error) then
            BPTF(1:N)%EXP(w_ch4) = 1d0 / (1d0 + (nf(1:N) / f_c(ch4))**2)
        else
            BPTF(:)%EXP(w_ch4) = 1d0
        end if
        if (f_c(gas4) /= error) then
            BPTF(1:N)%EXP(w_gas4) = 1d0 / (1d0 + (nf(1:N) / f_c(gas4))**2)
        else
            BPTF(:)%EXP(w_gas4) = 1d0
        end if

    !> experimental transfer function, see Aubinet et al. (2001, AFM)
    elseif (shape == 'sigma') then
        if (f_2(co2) /= error) then
            BPTF(1:N)%EXP(w_co2) = dexp(-0.346574d0 * (nf(1:N) / f_2(co2))**2)
        else
            BPTF(:)%EXP(w_co2) = 1d0
        end if
        if (f_2(h2o) /= error) then
            BPTF(1:N)%EXP(w_h2o) = dexp(-0.346574d0 * (nf(1:N) / f_2(h2o))**2)
        else
            BPTF(:)%EXP(w_h2o) = 1d0
        end if
        if (f_2(ch4) /= error) then
            BPTF(1:N)%EXP(w_ch4) = dexp(-0.346574d0 * (nf(1:N) / f_2(ch4))**2)
        else
            BPTF(:)%EXP(w_ch4) = 1d0
        end if
        if (f_2(gas4) /= error) then
            BPTF(1:N)%EXP(w_gas4) = dexp(-0.346574d0 * (nf(1:N) / f_2(gas4))**2)
        else
            BPTF(:)%EXP(w_gas4) = 1d0
        end if
    end if
end subroutine ExperimentalLPTF
