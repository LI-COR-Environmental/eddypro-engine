!***************************************************************************
! fit_models.f90
! --------------
!Copyright (C) 2011, LI-COR Biosciences
!
!This file is part of EddyPro (TM).
!
!EddyPro (TM) is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!EddyPro (TM) is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with EddyPro (TM).  If not, see <http://www.gnu.org/licenses/>.
!***************************************************************************

!***************************************************************************
! \file        src/fit_models.f90
! \brief       Fit TF model to actual ratio of gas to temperature spectra
! \version     3.0.0
! \date        2011-08-04
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FitTFModels(nbins)
    use m_fx_global_var
    use m_levenberg_marquardt
    implicit none

    !> interface to function fcn, defining the IIR and its jacobian
    interface
        subroutine fcn(m, npar, x, fvec, fjac, iflag)
            implicit none
            integer, parameter :: dbl   = kind(0.0d0)
            integer, intent(in)            :: m, npar
            real(kind = dbl), intent(in)    :: x(:)
            real(kind = dbl), intent(inout) :: fvec(:)
            real(kind = dbl), intent(out)   :: fjac(:,:)
            integer, intent(inout)         :: iflag
        end subroutine fcn
    end interface

    !> In/out variables
    integer, intent(in) :: nbins
    !> local variables
    integer :: nlong(GHGNumVar, MaxGasClasses)
    integer :: gas
    integer :: cls
    integer :: i
    integer :: m
    integer :: maxnlong
    real(kind = dbl), allocatable  :: fvec(:), fjac(:,:)
    integer, parameter :: npar_IIR = 2
    integer, parameter :: npar_sigma = 1
    integer :: info, ipvt_IIR(npar_IIR)  !, ipvt_sigma(npar_sigma)
    real(kind = dbl) :: IIRPar(npar_IIR),SigmaPar(npar_sigma)
    real(kind = dbl) :: tol = 1d-04
    type(LongSpectraType), allocatable :: lSpec(:, :)


    call log_msg(' inf=calculating cut-off frequencies.')
    write(*, '(a)', advance = 'no') ' Calculating cut-off frequencies..'

    !> Calculate length of un-binned spectra (lSpec), by looking at fnum for each bin
    call LongSpectraLength(nbins, maxnlong)

    !> Allocate lSpec based on detected nlong
    allocate(lSpec(maxnlong, MaxGasClasses))

    !> Create artificial dataset for regression
    call UnbinSpectra(lSpec, size(lSpec, 1), size(lSpec, 2), nlong, size(nlong, 1), size(nlong, 2), nbins)

    !> Allocate vectors for fit
    if (.not. allocated(xFit)) allocate(xFit(maxval(nlong)))
    if (.not. allocated(yFit)) allocate(yFit(maxval(nlong)))
    if (.not. allocated(zFit)) allocate(zFit(maxval(nlong)))
    if (.not. allocated(ddum)) allocate(ddum(maxval(nlong)))

    !> Substract high-frequencynoise if requested
    call SubtractHighFreqNoise(lSpec, size(lSpec, 1), size(lSpec, 2), nlong, size(nlong, 1), size(nlong, 2), nbins)

    !> regression for transfer functions, class-sorted
    do gas = co2, gas4
        do cls = 1, MaxGasClasses
            if (MeanBinSpecAvailable(cls, gas)) then
                !> Initialization of IIR and SIGMA filter parameters
                !> IIR  : Par(1) = Fn, Par(2)= f_co
                !> SIGMA: Par(1) = half-power frequency
                IIRPar(1) = 1d0
                IIRPar(2) = 0.02d0
                SigmaPar(1) = 0.5d0
                if (nlong(gas, cls) == 0) then
                    IIRPar(1:2) = error
                    SigmaPar(1) = error
                else
                    !> create dataset for regression, limiting natural frequencies
                    !> between fmin and fmax, user-defined for each gas
                    m = 0
                    do i = 1, nlong(gas, cls)
                        if (lSpec(i, cls)%fn(gas) > FCCsetup%SA%fmin(gas) &
                            .and. lSpec(i, cls)%fn(gas) < FCCsetup%SA%fmax(gas)) then
                            m = m + 1
                            !> makes regression using f*sp_x instead of sp_x
                            xFit(m) = lSpec(i, cls)%fn(gas)
                            yFit(m) = lSpec(i, cls)%of(gas) * lSpec(i, cls)%fn(gas)
                            zFit(m) = lSpec(i, cls)%ts(gas) * lSpec(i, cls)%fn(gas)
                        end if
                    end do

                    !> normalize to a reasonable factor, to avoid regression to
                    !> underflow or overflow
                    call NormalizeForRegression(m)

                    !> Compute regression parameters for IIR filter (Ibrom et al. 2007)
                    TFShape = 'iir'
                    allocate(fvec(m), fjac(m, npar_IIR))
!MAYBE WRONG        allocate(fvec(nlong(gas, cls)), fjac((nlong(gas, cls)), npar_IIR))
                    call lmder1(fcn, m, npar_IIR, IIRPar, fvec, fjac, tol, info, ipvt_IIR)
                    if ((IIRPar(1) == 1d0 .and. IIRPar(2) == 0.02d0) &
                        .or. info < 1 .or. info > 3) then
                        IIRPar(1:2) = error
                    end if
                    deallocate(fvec, fjac)

                    !> Compute regression parameters for sigma function (Aubinet et al. 2001)
!                    TFShape = 'sigma'
!                    allocate(fvec(nlong(gas, cls)), fjac((nlong(gas, cls)), npar_sigma))
!                    call lmder1(fcn, m, npar_sigma, SigmaPar, fvec, fjac, tol, info, ipvt_sigma)
!                    SigmaPar(1) = dabs(SigmaPar(1))
!                    if (SigmaPar(1) == 0.5d0 .or. info < 1 .or. info > 3) then
!                        SigmaPar(1) = error
!                    end if
!                    deallocate(fvec, fjac)
                end if
                !> Store regression params
                RegPar(gas, cls)%Fn = IIRPar(1)
                RegPar(gas, cls)%fc = IIRPar(2)
                RegPar(gas, cls)%f2 = SigmaPar(1)
            else
                RegPar(gas, cls)%Fn = error
                RegPar(gas, cls)%fc = error
                RegPar(gas, cls)%f2 = error
            end if
        end do
    end do

    if (allocated(xFit)) deallocate(xFit)
    if (allocated(yFit)) deallocate(yFit)
    if (allocated(zFit)) deallocate(zFit)
    if (allocated(ddum)) deallocate(ddum)

    write(*, '(a)') ' done.'
end subroutine FitTFModels

!***************************************************************************
! \file        src/fit_models.f90
! \brief       "Unbin" spectra, for the purposes of a weigthed regression
! \version     3.0.0
! \date        2011-08-04
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LongSpectraLength(nbins, maxnlong)
    use m_fx_global_var
    !> in/out variables
    integer, intent(in) :: nbins
    integer, intent(out) :: maxnlong
    !> local variables
    integer :: nlong(GHGNumVar, MaxGasClasses)
    integer :: gas
    integer :: cls
    integer :: bin


    nlong = 0
    do gas = co2, gas4
        do cls = 1, MaxGasClasses
            if (MeanBinSpecAvailable(cls, gas)) then
                !> If spectra are available, creates the long spectra
                do bin = 1, nbins
                    if (MeanBinSpec(bin, cls)%fnum(gas) > 0) &
                            nlong(gas, cls) = nlong(gas, cls) + MeanBinSpec(bin, cls)%fnum(gas)
                end do
            end if
        end do
    end do
    maxnlong = maxval(nlong)
end subroutine LongSpectraLength
!***************************************************************************
! \file        src/fit_models.f90
! \brief       "Unbin" spectra, for the purposes of a weigthed regression
! \version     3.0.0
! \date        2011-08-04
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine UnbinSpectra(lSpec, nrow, ncol, nlong, nnrow, nncol, nbins)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: nnrow, nncol
    integer, intent(in) :: nbins
    integer, intent(out) :: nlong(nnrow, nncol)
    type(LongSpectraType), intent(out) :: lSpec(nrow, ncol)
    !> local variables
    integer :: gas
    integer :: cls
    integer :: bin
    integer :: i
    integer :: err_cnt

    nlong = 0
    do gas = co2, gas4
        do cls = 1, nncol
            if (MeanBinSpecAvailable(cls, gas)) then
                !> Check if current-class spectra are available by checking if all rows of
                !> first column (fnum) are set to -9999
                err_cnt = 0
                do bin = 1, nbins
                    if (MeanBinSpec(bin, cls)%fnum(gas) == nint(error)) err_cnt = err_cnt + 1
                end do
                if (err_cnt == nbins) then
                    !> If current-class spectra are not available, set long spectra to -9999
                    lSpec(:, cls)%fn(gas) = error
                    lSpec(:, cls)%ts(gas) = error
                    lSpec(:, cls)%of(gas) = error
                else
                    !> If spectra are available, creates the long spectra
                    nlong(gas, cls) = 0
                    do bin = 1, nbins
                        if (MeanBinSpec(bin, cls)%fnum(gas) /= 0 .and. MeanBinSpec(bin, cls)%fnum(gas) /= nint(error)) then
                            do i = 1, MeanBinSpec(bin, cls)%fnum(gas)
                                nlong(gas, cls) = nlong(gas, cls) + 1
                                lSpec(nlong(gas, cls), cls)%fn(gas) = MeanBinSpec(bin, cls)%fn(gas)
                                lSpec(nlong(gas, cls), cls)%ts(gas) = MeanBinSpec(bin, cls)%ts(gas)
                                lSpec(nlong(gas, cls), cls)%of(gas) = MeanBinSpec(bin, cls)%of(gas)
                            end do
                        end if
                    end do
                end if
            end if
        end do
    end do
end subroutine UnbinSpectra

!***************************************************************************
! \file        src/fit_models.f90
! \brief       Take numbers to reasonable ranges prior to performing regressions
! \version     3.0.0
! \date        2011-08-04
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine NormalizeForRegression(N)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    !> local variables
    integer :: i
    real(kind=dbl) :: fit(N)
    real(kind=dbl) :: sq_diff(N)
    real(kind=dbl) :: sum_diff
    real(kind=dbl) :: gain

    !> Determine tentative sum of squared differences
    fit(1:N)     = zFit(1:N) / (1d0 + (xFit(1:N)/0.5d0)**2)
    sq_diff(1:N) = (fit(1:N) - yFit(1:N))**2
    sum_diff     = sum(sq_diff(1:N))

    !> detect the order of magnitude of sum_diff and derive the gain needed to
    !> take the sum to the order of O(10). This is given by the order of the sum (i),
    !> changed in sign and added with 1 (try it yourself..). Then this exponent is devided by
    !> 2 because it applies to the data, which in the sum appear as squared.
    gain = 1d0
    do i = -10, 10, 1
        if ((sum_diff/(10d0**i) > 1d0) .and. (sum_diff/(10d0**(i+1)) < 1d0)) then
            gain = 10d0 **((-i + 1)/2)
            exit
        end if
    end do

    !> Apply the gain to y and z
    yFit = yFit * gain
    zFit = zFit * gain
end subroutine NormalizeForRegression
