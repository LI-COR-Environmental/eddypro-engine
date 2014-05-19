!***************************************************************************
! fit_cospectral_model.f90
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
! \brief       Fit cospectral models (as from Runkle et al. 2012) to
!              stability sorted ensemble cospectra
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FitCospectralModel(nfit, dim1, dim2, FitStable, FitUnstable, fnrow)
    use m_fx_global_var
    use m_levenberg_marquardt
    implicit none

    !> interface to function fcn
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

    !> in/out variables
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    integer, intent(in) :: fnrow
    integer, intent(in) :: nfit(dim1, dim2)
    type(FitSpectraType), intent(in) :: FitUnstable(fnrow)
    type(FitSpectraType), intent(in) :: FitStable(fnrow)

    !> local variables
    integer, parameter :: npar = 3
    real(kind = dbl) :: InitPar(3)
    data InitPar(1) /1.d0   / &
         InitPar(2) /0.09d0 / &
         InitPar(3) /0.2d0  /

    integer :: m
    integer :: var
    integer :: info, ipvt_mass(npar)
    real(kind = dbl) :: lMassPar(npar)
    real(kind = dbl) :: tol = 1d-06
    real(kind = dbl), allocatable  :: fvec(:), fjac(:,:)
    character(32) :: cvar


    if (.not. allocated(xFit)) allocate (xFit(fnrow))
    if (.not. allocated(yFit)) allocate (yFit(fnrow))

    TFShape = 'cospectra_massman'
    !> Unstable case
    MassPar = MassParType(error, error, error)
    do var = w_ts, w_gas4
        xFit = 0d0
        yFit = 0d0
        m = nfit(var, unstable)
        if (m < 200) cycle !< includes variables that are not present
        select case (var)
            case (w_ts)
                cvar = 'w/T'
            case (w_co2)
                cvar = 'w/CO2'
            case (w_h2o)
                cvar = 'w/H2O'
            case (w_ch4)
                cvar = 'w/CH4'
            case (w_gas4)
                cvar = 'w/' // g4lab(1:g4l)
        end select

        write(*, '(a)', advance = 'no') &
            ' Fitting model cospectra to unstable ' // cvar(1:len_trim(cvar)) // ' cospectra..'

        !> create dataset for regression
        xFit(1:m) = FitUnstable(1:m)%fnorm(var)
        yFit(1:m) = dlog(FitUnstable(1:m)%of(var))
        !yFit(1:m) = FitUnstable(1:m)%of(var)

        !> Initialization of model parameter
        lMassPar = InitPar

        !> Regression parameters for model cospectrum, after Massman...
        allocate(fvec(m), fjac(m, npar))
        call lmder1(fcn, m, npar, lMassPar, fvec, fjac, tol, info, ipvt_mass)
        if ( (lMassPar(1) == InitPar(1) .and. lMassPar(2) == InitPar(2) .and. lMassPar(3) == InitPar(3)) &
            .or. info < 1 .or. info > 3) then
            MassPar(var, unstable) = MassParType(error, error, error)
        else
            MassPar(var, unstable)%a0    = lMassPar(1)
            MassPar(var, unstable)%fpeak = lMassPar(2)
            MassPar(var, unstable)%mu    = lMassPar(3)
        end if
        deallocate(fvec, fjac)
        write(*,'(a)') ' done'
    end do

    !> Stable case
    do var = w_ts, w_gas4
        xFit = 0d0
        yFit = 0d0
        m = nfit(var, stable)
        if (m < 100) cycle !< includes variables that are not present
        select case (var)
            case (w_ts)
                cvar = 'w/T'
            case (w_co2)
                cvar = 'w/CO2'
            case (w_h2o)
                cvar = 'w/H2O'
            case (w_ch4)
                cvar = 'w/CH4'
            case (w_gas4)
                cvar = 'w/' // g4lab(1:g4l)
        end select
        write(*, '(a)', advance = 'no') &
            ' Fitting model cospectra to stable ' // cvar(1:len_trim(cvar)) // ' cospectra..'

        !> create dataset for regression
        xFit(1:m) = FitStable(1:m)%fnorm(var)
        yFit(1:m) = dlog(FitStable(1:m)%of(var))

        !> Initialization of model parameter
        lMassPar = InitPar

        !> Regression parameters for model cospectrum, after Massman...
        allocate(fvec(m), fjac(m, npar))
        call lmder1(fcn, m, npar, lMassPar, fvec, fjac, tol, info, ipvt_mass)
        if ( (lMassPar(1) == InitPar(1) .and. lMassPar(2) == InitPar(2) .and. lMassPar(3) == InitPar(3)) &
            .or. info < 1 .or. info > 3) then
            MassPar(var, stable) = MassParType(error, error, error)
        else
            MassPar(var, stable)%a0    = lMassPar(1)
            MassPar(var, stable)%fpeak = lMassPar(2)
            MassPar(var, stable)%mu    = lMassPar(3)
        end if
        deallocate(fvec, fjac)
        write(*,'(a)') ' done'
    end do

    if (allocated(xFit)) deallocate(xFit)
    if (allocated(yFit)) deallocate(yFit)
    if (allocated(ddum)) deallocate(ddum)
end subroutine FitCospectralModel
