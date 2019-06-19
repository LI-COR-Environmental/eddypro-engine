!***************************************************************************
! fit_cospectral_models.f90
! -------------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2019, LI-COR Biosciences, Inc.  All Rights Reserved.
! Author: Gerardo Fratini
!
! This file is part of EddyPro®.
!
! NON-COMMERCIAL RESEARCH PURPOSES ONLY - EDDYPRO® is licensed for 
! non-commercial academic and government research purposes only, 
! as provided in the EDDYPRO® End User License Agreement. 
! EDDYPRO® may only be used as provided in the End User License Agreement
! and may not be used or accessed for any commercial purposes.
! You may view a copy of the End User License Agreement in the file
! EULA_NON_COMMERCIAL.rtf.
!
! Commercial companies that are LI-COR flux system customers 
! are encouraged to contact LI-COR directly for our commercial 
! EDDYPRO® End User License Agreement.
!
! EDDYPRO® contains Open Source Components (as defined in the 
! End User License Agreement). The licenses and/or notices for the 
! Open Source Components can be found in the file LIBRARIES-ENGINE.txt.
!
! EddyPro® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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
    include '../src_common/interfaces.inc'


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
            ' Fitting model cospectra to unstable ' // &
                cvar(1:len_trim(cvar)) // ' cospectra..'

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
        write(*,'(a)') ' Done'
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
        write(*,'(a)') ' Done'
    end do

    if (allocated(xFit)) deallocate(xFit)
    if (allocated(yFit)) deallocate(yFit)
    if (allocated(ddum)) deallocate(ddum)
end subroutine FitCospectralModel
