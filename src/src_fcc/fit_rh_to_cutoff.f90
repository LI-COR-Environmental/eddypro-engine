!***************************************************************************
! fit_rh_to_cutoff.f90
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
! \brief       Fit exponential curve to actual Fco vs. RH \n
!              after Ibrom et al. (2007, AFM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FitRh2Fco()
    use m_fx_global_var
    use m_levenberg_marquardt
    implicit none

    !> Interface to function fcn, defining the model and its jacobian
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

    !> Local variables
    integer :: RH
    integer :: m
    integer :: cls
    integer :: cls2
    integer :: cnt
    integer :: cnt2
    real(kind = dbl) :: avrg_fc = 0d0
    real(kind = dbl), allocatable  :: fvec(:), fjac(:,:)
    integer, parameter :: npar_EXP = 3
    integer :: info, ipvt_EXP(npar_EXP)
    real(kind = dbl) :: EXPPar(npar_EXP)
    real(kind = dbl) :: tol = 1d-04
    real(kind = dbl) :: mean_fc


    !> Allocate arrays for fits
    if (.not. allocated(xFit)) allocate(xFit(10))
    if (.not. allocated(yFit)) allocate(yFit(10))
    if (.not. allocated(ddum)) allocate(ddum(10))

    !> Preliminary validation of calculated cut-offs for water vapour
    where (RegPar(h2o, RH10:RH90)%fc > FCCMetadata%ac_freq / 2d0 .or. &
        RegPar(h2o, RH10:RH90)%fc < 0d0) &
        RegPar(h2o, RH10:RH90)%fc = error

    if (FCCMetadata%H2oPathType == 'open') then
        !> If the instrument associated to the first H2O reading is an open path
        !> fit exponential model analytically, such that the exponential function
        !> provides a constant value, equal to the mean value of f_cutoff among all
        !> RH classes
        call log_msg(' inf=fitting cut-off frequencies versus RH')
        write(*, '(a)', advance = 'no') ' Open-path H2O analyser: analytic fitting of cut-off frequencies vs. RH.. '
        mean_fc = 0d0
        cnt2 = 0
        do RH = RH10, RH90
            if (RegPar(h2o, RH)%fc /= error) then
                mean_fc = mean_fc + RegPar(h2o, RH)%fc
                cnt2 = cnt2 + 1
            end if
        end do
        mean_fc = mean_fc / cnt2
        RegPar(dum, dum)%e1 = 1d-15
        RegPar(dum, dum)%e2 = 1d-15
        RegPar(dum, dum)%e3 = dlog(mean_fc)
    else
        !> Fit exponential model by least squares minimization
        call log_msg(' inf=fitting cut-off frequencies versus RH')
        write(*, '(a)', advance = 'no') ' Fitting in-situ assessment of cut-off frequencies vs. RH.. '

        !> Initialization of function parameters (see Ibrom et al. 2007, AFM)
        !> EXP: Par(1) = A, Par(2)= B, Par(3)=C
        !> in function f(x) = exp**(A * x**2 + B * x + C)
        EXPPar(1) = -2.d0
        EXPPar(2) = -1.d0
        EXPPar(3) = -2.d0
        m = 0
        do cls = RH10, RH90
            if (RegPar(h2o, cls)%fc /= error) then
                m = m + 1
                xFit(m) = dfloat(cls) * 1d-1
                yFit(m) = RegPar(h2o, cls)%fc
            end if
        end do

        !> Extrapolate cut-off frequency with different policies, depending on how
        !> many have been correctly calculated
        if (m >= 3) then
            !> Perform regression if there are at least 3 RH/fc pairs
            TFShape = 'exponential'
            allocate(fvec(m), fjac(m, npar_EXP))
            call lmder1(fcn, m, npar_EXP, EXPPar, fvec, fjac, tol, info, ipvt_EXP)
            if ((EXPPar(1) == -2d0 .and. EXPPar(2) == -1.d0 .and. EXPPar(3) == -2.d0) &
                .or. info < 1 .or. info > 3) then
                EXPPar(1:3) = error
            end if
            deallocate(fvec, fjac)
            RegPar(dum, dum)%e1 = EXPPar(1)
            RegPar(dum, dum)%e2 = EXPPar(2)
            RegPar(dum, dum)%e3 = EXPPar(3)
        elseif(m == 2) then
            cnt = 0
            avrg_fc = 0d0
            do cls = RH10, RH90
                if (RegPar(h2o, cls)%fc /= error) then
                    cnt = cnt + 1
                    avrg_fc = avrg_fc + RegPar(h2o, cls)%fc
                    if (cnt == 2) then
                        avrg_fc = avrg_fc / cnt
                        RegPar(h2o, RH10:RH90)%fc = avrg_fc
                        exit
                    end if
                end if
            end do
        elseif(m == 1) then
            do cls = RH10, RH90
                if (RegPar(h2o, cls)%fc /= error) then
                    RegPar(h2o, RH10:RH90)%fc = RegPar(h2o, cls)%fc
                    exit
                end if
            end do
        elseif(m == 0) then
            RegPar(dum, dum)%e1 = error
            RegPar(dum, dum)%e2 = error
            RegPar(dum, dum)%e3 = error
        endif

        !> Now that the regession is done, sets Fc which are at -9999, to the closest valid value.
        if (m >= 3) then
            !> For low RH classes, uses higher ones
            do cls = RH50, RH10, - 1
                if (RegPar(h2o, cls)%fc == error) then
                    do cls2 = cls + 1, RH90
                        if (RegPar(h2o, cls2)%fc /= error) then
                            RegPar(h2o, cls)%fc = RegPar(h2o, cls2)%fc
                            RegPar(h2o, cls)%Fn = RegPar(h2o, cls2)%Fn
                            exit
                        end if
                    enddo
                end if
            end do
            !> For high RH classes, uses lower ones
            do cls = RH60, RH90
                if (RegPar(h2o, cls)%fc == error) then
                    do cls2 = cls - 1, RH10, - 1
                        if (RegPar(h2o, cls2)%fc /= error) then
                            RegPar(h2o, cls)%fc = RegPar(h2o, cls2)%fc
                            RegPar(h2o, cls)%Fn = RegPar(h2o, cls2)%Fn
                            exit
                        end if
                    enddo
                end if
            end do
        end if
    end if

    if (allocated(xFit)) deallocate(xFit)
    if (allocated(yFit)) deallocate(yFit)
    if (allocated(ddum)) deallocate(ddum)

    write(*, '(a)') 'Done.'
    call log_msg(' inf=fit calculated correctly.')
end subroutine FitRh2Fco
