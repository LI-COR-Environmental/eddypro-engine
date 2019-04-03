!***************************************************************************
! fit_rh_to_cutoff.f90
! --------------------
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

    !> Local variables
    integer :: RH
    integer :: m
    integer :: cls
    integer :: cls2
    integer :: cnt
    integer :: cnt2
    real(kind = dbl), allocatable  :: fvec(:), fjac(:,:)
    integer, parameter :: npar_EXP = 3
    integer :: info, ipvt_EXP(npar_EXP)
    real(kind = dbl) :: EXPPar(npar_EXP)
    real(kind = dbl) :: tol = 1d-04
    real(kind = dbl) :: mean_fc
    include '../src_common/interfaces.inc'


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
        write(*, '(a)', advance = 'no') ' Open-path H2O analyser: analytic &
            &fitting of cut-off frequencies vs. RH.. '
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
        write(*, '(a)', advance = 'no') ' Fitting in-situ assessment of &
            &cut-off frequencies vs. RH.. '

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
        if (m >= 4) then
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
        elseif(m == 3 .or. m == 2) then
            cnt = 0
            mean_fc = 0d0
            do cls = RH10, RH90
                if (RegPar(h2o, cls)%fc /= error) then
                    cnt = cnt + 1
                    mean_fc = mean_fc + RegPar(h2o, cls)%fc
                    if (cnt == m) then
                        mean_fc = mean_fc / cnt
                        RegPar(h2o, RH10:RH90)%fc = mean_fc
                        RegPar(dum, dum)%e1 = 1d-15
                        RegPar(dum, dum)%e2 = 1d-15
                        RegPar(dum, dum)%e3 = dlog(mean_fc)
                        exit
                    end if
                end if
            end do
        elseif(m == 1) then
            do cls = RH10, RH90
                if (RegPar(h2o, cls)%fc /= error) then
                    RegPar(h2o, RH10:RH90)%fc = RegPar(h2o, cls)%fc
                    RegPar(dum, dum)%e1 = 1d-15
                    RegPar(dum, dum)%e2 = 1d-15
                    RegPar(dum, dum)%e3 = dlog(RegPar(h2o, cls)%fc)
                    exit
                end if
            end do
        elseif(m == 0) then
            RegPar(dum, dum)%e1 = error
            RegPar(dum, dum)%e2 = error
            RegPar(dum, dum)%e3 = error
        endif

        !> Now that the regression is done, sets Fc which are at -9999, to the closest valid value.
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
end subroutine FitRh2Fco
