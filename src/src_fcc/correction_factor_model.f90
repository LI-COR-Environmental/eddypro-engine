!***************************************************************************
! correction_factor_model.f90
! ---------------------------
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
! \brief       Apply non linear regression to calculate correction factor \n
!              model parameters c1 and c2, after Ibrom et al. (2007, AFM) eq. 9 \n
!              using progressively degraded temperature time-series
! \author      Gerardo Fratini
! \sa
! \bug
! \deprecated
! \test
! \todo        Allow user to select F_low individual values to be either monotonic \n
!              (hard constrain), > 1 (soft constrain) or any value (no constrains) \n
!              So far, only "monotonic" option in available, but the code is present \n
!              and commented for the other two options.
!***************************************************************************
subroutine CorrectionFactorModel(ExFilename, NumExRecords)
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
    character(*), intent(in) :: ExFilename
    integer, intent(in) :: NumExRecords
    !> local variables
    integer :: i
    integer :: j
    integer :: m
    integer :: NumDegtRecords
    integer, parameter :: npar = 2
    integer :: info, ipvt(npar)
    real(kind = dbl), allocatable :: DegWT(:, :)
    real(kind = dbl), allocatable :: Fl(:, :)
    real(kind = dbl) :: tol = 1d-04
    real(kind = dbl), allocatable  :: fvec(:), fjac(:,:)
    real(kind = dbl) :: f_co(Nt)
    data f_co(1:Nt) /1.626d0, 0.614d0, 0.277d0, 0.133d0, 6.5d-2, 3.2d-2, 1.6d-2, 8d-3, 4d-3/


    write(*, '(a)') ' Fitting low-pass correction factor model as from Ibrom et al. (2007)..'
    allocate(DegWT(NumExRecords, Nt + 3))

    !> Retrieve "degraded" w'T' covariances from essentials file
    call ExtractColumnFromEssentials(ExFilename, NumExRecords, 'degraded_wT_covariances', &
        DegWT, size(DegWT, 1), size(DegWT, 2), NumDegtRecords)
    allocate(Fl(NumDegtRecords, Nt))

    if (.not. allocated(xFit)) allocate(xFit(NumDegtRecords * Nt))
    if (.not. allocated(yFit)) allocate(yFit(NumDegtRecords * Nt))
    if (.not. allocated(zFit)) allocate(zFit(NumDegtRecords * Nt))
    if (.not. allocated(zzFit)) allocate(zzFit(NumDegtRecords * Nt))
    if (.not. allocated(ddum)) allocate(ddum(NumDegtRecords * Nt))

    !> Perform a rough despiking of degraded covariances, to avoind unrealistic results
    do j = 1, size(DegWT, 2)
        call QuickDespiking(DegWT(1:NumDegtRecords, j), NumDegtRecords)
    end do

    !> calculate fc-dependent correction factors
    do j = 1, Nt
        where (DegWT(1: NumDegtRecords, j) /= 0d0 &
            .and. DegWT(1: NumDegtRecords, j) /= error .and. DegWT(1: NumDegtRecords, Nt + 1) /= error)
            Fl(1: NumDegtRecords, j) = DegWT(1: NumDegtRecords, Nt + 1) / DegWT(1: NumDegtRecords, j)
        elsewhere
            Fl(1: NumDegtRecords, j) = error
        end where
    end do

    TFShape = 'hyperbole'

    !> Unstable case
    !> Initialize regression parameters
    !> Par(1) = c1, Par(2) = c2 in eq. 9 in Ibrom et al. (2009, AFM)
    UnPar(1) = 2.5d-3
    UnPar(2) = 6.0d-4
    !> Create dataset for regression
    m = 0
    un_loop: do i = 1, NumDegtRecords
        if (DegWT(i, Nt + 1) /= error .and. DegWT(i, 1) /= error .and. DegWT(i, Nt) /= error &
            .and. DegWT(i, Nt + 2) /= error .and. DegWT(i, Nt + 3) /= error .and. DegWT(i, Nt + 3) < 0d0) then
            !> Code for monotonic individual F_low
            do j = 1, Nt
                if (j == 1) then
                    if (Fl(i, j) < 1d0) cycle un_loop
                else
                    if (Fl(i, j) < Fl(i, j - 1)) cycle un_loop
                end if
            end do

!            !> For unconstrained option: eliminate this filter
!            !> For individual F_low > 1
!            do j = 1, Nt
!                    if (degT(i)%Fl(j) < 1d0) cycle un_loop
!            end do

            !> Create dataset for regression
            do j = 1, Nt
                m = m + 1
                xFit(m)  = f_co(j)
                yFit(m)  = DegWT(i, Nt + 1)
                zzFit(m) = DegWT(i, j)
                zFit(m)  = DegWT(i, Nt + 2)
            end do
        end if
    end do un_loop

    !> Compute regression parameters for hyperbole function (Ibrom et al. 2007, eq. 9)
    allocate(fvec(m), fjac(m, npar))
    call lmder1(fcn, m, npar, UnPar, fvec, fjac, tol, info, ipvt)
    if ((UnPar(1) == 2.5d-3 .and. UnPar(2) == 6.0d-4) &
        .or. info < 1 .or. info > 3) then
        UnPar(1:2) = error
    end if
    deallocate(fvec, fjac)

    !> Stable case
    !> Initialize regression parameters
    !> Par(1) = c1, Par(2) = c2 in eq. 9 in Ibrom et al. (2009, AFM)
    StPar(1) = 6.5d-3
    StPar(2) = 3.0d-3
    !> Create dataset for regression
    m = 0
    st_loop: do i = 1, NumDegtRecords
        if (DegWT(i, Nt + 1) /= error .and. DegWT(i, 1) /= error .and. DegWT(i, Nt) /= error &
            .and. DegWT(i, Nt + 2) /= error .and. DegWT(i, Nt + 3) /= error .and. DegWT(i, Nt + 3) > 0d0) then
            !> Code for monotonic individual F_low
            do j = 1, Nt
                if (j == 1) then
                    if (Fl(i, j) < 1d0) cycle st_loop
                else
                    if (Fl(i, j) < Fl(i, j - 1)) cycle st_loop
                end if
            end do

!            !> For unconstrained option: eliminate this filter
!            !> For individual F_low > 1
!            do j = 1, Nt
!                    if (degT(i)%Fl(j) < 1d0) cycle st_loop
!            end do

            !> Create dataset for regression
            do j = 1, Nt
                m = m + 1
                xFit(m)  = f_co(j)
                yFit(m)  = DegWT(i, Nt + 1)
                zzFit(m) = DegWT(i, j)
                zFit(m)  = DegWT(i, Nt + 2)
            end do
        end if
    end do st_loop

    !> Compute regression parameters for hyperbole function (Ibrom et al. 2007, eq. 9)
    allocate(fvec(m), fjac(m, npar))
    call lmder1(fcn, m, npar, StPar, fvec, fjac, tol, info, ipvt)
    if ((StPar(1) == 1d0 .and. StPar(2) == 0.02d0) &
        .or. info < 1 .or. info > 3) then
        StPar(1:2) = error
    end if
    deallocate(fvec, fjac)

    if (allocated (DegWT)) deallocate (DegWT)
end subroutine CorrectionFactorModel

!***************************************************************************
!
! \brief       Despike dataset before regression, to limit the chance of \n
!              unrealistic results
! \author      Gerardo Fratini
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine QuickDespiking(vec, nrow)
    use m_fx_global_var
    implicit none
	!> in/out variables
    integer, intent(in) :: nrow
    real(kind = dbl), intent(inout) :: vec(nrow)
	!> local variables
    integer :: i
    real(kind = dbl) :: array(nrow, 1)
    real(kind = dbl) :: Mean(1)
    real(kind = dbl) :: StDev(1)


    !> Calculate mean and standard deviation
    array(:, 1) = vec(:)
    call AverageNoError(array, size(array, 1), size(array, 2), Mean, error)
    call CovarianceMatrixNoError(array, size(array, 1), size(array, 2), StDev, error)
    if (StDev(1) >= 0) then
        StDev(1) = dsqrt(StDev(1))
    else
        StDev(1) = error
    end if
    !> Now eliminate spikes, substituting with error codes
    do i = 1, nrow
        if (vec(i) <= Mean(1) - 3.5d0 * StDev(1) .or. vec(i) >= Mean(1) + 3.5d0 * StDev(1)) vec(i) = error
    end do
end subroutine QuickDespiking

