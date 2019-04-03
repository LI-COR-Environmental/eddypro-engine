!***************************************************************************
! bpcf_read_full_cos_wt.f90
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
! \brief       Read file with full cospectrum wT, for spectral correction
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ImportFullCospectra(CospFile, cospectra, nfreq, wanted, skip)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(FilelistType), intent(in) :: CospFile
    logical, intent(in) :: wanted(GHGNumVar)
    integer, intent(in) :: nfreq
    logical, intent(out) :: skip
    type(SpectraSetType), intent(out) :: cospectra(nfreq)
    !> local variables
    integer :: io_status
    integer :: read_status
    integer :: nvar
    integer :: ord(GHGNumVar)
    integer :: i
    integer :: j
    real(kind = dbl), allocatable :: aux(:)
    character(ShortInstringLen) :: dataline
    character(32) :: var
    character(11) :: covlabs(GHGNumVar)
    real(kind = dbl) :: cov(GHGNumVar)
    data covlabs / 'cov(w_u)', 'cov(w_v)', 'cov(w_w)', 'cov(w_ts)', &
                   'cov(w_co2)', 'cov(w_h2o)', 'cov(w_ch4)', 'cov(w_gas4)' /


    skip = .false.

    open(udf, file = CospFile%path, iostat = io_status)
    if (io_status /= 0) then
        cospectra = ErrSpec
        skip = .true.
        call ExceptionHandler(63)
        return
    end if

    !> Skip 10 header lines
    do i = 1, 10
        read(udf, *)
    end do

    !> Detects cospectra available in the file, if none of those wanted,
    !> exit with error code. The control is done on the presence
    !> of the label "cov(w_xx)"
    ord = 0
    nvar = 0
    read(udf, '(a)') dataline
    do
        if (index(dataline, ',') /= 0) then
            var = dataline(1:index(dataline, ',') - 1)
            nvar = nvar + 1
            do j = 1, GHGNumVar
                if (var == covlabs(j)) ord(j) = nvar
            end do
            dataline = dataline(index(dataline, ',') + 1:len_trim(dataline))
        else
            var = dataline(1:len_trim(dataline))
            nvar = nvar + 1
            do j = 1, GHGNumVar
                if (var == covlabs(j)) ord(j) = nvar
            end do
            exit
        end if
    end do

    if (nvar == 0) then
        cospectra = ErrSpec
        skip = .true.
        call ExceptionHandler(63)
        return
    end if

    if (any(ord>0 .and. wanted)) then
        !> Import covariances
        if (.not. allocated(aux)) allocate (aux(nvar))
        aux = 0d0
        cov = error
        read(udf, *, iostat = read_status) aux(1:nvar)
        where (ord > 0) cov = aux(ord)
        !> Skip one line
        read(udf, *)
        !> Read natural frequencies and full co-spectra
        i = 0
        do
            read(udf, *, iostat = read_status) aux(1:nvar)
            if (read_status > 0) cycle
            if (read_status < 0) exit
            i = i + 1
            if (i > nfreq) exit
            cospectra(i)%fn = aux(1)
            where (ord > 0) cospectra(i)%of(:) = aux(ord)
        end do
        close(udf)
        if (allocated(aux)) deallocate (aux)

        !> Un-normalize cospectra
        do j = w_u, w_gas4
            where (wanted(j) .and. cospectra(:)%fn /= 0d0 .and. cospectra(:)%fn /= error)
                cospectra(:)%of(j) = cospectra(:)%of(j) /  cospectra(:)%fn * cov(j)
            elsewhere
                cospectra(:)%of(j) = error
            end where
        end do

        !> Check that cospectra values are within reasonable values
        !> (not too high). Discard co-spectra if this is the case
        !> This is somewhat arbitrary, introduced to eliminate observed
        !> implausible cospectra It's very strict: one only outranged
        !> value will eliminate the whole cospectra set
        do j = w_u, w_gas4
            if (wanted(j) .and. any(dabs(cospectra(:)%of(j)) > MaxSpecValue)) then
                cospectra = ErrSpec
                skip = .true.
                exit
            end if
        end do
    else
        skip = .true.
    end if

end subroutine ImportFullCospectra

!***************************************************************************
! \brief       read file with full cospectrum wT, for spectral correction
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FullCospectraLength(Filepath, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: Filepath
    integer, intent(out) :: N
    !> local variables
    integer :: io_status
    integer :: i
    character(ShortInstringLen) :: dataline
    character(11) :: covlabs(GHGNumVar)
    data covlabs / 'cov(w_u)', 'cov(w_v)', 'cov(w_w)', 'cov(w_ts)', &
                   'cov(w_co2)', 'cov(w_h2o)', 'cov(w_ch4)', 'cov(w_gas4)' /


    open(udf, file = Filepath, iostat = io_status)
    if (io_status /= 0) then
        N = nint(error)
        return
    end if

    !> Skip 10 header lines
    do i = 1, 13
        read(udf, *, iostat = io_status)
        if (io_status /= 0) then
            N = nint(error)
            return
        end if
    end do

    N = 0
    do
        read(udf, '(a)', iostat = io_status) dataline
        if (io_status > 0) cycle
        if (io_status < 0) exit
        N = N + 1
    end do
end subroutine FullCospectraLength
