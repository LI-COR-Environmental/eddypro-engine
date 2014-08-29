!***************************************************************************
! bpcf_read_full_cosp_wt.f90
! --------------------------
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
    character(1024) :: string
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

    !> Detects cospectra available in the file, if none of those wanted, exit with error code
    !> The control is done on the presence of the label "cov(w_xx)"
    ord = 0
    nvar = 0
    read(udf, '(a)') string
    do
        if (index(string, ',') /= 0) then
            var = string(1:index(string, ',') - 1)
            nvar = nvar + 1
            do j = 1, GHGNumVar
                if (var == covlabs(j)) ord(j) = nvar
            end do
            string = string(index(string, ',') + 1:len_trim(string))
        else
            var = string(1:len_trim(string))
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

        !> Check that cospectra values are within reasonable values (not too high).
        !> Discard co-spectra if this is the case
        !> This is somewhat arbitrary, introduced to eliminate observed implausible cospectra
        !> It's very strict: one only outranged value will eliminate the whole cospectra set
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
    character(1024) :: string
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
        read(udf, '(a)', iostat = io_status) string
        if (io_status > 0) cycle
        if (io_status < 0) exit
        N = N + 1
    end do
end subroutine FullCospectraLength
