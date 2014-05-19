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
subroutine ReadFullCosWT(CospFile, FullCospectrum, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(FilelistType), intent(in) :: CospFile
    integer, intent(out) :: N
    type(FullCospType), intent(out) :: FullCospectrum
    !> local variables
    integer :: io_status
    integer :: read_status
    integer :: nvar
    integer :: ord
    integer :: i
    real(kind = dbl), allocatable :: aux(:)
    character(1024) :: string
    character(32) :: var
    real(kind = dbl) :: cov_w_ts

    open(udf, file = CospFile%path, iostat = io_status)
    if (io_status /= 0) then
        FullCospectrum%fn(:) = error
        FullCospectrum%wt(:) = error
        call ErrorHandle(2, 0, 14)
        return
    end if

    !> Skip 10 header lines
    do i = 1, 10
        read(udf, *)
    end do

    !> Detects if H cospectrum is available in the file, if not exit with error code
    !> The control is done on the presence of the label "cov(w_ts)"
    ord = 0
    nvar = 0
    read(udf, '(a)') string
    do
        if (index(string, ',') /= 0) then
            var = string(1:index(string, ',') - 1)
            nvar = nvar + 1
            if (var == 'cov(w_ts)') ord = nvar
            string = string(index(string, ',') + 1:len_trim(string))
        else
            var = string(1:len_trim(string))
            nvar = nvar + 1
            if (var == 'cov(w_ts)') ord = nvar
            exit
        end if
    end do

    if (ord > 0) then
        !> Read w/ts covariance from first line
        if (.not. allocated(aux)) allocate (aux(nvar))
        read(udf, *, iostat = read_status) aux(1:nvar)
        cov_w_ts = aux(ord)

        !> Read natural frequencies and full co-spectrum
        if (cov_w_ts /= error .and. cov_w_ts /= 0d0) then
            !> Skip one line
            read(udf, *)
            N = 0
            do
                N = N + 1
                read(udf, *, iostat = read_status) aux(1:nvar)
                if (read_status > 0) then
                    N = N - 1
                    cycle
                end if
                if (read_status < 0) exit
                FullCospectrum%fn(N) = aux(1)
                FullCospectrum%wt(N) = aux(ord)
            end do
            N = N - 1
            FullCospectrum%fn(N + 1: MaxNumRow) = error

            !> Un-normalize cospectrum
            where (FullCospectrum%fn(1:N) /= 0d0 .and. FullCospectrum%fn(1:N) /= error)
                FullCospectrum%wt(1:N) = &
                    FullCospectrum%wt(1:N) /  FullCospectrum%fn(1:N) * cov_w_ts
            elsewhere
                FullCospectrum%wt(1:N) = error
                FullCospectrum%fn(1:N) = error
            end where
        else
            N = nint(error)
        end if
    else
        N = nint(error)
    end if
    close(udf)
    if (allocated(aux)) deallocate (aux)


    !> Check that cospectrum values are within reasonable values (not too high).
    !> Discard spectrum if this is the case
    !> This is somewhat arbitrary, introduced to eliminate observed implausible cospectra
    !> It's very strict: one only outranged value will eliminate the whole cospectrum
    do i = 1, N
        if (dabs(FullCospectrum%wt(i)) > MaxSpecValue) then
            N = nint(error)
            exit
        end if
    end do
end subroutine ReadFullCosWT
