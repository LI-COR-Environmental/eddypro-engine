!***************************************************************************
! read_binned_file.f90
! --------------------
! Copyright (C) 2007-2011, Eco2s team, Gerardo Fratini
! Copyright (C) 2011-2015, LI-COR Biosciences
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
! \brief       Reads file containing binned co-spectra and import relevant ones
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadBinnedFile(InFile, BinSpec, BinCosp, nrow, nbins, skip)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(FileListType), intent(in) :: InFile
    type(SpectraSetType), intent(out) :: BinSpec(nrow)
    type(SpectraSetType), intent(out) :: BinCosp(nrow)
    logical, intent(out) :: skip
    integer, intent(out) :: nbins
    !> local variables
    integer :: i
    integer :: var
    integer :: open_status
    integer :: read_status
    real(kind = dbl) :: aux


    !> Open file
    open(udf, file = InFile%path, iostat = open_status)
    !> Control on error in file opening
    skip = .false.
    if (open_status /= 0) then
        skip = .true.
        call ExceptionHandler(62)
        return
    end if

    !> Read frequency and (co)spectra of w, Ts, co2, h2o, ch4, gas4
    BinSpec = ErrSpec
    BinCosp = ErrSpec
    i = 0
    do
        i = i + 1
        if (i > nrow) exit
        read(udf, *, iostat = read_status) BinSpec(i)%fnum, BinSpec(i)%fn, BinSpec(i)%fnorm, aux, aux, BinSpec(i)%of(w),  &
            BinSpec(i)%of(ts), BinSpec(i)%of(co2), BinSpec(i)%of(h2o), BinSpec(i)%of(ch4), BinSpec(i)%of(gas4), &
            aux, aux, BinCosp(i)%of(w_ts), BinCosp(i)%of(w_co2), BinCosp(i)%of(w_h2o), &
            BinCosp(i)%of(w_ch4), BinCosp(i)%of(w_gas4)
        if (read_status > 0) then
            i = i - 1
            cycle
        end if
        if (read_status < 0) exit
        BinCosp(i)%fnum = BinSpec(i)%fnum
        BinCosp(i)%fn   = BinSpec(i)%fn
        BinCosp(i)%fnorm = BinSpec(i)%fnorm
    end do
    nbins = i - 1

    !> Un-normalize binned spectra by dividing by the frequency
    do var = w, gas4
        where (BinSpec(1:nbins)%fn /= error &
            .and. BinSpec(1:nbins)%fn /= 0d0 .and. BinSpec(1:nbins)%of(var) /= error)
            BinSpec(1:nbins)%of(var) = BinSpec(1:nbins)%of(var) / BinSpec(1:nbins)%fn
        else where
            BinSpec(1:nbins)%of(var) = error
        end where
    end do

    !> Un-normalize binned cospectra by dividing by the frequency
    do var = w_ts, w_gas4
        where (BinCosp(1:nbins)%fn /= error &
            .and. BinCosp(1:nbins)%fn /= 0d0 .and. BinCosp(1:nbins)%of(var) /= error)
            BinCosp(1:nbins)%of(var) = BinCosp(1:nbins)%of(var) / BinCosp(1:nbins)%fn
        else where
            BinCosp(1:nbins)%of(var) = error
        end where
    end do

    !> Check that spectra values are within reasonable values (not too high).
    !> Individually discard spectra if this is not the case
    !> This is somewhat arbitrary, introduced to eliminate observed implausible spectra
    !> It's very strict: one only outranged value will eliminate the whole (co)spectrum
    ol: do var = w, gas4
        il: do i = 1, nbins
            if (BinSpec(i)%of(var) > MaxNormSpecValue) then
                BinSpec(:)%of(var) = error
                exit il
            end if
        end do il
    end do ol

    ol2: do var = w_ts, w_gas4
        il2: do i = 1, nbins
            if (dabs(BinCosp(i)%of(var)) > MaxNormSpecValue) then
                BinCosp(:)%of(var) = error
                exit il2
            end if
        end do il2
    end do ol2

    !> Similar filter as above, but now imposes that f*spectrum < 1 for each frequency
    ol3: do var = w, gas4
        il3: do i = 1, nbins
            if (BinSpec(i)%fn /= error .and. BinSpec(i)%of(var) /= error .and. BinSpec(i)%fn * BinSpec(i)%of(var) > 1d0) then
                BinSpec(:)%of(var) = error
                exit il3
            end if
        end do il3
    end do ol3

    ol4: do var = w_ts, w_gas4
        il4: do i = 1, nbins
            if (BinCosp(i)%fn /= error .and. BinCosp(i)%of(var) /= error .and. &
                dabs(BinCosp(i)%fn * BinCosp(i)%of(var)) > 10d0) then
                BinCosp(:)%of(var) = error
                exit il4
            end if
        end do il4
    end do ol4
end subroutine ReadBinnedFile
