!***************************************************************************
! spectral_analysis.f90
! ---------------------
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
! \brief       Calculate spectra and co-spectra \n
!              and write them on output file as requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SpectralAnalysis(date, time, bf, Set, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(in) :: Set(N, M)
    real(kind = dbl), intent(in) :: bf(Meth%spec%nbins + 1)
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    !> local variables
    integer :: i
    integer :: var
    integer :: bcnt(Meth%spec%nbins)
    real(kind = dbl) :: bnf(Meth%spec%nbins)
    real(kind = dbl) :: AuxSet(N, M)
    real(kind = dbl) :: nf(N/2)
    real(kind = dbl) :: sumw, auxsumw
    character(13) :: Datestring
    type (SpectralType) :: Spectrum(N/2 + 1)
    type (SpectralType) :: AuxSpectrum(N/2 + 1)
    type (SpectralType) :: BinnedSpectrum(Meth%spec%nbins)
    type (SpectralType) :: Cospectrum(N/2 + 1)
    type (SpectralType) :: AuxCospectrum(N/2 + 1)
    type (SpectralType) :: BinnedCospectrum(Meth%spec%nbins)
    type (SpectralType) :: Ogive(N/2 + 1)
    type (SpectralType) :: CoOgive(N/2 + 1)
    type (SpectralType) :: BinnedOgive(Meth%spec%nbins)
    type (SpectralType) :: BinnedCoOgive(Meth%spec%nbins)
    logical :: DoSpectrum(GHGNumVar)
    logical :: DoCospectrum(GHGNumVar)
    logical :: proceed

    !> Determine which spectra and cospectra can be calculated
    call DetectFeasibleSpectraAndCospectra(DoSpectrum, DoCospectrum)

    !> Determine if at least one (co)spectrum has to be written on output
    proceed = .false.
    do var = 1, GHGNumVar
        if (RPsetup%out_full_sp(var) .or. RPsetup%out_full_cosp(var)) then
            proceed = .true.
            exit
        end if
    end do

    !> If binned (co)spectra or at least one full (co)spectrum are requested,
    !> perform all related calculations
    write(*, '(a)') '  Calculating (co)spectra..'
    Datestring = date(1:4) // date(6:7) // date(9:10) &
               // '-' // time(1:2) // time(4:5)
    !> Define the "natural frequency" vector nf, extending
    !> from f_min = 1/N --> f_max = AcFreq/2 (= Nyquist frequency)
    !> it is period-dependent
    do i = 1, N/2
        nf(i) = dble(i) * Metadata%ac_freq / dble(N)
    end do
    !> Use "Aux" variables to calculate degraded covariances and to
    !> output full co-spectrum wT. "Aux" variables are not tapered
    AuxSet = Set
    !> (do not) taper dataset (squared window makes no tapering), but needed to
    !> calculate "window squared and summed"
    call Tapering('squared', AuxSet, N, M, auxsumw)

    !> Fft and calculate cospectra
    call FourierTransform(AuxSet, N, M)
    call AllCospectra(AuxSet, auxsumw, AuxSpectrum, AuxCospectrum, &
        DoSpectrum, DoCospectrum, N, M)

    !> Calculate ogives if requested
    if (RPsetup%out_bin_og) &
        call AllOgives(AuxSpectrum, AuxCospectrum, DoSpectrum, DoCospectrum, Ogive, CoOgive, N)

    !> Write out full co-spectra as requested
    if (proceed) &
        call WriteOutFullCoSpectra(Datestring, nf, &
            AuxSpectrum, AuxCospectrum, DoSpectrum, DoCospectrum, N)

    !> Output degraded covariances calculated by integrating filtered cospectra
    if (DoCospectrum(w_ts)) then
        call DegradedCovariances(nf, AuxCospectrum, N)
    else
        Essentials%degH(:) = error
    end if

    !> If binned (co-)spectra are requested, perform all related calculations
    if (RPsetup%out_bin_sp .or. RPsetup%out_bin_og) then
        !> From here, tapering is applied.
        !> Taper dataset
        call Tapering(RPsetup%tap_win, Set, N, M, sumw)
        !> Fft and calculate cospectra
        call FourierTransform(Set, N, M)
        call AllCospectra(Set, sumw, Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N, M)
        !> Normalization by variances and covariances
        call NormalizeCoSpectra(Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N)
        !> Binned cospectra session
        if (RPsetup%out_bin_sp) then
            !> Exponential binning of frequencies, spectra and co-spectra
            call ExpAvrgCospectra(bf, nf, Spectrum, Cospectrum, N/2 + 1, bnf &
                , BinnedSpectrum, BinnedCospectrum, bcnt)
            !> Write co-spectra on output file in csv format
            call WriteOutBinnedCoSpectra(Datestring, bnf, bcnt, BinnedSpectrum, BinnedCospectrum &
                , DoSpectrum, DoCospectrum)
        end if
        !> Ogive session
        if (RPsetup%out_bin_og) then
            !> Exponential binning of ogives
            call ExpAvrgOgives(bf, nf, Ogive, CoOgive, N/2 + 1, bnf &
                , BinnedOgive, BinnedCoOgive, bcnt)
            !> Write co-ogives on output file in csv format
            call WriteOutBinnedOgives(Datestring, bnf, bcnt, BinnedOgive, BinnedCoOgive &
                , DoSpectrum, DoCospectrum)
        end if
    end if
    write(*,'(a)') '  Done.'
end subroutine SpectralAnalysis

!***************************************************************************
!
! \brief       Determine which spectra and cospectra can be evaluated, given \n
!              available variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DetectFeasibleSpectraAndCospectra(DoSpectrum, DoCospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(inout) :: DoSpectrum(GHGNumVar)
    logical, intent(inout) :: DoCospectrum(GHGNumVar)
    !> local variables
    integer :: i
    character(9) :: hf_sr
    character(9) :: hf_do
    character(9) :: hf_sk, sf_sk
    character(9) :: hf_ds, sf_ds


    !> Spikes test (only if not despiked)
    if (RPsetup%filter_spectra_by_qc) then
        if (.not. RPsetup%filter_sr) then
            call int2char(IntHF%sr, hf_sr, 9)
            hf_sr(1:8) = hf_sr(2:9)
        else
            hf_sr = '000000000'
        end if

        !> Drop out test
        call int2char(IntHF%do, hf_do, 9)
        hf_do(1:8) = hf_do(2:9)

        !> Skewness/kurtosis test
        call int2char(IntHF%sk, hf_sk, 9)
        call int2char(IntSF%sk, sf_sk, 9)
        hf_sk(1:8) = hf_sk(2:9)
        sf_sk(1:8) = sf_sk(2:9)

        !> Discontinuities test
        call int2char(IntHF%ds, hf_ds, 9)
        call int2char(IntSF%ds, sf_ds, 9)
        hf_ds(1:8) = hf_ds(2:9)
        sf_ds(1:8) = sf_ds(2:9)
    else
        hf_sr = '000000000'
        hf_do = '000000000'
        hf_sk = '000000000'
        sf_sk = '000000000'
        hf_ds = '000000000'
        sf_ds = '000000000'
    end if

    !> Decision based on variables availability and on results of statistical tests (if the case)
    !> Spectra
    DoSpectrum = .false.
    do i = u, w
        if (SpecCol(i)%present &
            .and. hf_sr(i:i) /= '1' .and. hf_do(i:i) /= '1' &
            .and. hf_sk(i:i) /= '1' .and. sf_sk(i:i) /= '1' &
            .and. hf_ds(i:i) /= '1' .and. sf_ds(i:i) /= '1') &
            DoSpectrum(i) = .true.
    end do
    do i = ts, GHGNumVar
        !> NOTE: Control on soft flags eliminated for the scalars
        if (SpecCol(i)%present &
            .and. hf_sr(i:i) /= '1' .and. hf_do(i:i) /= '1' &
            .and. hf_sk(i:i) /= '1' .and. hf_ds(i:i) /= '1') &
            DoSpectrum(i) = .true.
    end do

    !> Cospectra
    DoCospectrum = .false.
    if (SpecCol(w)%present .and. SpecCol(u)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_do(w:w) /= '1' &
        .and. hf_sr(u:u) /= '1' .and. hf_do(u:u) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. sf_sk(w:w) /= '1' &
        .and. hf_sk(u:u) /= '1' .and. sf_sk(u:u) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. sf_ds(w:w) /= '1' &
        .and. hf_ds(u:u) /= '1' .and. sf_ds(u:u) /= '1') &
         DoCospectrum(w_u)   = .true.

    if (SpecCol(w)%present .and. SpecCol(v)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_do(w:w) /= '1' &
        .and. hf_sr(v:v) /= '1' .and. hf_do(v:v) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. sf_sk(w:w) /= '1' &
        .and. hf_sk(v:v) /= '1' .and. sf_sk(v:v) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. sf_ds(w:w) /= '1' &
        .and. hf_ds(v:v) /= '1' .and. sf_ds(v:v) /= '1') &
        DoCospectrum(w_v)   = .true.

    !> NOTE: Control on soft flags not used for scalars (Ts and gases)
    if (SpecCol(w)%present .and. SpecCol(ts)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_sr(ts:ts) /= '1' &
        .and. hf_do(w:w) /= '1' .and. hf_do(ts:ts) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. hf_sk(ts:ts) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. hf_ds(ts:ts) /= '1' &
        .and. sf_sk(w:w) /= '1' .and. sf_ds(w:w) /= '1') &
        DoCospectrum(w_ts) = .true.

    if (SpecCol(w)%present .and. SpecCol(co2)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_sr(co2:co2) /= '1' &
        .and. hf_do(w:w) /= '1' .and. hf_do(co2:co2) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. hf_sk(co2:co2) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. hf_ds(co2:co2) /= '1' &
        .and. sf_sk(w:w) /= '1' .and. sf_ds(w:w) /= '1') &
         DoCospectrum(w_co2) = .true.

    if (SpecCol(w)%present .and. SpecCol(h2o)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_sr(h2o:h2o) /= '1' &
        .and. hf_do(w:w) /= '1' .and. hf_do(h2o:h2o) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. hf_sk(h2o:h2o) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. hf_ds(h2o:h2o) /= '1' &
        .and. sf_sk(w:w) /= '1' .and. sf_ds(w:w) /= '1') &
        DoCospectrum(w_h2o) = .true.

    if (SpecCol(w)%present .and. SpecCol(ch4)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_sr(ch4:ch4) /= '1' &
        .and. hf_do(w:w) /= '1' .and. hf_do(ch4:ch4) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. hf_sk(ch4:ch4) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. hf_ds(ch4:ch4) /= '1' &
        .and. sf_sk(w:w) /= '1' .and. sf_ds(w:w) /= '1') &
        DoCospectrum(w_ch4) = .true.

    if (SpecCol(w)%present .and. SpecCol(gas4)%present &
        .and. hf_sr(w:w) /= '1' .and. hf_sr(gas4:gas4) /= '1' &
        .and. hf_do(w:w) /= '1' .and. hf_do(gas4:gas4) /= '1' &
        .and. hf_sk(w:w) /= '1' .and. hf_sk(gas4:gas4) /= '1' &
        .and. hf_ds(w:w) /= '1' .and. hf_ds(gas4:gas4) /= '1' &
        .and. sf_sk(w:w) /= '1' .and. sf_ds(w:w) /= '1') &
        DoCospectrum(w_n2o) = .true.

end subroutine DetectFeasibleSpectraAndCospectra

!***************************************************************************
!
! \brief       Calculate cospectral densities, starting from un-normalised Fourier \n
!              coefficients, as calculated by rfftf.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AllCospectra(Set, sumw, Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind=dbl), intent(in) :: sumw
    real(kind = dbl), intent(in) :: Set(N, M)
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    Type(SpectralType), intent(out) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(out) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: j
    real(kind = dbl) :: xx(N)
    real(kind = dbl) :: yy(N)

    write(*, '(a)', advance = 'no') '   Cospectral densities..'

    !> spectra
    do j = u, GHGNumVar
        if (DoSpectrum(j)) then
            xx(1:N) = Set(1:N, j)
            call SpectralDensity(xx, xx, Metadata%ac_freq, sumw, Spectrum%of(j), N)
        end if
    end do

    !> Cospectra
    do j = w_u, w_n2o
        if (j /= w_w) then
            if (DoCospectrum(j)) then
                xx(1:N) = Set(1:N, w)
                yy(1:N) = Set(1:N, j)
                call SpectralDensity(xx, yy, Metadata%ac_freq, sumw, Cospectrum%of(j), N)
            end if
        end if
    end do
    write(*,'(a)') ' done.'
end subroutine AllCospectra

!***************************************************************************
!
! \brief       Calculate cospectral densities, starting from un-normalised Fourier \n
!              coefficients, as calculated by rfftf.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AllOgives(Spectrum, Cospectrum, DoSpectrum, DoCospectrum, Ogive, CoOgive, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    Type(SpectralType), intent(in) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(in) :: Cospectrum(N/2 + 1)
    Type(SpectralType), intent(out) :: Ogive(N/2 + 1)
    Type(SpectralType), intent(out) :: CoOgive(N/2 + 1)
    !> local variables
    integer :: i
    integer :: j
    real (kind = dbl) :: df

    write(*, '(a)', advance = 'no') '   Ogives..'

    df = Metadata%ac_freq / dble(N)

    !> Ogive
    do j = u, GHGNumVar
        if (DoSpectrum(j)) then
            Ogive(N/2+1)%Of(j) = Spectrum(N/2+1)%of(j)
            do i = N/2, 1, -1
                Ogive(i)%of(j) = Ogive(i+1)%of(j) + Spectrum(i)%of(j) * df
            end do
        end if
    end do

    !> CoOgive
    do j = w_u, w_n2o
        if (j /= w_w) then
            if (DoCospectrum(j)) then
                CoOgive(N/2+1)%Of(j) = Cospectrum(N/2+1)%of(j)
                do i = N/2, 1, -1
                    CoOgive(i)%Of(j) = CoOgive(i+1)%of(j) + Cospectrum(i)%of(j)  * df
                end do
            end if
        end if
    end do

    !> Normalize
    do j = u, GHGNumVar
        if (DoSpectrum(j) .and. Stats%Cov(j, j) /= 0d0 .and. Stats%Cov(j, j) /= error) &
            Ogive(1:N/2 + 1)%of(j) = Ogive(1:N/2 + 1)%of(j) / Stats%Cov(j, j)
    end do

    do j = w_u, w_n2o
    if (DoCospectrum(j) .and. Stats%Cov(w, j) /= 0d0 .and. Stats%Cov(w, j) /= error) &
        CoOgive(1:N/2 + 1)%of(j) = CoOgive(1:N/2 + 1)%of(j) / Stats%Cov(w, j)
    end do

    write(*,'(a)') ' done.'
end subroutine AllOgives

!***************************************************************************
!
! \brief       Normalize calculated (co)spectra with the relevant (co)variance
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine NormalizeCoSpectra(Spectrum, Cospectrum, DoSpectrum, DoCospectrum, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: DoSpectrum(GHGNumVar)
    logical, intent(in) :: DoCospectrum(GHGNumVar)
    Type(SpectralType), intent(inout) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(inout) :: Cospectrum(N/2 + 1)
    !> local variables
    integer :: j

    do j = u, gas4
        if (DoSpectrum(j) .and. Stats%Cov(j, j) /= 0d0 .and. Stats%Cov(j, j) /= error) &
            Spectrum(1:N/2 + 1)%of(j) = Spectrum(1:N/2 + 1)%of(j) / Stats%Cov(j, j)
    end do

    do j = w_u, w_n2o
    if (DoCospectrum(j) .and. Stats%Cov(w, j) /= 0d0 .and. Stats%Cov(w, j) /= error) &
        Cospectrum(1:N/2 + 1)%of(j) = Cospectrum(1:N/2 + 1)%of(j) / Stats%Cov(w, j)
    end do
end subroutine NormalizeCoSpectra

!***************************************************************************
!
! \brief       Average (co)spectra in the binned frequency range
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExpAvrgCospectra(bf, nf, Spectrum, Cospectrum, N, bin_nf, &
    BinnedSpectrum, BinnedCospectrum, bin_cnt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: bf(Meth%spec%nbins + 1)
    real(kind = dbl), intent(in) :: nf(N)
    integer, intent(out) :: bin_cnt(Meth%spec%nbins)
    real(kind = dbl), intent(out) :: bin_nf(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedSpectrum(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedCospectrum(Meth%spec%nbins)
    type (SpectralType), intent(inout) :: Spectrum(N)
    type (SpectralType), intent(inout) :: Cospectrum(N)
    !> local variables
    integer :: i = 0
    integer :: j = 0

    !> Align co-spectra to frequencies before averaging
    do i = 1, N-1
        Spectrum(i) = Spectrum(i + 1)
        Cospectrum(i) = Cospectrum(i + 1)
    end do
    Spectrum(N)%of(:) = 0d0
    Cospectrum(N)%of(:) = 0d0

    write(*, '(a)', advance = 'no') '   Binning spectra and cospectra..'
    !> average variables in the exp-spaced ranges
    do i = 1, Meth%spec%nbins
        bin_cnt(i) = 0
        bin_nf(i) = 0.d0
        BinnedSpectrum(i)%of(u:gas4) = 0d0
        BinnedCospectrum(i)%of(w_u:w_n2o) = 0d0
        do j = 1, N
            if(nf(j) > bf(i) .and. nf(j) <= bf(i + 1)) then
                bin_nf(i) = bin_nf(i) + nf(j)
                BinnedSpectrum(i)%of(u:gas4) = &
                    BinnedSpectrum(i)%of(u:gas4) + Spectrum(j)%of(u:gas4)
                BinnedCospectrum(i)%of(w_u:w_n2o) = &
                    BinnedCospectrum(i)%of(w_u:w_n2o) + Cospectrum(j)%of(w_u:w_n2o)
                bin_cnt(i) = bin_cnt(i) + 1
            end if
        end do
        if(bin_cnt(i) /= 0) then
            bin_nf(i)    = bin_nf(i) / dble(bin_cnt(i))
            BinnedSpectrum(i)%of(u:gas4) = &
                BinnedSpectrum(i)%of(u:gas4) / dble(bin_cnt(i))
            BinnedCospectrum(i)%of(w_u:w_n2o) = &
                BinnedCospectrum(i)%of(w_u:w_n2o) / dble(bin_cnt(i))
        else
            bin_nf(i)    = error
            BinnedSpectrum(i)%of(u:gas4) = error
            BinnedCospectrum(i)%of(w_u:w_n2o) = error
        end if
    end do
    write(*,'(a)') ' done.'
end subroutine ExpAvrgCospectra

!***************************************************************************
!
! \brief       Average ogives in the binned frequency range
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExpAvrgOgives(bf, nf, Ogive, CoOgive, N, bin_nf, &
    BinnedOgive, BinnedCoOgive, bin_cnt)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: bf(Meth%spec%nbins + 1)
    real(kind = dbl), intent(in) :: nf(N)
    integer, intent(out) :: bin_cnt(Meth%spec%nbins)
    real(kind = dbl), intent(out) :: bin_nf(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedOgive(Meth%spec%nbins)
    type (SpectralType), intent(out) :: BinnedCoOgive(Meth%spec%nbins)
    type (SpectralType), intent(inout) :: Ogive(N)
    type (SpectralType), intent(inout) :: CoOgive(N)
    !> local variables
    integer :: i = 0
    integer :: j = 0

    !> Align co-ogives to frequencies before averaging
    do i = 1, N-1
        Ogive(i) = Ogive(i + 1)
        CoOgive(i) = CoOgive(i + 1)
    end do
    Ogive(N)%of(:) = 0d0
    CoOgive(N)%of(:) = 0d0

    write(*, '(a)', advance = 'no') '   Binning ogives..'
    !> average variables in the exp-spaced ranges
    do i = 1, Meth%spec%nbins
        bin_cnt(i) = 0
        bin_nf(i) = 0.d0
        BinnedOgive(i)%of(u:gas4) = 0d0
        BinnedCoOgive(i)%of(w_u:w_n2o) = 0d0
        do j = 1, N
            if(nf(j) >= bf(i) .and. nf(j) < bf(i + 1)) then
                bin_nf(i) = bin_nf(i) + nf(j)
                BinnedOgive(i)%of(u:gas4) = &
                    BinnedOgive(i)%of(u:gas4) + Ogive(j)%of(u:gas4)
                BinnedCoOgive(i)%of(w_u:w_n2o) = &
                    BinnedCoOgive(i)%of(w_u:w_n2o) + CoOgive(j)%of(w_u:w_n2o)
                bin_cnt(i) = bin_cnt(i) + 1
            end if
        end do
        if(bin_cnt(i) /= 0) then
            bin_nf(i)    = bin_nf(i)    / dble(bin_cnt(i))
            BinnedOgive(i)%of(u:gas4) = &
                BinnedOgive(i)%of(u:gas4) / dble(bin_cnt(i))
            BinnedCoOgive(i)%of(w_u:w_n2o) = &
                BinnedCoOgive(i)%of(w_u:w_n2o) / dble(bin_cnt(i))
        else
            bin_nf(i)    = error
            BinnedOgive(i)%of(u:gas4) = error
            BinnedCoOgive(i)%of(w_u:w_n2o) = error
        end if
    end do
    write(*,'(a)') ' done.'
end subroutine ExpAvrgOgives

!***************************************************************************
!
! \brief       Write binned (co)spectra on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBinnedCoSpectra(String, bnf, bcnt, BinnedSpectrum, BinnedCospectrum &
    , DoSpectrum, DoCospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bcnt(Meth%spec%nbins)
    real(kind = dbl), intent(in) :: bnf(Meth%spec%nbins)
    character(*), intent(in) :: String
    type (SpectralType), intent(in) :: BinnedSpectrum(Meth%spec%nbins)
    type (SpectralType), intent(in) :: BinnedCospectrum(Meth%spec%nbins)
    logical, intent(in)  :: DoSpectrum(GHGNumVar)
    logical, intent(in)  :: DoCospectrum(GHGNumVar)
    !> local variables
    integer :: i
    integer :: j
    character(64) :: e2sg(E2NumVar)
    character(256) :: BinCospectraPath
    character(10000) :: dataline = ''
    character(30) :: datum = ''

    e2sg(gas4) = SpecCol(gas4)%label(1:len_trim(SpecCol(gas4)%label))

    !> Open output file for binned co-spectra
    BinCospectraPath = BinCospectraDir(1:len_trim(BinCospectraDir)) // String &
                 // BinCospec_FilePadding // Timestamp_FilePadding // CsvExt
    open(udf, file = BinCospectraPath, encoding = 'utf-8')
    write(udf, *)           'normalised_and_exponentially_binned_(co)spectra'
    write(udf, *)           '-----------------------------------------------'
    write(udf, '(a)')       'natural_frequency_[Hz]_->_from_[1/_averaging_interval]_to_[acquisition_frequency_/_2]'
    write(udf, '(a)')       'normalized_frequency_[#]_->_natural_frequency_*_measuring_height_/_wind_speed'
    write(udf, '(a)')       'y-axis_->_natural_frequency_*_(co)spectrum_/_(co)variance'
    write(udf, '(a, f7.3)') 'acquisition_frequency_[Hz]_=_', Metadata%ac_freq
    write(udf, '(a, f7.3)') 'measuring_height_(z-d)_[m]_=_', (SpecCol(u)%Instr%height - Metadata%d)
    write(udf, '(a, f7.3)') 'wind_speed_[m+1s-1]_=_', LitePar%WS
    write(udf, '(a, i7)')   'averaging_interval_[min]_=_', RPsetup%avrg_len
    write(udf, '(a, i4)')   'number_of_bins_=_', Meth%spec%nbins
    write(udf, '(a, a)')    'tapering_window_=_', RPsetup%tap_win(1:len_trim(RPsetup%tap_win))
    write(udf, *) '#_freq,natural_frequency,normalized_frequency,f_nat*spec(u)/var(u),f_nat*spec(v)/var(v)&
        &,f_nat*spec(w)/var(w),f_nat*spec(ts)/var(ts),f_nat*spec(co2)/var(co2),f_nat*spec(h2o)/var(h2o),f_nat*spec(ch4)/var(ch4)&
        &,f_nat*spec(' // e2sg(gas4)(1:len_trim(e2sg(gas4)))// ')/var(' // e2sg(gas4)(1:len_trim(e2sg(gas4)))// ')&
        &,f_nat*cospec(w_u)/cov(w_u),f_nat*cospec(w_v)/cov(w_v),f_nat*cospec(w_ts)/cov(w_ts),f_nat*cospec(w_co2)/cov(w_co2)&
        &,f_nat*cospec(w_h2o)/cov(w_h2o),f_nat*cospec(w_ch4)/cov(w_ch4)&
        &,f_nat*cospec(w_' // e2sg(gas4)(1:len_trim(e2sg(gas4)))// ')/cov(w_' // e2sg(gas4)(1:len_trim(e2sg(gas4)))// ')'

    !> Write to output file in csv style
    do i = 1, Meth%spec%nbins
        call clearstr(dataline)
        call WriteDatumInt(bcnt(i), datum, '-9999')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(bnf(i), datum, '-9999.0')
        call AddDatum(dataline, datum, separator)
        if (bnf(i) /= error .and. LitePar%WS /= 0d0) then
            call WriteDatumFloat(bnf(i) * (SpecCol(u)%Instr%height - Metadata%d) / LitePar%WS, datum, '-9999.0')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, '-9999.0', separator)
        end if

        !> Spectra
        do j = u, gas4
            if (DoSpectrum(j) .and. bnf(i) /= error .and. BinnedSpectrum(i)%of(j) /= error) then
                call WriteDatumFloat(bnf(i) * BinnedSpectrum(i)%of(j), datum, '-9999.0')
                call AddDatum(dataline, datum, separator)
            else
                call AddDatum(dataline, '-9999.0', separator)
            end if
        end do
        !> Cospectra
        do j = w_u, w_n2o
            if (j /= w_w) then
                if (DoCospectrum(j) .and. bnf(i) /= error .and. BinnedCospectrum(i)%of(j) /= error) then
                    call WriteDatumFloat(bnf(i) * BinnedCospectrum(i)%of(j), datum, '-9999.0')
                    call AddDatum(dataline, datum, separator)
                else
                    call AddDatum(dataline, '-9999.0', separator)
                end if
            end if
        end do
        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
    end do
    close(udf)
end subroutine WriteOutBinnedCoSpectra

!***************************************************************************
!
! \brief       Write binned ogives on output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutBinnedOgives(String, bnf, bcnt, BinnedOgive, BinnedCoOgive &
    , DoSpectrum, DoCospectrum)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: bcnt(Meth%spec%nbins)
    real(kind = dbl), intent(in) :: bnf(Meth%spec%nbins)
    character(*), intent(in) :: String
    type (SpectralType), intent(in) :: BinnedOgive(Meth%spec%nbins)
    type (SpectralType), intent(in) :: BinnedCoOgive(Meth%spec%nbins)
    logical, intent(in)  :: DoSpectrum(GHGNumVar)
    logical, intent(in)  :: DoCospectrum(GHGNumVar)
    !> local variables
    integer :: i
    integer :: j
    character(64) :: e2sg(E2NumVar)
    character(256) :: BinOgivesPath
    character(10000) :: dataline = ''
    character(30) :: datum = ''

    e2sg(gas4) = SpecCol(gas4)%label(1:len_trim(SpecCol(gas4)%label))

    !> Open output file for binned co-spectra
    BinOgivesPath = BinOgivesDir(1:len_trim(BinOgivesDir)) // String &
                 // BinOgives_FilePadding // Timestamp_FilePadding // CsvExt
    open(udf, file = BinOgivesPath, encoding = 'utf-8')
    write(udf, '(a)')       'exponentially_binned_ogives'
    write(udf, '(a)')       '---------------------------'
    write(udf, '(a)')       'natural_frequency_[Hz]_->_from_[1/_averaging_interval]_to_[acquisition_frequency_/_2]'
    write(udf, '(a)')       'normalized_frequency_[#]_->_natural_frequency_*_measuring_height_/_wind_speed'
    write(udf, '(a)')       'y-axis_->_ogive'
    write(udf, '(a, f7.3)') 'acquisition_frequency_[Hz]_=_', Metadata%ac_freq
    write(udf, '(a, f7.3)') 'measuring_height_(z-d)_[m]_=_', (SpecCol(u)%Instr%height - Metadata%d)
    write(udf, '(a, f7.3)') 'wind_speed_[m+1s-1]_=_', LitePar%WS
    write(udf, '(a, i7)')   'averaging_interval_[min]_=_', RPsetup%avrg_len
    write(udf, '(a, i4)')   'number_of_bins_=_', Meth%spec%nbins
    write(udf, '(a, a)')    'tapering_window_=_', RPsetup%tap_win(1:len_trim(RPsetup%tap_win))
    write(udf, *) '#_freq,natural_frequency,normalized_frequency,og(u),og(v),og(w),og(ts),og(co2),og(h2o),og(ch4)&
        &,og(' // e2sg(gas4)(1:len_trim(e2sg(gas4)))// ')&
        &,og(w_u),og(w_v),og(w_ts),og(w_co2),og(w_h2o),og(w_ch4),og(w_' // e2sg(gas4)(1:len_trim(e2sg(gas4)))// ')'

    !> Write to output file in csv style
    do i = 1, Meth%spec%nbins
        call clearstr(dataline)

        call WriteDatumInt(bcnt(i), datum, '-9999')
        call AddDatum(dataline, datum, separator)
        call WriteDatumFloat(bnf(i), datum, '-9999.0')
        call AddDatum(dataline, datum, separator)
        if (bnf(i) /= error .and. LitePar%WS /= 0d0) then
            call WriteDatumFloat(bnf(i) * (SpecCol(u)%Instr%height - Metadata%d) / LitePar%WS, datum, '-9999.0')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, '-9999.0', separator)
        end if

        !> Ogives
        do j = u, gas4
            if (DoSpectrum(j)) then
                call WriteDatumFloat(BinnedOgive(i)%of(j), datum, '-9999.0')
                call AddDatum(dataline, datum, separator)
            else
                call AddDatum(dataline, '-9999.0', separator)
            end if
        end do
        !> Co-ogives
        do j = w_u, w_n2o
            if (j /= w_w) then
                if (DoCospectrum(j)) then
                    call WriteDatumFloat(BinnedCoOgive(i)%of(j), datum, '-9999.0')
                    call AddDatum(dataline, datum, separator)
                else
                    call AddDatum(dataline, '-9999.0', separator)
                end if
            end if
        end do
        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
    end do
    close(udf)
end subroutine WriteOutBinnedOgives

!***************************************************************************
!
! \brief       Write full (co)spectra on output file as requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFullCoSpectra(String, nf, Spectrum, Cospectrum, &
    DoSpectrum, DoCospectrum, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: nf(N/2)
    Type(SpectralType), intent(in) :: Spectrum(N/2 + 1)
    Type(SpectralType), intent(in) :: Cospectrum(N/2 + 1)
    logical, intent(in)  :: DoSpectrum(GHGNumVar)
    logical, intent(in)  :: DoCospectrum(GHGNumVar)
    character(*), intent(in) :: String
    !> local variables
    integer :: i
    integer :: var
    character(256) :: CospectraPath
    character(10000) :: dataline  = ''
    character(10000) :: dataline1 = ''
    character(10000) :: dataline2 = ''
    character(10000) :: dataline3 = ''
    character(30) :: datum = ''
    character(4) :: e2sg(GHGNumVar)

    write(*, '(a)', advance = 'no') '   Writing requested full (co)spectra on output file..'

    !> Convenient strings
    e2sg(u)   = 'u'
    e2sg(v)   = 'v'
    e2sg(w)   = 'w'
    e2sg(ts)  = 'ts'
    e2sg(co2) = 'co2'
    e2sg(h2o) = 'h2o'
    e2sg(ch4) = 'ch4'
    e2sg(gas4) = 'gas4'

    !> Open output file for binned co-spectra
    CospectraPath = CospectraDir(1:len_trim(CospectraDir)) // String &
                 // Cospec_FilePadding // Timestamp_FilePadding // CsvExt
    open(udf, file = CospectraPath, encoding = 'utf-8')
    write(udf, '(a)')       'Full_(co)spectra'
    write(udf, '(a)')       '----------------'
    write(udf, '(a)')       'natural_frequency_[Hz]_->_from_[1/_averaging_interval]_to_[acquisition_frequency_/_2]'
    write(udf, '(a)')       'normalized_frequency_[#]->_natural_frequency_*_measuring_height_/_wind_speed'
    write(udf, '(a)')       'y-axis_->_natural_frequency_*_(co)spectrum_/_(co)variance'
    write(udf, '(a, f7.3)') 'acquisition_frequency_[Hz]_=_', Metadata%ac_freq
    write(udf, '(a, f7.3)') 'measuring_height_(z-d)_[m]_=_', (SpecCol(u)%Instr%height - Metadata%d)
    write(udf, '(a, f7.3)') 'wind_speed_[m+1s-1]_=_', LitePar%WS
    write(udf, '(a, i7)')   'averaging_interval_[min]_=_', RPsetup%avrg_len
    write(udf, '(a)')       'tapering_window_=_SQUARED_(no_tapering_forced_by_EddyPro.&
                            &_Tapering_is_only_applied_for_binned_(co)spectra)'

    !> First writes (co)variances (header + numbers) and (co)spectra labels
    call clearstr(dataline1)
    call clearstr(dataline2)
    call clearstr(dataline3)
    call AddDatum(dataline1, ',', separator)
    call AddDatum(dataline2, ',', separator)
    call AddDatum(dataline3, 'natural_frequency,normalized_frequency', separator)
    do var = 1, GHGNumVar
        if (RPsetup%out_full_sp(var) .and. DoSpectrum(var)) then
            call AddDatum(dataline1, 'var(' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
            write(datum, *) Stats%cov(var, var)
            call AddDatum(dataline2, datum, separator)
            call AddDatum(dataline3, 'f_nat*spec(' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
        end if
    end do
    do var = 1, GHGNumVar
        if (RPsetup%out_full_cosp(var) .and. DoCospectrum(var)) then
            call AddDatum(dataline1, 'cov(w_' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
            write(datum, *) Stats%cov(w, var)
            call AddDatum(dataline2, datum, separator)
            call AddDatum(dataline3, 'f_nat*cospec(w_' //e2sg(var)(1:len_trim(e2sg(var))) // ')', separator)
        end if
    end do
    write(udf, '(a)') dataline1(1:len_trim(dataline1) - 1)
    write(udf, '(a)') dataline2(1:len_trim(dataline2) - 1)
    write(udf, '(a)') dataline3(1:len_trim(dataline3) - 1)

    !> Now write spectra and cospectra
    do i = 1, N/2
        call clearstr(dataline)
        !> Frequencies
        call WriteDatumFloat(nf(i), datum, '-9999.0')
        call AddDatum(dataline, datum, separator)
        if (nf(i) /= error .and. LitePar%WS /= 0d0) then
            call WriteDatumFloat(nf(i) * (SpecCol(u)%Instr%height - Metadata%d) / LitePar%WS, datum, '-9999.0')
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, '-9999.0', separator)
        end if
        !> spectra
        do var = 1, GHGNumVar
            if (RPsetup%out_full_sp(var) .and. DoSpectrum(var)) then
                if (nf(i) /= error .and. Stats%cov(var, var) /= 0d0 .and. Stats%cov(var, var) /= error) then
                    call WriteDatumFloat(nf(i) * Spectrum(i)%of(var) / Stats%cov(var, var), datum, '-9999.0')
                    call AddDatum(dataline, datum, separator)
                else
                    call AddDatum(dataline, '-9999.0', separator)
                end if
            end if
        end do
        !> cospectra
        do var = 1, GHGNumVar
            if (RPsetup%out_full_cosp(var) .and. DoCospectrum(var)) then
                if (nf(i) /= error .and. Stats%cov(w, var) /= 0d0 .and. Stats%cov(w, var) /= error) then
                    call WriteDatumFloat(nf(i) * Cospectrum(i)%of(var) / Stats%cov(w, var), datum, '-9999.0')
                    call AddDatum(dataline, datum, separator)
                else
                    call AddDatum(dataline, '-9999.0', separator)
                end if
            end if
        end do
        write(udf, '(a)') dataline(1:len_trim(dataline) - 1)
    end do
    close(udf)
    write(*,'(a)') '  done.'
end subroutine WriteOutFullCoSpectra
