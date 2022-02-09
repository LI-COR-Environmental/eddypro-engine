module mo_fftmax

    ! Written Feb 2022, Matthias Cuntz

    use m_numeric_kinds, only: i4, sp => sgl, dp => dbl

    implicit none

    private

    public :: correl  ! correlation function between two vectors
    public :: fftmax  ! index of maximum correlation and actual time lag in s
    public :: mean    ! 1st moment of an array

    integer, parameter :: spc = sp
    integer, parameter :: dpc = dp

    ! ------------------------------------------------------------------

contains

    ! ------------------------------------------------------------------

    ! Fast Fourier Transform of complex numbers
    function zfft(x)

        use fftpack, only: zffti, zfftf

        complex(dpc), dimension(:), intent(in) :: x
        complex(dpc), dimension(size(x, 1)) :: zfft

        integer :: n
        real(dpc), dimension(:), allocatable :: wsave

        n = size(x)
        zfft = x

        ! initialise
        allocate(wsave(4*n+15))
        call zffti(n, wsave)

        ! FFT
        call zfftf(n, zfft, wsave)

        deallocate(wsave)

        return

    end function zfft

    ! Inverse Fast Fourier Transform of complex numbers
    function izfft(x)

        use fftpack, only: zffti, zfftb

        complex(dpc), dimension(:), intent(in) :: x
        complex(dpc), dimension(size(x, 1)) :: izfft

        integer(i4) :: n
        real(dpc), dimension(:), allocatable :: wsave

        n = size(x)
        izfft = x

        ! initialise
        allocate(wsave(4*n+15))
        call zffti(n, wsave)

        ! inverse FFT
        call zfftb(n, izfft, wsave)
        izfft = izfft / real(n, dp)

        deallocate(wsave)

        return

    end function izfft

    ! ------------------------------------------------------------------

    ! correlation function between two vectors
    function correl(x, y, nadjust)

        real(dp), dimension(:), intent(in) :: x
        real(dp), dimension(:), intent(in) :: y
        integer(i4), intent(out), optional :: nadjust
        real(dp), dimension(size(x, 1)) :: correl

        integer(i4) :: n, nf

        n = size(x)
        nf = 2**floor(log(real(n, dp))/log(2.0_dp))
        correl(1:nf) = real(izfft(zfft(cmplx(x(1:nf), kind=dpc)) &
            * conjg(zfft(cmplx(y(1:nf), kind=dpc)))), dp)

        if (present(nadjust)) nadjust = nf

        return

    end function correl

    ! ------------------------------------------------------------------

    subroutine fftmax(col1, col2, error, &
        minlag, maxlag, ac_freq, &
        tlag, rlag, &
        ncorr, crosscorr)

        implicit none

        !> in/out variables
        real(dp), dimension(:), intent(in)  :: col1     ! data
        real(dp), dimension(:), intent(in)  :: col2
        real(dp),               intent(in)  :: error    ! error label of data
        integer,                intent(in)  :: minlag   ! minimum possible time lag
        integer,                intent(in)  :: maxlag   ! maximum possible time lag
        real(dp),               intent(in)  :: ac_freq  ! acquisition frequency
        real(dp),               intent(out) :: tlag     ! time lag [s]
        integer,                intent(out) :: rlag     ! time lag index
        integer,                            intent(out), optional :: ncorr      ! length of 2**n
        real(dp), dimension(size(col1, 1)), intent(out), optional :: crosscorr  ! cross-covariance series

        !> local variables
        integer(i4) :: i30
        real(dp) :: r1corr        ! 1/icorr
        integer  :: icorr         ! adjusted time series to length of 2**n
        integer  :: imax          ! index of maximum covariance
        real(dp), dimension(size(col1, 1)) :: acol1  ! anomalies
        real(dp), dimension(size(col2, 1)) :: acol2
        real(dp), dimension(size(col1, 1)) :: cov12  ! cross-covariance series
        real(dp), dimension(size(col1, 1)) :: corr12 ! cross-correlation series
        logical, dimension(size(col1, 1)) :: mask1  ! data /= error
        logical, dimension(size(col2, 1)) :: mask2  !
        real(dp), dimension(:), allocatable :: corr30 ! correlation series minlag to maxlag
        logical :: dominmax

        dominmax = .true.
        if ((minlag == 0) .and. (maxlag == 0)) then
            ! nothing to do if .not. present(crosscorr)
            dominmax = .false.
            if (.not. present(crosscorr)) then
                rlag = 0
                tlag = 0.0_dp
                return
            endif
        endif

        ! correlation is done with anomalies
        mask1 = col1 /= error
        mask2 = col2 /= error
        acol1 = merge(col1 - mean(col1, mask1), 0.0_dp, mask1)
        acol2 = merge(col2 - mean(col2, mask2), 0.0_dp, mask2)

        ! cross-correlation time series
        cov12 = correl(acol1, acol2, nadjust=icorr)

        ! cross-correlation time series: cov/n
        r1corr = 1.0_dp / real(icorr, dp)
        corr12(1:icorr) = cov12(1:icorr) * r1corr

        ! Sometimes corr is shifted on y-axis. I do not know why.
        ! It could be numeric because of the large amount of data.
        ! Workaround: shift corr around 0 in the tails of the time series,
        ! which are in the middle of corr.
        if (.not. ( (minval(corr12(1:icorr)) < 0.0_dp) .and. &
            (maxval(corr12(1:icorr)) > 0.0_dp)) ) &
            corr12(1:icorr) = corr12(1:icorr) - mean(corr12(icorr/3:2*icorr/3))

        ! Get maximum
        imax = maxloc(abs(corr12(1:icorr)), 1) - 1
        if (imax > icorr/2) imax = imax - icorr

        ! if time lag out of range
        if ((imax < minlag) .or. (imax > maxlag)) then
            if (dominmax) then
                i30 = maxlag + abs(minlag) + 1 ! minlag < 0 possible
                allocate(corr30(i30))
                corr30(1:maxlag+1) = corr12(1:maxlag+1)
                corr30(maxlag+2:i30) = corr12(icorr-abs(minlag)+1:icorr)
                ! maximum
                imax = maxloc(abs(corr30), 1) - 1
                if (imax > (maxlag+1)) imax = imax - i30
                ! for better plotting, uncomment next two lines
                icorr = i30
                corr12(1:icorr) = corr30(1:i30)
                deallocate(corr30)
            else
                imax = 0
            endif
        endif

        ! time lag index and actual time lag in s
        rlag = imax
        tlag = real(imax, kind=dp) / ac_freq

        if (present(ncorr)) ncorr = icorr

        if (present(crosscorr)) then
            crosscorr(1:icorr)  = corr12(1:icorr)
            crosscorr(icorr+1:) = 0.0_dp
        endif

    end subroutine fftmax

    ! ------------------------------------------------------------------

    function mean(dat, mask)
        ! sum(x)/n
        implicit none

        real(dp), dimension(:),           intent(IN)  :: dat
        logical,  dimension(:), optional, intent(IN)  :: mask
        real(dp)                                      :: mean

        real(dp) :: n

        logical,  dimension(size(dat)) :: maske

        maske(:)    = .true.
        if (present(mask)) then
            if (size(mask) /= size(dat)) stop 'Error mean: size(mask) /= size(dat)'
            maske = mask
        endif
        n  = real(count(maske),dp)
        if (n <= (1.0_dp+tiny(1.0_dp))) stop 'mean: n must be at least 2'

        ! Mean
        mean  = sum(dat(:), mask=maske)/n

    end function mean

end module mo_fftmax
