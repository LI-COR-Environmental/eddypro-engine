!***************************************************************************
! write_out_stats.f90
! -------------------
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
! \brief       Write statistics of sensitive variables on output
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutStats(unt, LocStats, string, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    type(StatsType), intent(in) :: LocStats
    integer, intent(in) :: N
    character(*), intent(in) :: string
    !local variables
    integer :: j = 0
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum = ''
    include '../src_common/interfaces.inc'


    call clearstr(dataline)
    !> add file info
    call AddDatum(dataline, string(1:len_trim(string)), separator)
    call WriteDatumInt(N, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> add mean values
    do j = u, pe
        if (j == ti1 .or. j == ti2) cycle
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Mean(j), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
    end do
    !> add wind direction
    call WriteDatumFloat(LocStats%wind_dir, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> add variances
    do j = u, pe
        if (j == ti1 .or. j == ti2) cycle
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Cov(j, j), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
    end do

    !> add u-covariances
    if (E2Col(u)%present .and. E2Col(v)%present) then
        call WriteDatumFloat(LocStats%Cov(u, v), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(w)%present) then
        call WriteDatumFloat(LocStats%Cov(u, w), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(ts)%present) then
        call WriteDatumFloat(LocStats%Cov(u, ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(co2)%present) then
        call WriteDatumFloat(LocStats%Cov(u, co2), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(h2o)%present) then
        call WriteDatumFloat(LocStats%Cov(u, h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(ch4)%present) then
        call WriteDatumFloat(LocStats%Cov(u, ch4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(gas4)%present) then
        call WriteDatumFloat(LocStats%Cov(u, gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if

    !> add remaining v-covariances
    if (E2Col(v)%present .and. E2Col(w)%present) then
        call WriteDatumFloat(LocStats%Cov(v, w), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(v)%present .and. E2Col(ts)%present) then
        call WriteDatumFloat(LocStats%Cov(v, ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(v)%present .and. E2Col(co2)%present) then
        call WriteDatumFloat(LocStats%Cov(v, co2), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(v)%present .and. E2Col(h2o)%present) then
        call WriteDatumFloat(LocStats%Cov(v, h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(v)%present .and. E2Col(ch4)%present) then
        call WriteDatumFloat(LocStats%Cov(v, ch4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(v)%present .and. E2Col(gas4)%present) then
        call WriteDatumFloat(LocStats%Cov(v, gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if

    !> add remaining w-covariances
    if (E2Col(w)%present .and. E2Col(ts)%present) then
        call WriteDatumFloat(LocStats%Cov(w, ts), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(co2)%present) then
        call WriteDatumFloat(LocStats%Cov(w, co2), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(h2o)%present) then
        call WriteDatumFloat(LocStats%Cov(w, h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(ch4)%present) then
        call WriteDatumFloat(LocStats%Cov(w, ch4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(gas4)%present) then
        call WriteDatumFloat(LocStats%Cov(w, gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(tc)%present) then
        call WriteDatumFloat(LocStats%Cov(w, tc), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(pi)%present) then
        call WriteDatumFloat(LocStats%Cov(w, pi), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(te)%present) then
        call WriteDatumFloat(LocStats%Cov(w, te), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(pe)%present) then
        call WriteDatumFloat(LocStats%Cov(w, pe), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
    end if

    !> add standard deviations
    do j = u, pe
        if (j == ti1 .or. j == ti2) cycle
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%StDev(j), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
    end do
    !> add skewness and kurtosis
    do j = u, pe
        if (j == ti1 .or. j == ti2) cycle
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Skw(j), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
    end do
    do j = u, pe
        if (j == ti1 .or. j == ti2) cycle
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Kur(j), datum, EddyProProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyProProj%err_label(1:len_trim(EddyProProj%err_label)), separator)
        end if
    end do
    write(unt, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutStats
