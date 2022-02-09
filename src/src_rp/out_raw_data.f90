!***************************************************************************
! out_raw_data.f90
! ----------------
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
! \brief       Output raw datasets upon request
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine OutRawData(date, time, Set, nrow, ncol, level)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: level
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    character(*), intent(in) :: date
    character(*), intent(in) :: time
    !> local variables
    integer :: i
    integer :: j
    integer :: num_var
    character(1) :: lev_char
    character(13) :: Datestring
    character(PathLen) :: RawPath
    character(256) :: string
    character(256) :: string_utf8
    real(kind = dbl) :: OutSet(nrow, ncol)
    integer :: hlen
    character(512) :: raw_out_corr

    !> If binned (co)spectra or at least one full (co)spectrum are requested,
    !> perform all related calculations
    write(lev_char, '(i1)') level
    write(*, '(a)', advance = 'no') '  Writing raw dataset (level ' // lev_char // ') on output..'

    Datestring = date(1:4) // date(6:7) // date(9:10) &
               // '-' // time(1:2) // time(4:5)

    !> Open output file for binned co-spectra
    RawPath = RawSubDir(level)(1:len_trim(RawSubDir(level))) &
        // Datestring // Raw_FilePadding // Timestamp_FilePadding // TxtExt

    open(udf, file = RawPath, encoding = 'utf-8')
    write(udf, '(a)') 'Raw dataset level ' // lev_char
    write(udf, '(a)') 'Legend:'
    write(udf, '(a)') 'u, v, w: wind components [m s-1]'
    write(udf, '(a)') 'ts: sonic temperature [K]'
    write(udf, '(a)') 'air_t: ambient temperature [K]'
    write(udf, '(a)') 'air_p: ambient pressure [Pa]'
    string =  'co2, ch4, 4th gas: molar density [mmol m-3],&
        & mole fraction [' // char(181) // 'mol/mol] or mixing ratio &
        &[' // char(181) // 'mol/mol], depending on raw data'
    call latin1_to_utf8(string, string_utf8)
    write(udf, '(a)') trim(adjustl(string_utf8))
    write(udf, '(a)') 'h2o: molar density [mmol m-3],&
        & mole fraction [mmol/mol] or mixing ratio [mmol/mol]&
        &, depending on raw data'
    write(udf, '(a)')
    if (level == 8) then
        raw_out_corr = '   '
        hlen = 3
        do j=5, ncol
            if (E2Col(j)%present) then
                raw_out_corr = raw_out_corr(1:hlen) // E2Col(j)%var
                hlen = hlen + 25
            end if
        end do
        write(udf, '(a)') raw_out_corr(1:hlen)
    else
        write(udf, '(a)') raw_out_header(1:len_trim(raw_out_header))
    end if

    !> Define output dataset
    num_var = 0
    if (level == 8) then
        do j = 5, ncol
            if (E2Col(j)%present) then
                num_var = num_var + 1
                OutSet(:, num_var) = Set(:, j)
            end if
        end do
    else
        do j = 1, ncol
            if (RPsetup%out_raw_var(j)) then
                num_var = num_var + 1
                OutSet(:, num_var) = Set(:, j)
            end if
        end do
    end if

    !> write output dataset
    do i = 1, nrow
        write(udf,*) OutSet(i, 1:num_var)
    end do

    close(udf)
    write(*, '(a)') ' Done.'

end subroutine OutRawData
