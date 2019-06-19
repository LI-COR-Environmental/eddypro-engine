!***************************************************************************
! writeout_qc_details.f90
! -----------------------
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
! \brief       Write results of stationarity and integral turbulence test on
!              file as requested
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutQCDetails(string, StDiff, DtDiff, STFlg, DTFlg)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: string
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    integer, intent(in) :: STFlg(GHGNumVar)
    integer, intent(in) :: DTFlg(GHGNumVar)
    !> local variables
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum

    call clearstr(dataline)

    !> File info
    call AddDatum(dataline, string(1:len_trim(string)), separator)
    !> Differences (%) of stationarity test
    call WriteDatumInt(STDiff%u, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%w, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (OutVarPresent(co2)) then
        call WriteDatumInt(STDiff%co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(h2o)) then
        call WriteDatumInt(STDiff%h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(ch4)) then
        call WriteDatumInt(STDiff%ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(gas4)) then
        call WriteDatumInt(STDiff%gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    call WriteDatumInt(STDiff%w_u, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STDiff%w_ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    if (OutVarPresent(co2)) then
        call WriteDatumInt(STDiff%w_co2, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(h2o)) then
        call WriteDatumInt(STDiff%w_h2o, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(ch4)) then
        call WriteDatumInt(STDiff%w_ch4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(gas4)) then
        call WriteDatumInt(STDiff%w_gas4, datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    !> Flags for the stationarity test
    call WriteDatumInt(STFlg(w_u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(STFlg(w_ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    if (OutVarPresent(co2)) then
        call WriteDatumInt(STFlg(w_co2), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(h2o)) then
        call WriteDatumInt(STFlg(w_h2o), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(ch4)) then
        call WriteDatumInt(STFlg(w_ch4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if
    if (OutVarPresent(gas4)) then
        call WriteDatumInt(STFlg(w_gas4), datum, EddyProProj%err_label)
        call AddDatum(dataline, datum, separator)
    end if

    !> Differences (%) of the well developed turbulence test
    call WriteDatumInt(DtDiff%u, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%w, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DtDiff%ts, datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    !> Flags for the well developed turbulence test
    call WriteDatumInt(DTFlg(u), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(w), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)
    call WriteDatumInt(DTFlg(ts), datum, EddyProProj%err_label)
    call AddDatum(dataline, datum, separator)

    write(uqc, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutQCDetails
