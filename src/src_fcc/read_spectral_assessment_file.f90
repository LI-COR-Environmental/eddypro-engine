!***************************************************************************
! read_spectral_assessment_file.f90
! ---------------------------------
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
! \brief       Read file containing results of spectral assessment to \n
!              retrieve cutoff frequencies
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadSpectralAssessmentFile()
    use m_fx_global_var
    implicit none
    !> local variables
    integer :: gas
    integer :: cls
    integer :: i
    integer :: open_status
    character(ShortInstringLen) :: dataline


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading spectral assessment file: '
    write(*,'(a)') '  ' // AuxFile%sa(1:len_trim(AuxFile%sa))
    open(udf, file = AuxFile%sa, status = 'old', iostat = open_status)

    RegPar%Fn = 0d0
    RegPar%fc = 0d0
    RegPar%f2 = 0d0
    if (open_status == 0) then
        write(*, '(a)') '  Spectral assessment file found, importing content..'
        !> skip 7 lines
        do i = 1, 7
            read(udf, *)
        end do
        !> Read H2O transfer functions for IIR filter
        do cls = RH10, RH90
            read(udf, '(a)') dataline
            dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
            read(dataline, *)  RegPar(h2o, cls)%Fn, RegPar(h2o, cls)%fc
        end do

        !> Read CO2, CH4 and GAS4 transfer functions for IIR filter
        do gas = co2, gas4
            if (gas == h2o) cycle
            !> Skip 2 lines
            read(udf, *)
            read(udf, *)
            do cls = JAN, DEC
                read(udf, '(a)') dataline
                dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
                read(dataline, *)  RegPar(gas, cls)%Fn, RegPar(gas, cls)%fc
            end do
        end do

        !> skip 4 lines
        do i = 1, 4
            read(udf, *)
        end do

        !> Read parameters of exponential fit fc vs. RH
        read(udf, *) RegPar(dum, dum)%e1, RegPar(dum, dum)%e2, RegPar(dum, dum)%e3

        !> skip 6 lines
        do i = 1, 6
            read(udf, *)
        end do
        !> Read parameters of Ibrom's model for spectral correction factor
        read(udf, '(a)') dataline
        dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
        read(dataline, *)  UnPar(1), UnPar(2)
        read(udf, '(a)') dataline
        dataline = dataline(index(dataline, '=') + 1: len_trim(dataline))
        read(dataline, *)  StPar(1), StPar(2)
        close(udf)
        write(*, '(a)') ' Done.'
    else
        !> If the specified file was not found or is empty,
        !> switches to an analytic method
        EddyProProj%hf_meth = 'moncrieff_97'
        call ExceptionHandler(65)
    end if
end subroutine ReadSpectralAssessmentFile
