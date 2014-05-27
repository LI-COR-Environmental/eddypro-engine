!***************************************************************************
! read_spectral_assessment_file.f90
! ---------------------------------
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
    character(512) :: string


    !> Open planar fit file and read rotation matrices
    write(*,'(a)') ' Reading spectral assessment file: '
    write(*,'(a)') '  ' // AuxFile%sa(1:len_trim(AuxFile%sa))
    open(udf, file = AuxFile%sa, status = 'old', iostat = open_status)

    RegPar%Fn = 0d0
    RegPar%fc = 0d0
    RegPar%f2 = 0d0
    if (open_status == 0) then
        write(*, '(a)', advance = 'no') ' Spectral assessment file found, importing content..'
        !> skip 7 lines
        do i = 1, 7
            read(udf, *)
        end do
        !> Read H2O transfer functions for IIR filter
        do cls = RH10, RH90
            read(udf, '(a)') string
            string = string(index(string, '=') + 1: len_trim(string))
            read(string, *)  RegPar(h2o, cls)%Fn, RegPar(h2o, cls)%fc
        end do

        !> Read CO2, CH4 and GAS4 transfer functions for IIR filter
        do gas = co2, gas4
            if (gas == h2o) cycle
            !> Skip 2 lines
            read(udf, *)
            read(udf, *)
            do cls = JAN, DEC
                read(udf, '(a)') string
                string = string(index(string, '=') + 1: len_trim(string))
                read(string, *)  RegPar(gas, cls)%Fn, RegPar(gas, cls)%fc
            end do
        end do

        !> skip 4 lines
        do i = 1, 4
            read(udf, *)
        end do

        !> Read parameters of exponential fit fc vs. RH
        read(udf, *)  RegPar(dum, dum)%e1, RegPar(dum, dum)%e2, RegPar(dum, dum)%e3

        !> skip 6 lines
        do i = 1, 6
            read(udf, *)
        end do
        !> Read parameters of Ibrom's model for spectral correction factor
        read(udf, '(a)') string
        string = string(index(string, '=') + 1: len_trim(string))
        read(string, *)  UnPar(1), UnPar(2)
        read(udf, '(a)') string
        string = string(index(string, '=') + 1: len_trim(string))
        read(string, *)  StPar(1), StPar(2)
        close(udf)
        write(*,*) ' Done.'
    else
       !> If the specified file was not found or is empty,m switches to an analytic method
        EddyProProj%hf_meth = 'moncrieff_97'
        call ExceptionHandler(65)
    end if
end subroutine ReadSpectralAssessmentFile
