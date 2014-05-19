!***************************************************************************
! error_handle.f90
! ----------------
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
! \brief       Manages error and Warning instances, possibly aborting execution
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ErrorHandle(prj, lev1, lev2)
    implicit none
    !> in/out variables
    integer, intent(in) :: prj
    integer, intent(in) :: lev1
    integer, intent(in) :: lev2
    !> local variables
    integer :: error_code

    error_code = prj*10000 + lev1*100 + lev2
    select case (error_code)
        !> Errors in common subroutines (0)
        case(0)
            write(*, *) ' Fatal error(0): An error occurred while retrieving user home path.'
            write(*, *) '                 EddyPro is not able to locate project files and must terminate.'
            write(*, *) '                 *Program execution aborted*'
            stop 1
        case(1)
            write(*, *) ' Fatal error(1): temporary file "flist.tmp" not created.'
            write(*, *) '                 most likely, no files with the selected &
                        &extension were found in the data folder.'
            write(*, *) '                 *Program execution aborted*'
            stop 1
        case(101)
            write(*, *) ' Error(201): an error occurred while reading dynamic metadata file. file not found or empty.'
            write(*, *) '             dynamic metadata not used.'
        case(2)
            write(*, *) ' Error(2): an error occurred while reading slow data file. file not found or empty.'
            write(*, *) '           biomet data of this file not used.'
        case(3)
            write(*, *) ' Error(3): error while reading metadata file in GHG archive. file not found or empty.'
            write(*, *) '           GHG archive appears corrupted or invalid.'
            write(*, *) '           GHG file skipped.'
        case(4)
            write(*, *) ' Error(4): error while reading raw data file in GHG archive. file not found or empty.'
            write(*, *) '           GHG archive appears corrupted or invalid.'
            write(*, *) '           GHG file skipped.'
        case(5)
            write(*, *) ' Error(5): error while reading biomet data file or corresponding metadata file. &
                        &file(s) not found or empty.'
            write(*, *) '           biomet data not used for this time period.'
        case(6)
            write(*, *) ' Error(6): error while opening current data file. file is empty.'
            write(*, *) '           file skipped'
        case(7)
            write(*, *)
            write(*, *) ' Error(7): error while opening INI-format file. looking for a solution..'
        case(8)
            write(*, *) ' Fatal error(8): temporary file "flist.tmp" is empty.'
            write(*, *) '                 no files with selected extension were &
                        &found in the raw data folder.'
            write(*, *) '                 *Program execution aborted*'
            stop 1
        case(9)
            write(*, *) ' Fatal error(9): error while opening statistics file.'
            write(*, *) '                 *Program execution aborted*'
            stop 1
        case(10)
            write(*, *) ' Fatal error(10): error while reading statistics file. &
                        &No valid data lines found in file.'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(11)
            write(*, *) ' Fatal error(11): no statistics data were found for the &
                        &selected time period.'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(12)
            write(*, *) ' Fatal error(12): error while opening provisional fluxes file.'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(13)
            write(*, *) ' Error(13): error while reading timelag optimization file. file not found or empty.'
            write(*, *) '            switching to "covariance maximization with default" method.'
        case(14)
            write(*, *) ' Error(14): error while unzipping raw archive. File skipped.'
        case(15)
            write(*, *) ' Warning(15): statistics file is too long. No more than 30000 stats can be&
                        &processed. Check results: most likely the last period will not be processed.'
        case(16)
            write(*, *) ' Fatal error(16): the Raw data directory does not containg any file&
                        & that matches the File name prototype.'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(17)
            write(*, *) '  Error(17): error while opening w/h2o covariances file.'
            write(*, *) '             EddyPro will use maximized w/h2o covariances instead. execution continues.'
        case(18)
            write(*, *) '  Error(18): error while reading w/h2o file. No valid data lines&
                        & were found in the file. EddyPro will use maximum w/h2o covariances instead. execution continues.'
        case(19)
            write(*, *) '  Error(19): number of covariances found in file larger than max supported (18000).&
            & at least one was not imported.'
        case(20)
            write(*, *) '  Fatal error(20): incorrect prototype for raw file names. please check "Raw file name format" entry'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(21)
            write(*, *) '  Fatal error(21): processing project file not found.'
            write(*, *) '                   *Program execution aborted*'
            stop 1
        case(22)
            write(*, *) '  Fatal error(22): alternative metadata file not found.'
            write(*, *) '                   *Program execution aborted*'
            stop 1
        case(23)
            write(*, *) '  Fatal error(23): error validating alternative metadata file.'
            write(*, *) '                   please check entries in the "Metadata file editor" and try again.'
            write(*, *) '                   *Program execution aborted*'
            stop 1
        case(24)
            write(*, *) '  Error(24): an error occurred while opening GHG file.'
            write(*, *) '             skipping to next file.'
        case(25)
            write(*, *) '  Error(25): an error occurred while validating embedded metadata file.'
            write(*, *) '             skipping to next file.'
        case(26)
            write(*, *) '  Error(26): an error occurred while opening ENE file.'
            write(*, *) '             skipping to next file.'
        case(27)
            write(*, *) '  Error(27): an error occurred while validating embedded INI file.'
            write(*, *) '             skipping to next file.'
        case(28)
            write(*, *) '  Error(28): an error occurred while opening raw file.'
            write(*, *) '             skipping to next file.'
        case(29)
            write(*, *) ' Error(29): error while reading planar-fit rotation matrices from file.'
            write(*, *) '            axis rotation method switched to double rotations.'
        case(30)
            write(*, *) ' Error(30): error while opening auxiliary file for planar-fit. &
                        &file not found or empty.'
            write(*, *) '            axis rotation method switched to double rotations.'
        case(31)
            write(*, *) ' Fatal error(31): TOB1 data format unrecognized. Please select either IEEE4 or FP2 data formats'
            write(*, *) '                  and run EddyPro again.'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(32)
            write(*, *) ' Fatal error(32): No valid GHG file was found in the selected data folder. The problem can &
                & be due to corrupted archives or invalid metadata files.'
            write(*, *) '                  *Program execution aborted*'
            stop 1
        case(33)
            write(*, *) ' Error(33): Number of wind records for this sector is less than requested.'
            write(*, *) '            Planar fit results set to -9999 for this sector.'
        case(34)
            write(*, *) ' Error(34): A problem occurred while calculating planar fit for this sector. singular matrix found.'
            write(*, *) '            Planar fit results set to -9999 for this sector.'
        case(35)
            write(*, *) ' Fatal error(35): Sorry, something went wrong. EddyPro was not able to process any raw file.'
            write(*, *) '                  Run will be closed without creating result files.'
            stop 1
        case(36)
            write(*, *) ' Fatal error(36): No output directory was selected.&
                &Select an output directory before running EddyPro.'
            write(*, *) '                       *Program execution aborted*'
            stop 1
        case(37)
            write(*, *) ' Error(37): No wind sector with valid planar fit. Switching to "double rotations"  method.'
        case(38)
            write(*, *) ' Warning(38): 0 sectors selected for planar fit. Forcing to 1 sector.'
        case(39)
            write(*, *) ' Error(39): Error while opening auxiliary file for time lag optimization. &
                        &file not found or empty.'
            write(*, *) '            time lag detection method switched to covariance maximization.'
        case(40)
            write(*, *) ' Error(40): An error occurred while retrieving wind sectors configuration &
                        &for planar fit. Forcing to 1 sector.'
        case(41)
            write(*, *) ' Warning(41): Sector excluded by user.'
        case(42)
            write(*, *) ' Error(42): Method for random uncertainty estimation not recognized. Random uncertainty not calculated.'
        case(43)
            write(*, *) ' Error(43): Time lag optimization failed. Switching to&
                & covariance maximization method for time lag detection.'
        case(44)
            write(*, *) ' Error(44): An error occurred while reading/interpreting biomet file. &
                & Biomet data not used for this period.'
        case(45)
            write(*, *) ' Warning(45): Not enough valid cospectra were found for fitting models, or fitting procedure failed.'
            write(*, *) '              Cospectra output file not created.'

        case(46)
            write(*, *) '  Fatal error(46): The dataset does not contain any raw file corresponding &
                &to the selected sub-period. Select another sub-period or another Raw data directory, or uncheck the option &
                &"Select subperiod" in the Dataset selection page.'
            write(*, *) '                   *Program execution aborted*'
            stop 1

        case(47)
            write(*, *) '  Fatal error(47): Length of raw files appears to be zero or less. Please review metadata &
                &file through the Metadata File Editor'
            write(*, *) '                   *Program execution aborted*'
            stop 1

        case(48)
            write(*, *) '  Fatal error(48): The dataset does not contain any raw file corresponding &
                &to the sub-period selected for planar fit. Select another sub-period &
                &or another Raw data directory, or uncheck the option &
                &"Select subperiod" in the planar fit dialogue.'
            write(*, *) '                   *Program execution aborted*'
            stop 1

        case(49)
            write(*, *) '  Fatal error(49): the dataset does not contain any raw file corresponding &
                &to the sub-period selected for time lag optimization. Select another sub-period &
                &or another Raw data directory, or uncheck the option &
                &"Select subperiod" in the time lag optimization dialogue.'
            write(*, *) '                   *Program execution aborted*'
            stop 1

        case(50)
            write(*, *) '  Fatal error(50): essential files does not contain any results corresponding &
                &to the selected sub-period. Select another sub-period, or uncheck the option &
                &"Select subperiod" in the Dataset selection page.'
            write(*, *) '                   *Program execution aborted*'
            stop 1

        case(51)
            write(*, *) '  Error(51): calculation of tube attenuation parameter (lambda) failed for at least one variable.'
            write(*, *) '             Processing continues but tube attenuation is not included in the spectral correction.'

        case(52)
            write(*, *) '  Fatal error(52): an internal error has occurred.'
            write(*, *) '                   *Program execution aborted*'
            stop 1

        case(53)
            write(*, *) '  Warning(53): no raw data file relevant to current averaging period was found.'
            write(*, *) '               Skipping to next averaging period.'



        case(10004)
            write(*, *) '  Error(10004): error while opening native binary file.'
            write(*, *) '                file skipped.'
        case(10005)
            write(*, *) '  Error(10005): error while reading native binary file.'
            write(*, *) '                file skipped.'
        case(10006)
            write(*, *) '  Error(10006): error while opening native ASCII file.'
            write(*, *) '                file skipped.'
        case(10007)
            write(*, *) '  Error(10007): error while reading file. dataset not created.'
        case(10008)
            write(*, *) ' Warning(10008): available samples not enough for an averaging period.'
            write(*, *) '                 Skipping to next averaging period.'
        case(10009)
            write(*, *) '  Error(10009): at least one wind component seems to be corrupted (too many implausible values).'
            write(*, *) '                Skipping to next averaging period.'

!*******************************************************************************************************************************************
!*******************************************************************************************************************************************

        !> Errors in FCC
        case(20001)
            write(*, *) ' Fatal error(20001): An error occurred while opening "essentials" file.'
            write(*, *) '                     *Program execution aborted*'
            stop 1
        case(20011)
            write(*, *) ' Fatal error(20011): No valid result records found in the "essentials" file.'
            write(*, *) '                     *Program execution aborted*'
            stop 1




        case(20002)
            write(*, *) ' Error(20002): An error occurred while reading binned (co)spectra file. File skipped.'
        case(200014)
            write(*, *) ' Error(20014): An error occurred while reading full (co)spectra file. File skipped.'
        case(20003)
            write(*, *) ' Error(20003): An error occurred while creating output file. Some spectral assessment results&
                & will not be written on output file.'
            write(*, *) '               Processing continues..'
        case(20004)
            write(*, *) ' Error(20004): An error occurred while reading spectral assessment file.'
            write(*, *) '               High-frequency spectral correction method switched to Moncrieff et al. (1997).'
        case(20005)
            write(*, *) ' Fatal error(20005): No valid records found in the "essentials" file.'
            write(*, *) '                      *Program execution aborted*'
            stop 1
        case(20006)
            write(*, *) ' Error(20006): An error occurred while reading spectral assessment file.'
            write(*, *) '               Switching to analytic spectral corrections (Moncrieff et al., 1997).'


        case(20007)
            write(*, *) ' Error(20007): unknown filter shape (neither "iir" nor &
                        &sigma") detected in lptf file.'
            write(*, *) '               lptf parameters set to -9999.'
        case(20008)
            write(*, *) ' Error(20008): error while opening fc=fc(RH) fit parameters file.'
            write(*, *) '               parameters set to -9999.'
        case(20009)
            write(*, *) ' Error(20009): error while creating file. fitting results &
                        &not reported on output.'
        case(20010)
            write(*, *)
            write(*, *) '  Warning(20010): footprint estimation not requested.'
        case(20012)
            write(*, *) ' Error(20012): error while opening file. quality flags not calculated.'
        case(20013)
            write(*, *) '  Warning(20013): quality flag calculation not requested.'


        !> Errors in SpectralAnalysis(4)
        case(40001)
            write(*, *) ' Fatal error(40001): No meteo data, either measured or estimated, &
                        &was imported. SpectralAnalysis cannot run.'
            write(*, *) '                     *Program execution aborted*'
            stop 1
        case(40002)
            write(*, *) ' Fatal error(40002): No provisional flux data was imported. &
                        &SpectralAnalysis cannot run.'
            write(*, *) '                     *Program execution aborted*'
            stop 1
        case(40004)
            write(*, *) ' Error(40004): error while creating file. output file not created.'

!********************************************************************************

        !> Obsolete, to be updated
        case(10901)
            write(*, *) ' Warning(10901): Not all variables (u,v,w,ts,co2,h2o) needed'
            write(*, *) '                 in ECCOCE are present in the current data file.'
            write(*, *) '                 File skipped and results set to -999.9.'
        case(10902)
            write(*, *) ' Warning(10902): Spikes number exceed 1% for at least one variable'
            write(*, *) '                 File skipped and results set to -999.9.'
        case(11001)
            write(*, *) ' Warning(11001): Error while opening auxiliary file for'
            write(*, *) '                 transfer function parameters.'
            write(*, *) '                 File not found or empty.'
            write(*, *) '                 Switching to analytic spectral correction.'
        case(11002)
            write(*, *) ' Warning(11002): At least one transfer function parameter results <=0.'
            write(*, *) '                 Experimental spectral correction not applied.'
            write(*, *) '                 Switching to analytic spectral correction.'
        case(12601)
            write(*, *) ' Warning(12601): Either air pressure or temperature is <= 0.'
            stop        '                 Program execution aborted'
        case(13001)
            write(*, *) ' Warning(13001): Impossible to run footprint model for this file.'
            write(*, *) '                 Footprint results set to -9999.'

        case(3003)
            write(*, *) 'Warning(3003): Anemometer not supported yet.'
        case(3004)
            write(*, *) 'Warning(3004): No correct data line detected in the ".raw" file'
            write(*, *) '               Check anemometer status or file ".ini" for a'
            write(*, *) '               proper setup.'
            write(*, *) '               Acquisition proceeds without field-level processing'
            write(*, *) '               and without writing ".ene" files.'
            stop        '               Only native anemometer streams will be saved.'

        case(40601)
            write(*, *) ' Warning(40601): Current data file too short or too many bad data lines.'
            write(*, *) '                 File skipped and results set to -999.9.'
        case(40501)
            write(*, *) ' Warning(40501): Something went wrong. Detected file duration is less than zero.'
            write(*, *) '                 Program execution aborted'
            stop 1
        case(40701)
            write(*, *) ' Warning(40701): Not all variables (u,v,w,ts,co2,h2o) needed in EddyPro:spec_corr are present&
                                          &in the current data file.'
            write(*, *) '                 File skipped and results set to -999.9.'

        case(5001)
            write(*, *) 'Warning(5001): Vertical wind data (W) missing. Processing impossible.'
            stop         '              Acquisition proceeds without field-level processing.'
        case(5002)
            write(*, *) 'Warning(5002): CO2 data missing. Processing incomplete.'
            write(*, *) '               Check results in "fieldfluxes.txt".'
        case(5003)
            write(*, *) 'Warning(5003): H2O data missing. Processing incomplete.'
            write(*, *)   '             Check results in "fieldfluxes.txt".'
        case(5004)
            write(*, *) 'Warning(5004): Temperature data (either Ts, SoS or Ta) missing.'
            write(*, *) '               Processing impossible.'
            stop         '              Acquisition proceeds without field-level processing.'
        case(5005)
            write(*, *) 'Warning(5005): Horizontal wind data (U) missing. Processing incomplete.'
            write(*, *)   '             Check results in "fieldfluxes.txt".'
        case(5006)
            write(*, *) 'Warning(5006): Horizontal wind data (V) missing. Processing incomplete.'
            write(*, *) '               Check results in "fieldfluxes.txt".'
        case default
            write(*, *) 'Warning: A general (unrecognized) error occurred.'
            stop         '        Program execution aborted.'
    end select
end subroutine ErrorHandle
