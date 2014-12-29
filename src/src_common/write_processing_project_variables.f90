!***************************************************************************
! write_processing_project_variables.f90
! --------------------------------------
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
! \brief       Read EddyPro configuration file in its preliminary
!              section, common to all sub-programs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteProcessingProjectVariables()
    use m_common_global_var
    implicit none
    !> local variables
    integer :: dot

    !> Initializations
    Auxfile%metadata   = 'none'
    Auxfile%biomet     = 'none'
    Dir%biomet         = 'none'
    Dir%main_out       = 'none'
    EddyProProj%fproto = 'none'

    !> Project general info
    select case (EPPrjCTags(16)%value(1:1))
        case('1')
            EddyProProj%run_mode =  'express'
        case('2')
            EddyProProj%run_mode =  'md_retrieval'
        case default
            EddyProProj%run_mode =  'advanced'
    end select

    EddyProProj%title  = EPPrjCTags(4)%value(1:len_trim(EPPrjCTags(4)%value))
    EddyProProj%id     = EPPrjCTags(5)%value(1:len_trim(EPPrjCTags(5)%value))
    if (EddyProProj%id(1:1) /= '_') then
        EddyProProj%id = 'eddypro_' // EddyProProj%id(1:len_trim(EddyProProj%id))
    else
        EddyProProj%id = 'eddypro_' // EddyProProj%id(1:len_trim(EddyProProj%id))
    end if

    !>  file type
    select case (EPPrjCTags(6)%value(1:1))
        case ('0')
        EddyProProj%ftype = 'licor_ghg'
        EddyProProj%fext = 'ghg'
        EddyProLog%iso_format = .true.
        EddyProProj%fproto = 'yyyy-mm-ddTHHMM'
        case ('1')
        EddyProProj%ftype = 'generic_ascii'
        case ('2')
        EddyProProj%ftype = 'tob1'
        case ('3')
        EddyProProj%ftype = 'eddymeas_bin'
        case ('4')
        EddyProProj%ftype = 'edisol_bin'
        case ('5')
        EddyProProj%ftype = 'generic_bin'
        case ('6')
        EddyProProj%ftype = 'alteddy_bin'
    end select

    !> If file type is different from GHG, metadata retrieval mode is not feasible
    !> so forces into advanced mode
    if (EddyProProj%ftype /= 'licor_ghg' .and. EddyProProj%run_mode ==  'md_retrieval') &
        EddyProProj%run_mode =  'advanced'

    !> File names prototype and related (for non-ENE, non-LICOR files)
    if (EddyProProj%ftype /= 'licor_ghg') then
        EddyProProj%fproto = EPPrjCTags(7)%value(1:len_trim(EPPrjCTags(7)%value))
        if (index(EddyProProj%fproto, '.') /= 0) then
            !> File extensions
            dot = index(EddyProProj%fproto, '.', .true.)
            EddyProProj%fext = EddyProProj%fproto(dot + 1:len_trim(EddyProProj%fproto))
            !> ISO format
            EddyProLog%iso_format = index(EddyProProj%fproto, 'mm') /= 0
        end if
    end if

    !> If file type is TOB1, check if user entered the format
    FileInterpreter%tob1_format = 'none'
    if (EddyProProj%ftype == 'tob1') then
        select case (EPPrjCTags(32)%value(1:len_trim(EPPrjCTags(32)%value)))
            case ('1')
                FileInterpreter%tob1_format = 'IEEE4'
            case ('2')
                FileInterpreter%tob1_format = 'FP2'
        end select
    end if

    !> Whether to use alternative metadata file
    EddyProProj%use_extmd_file = EPPrjCTags(9)%value(1:1) == '1'
    AuxFile%metadata ='none'
    if(EddyProProj%use_extmd_file) &
        AuxFile%metadata = EPPrjCTags(10)%value(1: len_trim(EPPrjCTags(10)%value))

    !> Whether to use dynamic metadata file
    EddyProProj%use_dynmd_file = EPPrjCTags(11)%value(1:1) == '1'
    if(EddyProProj%use_dynmd_file) &
        AuxFile%DynMD = EPPrjCTags(12)%value(1: len_trim(EPPrjCTags(12)%value))

    !> Settings for binary raw files
    if (EddyProProj%ftype(1:len_trim(EddyProProj%ftype)) == 'generic_bin') then
        !> select line terminator in ASCII header of binary files
        select case(EPPrjCTags(13)%value(1:1))
            case ('0')
            Binary%ascii_head_eol = 'cr/lf'
            case ('1')
            Binary%ascii_head_eol = 'lf'
            case ('2')
            Binary%ascii_head_eol = 'cr'
        end select
        !> select binary files endianess
        Binary%little_endian = EPPrjCTags(14)%value(1:1) == '1'
        !> Select number of bytes per variable
        Binary%nbytes = nint(EPPrjNTags(1)%value)
        !> Select number of ASCII header lines
        Binary%head_nlines = nint(EPPrjNTags(2)%value)
    end if

    !> Master sonic
    EddyProProj%master_sonic = EPPrjCTags(15)%value(1: len_trim(EPPrjCTags(15)%value))
    !> Variables to be used other than sonic ones
    EddyProProj%col(ts:pe) = nint(error)
    EddyProProj%col(ts)  = nint(EPPrjNTags(3)%value)
    EddyProProj%col(co2) = nint(EPPrjNTags(4)%value)
    EddyProProj%col(h2o) = nint(EPPrjNTags(5)%value)
    EddyProProj%col(ch4) = nint(EPPrjNTags(6)%value)
    EddyProProj%col(gas4) = nint(EPPrjNTags(7)%value)
    EddyProProj%col(tc)  = nint(EPPrjNTags(8)%value)
    EddyProProj%col(ti1) = nint(EPPrjNTags(9)%value)
    EddyProProj%col(ti2) = nint(EPPrjNTags(10)%value)
    EddyProProj%col(pi)  = nint(EPPrjNTags(11)%value)
    EddyProProj%col(te)  = nint(EPPrjNTags(12)%value)
    EddyProProj%col(pe)  = nint(EPPrjNTags(13)%value)
    EddyProProj%col(E2NumVar + diag72) = nint(EPPrjNTags(14)%value)
    EddyProProj%col(E2NumVar + diag75) = nint(EPPrjNTags(15)%value)
    EddyProProj%col(E2NumVar + diag77) = nint(EPPrjNTags(16)%value)

    !> if a column was selected for gas4, read diffusivity. If diffusivity is
    !> below zero, defaults to gas4 diffusivity
    if (EddyProProj%col(gas4) > 0) then
        Dc(gas4) = EPPrjNTags(17)%value * 1d-4 !< takes from cm+2s-1 to m+2s-1
        if (Dc(gas4) <= 0) Dc(gas4) = 0.00001436d0  !< default for N2O from Massman (1998, J. Atm. Env) Table 2.
        MW(gas4) = sngl(EPPrjNTags(18)%value) * 1e-3 !< takes from g+1mol-1 to kg+1mol-1
        if (MW(gas4) <= 0) MW(gas4) = 44.01e-3  !< default for N2O
    end if

    !> biomet measurements info
    select case (EPPrjCTags(17)%value(1:1))
        case ('1')
        EddyProProj%biomet_data = 'embedded'
        case ('2')
        EddyProProj%biomet_data = 'ext_file'
        case ('3')
        EddyProProj%biomet_data = 'ext_dir'
        case default
        EddyProProj%biomet_data = 'none'
    end select
    !> biomet files/folders as applicable
    if (EddyProProj%biomet_data == 'ext_file') &
        AuxFile%biomet = EPPrjCTags(18)%value(1: len_trim(EPPrjCTags(18)%value))
    if (EddyProProj%biomet_data == 'ext_dir') then
        Dir%biomet = EPPrjCTags(29)%value(1: len_trim(EPPrjCTags(29)%value))
        if (len_trim(Dir%biomet) == 0) then
            Dir%biomet = 'none'
        else
            EddyProProj%biomet_tail = EPPrjCTags(30)%value(1: len_trim(EPPrjCTags(30)%value))
            EddyProProj%biomet_recurse = EPPrjCTags(31)%value(1:1) == '1'
        end if
    end if

    !> select whether to output GHG-europe-formatted file
    EddyProProj%out_ghg_eu = EPPrjCTags(19)%value(1:1) == '1'
    !> select whether to output AmeriFlux-formatted file
    EddyProProj%out_amflux = EPPrjCTags(20)%value(1:1) == '1'
    !> select whether to output full output file
    EddyProProj%out_full = EPPrjCTags(21)%value(1:1) == '1'
    !> select whether to use fixed or dynamic output format
    EddyProProj%out_md = EPPrjCTags(39)%value(1:1) == '1'
    !> select whether to output average cospectra
    EddyProProj%out_avrg_cosp = EPPrjCTags(41)%value(1:1) == '1'
    !> select whether to output average spectra
    EddyProProj%out_avrg_spec = EPPrjCTags(43)%value(1:1) == '1'
    !> select whether to output biomet average values
    EddyProProj%out_biomet = EPPrjCTags(42)%value(1:1) == '1'
    !> select whether to use fixed or dynamic output format
    EddyProProj%fix_out_format = EPPrjCTags(37)%value(1:1) == '1'

    !> Select whether to apply high-pass theoretical spectral correction.
    !> It is independent from the choice of the low-pass method
    select case (EPPrjCTags(22)%value(1:1))
        case ('0')
            !> Do not apply low-frequency spectral correction
            EddyProProj%lf_meth = 'none'
        case ('1')
            EddyProProj%lf_meth = 'analytic'
    end select

    !> Select low-pass spectral correction method.
    select case (EPPrjCTags(23)%value(1:1))
        case ('0')
            !> Do not apply spectral correction (e.g. open-path)
            EddyProProj%hf_meth = 'none'
        case ('1')
            !> Correction after Moncrieff et al (1997, JH) fully analytical
            EddyProProj%hf_meth = 'moncrieff_97'
        case ('2')
            !> Correction after Horst (1997, BLM), in-situ/analytical
            EddyProProj%hf_meth = 'horst_97'
        case ('3')
            !> Correction after Ibrom et al (2007, AFM) fully in-situ
            EddyProProj%hf_meth = 'ibrom_07'
        case ('4')
            !> Correction after Fratini et al. 2010, fully in-situ
            EddyProProj%hf_meth = 'fratini_12'
        case ('5')
            !> Correction after Massman (2000, 2001), fully analytical
            EddyProProj%hf_meth = 'massman_00'
        case ('6')
            !> Custom correction, in-situ/analytical
            EddyProProj%hf_meth = 'custom'
        case default
            !> If not specified, set to none
            EddyProProj%hf_meth = 'none'
    end select

    !> select whether to fill gaps with error codes
    EddyProProj%make_dataset = EPPrjCTags(24)%value(1:1) == '1'

    !> start/end date and time of period to be processed
    if (EPPrjCTags(40)%value(1:1) == '1') then
        EddYProProj%start_date = EPPrjCTags(25)%value(1:len_trim(EPPrjCTags(25)%value))
        if (len_trim(EddYProProj%start_date) == 0) EddYProProj%start_date = 'none'
        EddYProProj%start_time = EPPrjCTags(26)%value(1:len_trim(EPPrjCTags(26)%value))
        if (len_trim(EddYProProj%start_time) == 0) EddYProProj%start_time = 'none'
        EddYProProj%end_date = EPPrjCTags(27)%value(1:len_trim(EPPrjCTags(27)%value))
        if (len_trim(EddYProProj%end_date) == 0) EddYProProj%end_date = 'none'
        EddYProProj%end_time = EPPrjCTags(28)%value(1:len_trim(EPPrjCTags(28)%value))
        if (len_trim(EddYProProj%end_time) == 0) EddYProProj%end_time = 'none'
    else
        EddYProProj%start_date = '1900-01-01'
        EddYProProj%start_time = '00:00'
        EddYProProj%end_date   = '2100-12-31'
        EddYProProj%end_time   = '23:59'
    end if
    if (EddYProProj%start_date == 'none') EddYProProj%start_date = '1900-01-01'
    if (EddYProProj%start_time == 'none') EddYProProj%start_time = '00:00'
    if (EddYProProj%end_date == 'none')   EddYProProj%end_date   = '2100-12-31'
    if (EddYProProj%end_time == 'none')   EddYProProj%end_time   = '23:59'

    !> select whether to apply WPL correction at all
    EddyProProj%wpl = EPPrjCTags(33)%value(1:1) /= '0'

    !> set error string
    EddyProProj%err_label = EPPrjCTags(36)%value(1:len_trim(EPPrjCTags(36)%value))
    if (len_trim(EddyProProj%err_label) == 0 .or. EddyProProj%err_label == 'none') &
        EddyProProj%err_label = '-9999.0'

    !> select footprint method
    select case (EPPrjCTags(34)%value(1:1))
        case ('0')
        Meth%foot = 'none'
        case ('1')
        Meth%foot = 'kljun_04'
        case ('2')
        Meth%foot = 'kormann_meixner_01'
        case ('3')
        Meth%foot = 'hsieh_00'
        case default
        Meth%foot = 'kljun_04'
    end select

    !> select quality-flagging method
    select case (EPPrjCTags(38)%value(1:1))
        case ('0')
        Meth%qcflag = 'none'
        case ('1')
        Meth%qcflag = 'mauder_foken_04'
        case ('2')
        Meth%qcflag = 'foken_03'
        case ('3')
        Meth%qcflag = 'goeckede_06'
        case default
        Meth%qcflag = 'mauder_foken_04'
    end select

    !> main output directory, only in Desktop mode
    if (EddyProProj%run_env /= 'embedded') then
        Dir%main_out = EPPrjCTags(35)%value
        if (len_trim(Dir%main_out) == 0) then
            write(*, *)
            call ExceptionHandler(36)
        end if
        call AdjDir(Dir%main_out, slash)
    end if

    !> Adjust paths
    call AdjFilePath(AuxFile%metadata, slash)
    call AdjFilePath(AuxFile%biomet, slash)
    call AdjDir(Dir%biomet, slash)
end subroutine WriteProcessingProjectVariables
