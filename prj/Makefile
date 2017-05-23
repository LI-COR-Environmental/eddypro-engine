#
# Makefile for EddyPro's engine executables
# Author: Gerardo Fratini, LI-COR Biosciences
#

#Portability
ifdef ComSpec
	#Windows
	RM = if exist FILE del /F /Q FILE
	MV = if exist FILE move FILE DESTINATION
	RP_EXE = eddypro_rp.exe
	FCC_EXE = eddypro_fcc.exe
	OS_DIR = win
	FC = mingw32-gfortran.exe
	MAKE_DIR = if not exist DIR mkdir DIR
	PATHSEP2 = \\
	LFLAGS = -static
else
	#Unix
	RM = rm -f FILE
	MV = mv FILE DESTINATION || true
	RP_EXE = eddypro_rp
	FCC_EXE = eddypro_fcc
	FC = gfortran
	MAKE_DIR = @mkdir -p DIR
	PATHSEP2 = /
	platform=$(shell uname)
	ifeq ($(platform), Linux)
		OS_DIR = linux
	else
		OS_DIR = mac
	endif
	LFLAGS =
endif
PATHSEP=$(strip $(PATHSEP2))

#Dirs as seen by Make
##RP
RP_SRC_DIR_f90d1 = ../src/src_rp/
RP_SRC_DIR_Fd2 = ../src/src_rp/fft4/
##FCC
FCC_SRC_DIR_f90d1 = ../src/src_fcc/
##Common
COMMON_SRC_DIR_f90d1 = ../src/src_common/
OBJS_DIR = ../obj/$(OS_DIR)/
EXE_DIR = ../bin/$(OS_DIR)/

#Dirs as seen by OS
##RP
OS_RP_SRC_DIR_f90d1 = ..$(PATHSEP)src$(PATHSEP)src_rp$(PATHSEP)
OS_RP_SRC_DIR_Fd2   = ..$(PATHSEP)src$(PATHSEP)src_rp$(PATHSEP)fft4$(PATHSEP)
##FCC
OS_FCC_SRC_DIR_f90d1 = ..$(PATHSEP)src$(PATHSEP)src_fcc$(PATHSEP)
##Common
OS_COMMON_SRC_DIR_f90d1 = ..$(PATHSEP)src$(PATHSEP)src_common$(PATHSEP)
OS_OBJS_DIR = ..$(PATHSEP)obj$(PATHSEP)$(OS_DIR)$(PATHSEP)
OS_EXE_DIR  = ..$(PATHSEP)bin$(PATHSEP)$(OS_DIR)$(PATHSEP)

#List of source file basenames (.f90, .F)
##RP
RP_SRCS_f90d1 = $(notdir $(wildcard $(RP_SRC_DIR_f90d1)*.f90))
RP_SRCS_Fd2 = $(notdir $(wildcard $(RP_SRC_DIR_Fd2)*.F))
##FCC
FCC_SRCS_f90d1 = $(notdir $(wildcard $(FCC_SRC_DIR_f90d1)*.f90))
##Common
COMMON_SRCS_f90d1 = $(notdir $(wildcard $(COMMON_SRC_DIR_f90d1)*.f90))

#List of object file basenames (.o)
##RP
RP_OBJS_f90d1=$(patsubst %.f90, %.o, $(RP_SRCS_f90d1))
RP_OBJS_Fd2=$(patsubst %.F, %.o, $(RP_SRCS_Fd2))
##FCC
FCC_OBJS_f90d1=$(patsubst %.f90, %.o, $(FCC_SRCS_f90d1))
##Common
COMMON_OBJS_f90d1=$(patsubst %.f90, %.o, $(COMMON_SRCS_f90d1))

#CFLAGS = -ffpe-trap=invalid,zero -Wall -pedantic -fall-intrinsics -fbounds-check -Wextra -Wcharacter-truncation -Wunderflow -Wno-compare-reals -std=f2003 -O3 -I $(OBJS_DIR)
#CFLAGS = -g -fcheck=all -Wall -pedantic -fall-intrinsics -fbounds-check -Wextra -Wcharacter-truncation -Wunderflow -Wno-compare-reals -std=f2003 -O3 -I $(OBJS_DIR)
CFLAGS = -Wall -pedantic -fall-intrinsics -fbounds-check -Wextra -Wcharacter-truncation -Wunderflow -Wno-compare-reals -std=f2003 -O3 -I $(OBJS_DIR)

VPATH = $(RP_SRC_DIR_f90d1):$(RP_SRC_DIR_Fd2):$(FCC_SRC_DIR_f90d1):$(COMMON_SRC_DIR_f90d1):$(OBJS_DIR)
RP_OBJS = $(addprefix $(OBJS_DIR), $(RP_OBJS_f90d1) $(RP_OBJS_Fd2)$(COMMON_OBJS_f90d1))
FCC_OBJS = $(addprefix $(OBJS_DIR), $(FCC_OBJS_f90d1) $(COMMON_OBJS_f90d1))

#Targets
all : rp fcc

rp : $(RP_OBJS_f90d1) $(RP_OBJS_Fd2) $(COMMON_OBJS_f90d1)
	@$(subst DIR,$(OS_EXE_DIR),$(MAKE_DIR))
	$(FC) $(LFLAGS) -o $(EXE_DIR)$(RP_EXE) $(RP_OBJS)
	@$(subst DESTINATION, $(OS_OBJS_DIR), $(subst FILE, *.mod, $(MV)))

fcc : $(FCC_OBJS_f90d1) $(COMMON_OBJS_f90d1)
	@$(subst DIR,$(OS_EXE_DIR),$(MAKE_DIR))
	$(FC) $(LFLAGS) -o $(EXE_DIR)$(FCC_EXE) $(FCC_OBJS)
	@$(subst DESTINATION, $(OS_OBJS_DIR), $(subst FILE, *.mod, $(MV)))

#Compilation
##RP
$(RP_OBJS_f90d1):
	@$(subst DIR,$(OS_OBJS_DIR),$(MAKE_DIR))
	$(FC) $(CFLAGS) -c $(RP_SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@
$(RP_OBJS_Fd2):
	@$(subst DIR,$(OS_OBJS_DIR),$(MAKE_DIR))
	$(FC) -c $(RP_SRC_DIR_Fd2)$(@:.o=.F) -o $(OBJS_DIR)$@
##FCC
$(FCC_OBJS_f90d1):
	@$(subst DIR,$(OS_OBJS_DIR),$(MAKE_DIR))
	$(FC) $(CFLAGS) -c $(FCC_SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@
##Common
$(COMMON_OBJS_f90d1):
	@$(subst DIR,$(OS_OBJS_DIR),$(MAKE_DIR))
	$(FC) $(CFLAGS) -c $(COMMON_SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean_all: clean_rp clean_fcc

clean_rp:
	@$(subst FILE, $(OS_OBJS_DIR)*, $(RM))
	@$(subst FILE, $(OS_EXE_DIR)$(RP_EXE), $(RM))

clean_fcc:
	@$(subst FILE, $(OS_OBJS_DIR)*, $(RM))
	@$(subst FILE, $(OS_EXE_DIR)$(FCC_EXE), $(RM))

#File dependencies
m_dates.o: \
	m_dates.f90
m_libdate.o: \
	m_libdate.f90
m_numeric_kinds.o: \
	m_numeric_kinds.f90
m_fp2_to_float.o: \
	m_fp2_to_float.f90 \
	m_numeric_kinds.o
m_index_parameters.o: \
	m_index_parameters.f90 \
	m_numeric_kinds.o
m_methane_tables.o: \
	m_methane_tables.f90 \
	m_numeric_kinds.o
m_typedef.o: \
	m_typedef.f90 \
	m_libdate.o \
	m_numeric_kinds.o
m_common_global_var.o: \
	m_common_global_var.f90 \
	m_libdate.o \
	m_dates.o \
	m_fp2_to_float.o \
	m_index_parameters.o \
	m_methane_tables.o \
	m_typedef.o
m_rp_global_var.o: \
	m_rp_global_var.f90 \
	m_common_global_var.o
tests_not_performed.o: \
	tests_not_performed.f90 \
	m_rp_global_var.o
first_file_of_current_period.o: \
	first_file_of_current_period.f90 \
	m_common_global_var.o
interpret_diagnostics.o: \
	interpret_diagnostics.f90 \
	m_rp_global_var.o
parse_ini_file.o: \
	parse_ini_file.f90 \
	m_common_global_var.o
basic_stats.o: \
	basic_stats.f90 \
	m_common_global_var.o
timelag_handle.o: \
	timelag_handle.f90 \
	m_rp_global_var.o
count_records_and_values.o: \
	count_records_and_values.f90 \
	m_common_global_var.o
wind_sector.o: \
	wind_sector.f90 \
	m_rp_global_var.o
footprint_handle.o: \
	footprint_handle.f90 \
	m_common_global_var.o
molefractions_and_mixingratios.o: \
	molefractions_and_mixingratios.f90 \
	m_rp_global_var.o
bpcf_Ibrom_07.o: \
	bpcf_Ibrom_07.f90 \
	m_common_global_var.o
planarfit_rotation_matrix.o: \
	planarfit_rotation_matrix.f90 \
	m_common_global_var.o
ts_extract_subperiod_indexes.o: \
	ts_extract_subperiod_indexes.f90 \
	m_common_global_var.o
write_out_biomet.o: \
	write_out_biomet.f90 \
	m_rp_global_var.o
write_icos_output_rp.o: \
	write_icos_output_rp.f90 \
	m_rp_global_var.o
write_out_icos_only_biomet.o: \
	write_out_icos_only_biomet.f90 \
	m_rp_global_var.o
out_raw_data.o: \
	out_raw_data.f90 \
	m_rp_global_var.o
bpcf_analytic_transfer_functions.o: \
	bpcf_analytic_transfer_functions.f90 \
	m_common_global_var.o
portion_of_file_in_current_period.o: \
	portion_of_file_in_current_period.f90 \
	m_common_global_var.o
unzip_archive.o: \
	unzip_archive.f90 \
	m_common_global_var.o
write_out_user_stats.o: \
	write_out_user_stats.f90 \
	m_rp_global_var.o
override_settings.o: \
	override_settings.f90 \
	m_rp_global_var.o
bpcf_aux_subs.o: \
	bpcf_aux_subs.f90 \
	m_common_global_var.o
qc_flags_subs.o: \
	qc_flags_subs.f90 \
	m_common_global_var.o
user_fluctuations.o: \
	user_fluctuations.f90 \
	m_common_global_var.o
writeout_qc_details.o: \
	writeout_qc_details.f90 \
	m_rp_global_var.o
random_error_handle.o: \
	random_error_handle.f90 \
	m_rp_global_var.o
bpcf_cospectral_models.o: \
	bpcf_cospectral_models.f90 \
	m_common_global_var.o
read_licor_ghg_archive.o: \
	read_licor_ghg_archive.f90 \
	m_rp_global_var.o
developed_turbulence_test.o: \
	developed_turbulence_test.f90 \
	m_rp_global_var.o
import_ascii.o: \
	import_ascii.f90 \
	m_common_global_var.o
wind_direction.o: \
	wind_direction.f90 \
	m_common_global_var.o
read_ini_rp.o: \
	read_ini_rp.f90 \
	m_rp_global_var.o
bpcf_li7550_analog_filters.o: \
	bpcf_li7550_analog_filters.f90 \
	m_common_global_var.o
rename_tmp_files_common.o: \
	rename_tmp_files_common.f90 \
	m_common_global_var.o
drift_retrieve_calibration_events.o: \
	drift_retrieve_calibration_events.f90 \
	m_rp_global_var.o
import_binary.o: \
	import_binary.f90 \
	m_common_global_var.o
write_processing_project_variables.o: \
	write_processing_project_variables.f90 \
	m_common_global_var.o
bpcf_moncrieff_97.o: \
	bpcf_moncrieff_97.f90 \
	m_common_global_var.o
eliminate_corrupted_variables.o: \
	eliminate_corrupted_variables.f90 \
	m_rp_global_var.o
import_native_data.o: \
	import_native_data.f90 \
	m_common_global_var.o
adjust_sonic_coordinates.o: \
	adjust_sonic_coordinates.f90 \
	m_rp_global_var.o
read_timelag_opt_file.o: \
	read_timelag_opt_file.f90 \
	m_rp_global_var.o
bpcf_read_full_cos_wt.o: \
	bpcf_read_full_cos_wt.f90 \
	m_common_global_var.o
rfftf.o: \
	rfftf.F
import_slt_edisol.o: \
	import_slt_edisol.f90 \
	m_common_global_var.o
air_and_cell_parameters.o: \
	air_and_cell_parameters.f90 \
	m_rp_global_var.o
retrieve_dynamic_metadata.o: \
	retrieve_dynamic_metadata.f90 \
	m_rp_global_var.o
change_24_into_00.o: \
	change_24_into_00.f90 \
	m_common_global_var.o
filter_data_for_diagnostics.o: \
	filter_data_for_diagnostics.f90 \
	m_rp_global_var.o
inform_of_metadata_problem.o: \
	inform_of_metadata_problem.f90 \
	m_common_global_var.o
aoa_calibration.o: \
	aoa_calibration.f90 \
	m_rp_global_var.o
set_licor_diagnostics.o: \
	set_licor_diagnostics.f90 \
	m_rp_global_var.o
copy_file.o: \
	copy_file.f90 \
	m_common_global_var.o
filter_raw_data_by_flags.o: \
	filter_raw_data_by_flags.f90 \
	m_rp_global_var.o
init_env.o: \
	init_env.f90 \
	m_common_global_var.o
binned_frequency_vector.o: \
	binned_frequency_vector.f90 \
	m_rp_global_var.o
sort_wind_by_sector.o: \
	sort_wind_by_sector.f90 \
	m_rp_global_var.o
create_master_timeseries.o: \
	create_master_timeseries.f90 \
	m_common_global_var.o
fix_planarfit_sectors.o: \
	fix_planarfit_sectors.f90 \
	m_rp_global_var.o
biomet_retrieve_embedded_data.o: \
	biomet_retrieve_embedded_data.f90 \
	m_rp_global_var.o
stationarity_test.o: \
	stationarity_test.f90 \
	m_rp_global_var.o
date_subs.o: \
	date_subs.f90 \
	m_common_global_var.o
flux_params.o: \
	flux_params.f90 \
	m_rp_global_var.o
biomet_subs.o: \
	biomet_subs.f90 \
	m_rp_global_var.o
storage.o: \
	storage.f90 \
	m_rp_global_var.o
def_sngl_indx_rng.o: \
	def_sngl_indx_rng.f90 \
	m_common_global_var.o
fluxes1_rp.o: \
	fluxes1_rp.f90 \
	m_rp_global_var.o
burba_terms.o: \
	burba_terms.f90 \
	m_rp_global_var.o
test_absolute_limits.o: \
	test_absolute_limits.f90 \
	m_rp_global_var.o
define_e2_set.o: \
	define_e2_set.f90 \
	m_common_global_var.o
gill_wm_w_boost.o: \
	gill_wm_w_boost.f90 \
	m_rp_global_var.o
configure_for_express.o: \
	configure_for_express.f90 \
	m_rp_global_var.o
test_attack_angle.o: \
	test_attack_angle.f90 \
	m_rp_global_var.o
define_used_variables.o: \
	define_used_variables.f90 \
	m_common_global_var.o
init_biomet_out.o: \
	init_biomet_out.f90 \
	m_rp_global_var.o
maths_subs.o: \
	maths_subs.f90 \
	m_common_global_var.o \
	m_numeric_kinds.o
create_datasets_rp.o: \
	create_datasets_rp.f90 \
	m_rp_global_var.o
retrieve_vars_selection.o: \
	retrieve_vars_selection.f90 \
	m_common_global_var.o
test_higher_moments.o: \
	test_higher_moments.f90 \
	m_rp_global_var.o
dir_sub.o: \
	dir_sub.f90 \
	m_common_global_var.o \
	m_typedef.o
init_external_biomet.o: \
	init_external_biomet.f90 \
	m_rp_global_var.o
max_wind_speed.o: \
	max_wind_speed.f90 \
	m_common_global_var.o
degraded_covariances.o: \
	degraded_covariances.f90 \
	m_rp_global_var.o
show_daily_advancement.o: \
	show_daily_advancement.f90 \
	m_common_global_var.o
test_spike_detection_mauder_13.o: \
	test_spike_detection_mauder_13.f90 \
	m_rp_global_var.o
exception_handler.o: \
	exception_handler.f90
init_user_outfiles.o: \
	init_user_outfiles.f90 \
	m_rp_global_var.o
metadata_file_validation.o: \
	metadata_file_validation.f90 \
	m_common_global_var.o
stats_operator_no_error.o: \
	stats_operator_no_error.f90 \
	m_common_global_var.o
test_timelag.o: \
	test_timelag.f90 \
	m_rp_global_var.o
files_in_chronological_order.o: \
	files_in_chronological_order.f90 \
	m_common_global_var.o
integral_turbulence_scale.o: \
	integral_turbulence_scale.f90 \
	m_rp_global_var.o
parse_file_name_with_prototype.o: \
	parse_file_name_with_prototype.f90 \
	m_common_global_var.o
CRA.o: \
	CRA.f90 \
	m_common_global_var.o
tag_run_mode.o: \
	tag_run_mode.f90 \
	m_common_global_var.o
tilt_correction.o: \
	tilt_correction.f90 \
	m_rp_global_var.o
fluctuations.o: \
	fluctuations.f90 \
	m_common_global_var.o
bpcf_Horst_97.o: \
	bpcf_Horst_97.f90 \
	m_common_global_var.o
planar_fit_aux_params.o: \
	planar_fit_aux_params.f90 \
	m_common_global_var.o
timestamp.o: \
	timestamp.f90 \
	m_numeric_kinds.o
write_out_ameriflux_rp.o: \
	write_out_ameriflux_rp.f90 \
	m_rp_global_var.o
fourier_transform.o: \
	fourier_transform.f90 \
	m_rp_global_var.o
optimize_timelags.o: \
	optimize_timelags.f90 \
	m_rp_global_var.o
bpcf_additional_horst_lenschow_09.o: \
	bpcf_additional_horst_lenschow_09.f90 \
	m_common_global_var.o
point_by_point_to_mixing_ratio.o: \
	point_by_point_to_mixing_ratio.f90 \
	m_common_global_var.o
ts_round_to_minute.o: \
	ts_round_to_minute.f90 \
	m_common_global_var.o
write_out_planar_fit.o: \
	write_out_planar_fit.f90 \
	m_rp_global_var.o
override_master_sonic_related_settings.o: \
	override_master_sonic_related_settings.f90 \
	m_rp_global_var.o
bpcf_anemometric_fluxes.o: \
	bpcf_anemometric_fluxes.f90 \
	m_common_global_var.o
potential_radiation.o: \
	potential_radiation.f90 \
	m_common_global_var.o
user_basic_stats.o: \
	user_basic_stats.f90 \
	m_common_global_var.o
write_output_files_rp.o: \
	write_output_files_rp.f90 \
	m_rp_global_var.o
power_of_two.o: \
	power_of_two.f90
bpcf_bandpass_spectral_corrections.o: \
	bpcf_bandpass_spectral_corrections.f90 \
	m_common_global_var.o
read_biomet_file.o: \
	read_biomet_file.f90 \
	m_rp_global_var.o
detect_master_sonic.o: \
	detect_master_sonic.f90 \
	m_rp_global_var.o
generate_t_cell.o: \
	generate_t_cell.f90 \
	m_common_global_var.o
user_timelag_handle.o: \
	user_timelag_handle.f90 \
	m_common_global_var.o
writeout_timelag_optimization.o: \
	writeout_timelag_optimization.f90 \
	m_rp_global_var.o
read_biomet_meta_file.o: \
	read_biomet_meta_file.f90 \
	m_rp_global_var.o
bpcf_fratini_12.o: \
	bpcf_fratini_12.f90 \
	m_common_global_var.o
read_metadata_file.o: \
	read_metadata_file.f90 \
	m_common_global_var.o
drift_correction.o: \
	drift_correction.f90 \
	m_rp_global_var.o
import_ascii_with_text.o: \
	import_ascii_with_text.f90 \
	m_common_global_var.o
write_out_stats.o: \
	write_out_stats.f90 \
	m_common_global_var.o
read_planar_fit_file.o: \
	read_planar_fit_file.f90 \
	m_rp_global_var.o
bpcf_massman_00.o: \
	bpcf_massman_00.f90 \
	m_common_global_var.o
replace_gap_with_linear_interpolation.o: \
	replace_gap_with_linear_interpolation.f90 \
	m_common_global_var.o
eddypro-rp_main.o: \
	eddypro-rp_main.f90 \
	m_rp_global_var.o
import_current_period.o: \
	import_current_period.f90 \
	m_rp_global_var.o
add_to_timelag_opt_dataset.o: \
	add_to_timelag_opt_dataset.f90 \
	m_rp_global_var.o
bpcf_only_lowfrequency_correction.o: \
	bpcf_only_lowfrequency_correction.f90 \
	m_common_global_var.o
rfftb.o: \
	rfftb.F
import_slt_eddysoft.o: \
	import_slt_eddysoft.f90 \
	m_common_global_var.o
adjust_timelag_opt_settings.o: \
	adjust_timelag_opt_settings.f90 \
	m_rp_global_var.o
rename_tmp_files_rp.o: \
	rename_tmp_files_rp.f90 \
	m_rp_global_var.o
calibrate_gas4.o: \
	calibrate_gas4.f90 \
	m_common_global_var.o
rffti.o: \
	rffti.F
import_tob1.o: \
	import_tob1.f90 \
	m_common_global_var.o
aoa_cal_nakai_2012.o: \
	aoa_cal_nakai_2012.f90 \
	m_rp_global_var.o
retrieve_sensor_params.o: \
	retrieve_sensor_params.f90 \
	m_rp_global_var.o
configure_for_embedded.o: \
	configure_for_embedded.f90 \
	m_common_global_var.o
filter_for_physical_thresholds.o: \
	filter_for_physical_thresholds.f90 \
	m_rp_global_var.o
init_continuous_dataset.o: \
	init_continuous_dataset.f90 \
	m_common_global_var.o
assess_daytime.o: \
	assess_daytime.f90 \
	m_rp_global_var.o
set_timelags.o: \
	set_timelags.f90 \
	m_rp_global_var.o
create_datasets_common.o: \
	create_datasets_common.f90 \
	m_common_global_var.o
fix_dataset_for_spectra.o: \
	fix_dataset_for_spectra.f90 \
	m_rp_global_var.o
biomet_enrich_vars_description.o: \
	biomet_enrich_vars_description.f90 \
	m_rp_global_var.o
spectral_analysis.o: \
	spectral_analysis.f90 \
	m_rp_global_var.o
cross_wind_corr.o: \
	cross_wind_corr.f90 \
	m_common_global_var.o
fix_timelag_opt_dataset.o: \
	fix_timelag_opt_dataset.f90 \
	m_rp_global_var.o
biomet_retrieve_external_data.o: \
	biomet_retrieve_external_data.f90 \
	m_rp_global_var.o
statistical_screening.o: \
	statistical_screening.f90 \
	m_rp_global_var.o
def_cross_indx_rng.o: \
	def_cross_indx_rng.f90 \
	m_common_global_var.o
fluxes0_rp.o: \
	fluxes0_rp.f90 \
	m_rp_global_var.o
biomet_units_conversions.o: \
	biomet_units_conversions.f90 \
	m_rp_global_var.o
tapering.o: \
	tapering.f90 \
	m_rp_global_var.o
define_all_var_set.o: \
	define_all_var_set.f90 \
	m_common_global_var.o
fluxes23_rp.o: \
	fluxes23_rp.f90 \
	m_rp_global_var.o
cleanup_e2_set.o: \
	cleanup_e2_set.f90 \
	m_rp_global_var.o
test_amp_res_drop_out.o: \
	test_amp_res_drop_out.f90 \
	m_rp_global_var.o
define_relative_separations.o: \
	define_relative_separations.f90 \
	m_common_global_var.o
infer_aoa_method.o: \
	infer_aoa_method.f90 \
	m_rp_global_var.o
make_dataset.o: \
	make_dataset.f90 \
	m_common_global_var.o
configure_for_md_retrieval.o: \
	configure_for_md_retrieval.f90 \
	m_rp_global_var.o
replace_sonic_temperature.o: \
	replace_sonic_temperature.f90 \
	m_common_global_var.o
test_discontinuities.o: \
	test_discontinuities.f90 \
	m_rp_global_var.o
define_user_set.o: \
	define_user_set.f90 \
	m_common_global_var.o
init_dynamic_medata.o: \
	init_dynamic_medata.f90 \
	m_rp_global_var.o
matrix_inversion.o: \
	matrix_inversion.f90 \
	m_common_global_var.o
default_vars_selection.o: \
	default_vars_selection.f90 \
	m_rp_global_var.o
set_os_environment.o: \
	set_os_environment.f90 \
	m_common_global_var.o
test_non_steady_wind.o: \
	test_non_steady_wind.f90 \
	m_rp_global_var.o
edit_ini_file.o: \
	edit_ini_file.f90 \
	m_common_global_var.o
init_outfiles_rp.o: \
	init_outfiles_rp.f90 \
	m_rp_global_var.o
median.o: \
	median.f90 \
	m_common_global_var.o
sort.o: \
	sort.f90 \
	m_common_global_var.o
despike_user_set.o: \
	despike_user_set.f90 \
	m_rp_global_var.o
spectral_density.o: \
	spectral_density.f90 \
	m_common_global_var.o
test_spike_detection_vickers_97.o: \
	test_spike_detection_vickers_97.f90 \
	m_rp_global_var.o
filelist_by_ext.o: \
	filelist_by_ext.f90 \
	m_common_global_var.o
int_2_flag.o: \
	int_2_flag.f90 \
	m_rp_global_var.o
multipliers_7700.o: \
	multipliers_7700.f90 \
	m_common_global_var.o
string_sub.o: \
	string_sub.f90 \
	m_common_global_var.o \
	m_numeric_kinds.o
add_to_cospectra_fit_dataset.o: \
	add_to_cospectra_fit_dataset.f90 \
	m_fx_global_var_mod.o
available_mean_spectra_cospectra.o: \
	available_mean_spectra_cospectra.f90 \
	m_fx_global_var_mod.o
correction_factor_model.o: \
	correction_factor_model.f90 \
	m_fx_global_var_mod.o \
	m_levenberg_marquardt.o
cospectra_qaqc.o: \
	cospectra_qaqc.f90 \
	m_fx_global_var_mod.o
bpcf_read_full_cos_wt.o: \
	bpcf_read_full_cos_wt.f90 \
	m_common_global_var.o
cospectra_sorting_and_averaging.o: \
	cospectra_sorting_and_averaging.f90 \
	m_fx_global_var_mod.o
eddypro-fcc_main.o: \
	eddypro-fcc_main.f90 \
	m_fx_global_var_mod.o
ensemble_cospectra_by_stability.o: \
	ensemble_cospectra_by_stability.f90 \
	m_fx_global_var_mod.o
fit_cospectral_models.o: \
	fit_cospectral_models.f90 \
	m_fx_global_var_mod.o \
	m_levenberg_marquardt.o
read_ex_record.o: \
	read_ex_record.f90 \
	m_common_global_var.o
read_ex2_record.o: \
	read_ex2_record.f90 \
	m_common_global_var.o
fit_rh_to_cutoff.o: \
	fit_rh_to_cutoff.f90 \
	m_fx_global_var_mod.o \
	m_levenberg_marquardt.o
fit_tf_models.o: \
	fit_tf_models.f90 \
	m_fx_global_var_mod.o \
	m_levenberg_marquardt.o
ridders_diff.o: \
	ridders_diff.f90 \
	m_common_global_var.o
fluxes1.o: \
	fluxes1.f90 \
	m_fx_global_var_mod.o
fluxes23.o: \
	fluxes23.f90 \
	m_fx_global_var_mod.o
extract_column_from_essentials.o: \
	extract_column_from_essentials.f90 \
	m_common_global_var.o
init_ex_vars.o: \
	init_ex_vars.f90 \
	m_fx_global_var_mod.o
fcn.o: \
	fcn.f90 \
	m_common_global_var.o
init_out_files.o: \
	init_out_files.f90 \
	m_fx_global_var_mod.o
m_fx_global_var_mod.o: \
	m_fx_global_var_mod.f90 \
	m_common_global_var.o
normalize_mean_spectra_cospectra.o: \
	normalize_mean_spectra_cospectra.f90 \
	m_fx_global_var_mod.o
output_spectral_assessment_results.o: \
	output_spectral_assessment_results.f90 \
	m_fx_global_var_mod.o
read_binned_file.o: \
	read_binned_file.f90 \
	m_fx_global_var_mod.o
ts_extract_subperiod_indexes_from_filelist.o: \
	ts_extract_subperiod_indexes_from_filelist.f90 \
	m_common_global_var.o
read_fx_ini.o: \
	read_fx_ini.f90 \
	m_fx_global_var_mod.o
kaimal_models.o: \
	kaimal_models.f90 \
	m_common_global_var.o
read_spectral_assessment_file.o: \
	read_spectral_assessment_file.f90 \
	m_fx_global_var_mod.o
retrieve_ex_var_by_timestamp.o: \
	retrieve_ex_var_by_timestamp.f90 \
	m_fx_global_var_mod.o
spectra_sorting_and_averaging.o: \
	spectra_sorting_and_averaging.f90 \
	m_fx_global_var_mod.o
subtract_high_freq_noise.o: \
	subtract_high_freq_noise.f90 \
	m_fx_global_var_mod.o
write_ameriflux_output.o: \
	write_ameriflux_output.f90 \
	m_fx_global_var_mod.o
report_imported_spectra.o: \
	report_imported_spectra.f90 \
	m_fx_global_var_mod.o
m_levenberg_marquardt.o: \
	m_levenberg_marquardt.f90
write_output_files_fcc.o: \
	write_output_files_fcc.f90 \
	m_fx_global_var_mod.o
