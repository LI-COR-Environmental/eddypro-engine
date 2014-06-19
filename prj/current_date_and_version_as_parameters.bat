set datetime=%date%, %time%
set datetime=%datetime:~0,17%
echo character(5), parameter  :: sw_ver = '5.2.0' > ..\src\src_common\version_and_date.inc
echo character(17), parameter :: build_date = '%datetime%' >> ..\src\src_common\version_and_date.inc
cd ..\src\src_common\
copy /b init_env.f90 +,, > nul
