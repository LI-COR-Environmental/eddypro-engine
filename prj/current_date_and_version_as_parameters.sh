#!/bin/bash
timestamp=`date "+%Y-%m-%d, %H:%M"`
version='7.0.0'
echo "character(5), parameter  :: sw_ver = '"$version"'" > ../src/src_common/version_and_date.inc
echo "character(17), parameter :: build_date = '"$timestamp"'" >> ../src/src_common/version_and_date.inc
touch ../src/src_common//init_env.f90