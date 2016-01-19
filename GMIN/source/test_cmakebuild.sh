#!/bin/bash

build_dir=$(readlink -f ../build_test)
source_dir=/scratch/khs26/svn/GMIN/source

mkdir ${build_dir} || exit 1
cd ${build_dir}

for compiler in gfortran pgf90 ifort; do
  # Check whether our compiler is loaded and try to load it, if not.
  ${compiler} --version > /dev/null 2>&1 || {
    module load ${compiler}
  }
  # If it still doesn't work, exit.
  ${compiler} --version > /dev/null 2>&1 || exit 1

  # Iterate through release and debug modes.
#  for build_type in Release Debug
  for build_type in Release
  do
    curdir=${build_dir}/${compiler}_${build_type}
    mkdir ${curdir}
    cd ${curdir}

    if [ ${compiler} == ifort ] || [ ${compiler} == pgf90 ]; then
      # If building with ifort or pgf90, we can build CHARMM35.
#      FC=${compiler} cmake -DCMAKE_BUILD_TYPE=${build_type} -DWITH_DMAGMIN=no -DWITH_CHARMM35=yes -DWITH_AMBER9=yes -DWITH_AMBER12=yes -DWITH_AMH=no ${source_dir} > /dev/null 2>&1
      FC=${compiler} cmake -DCMAKE_BUILD_TYPE=${build_type} -DWITH_AMBER9=yes -DWITH_AMBER12=yes -DWITH_AMH=yes -DWITH_OPEP=yes -DWITH_SPIN=yes -DWITH_OXDNA=yes ${source_dir} > /dev/null 2>&1
      for target in GMIN A9GMIN A12GMIN OPEPGMIN AMHGMIN SPINGMIN OXDNAGMIN; do
        echo -n "${compiler} ${build_type} ${target}: "
        make -j 8 ${target}  > /dev/null 2>&1
        if [ "$?" == "0" ]; then
          echo -e "\e[01;32mok\e[00m"
        else 
          echo -e "\e[00;31mfailed\e[00m"
        fi
      done
    else
      # Otherwise, don't build CHARMM35.
#      FC=${compiler} cmake -DCMAKE_BUILD_TYPE=${build_type} -DWITH_DMAGMIN=no -DWITH_CHARMM35=no -DWITH_AMBER9=yes -DWITH_AMBER12=yes -DWITH_AMH=no ${source_dir} > /dev/null 2>&1
      FC=${compiler} cmake -DCMAKE_BUILD_TYPE=${build_type} -DWITH_AMBER9=yes -DWITH_AMBER12=yes -DWITH_AMH=yes -DWITH_OPEP=yes -DWITH_SPIN=yes -DWITH_OXDNA=yes ${source_dir} > /dev/null 2>&1
#      FC=${compiler} cmake -DCMAKE_BUILD_TYPE=${build_type} -DWITH_DMAGMIN=no -DWITH_CHARMM35=no -DWITH_AMBER9=yes -DWITH_AMBER12=yes -DWITH_AMH=no ${source_dir}
      for target in GMIN A9GMIN A12GMIN OPEPGMIN AMHGMIN SPINGMIN OXDNAGMIN; do
#      for target in GMIN A9GMIN A12GMIN; do
        echo -n "${compiler} ${build_type} ${target}: "
        make -j 8 ${target}  > /dev/null 2>&1
        if [ "$?" == "0" ]; then
          echo -e "\e[01;32mok\e[00m"
        else 
          echo -e "\e[00;31mfailed\e[00m"
        fi
      done
    fi
  done
done
