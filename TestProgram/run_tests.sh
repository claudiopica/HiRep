#!/bin/bash -e

# if we are running inside a github action, change workdir
if [ ! -z "$GITHUB_WORKSPACE" ] ; then 
  cd $GITHUB_WORKSPACE/TestProgram
  export OMPI_ALLOW_RUN_AS_ROOT=1
  export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
fi

echo Working directory: `pwd`

while getopts ':h' c
do
  case $c in
    h) ../Make/Utils/write_mkflags.pl -h
    ;;
  esac
done

# if we run inside a docker container, remove openmpi weirdness
[ -f /.dockerenv ] && export OMPI_MCA_btl_vader_single_copy_mechanism=none

dirList=$1;
echo "Testing: $dirList"
echo "With flags:" "${@: 2}"

# write flags file
../Make/Utils/write_mkflags.pl -f ../Make/MkFlags.ini "${@: 2}" || exit 1

echo Building...
../Make/nj $dirList

for dir in $dirList ; do

[ ! -d "$dir" ] && echo Argument must be a subdirectory of TestProgram && exit 1

rm -f $dir/.test_failed $dir/.test_failed_*

echo Run Tests...
../Make/nj ${dir}_tests

done

./mk_summary.sh "${dirList[*]}"
exit $?
