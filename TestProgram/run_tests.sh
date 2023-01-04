#!/bin/bash -e
while getopts 'h' c
do
  case $c in
    h) ../Make/Utils/write_mkflags.pl -h
    ;;
  esac
done

# if we are running inside a github action, change workdir
[ ! -z "$GITHUB_WORKSPACE" ] && cd $GITHUB_WORKSPACE/TestProgram

[ ! -d "$1" ] && echo First argument must be a subdirectory of TestProgram && exit 1

../Make/Utils/write_mkflags.pl -f ../Make/MkFlags.ini "${@: 2}" || exit 1

echo Building...
../Make/nj ${1}

rm -f ${1}/.test_failed ${1}/.test_failed_*

echo Run Tests...
../Make/nj ${1}_tests

./mk_summary.sh ${1}

if compgen -G "${1}/.test_failed_*" >/dev/null ; then 
  touch ${1}/.test_failed
  exit 1 ;
else 
  exit 0 ;
fi
