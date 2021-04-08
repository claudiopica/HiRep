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

../Make/Utils/write_mkflags.pl -f ../Make/MkFlags ${@: 2} || exit 1

echo Cleaning...
( cd .. && make cleanall )

cd ./${1}

echo Building...
make -j1 

echo Run Tests...
make runtests

if [ -f .test_failed ]
then 
  exit 1 ;
else 
  exit 0 ;
fi
