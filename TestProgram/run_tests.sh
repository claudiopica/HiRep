#!/bin/bash -e
while getopts 'h' c
do
  case $c in
    h) echo heko ;
      ../Make/Utils/write_mkflags.pl -h
    ;;
  esac
done

[ ! -d "$1" ] && echo First argument must be a subdirectory of TestProgram && exit 1

# if we are running inside a github action, change workdir
[ ! -z "$GITHUB_WORKSPACE" ] && cd $GITHUB_WORKSPACE/TestProgram

../Make/Utils/write_mkflags.pl -f ../Make/MkFlags ${@: 2} || exit 1

echo Cleaning...
( cd .. && make cleanall )

cd ./${1}

echo Building...
make -j1 

echo Run Tests...
make runtests
