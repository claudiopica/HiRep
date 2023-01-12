#!/bin/bash

FILENAME=cinfo.c
MKDIR=$1
TOPDIR=$2
PWD=`pwd`

shift; shift
MACROS=$@

> ${FILENAME}

CMDOUT=`echo ${MACROS} | perl -pe 's/"/\\\"/g'`
echo "char MACROS[] = \""$CMDOUT"\\n\";" >> ${FILENAME}

CMDOUT=`perl -pe 's/\n/\\\n/g' ${MKDIR}/MkFlags.ini`
echo "char CI_mkflags[] = \""$CMDOUT"\\n\";" >> ${FILENAME}

if command -v lscpu >/dev/null 2>&1
then
    CMDOUT=`lscpu | perl -pe 's/\n/\\\n/g'`
else
if [[ -a /proc/cpuinfo ]] 
then
 CMDOUT=`perl -pe 's/\n/\\\n/g' /proc/cpuinfo`
else
 CMDOUT="No CPU info"
fi
fi
echo "char CI_cpuinfo[] = \""$CMDOUT"\\n\";" >> ${FILENAME}

if uname >/dev/null 2>/dev/null
then
 CMDOUT=`uname -a`
else
 CMDOUT="No VERSION info\n"
fi
echo "char CI_linux[] = \""$CMDOUT"\\n\";" >> ${FILENAME}

CMDOUT=`gcc -v 2>&1 | perl -pe 's/\n/\\\n/g'`
echo "char CI_gcc[] = \""$CMDOUT"\\n\";" >> ${FILENAME}

if command -v git >/dev/null 2>/dev/null
then
    echo "char CI_gitinfo[] = \""`git rev-parse --symbolic-full-name`"\";" >>${FILENAME}
    echo "char CI_gitrevision[] = \""`git log --format="%H" -n 1`"\";" >>${FILENAME}
else
    echo "char CI_gitinfo[] = \"\";" >> ${FILENAME}
    echo "char CI_gitrevision[11] = \""No version"\";" >>${FILENAME}
fi
