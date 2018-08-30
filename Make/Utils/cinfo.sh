#!/bin/bash

FILENAME=cinfo.c
MKDIR=$1
TOPDIR=$2
PWD=`pwd`

shift; shift
MACROS=$@

> ${FILENAME}


echo -n ${MACROS} | tr '"' '@' > cinfo.tmp
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char MACROS[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

awk '{printf "%s\\n",$0}' ${MKDIR}/MkFlags > cinfo.tmp
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char CI_mkflags[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

if [[ -a /proc/cpuinfo ]] 
then
 awk '{printf "%s\\n",$0}' /proc/cpuinfo > cinfo.tmp
else
 echo -n "No CPU info\n" > cinfo.tmp
fi
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char CI_cpuinfo[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

if [[ -a /proc/version ]] 
then
 awk '{printf "%s\\n",$0}' /proc/version > cinfo.tmp
else
 echo -n "No VERSION info\n" > cinfo.tmp
fi
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char CI_linux[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

gcc -v 2>&1 | awk '{printf "%s\\n",$0}' > cinfo.tmp
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char CI_gcc[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

if [ -x "`command -v git`" ]
    then
    len=`git rev-parse --symbolic-full-name|wc -c`+1
    echo -n "static char CI_gitinfo[${len}] = \"" `git rev-parse --symbolic-full-name` "\";" >>${FILENAME}
    echo "" >> ${FILENAME}
    
    len=`git log --format="%H" -n 1|wc -c`+1
    echo "static char CI_gitrevision[${len}] = \"" `git log --format="%H" -n 1` "\";" >>${FILENAME}
    echo "" >> ${FILENAME}
else
    echo -n "static char CI_gitinfo[1] = \"\";" >> ${FILENAME}
    echo "" >> ${FILENAME}

    echo "static char CI_gitrevisio[11] = \""No version"\";" >>${FILENAME}
    echo "" >> ${FILENAME}
fi

cat ${MKDIR}/Utils/${FILENAME}.tmpl >> ${FILENAME}