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

#svn --version >&- ; ret=$?
#if [ "${ret}" -eq "0" ] ; then
#  svn info ${TOPDIR} | awk '{printf "%s\\n",$0}' > cinfo.tmp
#  len=`cat cinfo.tmp | wc -c`+1
#  echo -n "static char CI_svninfo[${len}] = \"" >> ${FILENAME}
#  cat cinfo.tmp >> ${FILENAME}
#  echo "\";" >> ${FILENAME}
#  echo "" >> ${FILENAME}
#  rm cinfo.tmp

#  svn st -q ${TOPDIR} | awk '{printf "%s\\n",$0}' > cinfo.tmp
#  len=`cat cinfo.tmp | wc -c`+1
#  echo -n "static char CI_svnstatus[${len}] = \"" >> ${FILENAME}
#  cat cinfo.tmp >> ${FILENAME}
#  echo "\";" >> ${FILENAME}
#  echo "" >> ${FILENAME}
#  rm cinfo.tmp
#else
  echo -n "static char CI_svninfo[1] = \"\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
  echo -n "static char CI_svnstatus[1] = \"\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
#fi

#REV=$(svn info | grep Revision | awk '{ print $2 }')
#if [ -z "$REV" ]; then
  REV="0"
#fi
echo "static int CI_svnrevision = ${REV};" >> ${FILENAME}
echo "" >> ${FILENAME}

cat ${MKDIR}/Utils/${FILENAME}.tmpl >> ${FILENAME}
