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
echo -n "static char mkflags[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

awk '{printf "%s\\n",$0}' /proc/cpuinfo > cinfo.tmp
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char cpuinfo[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

awk '{printf "%s\\n",$0}' /proc/version > cinfo.tmp
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char linux[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

gcc -v 2>&1 | awk '{printf "%s\\n",$0}' > cinfo.tmp
len=`cat cinfo.tmp | wc -c`+1
echo -n "static char gcc[${len}] = \"" >> ${FILENAME}
cat cinfo.tmp >> ${FILENAME}
echo "\";" >> ${FILENAME}
echo "" >> ${FILENAME}
rm cinfo.tmp

svn --version >&- ; ret=$?
if [ "${ret}" -eq "0" ] ; then
  svn info ${TOPDIR} | awk '{printf "%s\\n",$0}' > cinfo.tmp
  len=`cat cinfo.tmp | wc -c`+1
  echo -n "static char svninfo[${len}] = \"" >> ${FILENAME}
  cat cinfo.tmp >> ${FILENAME}
  echo "\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
  rm cinfo.tmp

  svn st -q ${TOPDIR} | awk '{printf "%s\\n",$0}' > cinfo.tmp
  len=`cat cinfo.tmp | wc -c`+1
  echo -n "static char svnstatus[${len}] = \"" >> ${FILENAME}
  cat cinfo.tmp >> ${FILENAME}
  echo "\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
  rm cinfo.tmp
else
  echo -n "static char svninfo[1] = \"\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
  echo -n "static char svnstatus[1] = \"\";" >> ${FILENAME}
  echo "" >> ${FILENAME}
fi

cat ${MKDIR}/Utils/${FILENAME}.tmpl >> ${FILENAME}
