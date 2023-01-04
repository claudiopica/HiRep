#!/bin/bash
if [ $# -eq 0 ] ; then
    SUBDIRS=`find . -maxdepth 1 -mindepth 1 -type d `;
else
    SUBDIRS="$@";
fi

SEP="=============================\n"
printf "$SEP|       Test Summary        |\n$SEP";
rm -f .test_failed
for i in $SUBDIRS ; do
    failed=`compgen -G "$i/.test_failed_*" | awk -v RS=  '{$1=$1}1'`;
	if [ -z "$failed" ]; then
        printf "%19s  \e[32m%s\e[0m\n" `basename $i` "[ OK ]";
	else
        printf "\e[1;41m%19s\e[0m  \e[1;31m%s\e[0m\n" `basename $i` "[FAIL]";
        for f in $failed ; do
            printf "\e[1;31m%27s\e[0m\n" `basename ${f/.test_failed_/}`;
        done 
	    touch .test_failed ; 
	fi ;
done
if [ -f .test_failed ];
then exit 1 ; 
else exit 0 ; 
fi
