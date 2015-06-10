#!/bin/bash
help(){
        echo -e "Usage -c<channel> -T<Lt> -L<Ls> -i<input meson> -b<blocksize> -m<method> [-f]"
	echo -e "\tc) Channel name"
	echo -e "\tLt) Length in the T dir"
	echo -e "\tLs) Length in the Spatial dir"
	echo -e "\ti) mesonic input file"
	echo -e "\tb) block size"
	echo -e "\tm) method of analysis:\n\t   0 -> effective mass\n\t   1 -> Prony no excitation\n\t   2 -> Prony one excitation"
	echo -e "\tf) force recomputation"
        exit 0;
}



EXECDIR=`dirname $0`
EXEC="${EXECDIR}/bs_mesons"
OUTDIR="output"
RESDIR="results"
FINALCUTDIR="cut"
FORCEFLAG=""

while getopts "c:T:L:i:m:b:f" opt; do
    case $opt in
	h ) help ;;
	c ) CHANNEL=$OPTARG;;
	T ) LT=$OPTARG;;
	L ) LS=$OPTARG;;
	i ) INPUT=$OPTARG;;
	b ) BLKSIZE=$OPTARG;;
	m ) METHOD=$OPTARG;;
	f ) FORCEFLAG="-f";;
	* ) help ;;
    esac
done
shift $[ OPTIND - 1 ]

if [ -z "$CHANNEL" ] ||  [ -z "$LT" ] || [ -z "$LS" ] || [ -z "$INPUT" ]  || [ -z "$METHOD" ]
    then
    echo "$0: Missing parameter" 
    help
fi

if [ ! -f ${EXEC} ] 
    then
    echo "$0: Missing the executable file (${EXEC})" 
    exit 0
fi

if [ ! -f ${INPUT} ] 
    then
    echo "$0: Missing the input file (${INPUT})" 
    exit 0
fi

if [ ! "${METHOD}" -eq "0"  ] && [ ! "${METHOD}" -eq "1"  ] && [ ! "${METHOD}" -eq "2"  ]  
    then
    echo "$0: Method can take only values 0,1,2" 
    exit 0
fi


if [  "${CHANNEL}" != "mvmps" ] && [  "${CHANNEL}" != "mvkmps" ] && [  "${CHANNEL}" != "mpsfps" ] && [  "${CHANNEL}" != "mvfps" ] && [  "${CHANNEL}" != "mvkfps" ] && [  "${CHANNEL}" != "gmor" ] && [  "${CHANNEL}" != "mpsmpcac" ] && [  "${CHANNEL}" != "mps2mpcac" ] && [  "${CHANNEL}" != "fvfps" ] && [  "${CHANNEL}" != "fvkfps" ] && [  "${CHANNEL}" != "mamps" ] && [  "${CHANNEL}" != "makmps" ] && [  "${CHANNEL}" != "mamv" ] && [  "${CHANNEL}" != "makmvk" ] && [  "${CHANNEL}" != "id" ] && [  "${CHANNEL}" != "g0" ] && [  "${CHANNEL}" != "g5" ] && [  "${CHANNEL}" != "g0g5" ] && [  "${CHANNEL}" != "g1" ] && [  "${CHANNEL}" != "g2" ] && [  "${CHANNEL}" != "g3" ] && [  "${CHANNEL}" != "gk" ] && [  "${CHANNEL}" != "g0g1" ] && [  "${CHANNEL}" != "g0g2" ] && [  "${CHANNEL}" != "g0g3" ] && [  "${CHANNEL}" != "g0gk" ] && [  "${CHANNEL}" != "g5g1" ] && [  "${CHANNEL}" != "g5g2" ] && [  "${CHANNEL}" != "g5g3" ] && [  "${CHANNEL}" != "g5gk" ] && [  "${CHANNEL}" != "g0g5g1" ] && [  "${CHANNEL}" != "g0g5g2" ] && [  "${CHANNEL}" != "g0g5g3" ] && [  "${CHANNEL}" != "g0g5gk" ] && [ "${CHANNEL}" != "mpcac" ] && [ "${CHANNEL}" != "gps" ] && [ "${CHANNEL}" != "fps" ] && [ "${CHANNEL}" != "fv" ] && [ "${CHANNEL}" != "fvk" ]  && [ "${CHANNEL}" != "fa" ] && [ "${CHANNEL}" != "fak" ] && [ "${CHANNEL}" != "fid" ]
    then
    echo "$0: Wrong channel declaration, only primary or derived channels are allowed (${CHANNEL})"
    exit 0
fi

[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}
[ ! -d ${RESDIR} ] && mkdir ${RESDIR}

if [ -z "${BLKSIZE}" ]
    then
    BLKSIZE=30
fi


CHANTYPE=`${EXEC} info ${CHANNEL}|  grep -e"^2 EVAL_CTRL ${CHANNEL} " | tail -n1 | wc -w`
BASENAME=`basename $INPUT`
NAME="${BASENAME}_M${METHOD}"
OUTFILE="${OUTDIR}/${NAME}_TMP"
FINALCUTNAME="${FINALCUTDIR}/${BASENAME}_M${METHOD}"

[ -f ${RESDIR}/${NAME} ] && EVALUATED=`awk '$1=="'$CHANNEL'" {print $2}' ${RESDIR}/${NAME}`
if [ ! -z "$EVALUATED" ] && [ "$FORCEFLAG" != "-f" ] ; then
    echo " Channel $CHANNEL already evaluated use the flag -f to force re-evaluation."
    exit 0;
fi
echo " Evaluating channel ${CHANNEL}."


if [ ! -f ${FINALCUTNAME} ]
then
  echo "$0: ${FINALCUTNAME} missing."
  exit 0
fi

if [ ! -f ${RESDIR}/${NAME} ]
then
  > ${RESDIR}/${NAME}
fi

WANTED_PRIM_CHAN=`${EXEC} info ${CHANNEL}|  grep -e"^2 PRIMARY_CTRL .* P_EFF_B1 .* P_EFF_B2" | awk '$9==1 || $13==1 {print $3}' | tr '\n' ' '`
WANTED_DER_CHAN=`${EXEC} info ${CHANNEL}| grep -e"^2 EVAL_CTRL .* D_EFF_B1 .* D_EFF_B2" | awk '$5==1 || $7==1 {print $3}' | tr '\n' ' '`
WANTED_CHAN="$WANTED_PRIM_CHAN $WANTED_DER_CHAN"

for i in $WANTED_CHAN
  do
  LCUT=`awk '$1=="'$i'" {print $2}' ${FINALCUTNAME}`
  if [ -z "$LCUT" ] || [ "$LCUT" == "*" ] 
      then
      echo "$0: Right cut not defined for channel $i"
      exit 1
  fi
done



$EXEC ${CHANNEL} $METHOD $INPUT ${FINALCUTNAME} $LT $LS $BLKSIZE 1000 100 > $OUTFILE

FIT=`grep $OUTFILE -e "10 FIT "`

if [ ! -z "$FIT" ]
then
  awk '$1!="'$CHANNEL'" {print}' ${RESDIR}/${NAME} > ${RESDIR}/${NAME}.tmp
  mv ${RESDIR}/${NAME}.tmp ${RESDIR}/${NAME}
  grep $OUTFILE -e "10 FIT " | awk '{print $3, $4, $5, "#", $7, $9, $11}' >> ${RESDIR}/${NAME}
  rm $OUTFILE
else
  echo "$0: $EXEC unable to perform the fit for $CHANNEL"
  echo "$0: Output in $OUTFILE"
  exit 0
fi
