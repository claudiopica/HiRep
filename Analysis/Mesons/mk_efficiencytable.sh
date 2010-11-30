#!/bin/bash


help(){
        echo -e "Usage: $0 -c<channel> -T<Lt> -L<Ls> -i<input meson> [-e <min efficiency>] -b<blocksize> -m<method>"
	echo -e "\tc) Channel name"
	echo -e "\tLt) Length in the T dir"
	echo -e "\tLs) Length in the Spatial dir"
	echo -e "\ti) mesonic input file"
	echo -e "\te) minimal requested efficiency (default is 100%)"
	echo -e "\tb) block size"
	echo -e "\tm) method of analysis:\n\t   0 -> effective mass\n\t   1 -> Prony no excitation\n\t   2 -> Prony one excitation"
        exit 0;
}



EXECDIR=`dirname $0`
EXEC="${EXECDIR}/bs_mesons"
OUTDIR="output"
MEFFDIR="eff_mass"
EFFDIR="efficiency"
RCUTDIR="rightcut"


while getopts "c:T:L:i:e:m:b:" opt; do
    case $opt in
	h ) help ;;
	c ) CHANNEL=$OPTARG;;
	T ) LT=$OPTARG;;
	L ) LS=$OPTARG;;
	i ) INPUT=$OPTARG;;
	e ) MIN_EFFIC=$OPTARG;;
	b ) BLKSIZE=$OPTARG;;
	m ) METHOD=$OPTARG;;
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


if [  "${CHANNEL}" != "id" ] && [  "${CHANNEL}" != "g0" ] && [  "${CHANNEL}" != "g5" ] && [  "${CHANNEL}" != "g0g5" ] && [  "${CHANNEL}" != "g1" ] && [  "${CHANNEL}" != "g2" ] && [  "${CHANNEL}" != "g3" ] && [  "${CHANNEL}" != "gk" ] && [  "${CHANNEL}" != "g0g1" ] && [  "${CHANNEL}" != "g0g2" ] && [  "${CHANNEL}" != "g0g3" ] && [  "${CHANNEL}" != "g0gk" ] && [  "${CHANNEL}" != "g5g1" ] && [  "${CHANNEL}" != "g5g2" ] && [  "${CHANNEL}" != "g5g3" ] && [  "${CHANNEL}" != "g5gk" ] && [  "${CHANNEL}" != "g0g5g1" ] && [  "${CHANNEL}" != "g0g5g2" ] && [  "${CHANNEL}" != "g0g5g3" ] && [  "${CHANNEL}" != "g0g5gk" ]
    then
    echo "$0: Wrong channel declaration, only primary channels are allowed"
    exit 0
fi

if [ -z "$BLKSIZE" ]
    then
    BLKSIZE=30
fi

if [ -z "$MIN_EFFIC" ]
    then
    MIN_EFFIC=100
fi

NAME="${CHANNEL}_M${METHOD}_`basename $INPUT`"

[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}
[ ! -d ${MEFFDIR} ] && mkdir ${MEFFDIR}
[ ! -d ${EFFDIR} ] && mkdir ${EFFDIR}
[ ! -d ${RCUTDIR} ] && mkdir ${RCUTDIR}

CONTROL=1
LEFTCUT=4
RIGHTCUT=$((LT/2))

OUTFILE=${OUTDIR}/${NAME}_TMP

while [ "$CONTROL" -eq "1" ] 
  do
  
  if ((RIGHTCUT<LEFTCUT))
      then
      echo "Unable to obtain required Efficiency"
      exit 0
  fi
  
  CUTFILE=${OUTDIR}/${NAME}_rcut$RIGHTCUT
  echo "${CHANNEL} ${LEFTCUT} ${RIGHTCUT}" > $CUTFILE
  $EXEC ${CHANNEL}_eff $METHOD $INPUT $CUTFILE $LT $LS $BLKSIZE 1000 100 > $OUTFILE

  EFFIC=`grep $OUTFILE -e "20 EFFICIENCY_BS1" | awk '{print $3}'`
  
  [ ! -f "${EFFDIR}/`basename $INPUT`" ] && echo "# Channel Efficiency Method LeftCut RightCut" > ${EFFDIR}/`basename $INPUT`
  
  echo $CHANNEL $EFFIC $METHOD $LEFTCUT $RIGHTCUT >> ${EFFDIR}/`basename $INPUT`
  grep $OUTFILE -e "10 EFF " > ${MEFFDIR}/${NAME}_rcut${RIGHTCUT}
  echo -e "Efficiency: ${EFFIC}/${MIN_EFFIC}"
  if ((EFFIC>=MIN_EFFIC)) 
      then
      CONTROL=0
      grep $OUTFILE -e "10 EFF " > ${MEFFDIR}/${NAME} 
      echo $CHANNEL $EFFIC $METHOD $LEFTCUT $RIGHTCUT >> ${RCUTDIR}/"`basename $INPUT`"
  fi
  
  ((RIGHTCUT--))

  rm -f $OUTFILE
done