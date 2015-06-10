#!/bin/bash
help(){
        echo -e "Usage -c<channel> -T<Lt> -L<Ls> -i<input meson> [-e <min efficiency>] [-b <blocksize>] -m<method> [-f]"
	echo -e "\tc) Channel name"
	echo -e "\tLt) Length in the T dir"
	echo -e "\tLs) Length in the Spatial dir"
	echo -e "\ti) mesonic input file"
	echo -e "\te) minimal requested efficiency (default is 100%)"
	echo -e "\tb) block size"
	echo -e "\tm) method of analysis:\n\t   0 -> effective mass\n\t   1 -> Prony no excitation\n\t   2 -> Prony one excitation"
	echo -e "\tf) force recomputation"
        exit 0;
}



EXECDIR=`dirname $0`
EXEC="${EXECDIR}/bs_mesons"
OUTDIR="output"
MEFFDIR="eff_mass"
EFFDIR="efficiency"
RCUTDIR="rightcut"
MKEFFTABLE="${EXECDIR}/mk_efficiencytable.sh"
FINALCUTDIR="cut"
FORCEFLAG=""

while getopts "c:T:L:i:e:m:b:f" opt; do
    case $opt in
	h ) help ;;
	c ) CHANNEL=$OPTARG;;
	T ) LT=$OPTARG;;
	L ) LS=$OPTARG;;
	i ) INPUT=$OPTARG;;
	e ) MIN_EFFIC=$OPTARG;;
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

if [ ! -f ${MKEFFTABLE} ] 
    then
    echo "$0: Missing the script mk_efficiencytable.sh" 
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

[ ! -d ${EFFDIR} ] && mkdir ${EFFDIR}
[ ! -d ${OUTDIR} ] && mkdir ${OUTDIR}
[ ! -d ${FINALCUTDIR} ] && mkdir ${FINALCUTDIR}

if [ ! -z "${MIN_EFFIC}" ]
    then
    CMD_EFF="-e ${MIN_EFFIC}"
fi

if [ -z "${BLKSIZE}" ]
    then
    BLKSIZE=30
fi


CHANTYPE=`${EXEC} info ${CHANNEL}|  grep -e"^2 EVAL_CTRL ${CHANNEL} " | tail -n1 | wc -w`
BASENAME=`basename $INPUT`
NAME="${CHANNEL}_M${METHOD}_`basename $INPUT`"
CUTFILE="${OUTDIR}/${NAME}"
OUTFILE="${OUTDIR}/${NAME}_TMP"
FINALCUTNAME="${FINALCUTDIR}/${BASENAME}_M${METHOD}"

if [ ! -f ${FINALCUTNAME} ]
then
  > ${FINALCUTNAME}
fi

> ${CUTFILE}

echo -n "Evaluating $CHANNEL effective mass:"

RCUT=`awk '$1=="'$CHANNEL'" {print $3}' ${FINALCUTNAME}`
if [ ! -z "$RCUT" ] && [ "$FORCEFLAG" != "-f" ] ; then
    echo " Channel already evaluated."
    exit 0;
fi
echo " Evaluating channel."

if [ "$FORCEFLAG" == "-f" ] ; then
    
    WANTED_PRIM_CHAN=`${EXEC} info ${CHANNEL}|  grep -e"^2 PRIMARY_CTRL .* P_EFF_B1 .* P_EFF_B2" | awk '$9==1 || $13==1 {print $3}' | tr '\n' ' '`
    WANTED_DER_CHAN=`${EXEC} info ${CHANNEL}| grep -e"^2 EVAL_CTRL .* D_EFF_B1 .* D_EFF_B2" | awk '$5==1 || $7==1 {print $3}' | tr '\n' ' '`
    WANTED_CHAN="$WANTED_PRIM_CHAN $WANTED_DER_CHAN"

    for i in $WANTED_CHAN
      do
      awk '$1!="'$i'" {print}' ${FINALCUTNAME} > ${FINALCUTNAME}.tmp
      mv ${FINALCUTNAME}.tmp ${FINALCUTNAME}
    done
fi


if [ "$CHANTYPE" -eq "13"  ]
    then
    WANTED_PRIM_CHAN=`${EXEC} info ${CHANNEL}|  grep -e"^2 PRIMARY_CTRL .* P_EFF_B1 .* P_EFF_B2" | awk '$9==1 || $13==1 {print $3}' | tr '\n' ' '`
    WANTED_DER_CHAN=`${EXEC} info ${CHANNEL}| grep -e"^2 EVAL_CTRL .* D_EFF_B1 .* D_EFF_B2" | awk '$5==1 || $7==1 {print $3}' | tr '\n' ' '`

    for i in $WANTED_PRIM_CHAN 
      do
      $0 -c$i -T$LT -L$LS $CMD_EFF -i$INPUT -b$BLKSIZE -m$METHOD 
      RCUT=`awk '$1=="'$i'" {print $3}' $FINALCUTNAME`
      echo "$i 4 $RCUT" >> ${CUTFILE}
    done

    for i in $WANTED_DER_CHAN 
      do
      if [ "$i" != "${CHANNEL}" ] ; then
	  $0 -c$i -T$LT -L$LS $CMD_EFF -i$INPUT -b$BLKSIZE -m$METHOD 
	  RCUT=`awk '$1=="'$i'" {print $3}' $FINALCUTNAME`
	  echo "$i 4 $RCUT" >> ${CUTFILE}
      fi
    done

    echo "$CHANNEL 4 $((LT/2))" >> ${CUTFILE}
    $EXEC ${CHANNEL}_eff $METHOD $INPUT ${CUTFILE} $LT $LS $BLKSIZE 1000 100 > $OUTFILE

    awk '$1!="'$i'" {print}' ${FINALCUTNAME} > ${FINALCUTNAME}.tmp
    mv ${FINALCUTNAME}.tmp ${FINALCUTNAME}
    echo -e "$i\t*\t$((LT/2))" >> ${FINALCUTNAME}      
    grep $OUTFILE -e "10 EFF " > ${MEFFDIR}/${NAME}
    rm ${CUTFILE} $OUTFILE 
    
fi


if [ "$CHANTYPE" -eq "5"  ]
    then
    WANTED_PRIM_CHAN=`${EXEC} info ${CHANNEL}|  grep -e"^2 PRIMARY_CTRL .* P_EFF_B1 .* P_EFF_B2" | awk '$9==1 || $13==1 {print $3}' | tr '\n' ' '`
    WANTED_DER_CHAN=`${EXEC} info ${CHANNEL}| grep -e"^2 EVAL_CTRL .* D_EFF_B1 .* D_EFF_B2" | awk '$5==1 || $7==1 {print $3}' | tr '\n' ' '`

    for i in $WANTED_PRIM_CHAN 
      do
      $0 -c$i -T$LT -L$LS $CMD_EFF -i$INPUT -b$BLKSIZE -m$METHOD 
    done

    for i in $WANTED_DER_CHAN 
      do
      $0 -c$i -T$LT -L$LS $CMD_EFF -i$INPUT -b$BLKSIZE -m$METHOD 
    done
fi

if [ "$CHANTYPE" -eq "21"  ]
    then
    ${MKEFFTABLE} -c$CHANNEL -T$LT -L$LS $CMD_EFF -b$BLKSIZE -i$INPUT -m$METHOD
    [ -z "`cat $RCUTDIR/$BASENAME | grep -e$CHANNEL | tail -n1`" ] && echo "Unable to eval right cut for $CHANNEL" && exit 1

    RCUT=`cat $RCUTDIR/$BASENAME | grep -e"^$CHANNEL " | tail -n1| awk '{print $5}'`

    awk '$1!="'$CHANNEL'" {print}' ${FINALCUTNAME} > ${FINALCUTNAME}.tmp
    mv ${FINALCUTNAME}.tmp ${FINALCUTNAME}
    echo -e "$CHANNEL\t*\t$RCUT" >> ${FINALCUTNAME}
fi
