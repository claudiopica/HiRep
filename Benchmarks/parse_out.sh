#!/bin/bash

print_report_header () {
echo "#Number of Threads
#Gauge group
#MPI_size
#Global size (4x)
#Local size  (4x)
#linear algebra kernel (MACRO,AVX2,VECT)
#geometry type (NEW,OLD)
#RLXD
#Flop per site
#Byte per site
#Dirac data movement
#Massless Diracoperator (reps, data size in kb, time in msec, GFLOPS, BAND in GB/s)
#Massless fused Diracoperator (reps, data size in kb, time in msec, GFLOPS, BAND in GB/s)
#Job Output"
}


while getopts 'n:j:Hh' opt; do
  case "$opt" in
    n)
      outputfilename="$OPTARG"
      ;;

    j)
      jobname="$OPTARG"
      ;;

    H)
      print_report_header
      exit
      ;;
   
    ?|h)
      echo "Usage: $(basename $0) [-f < output filename>] [-j <jobreport filename>] [-H] [-h]"
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"


[ -z "$outputfilename" ] && exit 1
[ -z "$jobname" ] && exit 1




awk -v j=$jobname '{
    if(!support){support="MACRO";nofused=0;geometry="OLD"}
    if($0 ~ /\[SYSTEM\]\[0\]SIMD VECTORIZATION support enabled/){support="VECT" }
    if($0 ~ /\[SYSTEM\]\[0\]AVX2 support enabled/){support="AVX2"}
if($0 ~ /WITH_NEW_GEOMETRY/){geometry="NEW"}
if($0 ~ /\[OMP\]\[0\]Number of Threads requested =/) {printf "%s," ,$6}
if($0 ~ /\[SYSTEM\]\[0\]Gauge group:/) {split($3,a,"(");split(a[2],b,")"); printf "%s," ,b[1]}
if($0 ~ /\[SYSTEM\]\[0\]\[MPI_ID: 0\]\[MPI_size:/) {split($3,a,"]"); printf "%s," ,a[1]}
if($0 ~ /\[GEOMETRY_INIT\]\[0\]Global size is/) {split($4,a,"x"); printf "%s,%s,%s,%s," ,a[1],a[2],a[3],a[4]}
if($0 ~ /\[GEOMETRY_INIT\]\[0\]Local size is/) {split($4,a,"x"); printf "%s,%s,%s,%s," ,a[1],a[2],a[3],a[4]}
if($0 ~ /\[SETUP_RANDOM\]\[0\]RLXD/){split($2,a,","); split(a[2],b,"]") ; printf "%s,%s,%s," ,support,geometry,b[1]}
if($0 ~ /\[LA TEST\]\[0\]Flop per site =/){ printf "%s," ,$6}
if($0 ~ /\[LA TEST\]\[0\]Byte per site =/){ printf "%s," ,$6}
if($0 ~ /\[LA TEST\]\[0\]Dirac data movement =/){ printf "%s," ,$6}
if($0 ~ /\[LA TEST\]\[0\]Massless Diracoperator reps:/){ printf "%s,%s,%s,%s,%s," ,$5,$8,$11,$14,$17}
if($0 ~ /\[LA TEST\]\[0\]Massless fused Diracoperator reps:/){ printf "%s,%s,%s,%s,%s," ,$6,$9,$12,$15,$18;nofused=1}
if($0 ~ /\[SYSTEM\]\[0\]Process finalized./){ if(nofused==0){printf ",,,,,"};printf "%s\n",j}
    } ' $outputfilename
