#!/bin/bash

#test=(
#5 15 13 4 0 3 1 1 3 4 2 2 4 1 
#4 4 8 8 1 1 0 0 2 2 1 1 0 1
#4 4 8 8 1 1 0 0 2 2 1 1 0 0
#8 8 4 4 0 0 1 1 1 1 2 2 0 1
#)


ntests=100

for (( i=0; i<ntests; i++ )) ; do
#T X Y Z BT BX BY BZ npt npx npy npz myid myid_sign
   test[$((i*14+0))]=`expr $RANDOM % 10 + 2`; T=${test[$((i*14+0))]}
   test[$((i*14+1))]=`expr $RANDOM % 10 + 2`; X=${test[$((i*14+1))]}
   test[$((i*14+2))]=`expr $RANDOM % 10 + 2`; Y=${test[$((i*14+2))]}
   test[$((i*14+3))]=`expr $RANDOM % 10 + 2`; Z=${test[$((i*14+3))]}
   
   test[$((i*14+4))]=`expr $RANDOM % \( $T / 2 \)`; BT=${test[$((i*14+4))]}
   test[$((i*14+5))]=`expr $RANDOM % \( $X / 2 \)`; BX=${test[$((i*14+5))]}
   test[$((i*14+6))]=`expr $RANDOM % \( $Y / 2 \)`; BY=${test[$((i*14+6))]}
   test[$((i*14+7))]=`expr $RANDOM % \( $Z / 2 \)`; BZ=${test[$((i*14+7))]}
   
   if [ "${BT}" != "0" ] ; then test[$((i*14+8))]=`expr $RANDOM % 5 + 1` ; else test[$((i*14+8))]=1 ; fi ; npt=${test[$((i*14+8))]}
   if [ "${BX}" != "0" ] ; then test[$((i*14+9))]=`expr $RANDOM % 5 + 1` ; else test[$((i*14+9))]=1 ; fi ; npx=${test[$((i*14+9))]}
   if [ "${BY}" != "0" ] ; then test[$((i*14+10))]=`expr $RANDOM % 5 + 1` ; else test[$((i*14+10))]=1 ; fi ; npy=${test[$((i*14+10))]}
   if [ "${BZ}" != "0" ] ; then test[$((i*14+11))]=`expr $RANDOM % 5 + 1` ; else test[$((i*14+11))]=1 ; fi ; npz=${test[$((i*14+11))]}

   test[$((i*14+12))]=`expr $RANDOM % \( $npt \* $npx \* $npy \* $npz \)`
   test[$((i*14+13))]=`expr $RANDOM % 2`
done

echo ${test[@]}

#ntests=`expr ${#test[@]} / 14`

cp ../../Geometry/geometry_init.c.bak ../../Geometry/geometry_init.c
cp ../../Include/geometry.h.bak ../../Include/geometry.h

for (( i=0; i<ntests; i++ )) ; do
   awk '{par='${test[$((i*14+0))]}';if ($0 ~ /T=[0-9]+;/) {print "T="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+1))]}';if ($0 ~ /X=[0-9]+;/) {print "X="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+2))]}';if ($0 ~ /Y=[0-9]+;/) {print "Y="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+3))]}';if ($0 ~ /Z=[0-9]+;/) {print "Z="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+4))]}';if ($0 ~ /T_BORDER=[0-9]+;/) {print "T_BORDER="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+5))]}';if ($0 ~ /X_BORDER=[0-9]+;/) {print "X_BORDER="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+6))]}';if ($0 ~ /Y_BORDER=[0-9]+;/) {print "Y_BORDER="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c
   awk '{par='${test[$((i*14+7))]}';if ($0 ~ /Z_BORDER=[0-9]+;/) {print "Z_BORDER="par";"} else {print $0}}' ../../Geometry/geometry_init.c > tmp && mv tmp ../../Geometry/geometry_init.c

    awk '{par='${test[$((i*14+8))]}';if ($0 ~ /#define np_t *[0-9]+/) {print "#define np_t "par} else {print $0}}' ../../Include/geometry.h > tmp && mv tmp ../../Include/geometry.h
    awk '{par='${test[$((i*14+9))]}';if ($0 ~ /#define np_x *[0-9]+/) {print "#define np_x "par} else {print $0}}' ../../Include/geometry.h > tmp && mv tmp ../../Include/geometry.h
    awk '{par='${test[$((i*14+10))]}';if ($0 ~ /#define np_y *[0-9]+/) {print "#define np_y "par} else {print $0}}' ../../Include/geometry.h > tmp && mv tmp ../../Include/geometry.h
    awk '{par='${test[$((i*14+11))]}';if ($0 ~ /#define np_z *[0-9]+/) {print "#define np_z "par} else {print $0}}' ../../Include/geometry.h > tmp && mv tmp ../../Include/geometry.h
    awk '{par='${test[$((i*14+12))]}';if ($0 ~ /#define myid *[0-9]+/) {print "#define myid "par} else {print $0}}' ../../Include/geometry.h > tmp && mv tmp ../../Include/geometry.h
    awk '{par='${test[$((i*14+13))]}';if ($0 ~ /#define myid_sign *[0-9]+/) {print "#define myid_sign "par} else {print $0}}' ../../Include/geometry.h > tmp && mv tmp ../../Include/geometry.h

   make
   
   ./check_geometry_1
   pippo=`echo $?`;
   if [ "${pippo}" != "0" ] ; then exit ; fi
   
done
