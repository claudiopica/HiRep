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
   test[$((i*14+0))]=`expr $RANDOM % 3 + 3`; T=${test[$((i*14+0))]}
   test[$((i*14+1))]=`expr $RANDOM % 3 + 3`; X=${test[$((i*14+1))]}
   test[$((i*14+2))]=`expr $RANDOM % 3 + 3`; Y=${test[$((i*14+2))]}
   test[$((i*14+3))]=`expr $RANDOM % 3 + 3`; Z=${test[$((i*14+3))]}
   
   test[$((i*14+4))]=`expr $RANDOM % 3 + 1`; BT=${test[$((i*14+4))]} 
   test[$((i*14+5))]=`expr $RANDOM % 3 + 1`; BX=${test[$((i*14+5))]} 
   test[$((i*14+6))]=`expr $RANDOM % 3 + 1`; BY=${test[$((i*14+6))]} 
   test[$((i*14+7))]=`expr $RANDOM % 3 + 1`; BZ=${test[$((i*14+7))]} 

   test[$((i*14+8))]=`expr 2 \* \( $RANDOM % 3 \) `; dT=${test[$((i*14+8))]}
   test[$((i*14+9))]=`expr 2 \* \( $RANDOM % 3 \) `; dX=${test[$((i*14+9))]}
   test[$((i*14+10))]=`expr 2 \* \( $RANDOM % 3 \) `; dY=${test[$((i*14+10))]}
   test[$((i*14+11))]=`expr 2 \* \( $RANDOM % 3 \) `; dZ=${test[$((i*14+11))]}


   
   (( T=BT*(T + (BT*T)%2 )+dT  ))
   (( X=BX*(X + (BX*X)%2 )+dX  ))
   (( Y=BY*(Y + (BY*Y)%2 )+dY  ))
   (( Z=BZ*(Z + (BZ*Z)%2 )+dZ  ))
   

   test[$((i*14+0))]=$T
   test[$((i*14+1))]=$X
   test[$((i*14+2))]=$Y
   test[$((i*14+3))]=$Z

   test[$((i*14+4))]=$BT
   test[$((i*14+5))]=$BX
   test[$((i*14+6))]=$BY
   test[$((i*14+7))]=$BZ


   test[$((i*14+13))]=`expr  $BT \* $BX \* $BY \* $BZ `
done


#ntests=`expr ${#test[@]} / 14`


for (( i=0; i<ntests; i++ )) ; do

echo "GLB_T =  ${test[$((i*14+0))]} //Global T size">test_input
echo "GLB_X =  ${test[$((i*14+1))]}		">>test_input
echo "GLB_Y =  ${test[$((i*14+2))]}		">>test_input
echo "GLB_Z =  ${test[$((i*14+3))]}		">>test_input
echo "NP_T =   ${test[$((i*14+4))]}		">>test_input
echo "NP_X =   ${test[$((i*14+5))]}		">>test_input
echo "NP_Y =   ${test[$((i*14+6))]}      		">>test_input
echo "NP_Z =   ${test[$((i*14+7))]}		">>test_input
echo "level = 1			">>test_input
echo "seed = 12345              ">>test_input


(( ${test[$((i*14+13))]} > 18 )) && continue;
cat test_input

   openmpirun -np  ${test[$((i*14+13))]} ./check_geometry_1
   pippo=`echo $?`;
   if [ "${pippo}" != "0" ] ; then exit ; fi
   
   echo "Ok!"
   
done
