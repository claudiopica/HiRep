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

   echo "GLB_T =  ${test[$((i*14+0))]} //Global T size">test_input
   echo "GLB_X =  ${test[$((i*14+1))]}		">>test_input
   echo "GLB_Y =  ${test[$((i*14+2))]}		">>test_input
   echo "GLB_Z =  ${test[$((i*14+3))]}		">>test_input
   echo "NP_T =   ${test[$((i*14+4))]}		">>test_input
   echo "NP_X =   ${test[$((i*14+5))]}		">>test_input
   echo "NP_Y =   ${test[$((i*14+6))]}      		">>test_input
   echo "NP_Z =   ${test[$((i*14+7))]}		">>test_input
   echo "// Replicas">>test_input
   echo "N_REP = 1">>test_input
   echo "//Logger levels (default = -1)">>test_input
   echo "log:default = -1">>test_input
   echo "log:inverter = -1">>test_input
   echo "log:forcestat = 0">>test_input
   echo "// Random generator">>test_input
   echo "rlx_level = 1">>test_input
   echo "rlx_seed = 13813">>test_input
   echo "rlx_start = new">>test_input
   echo "rlx_state = rand_state">>test_input
   echo "rlx_store = 0">>test_input

   (( ${test[$((i*14+13))]} > 18 )) && ((i--)) && continue;
   #cat test_input
   export OMP_NUM_THREADS=1
   rm -f out_0
   mpirun-openmpi-gcc11 --mca shmem posix --oversubscribe -np ${test[$((i*14+13))]} ./check_geometry_1 -i ./test_input
   if [ "$?" != "0" ] ; then exit ; fi

   echo $i"/$ntests : Ok!"

done
