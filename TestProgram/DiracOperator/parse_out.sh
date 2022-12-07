#!/bin/bash
awk '{
if($0 ~ /\[OMP\]\[0\]Number of Threads requested =/) {printf "%s," ,$6}
if($0 ~ /\[SYSTEM\]\[0\]Gauge group:/) {split($3,a,"(");split(a[2],b,")"); printf "%s," ,b[1]}
if($0 ~ /\[SYSTEM\]\[0\]\[MPI_ID: 0\]\[MPI_size:/) {split($3,a,"]"); printf "%s," ,a[1]}
if($0 ~ /\[GEOMETRY_INIT\]\[0\]Global size is/) {split($4,a,"x"); printf "%s,%s,%s,%s," ,a[1],a[2],a[3],a[4]}
if($0 ~ /\[GEOMETRY_INIT\]\[0\]Local size is/) {split($4,a,"x"); printf "%s,%s,%s,%s," ,a[1],a[2],a[3],a[4]}
if($0 ~ /\[SETUP_RANDOM\]\[0\]RLXD/){split($2,a,","); split(a[2],b,"]") ; printf "%s," ,b[1]}
if($0 ~ /\[LA TEST\]\[0\]Flop per site =/){ printf "%s," ,$6}
if($0 ~ /\[LA TEST\]\[0\]Byte per site =/){ printf "%s," ,$6}
if($0 ~ /\[LA TEST\]\[0\]Dirac data movement =/){ printf "%s," ,$6}
if($0 ~ /\[LA TEST\]\[0\]Massless fused Diracoperator reps:/){ printf "%s,%s,%s,%s,%s," ,$6,$9,$12,$15,$18}
if($0 ~ /\[LA TEST\]\[0\]Massless Diracoperator reps:/){ printf "%s,%s,%s,%s,%s\n" ,$5,$8,$11,$14,$17}
    } ' $1