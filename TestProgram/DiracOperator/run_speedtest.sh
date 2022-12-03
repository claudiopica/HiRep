#!/bin/bash
MAXCORE=64
NSOCKETS=2
NODES=$1
localt=(4 6 8 10)
locall=6
MAXT=120

EXEC="speed_test_diracoperator"
outfile="out_0_n${NODES}"
reportfile="report_nodes${NODES}.csv"
jobmpi="job_n${NODES}.mpi"

((MAXPROC = NSOCKETS * MAXCORE * NODES))
((MAXPROCPERNODE = NSOCKETS * MAXCORE))

parse_out_script="./parse_out.sh"

rm -rf $jobmpi

echo "export  OMP_PROC_BIND=close" > $jobmpi
echo "export  OMP_PLACES=cores" >> $jobmpi
echo "scontrol show hostname \$SLURM_JOB_NODELIST | perl -ne 'chomb; print \"\$_\"x1'> myhostall_nodes${NODES}">> $jobmpi


cat <<EOF >$reportfile
#Number of Threads
#Gauge group
#MPI_size
#Global size (4x)
#Local size  (4x)
#RLXD
#Flop per site
#Byte per site
#Dirac data movement
#Massless fused Diracoperator (reps, data size in kb, time in msec, GFLOPS, BAND in GB/s)
#Massless Diracoperator (reps, data size in kb, time in msec, GFLOPS, BAND in GB/s)
EOF

for paral in $(seq 0 3); do
    if ((MAXPROC / 2 ** paral >= 1)); then
        if ((paral == 0)); then
            NPX=1
            NPY=1
            NPZ=1
        fi

        if ((paral == 1)); then
            NPX=2
            NPY=1
            NPZ=1
        fi

        if ((paral == 2)); then
            NPX=2
            NPY=2
            NPZ=1
        fi

        if ((paral == 3)); then
            NPX=2
            NPY=2
            NPZ=2
        fi

    else

        break
    fi

    ((maxnpt = MAXPROC / (2 ** paral) ))
    ((procpernode = MAXPROCPERNODE / (2 ** paral) ))

    if ((procpernode == 0)); then
        break
    fi

    ((freeprocpernode = MAXPROCPERNODE - procpernode))
    ((maxompproc = freeprocpernode + 1))

    NPTLIST=$(seq 1 $maxnpt)

    for npt in ${NPTLIST[*]}; do
        (( ( ( ( npt * NPX * NPY * NPZ ) / 2 ) * 2 ) != npt*NPX*NPY*NPZ )) && continue
        (( ( ( ( npt * NPX * NPY * NPZ ) / NODES ) * NODES ) != npt*NPX*NPY*NPZ )) && continue
        (( npt * NPX * NPY * NPZ  < NODES )) && continue        
        for lct in ${localt[*]}; do
            ((glbt = npt * lct ))
            ((glbt > MAXT)) && continue
            cat <<EOF >loc_speed_${NODES}_${npt}_${lct}_${paral}.in
// Global variables 
GLB_T = $glbt //Global T size
GLB_X = $((locall * NPX))
GLB_Y = $((locall * NPY))
GLB_Z = $((locall * NPZ))
NP_T = $npt
NP_X = ${NPX}
NP_Y = ${NPY}
NP_Z = ${NPZ}
rlx_level = 1
rlx_seed = 12345

// Replicas
N_REP = 1

//Logger levels (default = -1)
log:default = -1
log:inverter = -1
log:forcestat = 0
EOF
            omplist=(`seq 2 2 $maxompproc`)
            omplist=(1 ${omplist[@]})
            for ompproc in ${omplist[@]}; do
                (( ompproc * npt * NPX * NPY * NPZ >  MAXPROC )) && continue
                echo export OMP_NUM_THREADS=$ompproc >>$jobmpi
                echo rm -f ${outfile}_${npt}_${lct}_${ompproc} >>$jobmpi
                echo mpirun  -hostfile myhostall_nodes${NODES} -n $(( npt * NPX * NPY * NPZ )) --map-by node ./$EXEC -i loc_speed_${NODES}_${npt}_${lct}_${paral}.in -o ${outfile}_${NODES}_${npt}_${lct}_${paral}_${ompproc} >>$jobmpi
                echo -e "$parse_out_script ${outfile}_${NODES}_${npt}_${lct}_${paral}_${ompproc} >>$reportfile\n" >>$jobmpi
            done
        done
    done
done
