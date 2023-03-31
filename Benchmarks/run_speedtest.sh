#!/bin/bash
MAXCORE=64
NSOCKETS=2
NODES=$1

#SIZE of the lattice per MPI process in X,Y,Z
locall=8

#Largest local VOLUME per PROCESS (OMP*MPI) Change to modify the size of largest local volumes
(( LLVOL=8**4 ))
#Smallest local VOLUME per PROCESS (OMP*MPI) Change to modify the size of smallest local volumes
(( MLVOL=4**4 ))

#paral= 0 -> NPT>=1 NPX=NPY=NPZ=1
#paral= 1 -> NPT>=1 NPX=2 NPY=NPZ=1
#paral= 2 -> NPT>=1 NPX=NPY=2 NPZ=1
#paral= 3 -> NPT>=1 NPX=NPY=NPZ=2
paral=3 


(( MAXGLB_T = ( LLVOL * NODES * MAXCORE *NSOCKETS)/(2**paral * locall**3))) #MAX Global T
(( MINGLB_T = ( MLVOL * NODES * MAXCORE *NSOCKETS)/(2**paral * locall**3))) #MIN Global T


localt=`seq 4 $(( MAXGLB_T / NODES ))` #range of possible values of local T (4 is the smallest meaningful lattice size)

EXEC="speed_test_diracoperator"
outfile="out_n${NODES}"
reportfile="report_nodes${NODES}"

((MAXPROC = NSOCKETS * MAXCORE * NODES))#Maximal number of processes (MPI*OMPI) available
((MINPROC = NSOCKETS * (MAXCORE-3) * NODES )) #Minimal number of processes (MPI*OMPI) used
((MAXPROCPERNODE = NSOCKETS * MAXCORE))#Maximal number of processes (MPI*OMPI) available per node
((MAXNPT = MAXPROC / (2 ** paral))) #max NP_T


if (( MAXPROCPERNODE / (2 ** paral) == 0)); then
    echo "Not enough cores for the choosen parallelization"
    exit
fi


if (( (MAXPROCPERNODE / (2 ** paral)) * (2 ** paral) != MAXPROCPERNODE )); then
    echo "The number of cores per processor is not compatible with the choosen parallelization strategy"
    exit
fi

parse_out_script="./parse_out.sh"


MAXARRAYSIZE=`scontrol show config | grep -e MaxArraySize | awk '{print $3}'`

rm -rf job_${NODES}_*.mpi #report_nodes${NODES}_*.csv

print_slurm_header () {
    echo "#!/bin/bash
#SBATCH -J op_$1               # Job name
#SBATCH -o job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH -N $NODES                  # Total number of nodes requested
#SBATCH --time 0:02:00           # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --exclusive
#SBATCH --constraint=hm2
# Launch MPI-based executable

. ~/spack/share/spack/setup-env.sh 

#module load openmpi-4.1.4-clang-15.0.7-wuqcqaw
module load mpich-4.1-clang-15.0.7-o3jgla4
export LD_LIBRARY_PATH=:/home/srahman/spack/opt/spack/linux-almalinux8-zen/gcc-8.5.0/gcc-12.2.0-dplyzyl6twqs6bjthkfv3dcljl4u5hkp/lib64/

export I_MPI_DEBUG=5 
export OMP_DISPLAY_AFFINITY=TRUE 
unset MPI_DSM_OFF 
export MPI_DSM_VERBOSE=1 
export MPI_SHARED_VERBOSE=1 
export MPI_MEMMAP_VERBOSE=1 

export UCX_NET_DEVICES=mlx5_2:1
export UCX_RC_MLX5_TM_ENABLE=y

export OMP_AFFINITY_FORMAT=\"Thread Affinity: %.8i %.8n %.15A %.12H\"

export OMP_PROC_BIND=close # How I am going to bind omp threads on to those places 
export OMP_PLACES=cores    # Where I am going to place omp threads on the hardware, here is to cores
"
}


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


counter=0
ARRAYID=0

NPTLIST=$(seq 1 $MAXNPT)

for npt in ${NPTLIST[*]}; do
    (((((npt * NPX * NPY * NPZ) / NODES) * NODES) != npt * NPX * NPY * NPZ)) && continue
    ((npt * NPX * NPY * NPZ < NODES)) && continue
    for lct in ${localt[*]}; do
        ((glbt = npt * lct))
        ((glbt > MAXGLB_T)) && continue
        ((glbt < MINGLB_T)) && continue
        (( ( glbt / 2) * 2 != glbt )) && continue
        cat <<EOF >loc_speed_${NODES}_${npt}_${lct}_${locall}_${paral}.in
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
        ((maxompproc = (MAXPROC - npt * NPX * NPY * NPZ) / NODES + 1))
        omplist=($(seq 1 $maxompproc))
        
        for ompproc in ${omplist[@]}; do
            ((ompproc * npt * NPX * NPY * NPZ > MAXPROC)) && continue
            ((ompproc * npt * NPX * NPY * NPZ < MINPROC)) && continue
            (( counter+=1 ))
	    if (( counter == 1)) ; then
		(( ARRAYID += 1 ))
		jobmpi="job_${NODES}_${ARRAYID}.mpi"
		locreportfile="${reportfile}_${ARRAYID}.csv"

		print_slurm_header ${NODES}_${ARRAYID} > $jobmpi
		[ ! -f "$locreportfile" ] && $parse_out_script -H > $locreportfile
		ifelse="if"
            else
		ifelse="elif"
            fi
	    cat <<EOF >>$jobmpi
$ifelse (( SLURM_ARRAY_TASK_ID == $counter )) ; then
export OMP_NUM_THREADS=$ompproc
EXECONTROL=\`tail -n1 ${outfile}_${NODES}_${npt}_${lct}_${locall}_${paral}_${ompproc}_0 | grep -e "Process finalized" | wc -l\`
if (( EXECONTROL == 0 )) ; then
    srun -n $((npt * NPX * NPY * NPZ)) -c$ompproc  ./$EXEC -i loc_speed_${NODES}_${npt}_${lct}_${locall}_${paral}.in -o ${outfile}_${NODES}_${npt}_${lct}_${locall}_${paral}_${ompproc}
    $parse_out_script -n ${outfile}_${NODES}_${npt}_${lct}_${locall}_${paral}_${ompproc}_0 -j job.\${SLURM_JOB_ID}.out >>${locreportfile}
fi
EOF

	    if (( counter == MAXARRAYSIZE - 1 )) ; then
		counter=0
		echo "fi" >> $jobmpi      
	    fi

        done
    done
done
echo "fi" >> $jobmpi


echo "Jobarray should contain $(( counter + ( MAXARRAYSIZE - 1 )*(ARRAYID-1) )) elements in $ARRAYID submit scripts"
