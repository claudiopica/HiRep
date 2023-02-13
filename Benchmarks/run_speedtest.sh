#!/bin/bash
MAXCORE=16
NSOCKETS=2
NODES=$1
SMALLEST_T=4
locall=8
localt=`seq $SMALLEST_T $(( SMALLEST_T*MAXCORE*NSOCKETS ))`

EXEC="speed_test_diracoperator"
outfile="out_n${NODES}"
reportfile="report_nodes${NODES}"

((MAXPROC = NSOCKETS * MAXCORE * NODES))
((MAXPROCPERNODE = NSOCKETS * MAXCORE))
((MINPROC = NSOCKETS * ( MAXCORE - 3 ) * NODES))

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
##SBATCH --qos=scavenger
# Launch MPI-based executable

# module load GCC/11.2.0   
# module load  OpenMPI/4.1.1
source /opt/ohpc/pub/oneAPI/setvars.sh

export I_MPI_DEBUG=5 
export OMP_DISPLAY_AFFINITY=TRUE 
unset MPI_DSM_OFF 
export MPI_DSM_VERBOSE=1 
export MPI_SHARED_VERBOSE=1 
export MPI_MEMMAP_VERBOSE=1 


export  OMP_PROC_BIND=close # How I am going to bind omp threads on to those places 
export  OMP_PLACES=cores    # Where I am going to place omp threads on the hardware, here is to cores
"
}



counter=0
ARRAYID=0
for paral in $(seq 0 0); do
    if ((MAXPROC / 2 ** paral >= 1)); then
        if ((paral == 0)); then
            NPX=1
            NPY=1
            NPZ=1
            (( MAXT = SMALLEST_T * MAXCORE * NSOCKETS * NODES ))
        fi

        if ((paral == 1)); then
            NPX=2
            NPY=1
            NPZ=1
            (( MAXT = SMALLEST_T * MAXCORE * NSOCKETS * NODES / 2  ))
        fi

        if ((paral == 2)); then
            NPX=2
            NPY=2
            NPZ=1
            (( MAXT = SMALLEST_T * MAXCORE * NSOCKETS * NODES / 4 ))
        fi

        if ((paral == 3)); then
            NPX=2
            NPY=2
            NPZ=2
            (( MAXT = SMALLEST_T * MAXCORE * NSOCKETS * NODES / 8  ))
        fi

    else

        break
    fi
 
    ((maxnpt = MAXPROC / (2 ** paral)))
    ((procpernode = MAXPROCPERNODE / (2 ** paral)))

    if ((procpernode == 0)); then
        break
    fi

    NPTLIST=$(seq 1 $maxnpt)

    for npt in ${NPTLIST[*]}; do
        (((((npt * NPX * NPY * NPZ) / 2) * 2) != npt * NPX * NPY * NPZ)) && continue
        (((((npt * NPX * NPY * NPZ) / NODES) * NODES) != npt * NPX * NPY * NPZ)) && continue
        ((npt * NPX * NPY * NPZ < NODES)) && continue
        for lct in ${localt[*]}; do
            ((glbt = npt * lct))
            ((glbt > MAXT)) && continue
            (( ( glbt / 2) * 2 != glbt )) && continue
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
            (( maxompproc =  ( maxnpt - npt )/NODES + 1 ))
            #omplist=($(seq 2 2 $maxompproc))
            #omplist=(1 ${omplist[@]})
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
EXECONTROL=\`tail -n1 ${outfile}_${NODES}_${npt}_${lct}_${paral}_${ompproc}_0 | grep -e "Process finalized" | wc -l\`
if (( EXECONTROL == 0 )) ; then
    mpirun -n $((npt * NPX * NPY * NPZ)) --ppn $(((npt * NPX * NPY * NPZ) / NODES))  ./$EXEC -i loc_speed_${NODES}_${npt}_${lct}_${paral}.in -o ${outfile}_${NODES}_${npt}_${lct}_${paral}_${ompproc}
    $parse_out_script -n ${outfile}_${NODES}_${npt}_${lct}_${paral}_${ompproc}_0 -j job.\${SLURM_JOB_ID}.out >>${locreportfile}
fi
EOF

			if (( counter == MAXARRAYSIZE - 1 )) ; then
			    counter=0
			    echo "fi" >> $jobmpi      
			fi
                #sbatch $jobmpi
            done
        done
    done
done
echo "fi" >> $jobmpi


echo "Jobarray should contain $(( counter + ( MAXARRAYSIZE - 1 )*(ARRAYID-1) )) elements in $ARRAYID submit scripts"
