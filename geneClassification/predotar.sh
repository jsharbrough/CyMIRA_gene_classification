#!/bin/bash
#SBATCH --array 0-4
#SBATCH --time=24:00:00
#SBATCH --job-name=predotar
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --error=predotar_%A_%a.err
#SBATCH --output=predotar_%A_%a.out

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
fastaFile=$(python getLine.py fastas.fofn $SLURM_ARRAY_TASK_ID)
len="$((${#fastaFile}-6))"
fh=${fastaFile:0:$len}
echo $fh
module load jdk/1.8.0
java -jar predotar.jar p $fastaFile > $fh.predotar.out
