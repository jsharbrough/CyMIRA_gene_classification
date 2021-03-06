#!/bin/bash
#SBATCH --array 0-5
#SBATCH --time=24:00:00
#SBATCH --job-name=ipsort
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --error=ipsort_%A_%a.err
#SBATCH --output=ipsort_%A_%a.out
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

fastaFile=$(python getLine.py fastas.fofn $SLURM_ARRAY_TASK_ID)
len="$((${#fastaFile}-6))"
fh=${fastaFile:0:$len}
echo $fh
ipsort -F -i $fastaFile > $fh.ipsort.out
python ipsortOutput.py $fh.ipsort.out $fastaFile > $fh.processed.ipsort.out

