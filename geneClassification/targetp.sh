#!/bin/bash
#SBATCH --array=0-4
#SBATCH --time=24:00:00
#SBATCH --job-name=targetp
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --error=targetp_%A_%a.err
#SBATCH --output=targetp_%A_%a.out
echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

fastaFile=$(python getLine.py fastas.fofn $SLURM_ARRAY_TASK_ID)
len="$((${#fastaFile}-6))"
fh=${fastaFile:0:$len}
echo $fh
fofn=$(python splitFasta.py $fastaFile)
y=$(python lineNum.py $fofn)
x=0
for ((i=$x;i<=$y;i++));do arrayFile=$(python getLine.py $fofn $i);len="$((${#arrayFile}-6))";fileString=${arrayFile:0:$len};targetp -P $arrayFile > $fileString.targetp.out;rm $arrayFile;done
cat $fh*.targetp.out > $fh.targetp.out
rm $fh.*.targetp.out
python targetpOutput.py $fh.targetp.out $fastaFile > $fh.processed.targetp.out
