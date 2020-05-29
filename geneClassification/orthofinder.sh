#!/bin/bash
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --job-name=orthofinder
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --error=orthofinder.err
#SBATCH --output=orthofinder.out

orthofinder -t 24 -og -f /path/to/folder/containing/protein/fastas

#Need to include protein fastas (primary isoform only) from your species of interest as well as from Arabidopsis (Araport 11 fasta included in distribution)
