#!/bin/bash -l
#
# Run this file using 'qsub template.qsub argument'

# All lines starting with "#$" are SGE qsub commands

# Specify a project to run under
#$ -P wunderl

# Join standard output and error to a single file
#$ -j y

# Name the file where to redirect standard output and error
#$ -o logs/template.log

# Send an email when the job finishes
#$ -m eas

# Request a large memory node, this will affect your queue time, but it's better to overestimate
#$ -l mem_per_core=8G

# Request more cores, this will affect your queue time, make sure your program supports multithreading, or it's a waste
# -pe omp 1

# Now we write the script that the compute node will work on.

# load modules - MEME suite
module load perl/5.28.1
module load python3/3.10.12
module load openmpi/4.1.5
module load meme/5.5.5

# source the common functionality for logging (common.sh needs to be in the same directory as this qsub)
source ./common.sh

# prerun command, imported from common.sh
prerun


# INPUT PARAMETERS
cb_database="$(default $1 ${datafolder}/learn_starrmotifs/data/ict2022.cb)"

# OUTPUT PARAMETERS
meme_database="$(default $2 ${datafolder}/learn_starrmotifs/data/ict2022.meme)"
genome_fa=""
markov_model=""
tomtom_out=""

#####
# Step 1: Prepare motif databases
#####

# TF motif models were downloaded from https://resources.aertslab.org/papers/iregulon/motifColl-10k-all-public.tar.gz in cluster-buster format
# (see here for details on the format: https://aertslab.org/#data-resources-all, under *SUPPLEMENTARY MATERIAL TO THE IREGULON PAPER*)

# create Markov Background Model - order 3
fasta-get-markov -n m ${genome_fa} ${markov_model}

# convert PWM models to MEME format
chen2meme ${database} -bg ${markov_model} > ${meme_database}

#####
# Step 2: Compute pair-wise motif similarity
#####

tomtom \
	-dist kullback \
	-motif-pseudo 0.1 \
	-text \
	-min-overlap 1 \
	${meme_database} ${meme_database} \
> ${tomtom_out}

# remove last lines that have details
head -n -4 ${tomtom_out} > ${tomtom_out}

#####
# Step 3: Hierarchically cluster motifs by similarity in R (Motif_clustering.R)
#####

# postrun command, imported from common.sh
postrun

