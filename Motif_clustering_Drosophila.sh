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

# set the 'argument' variable equal to the first argument on the command line, else foo
# default function imported from common.sh
argument="$(default $1 foo)"

#####
# Step 1: Prepare motif databases
#####

# TF motif models were downloaded from https://resources.aertslab.org/papers/iregulon/motifColl-10k-all-public.tar.gz in cluster-buster format
# (see here for details on the format: https://aertslab.org/#data-resources-all, under *SUPPLEMENTARY MATERIAL TO THE IREGULON PAPER*)

# create Markov Background Model - order 3
fasta-get-markov -n m /groups/stark/genomes/dm3/dm3.fa dm3.3-order.markov

# convert PWM models to MEME format
chen2meme singletons/bergman*cb -bg dm3.3-order.markov > all.dbs.meme
chen2meme singletons/cisbp*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/flyfactorsurvey*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/homer*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/jaspar*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/stark*cb -bg dm3.3-order.markov >> all.dbs.meme
chen2meme singletons/idmmpmm*cb -bg dm3.3-order.markov >> all.dbs.meme

#####
# Step 2: Compute pair-wise motif similarity
#####

tomtom \
	-dist kullback \
	-motif-pseudo 0.1 \
	-text \
	-min-overlap 1 \
	all.dbs.meme all.dbs.meme \
> tomtom.all.txt

# remove last lines that have details
head -n -4 tomtom.all.txt > tomtom.all.treated.txt

#####
# Step 3: Hierarchically cluster motifs by similarity in R (Motif_clustering.R)
#####

# postrun command, imported from common.sh
postrun

