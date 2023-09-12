#!/usr/bin/bash -l
#SBATCH -p short -N 1 -c 48 -C xeon --mem 48gb --out logs/AMPtk_ITS.%A.log

CPU=$SLURM_CPUS_ON_NODE

if [ ! $CPU ]; then
 CPU=2
fi

#AMPtk needs to be loaded in miniconda2 for UCR HPCC
#We'll need to unload miniconda3 and load miniconda2 before load AMPtk
module unload miniconda2
module load amptk
module load usearch

#Set up basename for all the output that will be generated
BASE=ECDRE_LR_Y0_1_ITS

#Chnage this to match your data folder name
INPUT=illumina

#Pre-preocessing steps will use `amptk illumia` command for demultiplexed PE reads
if [ ! -f $BASE.demux.fq.gz ]; then
 amptk illumina -i $INPUT --merge_method vsearch -f ITS1-F -r ITS2 --require_primer off -o $BASE --usearch usearch9 --cpus $CPU --rescue_forward on --primer_mismatch 2 -l 300
fi

if [ ! -f $BASE.otu_table.txt ];  then
 amptk cluster -i $BASE.demux.fq.gz -o $BASE --uchime_ref ITS --usearch usearch9 --map_filtered -e 0.9 --cpus $CPU
fi

if [ ! -f $BASE.otu_table.taxonomy.txt ]; then
 amptk taxonomy -f $BASE.cluster.otus.fa -i $BASE.otu_table.txt -d ITS
fi
