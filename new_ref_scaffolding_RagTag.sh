#!/bin/bash
# RagTag workflow for scaffolding a contig-level assembly using a chromosome-scale reference

######################################
## Install Ragtag
######################################

# ml load conda
# conda env list #or: conda info --envs
### conda create -n ragtag python=3.8
# conda activate HE_conda_env
# conda install -c bioconda ragtag

######################################
## Scaffolding new ref assemblies
######################################

## Step 1: Estimate genomic divergence with Mash

ml load bioinfo-tools mash

ref="Ref.fa"
query_1="/pr_019_pool_herring.hifiasm0.19.7.default.bp.hap1.p_ctg.fasta"
query_2="/pr_019_pool_herring.hifiasm0.19.7.default.bp.hap2.p_ctg.fasta"

mash dist ${ref} ${query_1}
#${ref} ${query_1} 0.00590877	0	791/1000	#Mash distance is 0.0059 (Very low divergence (~0.59%))
mash dist ${ref} ${query_2}
#${ref} ${query_2} 0.00659353	0	771/1000

#-------------------------------------------------

## Step 2: Correct assembly errors in query contigs (optional)

#-------------------------------------------------

## Step 3: scaffold the contig-level assembly

ml load conda
conda activate HE_conda_env

REFERENCE="Ref.fa"
QUERY_1="/pr_019_pool_herring.hifiasm0.19.7.default.bp.hap1.p_ctg.fasta"
QUERY_2="/pr_019_pool_herring.hifiasm0.19.7.default.bp.hap2.p_ctg.fasta"
OUTPUT_DIR="/new_ref_haps_scaffolding"
THREADS=16

cd new_ref_haps_scaffolding

ragtag.py scaffold $REFERENCE $QUERY_1 \
    -o ${OUTPUT_DIR}/scaffolding_new.ref.hap1 \
    -t $THREADS \
    -r \
    --aligner minimap2 \
    --mm2-params '-x asm5'

mv ${OUTPUT_DIR}/scaffolding_new.ref.hap1/ragtag.scaffold.fasta ${OUTPUT_DIR}/scaffolding_new.ref.hap1/new.ref.hap1_scaffolded.fa
sed 's/^>Ref_\(.*\)_RagTag/>new_ref_hap1_\1/' ${OUTPUT_DIR}/scaffolding_new.ref.hap1/new.ref.hap1_scaffolded.fa > ${OUTPUT_DIR}/scaffolding_new.ref.hap1/new.ref.hap1_scaffolded_renamed.fa

ragtag.py scaffold $REFERENCE $QUERY_2 \
    -o ${OUTPUT_DIR}/scaffolding_new.ref.hap2 \
    -t $THREADS \
    -r \
    --aligner minimap2 \
    --mm2-params '-x asm5'

mv ${OUTPUT_DIR}/scaffolding_new.ref.hap2/ragtag.scaffold.fasta ${OUTPUT_DIR}/scaffolding_new.ref.hap2/new.ref.hap2_scaffolded.fa
sed 's/^>Ref_\(.*\)_RagTag/>new_ref_hap2_\1/' ${OUTPUT_DIR}/scaffolding_new.ref.hap2/new.ref.hap2_scaffolded.fa > ${OUTPUT_DIR}/scaffolding_new.ref.hap2/new.ref.hap2_scaffolded_renamed.fa

#-------------------------------------------------

## Step 4: Generate stats

ml load conda
conda activate HE_conda_env

OUTPUT_DIR="/new_ref_haps_scaffolding"

ragtag.py asmstats ${OUTPUT_DIR}/scaffolding_new.ref.hap1/new.ref.hap1_scaffolded_renamed.fa
ragtag.py asmstats ${OUTPUT_DIR}/scaffolding_new.ref.hap2/new.ref.hap2_scaffolded_renamed.fa

cat ${OUTPUT_DIR}/scaffolding_new.ref.hap1/ragtag.scaffold.stats
cat ${OUTPUT_DIR}/scaffolding_new.ref.hap2/ragtag.scaffold.stats

#-------------------------------------------------

## Step 5: Generate dotplots

#-------------------------------------------------

## Step 6: Run Spaln on the new scaffolded assemblies
## Inspect results in IGV

ml load bioinfo-tools spaln samtools

new_ref_h1="/new.ref.hap1_scaffolded_renamed.fa"
new_ref_h2="/new.ref.hap2_scaffolded_renamed.fa"

mkdir new_ref_h1_chr26
cd new_ref_h1_chr26

samtools faidx ${new_ref_h1} new_ref_hap1_chr26 > new_ref_h1_chr26.fa

time spaln -Wref.bkp -KP -t 8 new_ref_h1_chr26.fa
time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o new_ref_h1_chr26_HCE \
    ../HCE_query_protein_editName.fasta 2>&1 | tee new_ref_h1_chr26_HCE.log


cd /spaln_test
mkdir new_ref_h2_chr26
cd new_ref_h2_chr26

samtools faidx ${new_ref_h2} new_ref_hap2_chr26 > new_ref_h2_chr26.fa

time spaln -Wref.bkp -KP -t 8 new_ref_h2_chr26.fa
time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o new_ref_h2_chr26_HCE \
    ../HCE_query_protein_editName.fasta 2>&1 | tee new_ref_h2_chr26_HCE.log

