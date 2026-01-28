## Alignment, tree reconstruction, dN/dS analysis

##------------------
## Step 1: Alignment
##------------------

## Align the AA sequences with Mafft

module load bioinfo-tools MAFFT

mafft TGase.all.pep.v2_from_Andreas_n_herring_only_FXIII_modified.fasta  > alignment/AA_aln_input_order.fasta
mafft --reorder TGase.all.pep.v2_from_Andreas_n_herring_only_FXIII_modified.fasta  > alignment/AA_aln_aligned_order.fasta


## Align the DNA (CDS) sequences with PAL2NAL (along with AA alignment) to create codon-based DNA alignment.

PAL2NAL_DIR="/home/fahim/my_tools/pal2nal.v14"

${PAL2NAL_DIR}/pal2nal.pl alignment/AA_aln_input_order.fasta TGase.all.cds.v2_from_Andreas_n_herring_only_FXIII_modified.fasta -output fasta > alignment/CDS_aln.fasta
${PAL2NAL_DIR}/pal2nal.pl alignment/AA_aln_input_order.fasta TGase.all.cds.v2_from_Andreas_n_herring_only_FXIII_modified.fasta -output paml > alignment/CDS_aln.paml


##--------------------------
## Step 2: phylogenetic tree
##--------------------------

## Construct phylogenetic tree with IQtree

ml load bioinfo-tools iqtree

AA_aln="alignment/AA_aln_input_order.fasta"

#iqtree2 -s ${AA_aln} --seqtype AA --prefix IQtree/AA_based_ML_tree_BS -T AUTO --threads-max 8 -m MFP -B 1000
#iqtree2 -s ${AA_aln} --seqtype AA --prefix IQtree/AA_based_ML_tree -T 8 -m MFP


CDS_aln="alignment/CDS_aln.fasta"

#iqtree2 -s ${CDS_aln} --seqtype DNA --prefix IQtree/nt_based_ML_tree_BS -T AUTO --threads-max 8 -m MFP -B 1000

iqtree2 -s ${CDS_aln} --seqtype CODON --prefix IQtree/codon_based_ML_tree_BS -T AUTO --threads-max 8 -m MFP -B 1000


##----------------------
## Step 3: PAML (codeml)
##----------------------
## 1. Based on IQtree

cat > PAML_ML_tree_run_04.sh

#!/bin/bash

#SBATCH -A uppmax2025-2-114
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J PAML_ML_tree_run_04
#SBATCH --mail-type=all
#SBATCH --mail-user=fhmohamadnejad@gmail.com

ml load bioinfo-tools paml

codeml codeml_IQtree.ctl 2>&1 | tee codeml.log

##----------------------------
sbatch PAML_ML_tree_run_04.sh

