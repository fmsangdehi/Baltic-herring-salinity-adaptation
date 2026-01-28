#mkdir all_herring_hTGase_tree
cd all_herring_hTGase_tree

seq_dir="/Users/sprat_as_outgroup"
rsync -av ${seq_dir}/hTGase_AA_full-length.fa fahim@rackham.uppmax.uu.se:/home/fahim/all_herring_hTGase_tree
rsync -av ${seq_dir}/hTGase_CDS_full-length.fa fahim@rackham.uppmax.uu.se:/home/fahim/all_herring_hTGase_tree

##----------------------------
## Alignment
##----------------------------
## 1) Align the AA sequences with Mafft

module load bioinfo-tools MAFFT

mafft hTGase_AA_full-length.fa  > AA_aln_input_order.fasta
#mafft --reorder hTGase_AA_full-length.fa  > AA_aln_aligned_order.fasta

##----------------------------
## 2) Align the DNA (CDS) sequences with PAL2NAL (along with AA alignment) to create codon-based DNA alignment.

PAL2NAL_DIR="/home/fahim/my_tools/pal2nal.v14"

${PAL2NAL_DIR}/pal2nal.pl AA_aln_input_order.fasta hTGase_CDS_full-length.fa -output fasta > CDS_aln.fasta
${PAL2NAL_DIR}/pal2nal.pl AA_aln_input_order.fasta hTGase_CDS_full-length.fa -output paml > CDS_aln.paml

##----------------------------
## Tree construction
##----------------------------
## 1) Make an ML tree with IQ-TREE (codon-based phylogenetic tree)

ml load bioinfo-tools iqtree

CDS_aln="CDS_aln.fasta"

iqtree2 -s ${CDS_aln} --seqtype CODON1 --prefix herring_codon_tree -T AUTO --threads-max 8 -m MFP
iqtree2 -s ${CDS_aln} --seqtype CODON1 --prefix herring_codon_tree -T 6 -m MFP

##----------------------------
## 2) Make an ML tree with IQ-TREE (amino-acid-based phylogenetic tree)

ml load bioinfo-tools iqtree

AA_aln="AA_aln_input_order.fasta"

iqtree2 -s ${AA_aln} --seqtype AA --prefix herring_AA_tree -T AUTO --threads-max 8 -m MFP
iqtree2 -s ${AA_aln} --seqtype AA --prefix herring_AA_tree -T 4 -m MFP

