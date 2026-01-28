## Make AA-based ML phylogenetic trees with IQ-TREE (only_FXIII)

ml load bioinfo-tools iqtree

CDS_aln="/home/fahim/CDS_aln_concatenated_edited_noUnique_only_FXIII_renamed_subset.fasta"

iqtree2 -s ${CDS_aln} --seqtype NT2AA --prefix AA_based_ML_tree -T AUTO --threads-max 8 -m MFP -B 1000

