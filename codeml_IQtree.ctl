seqfile = ../alignment/CDS_aln.paml      * alignment file name
treefile = codon_based_ML_tree_BS_rerooted_ape_edited.nwk         * tree structure file name
outfile = results.txt        * output file name

noisy = 3                         * 0, 1, 2, 3, 9: how much rubbish on the screen
verbose = 1                       * 1: detailed output, 0: concise output
runmode = 0                       * 0: user-defined tree

seqtype = 1                       * 1: codons
ndata = 1
icode = 0                         * 0: universal genetic code
cleandata = 0

model = 1                         * 0: one dN/dS ratio
NSsites = 0                       * 0: no site-specific models
CodonFreq = 2                     * 2: F3x4 model for codon frequencies
estFreq = 0
clock = 0
fix_omega = 0                     * 1: fix omega, 0: estimate
omega = 0.5                         * initial or fixed omega
