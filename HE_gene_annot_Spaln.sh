## Gene annotation with Spaln

#########################
## Step1 - Extract chr 26
#########################

ml load bioinfo-tools samtools gnuparallel

export fa_files_dir="RagTag_genomes"
fa_files=(Ref_v2 Ref_v3 new_ref_hap{1,2} CS{2,4,5,7,8,10}_hap{1,2} BS{1,2,3,4,5,6}_hap{1,2} NSSH{2,10}_hap{1,2})

parallel -j 5 --env fa_files_dir '
  mkdir -p {}
  samtools faidx "${fa_files_dir}/{}.fa" {}_chr26 > {}/{}_chr26.mfa
' ::: "${fa_files[@]}"

#--------------
## Extract sprat chromosome 1 (OY735285.1)

sprat_genome_dir="/proj/snic2020-2-19/private/herring/users/fahime/ref_genomes_downloaded_original/Sprat_Darwins_tree_of_life_NCBI"
mkdir -p sprat
samtools faidx "${sprat_genome_dir}/GCA_963457725.1_fSprSpr1.1_genomic.fna" OY735285.1 > sprat/sprat_chr1.mfa


##########################
## Step2 - make database
## Step3 - gene annotation
##########################

ml load bioinfo-tools spaln

parent_dir="/home/fahim/annot"
cd ${parent_dir}

fa_files=(Ref_v2 Ref_v3 new_ref_hap{1,2} CS{2,4,5,7,8,10}_hap{1,2} BS{1,2,3,4,5,6}_hap{1,2} NSSH{2,10}_hap{1,2})

for i in ${fa_files[@]}; do
  cd ${i}
  time spaln -Wref.bkp -KP -t 8 ${i}_chr26.mfa
  
  #Annotate HCE
  time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o ${i}_chr26_HCE \
    ../query_protein_seqs/HCE_ENSCHAT00020074750_full-length_AA.fasta 2>&1 | tee ${i}_chr26_HCE.log
  mv ${i}_chr26_HCE.O0 ${i}_chr26_HCE.O0.gff3
  
  #Annotate hsp70
   time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o ${i}_hsp70 \
    ../query_protein_seqs/hsp70_ENSCHAT00020077596_full-length_AA.fasta 2>&1 | tee ${i}_hsp70.log
  mv ${i}_hsp70.O0 ${i}_hsp70.O0.gff3
  
  #Annotate dnaja3b
  time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o ${i}_dnaja3b \
    ../query_protein_seqs/dnaja3b_full-length_AA.fasta 2>&1 | tee ${i}_dnaja3b.log
  mv ${i}_dnaja3b.O0 ${i}_dnaja3b.O0.gff3
  
  #Annotate CCNB1
  time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o ${i}_CCNB1 \
    ../query_protein_seqs/CCNB1_full-length_AA.fasta 2>&1 | tee ${i}_CCNB1.log
  mv ${i}_CCNB1.O0 ${i}_CCNB1.O0.gff3
  
  cd ${parent_dir}
done

#--------------
## Gene annotation and mapping for sprat

cd sprat
time spaln -Wref.bkp -KP -t 8 sprat_chr1.mfa

#Annotate HCE
time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o sprat_chr1_HCE \
    ../query_protein_seqs/HCE_ENSCHAT00020074750_full-length_AA.fasta 2>&1 | tee sprat_chr1_HCE.log
mv sprat_chr1_HCE.O0 sprat_chr1_HCE.O0.gff3

#Annotate hsp70
time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o sprat_hsp70 \
    ../query_protein_seqs/hsp70_ENSCHAT00020077596_full-length_AA.fasta 2>&1 | tee sprat_hsp70.log
mv sprat_hsp70.O0 sprat_hsp70.O0.gff3

#Annotate dnaja3b
time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o sprat_dnaja3b \
    ../query_protein_seqs/dnaja3b_full-length_AA.fasta 2>&1 | tee sprat_dnaja3b.log
mv sprat_dnaja3b.O0 sprat_dnaja3b.O0.gff3

#Annotate CCNB1
time spaln -Q7 -LS -M500 -dref -O0,1,2,3,4,6,7 -t 8 -o sprat_CCNB1 \
    ../query_protein_seqs/CCNB1_full-length_AA.fasta 2>&1 | tee sprat_CCNB1.log
mv sprat_CCNB1.O0 sprat_CCNB1.O0.gff3


##########################
## Step4 - seq extraction
##########################

## Extraction of CDS using GFF3/GTF file by gffread (full-length CDS)
## and writing a protein fasta file with the translation of CDS for each record

ml load bioinfo-tools gffread

parent_dir="/home/fahim/annot"
cd ${parent_dir}

mkdir all_fasta_output

fa_files=(Ref_v2 Ref_v3 new_ref_hap{1,2} CS{2,4,5,7,8,10}_hap{1,2} BS{1,2,3,4,5,6}_hap{1,2} NSSH{2,10}_hap{1,2})

for i in ${fa_files[@]}; do
    gffread ${i}/${i}_chr26_HCE.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -x - | sed "s/^>/>${i}_chr26_/" >> all_fasta_output/HCE_chr26_CDS_full-length.fa
    gffread ${i}/${i}_chr26_HCE.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -y - | sed "s/^>/>${i}_chr26_/" >> all_fasta_output/HCE_chr26_AA_full-length.fa
    
    gffread ${i}/${i}_hsp70.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -x - | sed "s/^>/>${i}_/" >> all_fasta_output/hsp70_CDS_full-length.fa
    gffread ${i}/${i}_hsp70.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -y - | sed "s/^>/>${i}_/" >> all_fasta_output/hsp70_AA_full-length.fa
    
    gffread ${i}/${i}_dnaja3b.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -x - | sed "s/^>/>${i}_/" >> all_fasta_output/dnaja3b_CDS_full-length.fa
    gffread ${i}/${i}_dnaja3b.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -y - | sed "s/^>/>${i}_/" >> all_fasta_output/dnaja3b_AA_full-length.fa
    
    gffread ${i}/${i}_CCNB1.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -x - | sed "s/^>/>${i}_/" >> all_fasta_output/CCNB1_CDS_full-length.fa
    gffread ${i}/${i}_CCNB1.O0.gff3 -g ${i}/${i}_chr26.mfa \
    -y - | sed "s/^>/>${i}_/" >> all_fasta_output/CCNB1_AA_full-length.fa
done

#--------------
## Extraction of CDS and protein from sprat

gffread sprat/sprat_chr1_HCE.O0.gff3 -g sprat/sprat_chr1.mfa \
    -x - | sed "s/^>/>sprat_chr1_/" >> all_fasta_output/HCE_chr26_CDS_full-length.fa
gffread sprat/sprat_chr1_HCE.O0.gff3 -g sprat/sprat_chr1.mfa \
    -y - | sed "s/^>/>sprat_chr1_/" >> all_fasta_output/HCE_chr26_AA_full-length.fa

gffread sprat/sprat_hsp70.O0.gff3 -g sprat/sprat_chr1.mfa \
    -x - | sed "s/^>/>sprat_/" >> all_fasta_output/hsp70_CDS_full-length.fa
gffread sprat/sprat_hsp70.O0.gff3 -g sprat/sprat_chr1.mfa \
    -y - | sed "s/^>/>sprat_/" >> all_fasta_output/hsp70_AA_full-length.fa

gffread sprat/sprat_dnaja3b.O0.gff3 -g sprat/sprat_chr1.mfa \
    -x - | sed "s/^>/>sprat_/" >> all_fasta_output/dnaja3b_CDS_full-length.fa
gffread sprat/sprat_dnaja3b.O0.gff3 -g sprat/sprat_chr1.mfa \
    -y - | sed "s/^>/>sprat_/" >> all_fasta_output/dnaja3b_AA_full-length.fa

gffread sprat/sprat_CCNB1.O0.gff3 -g sprat/sprat_chr1.mfa \
    -x - | sed "s/^>/>sprat_/" >> all_fasta_output/CCNB1_CDS_full-length.fa
gffread sprat/sprat_CCNB1.O0.gff3 -g sprat/sprat_chr1.mfa \
    -y - | sed "s/^>/>sprat_/" >> all_fasta_output/CCNB1_AA_full-length.fa


#################################
## Step5 - Collect all gff3 files
#################################
## Copy all gff3 files into a single directory.

parent_dir="/home/fahim/annot"
cd ${parent_dir}

mkdir all_gff3_files_HCE_chr26 all_gff3_files_hsp70 all_gff3_files_dnaja3b all_gff3_files_CCNB1

fa_files=(Ref_v2 Ref_v3 new_ref_hap{1,2} CS{2,4,5,7,8,10}_hap{1,2} BS{1,2,3,4,5,6}_hap{1,2} NSSH{2,10}_hap{1,2})

for i in ${fa_files[@]}; do
  cp ${i}/${i}_chr26_HCE.O0.gff3 all_gff3_files_HCE_chr26
  cp ${i}/${i}_hsp70.O0.gff3 all_gff3_files_hsp70
  cp ${i}/${i}_dnaja3b.O0.gff3 all_gff3_files_dnaja3b
  cp ${i}/${i}_CCNB1.O0.gff3 all_gff3_files_CCNB1
done

#--------------
## Copy gff3 file of sprat into the same directory.

cp sprat/sprat_chr1_HCE.O0.gff3 all_gff3_files_HCE_chr26
cp sprat/sprat_hsp70.O0.gff3 all_gff3_files_hsp70
cp sprat/sprat_dnaja3b.O0.gff3 all_gff3_files_dnaja3b
cp sprat/sprat_CCNB1.O0.gff3 all_gff3_files_CCNB1

