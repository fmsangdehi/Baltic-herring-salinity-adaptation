## Install and run liftoff

#-------------------------------------------
## Create a conda environment and install liftoff

cd /crex/proj/snic2020-2-19/private/herring/users/fahime/
#mkdir envs
cd envs

ml load conda
conda create --prefix /crex/proj/snic2020-2-19/private/herring/users/fahime/envs/HE_conda_env
conda env list
conda activate /crex/proj/snic2020-2-19/private/herring/users/fahime/envs/HE_conda_env
conda list
conda install -c bioconda liftoff

#-------------------------------------------
## Create a directory and upload the gtf/gff3 file

#mkdir liftoff_hTGase_chr17
cd liftoff_hTGase_chr17

rsync -av /Users/FTG_chr17_Ensembl_corrected.gff3 fahim@rackham.uppmax.uu.se:/home/fahim/liftoff_hTGase_chr17

#-------------------------------------------
## Extract selected chromosomes from ref assembly

ref_assembly=/proj/snic2020-2-19/private/herring/assembly/Ch_v2.0.2.fasta

module load bioinfo-tools samtools

samtools faidx ${ref_assembly} chr17 > ref_chr17.fa
sed -i 's/^>chr17$/>17/' ref_chr17.fa

#-------------------------------------------
## Extract the target contigs from PacBio assemblies

cat >hap_list.txt

/proj/snic2020-2-19/private/herring/users/fahime/new_ref/pr_019_pool_herring.hifiasm0.19.7.default.bp.hap1.p_ctg.fasta
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/F1.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/F2.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/F3.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/F4.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/F5.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/F6.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/CS2.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/CS4.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/CS5.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/CS7.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/CS8.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/CS10.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/NSSH2.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap1/NSSH10.hap1.fa
/proj/snic2020-2-19/private/herring/users/fahime/new_ref/pr_019_pool_herring.hifiasm0.19.7.default.bp.hap2.p_ctg.fasta
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/F1.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/F2.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/F3.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/F4.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/F5.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/F6.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/CS2.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/CS4.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/CS5.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/CS7.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/CS8.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/CS10.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/NSSH2.hap2.fa
/proj/snic2020-2-19/private/herring/users/fahime/hifiasm_assemblies_renamed_new/hap2/NSSH10.hap2.fa

#-------------------------------------------

module load bioinfo-tools samtools seqtk

WORKDIR="/home/fahim/liftoff_hTGase_chr17"
paf_table="/home/fahim/hTGase_chr17/1_dotplot_of_target_region/2_obtain_target_seq/paf_table_simplified_noHeader_filtered.txt"
ID_list="${WORKDIR}/ID_list.txt"
hap_list="${WORKDIR}/hap_list.txt"

num_lines=$(wc -l < "${paf_table}")

for i in $(seq 1 ${num_lines}); do
    ID=$(awk -v i="${i}" 'NR==i {split($1, arr, "."); if(length(arr[1])==0) print $1; else print arr[1]}' "${paf_table}")
    line=$(grep -n ${ID} ${ID_list} | cut -d ':' -f 1)
    assembly=$(sed -n "${line}p" ${hap_list})
    contig=$(awk -F'\t' -v i="${i}" 'NR==i {print $7}' ${paf_table})
#    strand=$(awk "NR==$i {print \$6}" ${paf_table})

    mkdir ${ID}
    samtools faidx ${assembly} ${contig} > ${ID}/${ID}_FTG_chr17_contig.fa
done

#-------------------------------------------
## Run liftoff

ml load conda
conda activate /crex/proj/snic2020-2-19/private/herring/users/fahime/envs/HE_conda_env

WORKDIR="/home/fahim/liftoff_hTGase_chr17"
paf_table="/home/fahim/hTGase_chr17/1_dotplot_of_target_region/2_obtain_target_seq/paf_table_simplified_noHeader_filtered.txt"

num_lines=$(wc -l < "${paf_table}")

cd ${WORKDIR}

for i in $(seq 1 ${num_lines}); do
    ID=$(awk -v i="${i}" 'NR==i {split($1, arr, "."); if(length(arr[1])==0) print $1; else print arr[1]}' "${paf_table}")
    
    liftoff -g FTG_chr17_Ensembl_corrected.gff3 \
    -u ${ID}/unmapped_features.txt \
    -dir ${ID}/intermediate_files \
    -o ${ID}/${ID}_FTG_chr17.gff3 \
    ${ID}/${ID}_FTG_chr17_contig.fa \
    ref_chr17.fa
done

#-------------------------------------------
## Extraction of CDS using GFF3/GTF file by gffread (full-length CDS)
## and writing a protein fasta file with the translation of CDS for each record

module load bioinfo-tools gffread

WORKDIR="/home/fahim/liftoff_hTGase_chr17"
paf_table="/home/fahim/hTGase_chr17/1_dotplot_of_target_region/2_obtain_target_seq/paf_table_simplified_noHeader_filtered.txt"

num_lines=$(wc -l < "${paf_table}")

cd ${WORKDIR}

for i in $(seq 1 ${num_lines}); do
    ID=$(awk -v i="${i}" 'NR==i {split($1, arr, "."); if(length(arr[1])==0) print $1; else print arr[1]}' "${paf_table}")
    
    gffread ${ID}/${ID}_FTG_chr17.gff3 -g ${ID}/${ID}_FTG_chr17_contig.fa \
    -x - | sed "s/^>transcript:/>${ID}_/" >> FTG_chr17_CDS_full-length.fa
    
    gffread ${ID}/${ID}_FTG_chr17.gff3 -g ${ID}/${ID}_FTG_chr17_contig.fa \
    -y - | sed "s/^>transcript:/>${ID}_/" >> FTG_chr17_AA_full-length.fa
done

gffread FTG_chr17_Ensembl_corrected.gff3 -g ref_chr17.fa \
    -x - | sed "s/^>transcript:/>/" >> FTG_chr17_CDS_full-length.fa

gffread FTG_chr17_Ensembl_corrected.gff3 -g ref_chr17.fa \
    -y - | sed "s/^>transcript:/>/" >> FTG_chr17_AA_full-length.fa

#-------------------------------------------
## Split the multi-FASTA files into three files based on transcript IDs.

module load bioinfo-tools SeqKit

seqkit grep -r -p "ENSCHAT00020026568" FTG_chr17_CDS_full-length.fa -o FTG_chr17_CDS_gene1.fa
seqkit grep -r -p "ENSCHAT00020026651" FTG_chr17_CDS_full-length.fa -o FTG_chr17_CDS_gene2.fa
seqkit grep -r -p "ENSCHAT00020026742" FTG_chr17_CDS_full-length.fa -o FTG_chr17_CDS_gene3.fa

seqkit grep -r -p "ENSCHAT00020026568" FTG_chr17_AA_full-length.fa -o FTG_chr17_AA_gene1.fa
seqkit grep -r -p "ENSCHAT00020026651" FTG_chr17_AA_full-length.fa -o FTG_chr17_AA_gene2.fa
seqkit grep -r -p "ENSCHAT00020026742" FTG_chr17_AA_full-length.fa -o FTG_chr17_AA_gene3.fa

