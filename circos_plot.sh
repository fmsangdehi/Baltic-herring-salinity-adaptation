## Circular synteny plot with Circos based on NUCmer show-coords output, along with HE gene highlight

##-----------------------
## Create Karyotype Files
##-----------------------

## The columns are: chr - ID LABEL START END COLOR

herring_ref_v3_fasta_index="Ref_v3.fa.fai"
sprat_ref_fasta_index="GCA_963457725.1_fSprSpr1.1_genomic.fna.fai"

grep "^Ref_v3_chr" ${herring_ref_v3_fasta_index} | awk -v OFS="\t" '{print "chr", "-", $1, $1, 0, $2, $1}' >> karyotype_herring.txt
grep "^OY" ${sprat_ref_fasta_index} | grep -v "OY735306.1" | awk -v OFS="\t" '{print "chr", "-", $1, $1, 0, $2, $1}' >> karyotype_sprat.txt


##-----------------------------------------
## Convert show-coords to Circos links file
##-----------------------------------------

show_coords_dir="1.2_nucmer_default_settings"
rsync -av ${show_coords_dir}/show_coords.coords fahim@rackham.uppmax.uu.se:/home/fahim/Circos

awk -v OFS="\t" 'NR>5 {print $14, $1, $2, $15, $3, $4}' show_coords.coords > links.txt


##--------------------------------------
## Filter only chromosomes in links file
##--------------------------------------

## This is to include only chromosomes in karyotype.txt

cut -f3 karyotype_herring.txt | sort | uniq >> valid_chromosomes.txt
cut -f3 karyotype_sprat.txt | sort | uniq >> valid_chromosomes.txt

awk 'NR==FNR {valid[$1]=1; valid[$4]=1; next} valid[$1] && valid[$4]' valid_chromosomes.txt links.txt > links.filtered.txt


##-----------------------------
## Fix and reduce the links file
##-----------------------------

## This is to:
## 1) Sort the coordinates to ensure start < end for each side of the link

awk '{
    if ($2 > $3) { temp=$2; $2=$3; $3=temp }
    if ($5 > $6) { temp=$5; $5=$6; $6=temp }
    print $1, $2, $3, $4, $5, $6
}' OFS="\t" links.filtered.txt > links.fixed.txt


## 2) Reduce the number of lines in the links file as Circos accepts maximum of [25000].

awk '$3-$2 > 4000 && $6-$5 > 4000' links.fixed.txt > links.fixed.reduced.txt

##--------------------------
## Highlight and label genes
##--------------------------

## Format: chromosome start end [options]

#OY735286.1	40888000	40938000	fill_color=128,128,128 #OY735286.1 gene was not captured

cat > genes.txt

Ref_v3_chr5	18266000	18316000	fill_color=0,204,153
Ref_v3_chr20	1780000	2090000	fill_color=0,102,102
Ref_v3_chr22	25820000	26140000	fill_color=51,153,255
Ref_v3_chr22	27069000	27099000	fill_color=51,153,255
Ref_v3_chr24	6412000	6452000	fill_color=0,102,204
Ref_v3_chr24	7836000	7866000	fill_color=0,102,204
Ref_v3_chr26	5870000	6770000	fill_color=0,0,204
OY735285.1	24981000	25051000	fill_color=128,128,128
OY735285.1	27493000	27523000	fill_color=128,128,128
OY735287.1	15815000	15915000	fill_color=128,128,128
OY735290.1	25200000	25900000	fill_color=128,128,128
OY735290.1	27430000	27460000	fill_color=128,128,128
OY735291.1	26564000	26594000	fill_color=128,128,128
OY735293.1	17520000	17610000	fill_color=128,128,128
OY735298.1	10098000	10128000	fill_color=128,128,128


##--------------------------
## Circos configuration file
##--------------------------

cat > circos.conf


karyotype = karyotype_herring.txt,karyotype_sprat.txt
chromosomes_units = 1000000
chromosomes_order = Ref_v3_chr1,Ref_v3_chr18,Ref_v3_chr2,Ref_v3_chr3,Ref_v3_chr12,Ref_v3_chr4,Ref_v3_chr15,Ref_v3_chr10,Ref_v3_chr5,\
                    Ref_v3_chr6,Ref_v3_chr7,Ref_v3_chr24,Ref_v3_chr26,Ref_v3_chr8,Ref_v3_chr9,Ref_v3_chr11,\
                    Ref_v3_chr13,Ref_v3_chr14,Ref_v3_chr16,Ref_v3_chr17,Ref_v3_chr19,Ref_v3_chr22,Ref_v3_chr25,\
                    Ref_v3_chr20,Ref_v3_chr21,Ref_v3_chr23,\
                    OY735287.1,OY735290.1,OY735288.1,OY735304.1,OY735301.1,OY735298.1,\
                    OY735297.1,OY735296.1,OY735303.1,OY735293.1,OY735285.1,\
                    OY735299.1,OY735292.1,OY735294.1,OY735286.1,OY735300.1,OY735295.1,OY735291.1,OY735289.1,OY735305.1,OY735302.1

chromosomes_reverse = OY735285.1,OY735286.1,OY735287.1,OY735290.1,OY735291.1,OY735292.1,OY735293.1,OY735294.1,OY735295.1,OY735296.1,OY735297.1,OY735299.1,OY735300.1,OY735301.1,OY735302.1,OY735305.1
#chromosomes_reverse = OY735289.1,OY735288.1,OY735304.1,OY735303.1,OY735298.1

#chromosomes_color   = Ref_v3_chr1=red,Ref_v3_chr2=orange,Ref_v3_chr3=white

<image>
dir = .
file = circos_synteny.png
png = yes
svg = yes
radius = 1500p
angle_offset = -90
</image>

<ideogram>
<spacing>
default = 0.005r
</spacing>
radius = 0.8r
thickness = 20p
fill = yes
stroke_thickness = 2p
stroke_color = black
show_label = yes
label_font = default
label_radius = dims(ideogram,radius) + 0.02r
label_size = 36
label_parallel = no
</ideogram>

<ticks>
show_ticks = no
show_tick_labels = no
</ticks>

<links>
<link>
file = links.fixed.reduced.txt
ribbon = yes
thickness = 2
radius = 0.92r
bezier_radius = 0r

<rules>
<rule>
condition = 1
color = eval(substr(var(chr1), 6))
</rule>
</rules>

</link>
</links>

<highlights>
<highlight>
type = highlight
file = genes.txt
r0   = 0.93r
r1   = 0.99r
z    = 0
#stroke_thickness = 2
</highlight>
</highlights>

<colors>
_chr1  = 255,102,102     # Red  
_chr18 = 255,153,102     # Orange  
_chr2  = 255,204,102     # Amber  
_chr3  = 255,255,102     # Yellow  
_chr12 = 204,255,102     # Yellow-green  
_chr4  = 153,255,102     # Light green  
_chr15 = 102,255,102     # Green  
_chr10 = 102,255,153     # Greenish-cyan  
_chr5  = 102,255,204     # Aqua  
_chr6  = 102,255,255     # Cyan  
_chr7  = 102,204,255     # Light sky blue  
_chr24 = 102,153,255     # Sky blue  
_chr26 = 102,102,255     # Blue  
_chr8  = 153,102,255     # Violet  
_chr9  = 204,102,255     # Purple  
_chr11 = 255,102,255     # Magenta  
_chr13 = 255,102,204     # Rose  
_chr14 = 255,102,153     # Deep pink  
_chr16 = 255,102,128     # Coral pink  
_chr17 = 255,102,178     # Hot pink  
_chr19 = 204,153,255     # Lavender  
_chr22 = 153,204,255     # Baby blue  
_chr25 = 102,204,204     # Teal  
_chr20 = 102,153,153     # Muted teal  
_chr21 = 153,153,102     # Olive khaki  
_chr23 = 192,128,128     # Warm beige
gray = 192,192,192
</colors>

<<include /sw/apps/circos/0.69-9/rackham/etc/housekeeping.conf>>
<<include /sw/apps/circos/0.69-9/rackham/etc/colors_fonts_patterns.conf>>

max_links* = 150000


##-----------
## Run Circos
##-----------

ml load bioinfo-tools circos

circos -conf circos.conf
