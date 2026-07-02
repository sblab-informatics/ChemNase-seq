# CnT data processing for the G4 ligand-PROTAC projects of paper "Targeted Degradation of G-Quadruplex Binding Cellular Proteins Using G4-Ligand-Based PROTACs"
sn=$1 # sample name
$out_dir="raw_data"
R1="unspecified/*.UnspecifiedIndex.*.s_1.r_1.fq.gz"
R2="unspecified/*.UnspecifiedIndex.*.s_1.r_2.fq.gz"

# Step 1: Demultiplex the sequencing data
demuxFQ -b unmatched_barcodes.fq -c -d -i -e -t 1 -r 0.01 -R -l 9 -o $out_dir -s demul_summary_r1.txt demul_index_r1 $R1
demuxFQ -b unmatched_barcodes.fq -c -d -i -e -t 1 -r 0.01 -R -l 9 -o $out_dir -s demul_summary_r2.txt demul_index_r2 $R2

# Step 2: quality control
fastqc -o fastqc_results raw_data/${sn}r1.fq.gz raw_data/${sn}r2.fq.gz
multiqc fastqc_results

# Step 3: trim the sequencing data with fastp
fastp -w 32 -p -q 25 \
    -i raw_data/${sn}r1.fq.gz -I raw_data/${sn}r2.fq.gz \
    -o trimmed/${sn}r1.fq.gz -O trimmed/${sn}r2.fq.gz \
    -h trim_logs/${sn}html \
    -j trim_logs/${sn}json

# Step4: alignment with bowtie2
R1="trimmed/${sn}r1.fq.gz"
R2="trimmed/${sn}r2.fq.gz"

idx="bowtie2_index/hg38"
log="bowtie2_logs"

bowtie2 -p 32 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 \
    -X 1500 -x $idx -1 $R1 -2 $R2 2> $log/${sn}metrics.txt \
    | samtools view -bS -@ 32 -o raw/${sn}bam

# Step 5: remove the reads from mitochondria
samtools view -@ 32 -h raw/${sn}bam | grep -v chrM | samtools sort -@ 32 -o rmmt/${sn}rmMT.bam

# Step 6: mark duplicates with Picard
picard="java -jar picard.jar"
${picard} MarkDuplicates -I rmmt/${sn}rmMT.bam -O dupmarked/${sn}markdup.bam -M dupmarked/${sn}dupMark.txt

# Step 7: reads filtering
samtools view -@ 32 -h -b -F 1804 -f 2 -q 10 dupmarked/${sn}markdup.bam > filtered/${sn}bam

# Step 8: generate count per million reads (CPM) normalized profile
samtools index -@ 32 filtered/${sn}bam
bamCoverage -p 32 -bs 10 --normalizeUsing CPM --effectiveGenomeSize 2913022398 \
    --extendReads --minFragmentLength 10 --maxFragmentLength 1000 \
    -b filtered/${sn}bam -o bigwigs/${sn}cpm.bw

# Step 9: generate mapping profile based on fragments
samtools sort -@ 32 -n -O bam -o sorted/${sn}sortName.bam filtered/${sn}bam

bedtools bamtobed -i sorted/${sn}sortName.bam -bedpe > beds/${sn}bed
awk '$1==$4 && $6-$2 < 1500 && $6-$2 > 10 {print $0}' beds/${sn}bed > beds/${sn}clean.bed
cut -f 1,2,6 beds/${sn}clean.bed | sort -k1,1 -k2,2n -k3,3n  > beds/${sn}fragments.bed

fcount=$(cat beds/${sn}fragments.bed | wc -l)
scale_factor=$(awk "BEGIN {print "1000000"/"${fcount}"}")

gsize="hg38.chrom.sizes"
bedtools genomecov -bg -scale $scale_factor -i beds/${sn}fragments.bed -g $gsize > bedgraphs/${sn}fragments.bedgraph
bedGraphToBigWig bedgraphs/${sn}fragments.bedgraph $gsize bigwigs/${sn}fragment.bw

# Step 10: peak calling with SEACR
seacr="SEACR-1.3/SEACR_1.3.sh"
$seacr bedgraphs/${sn}fragments.bedgraph 0.01 non stringent peaks/${sn}fragments.top0.01

