#!/bin/bash
set -e

LOGFILE="variant_calling.log"
exec > >(tee -a "$LOGFILE") 2>&1

# Define a logging function with a timestamp.
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log "Starting viral variant calling pipeline with BWA-MEM..."

# Download reference genome
log "Downloading reference genome..."
wget --quiet --show-progress \
     --output-document=reference.fa \
     "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001802.1&rettype=fasta"

log "Reference genome downloaded."

# Prefetch and fastq-dump steps.
log "Prefetching SRR33784444..."
prefetch SRR33784444 --progress
log "Prefetch complete."

log "Running fastq-dump to split FASTQ files..."
fastq-dump --split-files SRR33784444/
log "fastq-dump complete."

# Run FastQC on raw FASTQ files.
log "Running FastQC on raw FASTQ files..."
fastqc SRR33784444_1.fastq SRR33784444_2.fastq
log "FastQC on raw FASTQ files complete."

# Create directory for trimmed reads and run Trimmomatic.
log "Creating directory for trimmed reads and running Trimmomatic..."
mkdir -p Trimmed
cd Trimmed

# Copy reference to working directory
cp ../reference.fa .

trimmomatic PE ../SRR33784444_1.fastq ../SRR33784444_2.fastq \
    SRR33784444_1_paired.fastq SRR33784444_1_unpaired.fastq \
    SRR33784444_2_paired.fastq SRR33784444_2_unpaired.fastq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
    SLIDINGWINDOW:4:20 MINLEN:50
log "Trimmomatic processing complete."

# Run FastQC on trimmed files.
log "Running FastQC on trimmed FASTQ files..."
fastqc SRR33784444_1_paired.fastq SRR33784444_2_paired.fastq SRR33784444_1_unpaired.fastq SRR33784444_2_unpaired.fastq
log "FastQC on trimmed FASTQ files complete."

# Build BWA index for reference genome
log "Building BWA index for reference genome..."
bwa index reference.fa
log "BWA index building complete."

# Alignment with BWA-MEM (better for viral sequences)
log "Running BWA-MEM alignment..."
bwa mem -t 4 -M -R "@RG\tID:Sample\tSM:Sample\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    reference.fa SRR33784444_1_paired.fastq SRR33784444_2_paired.fastq > alignment_bwa.sam
log "BWA-MEM alignment complete."

# Check alignment statistics
log "Checking alignment statistics..."
samtools flagstat alignment_bwa.sam | tee alignment_stats.txt
log "Alignment statistics saved to alignment_stats.txt"

cd ..
mv Trimmed/alignment_bwa.sam .
mv Trimmed/alignment_stats.txt .

# Convert SAM to BAM, sort, and index using SAMtools
log "Converting SAM to BAM..."
samtools view -bS alignment_bwa.sam > temp.bam
log "Sorting BAM file..."
samtools sort temp.bam -o sorted.bam
log "Removing temporary BAM file..."
rm temp.bam
log "Indexing sorted BAM file..."
samtools index sorted.bam

# Verify read groups are present
log "Verifying read groups in sorted BAM file..."
samtools view -H sorted.bam | grep "@RG" || log "WARNING: No read groups found in BAM file"

# Check coverage statistics
log "Calculating coverage statistics..."
samtools depth sorted.bam > coverage.txt
awk '{sum+=$3; count++} END {print "Average coverage:", sum/count}' coverage.txt | tee coverage_summary.txt
log "Coverage statistics saved."

# Mark duplicates using GATK MarkDuplicates.
log "Marking duplicates with GATK..."
gatk MarkDuplicates \
    -I sorted.bam \
    -M metrics.txt \
    -O temp_unique_reads.bam \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true
log "Duplicates marked."

# Add/fix read groups to ensure GATK compatibility
log "Adding/fixing read groups for GATK compatibility..."
gatk AddOrReplaceReadGroups \
    -I temp_unique_reads.bam \
    -O unique_reads.bam \
    -RGID Sample \
    -RGSM Sample \
    -RGPL ILLUMINA \
    -RGLB lib1 \
    -RGPU unit1 \
    --CREATE_INDEX true
log "Read groups added/fixed."

# Clean up temporary file
rm temp_unique_reads.bam temp_unique_reads.bai

# Verify read groups in final BAM file
log "Verifying read groups in final BAM file..."
samtools view -H unique_reads.bam | grep "@RG"
log "Read group verification complete."

# Create reference files needed for GATK
log "Creating sequence dictionary for reference..."
gatk CreateSequenceDictionary -R reference.fa -O reference.dict
log "Indexing reference with samtools..."
samtools faidx reference.fa

# Call variants using GATK HaplotypeCaller with viral-specific parameters
log "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R reference.fa \
    -I unique_reads.bam \
    -O output.vcf \
    --native-pair-hmm-threads 4 \
    --standard-min-confidence-threshold-for-calling 10 \
    --minimum-mapping-quality 20 \
    --max-reads-per-alignment-start 0
log "Variant calling complete."

# Filter variants using GATK VariantFiltration with viral-specific filters
log "Filtering variants with GATK VariantFiltration..."
gatk VariantFiltration \
    -R reference.fa \
    -V output.vcf \
    -O filtered_output.vcf \
    --filter-expression "QD < 2.0" \
    --filter-name "LowQD" \
    --filter-expression "FS > 60.0" \
    --filter-name "HighFS" \
    --filter-expression "MQ < 30.0" \
    --filter-name "LowMQ" \
    --filter-expression "DP < 10" \
    --filter-name "LowDepth" \
    --filter-expression "QUAL < 30.0" \
    --filter-name "LowQual"
log "Variant filtration complete."

# Create a high-confidence variant set
log "Creating high-confidence variant set..."
gatk SelectVariants \
    -R reference.fa \
    -V filtered_output.vcf \
    -O high_confidence_variants.vcf \
    --exclude-filtered true
log "High-confidence variants selected."

# Additional viral-specific analysis
log "Generating variant summary statistics..."
bcftools stats filtered_output.vcf > variant_stats.txt
log "Variant statistics saved to variant_stats.txt"

# Create consensus sequence
log "Creating consensus sequence..."
bcftools consensus -f reference.fa filtered_output.vcf > consensus_sequence.fa
sed -i 's/>.*/>Patient_Sample_Consensus/' consensus_sequence.fa
log "Consensus sequence created as consensus_sequence.fa"

# Generate coverage plot data
log "Generating coverage data for visualization..."
samtools depth sorted.bam | awk '{print $2"\t"$3}' > coverage_plot.txt
log "Coverage plot data saved to coverage_plot.txt"

# Final summary
log "Generating final summary..."
echo "=== Viral Variant Calling Pipeline Summary ===" > pipeline_summary.txt
echo "Date: $(date)" >> pipeline_summary.txt
echo "" >> pipeline_summary.txt
echo "Alignment Statistics:" >> pipeline_summary.txt
cat alignment_stats.txt >> pipeline_summary.txt
echo "" >> pipeline_summary.txt
echo "Coverage Statistics:" >> pipeline_summary.txt
cat coverage_summary.txt >> pipeline_summary.txt
echo "" >> pipeline_summary.txt
echo "Variant Statistics:" >> pipeline_summary.txt
grep -E "^SN" variant_stats.txt >> pipeline_summary.txt
echo "" >> pipeline_summary.txt
echo "Output Files Generated:" >> pipeline_summary.txt
echo "- sorted.bam: Aligned reads" >> pipeline_summary.txt
echo "- unique_reads.bam: Deduplicated reads with proper read groups" >> pipeline_summary.txt
echo "- output.vcf: Raw variants" >> pipeline_summary.txt
echo "- filtered_output.vcf: Filtered variants" >> pipeline_summary.txt
echo "- high_confidence_variants.vcf: High-confidence variants only" >> pipeline_summary.txt
echo "- consensus_sequence.fa: Consensus genome sequence" >> pipeline_summary.txt
echo "- coverage_plot.txt: Coverage data for plotting" >> pipeline_summary.txt

log "Viral variant calling pipeline finished successfully."
log "Check pipeline_summary.txt for detailed results."