#! /bin/bash

# Gabriel Hoffman
# March 2, 2015
#
# countVirusReads.sh
#
# IN: FASTQ of read that don't map to human reference
# OUT: Count number of reads mapping to each viral genome

# Method: Assemble contigs with Trinity's Inchworm, blast contigs against NCBI viral genomes

# Invoke:
# countVirusReads.sh [unmapped fastq (text or gzip)] [SEQ_DATABASE] [NCPU] [Work folder]

# Prereq's:

#########################
# Create virus database #
#########################

# NCBI:
# http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&opt=Virus
# Click on "Sequence Info" (top right)
# Send to: file : FASTA
# Build database
# makeblastdb -in /hpc/users/hoffmg01/annotationData/NCBI_virus_complete_genome.fasta -dbtype nucl

################
# Load modules #
################

module purge
module load trinity/2013-08-14
module load bowtie2/2.1.0
module load tophat/2.0.6
module load samtools/0.1.19 # NOTE: using samtools/1.* causes error in prep_reads
module load bamutil/1.0.11
module load blast/2.2.26
module load R/3.1.1

# Arguments 
###########

FASTQ1=$1
FASTQ2=$2
SEQ_DATABASE=$3
NCPU=$4
WORK=$5

# Print arguments
#################

echo "Arguments:"
echo "FASTQ1: $FASTQ1"
echo "FASTQ2: $FASTQ2"
echo "SEQ_DATABASE: $SEQ_DATABASE"
echo "NCPU: $NCPU"
echo -e "WORK: $WORK\n"

# Check arguments
#################

if [ ! -e $FASTQ ]; 
then
	echo "FASTQ does not exist";
	exit 1;
fi

if [ ! -e $SEQ_DATABASE ]; 
then
	echo "SEQ_DATABASE does not exist";
	exit 1;
fi

if [ ! -e $WORK ]; 
then
	echo "WORK does not exist";
	exit 1;
fi


# Make folder
mkdir -p $WORK $WORK/scripts

# Analysis 
###########

# Assemble contigs
Trinity.pl --seqType fq --JM 3G --CPU $NCPU --no_run_chrysalis --no_run_butterfly --left $FASTQ1 --right $FASTQ2 --output $WORK/trinity_inchworm

# make output directory
mkdir -p $WORK/trinity_inchworm_tophat_unmapped

# # Build Bowtie index on Inchworm contigs
bowtie2-build $WORK/trinity_inchworm/inchworm.K25.L25.DS.fa $WORK/trinity_inchworm/inchworm.K25.L25.DS

# # run TopHat to map reads from FASTQ to Inchworm contigs
# tophat -p $NCPU --bowtie1 -o $WORK/trinity_inchworm_tophat_unmapped $WORK/trinity_inchworm/inchworm.K25.L25.DS $FASTQ
tophat -p $NCPU -o $WORK/trinity_inchworm_tophat_unmapped $WORK/trinity_inchworm/inchworm.K25.L25.DS $FASTQ1 $FASTQ2

# # index BAM file of mapped reads
# samtools index $WORK/trinity_inchworm_tophat_unmapped/accepted_hits.bam

# Write header to BALSTN results file
echo -e 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' | tr ' ' '\t' > $WORK/blastn_trinity_on_seqdb.tsv

# Blast against virus database
blastn -query $WORK/trinity_inchworm/inchworm.K25.L25.DS.fa -db ${SEQ_DATABASE} -evalue .03 -outfmt 6 | sort -k1 >> $WORK/blastn_trinity_on_seqdb.tsv

# Write R script that combines statistics from Inchwor, BLAST and TopHat
echo "
data = read.table( '$WORK/blastn_trinity_on_seqdb.tsv', header=TRUE, sep='\t', stringsAsFactors=FALSE)

seqAnnot = c()
nreads = c()
for(i in 1:nrow(data) ){
	cat('\r', i, ' / ', nrow(data), sep='' )
	cmd = paste( \"grep '\", data\$sseqid[i], \"' ${SEQ_DATABASE}\", sep='')
	
	id = system( cmd, intern=TRUE)
	id = paste( unlist(strsplit(id, ' '))[-1], collapse=' ')
	id = gsub(', complete genome', '', id)
	id = gsub(' complete genome', '', id)

	seqAnnot[i] = id

	cmd = paste(\"grep '\", data\$qseqid[i], \"' $WORK/trinity_inchworm/inchworm.K25.L25.DS.fa | cut -d' ' -f 3\", sep='')
	res = system( cmd, intern=TRUE)

	nreads[i] = as.numeric( res )
}
data\$seqAnnot = seqAnnot
data\$nreads_inchworm = nreads

# Add counts from TopHat
mappedReads = system('samtools view $WORK/trinity_inchworm_tophat_unmapped/accepted_hits.bam | cut -f3',intern=TRUE)
mappedReads = table(mappedReads)

idx = match( data\$qseqid, names(mappedReads))
data\$nreads_tophat = mappedReads[idx]

write.table( data, '$WORK/quantify_virus.tsv', row.names=FALSE, quote=FALSE, sep='\t')
" > $WORK/scripts/combine_results.R

# Run script
Rscript $WORK/scripts/combine_results.R

# Compare Inchworm counts with TopHat counts
# Note that a lot of the counts are rRNA contamination
# Mapping to the human reference should absorb these reads
#	so they won't be in the unmapped reads files
# plot( dataUniq$nreads, mappedReads[idx], ylab="# TopHat reads", xlab="# Inchworm reads")

# fit = lm( mappedReads[idx] ~ dataUniq$nreads )
# abline(coef(fit)[1], coef(fit)[2], lwd=2, col="grey80", lty=2)
# value = format( summary(fit)$r.squared, digits=2)
# legend("topleft", legend=bquote(R^2 == .(value)), bty='n', cex=1.4)


