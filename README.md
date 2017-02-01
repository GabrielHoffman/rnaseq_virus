# Documentation RNASeq repository #

A collection of scripts implementing analyses for RNA-Seq data, created by Gabriel Hoffman at the Icahn School of Medcine at Mount Sinai

# Detecting viral reads in FASTQ files #

In addition to capturing the expression of human transcripts, RNA-Seq FASTQ files can also contain reads from viral genomes.  Here I present a simple pipeline to quantify the expression of viral transcripts.

1) Get FASTQ of reads that don't map to human reference.  This subtract human reads that we don't want to consider in the virus analysis

    # write to this file prefix
    OUTFILE=~/test
    
    module load star/2.4.0g1
    NCPU=20
    REF=/sc/orga/projects/PBG/REFERENCES/hg19
    GENOME=$REF/star/2.4.0d/Homo_sapiens.GRCh37.ensembl70.overhang75bp
    GTF=$REF/ensembl/Homo_sapiens.GRCh37.70.processed.gtf
    
    # Paired end reads
	STAR --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --genomeDir $GENOME \ 
	   --sjdbGTFfile $GTF --runThreadN $NCPU --outReadsUnmapped Fastx --outStd SAM --outSAMmode Full \
	   --outFileNamePrefix $OUTFILE --readFilesCommand zcat --readFilesIn $FASTQ1 $FASTQ2 > $OUTFILE.sam
	
	# Single end
	# You can use a single FASTQ or 2 if you are combining multiple runs
	STAR --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outSAMstrandField intronMotif --genomeDir $GENOME \
	   --sjdbGTFfile $GTF --runThreadN $NCPU --outReadsUnmapped Fastx --outStd SAM --outSAMmode Full \
	   --outFileNamePrefix $OUTFILE --readFilesCommand cat --readFilesIn <(zcat $FASTQ1 $FASTQ2) > $OUTFILE.sam
	
	# The FASTQ of unmapped reads is called ${OUTFILE}Unmapped.out.mate1

2) Obtain FASTA files of virus sequence from [NCBI](http://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?taxid=10239&opt=Virus).  Or this can be any FASTA file of any sequences.  I have made a public file of all virus genomes from NCBI:

    /sc/orga/projects/ngs/resources/NCBI_virus_db/NCBI_virus_complete_genome.fasta
and I have  created the BLAST database files.    
    
Or create your own:

    # Create a BLAST database from this FASTA file
    module load blast/2.2.26
    makeblastdb -in $SEQDB -dbtype nucl
    
3) Run countVirusReads.sh on FASTQ of unmapped reads.  This script uses Trinity/Inchworm to assemble reads into contigs, BLAST to identify sequences in database matching these contigs, and TopHat to count the number of reads mapping to each contig

    OUTFOLDER=~/test/
    NCPU=3
    
    # Arguments: [FASTQ of unmapped reads] [FASTA file for BLAST database] [NCPU] [folder to save results]
    countVirusReads.sh ${OUTFILE}Unmapped.out.mate1 $SEQDB $NCPU $OUTFOLDER

    # The file $OUTFOLDER/quantify_virus.tsv stores the results of Inchworm, BLAST and TopHat.
    # For each Inchworm contig, it shows the BLAST statistics and the number of reads mapped to the contig by TopHat


