#######################################################################################################################
#       A script for unzipping, aligning, deduplicating and converting fastq files into a bigwig and bedgraph
#
#                   To run, place in folder with zipped single end R1.fastq.gz file and type:
#
#                                nohup perl ChIP_pipe_PEproperpairs.pl <$genome> & 
#
#                             Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################
#!/usr/bin/perl -w
use strict;

#load modules required for script
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 bowtie/1.1.2 ucsctools/1.0 bedtools/2.25.0));
print "samtools/0.1.19, bowtie/1.1.2, ucsctools/1.0 and bedtools/2.25.0 loaded\n\n";

#input the genome required (e.g. mm9 or hg19)
my $genome = $ARGV[0];

# concatenate lane fastq.gz's into R1 and R2 files
system ("nohup cat *L001_R1* *L002_R1* *L003_R1* *L004_R1* > R1.fastq.gz");
system ("nohup cat *L001_R2* *L002_R2* *L003_R2* *L004_R2* > R2.fastq.gz");
print "R1 and R2 lane fastq.gz files combined to R1.fastq.gz and R2.fastq.gz respectively\n\n";

# unzip R1 and R2 files
system ("nohup zcat R1.fastq.gz > R1.fastq");
system ("nohup zcat R2.fastq.gz > R2.fastq");
print "R1 and R2 fastq.gz files unzipped\n\n";

# align R1 and R2 files
system ("nohup bowtie -m 2 -p 4 --maxins 1000 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/$genome/$genome alignment.sam -1 R1.fastq -2 R2.fastq");
print "alignment complete\n\n";

# unzip sam
#system ("nohup zcat alignment.sam.gz > alignment.sam");
#print "alignment file unzipped\n\n";

# filter properly paired reads
system ("nohup samtools view -bS -f 3 -o properpairs.bam alignment.sam");
print "mapped, properly paired pairs filtered\n\n";

# count mapped proper pairs
system ("nohup samtools flagstat properpairs.bam > properpairs.flagstat");
print "properly paired mapped reads counted\n\n";

# sort by coordinate
system ("nohup samtools sort properpairs.bam properpairs.srtC");
print "sort by coordinate complete\n\n";

# remove PCR duplicates
system ("nohup samtools rmdup properpairs.srtC.bam filteredpairs.bam");
print "duplicates removed\n\n";

# COMMAND 5; sort the merged bam file by name
system ("nohup samtools sort -n filteredpairs.bam filteredpairs.srtN");
print "sort by coordinate complete\n\n";

# COMMAND 6; convert properpairs.bam into a bedpe file
system ("nohup bedtools bamtobed -bedpe -mate1 -i filteredpairs.srtN.bam | grep chr | cut -f 1,2,3,5,6,9 > TEMP_bedpe_fragments.txt");
print "temporary fragments file made from properpairs.bam\n\n";

# COMMAND 7; make frags.bed
my $awk1= '{if ($6=="+") {print $1"\t"$2"\t"$5"\t""0""\t""0""\t""+" > "frags.bed"} else if ($6=="-") {print $1"\t"$4"\t"$3"\t""0""\t""0""\t""-" > "frags.bed"}}';
system ("nohup awk '$awk1' TEMP_bedpe_fragments.txt");
print "temporary fragments file converted into frags.bed\n\n";

# COMMAND 9; sort frags file
system ("nohup sort -k1,1 -k2,2n frags.bed > frags.srt.bed");
print "fwd and rev bed6 files combined\n\n";

# get genome file
system ("fetchChromSizes $genome > $genome.genome");
print "$genome genome file retrieved\n\n";

# make bedraph
system ("nohup bedtools genomecov -i frags.srt.bed -g $genome.genome -bg > fullfrag.bg");
print "bedgraph sorted\n\n";

# sort bedgraph
system ("nohup sort -k1,1 -k2,2n fullfrag.bg > fullfrag.srt.bg");
print "bedgraph made\n\n";

# make bigWig
system ("nohup bedGraphToBigWig fullfrag.srt.bg $genome.genome fullfrag.bw");
print "bigWig made\n\n";

# cleanup temp files
system ("rm -rf *L00* *.fastq *.sam *TEMP* properpairs.bam filteredpairs.bam frags.bed fullfrag.bg");
print "temp files removed\n\n";

exit;
