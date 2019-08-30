#######################################################################################################################
#       A script to align paired RNA-seq data, extract reads mapping to introns, remove reads mapping to exonic or
#							non coding RNA and to produce stranded bigwigs of intronic reads (filtered) and all reads (allreads.bw)
#
#                   To run, place in folder with zipped paired end R1.fastq.gz and R2.fastq.gz file and type:
#
#                                nohup perl alignment2IntronicBamBgBw.pl <$genome> & 
#
#                              Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w
use strict;
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 bowtie/1.1.2 bedtools/2.25.0 ucsctools/1.0));
print "samtools/0.1.19, bowtie/1.1.2, bedtools/2.25.0 and ucsctools/1.0 loaded\n\n";
my $start_run = time();

my $genome = $ARGV[0];
chomp $genome;

# unzip R1 and R2 files
system ("nohup zcat R1.fastq.gz > R1.fastq");
system ("nohup zcat R2.fastq.gz > R2.fastq");
print "R1 and R2 fastq.gz files unzipped\n\n";

# align the fastq files
system "(nohup bowtie -m 2 -p 4 --maxins 1000 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/$genome/$genome alignment.sam -1 R1.fastq -2 R2.fastq)";
print "reads aligned\n\n";

# filter mapped reads
system "(nohup samtools view -bSh -F 4 -o alignment.bam alignment.sam)";
print "alignment.sam -mapped-> alignment.bam\n\n";

# make a symbolic links to the ncRNA, introns, and exons files
if ($genome =~ /mm9/)
{
	system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Mouse/nonCodingRNALists_mm9/Final/combined.srt.slop100r100.bed)";
    system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Mouse/UCSCgenes/ucsc.exons.mm9.bed)";
    system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Mouse/UCSCgenes/ucsc.introns.mm9.bed)";
    print "ln -s --> .beds made\n\n";
}
elsif ($genome =~ /hg19/)
{
    system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Human/nonCodingRNALists_hg19/ncRNA_Final/combined.srt.slopl100r100.bed)";
    system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Human/UCSC_final/ucsc.exons.hg19.bed)";
    system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Human/UCSC_final/ucsc.introns.hg19.bed)";
    print "ln -s --> .beds made\n\n";
}
else
{
    print "ERROR: genome does not equal mm9 or hg19, exiting script\n";
    exit;
}


# get reads that map 100% to an intron (regardless of overlap between gene isoforms)
system "(nohup fetchChromSizes $genome > $genome.genome)";
print "genome sizes for $genome retrieved\n\n";
system "(nohup bedtools intersect -a alignment.bam -b ucsc.introns.$genome.bed -wa -s -split -f 1 > alignment.intronic.bam)";
print "reads mapping 100% to introns extracted\n\n";

# remove reads mapping to exons and ncRNA (even by 1bp)
system "(nohup bedtools intersect -a alignment.intronic.bam -b ucsc.exons.$genome.bed combined.srt.slopl100r100.bed -v -wa -s -split > alignment.intronic.noexons.ncRNAremoved.bam)";
print "reads mapping by even 1bp to ncRNA or exons removed\n\n";

# sort intronic.noexons.ncRNAremoved.bam by coordinate --> filtered.bam
system "(nohup samtools sort alignment.intronic.noexons.ncRNAremoved.bam filtered)";
print "intronic.noexons.ncRNAremoved.bam -srtC > filtered.bam\n\n";

# get genome file
system ("fetchChromSizes $genome > $genome.genome");
print "$genome genome file retrieved\n\n";

# make coverage files
system "(nohup bedtools genomecov -g $genome.genome -ibam filtered.bam -bg -strand + > filtered.plus.bg)";
system "(nohup bedtools genomecov -g $genome.genome -ibam filtered.bam -bg -strand - > filtered.minus.bg)";
print "filtered.bam --> (+/-) .bgs \n\n";

# sort coverage files
system "(nohup sort -k1,1 -k2,2n filtered.plus.bg > filtered.plus.srt.bg)";
system "(nohup sort -k1,1 -k2,2n filtered.minus.bg > filtered.minus.srt.bg)";
print ".bgs sorted\n\n";

#make bigwigs
system "(nohup bedGraphToBigWig filtered.plus.srt.bg $genome.genome filtered.plus.bw)";
system "(nohup bedGraphToBigWig filtered.minus.srt.bg $genome.genome filtered.minus.bw)";
print "stranded bigWigs from bedGraphs\n\n";

# get millions of mapped reads
system "(nohup samtools flagstat filtered.bam > filtered.flagstat)";
print "flagstats calculated on filtered.bam (intronic (100%), noExons (even 1E7%))\n\n";

# cp hubber files
system "(nohup cp /t1-data/user/hugheslab/mlarke/scripts/Hubber/Hubber4IntronicRNAseq.pl .)";
system "(nohup cp /t1-data/user/hugheslab/mlarke/scripts/Hubber/info.txt .)";
print "hubber files copied to directory\\n";

# zip the alignment.sam file
system ("nohup gzip --best alignment.sam");
print "alignment.sam zipped\n\n";

# cleanup temp files
system ("rm -rf R1.fastq R2.fastq *.sam all.plus.bg all.minus.bg filtered.plus.bg filtered.minus.bg");
print "temp files removed\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;

