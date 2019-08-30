#######################################################################################################################
#       A script for unzipping, aligning, deduplicating, removal of reads mapping to noncodingRNA for GRO-seq data
#                    (requires lists of annotated non coding RNA- see methods for generation)
#
#                     To run, place in folder with zipped single end R1.fastq.gz file and type:
#
#                                    nohup perl GROseq_pipe.pl <$genome> & 
#
#                          Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################
#!/usr/bin/perl -w
use strict;
my $start_run = time();
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 bowtie/1.1.2 bedtools/2.25.0 ucsctools/1.0));
print "samtools/1.3, bowtie/1.1.2 and bedtools/2.25.0 loaded\n\n";
my $genome = $ARGV[0];

# COMMAND 0; *get $path
my $path; # (non coding RNA)

if ($genome =~ /mm9/)
{
    $path = "/t1-data/user/hugheslab/mlarke/Mouse/nonCodingRNALists/combined.srt.slop100r100.bed";
}
elsif ($genome =~ /hg19/)
{
    $path = "/t1-data/user/hugheslab/mlarke/Human/nonCodingRNALists/combined.srt.slop100r100.bed";
}
else
{
    print "ERROR: genome does not equal mm9 or hg19, exiting script\n";
    exit;
}

# COMMAND 1: generate unzipped copies of master R1 file
system ("nohup zcat R1.fastq.gz > R1.fastq");
print "R1.fastq.gz unzipped\n\n";

# COMMAND 2: align the R1 file
system ("nohup bowtie -m 2 -p 1 --maxins 1000 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/$genome/$genome R1.fastq alignment.sam");
print "R1.fastq aligned to $genome\n\n";

# COMMAND 3; filter properly paired reads
system ("nohup samtools view -bS -o mapped.bam alignment.sam");
print "mapped, properly paired pairs filtered\n\n";

# COMMAND 4; sort by coordinate
system ("nohup samtools sort mapped.bam mapped.srtC");
print "sort by coordinate complete\n\n";

# COMMAND 5; remove duplicates
system ("nohup samtools rmdup -s mapped.srtC.bam filtered.bam");
print "duplicates removed\n\n";

# COMMAND 6; convert bam to bed
system ("nohup bedtools bamtobed -i filtered.bam > filtered.bed");
print "duplicates removed\n\n";

# COMMAND 7; intersect the bam file withe the snRNA file
system ("nohup bedtools intersect -a filtered.bed -b $path -v -s > frags.snRNAremoved.bed");
print "snRNA etc removed from mapped.srtC.bam\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;
