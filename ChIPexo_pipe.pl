
#######################################################################################################################
#  A script for unzipping, aligning, deduplicating and converting fastq files into a bam fie for further proscessing
#
#                   To run, place in folder with zipped paired end R1.fastq.gz and R2.fastq.gz file and type:
#
#                                           nohup perl ChIPexo_pipe.pl <$genome> & 
#
#                                   Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w

use strict;
my $start_run = time();
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 bowtie/1.1.2));
print "samtools/0.1.19, deeptools/2.2.2\n\n";

# align the fastq file
system ("nohup bowtie -m 2 -p 4 --maxins 1000 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/hg19/hg19 R1.fastq > alignment.sam");

# filter mapped reads
system ("nohup samtools view -bS -F 4 -o mapped.bam alignment.sam");
print "mapped, properly paired pairs filtered\n\n";

# sort by coordinate
system ("nohup samtools sort mapped.bam mapped.srtC");
print "sort by coordinate complete\n\n";

# remove duplicates
system ("nohup samtools rmdup -s mapped.srtC.bam filtered.bam");
print "duplicates removed\n\n";

# zip the alignment.sam file
system ("nohup gzip --best alignment.sam");
print "alignment.sam zipped\n\n";

# cleanup temp files
system ("rm -rf alignment.bam");
print "temp files removed\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;
