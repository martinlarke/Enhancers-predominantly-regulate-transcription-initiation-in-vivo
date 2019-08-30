#######################################################################################################################
#       A script for unzipping, aligning, deduplicating, and converting Mnase read pairs into reconstructed DNA fragments

#                   To run, place in folder with paired end R1.fastq.gz and R2.fastq.gz file and type:
#
#                                nohup perl MNase_pipe.pl <$genome (mm9/hg19)> &
#
#                          Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w
use strict;
my $start_run = time();
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 bowtie/1.1.2 bedtools/2.25.0 ucsctools/1.0));
print "samtools/0.1.19, bowtie/1.1.2 uscsctools/1.0 and and bedtools/2.25.0 loaded\n\n";

my $genome = $ARGV[0];
chomp $genome;

# generate unzipped copies of master R1 file
system ("nohup zcat R1.fastq.gz > R1.fastq");
system ("nohup zcat R2.fastq.gz > R2.fastq");
print "R1.fastq.gz and R2.fastq.gz unzipped\n\n";

# align the R1 and R2 files
system ("nohup bowtie -m 1 -p 4 --maxins 250 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/$genome/$genome alignment.sam -1 R1.fastq -2 R2.fastq");
print "fastqs mapped\n\n";

# filter properly paired reads
system ("nohup samtools view -bS -f 3 -o properpairs.bam alignment.sam");
print "mapped, properly paired pairs filtered\n\n";

# sort by coordinate
system ("nohup samtools sort properpairs.bam properpairs.srtC");
print "sort by coordinate complete\n\n";

# remove duplicates
system ("nohup samtools rmdup properpairs.srtC.bam filtered.bam");
print "duplicates removed\n\n";

system ("nohup samtools sort -n filtered.bam filtered.srtN");
print "sort by coordinate complete\n\n";

# convert into a bedpe file
system ("nohup bedtools bamtobed -bedpe -mate1 -i filtered.srtN.bam | grep chr | cut -f 1,2,3,5,6,9 > TEMP_bedpe_fragments.txt");
print "temporary fragments file made from properpairs.bam\n\n";

# make plus and minus strand fragments
my $awk1= '{if ($6=="+") print $1"\t"$2"\t"$5"\t""0""\t""0""\t""+"}';
my $awk2= '{if ($6=="-") print $1"\t"$4"\t"$3"\t""0""\t""0""\t""-"}';
system ("nohup cat TEMP_bedpe_fragments.txt | awk '$awk1' > plus_frags.bedpe");
system ("nohup cat TEMP_bedpe_fragments.txt | awk '$awk2' > minus_frags.bedpe");
print "temporary fragments file converted into plus and minus strand fragment files\n\n";

# combine stranded frags
system ("nohup cat plus_frags.bedpe minus_frags.bedpe > frags.bedpe");
print "fwd and rev bed6 files combined to frags.bedpe\n\n";

# zip the alignment.sam file
system ("nohup gzip --best alignment.sam");
print "alignment.sam zipped\n\n";

# cleanup temp files
system ("rm -rf *L00* *.fastq *.sam");
print "temp files removed\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;
