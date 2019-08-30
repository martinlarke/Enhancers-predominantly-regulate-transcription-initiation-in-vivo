#######################################################################################################################
#       A script to align paired RNA-seq data, remove reads mapping non coding RNA and to produce stranded bigwigs
#
#                   To run, place in folder with zipped paired end R1.fastq.gz R2.fastq.gz file and type:
#
#                                nohup perl alignment2allreads.pl <$genome> & 
#
#                           Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w
use strict;
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 rna-star/2.4.2a bedtools/2.25.0 ucsctools/1.0));
print "samtools/0.1.19, rna-star/2.4.2a, bedtools/2.25.0 and ucsctools/1.0 loaded\n\n";
my $start_run = time();

my $genome = $ARGV[0];
chomp $genome;

# align the R1 and R2 files
system ("nohup STAR --genomeDir /databank/indices/star/$genome --readFilesCommand zcat --readFilesIn R1.fastq.gz R2.fastq.gz --outFileNamePrefix 1_ --outFilterMultimapNmax 2");
print "reads aligned to $genome\n\n";

# filter mapped reads
system "(nohup samtools view -bSh -F 4 -o alignment.bam 1_Aligned.out.sam)";
print "alignment.sam -mapped-> alignment.bam\n\n";

# make a symbolic links to the ncRNA
if ($genome =~ /hg19/)
{
    system "(nohup ln -s /t1-data/user/hugheslab/mlarke/Human/nonCodingRNALists_hg19/ncRNA_Final/combined.srt.slopl100r100.bed)";
    print "ln -s --> .beds made\n\n";
}
else
{
    print "ERROR: genome does not equal hg19, exiting script\n";
    exit;
}

# remove reads mapping ncRNA (even by 1bp)
system "(nohup bedtools intersect -a alignment.bam -b combined.srt.slopl100r100.bed -v -wa -s -split > alignment.ncRNAremoved.bam)";
print "reads mapping by even 1bp to ncRNA removed\n\n";

# sort by coordinate
system ("nohup samtools sort alignment.ncRNAremoved.bam alignment.srtC");
print "sort by coordinate complete\n\n";

# get genome file
system ("fetchChromSizes $genome > $genome.genome");
print "$genome genome file retrieved\n\n";

# make coverage files
system "(nohup bedtools genomecov -g $genome.genome -ibam alignment.srtC.bam -bg -strand + > all.plus.bg)";
system "(nohup bedtools genomecov -g $genome.genome -ibam alignment.srtC.bam -bg -strand - > all.minus.bg)";
print "filtered.bam --> (+/-) .bgs \n\n";

# sort coverage files
system "(nohup sort -k1,1 -k2,2n all.plus.bg > all.plus.srt.bg)";
system "(nohup sort -k1,1 -k2,2n all.minus.bg > all.minus.srt.bg)";
print ".bgs sorted\n\n";

# make bigwigs
system "(nohup bedGraphToBigWig all.plus.srt.bg $genome.genome all.plus.bw)";
system "(nohup bedGraphToBigWig all.minus.srt.bg $genome.genome all.minus.bw)";
print "(+/-) .bgs --> (+/-) .bws \n\n";

# get millions of mapped reads
system "(nohup samtools flagstat alignment.srtC.bam > alignment.flagstat)";
print "flagstats calculated on filtered.bam (intronic (100%), noExons (even 1E7%))\n\n";

# cp hubber files
system "(nohup cp /t1-data/user/hugheslab/mlarke/scripts/Hubber/Hubber4IntronicRNAseq.pl .)";
system "(nohup cp /t1-data/user/hugheslab/mlarke/scripts/Hubber/info.txt .)";
print "hubber files copied to directory\\n";

# zip the alignment.sam file
system ("nohup gzip --best 1_Aligned.out.sam");
print "alignment.sam zipped\n\n";

# cleanup temp files
system ("rm -rf R1.fastq R2.fastq *.sam all.plus.bg all.minus.bg");
print "temp files removed\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;

