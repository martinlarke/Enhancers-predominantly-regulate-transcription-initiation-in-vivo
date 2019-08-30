#!/usr/bin/perl -w

#######################################################################################################################
#       A script for unzipping, aligning, deduplicating and converting fastq files into a strand specific bigwig
#                                           of 5' and 3' RNA molecule ends
#
#                   To run, place in folder with zipped paired end R1.fastq.gz and R2.fastq.gz file and type:
#
#                                nohup perl scaRNAseq_pipe.pl <$genome> & 
#
#                          Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################
use strict;
my $start_run = time();
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/1.3 bowtie/1.1.2 bedtools/2.25.0 ucsctools/1.0));
print "samtools/1.3, bowtie/1.1.2 uscsctools/1.0 and and bedtools/2.25.0 loaded\n\n";

my $genome = $ARGV[0];
chomp $genome;

# COMMAND 0; *get $path
my $path; # (non coding RNA)
my $path2; # (splice sites)

if ($genome =~ /mm9/)
{
    $path = "/t1-data/user/hugheslab/mlarke/Mouse/nonCodingRNALists_mm9/Final/combined.srt.slop100r100.bed";
    $path2 = "/t1-data/user/hugheslab/mlarke/Mouse/UCSCgenes/SpliceSites/ucsc.mm9.codingexons.both.3primeSS.bed";
}
elsif ($genome =~ /hg19/)
{
    $path = "/t1-data/user/hugheslab/mlarke/Human/nonCodingRNALists_hg19/combined.slop100.srt.bed6";
    $path2 = "/t1-data/user/hugheslab/mlarke/Human/UCSC_final/SpliceSites/ucsc.hg19.codingexons.both.3primeSS.bed";
}
else
{
    print "ERROR: genome does not equal mm9 or hg19, exiting script\n";
    exit;
}

my $hubber1 = "/t1-data/user/hugheslab/mlarke/scripts/Hubber/Hubber4scaRNAseq.pl";
my $hubber2 = "/t1-data/user/hugheslab/mlarke/scripts/Hubber/Hubber4FwdRevBws.pl";
my $info1   = "/t1-data/user/hugheslab/mlarke/scripts/Hubber/info.txt";

# COMMAND 1: concatenate lane fastq.gz files to master R1 and R1 files
system ("nohup cat *L001_R1* *L002_R1* *L003_R1* *L004_R1* > R1.fastq.gz");
system ("nohup cat *L001_R2* *L002_R2* *L003_R2* *L004_R2* > R2.fastq.gz");
print "R1 and R2 lane fastq.gz files combined to R1.fastq.gz and R2.fastq.gz respectively\n\n";

# COMMAND 2: generate unzipped copies of master R1 file
system ("nohup zcat R1.fastq.gz > R1.fastq");
system ("nohup zcat R2.fastq.gz > R2.fastq");
print "R1.fastq.gz and R2.fastq.gz unzipped\n\n";

# COMMAND 3: align the R1 and R2 files
system ("nohup bowtie -m 2 -p 3 --maxins 1000 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/$genome/$genome alignment.sam -1 R1.fastq -2 R2.fastq");

# COMMAND 4; filter properly paired reads
system ("nohup samtools view -bS -f 3 -o properpairs.bam alignment.sam");
print "mapped, properly paired pairs filtered\n\n";

# COMMAND 5; sort the merged bam file by name
system ("nohup samtools sort -n -T temp -o  filtered.srtN.bam properpairs.bam");
print "sort by coordinate complete\n\n";

# COMMAND 6; convert properpairs.bam into a bedpe file
system ("nohup bedtools bamtobed -bedpe -mate1 -i filtered.srtN.bam | grep chr | cut -f 1,2,3,5,6,9 > TEMP_bedpe_fragments.txt");
print "temporary fragments file made from properpairs.bam\n\n";

# COMMAND 7; make plus and minus strand fragments
my $awk1= '{if ($6=="+") print $1"\t"$2"\t"$5"\t""0""\t""0""\t""+"}';
my $awk2= '{if ($6=="-") print $1"\t"$4"\t"$3"\t""0""\t""0""\t""-"}';
system ("nohup cat TEMP_bedpe_fragments.txt | awk '$awk1' > plus_frags.bedpe");
system ("nohup cat TEMP_bedpe_fragments.txt | awk '$awk2' > minus_frags.bedpe");
print "temporary fragments file converted into plus and minus strand fragment files\n\n";

# COMMAND 8; combine stranded frags
system ("nohup cat plus_frags.bedpe minus_frags.bedpe > frags.bedpe");
print "fwd and rev bed6 files combined to frags.bedpe\n\n";

# COMMAND 9; sort frags file
system ("nohup sort -k1,1 -k2,2n frags.bedpe > frags.srt.bedpe");
print "fwd and rev bed6 files combined\n\n";

# COMMAND 10; intersect the sorted frags file with the snRNA file
system ("nohup bedtools intersect -a frags.srt.bedpe -b $path -v -s > frags.snRNAremoved.bed");
print "snRNA etc removed from frags.srt.bedpe\n\n";

# COMMAND 11; intersect the sorted frags file with the splice sites file
system ("nohup bedtools intersect -a frags.snRNAremoved.bed -b $path2 -v -s > frags.snRNAremoved.ssRemoved.bed");
print "reads overlapping splice sites removed from frags.snRNAremoved.bed\n\n";

# COMMAND 12: get genome file
system ("fetchChromSizes $genome > $genome.genome");
print "$genome genome file retrieved\n\n";

#### VISUALISATION

# 5' and 3' most bases
system ("nohup mkdir standard_analysis");
system ("nohup bedtools genomecov -i frags.snRNAremoved.ssRemoved.bed -g $genome.genome -strand + -5 -bg > ./standard_analysis/fwd.tss.bg");
system ("nohup bedtools genomecov -i frags.snRNAremoved.ssRemoved.bed -g $genome.genome -strand + -3 -bg > ./standard_analysis/fwd.tes.bg");
system ("nohup bedtools genomecov -i frags.snRNAremoved.ssRemoved.bed -g $genome.genome -strand - -5 -bg > ./standard_analysis/rev.tss.bg");
system ("nohup bedtools genomecov -i frags.snRNAremoved.ssRemoved.bed -g $genome.genome -strand - -3 -bg > ./standard_analysis/rev.tes.bg");
system ("nohup bedGraphToBigWig ./standard_analysis/fwd.tss.bg $genome.genome ./standard_analysis/fwd.tss.bw");
system ("nohup bedGraphToBigWig ./standard_analysis/fwd.tes.bg $genome.genome ./standard_analysis/fwd.tes.bw");
system ("nohup bedGraphToBigWig ./standard_analysis/rev.tss.bg $genome.genome ./standard_analysis/rev.tss.bw");
system ("nohup bedGraphToBigWig ./standard_analysis/rev.tes.bg $genome.genome ./standard_analysis/rev.tes.bw");
print "stranded bedgraphs of 5' and 3' ends converted to bigWigs\n\n";
system ("nohup cp $hubber1 $info1 standard_analysis/");
print "stranded bdg and bw made for 5' and 3' ends\n\n";

# 5' and 3' bases with 4bp window
system ("nohup mkdir r4_analysis");
system ("nohup cp $hubber1 $info1  r4_analysis/");
my $awk5='{if (($6=="+") && ($2 > 0)) print $1"\t"($2-2)"\t"($2+2)"\t"$4"\t"$5"\t"$6}';
my $awk6='{if (($6=="+") && ($2 > 0)) print $1"\t"($3-2)"\t"($3+2)"\t"$4"\t"$5"\t"$6}';
my $awk7='{if (($6=="-") && ($2 > 0)) print $1"\t"($3-2)"\t"($3+2)"\t"$4"\t"$5"\t"$6}';
my $awk8='{if (($6=="-") && ($2 > 0)) print $1"\t"($2-2)"\t"($2+2)"\t"$4"\t"$5"\t"$6}';
system ("nohup cat frags.snRNAremoved.ssRemoved.bed | awk '$awk5' > ./r4_analysis/plus_frags.snRNAremoved.5prime.r4.bedpe");
system ("nohup cat frags.snRNAremoved.ssRemoved.bed | awk '$awk6' > ./r4_analysis/plus_frags.snRNAremoved.3prime.r4.bedpe");
system ("nohup cat frags.snRNAremoved.ssRemoved.bed | awk '$awk7' > ./r4_analysis/minus_frags.snRNAremoved.5prime.r4.bedpe");
system ("nohup cat frags.snRNAremoved.ssRemoved.bed | awk '$awk8' > ./r4_analysis/minus_frags.snRNAremoved.3prime.r4.bedpe");
system ("nohup cat ./r4_analysis/plus_frags.snRNAremoved.3prime.r4.bedpe ./r4_analysis/minus_frags.snRNAremoved.3prime.r4.bedpe > ./r4_analysis/both_frags.snRNAremoved.SSremoved.3prime.r4.bedpe");
system ("nohup cat ./r4_analysis/plus_frags.snRNAremoved.5prime.r4.bedpe ./r4_analysis/minus_frags.snRNAremoved.5prime.r4.bedpe > ./r4_analysis/both_frags.snRNAremoved.SSremoved.5prime.r4.bedpe");
system ("nohup sort -k1,1 -k2,2n ./r4_analysis/both_frags.snRNAremoved.SSremoved.5prime.r4.bedpe > ./r4_analysis/both_frags.snRNAremoved.SSremoved.5prime.r4.srt.bedpe");
system ("nohup sort -k1,1 -k2,2n ./r4_analysis/both_frags.snRNAremoved.SSremoved.3prime.r4.bedpe > ./r4_analysis/both_frags.snRNAremoved.SSremoved.3prime.r4.srt.bedpe");
system ("nohup bedtools genomecov -i ./r4_analysis/both_frags.snRNAremoved.SSremoved.3prime.r4.srt.bedpe -g $genome.genome -strand + -bg > ./r4_analysis/fwd.tes.bg");
system ("nohup bedtools genomecov -i ./r4_analysis/both_frags.snRNAremoved.SSremoved.3prime.r4.srt.bedpe -g $genome.genome -strand - -bg > ./r4_analysis/rev.tes.bg");
system ("nohup bedtools genomecov -i ./r4_analysis/both_frags.snRNAremoved.SSremoved.5prime.r4.srt.bedpe -g $genome.genome -strand + -bg > ./r4_analysis/fwd.tss.bg");
system ("nohup bedtools genomecov -i ./r4_analysis/both_frags.snRNAremoved.SSremoved.5prime.r4.srt.bedpe -g $genome.genome -strand - -bg > ./r4_analysis/rev.tss.bg");
system ("nohup bedGraphToBigWig ./r4_analysis/fwd.tes.bg $genome.genome ./r4_analysis/fwd.tes.bw");
system ("nohup bedGraphToBigWig ./r4_analysis/rev.tes.bg $genome.genome ./r4_analysis/rev.tes.bw");
system ("nohup bedGraphToBigWig ./r4_analysis/fwd.tss.bg $genome.genome ./r4_analysis/fwd.tss.bw");
system ("nohup bedGraphToBigWig ./r4_analysis/rev.tss.bg $genome.genome ./r4_analysis/rev.tss.bw");
print "stranded bdg and bw made for 4bp windowed 5' and 3' ends\n\n";

# full fragments
system ("nohup mkdir fullfrag_analysis");
system ("nohup bedtools genomecov -i frags.snRNAremoved.ssRemoved.bed -g $genome.genome -strand + -bg > ./fullfrag_analysis/fwd.bg");
system ("nohup bedtools genomecov -i frags.snRNAremoved.ssRemoved.bed -g $genome.genome -strand - -bg > ./fullfrag_analysis/rev.bg");
system ("nohup bedGraphToBigWig ./fullfrag_analysis/fwd.bg $genome.genome ./fullfrag_analysis/fwd.fullfrag.bw");
system ("nohup bedGraphToBigWig ./fullfrag_analysis/rev.bg $genome.genome ./fullfrag_analysis/rev.fullfrag.bw");
system ("nohup cp $hubber2 $info1  fullfrag_analysis/");
print "stranded bdg and bw made for fwd and rev frags\n\n";

# COMMAND 15: zip the alignment.sam file
system ("nohup gzip --best alignment.sam");
print "alignment.sam zipped\n\n";

# COMMAND 16: cleanup temp files
system ("rm -rf *L00* *.fastq *.sam TEMP_bedpe_fragments.txt frags.bedpe frags.srt.bedpe");
print "temp files removed\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;
