#######################################################################################################################
#               A script for analysing mNETseq, which can be used to produce scaled bigWigs for visualisation
#       Analysis stages: unzipping, alignment, PCR duplicate removal, make stranded bw of fullfrag and 2bp windowed tes
#
#               To run, place in folder with zipped paired end R1.fastq.gz and R2.fastq.gz file and type:
#
#                   nohup perl mNETseq_pipe_withDupXYMandSSremoval.pl <$genome> <numeric scaling factor> & 
#
#                               Coded by Martin Larke, University of Oxford 12-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w
use strict;
my $start_run = time();
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/1.3 bowtie/1.1.2 bedtools/2.25.0 ucsctools/1.0));
print "samtools/1.3, bowtie/1.1.2 uscsctools/1.0 and and bedtools/2.25.0 loaded\n\n";

my $genome = $ARGV[0];
my $scalingfactor = $ARGV[1];
chomp $genome;

# get $path
my $path2; # (splice sites)

if ($genome =~ /mm9/)
{
    $path2 = "/t1-data/user/hugheslab/mlarke/Mouse/UCSCgenes/SpliceSites/ucsc.mm9.codingexons.both.3primeSS.bed";
}
elsif ($genome =~ /hg19/)
{
    $path2 = "/t1-data/user/hugheslab/mlarke/Human/UCSC_final/SpliceSites/ucsc.hg19.codingexons.both.3primeSS.bed";
}
else
{
    print "ERROR: genome does not equal mm9 or hg19, exiting script\n";
    exit;
}

my $hubber1 = "/t1-data/user/hugheslab/mlarke/scripts/Hubber/Hubber4mNETseq.pl";
my $info1   = "/t1-data/user/hugheslab/mlarke/scripts/Hubber/info.txt";
my $hubber2 = "/t1-data/user/hugheslab/mlarke/scripts/Hubber/Hubber4mNETseqFullfrag.pl";

# concatenate lane fastq.gz files to master R1 and R1 files
#system ("nohup cat *L001_R1* *L002_R1* *L003_R1* *L004_R1* > R1.fastq.gz");
#system ("nohup cat *L001_R2* *L002_R2* *L003_R2* *L004_R2* > R2.fastq.gz");
#print "R1 and R2 lane fastq.gz files combined to R1.fastq.gz and R2.fastq.gz respectively\n\n";

# generate unzipped copies of master R1 and R2 fileS
system ("nohup zcat R1.fastq.gz > R1.fastq");
system ("nohup zcat R2.fastq.gz > R2.fastq");
print "R1.fastq.gz and R2.fastq.gz unzipped\n\n";

# align the R1 and R2 files
system ("nohup bowtie -m 2 -p 3 --maxins 1000 --best --strata --chunkmb 256 --sam /databank/indices/bowtie/$genome/$genome alignment.sam -1 R1.fastq -2 R2.fastq");
print "alignment complete\n\n";

# filter properly paired reads
system ("nohup samtools view -bS -f 3 -o properpairs.bam alignment.sam");
print "mapped, properly paired pairs filtered\n\n";

# sort by coordinate
system ("nohup samtools sort -T temp -o properpairs.srtC.bam properpairs.bam");
print "sort by coordinate complete\n\n";

# remove PCR duplicates
system ("nohup samtools rmdup properpairs.srtC.bam filteredpairs.bam");
print "duplicates removed\n\n";

# intersect the filteredpairs.bam file with the splice sites file
system ("nohup bedtools intersect -abam filteredpairs.bam -b $path2 -v -s > filtered.ssRemoved.bam");
print "reads overlapping splice sites removed from filterd.bam \n\n";

# convert into sam
system ("nohup samtools view -h -@ 8 -o filtered.ssRemoved.sam filtered.ssRemoved.bam");
print "file converted into .sam format\n\n";

# remove x, y and M chr peaks
my $grep = "chr[a-zA-Z]";
system ("nohup grep -v '$grep' filtered.ssRemoved.sam > filtered.ssRemoved.noXYM.sam");
print "filtered.ssRemoved.sam parsed to remove non numberical and random chromosomes\n\n";

# filter properly paired reads and convert back into to bam
system ("nohup samtools view -bS -f 3 -o filtered.ssRemoved.noXYM.bam filtered.ssRemoved.noXYM.sam");
print "mapped, properly paired pairs filtered\n\n";

# sort the merged bam file by name
system ("nohup samtools sort -n -T temp -o filtered.ssRemoved.noXYM.srtN.bam filtered.ssRemoved.noXYM.bam");
print "sort by coordinate complete\n\n";

# get aligned read stats
system ("nohup samtools flagstat filtered.ssRemoved.noXYM.srtN.bam > flagstat.filtered.ssRemoved.noXYM.srtN.txt");
print "flagstat analysis complete\n\n";

# convert filtered.ssRemoved.noXYM.srtN.bam into a bedpe file
system ("nohup bedtools bamtobed -bedpe -mate1 -i filtered.ssRemoved.noXYM.srtN.bam | grep chr | cut -f 1,2,3,5,6,9 > TEMP_bedpe_fragments.txt");
print "temporary fragments file made from filteredpairs.srtN.bam\n\n";

# make plus and minus strand fragments
my $awk1= '{if ($6=="+") print $1"\t"$2"\t"$5"\t""0""\t""0""\t""+"}';
my $awk2= '{if ($6=="-") print $1"\t"$4"\t"$3"\t""0""\t""0""\t""-"}';
system ("nohup cat TEMP_bedpe_fragments.txt | awk '$awk1' > plus_frags.bedpe");
system ("nohup cat TEMP_bedpe_fragments.txt | awk '$awk2' > minus_frags.bedpe");
print "temporary fragments file converted into plus and minus strand fragment files\n\n";

# combine stranded frags
system ("nohup cat plus_frags.bedpe minus_frags.bedpe > frags.bedpe");
print "fwd and rev bed6 files combined to frags.bedpe\n\n";

# sort frags file
system ("nohup sort -k1,1 -k2,2n frags.bedpe > frags.srt.bedpe");
print "fwd and rev bed6 files combined\n\n";

# get genome file
system ("fetchChromSizes $genome > $genome.genome");
print "$genome genome file retrieved\n\n";

## make tes file
my $awk3='{if ($6~/\+/) {print $1"\t"$3"\t"$3"\t""0""\t""0""\t"$6 > "frags.tes.bed"} else if ($6~/\-/) {print $1"\t"$2"\t"$2"\t""0""\t""0""\t"$6 > "frags.tes.bed"}}';
system ("nohup awk '$awk3' frags.srt.bedpe");
print "transcription end site (tes) file made\n\n";

## window tes file
system ("nohup bedtools slop -i frags.tes.bed -b 2 -g $genome.genome > frags.tes.windowed.bed");
print "transcription end sites windowed\n\n";

# make stranded tes bedgraphs
system ("nohup mkdir r4_analysis");
system ("nohup cp $hubber1 r4_analysis/");
system ("nohup cp $info1 r4_analysis/");
system ("nohup bedtools genomecov -i frags.tes.windowed.bed -g $genome.genome -bg -strand + -scale $scalingfactor > ./r4_analysis/fwd.tes.bg");
system ("nohup bedtools genomecov -i frags.tes.windowed.bed -g $genome.genome -bg -strand - -scale $scalingfactor > ./r4_analysis/rev.tes.bg");
print "stranded bg made for full frag and 4bp windowed tes\n\n";

# make fullfrag bedgraphs
system ("nohup mkdir fullfrag_analysis");
system ("nohup cp $hubber2 fullfrag_analysis/");
system ("nohup cp $info1 fullfrag_analysis/");
system ("nohup bedtools genomecov -i frags.srt.bedpe -g $genome.genome -bg -strand + -scale $scalingfactor > ./fullfrag_analysis/fwd.fullfrag.bg");
system ("nohup bedtools genomecov -i frags.srt.bedpe -g $genome.genome -bg -strand - -scale $scalingfactor > ./fullfrag_analysis/rev.fullfrag.bg");
print "stranded fullfrag bgs made\n\n";

# convert to bigwigs
system ("nohup bedGraphToBigWig ./r4_analysis/fwd.tes.bg $genome.genome ./r4_analysis/fwd.tes.bw");
system ("nohup bedGraphToBigWig ./r4_analysis/rev.tes.bg $genome.genome ./r4_analysis/rev.tes.bw");
system ("nohup bedGraphToBigWig ./fullfrag_analysis/fwd.fullfrag.bg $genome.genome ./fullfrag_analysis/fwd.fullfrag.bw");
system ("nohup bedGraphToBigWig ./fullfrag_analysis/rev.fullfrag.bg $genome.genome ./fullfrag_analysis/rev.fullfrag.bw");
print "bgs converted into bigwigs\n\n";

# cleanup temp files
system ("rm -rf *L00* *.fastq properpairs.bam filtered.ssRemoved.sam filtered.ssRemoved.noXYM.bam TEMP_bedpe_fragments.txt plus_frags.bedpe minus_frags.bedpe frags.bedpe");
print "temp files removed\n\n";

# zip the alignment.sam file
system ("nohup gzip --best alignment.sam *.bg *.bed *.bg");
print "alignment.sam zipped\n\n";

my $end_run = time();
my $run_time = ($end_run - $start_run);
my $minuteruntime = ($run_time/60);
if ($minuteruntime < 1)	{print "Job took $run_time second(s)\n";}
else			{print "Job took $minuteruntime minute(s)\n";}

exit;
