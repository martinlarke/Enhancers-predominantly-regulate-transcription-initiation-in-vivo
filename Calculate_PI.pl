
#######################################################################################################################
#  A script to calculate pausing indices using a promoter window of -30-->+300 relative to TSS and gene body window as
#  the rest of the annotated gene. Calculated denstiy of Pol2 reads (bam) or frags (bed) over windows and divides promoter
#  by gene body to generate pausing indices.
#
#  To run, place in folder with .bam or .bed of reads/frags and a .bed (whole gene) format download of the ucsc genes
#  from the table browser:
#
#                                nohup perl Calculate_PI.pl <$genes.bed6> <reads.bam/frags.bed> & 
#
#                                      Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Env::Modulecmd;
Env::Modulecmd::load (qw(samtools/0.1.19 bedtools/2.25.0));
print "samtools/0.1.19 and bedtools/2.25.0 loaded\n\n";

my $genes = $ARGV[0]; # list of annotated ucsc genes in whole bed format (from table browser)
my $frags = $ARGV[1]; # bam or bed file of Pol2 ChIP-seq reads or fragments (PCR duplicates removed)

chomp $genes;
chomp $frags;

# make sure genes are at least 301bp long!
my $awk= '{if (($3-$2)>601) print $0}';
system ("nohup awk '$awk' $genes > genes.over301.bed");
print "genes < 300bp parsed out\n\n";

#make windows -30 +300 (prom), +301 --> end of gene (body)
my $awk1= '{if ($6=="+") {print $1"\t"($2-30)"\t"($2+300)"\t"$4"\t"$5"\t"$6 > "prom.bed"} else if ($6=="-") {print $1"\t"($3-300)"\t"($3+30)"\t"$4"\t"$5"\t"$6 > "prom.bed"}}';
my $awk2= '{if ($6=="+") {print $1"\t"($2+301)"\t"$3"\t"$4"\t"$5"\t"$6 > "body.bed"} else if ($6=="-") {print $1"\t"$2"\t"($3-301)"\t"$4"\t"$5"\t"$6 > "body.bed"}}';
system ("nohup awk '$awk1' genes.over301.bed");
system ("nohup awk '$awk2' genes.over301.bed");
print "promoter and genebody windows made\n\n";

#sort promoter and gene body windowed files
system ("nohup sort -k1,1 -k2,2n prom.bed > prom.srt.bed");
system ("nohup sort -k1,1 -k2,2n body.bed > body.srt.bed");
print "windowed files sorted\n\n";

# calculate coverage over windows
system ("nohup bedtools coverage -a prom.srt.bed -b $frags -counts > prom.count");
system ("nohup bedtools coverage -a body.srt.bed -b $frags -counts > body.count");
print "read coverage over windows calculated\n\n";

# combine the coverage files
open PROMCOV, "prom.count" or die "cant open prom.count";
my %prom_cov;

while (my $line1 =<PROMCOV>)
{
    chomp $line1;
    my ($chr, $start, $stop, $name, $cigar, $strand, $cov)=split(/\t/, $line1);
    $prom_cov{$name}= $chr;
    $prom_cov{$name} .= "&$start";
    $prom_cov{$name} .= "&$stop";
    $prom_cov{$name} .= "&$name";
    $prom_cov{$name} .= "&$cigar";
    $prom_cov{$name} .= "&$strand";
    $prom_cov{$name} .= "&$cov";
 
}
close PROMCOV;

# open an output file for the pausing indices and print the header
open BED, '>', "PIs.txt" or die "PIs.txt";
print BED "prom_chr\tprom_start\tprom_stop\tprom_strand\tprom_name\tprom_cov\tbody_name\tbody_start\tbody_stop\tbody_strand\tbody_cov\n";

# open gene body coverage files and pair to promoter coverage based on IDS
open GENECOV, "body.count" or die "cant open body.count";
my %bodycov;
while (my $line2 =<GENECOV>)
{
    chomp $line2;
    my ($body_chr, $body_start, $body_stop, $body_name, $body_cigar, $body_strand, $body_cov)=split(/\t/, $line2);
    
    if (exists $prom_cov{$body_name})
    {
        my $storedvalue = $prom_cov{$body_name};
        my ($prom_chr, $prom_start, $prom_stop, $prom_name, $prom_cigar, $prom_strand, $prom_cov)=split(/&/, $storedvalue);
        print BED "$prom_chr\t$prom_start\t$prom_stop\t$prom_strand\t$prom_name\t$prom_cov\t$body_name\t$body_start\t$body_stop\t$body_strand\t$body_cov\n";
        
        if (($prom_cov > 0) && ($body_cov > 0))
        {
            #calculate pi's
            my $prom_density = ($prom_cov/330);
            my $body_size = $body_stop-$body_start;
            my $body_density = ($body_cov/$body_size);
            my $pi = ($prom_density/$body_density);
            
            # print pi's
            print BED "$prom_chr\t$prom_start\t$prom_stop\t$prom_name\t$pi\n";
        }
        else{next;}
    }
    else{next;}
}
close BED;
exit;
