#######################################################################################################################
#    Extracts RNA molecules from scaRNA frags file, converts to colour coded density track for UCSC
#
#    To run, place in folder with region to analyse in .bed for mate e.g gene sof interest:
#
# nohup perl frags2bed.pl <frags.bed, e.g. frags.snRNAremoved.bed>  <.bed of regions to analyse>  <minimum_coverage_depth>  <genomebuild> &
#
#                 e.g:  nohup perl frags2bed.pl frags.snRNAremoved.bed bed alphaglobin.bed 5 mm9 &
#
#                             Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

use strict;
use Env::Modulecmd;
Env::Modulecmd::load (qw(ucsctools/1.0 bedtools/2.25.0));

# open inputs and outputs, parse .sam filename
my $frags = $ARGV[0];
my $region = $ARGV[1];
my $min_depth = $ARGV[2];
my $genome = $ARGV[3];

# parse out frags matching region and add to hash (with a count for each uniqe fragment)
system ("bedtools intersect -wa -a $frags -b $region > regionfrags.bed");
open REGIONFRAGS, "regionfrags.bed" or die "Can't open regionfrags.bed";
my %hash;

while (my $line =<REGIONFRAGS>)
{
    chomp $line;
    my @array = split(/\t/, $line);
    my $chr = $array[0];
    my $start = $array[1];
    my $stop = $array[2];
    my $strand = $array[5];
    
    $hash{$chr}{$start}{$stop}{$strand}++;
}

#print output bed file with scores (log_10) for each fragment
open BED, '>',"output.bed" or die "Can't create output.bed";
print BED "track name=\"reconstructed reads.bed\" description=\"reconstructed reads displayed as a bed file, minimum read depth = $min_depth.\" useScore=1 spectrum=on scoreMin=1 scoreMax=25\n";  #use score colours reads according to frequency (set Min and Max to represnet upper and lower limits of reads for a given position in the genome (visualise bed to estrimate max fragment depth
my $log_10_score = 0;
my $rounded = 0;

foreach my $chr (sort keys %hash)
{   
    foreach my $start (sort keys %{$hash{$chr}})
    {
        foreach my $stop (sort keys %{$hash{$chr}{$start}})
        {
            foreach my $strand (sort keys %{$hash{$chr}{$start}{$stop}})
            {
                if ($hash{$chr}{$start}{$stop}{$strand} < $min_depth){next;}     # skip low coverage fragments
                else
                {
                    print BED "$chr\t$start\t$stop\t$hash{$chr}{$start}{$stop}{$strand}\t$hash{$chr}{$start}{$stop}{$strand}\t$strand\n";
                }
            }
        }
    }
}

close REGIONFRAGS;
close BED;
exit;
