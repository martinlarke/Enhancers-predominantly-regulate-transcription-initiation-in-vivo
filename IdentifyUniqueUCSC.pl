#######################################################################################################################
#  A script load in all UCSC genes from the table browser and concatenate gene ids for genes with identical coordinates
#
#                   To run, place in folder with the ucsc annotated whole genes in .bed6 format
#
#                                nohup perl IdentifyUniqueUCSC.pl <$genes.bed6> & 
#
#                          Coded by Martin Larke, University of Oxford 30-Aug-2019
#######################################################################################################################

#!/usr/bin/perl -w
use strict;
use Data::Dumper;
$Data::Dumper::Sortkeys  = 1;

my $genes= $ARGV[0];
my $filename = $genes;
$filename =~ s/.bed//gi;

open GENES, "$genes" or die "cant open $genes";
my $chr;
my $start;
my $stop;
my $id;
my $strand;
my $count;
my %genes;

while (my $line =<GENES>)
{
    #$count++;
    #if ($count==5){last;}
    
    chomp $line;
    my @columns = split(/\t/, $line);
    $chr    = $columns[0];
    $start  = $columns[1];
    $stop   = $columns[2];
    $id     = $columns[3];
    $strand = $columns[5];
    
    if (exists $genes{$strand}{$chr}{$start}{$stop})
    {
         $genes{$strand}{$chr}{$start}{$stop}.="__$id"; 
    }
    else
    {
        $genes{$strand}{$chr}{$start}{$stop}.=$id;
    }
    
}
#print Dumper (\%genes);
close GENES;

open BED,  '>',"$filename.unique.bed" or die "$filename.unique.bed";
open GFF,  '>',"$filename.unique.gff" or die "$filename.unique.gff";

foreach my $strand (sort keys %genes)
{
    foreach my $chr (sort keys %{$genes{$strand}})
    {
        foreach my $start (sort keys %{$genes{$strand}{$chr}})
        {
            foreach my $stop (sort keys %{$genes{$strand}{$chr}{$start}})
            {
                my $id = $genes{$strand}{$chr}{$start}{$stop};
                print BED "$chr\t$start\t$stop\t$id\t0\t$strand\n";
                print GFF "$chr\tucsc\thg19\tgene\t$start\t$stop\t0\t$strand\t0\t0\n";
            }
        }
    }
}
close BED;
close GFF;

exit;
