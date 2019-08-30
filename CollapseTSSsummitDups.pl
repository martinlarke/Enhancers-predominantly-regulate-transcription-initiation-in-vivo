#!/usr/bin/perl -w
use strict;
use Data::Dumper;
$Data::Dumper::Sortkeys  = 1;

########################################### Coded by Martin Larke 21-July-2018 ########################################
#
#   load in all TSS form the combined TSS summits files  collapse exact duplicates and count number of occurances in file
#
#
#   - Run command:  nohup perl  CollapseTSSsummitDups.pl combinedTSS.summit.bed &
#
#####################################################################################################################

my $tss= $ARGV[0];
my $filename = $tss;
$filename =~ s/.bed//gi;

open TSS, "$tss" or die "cant open $tss";

my $chr;
my $start;
my $id;
my $strand;

#my $count;
my %tss;

while (my $line =<TSS>)
{
    #$count++;
    #if ($count==100){last;}
    
    chomp $line;
    my @columns = split(/\t/, $line);
    $chr    = $columns[0];
    $start  = $columns[1];
    $id     = $columns[3];
    $strand = $columns[5];
    
    #add to hash and collapse duplicates by adding a count
    $tss{$strand}{$chr}{$start}++; 

}
#print Dumper (\%tss);
close TSS;

open BED,  '>',"$filename.nodups.bed" or die "$filename.nodups.bed";

foreach my $strand (sort keys %tss)
{
    foreach my $chr (sort keys %{$tss{$strand}})
    {
        foreach my $start (sort keys %{$tss{$strand}{$chr}})
        {
            my $count = $tss{$strand}{$chr}{$start};
            print BED "$chr\t$start\t$start\t$count\t0\t$strand\n";
        }
    }
}
close BED;

exit;
