#!/usr/bin/perl -w
use strict;
#use Data::Dumper;
#$Data::Dumper::Sortkeys  = 1;

### run command
#
#   nohup perl PlotObsVsAnn.pl distancetoTSS.txt &
#
#   N.B distancetoTSS.txt must be present (make as follows):
#
#   nohup cut -f 13 closestUCSCtss.500.txt | sort | uniq -c | awk -F" " '{print $1"\t"$2}' > distancetoTSS.txt

my $distancetoTSS = $ARGV[0];

### make a hash from -500 to + 500
my %roll;
my $rollstart   = -500;
my $rollstop    = +500;
for (my $i = $rollstart; $i<=$rollstop; +1)
{
    $roll{$i++}=0;
}
#print Dumper (\%roll);

### read in the distancetoTSS.sum.bed

open TSS, "$distancetoTSS" or die "cant open $distancetoTSS";

my %plot;
while (my $line =<TSS>)
{
    chomp $line;
    my ($number, $distance)=split(/\t/, $line);
    $roll{$distance}=$number;
}
close TSS;
#print Dumper (\%roll);

#dereference and print hash
open PLOT,  '>',"TSS_ObsVsAnn.plot" or die "TSS_ObsVsAnn.plot";
foreach my $position (sort by_number keys %roll)
{
    my $count = $roll{$position};
    print PLOT "$position\t$count\n";
}
close PLOT;

# ascii sort
sub by_number {($a <=> $b);}

exit;
