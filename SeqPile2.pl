#!/usr/bin/perl -w


=head1 NAME

   SeqPile2.pl 

=head1 SYNOPSIS

   SeqPile2.pl -coor <Coordinates of required regions in 3 or 6 column file> -build <mm9, hg18 or hg19> -mot <sequence motif eg [CG]TG[AC].GG> -win <length of window upstream and downstream, default 2000> -bin <size of bin, default 10bp> -name Fileout -rev <n, default in y>

=head1 EXAMPLE

   SeqPile2.pl -coor peaks.bed -build mm9 -win 3000 -bin 50 -mot [CG]TG[AC].GG -name MyTest -len 6

=head1 OPTIONS

   -coor 	: input file: coordinate in as chr start stop (bed format)
   -build 	: (mm9, hg18 and hg19 at moment)
   -win 	: window upstream and downstream of midpoint Default = 2000
   -bin 	: window into which to bin reads Default = 10
   -mot 	: Motif to quantify usig basic regex, [CG]TG[AC].GG,which equals C or G, TG, A or C, any base, GG
   -rev 	: to also count the reverse and complement of motif, y or n, default is y
   -name 	: base name for output files, variables such as motif and bind added to file name
   -length	: length of the motif

 
=head1 DESCRIPTION

    Quantifies the density of a sequence motif in a window around the centre of peaks or genomic coordinates.  If 6 column bed used with strand as 6th column [+/-] then will reverse and complement input sequence.
    Outputs a png from GNUpot of distribution and text file of the position of bins in base pairs in column 1 and the number of occurence of motifs in column 2. Minimum input requirment is 3 column; chr, start and stop - tabseparated. 

=head1 AUTHOR

   Jim Hughes (University of Oxford)

=head1 Update Record

    21.6.11 JHughes, recoded by JHughes  1.2.2018

=cut

#-------------------------------------------------------------------------------

	my $help=0;
	my $man=0;	

use strict;
use Pod::Usage;
use Getopt::Long;
use Bio::DB::Sam;
use Chart::Gnuplot;
use Data::Dumper;
$Data::Dumper::Sortkeys  = 1;

my %coordinates;
my %Dens1;
my $bin = 10;
my $windowLength = 2000;
my $reverse = "y";

# FAST retieval of sub-sequence from (indexed) fasta file:   (to index:  samtools faidx /full/path/to/fasta.file ) 8.0

&GetOptions ("coor=s"=>\my $coors,
	     "build=s" => \my $build,
		 "name=s" => \my $file_name,
	     "mot=s" => \my $mot,
	     "bin=s" => \$bin,
	     "win=i" => \$windowLength,
	     "rev=s" => \$reverse,
		 "len=i" => \my $motiflength,
	     "h|help"=>\$help,
	     "man"=>\$man);


	pod2usage(1) if $help;
	pod2usage(-verbose=>2) if $man;
	pod2usage(2) unless ($coors);
	
  open (OUTPUT, ">$file_name\_$bin\_$mot\_$windowLength\_Rev\_$reverse.txt");
	
  $mot =~ tr/tgca/TGCA/;

#print "motif entered = $mot\n\n";
#my $motifsize = length($mot);
#print "motif size = $motifsize bp\n\n";

  my $fai;
  my $revMot = reverse(uc($mot));
  $revMot =~ tr/ACGT/TGCA/;
  $revMot =~ tr/][/[]/;

  
    if ($build eq 'mm9'){
	$fai = Bio::DB::Sam::Fai->load('/databank/raw/mm9/mm9.fa');  
    }
    if ($build eq 'hg18'){
	$fai = Bio::DB::Sam::Fai->load('/databank/raw/ens_human_ARCHIVE/NCBI36/ens_human_chrs/all_chrs.fasta');  
    }
    if ($build eq 'hg19'){    
	$fai = Bio::DB::Sam::Fai->load('/databank/raw/hg19_full/hg19_full.fa');
	}
	
my $totalpos = ($windowLength * 2) + 1;
my $posCount=1;

until ($posCount == $totalpos){
	
	$Dens1{$posCount}=0;
	$posCount++;
	}


open(INFO1, $coors) || die $!;

while (<INFO1>) {
   chomp;
   unless ($_ =~ /track/g){
	my ($chr, $start, $stop, $Name, $Score, $Strand)=split(/\t/);
	$coordinates{$chr}{$start}{$stop}=$Strand;

   }
}
	     

## set a hash to catch the cht start stops and relative position of the motif
my %sasquatch;

my $counteinput = 0;

foreach my $Schr ( sort keys %coordinates){
    foreach my $Sstart (sort by_number keys %{$coordinates{$Schr}}){
		foreach my $Sstop (keys %{$coordinates{$Schr}{$Sstart}}){
		
		my $storedStrand = $coordinates{$Schr}{$Sstart}{$Sstop};
	    $counteinput++;
	    
		my $midcoor;
		
		if ($Sstart == $Sstop){
			$midcoor = $Sstart;
		} elsif ($Sstop > $Sstart) {
			my $mid = int($Sstop - $Sstart)/2;
			$midcoor = $Sstart + $mid;
		} else {
			die "Dodgy input file format\n";
		}
		
		
		
		
	    my $windowstart = int($midcoor - $windowLength);
	    my $windowstop = int($midcoor + $windowLength);
	    
		#print "$Schr $windowstart $windowstop\n";
   
	    my $fai_location = $Schr . ':' . $windowstart . '-' . $windowstop;
	    my $subseq = $fai->fetch($fai_location);
    
	    $subseq =~ tr/tgca/TGCA/;
		
		if ($storedStrand =~ /\-/){
			
			$subseq = reverse(uc($subseq));
			$subseq =~ tr/ACGT/TGCA/;
						
		}
    
	    pos($subseq) = 0;
		
		## match the motif to the sequence rev
		
	    while ($subseq =~ m/($mot)/g) {
			
			## find the end of the match position in the sequence
			my $end_offset = pos($subseq);
			
			## find the start by subtracting the length of the motif
			#my $start_offset = ($end_offset - (length($mot)));
			my $start_offset = ($end_offset - $motiflength);
			$Dens1{$start_offset}++;
			
			$sasquatch{$Schr}{$Sstart}{$Sstop}{'fwdmotifs'}= $start_offset;   # assign the coordinates and position of the motif to the sasquatch hash
	    }
		if ($reverse eq "y"){
			unless ($mot eq $revMot){
	    
		## match the motif to the sequence rev
			while ($subseq =~ m/($revMot)/g) {
				my $end_offset = pos($subseq);
				$Dens1{$end_offset}++;
				
			#$sasquatch{$Schr}{$Sstart}{$Sstop}{'revmotifs'}= $end_offset;   # assign the coordinates and position of the motif to the sasquatch hash (N.B this is the end of the motif upstream of the TSS)
		}
	    }
	}
}}}

my $bincounter = 1;
my %binscore;

my $endbin = $bin;

my $Position = $windowLength * -1;


until ($bincounter == ($windowLength * 2) + 1 ){
	
	if ($bincounter == ($endbin) + 1) {
	
		$Position = $endbin - $windowLength;
		$endbin = $endbin + $bin;
	
	}

	if (exists $Dens1{$bincounter}) {
	    $binscore{$Position} +=  $Dens1{$bincounter};
	}
	
	$bincounter++;
}

#print Dumper (\%sasquatch);

    
my @x;
my @y;   
    
foreach my $BinPosition ( sort by_number keys %binscore){
    
    
    my $finalscore = $binscore{$BinPosition};
    my $adjFinalScore = $finalscore/$counteinput;
    print OUTPUT "$BinPosition\t$finalscore\n";
	
	push @x, $BinPosition;
	push @y, $finalscore;
	
}


# Chart object
my $chart = Chart::Gnuplot->new(
    output => "$file_name\_$bin\_$mot\_$windowLength\_Rev\_$reverse.png",
);

# Data set object
my $dataSet = Chart::Gnuplot::DataSet->new(
    xdata => \@x,
    ydata => \@y,
	style => "lines",
    title => "$file_name",
);

$chart->plot2d($dataSet);

system ("display -alpha Deactivate $file_name\_$bin\_$mot\_$windowLength\_Rev\_$reverse.png &") == 0 or die "couldn't display $file_name\_$bin\_$mot\_$windowLength\_Rev\_$reverse.png\n";

    
    
##################subs###################


sub by_number {
	($a <=> $b);
	}

exit;
