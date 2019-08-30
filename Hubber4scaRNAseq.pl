#!/usr/bin/perl -w
use strict;

#######################################################################################################################
#       A script for ceating a UCSC hub containing fwd and rev strand coverage of scaRNAseq full RNA molecules
#
#                           To run, place in folder with *info.txt and type:
#
#                                nohup perl Hubber4scaRNAseq.pl info.txt & 
#
#                          Coded by Martin Larke, University of Oxford 30-Aug-2019
#
#
#                   *info.txt is a tab delimited text file with the following values:
#
# 1) $hubname = scaRNAseq_sample1
# 2) $description = scaRNAseq_sample1
# 3) $publicfolder = /public/user/Mouse/scaRNAseq/  ###N.B THIS FOLDER MUST EXIST
# 4) $genome = mm9
#
#
#######################################################################################################################


# OPEN INPUT/OUTPUTS
my $hubinfo = $ARGV[0];
open INFO, "$hubinfo" or die "Can't open $hubinfo";

#EXTRACT HUB INFO
my $hubname;
my $description;
my $publicfolder;
my $genome;

while (my $line =<INFO>)
{
    chomp $line;
    ($hubname, $description, $publicfolder, $genome) = split(/\t/, $line);
}
close INFO;

#MAKE FOLDERS
system ("nohup mkdir $publicfolder/$hubname");                      #main hub directory
system ("nohup mkdir $publicfolder/$hubname/BigWigFolder");         #bw folder
system ("nohup mkdir $publicfolder/$hubname/BigWigFolder/$genome"); #bw sub directories
system ("nohup mkdir $publicfolder/$hubname/HubFolder/");           #hub folder
system ("nohup mkdir $publicfolder/$hubname/HubFolder/$genome");    #hub sub directories

#COPY THE BIGWIGS
system ("nohup cp *.bw* $publicfolder/$hubname/BigWigFolder/$genome");

#MAKE hub.txt
open HUB, '>', "$publicfolder/$hubname/HubFolder/hub.txt" or die "Can't create hub.txt";
print HUB "hub $hubname\n";
print HUB "shortLabel $hubname\n";
print HUB "longLabel $description\n";
print HUB "genomesFile genomes.txt\n";
print HUB "email martin.larke\@gmail.com\n";
close HUB;

#MAKE genomes.txt
open GENOME, '>', "$publicfolder/$hubname/HubFolder/genomes.txt" or die "Can't create genomes.txt";
print GENOME "genome $genome\n";
print GENOME "trackDb $genome/tracks.txt\n";
close GENOME;

#MAKE trackDb.txt
open TRACKDB, '>', "$publicfolder/$hubname/HubFolder/$genome/trackDb.txt" or die "Can't create trackDb.txt";
print TRACKDB "include tracks.txt\n";
close TRACKDB;

#MAKE tracks.txt
#
# header
open TRACKS, '>', "$publicfolder/$hubname/HubFolder/$genome/tracks.txt" or die "Can't create trackDb.txt";
print TRACKS "track $hubname\n";
print TRACKS "container multiWig\n";
print TRACKS "configurable on\n";
print TRACKS "shortLabel $hubname\n";
print TRACKS "longLabel $description\n";
print TRACKS "visibility full\n";
print TRACKS "type bigWig\n";
print TRACKS "autoScale on\n";
print TRACKS "aggregate transparentOverlay\n";
print TRACKS "windowingFunction maximum\n";
print TRACKS "showSubtrackColorOnUi on\n";
print TRACKS "alwaysZero on\n";
print TRACKS "dragAndDrop subtracks\n\n"; # extra new line here for spacing
#
# track 1
print TRACKS "track TSS_Plus\n";
print TRACKS "parent $hubname\n";
print TRACKS "bigDataUrl http://userweb.molbiol.ox.ac.uk$publicfolder/$hubname/BigWigFolder/$genome/fwd.tss.bw\n";
print TRACKS "shortLabel TSS_Plus\n";
print TRACKS "longLabel TSS_Plus\n";
print TRACKS "type bigWig\n";
print TRACKS "color 153,0,204\n";
print TRACKS "priority 100\n\n";
#
# track 2
print TRACKS "track PauseSites_Plus\n";
print TRACKS "parent $hubname\n";
print TRACKS "bigDataUrl http://userweb.molbiol.ox.ac.uk$publicfolder/$hubname/BigWigFolder/$genome/fwd.tes.bw\n";
print TRACKS "shortLabel PauseSites_Plus\n";
print TRACKS "longLabel PauseSites_Plus\n";
print TRACKS "type bigWig\n";
print TRACKS "color 230,159,0\n";
print TRACKS "priority 100\n\n";
#
# track 3
print TRACKS "track TSS_Minus\n";
print TRACKS "parent $hubname\n";
print TRACKS "bigDataUrl http://userweb.molbiol.ox.ac.uk$publicfolder/$hubname/BigWigFolder/$genome/rev.tss.bw\n";
print TRACKS "negateValues on\n";
print TRACKS "shortLabel TSS_Minus\n";
print TRACKS "longLabel TSS_Minus\n";
print TRACKS "type bigWig\n";
print TRACKS "color 153,0,204\n";
print TRACKS "priority 100\n\n";
#
# track 4
print TRACKS "track PauseSites_Minus\n";
print TRACKS "parent $hubname\n";
print TRACKS "bigDataUrl http://userweb.molbiol.ox.ac.uk$publicfolder/$hubname/BigWigFolder/$genome/rev.tes.bw\n";
print TRACKS "negateValues on\n";
print TRACKS "shortLabel PauseSites_Minus\n";
print TRACKS "longLabel PauseSites_Minus\n";
print TRACKS "type bigWig\n";
print TRACKS "color 230,159,0\n";
print TRACKS "priority 100\n\n";

close TRACKS;

#MAKE HUB ADDRESS
open URL, '>', "$hubname.url.txt" or die "Can't open $hubname.url.txt";
print URL "Hub address is:\n";
print URL "http://userweb.molbiol.ox.ac.uk$publicfolder$hubname/HubFolder/hub.txt\n";
close URL;

exit;
