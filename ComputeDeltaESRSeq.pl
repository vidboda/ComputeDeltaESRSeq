#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use Net::Ping;
use modules::sequence;

$Getopt::Std::STANDARD_HELP_VERSION = 1;

###############################################################################################################
##	Script to retrieve all info necessary to compute a delta ESR score for a series of genomic variants  ##
##	david baux 09/2017										     ##
##	david.baux@inserm.fr										     ##
###############################################################################################################


#######
# the idea is to:
#	-retrieve the surrouding sequence of the substitution using http://togows.org/api/ucsc/hg19/$chr:$x-$y" (we need 5 before and 5 after the variant)
#	-compute 6 ESR scores for all hexamers overlapping the variant (wt and mut) using the hexamers.csv file adapted from Ke et al., 2011
#	-done...
#######

if (! -f 'data/hexamers.txt') {die "\nNo hexamer file in the data directory\n"}

#check togows.org availability
my $p = Net::Ping->new();
if (!defined($p->ping("togows.org", 1))) {die "\n togows.org is not reachable, please check your internet connection\n"}

#deal with options
my (%opts, $list, $genome,$path);
getopts('l:g:', \%opts);

if ((not exists $opts{'l'}) || ($opts{'l'} !~ /\.txt$/o) || (not exists $opts{'g'})  || $opts{'g'} !~ /hg(19|38)$/o) {
	&HELP_MESSAGE();
	exit
}

if ($opts{'l'} =~ /(.+)([^\/]+)\.txt$/o) {$path = $1.$2; $list = $2} #get file path and prefix
elsif ($opts{'l'} =~ /([^\/]+)\.txt$/o) {$list = $1; $path = $1}
if ($opts{'g'} =~ /hg(19|38)/) {$genome = "hg$1"}

#two possible inputs
#chr	pos	ref	alt	strand
#hgvs	strand
my %ESR;
&main();

sub main {
	open F, "$path.txt" or die $!;
	my ($chr, $pos, $start, $end, $strand, $sequence_obj);
	while (my $line = <F>) {
		if ($line !~ /^#/o) {#ignore comments
			chomp($line);
			my @input = split(/\t/, $line);
			if ($#input == 1) {#hgvs case
				if ($input[1] =~ /^[\+-]$/o && $input[0] =~ /^(chr[\dXY]{1,2}):g\.(\d+)([ATGCatgc])>([ATCGatgc])$/o) {
					$sequence_obj = sequence->new($1, $2, uc($3), uc($4), $input[1]);
				}
				else {die "ERROR: bad format for input $line in $list.txt $input[0] - $input[1]\n"}
			}
			elsif ($#input == 4) {#chr pos case
				if (($input[0] =~ /^chr[\dXY]{1,2}$/o || $input[0] =~ /^[\dXY]{1,2}$/o) && $input[1] =~ /^\d+$/o && $input[2] =~ /^[ATCGatgc]$/o && $input[3] =~ /^[ATCGatgc]$/o && $input[4] =~ /^[\+-]$/o) {
					$sequence_obj = sequence->new($input[0], $input[1], $input[2], $input[3], $input[4]);
				}
				else {die "ERROR: bad format for input $line in $list.txt $input[0] $input[1] $input[2] $input[3] $input[4]\n"}
			}
			else {die "ERROR: bad format for input $line in $list.txt $#input Don't forget the strand!!!!\n"}
			$sequence_obj->getSurroundings();			
			print "Computing\t";
			$sequence_obj->toPrint();
			my ($wtESR, $mtESR, $ESR) = $sequence_obj->ESR();
			#$objects{$sequence_obj->getChr()."-".$sequence_obj->getPos()."-".$sequence_obj->getRef()."-".$sequence_obj->getAlt()} = $sequence_obj;
			$ESR{$sequence_obj->getChr()."-".$sequence_obj->getPos()."-".$sequence_obj->getRef()."-".$sequence_obj->getAlt()} = [$sequence_obj, $wtESR, $mtESR, $ESR];
		}
	}
	open G, '>results/'.$list.'_ESR.txt' or die $!;
	print G "#Chr\tPosition\tRef\tAlt\tHGVS\twtESR\tmtESR\tDeltaESR\n";
	foreach my $key (sort(keys %ESR)) {
		my $seq_obj = $ESR{$key}->[0];
		my $hgvs = $seq_obj->getChr().":g.".$seq_obj->getPos().$seq_obj->getRef().">".$seq_obj->getAlt();
		print G $seq_obj->getChr()."\t".$seq_obj->getPos()."\t".$seq_obj->getRef()."\t".$seq_obj->getAlt()."\t".$hgvs."\t".$ESR{$key}->[1]."\t".$ESR{$key}->[2]."\t".$ESR{$key}->[3]."\n";
	}
	close G;
}


exit;

sub HELP_MESSAGE {
	print "\nUsage: perl -T ComputeDeltaESRSeq.pl -l path/to/variant_list.txt -g genome_version \nSupports --help or --version\n\n
### This script computes the Delta ESR Seq score for exonic substitutions
### input format must be a file of format:
### either tabulated:
### Chr	Pos	Strand
### chr1	87466579	+
### or variants in HGVS genomic format plus the strand of interest, e.g.:
### Var	Strand
### chr1:g.87466579C>T	+
### -l txt file, list of variants
### -g genome version, hg19/hg38
### output is a tabulated file with variants (sorted), ESR WT, ESR mutant, Delta ESR
### see README.md for installation instructions
### WARNING: providing exonic positions is your responsibility!!!
### the method is adapted from di Giacomo et al., 2016 and the raw data (hexamers scores) from Ke et al., 2011.
### contact: david.baux\@inserm.fr\n\n"
}

sub VERSION_MESSAGE {
	print "\nVersion 1.0 13/09/2017\n"
}