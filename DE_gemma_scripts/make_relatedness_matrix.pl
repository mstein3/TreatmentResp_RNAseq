#!/usr/bin/perl

use strict;
use warnings;

# Adapt this script from Katie Igartua to use with new imputation data: I've copied over the additve coef file to my current scratch/mstein3/cytoQTL_2 dir and zcat'd it.
# File: /group/nicolae-lab/users/mstein3/additive.coef.2015-12-15
#Command : perl make_relatedness_matrix.pl INPUTFILE OUTPUTFILE KINSHIP
#INPUTFUILE - order of findivs for matrix file. Must have a header column
#OUPUTFILe - name of file for output matrix
#KINSHIP - input file with pairwise kiship, one line per pair of individuals (/group/nicolae-lab/users/mstein3/additive.coef.2015-12-15)
# Check if the diagonal is around 1, or around 0.5. If around 0.5, multiply by 2 to get the GRM. 

# In this case: perl make_relatedness_matrix.pl /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_InputFiles/FINDIV_only_224HT_list.txt /group/ober-resources/users/mstein3/rna.seq.2017/DE_gemma_InputFiles/NbyN_LPSNull_Kinship_224.txt /group/nicolae-lab/users/mstein3/additive.coef.2015-12-15

my $file1=$ARGV[2]; #KINSHIP

my $outfile=$ARGV[1]; #OUTPUT MATRIX
open (OUT, ">$outfile") or die "Cannot open output file\n";


open(IN,$file1) || die "can't read $file1\n";


#MAKE DICTIONARY

my $findiv_combinations;
<IN>;
while (my $line = <IN>)
{
	next if $.==1;
	chomp $line;
	my @data = split /\t/, $line;
	$findiv_combinations->{$data[0]."_".$data[1]}=$data[2];
	$findiv_combinations->{$data[1]."_".$data[0]}=$data[2];
	}
close IN;

my $file2=$ARGV[0];


#MAKE ARRAY OF SAMPLE ORDER 

open(IN2,$file2) || die "can't read $file2\n";
my @findiv;
<IN2>;
while (my $line = <IN2>)
{
	chomp $line;
	my @data = split /\t/, $line;
	push(@findiv, $data[0]);
	}
close IN2;
	
#PRINT ORDER MATRIX IN ORDER OF ARRAY

foreach (my $findiv1=0; $findiv1<scalar(@findiv);$findiv1++){
	foreach (my $findiv2=0; $findiv2<scalar(@findiv);$findiv2++){ 
		if (defined $findiv_combinations->{$findiv[$findiv1]."_".$findiv[$findiv2]}){
			print OUT $findiv_combinations->{$findiv[$findiv1]."_".$findiv[$findiv2]};
			print OUT "\t";
		}
		else {
			print OUT "NA";
			#print $findiv[$findiv1]."_".$findiv[$findiv2];
			print "$findiv[$findiv1]";
			print "\n";
			print OUT "\t";
		}	
	}
	print OUT "\n";
	}

	
