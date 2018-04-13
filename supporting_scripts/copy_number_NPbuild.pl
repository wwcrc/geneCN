#! /usr/bin/perl

# copy_number_NPbuild.pl - A modification of depth-for-geneCN.pl to generate coverage data for
# a panel of normal (non-cancer) samples and add the median depth for each genomic CN 
# window to the static file that will be used by depth-for-geneCN.pl.
# Copyright (C) 2015 14MG
# Copyright (C) 2017-2018 University of Glasgow
# Author (s): Susie Cooke, John Marshall
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Required inputs are a list of bam files, a file of the sample genders in each bam file 
# and a design-specific bed file.

use strict;
use warnings;

use Getopt::Long;
use Bio::DB::HTS;

my $version = '1.1';


# Read in command line and options

my $bed_file;
my $outfile;
my $gender_file;
my $ignore_gender = 0;

GetOptions ("b=s" => \$bed_file,
			"s=s" => \$gender_file,
			"x" => \$ignore_gender,
			"o=s" => \$outfile,
			"version" => \&version,
			"help" => \&help )
	or die "Error in command line arguments\n";

unless (scalar @ARGV >= 1 && defined $bed_file && (defined $gender_file || $ignore_gender)){
	die "Error in command line. Try --help\n";
}
if (defined $gender_file && $ignore_gender) {
	die "-s and -x cannot be used together\n";
}

my @bam_files = @ARGV;
my %bams;

foreach (@bam_files) {
	$bams{$_}++;
}

# Parse file of genders

my @females;
my @males;

if (defined $gender_file) {
	open (FILE1, '<', $gender_file)
		or die "Cannot open file: $gender_file: $!\n";
	
	while (<FILE1>) {
		chomp;
		my ($sample, $gender) = split (/\t/, $_);
		if (lc($gender) eq 'female') {
			if (exists $bams{$sample}) {
				push (@females, $sample);
			}
		}
		elsif (lc($gender) eq 'male') {
			if (exists $bams{$sample}) {
				push (@males, $sample);
			}
		}
		else {die "Unrecognised gender of $gender for sample: $sample\n";}
	}
	close FILE1;

	if (scalar @females < 3 || scalar @males < 3) {
		die "Insufficient samples to generate median values for chrX/chrY"
	}
}

# Prepare output file and add header

open (OUTFILE, '>', $outfile)
	or die "Cannot open file for writing: $outfile: $!\n";
print OUTFILE "#Chrom\tStart\tEnd\tGC_proportion\tFeature\tNP_median_depth\n";

# Main processing

my %feature_counts;	# A counter for how many windows each feature has
open (BEDFILE, '<', $bed_file)
	or die "Cannot open file: $bed_file: $!\n";

while (<BEDFILE>) { # iterate through bed file
	next if $_ =~ /^#/; # Skip header line(s) if present
	chomp;
	my ($chr, $start, $end, $GC, $feature) = split ('\t', $_);
	$start++; # Shift the first coordinate to be 1-based
	$feature_counts{$feature}++;
	my $median;
	if ($chr =~ /([cC][hH][rR])?[Xx]$/) {
		my @depths;
        for my $file (@females) {# Get depths from bam files
            push @depths, average_depth($file, $chr, $start, $end);
        }
		$median = median(@depths);
		if ($feature eq 'background') {
			$feature = 'other';
		}
	}
	elsif ($chr =~ /([cC][hH][rR])?[Yy]$/) {
        my @depths;
        for my $file (@males) {# Get depths from bam files
            push @depths, average_depth($file, $chr, $start, $end);
        }
        $median = 2*median(@depths); # double to prevent chrY plotting at the same level as autosomes
		if ($feature eq 'background') {
			$feature = 'other';
		}
	}
	else {
        my @depths;
        for my $file (@bam_files) {# Get depths from bam files
            push @depths, average_depth($file, $chr, $start, $end);
        }
		$median = median(@depths);
	}
	next if $median == 0;
    $start = $start-1; # change back to zero-based
	print OUTFILE "$chr\t$start\t$end\t$GC\t$feature\t$median\n";
}
close BEDFILE;
close OUTFILE;

for my $key (keys %feature_counts) {
	next if $feature_counts{$key} > 3;
	print STDERR "Warning: $key has only $feature_counts{$key} windows\n";
}

### SUBROUTINES ###

# calculate average depth in region

sub average_depth {
    my ($file, $chr, $start, $end) = @_;
    
    my $sam = Bio::DB::HTS->new(-bam => $file) or die "Cannot find bamfile $file: $!\n";
    eval {$sam->hts_index}; die "Cannot find index for $file: $!\n" if $@;
    #  calculate depth ignoring MQ zero, duplicates, non-primary, QC fail and supplementary
    my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$chr,-start=>$start,-end=>$end,-filter=>sub { my $a = shift; return ($a->flag & 3840) == 0 && $a->qual > 0 });
    if ($coverage) {
        my @data = $coverage->coverage;
        my $total = 0;
        for (@data) {
            $total += $_; # Add together the depths for all positions in the region
        }
        return $total/@data; # Divide by the size of the window (in bp)
    }
    else {
        return 0;
    }
}

# calculate median of array

sub median {
	my @numbers = @_;
	my @sorted = sort {$a <=> $b} @numbers;
	my $length = $#sorted+1;
	if ($length % 2) {
		return $sorted[($length-1)/2];
	}
	else {
		return ($sorted[$length/2]+$sorted[($length/2)-1])/2;
	}
}

# Version information

sub version {
print "Version $version\n";
exit 0;
}

# Help information

sub help {
print "Usage: copy_number_NPbuild.pl [-s FILE] [-x] -o FILENAME -b BEDFILE BAM1 BAM2 .. BAMn
Options:\t-b BED\tBed file of regions for copy number
\t\t-s FILE\tTab-delimited text file of bam file name and sample gender
\t\t-x\tUse instead of -s if bed file does not contain any regions on chrX/Y
\t\t-o STR\tName of output file
\t\t--version
\t\t--help\n";
exit 0;
}


