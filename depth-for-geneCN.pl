#! /usr/bin/perl

# depth-for-geneCN.pl - A script to generate coverage data for a sample of interest in preparation for copy-number calling.
# Copyright (C) 2015 14MG
# Copyright (C) 2017-2018 University of Glasgow
# Author(s): Susie Cooke, John Marshall
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

# Required inputs are the sample bam file and a design-specific bed file.

use strict;
use warnings;

use Getopt::Long;
use Bio::DB::HTS;

my $version = '2.0';

# Read in command line and options

my $bed_file;
my $output_file = 'coverage.tmp';

GetOptions (	"b=s" => \$bed_file,
		"o=s" => \$output_file,
		"version" => \&version,
		"help" => \&help )
or die "Error in command line arguments\n";

unless (scalar @ARGV == 1 && defined $bed_file){
	die "Error in command line. Try --help\n";
}

my $bam_file = $ARGV[0];

# Read in bed file

open (BEDFILE, '<', $bed_file)
	or die "Cannot open file: $bed_file: $!\n";

my @bed_lines;

while (<BEDFILE>) {
	next if $_ =~ /^#/; # Skip header line(s) if present
	push (@bed_lines, $_);
}

close BEDFILE;

# Get depths from bam files

my $sam = Bio::DB::HTS->new(-bam  =>$bam_file)
	or die "Cannot find bamfile $bam_file: $!\n";
eval {$sam->hts_index}; die "Cannot find index for $bam_file: $!\n" if $@;

open OUTFILE, ">$output_file"  # 2-argument open so that -o - works
	or die "Cannot open file for writing: $output_file: $!\n";

my $counter = 0;

foreach my $line (@bed_lines) {
	chomp $line;
	$counter++;
	my @position = split ('\t', $line);
	if ($position[5] == 0) {
		die "Column 6 of the input bedfile cannot be zero, zero found at line $counter\n";
	}
	$position[1]++; # Shift the first coordinate to be 1-based
	#  calculate depth ignoring MQ zero, duplicates, non-primary, QC fail and supplementary
	my ($coverage) = $sam->features(-type=>'coverage',-seq_id=>$position[0],-start=>$position[1],-end=>$position[2],-filter=>sub { return ($_[0]->flag & 3840) == 0 });
	my $ave_depth = 0;
	if ($coverage) {
		my @data = $coverage->coverage;
		my $total = 0;
		for (@data) {
			$total += $_; # Add together the depths for all positions in the region
		}
		$ave_depth = $total/@data; # Divide by the size of the window (in bp)
	}
	my $zero_based = $position[1]-1; # change back to zero-based
	print OUTFILE "$position[0],$zero_based,$position[2],$ave_depth\n";
}
close OUTFILE or die "Cannot close $output_file: $!\n";
exit 0;

# Version information

sub version {
print "Version $version\n";
exit 0;
}

# Help information

sub help {
print "Usage: depth-for-geneCN.pl -b BEDFILE -o FILENAME BAMFILE
Options:\t-b BED\tNormal panel bed file of regions for copy number
\t\t-o FILE\tWrite output to FILE (defaults to coverage.tmp)
\t\t--version
\t\t--help\n";
exit 0;
}
