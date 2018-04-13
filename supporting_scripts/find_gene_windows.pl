#! /usr/bin/perl

# find_gene_windows.pl
# Copyright (C) 2017-2018 University of Glasgow
# 01 JUN 2017
# Author: Susie Cooke
# This script takes a RefSeq download of chrom, txStart, txEnd, name2 and works out the 
# maximum span of all transcripts for each gene. This window can then be expanded by a 
# specified number of bases.

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

use strict;
use warnings;
use Getopt::Long;

my $version = '1.0';

my $outfile;
my $transcript_boundaries; # refSeq/equivalent download from UCSC TableBrowser
my $gene_list;
my $extend_by = 0; # number of bases to add upstream and downstream

# options
GetOptions ("o=s" => \$outfile,
			"t=s" => \$transcript_boundaries,
			"g=s" => \$gene_list,
			"e=i" => \$extend_by,
			"version" => \&version,
			"help" => \&help )
or die "Error in command line arguments\n";

unless (defined $transcript_boundaries && defined $gene_list) {
	die "Missing command line parameter. Try --help\n";
}

# Main processing

my %starts;
my %ends;
my %chrs;

open INFILE1, '<', $transcript_boundaries
	or die "Cannot open file $transcript_boundaries: $!";
	
while (<INFILE1>) {
	chomp;
	next if $_ =~ /^#/; # Skip header line(s) if present
	my @data = split (/\t/, $_);
	unless ($#data == 3) {
		die "Wrong number of columns in transcript file\n";
	}
	next if $data[0] =~ /_/; # ignore non-standard contigs
	if (defined $starts{$data[3]} && $starts{$data[3]} <= $data[1]) {} 
	else { 
		$starts{$data[3]} = $data[1]; # extend window if new transcript spans beyond current one
		$chrs{$data[3]} = $data[0];
	}
	if (defined $ends{$data[3]} && $ends{$data[3]} >= $data[2]) {} 
	else { $ends{$data[3]} = $data[2]; } # extend window if new transcript spans beyond current one
}
close INFILE1;

open INFILE2, '<', $gene_list
	or die "Cannot open file $gene_list: $!";
open OUTFILE, '>', $outfile
	or die "Cannot open file $outfile: $!";
	
while (<INFILE2>) {
	chomp;
	if (defined $starts{$_} && defined $ends{$_}) {
		print OUTFILE $chrs{$_}, "\t", $starts{$_}-$extend_by, "\t", $ends{$_}+$extend_by, "\t", $_, "\n";
	}
	else { print "Cannot find coordinates for $_ in transcript file. Please check gene name is valid\n"; }
}
close INFILE2;
close OUTFILE;

# version
sub version {
print "Version: $version\n";
exit 0;
}

# help
sub help {
print "Usage: find_gene_windows.pl -o FILENAME -t FILE -g FILE -e NUMBER
Options:\t-t TXT\tFile containing table of transcripts
\t\t-g TXT\tFile containing one gene name per line
\t\t-o STR\tName of output file
\t\t--version
\t\t--help\n";
exit 0;
}
