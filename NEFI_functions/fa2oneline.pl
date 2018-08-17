#!/usr/bin/perl -w
#this will convert a fasta file where the sequence is written on multiple line to a fasta file where the sequence is written on one line.
#example: perl fa2oneline.pl sample.fna > out.fna
use strict;

my $input_fasta=$ARGV[0];
open(IN,"<$input_fasta") || die ("Error opening $input_fasta $!");

perl -pe '/^>/ ? print "\n" : chomp' input_fasta | tail -n +2
