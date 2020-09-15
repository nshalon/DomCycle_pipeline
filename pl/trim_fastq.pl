#!/usr/bin/env perl

use strict;
use warnings FATAL => qw(all);
use File::Basename;

if ($#ARGV == -1) {
	print "usage: $0 <ifn> <offset> <length> <ofn>\n";
	exit 1;
}

my $ifn = $ARGV[0];
my $offset = $ARGV[1];
my $rlen = $ARGV[2];
my $ofn = $ARGV[3];

my $iindex = 0;
open(OUT, ">", $ofn) or die;
open(IN, $ifn) or die;
while (my $line = <IN>) {
    chomp($line);
    $line = substr($line, $offset, $rlen) if ($iindex++ % 2 == 1);
    print OUT $line, "\n";
}

close(OUT);
close(IN);
