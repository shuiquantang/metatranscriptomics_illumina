#!/usr/bin/perl
#use strict;
#use warnings;

use vars qw($opt_f);
&Getopts('f:');         

# create a log file

opendir(my $dir, $opt_f) || die "can't opendir $opt_f: $!";
my @genome_id = grep { !/^\./ } readdir($dir);
open(my $f2, ">acc2genome.tsv") or die;
foreach my $i (@genome_id){
    system("gunzip -c $opt_f/$i > genome.fna");
    open(my $f1, "<genome.fna") or die;
    while (my $line = <$f1>) {
        chomp $line;
        if ($line =~ /^>/) {
            my @title = split(/\ /, $line);
            my $acc = substr($title[0], 1, length($title[0])-1);
            print($f2 "$acc\t$i\n");
        }
        
    }
    close $f1;
    
}
close $f2;

sub Getopts {
    local($argumentative) = @_;
    local(@args,$_,$first,$rest);
    local($errs) = 0;
    @args = split( / */, $argumentative );
    while(@ARGV && ($_ = $ARGV[0]) =~ /^-(.)(.*)/) {
		($first,$rest) = ($1,$2);
		$pos = index($argumentative,$first);
		if($pos >= 0) {
			if($args[$pos+1] eq ':') {
				shift(@ARGV);
				if($rest eq '') {
					++$errs unless(@ARGV);
					$rest = shift(@ARGV);
				}
				eval "
				push(\@opt_$first, \$rest);
				if (!defined \$opt_$first or \$opt_$first eq '') {
					\$opt_$first = \$rest;
				}
				else {
					\$opt_$first .= ' ' . \$rest;
				}
				";
			}
			else {
				eval "\$opt_$first = 1";
				if($rest eq '') {
					shift(@ARGV);
				}
				else {
					$ARGV[0] = "-$rest";
				}
			}
		}
		else {
			print STDERR "Unknown option: $first\n";
			++$errs;
			if($rest ne '') {
				$ARGV[0] = "-$rest";
			}
			else {
				shift(@ARGV);
			}
		}
	}
    $errs == 0;
}