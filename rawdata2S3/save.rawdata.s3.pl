#!/usr/bin/perl
#use strict;
#use warnings;
use vars qw($opt_i $opt_w);
&Getopts('i:w:');
open(my $table, "<$opt_i") or die;
if (!$opt_w) {
    $opt_w = 0;
}

my %s3;
while (my $line = <$table>) {
    chomp $line;
    my @coln = split(/,/, $line);
    my $random_string = $coln[1];
    my $rawdata = $coln[2];
    my $report = $coln[3];
    my $rawdata_s3 = "epiquest/epiquest_$coln[0]/$random_string/rawdata/$rawdata";
    my $report_s3 = "epiquest/epiquest_$coln[0]/$random_string/analysis/$report";
    #check if the file exists locally
    if (-e $rawdata and -e $report) {
        @{$s3{$coln[0]}{'rawdata'}}=($rawdata, $rawdata_s3);
        @{$s3{$coln[0]}{'report'}}=($report, $report_s3);
    }else{
        print ("$rawdata or $report does not exist!!\n");
        exit;
    }
    #check if the file exists on s3
    my $rawdata_s3_exist = `aws s3 ls s3://$rawdata_s3`;
    chomp $rawdata_s3_exist;
    my $report_s3_exist = `aws s3 ls s3://$report_s3`;
    chomp $report_s3_exist;
    if ($opt_w == 1) {
    }else{
        if ($report_s3_exist or $rawdata_s3_exist) {
            print ("$rawdata or $report exists on s3!!\n");
            exit;
        }
    }
       
}

close $table;

open(my $f1, ">download_links.txt") or die ;

foreach my $i (keys%s3){
    my @rawdata = @{$s3{$i}{'rawdata'}};
    my @report = @{$s3{$i}{'report'}};
    system("aws s3 cp ./$rawdata[0] s3://$rawdata[1] --acl public-read-write");
    system("aws s3 cp ./$report[0] s3://$report[1] --acl public-read-write");
    print($f1 "https://s3.amazonaws.com/$rawdata[1]\n");
    print($f1 "https://s3.amazonaws.com/$report[1]\n");
}

close $f1;

# subroutine to read input parameters, obtained on-line
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