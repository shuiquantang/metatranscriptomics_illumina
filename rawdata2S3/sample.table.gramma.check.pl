#!/usr/bin/perl
#use strict;
use warnings;
use File::Basename;
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);

# example 
# perl sample.table.gramma.check.pl -i metadata/

use vars qw($opt_i);
&Getopts('i:');
my $input_folder = $opt_i;
opendir(my $dir, $input_folder) || die "can't opendir $input_folder: $!";
my @tables = grep { !/^\./ } readdir($dir);

open(my $report, ">sample.table.gramma.report.txt") or die;
chdir($input_folder);

foreach my $i (@tables){
    if ($i =~ /\~$/) {
	next;
    }
    print ($report "\n------------Checking $i-----------\n");
    check_sample_table_format($i, $report);
}

close $report;


sub check_sample_table_format{
    my $sample_metadata_table = $_[0];
    my $report=$_[1];
    open(my $metadata, "<$sample_metadata_table") or die ("sample metadata table file, $sample_metadata_table, does not exist in work directory\n");
    my %groups;
    my %samples;
    my $die=0;
    while (my $line = <$metadata>) {
	chomp $line;
	$line =~ s/\"//g;
	if ($line eq ''){ next;}
	if ($line =~ /^#/) {next;}
	my @coln = split(/\,/, $line);
	my $group_id = $coln[3].'_'.$coln[4];
	my $sample_id = $coln[1].'_'.$coln[0];
	my $label = $coln[5];
	my $project_id = $coln[1];
	#count colns;
	my $coln_no = 0;
	foreach my $i (@coln){
	    $i =~ s/^\s+|\s+$//g; # remove white space
	    if (length($i)>0) {
		$coln_no++;
	    }
	}
	$groups{$group_id}{'id'}{$sample_id}++;
	$groups{$group_id}{'label'}{$label}++;
	$groups{$group_id}{'coln'}{$coln_no}++;
	$groups{$group_id}{'project'}{$project_id}++;
	#check the format of sample types
	my @seq_type = ('Bac16Sv13', 'Bac16Sv34', 'Bac16Sv4', 'Bac16Sv68', 'FungiITS', 'Euk18S');
	if (! grep {/$coln[4]/} @seq_type) {
	    print($report "sample $sample_id, has an unknown seq_type, $coln[4]!\n");
	    $die=1;
	}
	#check illegal letters
	foreach my $i (@coln){
            chomp($i);
            if ($i ne ''){
                if ($i =~ /[^a-zA-Z0-9\.]/) {
		    print($report "$i in the row of sample $sample_id contains illegal letters!\n");
		    $die=1;
		}
            }
		
	}
	#check total columns
	if ((scalar@coln) < 5) {
	    print($report "sample $sample_id has at least five columns separated by \',\'\n");
	    $die=1;
	}
    }
    foreach my $i (keys%groups){
	my %id = %{$groups{$i}{'id'}};
	my %label = %{$groups{$i}{'label'}};
	my @coln_counts = keys(%{$groups{$i}{'coln'}});
	my @projects = keys(%{$groups{$i}{'project'}});
	foreach my $j (keys%id){
	    if ($id{$j}>1) {
		print ($report "Group $i, sample $j, is not unique!\n");
		$die=1;
	    }
	    
	}
	foreach my $j (keys%label){
	    if ($label{$j}>1) {
		print($report "Group $i, label $j, is not unique!\n");
		$die=1;
	    }
	    
	}
	if (scalar(@coln_counts)>1) {
	    print ($report "Group $i has unequal columns!\n");
	    $die=1;
	}
	if (scalar(@projects)>1) {
	    if ($i =~ /^QC_/) {
		
	    }else{
		print ($report "warning: Group $i has more than one project_id!\n");
	    }
	}
	
    }
    close $metadata;
    
    if ($die) {
	#die("Sample info table fails sanity check\n");
    }
}

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
