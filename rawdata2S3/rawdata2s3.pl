#!/usr/bin/perl
use strict;
use warnings;
#inputs -s rawdata_folder -m metadata_folder -t manual_template_file -l project_info_file -d run_date -x sequencer -w controller to overwrite
use File::Basename;
use File::Copy qw(copy move);
use File::Path qw(make_path remove_tree);

use vars qw($opt_s $opt_l $opt_d $opt_m $opt_t $opt_w $opt_x);
&Getopts('s:l:d:m:t:w:x:');         

# machine code
my $sequencer = 'MiniSeq1';
if ($opt_x) {$sequencer = $opt_x;}

my $WorkDir = `pwd`;
chomp $WorkDir;

# create a log file
open(my $log, ">log.txt") or die;

# rename the raw seqs
my $run_date = $opt_d;
my $new_folder = "new_data";
my $meta_data_folder = $opt_m;
rename_raw_reads($opt_s, $new_folder);
my $manual_template_file = $opt_t;
my %projects;
my $overwrite=0;
if ($opt_w){$overwrite=$opt_w};

# read project information table
open(my $info_table, "$opt_l") or die;
my %str;
my %metadata;
my %sample_info_table;
my @order;
while (my $line = <$info_table>) {
    chomp $line;
    $line =~ s/"//g;
    if ($line =~ /^#/){
        next;
    }
    my @coln = split(/\,/, $line);
    if (scalar@coln>3) {
        $str{$coln[0]}=$coln[1];
        $metadata{$coln[0]}=$coln[3];
	$sample_info_table{$coln[0]}=$coln[2];
        push (@order, $coln[0]);
    }
}

# check if all samples have corresponding raw sequencing data files
chdir ("$WorkDir/$meta_data_folder");
foreach my $i (keys%metadata){
    my %samples;
    my %sample_label;
    my $file = $metadata{$i};
    open(my $table, "<$file") or die;
    while (my $line = <$table>) {
	chomp $line;
	$line =~ s/"//g;
	if ($line =~ /^#/) {
	    next;
	}else{
	    my @coln = split(/\,/, $line);
	    my $sample_id = $coln[1].'_'.$coln[0];
	    $samples{$sample_id}++;
	    $sample_label{$coln[5]}++;
	}
	
    }
    $file = $sample_info_table{$i};
    open(my $info_table, "<$file") or die;
    while (my $line = <$info_table>) {
	chomp $line;
	$line =~ s/"//g;
	if ($line =~ /^#/) {
	    next;
	}
	my @coln = split(/\,/, $line);
	if (scalar(@coln)<2) {
	    next;
	}
	
	if (exists($sample_label{$coln[1]})) {
	    delete($sample_label{$coln[1]});
	}else{
	    print ($log "sample $coln[1] is not included in the master table\n");
	}
	
    }
    foreach my $i (keys%sample_label){
	print($log "sample $i is not included in the sample information table\n");
    }
    
    foreach my $i (keys%samples){
	my $R1 = "$i\_R1.fastq.gz";
	my $R1_path = "$WorkDir/$new_folder/$R1";
	my $R2 = "$i\_R2.fastq.gz";
	my $R2_path = "$WorkDir/$new_folder/$R2";
	if (! -e $R1_path) {
	    print ($log "Warning:$R1 does not exist!\n");
	}
	if (! -e $R2_path) {
	    print ($log "Warning:$R2 does not exist!\n");
	}
	
    }
    
}


# zip raw sequence files
chdir ("$WorkDir/$new_folder");
my $zip_folder = "$WorkDir/zip_folder";
make_path($zip_folder);

my %rawdata_s3;
foreach my $i (keys%str){
    my @files = @{$projects{$i}};
    my $files_str = join(' ', @files);
    my $zip_file = "$i.rawdata.zip";
    #zip files
    if (scalar@files>1) {
        system("zip $zip_file $files_str");
        system("rm $files_str");
        system("mv $zip_file $zip_folder")
    }
}



# create manuals for each project from template and save to s3
foreach my $i (@order){
    chdir("$WorkDir");
    mkdir("manuals");
    my $projectID = $i;
    my $random_string = $str{$i};
    my $rawdata = "$i.rawdata.$run_date.zip";
    my $metadata = "$i.master.table.$run_date.txt";
    my $sample_info = $sample_info_table{$i};
    if (length($sample_info)>0) {
	$sample_info = "$i.sample.info.$run_date.csv";
    }
    open(my $manual, ">manuals/$i.manual") or die;
    open(my $template, "<$manual_template_file") or die;
    while (my $line = <$template>) {
	chomp $line;
	if ($line eq 'projectID=') {
	    $line.="\"$projectID\"";
	}elsif($line eq 'runID='){
	    $line.="\"$sequencer.$run_date\"";
	}elsif($line eq 'random_str='){
	    $line.="\"$random_string\"";
	}elsif($line eq 'rawdata='){
	    $line.="\"$rawdata\"";
	}elsif($line eq 'metadata='){
	    $line.="\"$metadata\"";
	}elsif($line eq 'tag='){
	    $line.="\"$run_date\"";
	}elsif($line eq 'sampleinfo='){
	    $line.="\"$sample_info\"";
	}
	print ($manual "$line\n");
    }
    close $manual;
    close $template;
    
}

#upload rawdata zip file, metadata table, and manual file to s3
foreach my $i (keys%str){
    chdir("$WorkDir");
    my $rawdata_path;
    my $metadata_path;
    my $manual_path;
    my $permission;
    if ($i =~ /^md[0-9]+/) {
	$rawdata_path = "midog/Projects/$i/rawdata";
	$metadata_path = "midog/Projects/$i/metadata";
	$manual_path = "midog/Projects/$i/manual";
	$permission = 0;
    }else{
	$rawdata_path = "epiquest/epiquest_$i/$str{$i}/rawdata";
	$metadata_path = "epiquest/epiquest_$i/$str{$i}/metadata";
	$manual_path = "epiquest/epiquest_$i/$str{$i}/manual";
	$permission = 1;
    }
    my $s3_rawdata_file = "$i.rawdata.$run_date.zip";
    my $local_rawdata_file = "$zip_folder/$i.rawdata.zip";
    my $s3_metadata_file = "$i.master.table.$run_date.txt";
    my $local_metadata_file = "$meta_data_folder/$metadata{$i}";
    my $s3_manual_file = "$i.manual";
    my $local_manual_file = "manuals/$i.manual";
    my $s3_sample_info_file = "$i.sample.info.$run_date.csv";
    my $local_sample_info_file = "$meta_data_folder/$sample_info_table{$i}";
    upload2S3($local_rawdata_file, $rawdata_path, $s3_rawdata_file, $permission);
    upload2S3($local_metadata_file, $metadata_path, $s3_metadata_file, $permission);
    upload2S3($local_sample_info_file, $metadata_path, $s3_sample_info_file, $permission);
    upload2S3($local_manual_file, $manual_path, $s3_manual_file, $permission);
}

# create summary file
chdir ("$WorkDir");
open(my $f1, ">summary.txt") or die;
foreach my $i (@order){
    my $projectID = $i;
    my $random_string = $str{$i};
    my $rawdata = "$i.rawdata.$run_date.zip";
    my $metadata = "$i.master.table.$run_date.txt";
    print ($f1 "$projectID\t$random_string\t$rawdata\t$metadata\n");
}
close $f1;

#save QC files to s3 by the date and the sequencer
chdir ("$WorkDir/$new_folder");
system("zip QC.zip *");
if (-e "QC.zip") {
    system("mv QC.zip $zip_folder");
    system("aws s3 cp $zip_folder/QC.zip s3://zymo-files/Shuiquan/MiSeq_QC/$sequencer.$opt_d/QC.zip");
}

# upload project information table, log file, summary file onto S3 QC folder.
chdir("$WorkDir");
system("aws s3 cp ./$opt_l s3://zymo-files/Shuiquan/MiSeq_QC/$sequencer.$opt_d/$opt_l");
system("aws s3 cp ./log.txt  s3://zymo-files/Shuiquan/MiSeq_QC/$sequencer.$opt_d/log.txt");
system("aws s3 cp ./summary.txt  s3://zymo-files/Shuiquan/MiSeq_QC/$sequencer.$opt_d/summary.txt");

# subroutines
sub rename_raw_reads{
    my $old_folder = shift;
    my $new_folder = shift;
    opendir(my $dir, $old_folder) || die "can't opendir $opt_s: $!";
    my @read_files = grep { !/^\./ } readdir($dir);
    #my @read_files = readdir($dir);
    make_path($new_folder);
    foreach my $i (@read_files){
        my @dots = split(/\./, $i);
        my @underscores = split(/\_/, $dots[0]);
        if (scalar@underscores <4) {
            print ($log "$i is not formated properly\n");
            next;
        }
        my @dashes = split(/\-/, $underscores[0]);
	# if the file is a paired-end file
	if ($i =~ /_R1_|_R2_/) {
	    push (@dashes, $underscores[3]);
	}
	
        my $id = join('_', @dashes);
        my @labels = ($id, $dots[1],$dots[2]);
        my $new_name = join('.', @labels);
        copy ("$old_folder/$i", "$new_folder/$new_name");
    }
}

sub upload2S3{
    my $local_file = shift;
    my $s3_path = shift;
    my $s3_file = shift;
    my $permission = shift;
    my $s3_exist = `aws s3 ls s3://$s3_path/$s3_file`;
    chomp $s3_exist;
    if ($s3_exist) {
        
        if ($overwrite == 1) {
	    if ($permission == 1) {
		system("aws s3 cp $local_file s3://$s3_path/$s3_file --acl public-read-write");
	    }else{
		system("aws s3 cp $local_file s3://$s3_path/$s3_file");
	    }
        }else{
	    print ($log "$s3_file exists on s3; to overwrite: -w 1\n");
            exit;
        }
    }else{
        if ($permission == 1) {
	    system("aws s3 cp $local_file s3://$s3_path/$s3_file --acl public-read-write");
	}else{
	    system("aws s3 cp $local_file s3://$s3_path/$s3_file");
	}
    }
    
}

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
