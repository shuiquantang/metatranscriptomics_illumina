#!/usr/bin/perl -I /home/stdell/Desktop/VirtualBox/scripts/metatrans/scripts/bin
use strict;
use warnings;
use File::Basename;

#introduction
#fastqc -> trimmomatic -> centrifuge -> diamond
#inputs: 1.raw sequence folder, 2. sample info table,

#---------------------------------------Preparation--------------------------------------------
# define parameters
use Inputs;
use SanityCheck;
use IlluminaAnalysis;
use Getopts qw(Getopts);
use File::Spec::Functions;
use Cwd 'abs_path';
use vars qw($opt_o $opt_t $opt_r $opt_k $opt_R $opt_T $opt_P);
&Getopts('t:o:r:k:R:T:P:');
my $work_dir = `pwd`;
chomp $work_dir;
my $script_dir = dirname(__FILE__);
$script_dir = abs_path($script_dir);
my $bin_path = catdir($script_dir,'bin');
my $phylogeny = Inputs::ranks_to_use($script_dir);
my $read_processing_dir = "$work_dir/$opt_o";
#raw sequence folder
my $rawdata='rawdata';

# sample information table
my $sample_info_table; if (!$opt_t) { print('input the file of master table with -t\n');die}else{$sample_info_table=$opt_t;}

#Sanity checks
SanityCheck::check_sample_table_format($sample_info_table);
SanityCheck::check_rawdata_files($sample_info_table, $rawdata);

chdir($work_dir);

# Extract Sequence File information from the sample info table
my ($analysis_groups, $analysis_group_info, $sample_info, $cat_titles) = Inputs::read_sample_table($sample_info_table);
mkdir ($read_processing_dir);

# fastqc
my $QC_folder = "$read_processing_dir/FastQC";
mkdir ($QC_folder);
IlluminaAnalysis::FastQC($rawdata, $QC_folder, $sample_info, $bin_path);


# trim the reads and SortMeRNA
my $trim_folder = "$read_processing_dir/Trimmomatic";
mkdir ($trim_folder);
my $max_memory_per_process_trimmomatic = 10000;
my $threads_trimmomatic = Inputs::parallel_process_allocation($max_memory_per_process_trimmomatic);
my $max_memory_per_process_ribodetector = 30000;
my $threads_SortMeRNA = Inputs::parallel_process_allocation($max_memory_per_process_ribodetector);

IlluminaAnalysis::trim_and_filter($rawdata, $trim_folder, $sample_info, $bin_path, $threads_trimmomatic, $threads_SortMeRNA, $QC_folder);

