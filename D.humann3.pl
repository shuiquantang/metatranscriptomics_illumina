#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin

use warnings;
use File::Basename;

#introduction
#fastqc -> trimmomatic -> centrifuge -> diamond
#inputs: 1.raw sequence folder, 2. sample info table,

#---------------------------------------Preparation--------------------------------------------
# define parameters
use Inputs;
use SanityCheck;
use HumannAnalysis;
use Getopts qw(Getopts);
use File::Spec::Functions;
use Cwd 'abs_path';
use Getopts qw(Getopts);
use vars qw($opt_o $opt_t $opt_r $opt_k $opt_R $opt_T $opt_P);
&Getopts('t:o:r:k:R:T:P:');

my $work_dir = `pwd`;
chomp $work_dir;
my $script_dir = dirname(__FILE__);
$script_dir = abs_path($script_dir);
my $bin_path = catdir($script_dir,'bin');
my $phylogeny = Inputs::ranks_to_use($script_dir);

#a controller to control whether or not to bypass uniref search

# sample information table
my $sample_info_table; if (!$opt_t) { print('input the file of master table with -t\n');die}else{$sample_info_table=$opt_t;}
chdir($work_dir);

# Extract Sequence File information from the sample info table
my ($analysis_groups, $analysis_group_info, $sample_info, $cat_titles) = Inputs::read_sample_table($sample_info_table);

# Run humann2
my $humann_folder = $opt_o;
my $trim_folder = "$work_dir/read_processing/Trimmomatic";
mkdir ($humann_folder);

my $max_memory_per_process = 25000;
my $humann_subset_size = 100000000; # 2M paired-end reads, equivalent to 4M single reads
my $threads = Inputs::parallel_process_allocation($max_memory_per_process);
HumannAnalysis::humann3($trim_folder, $humann_folder, $sample_info, $threads, $humann_subset_size);
