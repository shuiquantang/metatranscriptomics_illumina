#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin
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
use Getopts qw(Getopts);
use File::Spec::Functions;
use DiamondAnalysis;
use Cwd 'abs_path';
use vars qw($opt_o $opt_t $opt_r $opt_k $opt_R $opt_T $opt_P);
&Getopts('t:o:r:k:R:T:P:');
my $work_dir = `pwd`;
chomp $work_dir;
my $script_dir = dirname(__FILE__);
$script_dir = abs_path($script_dir);
my $phylogeny = Inputs::ranks_to_use($script_dir);

# in the reference database folder, there is three folders: centrifuge, AMR_genes, virulence_genes

# sample information table
my $sample_info_table; if (!$opt_t) { print('input the file of master table with -t\n');die}else{$sample_info_table=$opt_t;}

#Sanity checks
SanityCheck::check_sample_table_format($sample_info_table);
#SanityCheck::check_rawdata_files($sample_info_table,$rawdata);

chdir($work_dir);

# Extract Sequence File information from the sample info table
my ($analysis_groups, $analysis_group_info, $sample_info, $cat_titles) = Inputs::read_sample_table($sample_info_table);

# Run AMR_diamond
my $diamond_dir = $opt_o;
mkdir($diamond_dir);

my $trim_folder = "$work_dir/read_processing/Trimmomatic";
my $AMR_folder = "$diamond_dir/AMR";
mkdir ($AMR_folder);
my $taxa_abun_table = "$work_dir/sourmash/abun_table/all_abun_table.tsv";
my %species_abun;
Inputs::read_abun_table($taxa_abun_table, \%species_abun);

my $max_memory_per_process = 30000;
my $threads = Inputs::parallel_process_allocation($max_memory_per_process);
DiamondAnalysis::AMR_Diamond($trim_folder, $AMR_folder, $sample_info, \%species_abun, $threads);

# Run virulence_diamond
my $vir_folder = "$diamond_dir/virulence_genes";
mkdir ($vir_folder);
DiamondAnalysis::viru_Diamond($trim_folder, $vir_folder, $sample_info, \%species_abun, $threads);


