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
use StrainScan;
use DiamondAnalysis;
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

chdir($work_dir);

#Centrifuge Taxonomy Polisher
my $max_memory_per_process = 16000;
my $subset_size = 10000000;
my $threads = Inputs::parallel_process_allocation($max_memory_per_process);
my $trim_folder = 'read_processing/Trimmomatic';
my $abun_table_folder = 'sourmash/abun_table';
my $reads_cutoff = 10000;
my $output_dir = $opt_o;
my $top_species_number = 10; #top 10 species to keep
my %chosen_species;

my $StrainScan_ref_db_loci_s3 = "s3://zymo-files/WGS_Pipeline/shotgun_database/StrainScan_221101";
my $StrainScan_ref_dict_file = "$StrainScan_ref_db_loci_s3/StrainScan_ref_dict.tsv";

system("mkdir $output_dir");
my @species_wishlist = ();

StrainScan::select_species($subset_size, $abun_table_folder, $output_dir, $top_species_number, $reads_cutoff, \%chosen_species, $trim_folder, \@species_wishlist, $StrainScan_ref_dict_file);

StrainScan::scan_strains($trim_folder, $output_dir, \%chosen_species, $StrainScan_ref_dict_file, $StrainScan_ref_db_loci_s3);
