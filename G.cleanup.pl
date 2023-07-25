#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin
use strict;
use warnings;
use File::Basename;

#---------------------------------------Preparation--------------------------------------------
# define parameters
use Inputs;
use SanityCheck;
use Getopts qw(Getopts);
use File::Spec::Functions;
use Cwd 'abs_path';
use vars qw($opt_t $opt_o);
&Getopts('t:o:');
my $work_dir = `pwd`;
chomp $work_dir;
my $output_dir = $opt_o;
my $script_dir = dirname(__FILE__);
$script_dir = abs_path($script_dir);
my $bin_path = catdir($script_dir,'bin');
my $phylogeny = Inputs::ranks_to_use($script_dir);

# sample information table
my $sample_info_table; if (!$opt_t) { print('input the file of master table with -t\n');die}else{$sample_info_table=$opt_t;}

#Sanity checks
SanityCheck::check_sample_table_format($sample_info_table);
#SanityCheck::check_abundance_table($sample_info_table,"$work_dir/abun_table/all_abun_table.tsv");

chdir($work_dir);

# Extract Sequence File information from the sample info table
my ($analysis_groups, $analysis_group_info, $sample_info, $cat_titles) = Inputs::read_sample_table($sample_info_table);

my @taxa_domains = ('All', 'Prokaryote', 'Eukaryote', 'Virus');
my @ranks = ("1.Superkingdom","2.Phylum","3.Order","4.Family","5.Genus","6.Species","7.Strain");


my @ana_groups = keys%{$analysis_groups};

foreach my $i (@ana_groups){
    system("rm $output_dir/$i/AbundanceTables/AbundanceTable.csv");
    foreach my $j (@taxa_domains){
        system("rm $output_dir/$i/$j/AbundanceTables/ReadAbundance.tsv");
        foreach my $z (@ranks){
            system("rm $output_dir/$i/$j/AbundanceTables/$z/abun_table*");
        }
    }
}