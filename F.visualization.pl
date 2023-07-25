#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin
use strict;
use warnings;
use File::Basename;
use QiimeAnalysis;
use TaxaHeatmap;
use FuncHeatmap;
use AlphaDiversity;
use BetaDiversity;
use LefseAnalysis;
use AnalysisPreparation;
use Parallel::Loops;

#---------------------------------------Preparation--------------------------------------------
# define parameters
use Inputs;
use SanityCheck;
use Getopts qw(Getopts);
use File::Spec::Functions;
use Cwd 'abs_path';
use vars qw($opt_o $opt_t $opt_r $opt_k $opt_R $opt_T $opt_P);
&Getopts('t:o:r:k:R:T:P:');
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

chdir($work_dir);

# Extract Sequence File information from the sample info table
my ($analysis_groups, $analysis_group_info, $sample_info, $cat_titles) = Inputs::read_sample_table($sample_info_table);

my @taxa_domains = ('all', 'prokaryote', 'eukaryote', 'virus');


#if (-e $output_dir) {system ("rm -r $output_dir");} 
mkdir($output_dir);
chdir($output_dir);

my $internal_project_id = $opt_P;
my $random_str = $opt_R;
my $rawdata_tag = $opt_T;
create_rawdata_link_table($internal_project_id, $random_str, $rawdata_tag, $sample_info);


my $max_memory_per_process = 40000;
my $threads = Inputs::parallel_process_allocation($max_memory_per_process);

my $pl = Parallel::Loops->new($threads);
my @ana_groups = keys%{$analysis_groups};

#foreach my $i (keys%{$analysis_groups}){
$pl -> foreach (\@ana_groups, sub{
    my $i=$_;
    mkdir($i);
    # make relevant copies of FastQC, Trimmomatic 
    AnalysisPreparation::rename_sample_and_prepare_inputs($work_dir, $output_dir, $i, $analysis_groups, $analysis_group_info, $sample_info, \@taxa_domains);
    chdir($i);
    my $category_number = Inputs::create_mapping_file($i, $analysis_group_info, $analysis_groups, $cat_titles);
    foreach my $j (@taxa_domains){
        mkdir ("$j");
        # cp all mapping files and sample list file to the subfolders
        system("cp *.mapping.file.txt $j/");
        system("cp sample.list.txt $j/");
        chdir ("$j");
        my $abun_table = "$work_dir/sourmash/abun_table/$j\_abun_table.tsv";
        my $new_abun_table = "abun_table.tsv";
        my $taxa_no = QiimeAnalysis::rename_sample_in_abun_table($abun_table, $analysis_groups, $analysis_group_info, $new_abun_table, $i);
        if ($taxa_no==0) {
            print("No taxa is avaliable in this group!\n");
            chdir("..");
            next;
        }
        
        # save a copy of the abundance file into the group folder, so that we can count the host DNA abundance
        if ($j eq 'all') {
            system("cp abun_table.tsv ../");
        }
        
        # biom format conversion
        my $abun_table_folder = 'abun_tables';
        QiimeAnalysis::biom_conversion($new_abun_table, $abun_table_folder, $phylogeny);
        
        # barplots
        my $barplot_folder = 'barplots';
        QiimeAnalysis::taxa_composition_barplot($abun_table_folder, $barplot_folder, $phylogeny);
        
        #heatmaps
        my $heatmap_folder = 'heatmaps';
        TaxaHeatmap::plot_heatmap_by_category($abun_table_folder, $category_number, $cat_titles, $bin_path, $heatmap_folder, $phylogeny);
        
        #alpha_diversity
        my $alpha_div_folder = 'alpha_div';
        AlphaDiversity::alpha_diversity($abun_table_folder, $analysis_groups, $category_number, $cat_titles, $bin_path, $alpha_div_folder, $phylogeny, $i);
        
        #beta_diversity
        my $beta_div_folder = 'beta_div';
        BetaDiversity::beta_diversity($abun_table_folder, $analysis_groups, $category_number, $cat_titles, $bin_path, $beta_div_folder, $phylogeny, $analysis_group_info, $i);
        
        #lefse analysis
        LefseAnalysis::plot_lefse_by_category($bin_path, $category_number, $cat_titles, $i, $analysis_group_info, $abun_table_folder, $phylogeny);
        
        system("rm *.mapping.file.txt sample.list.txt");
        chdir("..");
        
    }
    my $host_abun_table = "$work_dir/sourmash/abun_table/host_dna_abun.csv";
    my $new_abun_table = "host_dna_abun.csv";
    if (-e $host_abun_table) {
        my $taxa = QiimeAnalysis::rename_sample_in_abun_table($host_abun_table, $analysis_groups, $analysis_group_info, $new_abun_table, $i);
    }
    
    # Functional analysis vasualization
    
    if (-d "humann") {
        system("cp *.mapping.file.txt humann/"); system("cp sample.list.txt humann/");
        chdir ("humann");
        #gene_family abundance heatmap, only count the gene family with species designation
        my $abun_file = 'species_gene_fam_cpm_filt.tsv';
        my $tag = 'gene_fam_cpm';
        FuncHeatmap::plot_heatmap_by_category($abun_file, $category_number, $cat_titles, $bin_path, $tag);
        #LefseAnalysis::plot_humann_lefse_by_category($bin_path, $category_number, $cat_titles, $i, $analysis_group_info, $abun_file, $tag);
        #pathway abundance lefse, lefse can deal with the file that has pathway and pathway stratification combined.
        $abun_file = 'combined_pathway_abun_filt.tsv';
        $tag = 'pathway_abun';
        LefseAnalysis::plot_humann_lefse_by_category($bin_path, $category_number, $cat_titles, $i, $analysis_group_info, $abun_file, $tag);
        #pathway abundance heatmap without species stratification.
        $abun_file = 'pathway_abun_filt.tsv';
        $tag = 'pathway_abun';
        FuncHeatmap::plot_heatmap_by_category($abun_file, $category_number, $cat_titles, $bin_path, $tag);
        #pathway abundance heatmap with species stratification alone.
        $abun_file = 'species_pathway_abun_filt.tsv';
        $tag = 'species_pathway_abun';
        FuncHeatmap::plot_heatmap_by_category($abun_file, $category_number, $cat_titles, $bin_path, $tag);
        system("rm *.mapping.file.txt sample.list.txt");
        chdir("..");
    }
    chdir("..");
#}
});

sub create_rawdata_link_table{
    my $project_id = shift;
    my $random_str = shift;
    my $run_id = shift;
    my $sample_info =shift;
    my @samples = sort {$a cmp $b} keys%{$sample_info};
    my $rawdata_link_table = "$project_id\_rawdata_links.tsv";
    open(my $f1, ">$rawdata_link_table") or die;
    print($f1 "Internal_id\tCustomer_label\tread1\tread2\n");
    foreach my $i (@samples){
        my $label = $sample_info->{$i}{'label'};
        my $R1 = $sample_info->{$i}->{'R1'};
        if ($R1 =~ /$project_id/) {
            my $http_link_R1 = "https://zymo-microbiomics-service.s3.amazonaws.com/epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$R1";
            my $s3_link_R1 = "s3://zymo-microbiomics-service/epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$R1";
            system("aws s3api put-object-acl --bucket zymo-microbiomics-service --key epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$R1 --acl public-read");
            my $http_link_R2 ='';
            if (exists($sample_info->{$i}->{'R2'})) {
                my $R2=$sample_info->{$i}->{'R2'};
                $http_link_R2 = "https://zymo-microbiomics-service.s3.amazonaws.com/epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$R2";
                my $s3_link_R2 = "s3://zymo-microbiomics-service/epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$R2";
                system("aws s3api put-object-acl --bucket zymo-microbiomics-service --key epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$R2 --acl public-read");
            }
            print($f1 "$i\t$label\t$http_link_R1\t$http_link_R2\n");
        }
    }
    close $f1;
    my $s3_path = "s3://zymo-microbiomics-service/epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$rawdata_link_table";
    my $s3_https_path = "https://zymo-microbiomics-service.s3.amazonaws.com/epiquest/epiquest_$project_id/$random_str/rawdata/$run_id/$rawdata_link_table";
    system("aws s3 cp $rawdata_link_table $s3_path --acl public-read-write");
    system("rm -r $rawdata_link_table");
    open(my $f2, ">rawdata_link.txt") or die;
    print($f2 "$s3_https_path");
    close $f2;
    
}
