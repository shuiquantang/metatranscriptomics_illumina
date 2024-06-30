#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin
#use strict;
#use warnings;
#---------------------------------------Preparation--------------------------------------------
# define parameters
use Inputs;
use SanityCheck;
use RefineSourmash;
use Parallel::Loops;
use File::Basename;
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
my $output_dir = "$work_dir/$opt_o";
my $reference_db = $opt_r;
#raw sequence folder
my $kmer_size = $opt_k;

# sample information table
my $sample_info_table; if (!$opt_t) { print('input the file of master table with -t\n');die}else{$sample_info_table=$opt_t;}

#Sanity checks
SanityCheck::check_sample_table_format($sample_info_table);


# Extract Sequence File information from the sample info table
my ($analysis_groups, $analysis_group_info, $sample_info, $cat_titles) = Inputs::read_sample_table($sample_info_table);

mkdir($output_dir);
chdir($output_dir);

# Use a subset of each file for the analysis
my @samples = sort(keys%{$sample_info});
my $rawdata_dir = "$work_dir/read_processing/Trimmomatic";
my $cmd = '';
my $cpu = `nproc`;
chomp $cpu;

my $max_memory_per_process = 4000;
my $threads = Inputs::parallel_process_allocation($max_memory_per_process);
my $pl = Parallel::Loops->new($threads);


#foreach my $i (@samples){
$pl -> foreach (\@samples, sub{
    my $i=$_;
    my $sample_output_dir = "$output_dir/$i";
    mkdir($sample_output_dir);
    chdir($sample_output_dir);
    my $path = `pwd`;
    my $log = "$sample_output_dir/sourmash.log.txt";
    system("touch $log");
    #Build signature of input fastq files
    my $j = $kmer_size;
    my $fastq = "$rawdata_dir/$i.fastq.gz";
    my $fastq_sig = "$i.k$j.sig.gz";
    my $cmd = "sourmash sketch dna -p k=$j,scaled=1000,abund $fastq -o $fastq_sig --name $i 2>&1 > /dev/null";
    if (-e $fastq) {
	Inputs::print_and_execute($cmd, $log);
    }else{
	return;
    }
    #Build match with reference databases
    my %ref_sourmash = (
            'virus'=>"$reference_db/sourmash/genbank-2022.03-viral-k$j.zip", #original sourmash genbank virus database
            'protozoa'=>"$reference_db/sourmash/genbank-2022.03-protozoa-k$j.zip", #original sourmash genbank protozoa database
            'fungi'=>"$reference_db/sourmash/genbank-2022.03-fungi-k$j.zip", #original sourmash genbank fungi database
            'GTDB'=>"$reference_db/sourmash/gtdb-rs207.genomic.k$j.zip", #original sourmash gtdb rs207 database
            'GTDB_rep'=>"$reference_db/sourmash/gtdb-rs207.genomic-reps.k$j.zip ", # original sourmash gtdb rs207 rep database
            'zymo'=>"$reference_db/sourmash/zymo_genomes.k$j.zip", # zymo genomes
               );
    my %lineage_sourmash = (
            'virus'=>"$reference_db/sourmash/lineage/genbank-2022.03-viral.lineages.csv", 
            'protozoa'=>"$reference_db/sourmash/lineage/genbank-2022.03-protozoa.lineages.csv", 
            'fungi'=>"$reference_db/sourmash/lineage/genbank-2022.03-fungi.lineages.csv", 
            'GTDB'=>"$reference_db/sourmash/lineage/gtdb-rs207.taxonomy.with-strain.csv", #
            'GTDB_rep'=>"$reference_db/sourmash/lineage/gtdb-rs207.taxonomy.reps.csv",
            'zymo'=>"$reference_db/sourmash/lineage/zymo_genomes_tax.csv",
               );
    
    #Build match with reference databases
    my %ref_sp = (
            'virus'=>"$reference_db/pathogen_id/virus.k$j.zip", # # Extracted all virus from assembly_summary_refseq.txt, and remove all phages
            'protozoa'=>"$reference_db/pathogen_id/genbank-2022.03-protozoa-k$j.zip", #original sourmash genbank protozoa database
            'fungi'=>"$reference_db/pathogen_id/fungi.k$j.zip", # Curated from assembly_summary_genbank.txt, keep one representative genome for every species
            'GTDB'=>"$reference_db/pathogen_id/gtdb-rs207.genomic-reps.dna.k$j.zip", # GTDB rep database and keep only genomes belonging to human microbiome and midog related bacteria genus.
            'zymo'=>"$reference_db/sourmash/zymo_genomes.k$j.zip", # zymo genomes.
               );
    my %lineage_sp = (
            'virus'=>"$reference_db/pathogen_id/lineage/virus_lineages.csv",
            'protozoa'=>"$reference_db/pathogen_id/lineage/genbank-2022.03-protozoa.lineages.csv",
            'fungi'=>"$reference_db/pathogen_id/lineage/fungi_lineages.csv",
            'GTDB'=>"$reference_db/pathogen_id/lineage/gtdb-rs207.taxonomy.reps.refined.csv",
            'zymo'=>"$reference_db/sourmash/lineage/zymo_genomes_tax.csv",
               );
    my $threshold_bp_1 = 1000; #50000 by default
    my $sketch_output = "$i.k$j.$threshold_bp_1.sketch.csv";
    #$cmd = "sourmash gather $fastq_sig  $ref_sourmash{'zymo'} --dna --ksize $j --threshold-bp $threshold_bp_1 -o $sketch_output 2>> $log"; #for debug purpose
    $cmd = "sourmash gather $fastq_sig $ref_sp{'virus'} $ref_sourmash{'protozoa'}  $ref_sp{'fungi'}  $ref_sourmash{'GTDB_rep'}  $ref_sourmash{'zymo'} --dna --ksize $j --threshold-bp $threshold_bp_1 -o $sketch_output 2>> $log"; #
    Inputs::print_and_execute($cmd, $log);
    
    if (-e $sketch_output) {
        #$cmd = "sourmash tax annotate -g $sketch_output -t $lineage_sourmash{'zymo'} 2>> $log";  #for debug purpose
        $cmd = "sourmash tax annotate -g $sketch_output -t $lineage_sp{'virus'} $lineage_sourmash{'protozoa'}  $lineage_sp{'fungi'}  $lineage_sourmash{'GTDB_rep'}  $lineage_sourmash{'zymo'} 2>> $log";
        Inputs::print_and_execute($cmd, $log);
        my $abun_file_old = "$i.k$j.$threshold_bp_1.sketch.with-lineages.csv";
        if (-e $abun_file_old) {
            my $abun_file = "$i.k$j.$threshold_bp_1.filtered.csv";
            $cmd ="perl $script_dir/bin/sourmash_hits_filter.pl -i $abun_file_old -o $abun_file -u 5 -l 1000 -L 50000 -S 0.05 -s 5 -G 0.01 -g 15";
            Inputs::print_and_execute($cmd, $log);
        }
    }
    
    #Tax annotate
    #use 1kb $threshold_bp and a species-level ref database
    
    my $threshold_bp_2 = 500; #50000 by default
    my $pathogen_id_controller = 0;
    if ($pathogen_id_controller == 1) {
        $sketch_output = "$i.k$j.$threshold_bp_2.sketch.csv";
        #$cmd = "sourmash gather $fastq_sig $ref_sp{'zymo'} --dna --ksize $j --threshold-bp $threshold_bp_2 -o $sketch_output 2>> $log"; #for debug purpose
        $cmd = "sourmash gather $fastq_sig $ref_sp{'virus'} $ref_sp{'protozoa'}  $ref_sp{'fungi'}  $ref_sp{'GTDB'}  $ref_sp{'zymo'} --dna --ksize $j --threshold-bp $threshold_bp_2 -o $sketch_output 2>> $log"; #
        Inputs::print_and_execute($cmd, $log);
        #Tax annotate
        if (-e $sketch_output) {
            #$cmd = "sourmash tax annotate -g $sketch_output -t $lineage_sp{'zymo'} 2>> $log";  #for debug purpose
            $cmd = "sourmash tax annotate -g $sketch_output -t $lineage_sp{'virus'} $lineage_sp{'protozoa'}  $lineage_sp{'fungi'}  $lineage_sp{'GTDB'}  $lineage_sp{'zymo'} 2>> $log";
            Inputs::print_and_execute($cmd, $log);
            my $abun_file_old = "$i.k$j.$threshold_bp_2.sketch.with-lineages.csv";
            if (-e $abun_file_old) {
                my $abun_file = "$i.k$j.$threshold_bp_2.filtered.csv";
                $cmd = "perl $script_dir/bin/sourmash_hits_filter.pl -i $abun_file_old -o $abun_file -u 5 -l 1000 -L 50000 -S 0.05 -s 5 -G 0.01 -g 15";
                Inputs::print_and_execute($cmd, $log);
            }
        }
    }
    my $output1 = "$i.k$j.$threshold_bp_1.filtered.csv";
    my $output2 = "$i.k$j.$threshold_bp_2.filtered.csv";
    my $merged = "$i.merged.csv";
    if ((-e $output2) || ($output1 eq $output2)) {
        merge_sourmash_tables($output1, '', $merged, 1);
    }else{
        merge_sourmash_tables($output1, $output2, $merged, 0.5);
    }
    chdir("..");
});
#}


sub merge_sourmash_tables{
    my $table1 = shift;
    my $table2 = shift;
    my $new_table = shift;
    my $weight = shift;
    my %genome_info1;
    my %hit_abun1;
    my @headers1 = RefineSourmash::read_sourmash_hit_table($table1, \%hit_abun1, \%genome_info1, $weight);
    my %genome_info2;
    my %hit_abun2;
    my @headers2 = RefineSourmash::read_sourmash_hit_table($table2, \%hit_abun2, \%genome_info2, $weight);
    my @headers;
    if (scalar(@headers1)>1) {
        @headers = @headers1;
    }else{
        @headers = @headers2;
    }
    my %genome_info=%genome_info1;
    foreach my $i (keys(%genome_info2)){
        if (exists($genome_info{$i})) {
            $genome_info2{$i}{'f_unique_weighted'}*=2;
        }
        my $species = $genome_info2{$i}{'species_name'};
        if (is_new_species($species, \%genome_info)) {
            $genome_info{$i} = $genome_info2{$i};
        }
    }
    #keep top 100 genome hits for prokaryote, eukaryote and virus.
    my @genomes_picked = pick_top_hits(\%genome_info);
    my @genome_order = sort{$genome_info{$b}{'f_unique_weighted'} <=> $genome_info{$a}{'f_unique_weighted'}}@genomes_picked;
    open(my $f1, ">$new_table") or die;
    my $header_line = join(",", @headers);
    print($f1 "$header_line\n");
    foreach my $i (@genome_order){
        my @coln = ();
        for(my $j=0; $j<scalar(@headers); $j++){
            my $value = $genome_info{$i}{$headers[$j]};
            if ($value =~ /\,/) {
                $value = "\"$value\"";
            }
            
            push(@coln, $value);
        }
        my $line = join(",", @coln);
        print($f1 "$line\n");
    }
    close $f1;
    
}

sub pick_top_hits{
    my $genome_info = shift;
    my %genome_info_by_domain;
    my @genomes_picked;
    foreach my $genome (keys(%{$genome_info})){
        my $lineage = $genome_info->{$genome}->{'lineage'};
        my @ranks = split(/\;/, $lineage);
        my $domain = substr($ranks[0],3,length($ranks[0])-3);
        $genome_info_by_domain{$domain}{$genome}=$genome_info->{$genome}->{'f_unique_weighted'};
    }
    foreach my $i (keys(%genome_info_by_domain)){
        my %abun = %{$genome_info_by_domain{$i}};
        my @genomes = sort{$abun{$b} <=> $abun{$a}}keys%abun;
        my $cutoff = 100;
        if ($i eq 'Bacteria') {
            $cutoff = 100;
        }
        for (my $j=0; $j<scalar(@genomes); $j++){
            if ($j<$cutoff) {
                push(@genomes_picked, $genomes[$j]);
            }
        }
    }
    return (@genomes_picked);
}

sub is_new_species{
    my $species = shift;
    my $genome_info = shift;
    foreach my $i (keys(%{$genome_info})){
        my $ref_sp = $genome_info->{$i}->{'species_name'};
        if ($species eq $ref_sp) {
            return (0); #not a new species
        }
        
    }
    return(1);# a new species
}


