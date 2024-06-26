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
my $sourmash_dir = "$work_dir/sourmash";
my $cpu = `nproc`;
chomp $cpu;

my $max_memory_per_process = 2000;
my $threads = Inputs::parallel_process_allocation($max_memory_per_process);
my $phylogeny = Inputs::ranks_to_use($script_dir);
my $cap_size = 10000000;
RefineSourmash::TaxaPolisher($rawdata_dir, $sourmash_dir, $output_dir, $sample_info, $bin_path, $threads, $reference_db, $phylogeny, $script_dir, $cap_size);
chdir($output_dir);
my $abun_dir = "$output_dir/abun_table";
mkdir($abun_dir);
merge_sourmash_outputs($abun_dir,\@samples, $script_dir);

sub merge_sourmash_outputs{ 
    my $abun_dir = shift;
    my $samples = shift;
    my $script_dir = shift;
    my %strain_abun;
    collect_abun($abun_dir, $samples, \%strain_abun, $script_dir);
    if (scalar(keys%strain_abun)>0) {
        write_abun_files(\%strain_abun, $abun_dir, $samples);
    }
    
}

sub write_abun_files{
    my $abun = shift;
    my $output_dir = shift;
    my $samples = shift;
    my %microbial_abun;
    my %host_abun;
    my $base = 100000;# pseudo total abundance assumed
    separate_host_dna($abun, $samples, \%microbial_abun, \%host_abun);
    open(my $all, ">$output_dir/all_abun_table.tsv") or die;
    open(my $euk, ">$output_dir/eukaryote_abun_table.tsv") or die;
    open(my $prok, ">$output_dir/prokaryote_abun_table.tsv") or die;
    open(my $virus, ">$output_dir/virus_abun_table.tsv") or die;
    open(my $host, ">$output_dir/host_dna_abun.csv") or die;
    my $header = join("\t", @{$samples});
    $header = "#taxon_id\t$header";
    print($all "$header\n");
    print($euk "$header\n");
    print($prok "$header\n");
    print($virus "$header\n");
    print($host "$header\n");
    my @strains = sort{$a cmp $b}keys(%microbial_abun);
    my $host_abun =0;
    my @host_DNA;
    foreach my $i (@strains){
        my @coln = ($i);
        foreach my $j (sort@{$samples}){
            my $value =0;
            if (exists($microbial_abun{$i}{$j})) {
                $value = $microbial_abun{$i}{$j}*$base;
            }
            push(@coln, $value);
        }
        my $line = join("\t", @coln);
        print($all "$line\n");
        if ($line =~ /^d__Eukaryota;/) {
	    print($euk "$line\n");
	}elsif($line =~ /^d__Virus/){
	    print($virus "$line\n");
	}else{
	    print($prok "$line\n");
	}
        
    }
    foreach my $i (keys(%host_abun)){
        my @coln = ($i);
        foreach my $j (@{$samples}){
            my $value =0;
            if (exists($host_abun{$i}{$j})) {
                $value = $host_abun{$i}{$j};
            }
            push(@coln, $value);
        }
        my $line = join("\t", @coln);
        print($host "$line\n");
        
    }
    
    close $all;
    close $euk;
    close $prok;
    close $virus;
    close $host;
}


sub separate_host_dna{
    #($abun, \%microbial_abun, \%host_abun);
    my $abun = shift;
    my $samples = shift;
    my $microbial_abun = shift;
    my $host_abun = shift;
    foreach my $i (keys%{$abun}){
        foreach my $j (keys%{$abun->{$i}}){
            $microbial_abun->{$i}->{$j} = int($abun->{$i}->{$j}*100000)/100000;
        }
    }
    open(my $f1, "<$rawdata_dir/summary.tsv") or die;
    my $header = <$f1>;
    while (my $line = <$f1>) {
        chomp $line;
        my @coln = split(/\t/, $line);
        my $abun = 0;
        if ($coln[7] =~ /\%/) {
            $abun = substr($coln[7],0,length($coln[7])-1);
        }else{
            $abun = $coln[7];
        }
        
        $host_abun->{'Host_reads'}->{$coln[0]}=$abun/100;
        $host_abun->{'Microbial_reads'}->{$coln[0]} = 1-$host_abun->{'Host_reads'}->{$coln[0]};
    }
}
sub collect_abun{
    my $abun_dir = shift;
    my $samples = shift;
    my $strain_abun = shift;
    my $script_dir = shift;
    foreach my $i (@{$samples}){
	my $abun_file = "$i/$i.after.mapping.csv";
        my %info;
        my @strains = read_csv_table($abun_file, \%info, "name");
        my $sum = 0;
        foreach my $i (@strains){
            $sum+=$info{$i}{'read_counts'};
        }
        foreach my $j (@strains){
            my $ranks = $info{$j}{'lineage'};
            my $lineage = refine_lineage($ranks, $j);#remove the level of class
            my $abun = $info{$j}{'read_counts'};
            if ($sum>0) {
                $strain_abun->{$lineage}{$i}=int($abun/$sum*1000000)/1000000;
            }
        }
    }
}

sub refine_lineage{
    my $lineage = shift;
    my $strain = shift;
    my @ranks = split(/\;/, $lineage);
    my @abbr = qw(d__ p__ c__ o__ f__ g__ s__ z__);
    my %rank;
    for (my $i=0; $i<scalar(@ranks); $i++){
        if ($ranks[$i] =~ /^$abbr[$i]/) {
            $rank{$abbr[$i]}=$ranks[$i];
        }else{
            $rank{$abbr[$i]}=$abbr[$i].$ranks[$i];
        }
    }
    my @new_ranks=();
    foreach my $i (@abbr){
        if ($i eq 'c__') {
            next;
        }elsif($i eq 'z__'){
            $strain =~ s/;/,/g;
            push(@new_ranks, "z__$strain");
        }else{
            $rank{$i} =~ s/;/,/g;
            push(@new_ranks, $rank{$i});
        }
        
    }
    
    my $new_lineage = join(";", @new_ranks);
}

sub read_csv_table{
    use Text::CSV;
    use Encode qw/encode decode/;
    my $table = shift;
    my $hash = shift;
    my $key_coln_name = shift;
    my @row_order;
    open(my $f1, "<$table") or return (@row_order);
    my $header = <$f1>;
    $header = decode("UTF-8", $header);
    chomp $header; $header =~ s/\r//g;
    my $csv = Text::CSV->new();
    $csv->parse($header);
    my @header = $csv->fields();
    my $key_coln = 0;
    for (my $i=0; $i<scalar(@header); $i++){
	if ($header[$i] eq $key_coln_name) {
	    $key_coln=$i;
	}
    }
   
    my @error;
    while (my $line = <$f1>) {
	$line = decode("UTF-8", $line);
	chomp $line; $line =~ s/\r//g;
	my $status = $csv->parse($line);
	if (!$status) {
	    push(@error, $line);
	}
	
	my @coln = $csv->fields();
	$coln[$key_coln] =~ s/^\s+|\s+$//g;
	if (length($coln[$key_coln])==0) {
	    next;
	}
	    
	push(@row_order, $coln[$key_coln]);
	for (my $i=0; $i<scalar(@header); $i++){
	    $hash->{$coln[$key_coln]}{$header[$i]}=$coln[$i];
	}
    }
    return(@row_order);
}
