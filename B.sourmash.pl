#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin
#use strict;
#use warnings;
#---------------------------------------Preparation--------------------------------------------
# define parameters
use Inputs;
use SanityCheck;
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
my $sourmash_database = "$reference_db/sourmash";
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
my @samples = keys%{$sample_info};
my $rawdata_dir = "$work_dir/read_processing/Trimmomatic";
my $cmd = '';
my $cpu = `nproc`;
chomp $cpu;

my $max_memory_per_process = 24000;
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
    my $cmd = "sourmash sketch dna -p k=$j,scaled=1000,abund $fastq -o $fastq_sig --name $i 2>> $log";
    if (-e $fastq) {
	Inputs::print_and_execute($cmd, $log);
    }else{
	return;
    }
    #Build match with reference databases
    my $genbank_v = "$sourmash_database/genbank-2022.03-viral-k$j.zip";
    my $genbank_p = "$sourmash_database/genbank-2022.03-protozoa-k$j.zip";
    my $genbank_f = "$sourmash_database/genbank-2022.03-fungi-k$j.zip";
    my $GTDB_all = "$sourmash_database/gtdb-rs207.genomic.k$j.zip";
    my $zymo_genomes = "$sourmash_database/zymo_genomes.k$j.zip";
    my $host_genomes = "$sourmash_database/host_genomes.k$j.zip";
    my $threshold_bp = 50000; #50000 by default
    my $sketch_output = "$i.k$j.sketch.csv";
    #$cmd = "sourmash gather $fastq_sig $zymo_genomes --dna --ksize $j --threshold-bp $threshold_bp -o $sketch_output 2>> $log"; #
    $cmd = "sourmash gather $fastq_sig $genbank_v $genbank_p $genbank_f $GTDB_all $zymo_genomes $host_genomes --dna --ksize $j --threshold-bp $threshold_bp -o $sketch_output 2>> $log"; #
    Inputs::print_and_execute($cmd, $log);

    #Tax metagenome
    my $gb_lineage_v = "$sourmash_database/lineage/genbank-2022.03-viral.lineages.csv";
    my $gb_lineage_p = "$sourmash_database/lineage/genbank-2022.03-protozoa.lineages.csv";
    my $gb_lineage_f = "$sourmash_database/lineage/genbank-2022.03-fungi.lineages.csv";
    my $GTDB_lineage = "$sourmash_database/lineage/gtdb-rs207.taxonomy.with-strain.csv";
    my $zymo_genome_tax = "$sourmash_database/lineage/zymo_genomes_tax.csv";
    my $host_genome_tax = "$sourmash_database/lineage/host_genomes_tax.csv";
    $cmd = "sourmash tax metagenome -g $sketch_output -t $gb_lineage_v $gb_lineage_p $gb_lineage_f $GTDB_lineage $MAGs_lineage $zymo_genome_tax $host_genome_tax --output-format krona csv_summary kreport -r order -o $i.k$j.metagenome 2>>$log";
    Inputs::print_and_execute($cmd, $log);
        
    #Tax annotate
    $cmd = "sourmash tax annotate -g $sketch_output -t $gb_lineage_v $gb_lineage_p $gb_lineage_f $GTDB_lineage $MAGs_lineage $zymo_genome_tax $host_genome_tax 2>> $log";
    Inputs::print_and_execute($cmd, $log);
    chdir("..");
});
#}

my $abun_dir = "$output_dir/abun_table";
mkdir($abun_dir);
merge_sourmash_outputs($abun_dir,\@samples, $kmer_size);


sub merge_sourmash_outputs{
    my $abun_dir = shift;
    my $samples = shift;
    my $kmer_size = shift;
    my %strain_abun;
    collect_abun($abun_dir, $samples, \%strain_abun, $kmer_size);
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
        foreach my $j (@{$samples}){
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
            if ($i =~ /o__Primates|o__Rodentia|o__Carnivora/){
                $host_abun->{'Host_reads'}->{$j} += $abun->{$i}->{$j};
                
            }else{
                $microbial_abun->{$i}->{$j} = $abun->{$i}->{$j}
            }
        }
    }
    
    foreach my $j (@{$samples}){
        if (! exists($host_abun->{'Host_reads'}->{$j})) {
            $host_abun->{'Host_reads'}->{$j}=0;
        }
        $host_abun->{'Microbial_reads'}->{$j} = 1-$host_abun->{'Host_reads'}->{$j};
    }
    foreach my $i (keys%{$microbial_abun}){
        foreach my $j (keys%{$microbial_abun->{$i}}){
            $microbial_abun->{$i}->{$j} = int($microbial_abun->{$i}->{$j}/$host_abun->{'Microbial_reads'}->{$j}*100000)/100000;
        }
    }
}
sub collect_abun{
    my $abun_dir = shift;
    my $samples = shift;
    my $strain_abun = shift;
    my $kmer_size = shift;
    foreach my $i (@{$samples}){
        my $abun_file = "$i/$i.k$kmer_size.sketch.with-lineages.csv";
        my %info;
        if (! -e $abun_file) {
            next;
        }
        
        my @strains = read_csv_table($abun_file, \%info, "name");
        my @reliable_strains;
        foreach my $i (@strains){
            my $identity = $info{$i}{'match_containment_ani'};
            my $intersect_bp = $info{$i}{'unique_intersect_bp'};
            if (($identity >= 0.935)||($intersect_bp>1000000)) {
                push(@reliable_strains, $i);
            }
        }
        @strains = @reliable_strains;
        my $sum = 0;
        foreach my $i (@strains){
            $sum+=$info{$i}{'f_unique_weighted'};
        }
        foreach my $j (@strains){
            my $ranks = $info{$j}{'lineage'};
            my $lineage = refine_lineage($ranks, $j);
            my $abun = $info{$j}{'f_unique_weighted'};
            if ($sum>0) {
                $strain_abun->{$lineage}{$i}=int($abun/$sum*100000)/100000;
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
    open(my $f1, "<$table") or die ("$table\n");
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
    my @row_order;
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
