#!/usr/bin/perl
#use strict;
#use warnings;
use vars qw($opt_s $opt_r $opt_l $opt_o $opt_t $opt_C $opt_T $opt_S);
&Getopts('s:r:l:o:C');
# -s raw read sequences
# -r reference genome folder
# -l initial list of hits
# -o output_dir

my $fastq_file = $opt_s;
my $ref_genome_dir = $opt_r;
my $sourmash_hit_file = $opt_l;
my $output_dir = $opt_o;
my $sample_id = $opt_t;
my $cap_size = $opt_C;
my $threads = $opt_T;
my $script_dir = $opt_S;

my $fastq = "$output_dir/$sample_id.fastq";
my $fastq_bk = "$output_dir/$sample_id.bk.fastq";
# need to rename the contig names to simplify classification of reads
rename_and_bk_fastq($fastq_file, $fastq, $fastq_bk);


my %sourmash_hit_info; # store genome size information
my %sourmash_hit_fate;
read_sourmash_hit_table($sourmash_hit_file, \%sourmash_hit_info);

my $group_bp_cap = 100;# 100Mb, each group of genomes has the maximum of total genome size of 100M bp
my @ref_genome_groups;# an array of genome_indexes of each group sorted by abundance
download_genomes($output_dir, \%genome_info, $ref_genome_dir);
group_and_index_genomes($ref_genome_dir, \%sourmash_hit_info, $group_bp_cap, @ref_genome_groups, $output_dir, $threads);


my %mapping_hit_size;
my $read_cutoff=$cap_size*0.000001;
if ($read_cutoff<10) {
    $read_cutoff=10;
}
my $SAM = "$output_dir/$sample_id.sam";
foreach my $genome_group_index (@ref_genome_groups){    
    map_reads($fastq, $genome_group_index, $SAM, $threads);
    #parse a subset of reads in the bam file
    my %unmapped_reads;
    my $subset_size = 200000;
    my ($hit_info, $read_info, $read_counts) = parse_sam_subset($SAM, $subset_size, \%unmapped_reads); # do not enrol read information into memory -> out of memory
    my $survived_hit_size = bam_hier_clustering($hit_info, $read_info, $read_cutoff); 
    my $left_hit_size=();
    my ($left_hit_info, $left_read_info);
    if ($read_counts == $subset_size) {
	# parse the sam file again, reduce the abundance of a hit by 1 if a assigned read was reassigned
	($left_hit_info, $left_read_info, $survived_hit_size) = substract_sam($SAM, $survived_hit_size, \%unmapped_reads); #update the survived hit size as well
	$left_hit_size = bam_hier_clustering($left_hit_info, $left_read_info, $read_cutoff);
    }
    save_hits(\%mapping_hit_size, $survived_hit_size, $left_hit_size);# keep checking what hits survived
    update_reads($output_dir, $fastq_bk, $fastq, \%unmapped_reads);
}

#redo mapping with curated hits to get accurate abundance and perform filtration by coverage.
my @survival =  keys%mapping_hit_size;
my $survival_genome_index = "$output_dir/survival_genome.idx";
format_ref_db(\@survival, $output_dir, $survival_genome_index);
map_reads($fastq_bk, $survival_genome_index, $SAM, $threads);
my ($hit_size, $unique_hit_size) = parse_sam_file($SAM, \%mapping_hit_size);
my ($NR_hit_size, $NR_hit_info) = divide_sam_file($SAM, $unique_hit_size, $hit_size);
system("echo step6 start: genome coverage filtering\n");
my $ave_read_length = 150; 
my $optimal_hit_size = filter_hits_by_genome_coverage($NR_hit_size, $output_dir, $genome_info, $ave_read_length, $lineage, $script_dir);
    
#summarize read fate into a table
system("echo step7 start: print results\n");
my %taxa2abun;
print_summary_results($optimal_hit_size, $phylogeny, $lineage, $NR_hit_info, \%taxa2abun, $genome_info, $output_dir);
print_hit_fate(\%sourmash_hit_fate, \%sourmash_hit_info, \%mapping_hit_size, $optimal_hit_size, $output_dir, $lineage, $genome_info, $phylogeny);

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