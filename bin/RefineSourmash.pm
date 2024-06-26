#!/usr/bin/perl
use strict;
use warnings;
#use Inputs;
use Storable 'dclone';
use Parallel::Loops;
use Text::CSV;
use Inputs;
package RefineSourmash;

#--------------------

sub TaxaPolisher{
    my $input_dir = $_[0];
    my $sourmash_dir = $_[1];
    my $output_dir = $_[2];
    my $sample_info = $_[3];
    my $bin_path = $_[4];
    my $threads = $_[5];
    my $database_folder = $_[6];
    my $phylogeny = $_[7];
    my $script_dir = $_[8];
    my $cap_size =$_[9];
    my @samples = sort(keys%{$sample_info});
    my %sample2taxa_abun_table;    
    #summarize sourmash results
    #foreach my $i (@samples){
    my $cpu = `nproc`;
    chomp $cpu;
    my $cpu_per_thread=1;
    if ($threads>=scalar(@samples)) {
	$threads = scalar(@samples);
    }
    
    $cpu_per_thread = sprintf("%0.f", $cpu/$threads);
    if ($cpu_per_thread == 0) {
	$cpu_per_thread = 1;
    }
    my %fastq_info;
    my $trim_summary = "$input_dir/summary.tsv";
    read_fastq_summary($trim_summary, \%fastq_info);
    
    my $pl = Parallel::Loops->new($threads);
    $pl -> share(\%sample2taxa_abun_table);
    #foreach my $i (@samples){
    $pl -> foreach (\@samples, sub{
        my $i = $_;
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1= "$input_dir/$i.fastq.gz";
	mkdir("$output_dir/$i");
	my $log = "$output_dir/$i/polish.log.txt";
	system("touch $log");
	if (! -e $R1) { # no readsq
	    Inputs::log_only("$R1 does not exists!", "$output_dir/$i/polish.log.txt");
	    my %taxa2abun=();
	    $sample2taxa_abun_table{$i}=\%taxa2abun;#in709_1_class_filtered
	}elsif(! -e "$sourmash_dir/$i/$i.merged.csv"){ # no reads mapped to reference genome
	    Inputs::log_only("$sourmash_dir/$i/$i.merged.csv does not exists!", "$output_dir/$i/polish.log.txt");
	    my %taxa2abun=();
	    $sample2taxa_abun_table{$i}=\%taxa2abun;
	}else{
	    my $ave_read_length = $fastq_info{$i}{'read_length(bp)'};
	    my $taxa2abun = refine_sourmash_results("$output_dir", "$sourmash_dir/$i/$i.merged.csv", $phylogeny, $R1, $database_folder, $cpu_per_thread, $i, $script_dir, $cap_size, $bin_path, $log, $ave_read_length);
	    $sample2taxa_abun_table{$i}=$taxa2abun;
	    cleanup_temp_files($input_dir, $sourmash_dir, $output_dir, $i);
	}
	
    });
    #}
}


sub read_fastq_summary{
    my $file = shift;
    my $fastq_info = shift;
    open(my $f1, "<$file") or die;
    my @titles;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^internal_id/) {
	    @titles = split(/\t/, $line);
	}else{
	    my @coln = split(/\t/, $line);
	    for (my $i=1; $i<scalar(@coln); $i++){
		$fastq_info->{$coln[0]}->{$titles[$i]}= $coln[$i];
	    }
	}
	
    }
    close $f1;
}



sub cleanup_temp_files{
    my $trim_dir = shift;
    my $sourmash_dir = shift;
    my $dir = shift;
    my $sample_id = shift;
    system("rm -r $dir/$sample_id/genom_cov/*.bam $dir/$sample_id/genom_cov/*.bai $dir/$sample_id/genom_cov/*.cov.txt");
    #system("rm -r $dir/$sample_id/genom_cov/*.bam $dir/$sample_id/genom_cov/*.bai");
    system("rm -r $dir/$sample_id/ref_genomes $dir/$sample_id/sam_files");
    system("rm -r $dir/$sample_id/results.sam");
    system("rm -r $sourmash_dir/$sample_id/*.sketch*.csv $sourmash_dir/$sample_id/*.sig.gz $dir/$sample_id/*.fastq");
    #system("rm -r $trim_dir/$sample_id.fastq.gz");
}

sub refine_sourmash_results{
    my $output_dir=shift;
    my $hit_table=shift;
    my $phylogeny = shift;
    my $R1 = shift;
    my $database = shift;
    my $cpu_per_thread = shift;
    my $sample_id = shift;
    my $script_dir = shift;
    my $cap_size = shift;
    my $bin_path = shift;
    my $log=shift;
    my $ave_read_length = shift;
    my $hit_length_cutoff= 80; # hits less than 80 bp match 

    #system("mkdir $output_dir/$sample_id");
    
    $output_dir = "$output_dir/$sample_id";
    my $new_R1 = 'R1.fastq';
    my $bk_R1 = 'R1_bk.fastq';
    system("gunzip -c $R1 > $output_dir/$new_R1");
    system("cp $output_dir/$new_R1 $output_dir/$bk_R1");
    my $subset_size = 10000; #200000 default
    my $genome_bin_size = 500; # 100Mb genome bp
    my %mapping_hier_hit_size=();
    my %unmapped_reads=();
    my @sourmash_hit_no=();
    my $read_cutoff=$cap_size*0.000001;
    if ($read_cutoff<10) {
	$read_cutoff=10;
    }
    my $genome_info = get_genome_info("$database/genome_info.tsv");
    my %sourmash_hit_size;
    my @headers = read_sourmash_hit_table("$hit_table", \%sourmash_hit_size, $genome_info, 1);
    if (scalar(keys(%sourmash_hit_size))==0) {
	Inputs::log_only("No useful sourmash hits found for sample of $sample_id!", $log);
	return ();
    }
    my $s3_dir = "s3://zymo-files/WGS_Pipeline/shotgun_database/shotgun_genomes_240521";
    my @ref_genome_bins;# an array of genome_indexes of each group sorted by abundance
    mkdir("$output_dir/ref_genomes");
    Inputs::log_only("Download and index genome bins", $log);
    download_genomes($output_dir, \%sourmash_hit_size, $genome_info, $s3_dir, $log, $cpu_per_thread);
    group_and_format_genomes($output_dir, \%sourmash_hit_size, $genome_bin_size, \@ref_genome_bins, $cpu_per_thread, $genome_info);
    
    foreach my $genome_index_per_bin (@ref_genome_bins) {
	Inputs::log_only("Processing $genome_index_per_bin", $log);
	read_mapping("$output_dir/$new_R1", $genome_index_per_bin, $output_dir, $cpu_per_thread, $hit_length_cutoff);
	my ($hit_info, $read_info, $reads) = parse_sam_subset($output_dir, $subset_size, \%unmapped_reads); # do not enrol read information into memory -> out of memory
	my $survived_hit_size = mapping_hier_clustering($hit_info, $read_info, $read_cutoff); # hirachical clustering to keep dominant hits, 
	my $left_hit_size=();
	my ($left_hit_info, $left_read_info);
	if ($reads == $subset_size) {
	    #update the survived hit size as well, the hits that are removed from the subset-screening will be enrolled and evaluated thoroughly agains
	    # the survived hits will have their size updated again.
	    ($left_hit_info, $left_read_info, $survived_hit_size) = substract_sam($output_dir, $survived_hit_size, \%unmapped_reads);
	    $left_hit_size = mapping_hier_clustering($left_hit_info, $left_read_info, $read_cutoff);
	}
	#combine all survived hits.
	save_hits(\%mapping_hier_hit_size, $survived_hit_size, $left_hit_size);
	#update $new_R1 to keep unmapped reads only
	update_reads($output_dir, $bk_R1, $new_R1, \%unmapped_reads);
	#empty unmapped reads
	%unmapped_reads=();
    }
    
    #redo mapping with curated hits to get accurate abundance and perform filtration by coverage.
    my $unmapped_read_no = scalar(keys%unmapped_reads);
    my @survival =  keys%mapping_hier_hit_size;
    my $ref_genome_index = 'all.fasta';
    Inputs::log_only("index and mapping surviving genomes", $log);
    format_ref_db(\@survival, $output_dir, $ref_genome_index);
    read_mapping("$output_dir/$bk_R1", $ref_genome_index, $output_dir, $cpu_per_thread, $hit_length_cutoff);
    my ($hit_size, $unique_hit_size) = parse_sam_file($output_dir, \%mapping_hier_hit_size);
    my ($NR_hit_size, $NR_hit_info) = divide_sam_file($output_dir, $unique_hit_size, $hit_size);
    system("echo step6 start: genome coverage filtering\n");
    Inputs::log_only("filter hits by genome coverage", $log);
    my $optimal_hit_size = filter_hits_by_genome_coverage($NR_hit_size, $output_dir, $genome_info, $ave_read_length, $script_dir, $cpu_per_thread);
    
    #summarize read fate into a table
    system("echo step7 start: print results\n");
    my $new_sourmash_table = "$output_dir/$sample_id.after.mapping.csv";
    Inputs::log_only("print taxa polish results", $log);
    print_new_sourmash_table(\@headers, $optimal_hit_size, $genome_info, $new_sourmash_table);
    my %taxa2abun;
    print_summary_results($optimal_hit_size, $phylogeny, $NR_hit_info, \%taxa2abun, $genome_info, $output_dir);
    print_hit_fate(\%sourmash_hit_size, \%mapping_hier_hit_size, $optimal_hit_size, $output_dir, $genome_info, $phylogeny);
    
}

sub print_new_sourmash_table{
    my $headers = shift;
    my $optimal_hit_size = shift;
    my $genome_info = shift;
    my $new_sourmash_table = shift;
    if (scalar(keys(%{$optimal_hit_size}))==0) {
	return;
    }
    open(my $f1, ">$new_sourmash_table") or die;
    $headers->[0]= 'hits';
    $headers->[1]= 'read_counts';
    my $header = join(",", @{$headers});
    print($f1 "$header\n");
    my @hits = sort{$optimal_hit_size->{$b} <=> $optimal_hit_size->{$a}}keys(%{$optimal_hit_size});
    foreach my $i (@hits){
	my $read_counts = $optimal_hit_size->{$i};
	my @coln = ($i, $read_counts);
	for(my $j=2; $j<scalar(@{$headers}); $j++){
	    my $value = $genome_info->{$i}->{$headers->[$j]};
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

#group_and_format_genomes($ref_dir, \%sourmash_hit_size, $genome_bin_size, \@ref_genome_bins, $cpu_per_thread)


sub group_and_format_genomes{
    my $dir = shift;
    my $sourmash_hit_size = shift;
    my $genome_bin_size = shift;
    my $ref_genome_bins = shift;
    my $threads = shift;
    my $genome_info = shift;
    my @genome_order = sort{$sourmash_hit_size->{$b} <=> $sourmash_hit_size->{$a}}keys(%{$sourmash_hit_size});
    if (scalar(@genome_order)==0) {
	return;
    }
    
    my $total_genome_size = 0;
    my @genome_bins;
    my @genomes;
    foreach my $i (@genome_order){
	if (! exists($genome_info->{$i})) {
	    next;
	}
	my $size = $genome_info->{$i}->{'genome_size'};
	$total_genome_size+=$size;
	push(@genomes, "$i");
	if ($total_genome_size>$genome_bin_size*1000000) {
	    my @genomes_to_enroll = @genomes;
	    push(@genome_bins, \@genomes_to_enroll);
	    $total_genome_size = 0;
	    @genomes = ();
	}
	
    }
    if (scalar(@genomes)>0) { # collect the last batch if exists.
	push(@genome_bins, \@genomes);
    }
    my $pl = Parallel::Loops->new($threads);
    my $bins = scalar(@genome_bins);
    my @index = (0..$bins-1);
    #foreach my $i (@index){
    $pl -> foreach (\@index, sub{
        my $i = $_;
	my $genome_index = "bin$i.fasta";
	format_ref_db($genome_bins[$i], $dir, $genome_index);
    });
    #}
    for (my $i =0; $i<scalar(@genome_bins); $i++){
	my $genome_index = "bin$i.fasta";
	push(@{$ref_genome_bins}, $genome_index);
    }
}


#download_genomes($ref_dir, \%sourmash_hit_size, $genome_info, $s3_dir);
sub download_genomes{
    my $output_dir = shift;
    my $sourmash_hit_size = shift;
    my $genome_info = shift;
    my $s3_dir = shift;
    my $log = shift;
    my $threads=shift;
    my @genomes = sort(keys(%{$sourmash_hit_size}));
    if ($threads>5) {
	$threads = 5;
    }
    
    my $pl = Parallel::Loops->new($threads);

    #foreach my $i (@genomes){
    $pl -> foreach (\@genomes, sub{
        my $i = $_;
	if (exists($genome_info->{$i})) {
	    my $path = $genome_info->{$i}->{'path'};
	    my $s3_path = "$s3_dir/$path";
	    my $cmd = "aws s3 cp $s3_path $output_dir/ref_genomes/$i.fna.gz";
	    Inputs::print_and_execute($cmd, $log);
	    concatenate_contigs("$output_dir/ref_genomes", $i, $genome_info);
	}else{
	    Inputs::log_only("The genome of $i was not found", $log);
	}
    });
    #}
}

sub concatenate_contigs{
    my $dir = shift;
    my $id = shift;
    my $genome_info = shift;
    open(my $f1, " gunzip -c $dir/$id.fna.gz |") or die;
    open(my $f2, ">$dir/$id.fasta") or die;
    print($f2 ">$id\n");
    my $polyN='';
    while (my $line = <$f1>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		if (length($polyN)==0){
		    $polyN = 'N' x 100;
		}else{
		    print($f2 "$polyN\n");
		    $genome_info->{$id}->{'genome_size'}+=100;
		}
	    }else{
		print($f2 "$line\n");
	    }
    }
    close $f1;
    close $f2;
    
}

#read_sourmash_hit_table
sub read_sourmash_hit_table{
    my $hit_table = shift;
    my $hit_size = shift;
    my $hit_info=shift;
    my $weight = shift;
    open(my $f1, "<$hit_table") or return;
    my $title = <$f1>;chomp $title;
    my $csv = Text::CSV->new();
    $csv->parse($title);
    my @title = $csv->fields();
    my $index_ID=0;
    for(my $i=0; $i<scalar(@title); $i++){
	if ($title[$i] eq 'name') {
	    $index_ID=$i;
	}
	
    }
    my %hits_of_interest;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line=~/^Y/) {
	    $csv->parse($line);
	    my @coln = $csv->fields();
	    my $name = $coln[$index_ID];
	    my @id = split(/\ /, $name);
	    if (exists($hit_info->{$id[0]}->{'genome_size'})) {# exclude genomes that are too big
		if ($hit_info->{$id[0]}->{'genome_size'}>100000000) {
		    next;
		}
	    }
	    for(my $i = 0; $i<scalar(@title); $i++){
		$hit_info->{$id[0]}->{$title[$i]}=$coln[$i];
	    }
	    $hit_info->{$id[0]}->{'f_unique_weighted'}=($hit_info->{$id[0]}->{'f_unique_weighted'})*$weight;
	    $hit_size->{$id[0]}=$hit_info->{$id[0]}->{'f_unique_weighted'};
	    $hits_of_interest{$id[0]}++;
	    my $lineage = $hit_info->{$id[0]}->{'lineage'};
	    my $species = extract_species_name($lineage);
	    $hit_info->{$id[0]}->{'species_name'} = $species;
	}
    }
    foreach my $i (keys(%{$hit_info})){
	if (!exists($hits_of_interest{$i})) {
	    delete($hit_info->{$i});
	}
	
    }
    close $f1;
    return(@title);
}

sub update_hits{
    my $hit2reads=$_[0];
    my $hits_to_skip=$_[1];
    foreach my $i (@{$hits_to_skip}){
	delete($hit2reads->{$i});
    }
}

sub get_top_hits{
    my $hit_size = $_[0];
    my $hit_limit = $_[1];
    my @hits = sort{$hit_size->{$b} <=> $hit_size->{$a}} keys(%{$hit_size});
    if (scalar@hits>$hit_limit) {
	@hits = @hits[0..$hit_limit-1];
    }
    return(\@hits);  
}

sub copy_and_filter_host{
    my $R1=$_[0];
    my $new_R1=$_[1];
    my $host_reads = $_[2];
    my $dir = $_[3];
    my $cap_size = $_[4];
    my $bin_path = $_[5];
    my @read_order = sort{$a <=> $b} keys(%{$host_reads});
    open(my $host, ">$dir/host_reads.txt") or die;
    my $host_read_count = scalar(@read_order);
    my $host_read_str = join(",", @read_order);
    print($host ">host_reads $host_read_count\n$host_read_str\n");
    open(my $f1, "<$R1") or die;
    open(my $f3, ">$dir/$new_R1") or die;
    open(my $f5, ">$dir/R1.fasta") or die;
    my $counter=0;
    my $added=0;
    my $bp =0;
    while (my $r1 = <$f1>) {
	my $r2 = <$f1>; my $r3 = <$f1>; my $r4 = <$f1>;
	if (scalar(@read_order)>0) {
	    if ($counter==$read_order[0]) { #skip host reads
		shift(@read_order); 
	    }else{
		print($f3 "$r1$r2$r3$r4");
		$r1 = substr($r1, 1, length($r1)-1);
		print($f5 ">$r1$r2");
		$added++;		
	    }
	}else{
	    print($f3 "$r1$r2$r3$r4");
	    $r1 = substr($r1, 1, length($r1)-1);
	    print($f5 ">$r1$r2");
	    $added++;
		
	}
	$counter++;
	$bp+=length($r2)-1;
	
	if ($added == $cap_size) { # keep a cap size of the reads 
	    last;
	}
    }
    close $f1;
    close $f3;
    close $f5;
    
    system("$bin_path/sdust/sdust $dir/R1.fasta > $dir/R1_low_complexity_reads.txt");
    my %low_complex_reads;
    read_sdust_output("$dir/R1_low_complexity_reads.txt", \%low_complex_reads);
    
    my @low_complex_reads = sort {$a <=> $b} keys(%low_complex_reads);
    open(my $f0, ">$dir/low_complexity_reads.txt") or die;
    my $content = join("\n", @low_complex_reads);
    print($f0 "$content\n");
    close $f0;
    system("rm $dir/R1.fasta $dir/R1_low_complexity_reads.txt");
    
    open(my $F1, "<$dir/$new_R1") or die;
    open(my $F3, ">$dir/R1.N.fastq") or die;
    while (my $r1 = <$F1>) {
	my $r2 = <$F1>; my $r3 = <$F1>; my $r4 = <$F1>;
	
	if (scalar(@low_complex_reads)>0) {
	    my $query = $low_complex_reads[0];
	    if ($r1=~/^\@$query /) { #skip low diversity reads
		shift(@low_complex_reads);
	    }else{
		print($F3 "$r1$r2$r3$r4");	
	    }
	}else{
	    print($F3 "$r1$r2$r3$r4");
	}
    }
    close $F1;
    close $F3;
    system("mv $dir/R1.N.fastq $dir/$new_R1");
    
    my $ave_read_length = int($bp/$counter/2);
    return($ave_read_length);
    
}

sub read_sdust_output{
    my $file = shift;
    my $reads = shift;
    open(my $f1, "<$file") or die;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (($coln[2]-$coln[1])>30) {
	    $reads->{$coln[0]}++;
	}
    }
    close $f1;
    
}

sub not_low_complexity{
    my $seq1 = shift;
    my $seq2 = shift;
    open(my $f1, ">tmp.fasta") or die;
    print($f1 ">1\n$seq1\n>2\n$seq2\n");
    close $f1;
    my $complexity = `sdust/sdust tmp.fasta`;
    if ($complexity =~ /[12]/) {
	return (0);# low complexity
    }else{
	return (1);# not low complexity
    }
    
}


# print_hit_fate(\%sourmash_hit_size, \%mapping_hier_hit_size, $optimal_hit_size, $output_dir, $genome_info, $phylogeny);
sub print_hit_fate{
    my $sourmash_hit_fate=shift;
    my $mapping_hier_hit_size = shift;
    my $optimal_hit_size = shift;
    my $dir = shift;
    my $genome_info=shift;
    my $phylogeny=shift;
    open(my $f1, ">$dir/hit_fate.tsv") or die;
    print($f1 "genome_id\ttaxid\ttaxa\tsourmash_abun\tmapping_reads\tlast_refinery_reads\tseq_identity(%)\tgenome_coverage(%)\tgenome_ave_depth\tread_depth(>1x)\ttotal_mapped_reads\tbp_distribution_ratio\n");
    my @hits = sort {$sourmash_hit_fate->{$b} cmp $sourmash_hit_fate->{$a}} keys%{$sourmash_hit_fate};
    foreach my $i (@hits){
	my $lineage = $genome_info->{$i}->{'lineage'};
	my $cent_fate = $sourmash_hit_fate ->{$i};
	my $mapping_fate ='';
	my $opt_fate='';
	if (exists($mapping_hier_hit_size->{$i})) {
	        $mapping_fate = "$mapping_hier_hit_size->{$i}";
		if (exists($optimal_hit_size->{$i})) {
		    $opt_fate = "$optimal_hit_size->{$i}";
		}else{
		    $opt_fate = 'removed';
		}
		
	}else{
		$mapping_fate = 'removed'    
	}
	if (-e "$dir/genom_cov/$i.cov.summary.txt") {
	    open(my $f2, "<$dir/genom_cov/$i.cov.summary.txt") or die;
	    my %hit_mapping_info;
	    while (my $line = <$f2>) {
		chomp $line;
		if ($line =~ /\:/) {
		    my @coln = split(/\:/, $line);
		    $hit_mapping_info{$coln[0]}=$coln[1];
		}
	    }
	    my $seq_iden = $hit_mapping_info{'seq_identity(%)'};
	    my $genome_cov = $hit_mapping_info{'genome_coverage(%)'};
	    my $genome_ave_depth = $hit_mapping_info{'genome_ave_depth'};
	    my $read_depth = $hit_mapping_info{'read_depth'};
	    my $total_mapped_reads = $hit_mapping_info{'total_mapped_reads'};
	    my $bp_dist_ratio = sprintf("%.4f", $hit_mapping_info{'distribution_ratio'});
	    
	    print ($f1 "$i\tNA\t$lineage\t$cent_fate\t$mapping_fate\t$opt_fate\t$seq_iden\t$genome_cov\t$genome_ave_depth\t$read_depth\t$total_mapped_reads\t$bp_dist_ratio\n");
	}else{
	    print ($f1 "$i\tNA\t$lineage\t$cent_fate\t$mapping_fate\t$opt_fate\t\t\t\t\t\t\n");
	}
	
	
	
    }
    close $f1;
}

#print_summary_results($optimal_hit_size, $phylogeny, $NR_hit_info, \%taxa2abun, $genome_info, $output_dir);
sub print_summary_results{
    my $optimal_hit_size=$_[0];
    my $phylogeny = $_[1];
    my $NR_hit_info = $_[2];
    my $taxa2abun = $_[3];
    my $genome_info = $_[4];
    my $dir = $_[5];
    open(my $f1, ">$dir/hit2read.txt") or die;
    open(my $f3, ">$dir/hit2tax.txt") or die;
    print($f3 "hits\ttaxonomy\tread_counts\n");
    open(my $f2, ">$dir/read2hit.tsv") or die;
    print($f2 "read_id\tgenome_id\tspecies\tlineage\n");
    foreach my $i (keys%{$optimal_hit_size}){
	my $lineage = $genome_info->{$i}->{'lineage'};
	my @ranks = split(/\;/, $lineage);
	my $species = $ranks[6];
	if ($species =~ /__/){
	    $species = substr($species, 3, length($species)-3);
	}
	my @reads = @{$NR_hit_info->{$i}};
	my $read_str = join(",", @reads);
        my $read_count = scalar(@reads);
	print($f1 ">$i\_$lineage\n$read_str\n");
        print($f3 "$i\t$lineage\t$read_count\n");
	foreach my $j (@reads){
	    print($f2 "$j\t$i\t$species\t$lineage\n");
	}
	$taxa2abun->{$lineage} += $optimal_hit_size->{$i};
    }
    close $f1;
    close $f2;
    close $f3;
}


sub format_ref_db {
    my $hit_list=$_[0];
    my $dir = $_[1];
    my $genome_index= $_[2];
    my @genomes='';
    if (-e "$dir/ref_genomes/$genome_index") {
	system("rm $dir/ref_genomes/$genome_index");
    }    
    system("touch $dir/ref_genomes/$genome_index");
    foreach my $i (@{$hit_list}){
	my $file = "$dir/ref_genomes/$i.fasta";
	if (-e $file) {
	    system("cat $file >>$dir/ref_genomes/$genome_index");
	}
    }
    my $cmd = "bwa index -a bwtsw $dir/ref_genomes/$genome_index > /dev/null 2>&1";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt");
}

sub read_mapping{
    my $R1 = shift;
    my $genome_index = shift;
    my $dir = shift;
    my $threads = shift;
    my $hit_length_cutoff = shift;
    my $cpu = `nproc`;
    chomp $cpu;
    $threads = $cpu;
    my $cmd = "bwa mem -a -t $threads $dir/ref_genomes/$genome_index $R1 -o $dir/results1.sam > /dev/null 2>&1";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt");
    #remove the hits of short length from the sam file.
    open(my $f1, "<$dir/results1.sam") or die;
    open(my $f2, ">$dir/results.sam") or die;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~/^@/) {
	    my @coln = split(/\t/, $line);
	    my $cigar = $coln[5];
	    my @pattern = $cigar =~ m/(\d+)[MI]/g;
	    my $size = 0;
	    foreach my $i (@pattern){
		$size+=$i;
	    }
	    if ($size>0 && $size<$hit_length_cutoff) {
		next;
	    }
	}
	
	print($f2 "$line\n");
	
    }
    close $f1;
    close $f2;
    system("rm $dir/results1.sam");

}

sub parse_sam_subset{
    my $dir=$_[0];
    my $subset_size = $_[1];
    my $unmapped_reads = $_[2];
    open(my $f1, "<$dir/results.sam") or die;
    my %hit_info;
    my %read_info;
    my $last_read='';
    my @records;
    my $counter=0;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
		if (scalar(@records)>0) {
		    if ($coln[0] eq $last_read) {
			push(@records, $line);
		    }else{
			my $no_hit = 0;
			foreach my $line (@records){
			    my @coln = split(/\t/, $line);
			    if ($coln[2] ne '*') {
				$hit_info{$coln[2]}{$coln[0]}++;
				$read_info{$coln[0]}{$coln[2]}++;
				$no_hit=1;
			    }
			}
			if ($no_hit == 0) {
			    $unmapped_reads->{$last_read}++;
			}
			
			@records = ();
			push(@records, $line);
			$last_read = $coln[0];
			$counter++;
			if ($counter == $subset_size) {
			    last;
			}
			
		    }
		}else{
		    push(@records, $line);
		    $last_read = $coln[0];
		}
	}
	
    }
    close $f1;
    
    return(\%hit_info, \%read_info, $counter);
    
}

#divide_sam_file($output_dir, $unique_hit_size);
sub divide_sam_file{
    my $dir=shift;
    my $unique_hit_size = shift;
    my $hit_size = shift;
    open(my $f1, "<$dir/results.sam") or die;
    my %hit_info;
    my %read_info;
    my $last_read='';
    my %hits;
    my %NR_hit_size;
    my %NR_hit_info;
    system("mkdir $dir/sam_files");
    my @hits = keys%{$hit_size};
    my %sam_fhs;
    foreach my $i (@hits){
	open($sam_fhs{$i}, ">$dir/sam_files/$i.sam") or die;
    }
    
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
	    if ($coln[2] ne '*') {
		if (scalar(keys%hits)>0) {
		    if ($coln[0] eq $last_read) {
			push(@{$hits{$coln[2]}}, $line);
		    }else{
			my $assigned_hit = assign_to_dominant($unique_hit_size, \%hits);
			my $size = scalar(@{$hits{$assigned_hit}});
			my @reads;
			if ($size>1) {
			    @reads = resolve_repeats(\@{$hits{$assigned_hit}});
			}else{
			    @reads = @{$hits{$assigned_hit}};
			}
			$size = scalar(@reads);
			$NR_hit_size{$assigned_hit}+=$size;
			push(@{$NR_hit_info{$assigned_hit}}, $last_read);
			my $str = join("\n", @reads);
			if (exists($sam_fhs{$assigned_hit})) {
			    my $fh = $sam_fhs{$assigned_hit};
			    print($fh "$str\n");
			}
			%hits = ();
			push(@{$hits{$coln[2]}}, $line);
			$last_read = $coln[0];
			
		    }
		}else{
		    push(@{$hits{$coln[2]}}, $line);
		    $last_read = $coln[0];
		}
		
	    }
	}else{
	    my $found = 0;
	    foreach my $i (keys%sam_fhs){
		if ($line =~ /SN:$i/) {
		    my $fh = $sam_fhs{$i};
		    print($fh "$line\n");
		    $found =1;
		}
	    }
	    if ($found == 0) {
		foreach my $i (keys%sam_fhs){
		    my $fh = $sam_fhs{$i};
		    print($fh "$line\n");
		}
	    }
	    
	}
	
    }
    close $f1;
    
    foreach my $i (keys%sam_fhs){
	close $sam_fhs{$i};
    }
    
    return(\%NR_hit_size, \%NR_hit_info);
    
}


sub mapping_hier_clustering{
    my $hit_info = $_[0];
    my $read_info = $_[1];
    my $reads_cutoff= $_[2];
    my %opt_hit_info;
    my $rank = 1;
    while (scalar(keys(%{$hit_info}))>0) {
	my $highest = get_highest_hit($hit_info);
	$rank++;
	my %reads = %{$hit_info->{$highest}};
	$opt_hit_info{$highest} = $hit_info->{$highest};
	delete($hit_info->{$highest});
	foreach my $j (keys%reads){
	    my %hits = %{$read_info->{$j}};
	    foreach my $i (keys%hits){
		if ($i ne $highest) {
		    delete($hit_info->{$i}->{$j});
		}
	    }
	    delete($read_info->{$j});
	}
	foreach my $x (keys%{$hit_info}){
	    my %reads = %{$hit_info->{$x}};	    
	    if (scalar(keys%reads)<$reads_cutoff) {
		delete($hit_info->{$x});
	    }
	    
	}
	
    }
    my %opt_hit_size;
    foreach my $i (keys%opt_hit_info){
	my %reads = %{$opt_hit_info{$i}};
	$opt_hit_size{$i} = scalar(keys%reads);
    }
    
    return (\%opt_hit_size);
}

sub substract_sam{ #$output_dir, $survived_hit_size, $subset_size, $cap_size, %unmapped_reads
    my $dir=$_[0];
    my $survived_hits = $_[1];
    my $unmapped_reads = $_[2];
    my %new_survived_hits;
    my %hit_info;
    my %read_info;
    open(my $f1, "<$dir/results.sam") or die;
    my $last_read='';
    my @records;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
		if (scalar(@records)>0) {
		    if ($coln[0] eq $last_read) {
			push(@records, $line);
		    }else{		
			my $chosen = 0;
			my $highest_hit='';
			my $highest_abun = 0;
			foreach my $line (@records){
			    my @coln = split(/\t/, $line);
			    if (exists($survived_hits->{$coln[2]})) {
				$chosen = 1;
				if ($survived_hits->{$coln[2]}>=$highest_abun) {
				    $highest_abun = $survived_hits->{$coln[2]};
				    $highest_hit = $coln[2];
				}
				
			    }
			    
			}
			my $read_mapped = '0';
			if ($chosen == 0) {
			    foreach my $record (@records){
				my @old_coln = split(/\t/, $record);
				if ($old_coln[2] eq '*') {
				    #code
				}else{
				    
				    $hit_info{$old_coln[2]}{$old_coln[0]}++;
				    $read_info{$old_coln[0]}{$old_coln[2]}++;
				    $read_mapped = 1;
				}
				
			        
			    }
			    if ($read_mapped == 0) { # save a list of unmapped reads
			        $unmapped_reads -> {$last_read}++;
			    }
			}else{
			    $new_survived_hits{$highest_hit}++;
			}
			
			
			
			@records = ();
			push(@records, $line);
			$last_read = $coln[0];
			
		    }
		}else{
		    push(@records, $line);
		    $last_read = $coln[0];
		}
		
	}
	
    }
    close $f1; # update the read counts for the pre-selected ones.
    
    return(\%hit_info, \%read_info, \%new_survived_hits);
    
}


# save_hits(\%hier_hit_size, $survived_hit_size, $left_hit_size);
sub save_hits {
    my $hier_hit_size = $_[0];
    my $survived_hit_size = $_[1];
    my $left_hit_siza = $_[2];
    foreach my $i (keys(%{$survived_hit_size})){
	$hier_hit_size->{$i} = $survived_hit_size ->{$i};
    }
    foreach my $i (keys(%{$left_hit_siza})){
	$hier_hit_size->{$i} = $left_hit_siza ->{$i};
    }
}

# update_reads($output_dir, $new_R1, $new_R2, %unmapped_reads);
sub update_reads{
    my $dir = $_[0];
    my $bk_R1 = $_[1];
    my $R1 =$_[2];
    my $unmapped_reads = $_[3];
    my @read_order = sort{$a <=> $b} keys(%{$unmapped_reads});
    open(my $f1, "<$dir/$bk_R1") or die;
    open(my $f3, ">$dir/$R1") or die;
    my $counter=0;
    while (my $r1 = <$f1>) {
	my $r2 = <$f1>; my $r3 = <$f1>; my $r4 = <$f1>;
	if ($counter == $read_order[0]) {
	    print($f3 "$r1$r2$r3$r4");
	    shift(@read_order);
	    if (scalar(@read_order)==0) {
	        last;
	    }
	}
	$counter++;
	
    }
    close $f1;
    close $f3;    
}

sub get_highest_hit{
    my $hit_info = shift;
    my %hit_size;
    foreach my $i (keys%{$hit_info}){
	my %reads = %{$hit_info->{$i}};
	$hit_size {$i} = scalar(keys%reads);
    }
    my @order = sort {$hit_size{$a} <=> $hit_size{$b}} keys%hit_size;
    my $highest = pop@order;
    return $highest;
}


sub parse_sam_file{
    my $dir = shift;
    my $selected_hits = shift;
    open(my $f1, "<$dir/results.sam") or die;
    my %hit_size;
    my %unique_hit_size;
    my $last_read='';
    my @records;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
	    if ($coln[2] ne '*') {
		if (scalar(@records)>0) {
		    if ($coln[0] eq $last_read) {
			push(@records, $line);
		    }else{
			add_and_count(\@records, \%hit_size, \%unique_hit_size, $selected_hits);
			@records = ();
			push(@records, $line);
			$last_read = $coln[0];
		    }
		}else{
		    push(@records, $line);
		    $last_read = $coln[0];
		}
		
	    }
	}
	
    }
    close $f1;
    
    return(\%hit_size, \%unique_hit_size);
    
}

sub add_and_count{
    my $records = $_[0];
    my $hit_size = $_[1];
    my $unique_hit_size = $_[2];
    my $selected_hits= $_[3];
    my %hits;
    foreach my $line (@{$records}){
	my @coln = split(/\t/, $line);
	if ($coln[2] !~ /^[A-Z]+[_0-9.]+$|^zymo/) {
	    return;
	}else{
	    if (exists($selected_hits->{$coln[2]})) {
		$hits{$coln[2]}++;
	    }
	}
	
    }
    
    foreach my $i (keys%hits){
	my $counts = $hits{$i};
	if ($counts>2) {
	    $hit_size->{$i}+=2;
	}else{
	    $hit_size->{$i}+=$counts;
	}
    }
    if (scalar(keys%hits)==1) {
	my @list = keys(%hits);
	my $counts = $hits{$list[0]};
	if ($counts>2) {
	    $unique_hit_size->{$list[0]}+=2;
	}else{
	    $unique_hit_size->{$list[0]}+=$counts;
	}
    }
    
}


sub print_sourmash_clustering_results{
    my $dir = shift;
    my $fate = shift;
    open(my $f1, ">$dir/sourmash_fate.txt") or die;
    foreach my $i (keys%{$fate}){
	my $note = $fate->{$i};
	print($f1 "$i\t$note\n");
    }
    close $f1;
}


sub assign_by_abun{
    my $hit_abun = $_[0];
    my $hits = $_[1];
    my $sum=0;
    my @hits = keys%{$hits};
    my %hit_sum;
    if (scalar(@hits)==1) {
	return($hits[0]);
    }
    
    foreach my $i (@hits){
	my $abun=1;
	if (exists($hit_abun->{$i})) {
	    $abun = $hit_abun->{$i};
	}
	$sum+=$abun;
	$hit_sum{$i} = $sum;
    }
    
    my $random_pick = rand($sum);
    foreach my $i (@hits){
	if ($random_pick<=$hit_sum{$i}) {
	    return ($i);
	}
	
    }
    
}

sub assign_to_dominant{
    my $unique_hit_abun = $_[0];
    my $hits = $_[1];
    my $sum=0;
    my @hits = keys%{$hits};
    if (scalar(@hits)==1) {
	return($hits[0]);
    }
    my @new_hits;
    foreach my $i (keys%{$hits}){
	if (exists($unique_hit_abun->{$i})) {
	   push(@new_hits, $i);
	}
    }
    @new_hits = sort{$unique_hit_abun->{$b}<=>$unique_hit_abun->{$a}}@new_hits;
    return($hits[0]);
    
}

sub resolve_repeats{ 
    # this function randomly assign the read to one position
    my $records = shift;
    my @fwd;
    foreach my $i (@{$records}){
	push(@fwd, $i);
    }   
    my @reads;
    if (scalar(@fwd)>0) {
	my $index = rand(@fwd);
	my $fwd_choice = $fwd[$index];
	push(@reads, $fwd_choice);
    }
    return(@reads);
}


sub create_hit_samfile{
    my $hits = shift;
    my $dir = shift;
    my %fhs;
    foreach my $i (@{$hits}){
	open(my $f1, ">$dir/sam_files/$i.sam") or die;
	$fhs{$i} = $f1;
    }
    return(\%fhs);
}


sub get_hit_size{
    my $hit2reads=$_[0];
    my %hit_size;
    foreach my $i (keys(%{$hit2reads})){
	my @reads = keys(%{$hit2reads->{$i}});
	$hit_size{$i}+=scalar(@reads);
    }
    return (\%hit_size);
}



# my ($remain_read_info, $remain_hit_info, $sourmash_hier_hit_info_subset) = substract_hits_and_reads($sourmash_read_info, $sourmash_hier_hit_info_subset);
sub substract_hits_and_reads{
    my $read_info=$_[0];
    my $hit2remove_size = $_[1];
    my @hits = keys%{$hit2remove_size};
    @hits = sort{$hit2remove_size->{$b}<=>$hit2remove_size->{$a}} @hits;
    my %hit_info_subset;
    my %hit_info_remain;
    
    foreach my $i (keys%{$read_info}){
	my $hits = $read_info->{$i};
	foreach my $j (@hits){
	    if (exists($hits->{$j})) {
		delete($read_info->{$i});
		$hit_info_subset{$j}{$i}++;
		last;
	    }
	    
	}
	
    }
    
    foreach my $i (keys%{$read_info}){
	my $hits = $read_info->{$i};
	foreach my $j (keys%{$hits}){
	    $hit_info_remain{$j}{$i}++;
	}
    }
    return($read_info, \%hit_info_remain, \%hit_info_subset);
}


sub taxid2tax{
    my $taxid = shift;
    my $phylogeny = shift;
    my $lineage = shift;
    my @ranks = @{$lineage->{$taxid}};
    my @rank_index=@{$phylogeny -> {'ranks_to_use'}};
    my @new_ranks;
    foreach my $i (@rank_index){
	push (@new_ranks, $ranks[$i-1]);
    }
    return (\@new_ranks);
}


sub filter_hits_by_genome_coverage{
    my $NR_hit_size = $_[0];
    my $dir = $_[1];
    my $genome_info = $_[2];
    my $read_size = $_[3];
    my $script_dir = $_[4];
    my $threads = $_[5];
    my %optimal_hit_size;
    my @genomes = sort(keys%{$NR_hit_size});
    my $cpu = `nproc`;
    chomp $cpu;
    my $pl = Parallel::Loops->new($cpu);
    $pl -> share(\%optimal_hit_size);
    $pl -> foreach (\@genomes, sub{
	my $i=$_;
    #foreach my $i (@genomes){
	my $sam_file = "$dir/sam_files/$i.sam";
	my $ref_genome = "$dir/ref_genomes/$i.fasta";
	my ($tot_mapped_reads, $cov_table) = sam2cov($dir, $i, $sam_file, $ref_genome, $script_dir, $read_size, $genome_info);
	if (removed_by_genome_cov($dir, $cov_table, $ref_genome, $tot_mapped_reads, $genome_info, $i, $read_size)) {
	}else{
	    $optimal_hit_size{$i} = $NR_hit_size->{$i};
	}
	
    #}
    });
    return(\%optimal_hit_size);
}

sub sam2cov{
    my $dir=shift;
    my $assemblyid = shift;
    my $sam_file = shift;
    my $ref_genome = shift;
    my $script_dir = shift;
    my $read_size = shift;
    my $genome_info = shift;
    system("mkdir $dir/genom_cov");
    my $seq_identity = get_mapping_identity($sam_file);
    $genome_info->{$assemblyid}->{'seq_identity'} = $seq_identity;
    my $cmd = "samtools view -S -b -F 256 $sam_file >$dir/genom_cov/$assemblyid.bam";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt"); # only keep the primary alignments
    $cmd = "samtools sort $dir/genom_cov/$assemblyid.bam -o $dir/genom_cov/$assemblyid.sorted.bam";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt"); 
    $cmd = "samtools rmdup $dir/genom_cov/$assemblyid.sorted.bam $dir/genom_cov/$assemblyid.rd.bam";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt");#remove duplicate
    $cmd = "samtools index $dir/genom_cov/$assemblyid.rd.bam";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt");
    my $lineage = $genome_info->{$assemblyid}->{'lineage'};
    my @ranks = split(/\;/,$lineage);
    my $strain = pop@ranks;
    if ($strain =~ /__/) {
	$strain = substr($strain,3, length($strain)-3);
    }
    $strain =~ s/ /\\ /g;
    my $genome_size = $genome_info->{$assemblyid}->{'genome_size'};
    #system("Rscript $script_dir/bin/read_depth_plot.r -l $assemblyid -n $strain -r $read_size -b $dir/genom_cov/$assemblyid.sorted.bam -s $genome_size -o $dir/genom_cov/$assemblyid\_cov.pdf");
    #plot_genome_coverage("$dir/genom_cov/$taxid.sorted.bam", $ref_genome, $taxid);
    $cmd = "bedtools genomecov -bga -split -ibam $dir/genom_cov/$assemblyid.sorted.bam > $dir/genom_cov/$assemblyid.genom.cov.txt"; #-g $ref_genome > $dir/genom_cov/$assemblyid.genom.cov.txt";
    Inputs::print_and_execute($cmd, "$dir/polish.log.txt");
    my $tot_reads = `samtools view -c -F 260 $dir/genom_cov/$assemblyid.sorted.bam`;
    chomp $tot_reads;
    return($tot_reads, "$dir/genom_cov/$assemblyid.genom.cov.txt");
}

sub removed_by_genome_cov{
    my $dir=shift;
    my $cov_table = shift;
    my $ref_genome = shift;
    my $tot_mapped_reads = shift;
    my $genome_info = shift;
    my $assemblyid = shift;
    my $read_size = shift;
    my $unit_size = 500;
    open(my $f2, ">$dir/genom_cov/$assemblyid.cov.summary.txt") or die;
    
    my $no_fragments = int($genome_info->{$assemblyid}->{'genome_size'}/$unit_size);
    my $maximum_fragments;
    if ($tot_mapped_reads>=$no_fragments) {
	$maximum_fragments = $no_fragments;
    }elsif($tot_mapped_reads>=10){
	$maximum_fragments = $tot_mapped_reads;
    }else{
	return(1); # remove the ones less than 10 reads;
    }
    open(my $f1, "<$cov_table") or die;
    my $upper = $unit_size;
    my $bp = 0;
    my @read_dist;
    my $genome_bp = 0;
    my $covered_genome_bp = 0;
    my $total_bp = 0;
    
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^$assemblyid/) {
	    next;
	}
	my @coln = split(/\t/, $line);
	my $start = $coln[1];
	my $end = $coln[2];
	my $depth = $coln[3];
	$genome_bp = $coln[2];
	if ($depth>0) {
	    $covered_genome_bp += $end-$start;
	    $total_bp += ($end-$start)*$depth;
	}
	
	if ($end>=$upper &&$start<$upper) {
	    while ($upper<=$end) {
		$bp+=($upper-$start)*$depth;
		my $read_no = int($bp*10)/10;
		push(@read_dist, $read_no);
		$start = $upper;
		$upper+=$unit_size;
		$bp=0;
		if ($upper>$end) {
		    $bp = ($end-$start)*$depth;
		}
	    }
	    
	}elsif($end<$upper){
	    $bp+=($end-$start)*$depth;
	    if ($end >= $genome_info->{$assemblyid}->{'genome_size'}) {
		my $read_no = int($bp*10)/10;
		push(@read_dist, $read_no);
	    }
	    
	}
    }
    # false postive: 90% reads located in <= 10% of the genome
    my $total_reads=0;
    foreach my $i (@read_dist){
	$total_reads += $i;
    }
    my $seq_identity = $genome_info->{$assemblyid}->{'seq_identity'};
    my $genome_perc_coverage =  sprintf("%.2f", $covered_genome_bp/$genome_bp*100);
    my $genome_ave_depth =  sprintf("%.1f", $total_bp/$genome_bp);
    my $read_depth =  sprintf("%.1f", $total_bp/$covered_genome_bp);
    print($f2 "seq_identity(%):$seq_identity\ngenome_coverage(%):$genome_perc_coverage\ngenome_ave_depth:$genome_ave_depth\nread_depth:$read_depth\n");
    print($f2 "total_mapped_reads:$tot_mapped_reads\ntotal_fragments:$no_fragments\nmaximum_fragments:$maximum_fragments\n");
    @read_dist = sort{$b<=>$a}@read_dist;
    my $read_distribution = join("\ ", @read_dist);
    my $accu_reads = 0;
    for (my $i=0; $i<=0.15*$maximum_fragments;$i++){ #10% regions of maximum possible regions
	$accu_reads+=$read_dist[$i];
    }
    my $ratio = $accu_reads/$total_reads;
    print($f2 "accumulation_bp:$accu_reads\ntotal_bp:$total_reads\ndistribution_ratio:$ratio\n");
    print($f2 "bp_distribution:\n$read_distribution\n");
    close $f2;
    
    if ($ratio>=0.9) { # if 10% of maximum possible mappable regions contains 90% of all mapped reads, it is a false positive
	return (1);
    }else{
	return (0);
    }
    
}


sub get_genome_info {
    my $genome_info_table = shift;
    my %genome_info;
    open(my $f1, "<$genome_info_table") or die;
    my @title;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if ($line =~ /^assembly_id/) {
	    @title = @coln;
	}else{
	    if (scalar(@coln)<3) {
		next;
	    }
	    for (my $i=1; $i<scalar(@title); $i++){
		$genome_info{$coln[0]}{$title[$i]}=$coln[$i];
	    }
	}
    }
    return \%genome_info;
}

sub sourmash_hier_clustering{
    my $hit_info = $_[0];# hash reference to a copy of the hash
    my $read_info = $_[1];
    my $cutoff = $_[2];
    my $hit_fate;
    %{$hit_fate}=();
    my $new_hit_info;
    %{$new_hit_info}=();
    #my $new_read_info;
    while (scalar(keys(%{$hit_info}))>0) {
	my $highest = get_highest_hit($hit_info);
	my %reads = %{$hit_info->{$highest}};
	$hit_fate -> {$highest} = 'best_hits';
	$new_hit_info->{$highest} = $hit_info->{$highest};
	delete($hit_info->{$highest});
	foreach my $j (keys%reads){
	    #$new_read_info->{$j}->{$highest}++;
	    my %hits = %{$read_info->{$j}};
	    foreach my $i (keys%hits){
		if ($i ne $highest) {
		    delete($hit_info->{$i}->{$j});
		}
	    }
	}
	foreach my $x (keys%{$hit_info}){
	    my %reads = %{$hit_info->{$x}};
	    if (scalar(keys%reads)==0) {
		delete($hit_info->{$x});
		$hit_fate->{$x} = 'reads_reassigned';
	    }elsif (scalar(keys%reads)<$cutoff) {
		delete($hit_info->{$x});
		$hit_fate->{$x} = 'reads_below_cutoff';
	    }
	    
	}
	
    }
    
    foreach my $x (keys%{$new_hit_info}){
	my %reads = %{$new_hit_info->{$x}};
	if (scalar(keys%reads)<$cutoff) {
	    delete($new_hit_info->{$x});
	    $hit_fate->{$x} = 'reads_below_cutoff';
	}
    }
    
    return ($new_hit_info, $hit_fate);
}

sub get_lineage{
    my $node_id =shift;
    my $nodes=shift;
    my $names = shift;
    my $phylogeny = shift;
    my $base_name = $names->{$node_id}{'name'};
    my %rank2name;
    my @lineage;
    my @rank_lvl= @{$phylogeny -> {'ranks'}};
    my @rank_prefix= @{$phylogeny -> {'prefix'}};
    my $num=0;
    my $rank='';
    
    if (length($node_id)==0) {
	return(\@lineage);
    }
    
    while ($node_id != 1 and $num<50) {
	my $name='';
	if (exists($nodes->{$node_id}->{'rank'})) {
	    $rank=$nodes->{$node_id}->{'rank'};
	    $name=$names->{$node_id}->{'name'};
	}else{
	    return(\@lineage);
	}
	
	$rank2name{$rank}=$name;
	my $parent_id = $nodes -> {$node_id} -> {'parent'};
	$node_id = $parent_id;
	$num++;
    }
    
    if (exists($rank2name{'species'})) {
	extract_species_info($base_name, \%rank2name);
    }
    
    for (my $i=0; $i<scalar(@rank_prefix);$i++){
	my $id = $rank_lvl[$i];
	my $rank_id ="unknown";
	if (exists($rank2name{$id})) {
	    $rank_id=$rank2name{$id};
	}
	my $prefix = $rank_prefix[$i];
	my $new_id = $prefix.$rank_id;
	push(@lineage, $new_id);
    }
    
    
    return \@lineage;
    
}

sub extract_species_info{
    my $base_name = shift;
    my $rank2name = shift;
    if (exists($rank2name->{'subspecies'})) {
	if ($base_name =~ / serotype | serovar /) {
	    my @name = split(/ serotype | serovar | str. | str | strain /, $base_name);
	    $rank2name->{'serotype'} = $name[1];
	}else{
	    $rank2name->{'serotype'} = 'unknown'
	}
	if ($base_name =~ / str. | str | strain /) {
	    my @name = split(/ str. | str | strain /, $base_name);
	    $rank2name->{'strain'} = $name[1];
	}else{
	    $rank2name->{'strain'} = 'unknown';
	}
	
    }elsif ($rank2name->{'species'} ne $base_name) {
	$rank2name->{'subspecies'} ='unknown'; # if subspecies are not recorded before hand it is unknown.
	if ($base_name =~ / serotype | serovar /) {
	    my @name = split(/ serotype | serovar | str. | str | strain /, $base_name);
	    $rank2name->{'serotype'} = $name[1];
	}else{
	    $rank2name->{'serotype'} = 'unknown'
	}
	if ($base_name =~ / str. | str | strain /) {
	    my @name = split(/ str. | str | strain /, $base_name);
	    $rank2name->{'strain'} = $name[1];
	}else{
	    $rank2name->{'strain'} = $base_name;
	}
    }else{
	$rank2name->{'subspecies'} ='unknown';
	$rank2name->{'serotype'} = 'unknown';
	if ($base_name =~ / sp. /) {
	    $rank2name->{'species'} = 'unknown'; 
	    my @name = split(/ sp. | str. | str | strain /, $base_name); # Dehalobacter sp. CF, Dehalobacter sp. str. CF, Limnodrilus sp. strain 12
	    $rank2name->{'strain'} = pop(@name);
	}else{
	    $rank2name->{'strain'} = 'unknown';
	}
	
    }
}

sub ranks_to_use {
    my $script_dir = shift;
    my %phylogeny;
    open(my $f1, "<$script_dir/ranks_to_use.csv") or die;
    while (my $line = <$f1>) {
	chomp $line;
	$line =~ s/\"//g;
	$line =~ s/\r//g;
	my @coln = split(/\,/, $line);
	if (scalar(@coln)<1) {
	    next;
	}else{
	    my $key = shift(@coln);
	    $phylogeny{$key}=\@coln;
	}
    }
    return(\%phylogeny);
}


sub extract_species_name{
    my $lineage = shift;
    my @ranks = split(/;/, $lineage);
    my $genus='';
    my $species='';
    foreach my $i (reverse(@ranks)){
	if ($i =~ /^z__/) {
	    next;
	}elsif ($i =~ /^s__/) {
	    $species = substr($i, 3, length($i)-3);
	    if ($species =~ /[Uu]known/){
		$species = 'sp.';
	    }
	    
	}elsif ($i =~ /^g__/) {
	    $genus = substr($i, 3, length($i)-3);
	    if ($genus =~ /[Uu]known/) {
		$genus='';
	    }elsif($genus =~ /\-/){
		$genus='';
	    }else{
		return("$genus $species");
	    }
	}else{
	    $genus = substr($i, 3, length($i)-3);
	    my $rank = substr($i,0,1);
	    if ($genus =~ /[Uu]known/) {
		$genus = '';
	    }else{
		return("($rank)$genus $species");
	    }
	    
	}
    }
}

sub get_mapping_identity{
    my $sam_file = shift;
    open(my $f1, "<$sam_file") or die;
    my $matches = 0;
    my $mismatches = 0;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^@/) {
	    next;
	}
	my @coln = split(/\t/, $line);
	my $MD_field = $coln[12];
	my @info = split(/\:/, $MD_field);
	my @numbers = split(/[^\d]+/, $info[2]);
	
	foreach my $i (@numbers){
	    $matches+=$i;
	}
	my @letters = split(/[^A-Z]+/, $info[2]);
	
	foreach my $i (@letters){
	    $mismatches += length($i);
	}
    }
    my $tot = $matches+$mismatches;
    my $seq_iden = sprintf ("%.4f", $matches/$tot*100);
    return($seq_iden);
}


1;
