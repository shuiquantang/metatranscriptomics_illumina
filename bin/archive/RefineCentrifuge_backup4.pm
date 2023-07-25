#!/usr/bin/perl
use strict;
use warnings;
#use Inputs;
use Storable 'dclone';
use Parallel::Loops;
package RefineCentrifuge;

#--------------------

sub TaxaPolisher{
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $bin_path = $_[3];
    my $threads = $_[4];
    my $database_folder = $_[5];
    my $phylogeny = $_[6];
    my $script_dir = $_[7];
    my $database = "$database_folder/zymo_centrifuge";
    my $ncbi_tax_database = "$database_folder/NCBI_taxonomy";
    my $min_hit_length = 80;
    my @samples = keys%{$sample_info};
    my %sample2taxa_abun_table;
    my ($nodes, $names)=read_ncbi_taxa_database($ncbi_tax_database);
    
    #summarize centrifuge results
    #foreach my $i (@samples){
	
    my $pl = Parallel::Loops->new(1);
    $pl -> share(\%sample2taxa_abun_table);
    my $threads_per_samples=1;
    if ($threads/scalar(@samples)>2) {
	$threads_per_samples = int($threads/scalar(@samples));
    }
    
    $pl -> foreach (\@samples, sub{
        my $i = $_;
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1;
	my $R2;
	my $sample_type;
	if ($seq_type =~ /\.pe/) {
	    $sample_type = 2;
	    $R1 = "$input_dir/$i\_R1.paired.fastq.gz";
	    $R2 = "$input_dir/$i\_R2.paired.fastq.gz";
	}else{
	    $sample_type = 1;
	    $R1 = "$input_dir/$i.fastq.gz";
	    $R2 = "";
	}
	
	my $taxa2abun = refine_centrifuge_results("$output_dir", "$i\_class.tsv", $nodes, $names, $phylogeny, $seq_type, $R1, $R2, $database, 15, $i, $script_dir);
	$sample2taxa_abun_table{$i}=$taxa2abun;
	
    });
    #}
    
    my $abun_table_folder='abun_table';
    break_and_print_taxa_abun_table($abun_table_folder, \@samples, \%sample2taxa_abun_table);
}

sub break_and_print_taxa_abun_table{
    my $output_dir = $_[0];
    my @samples = @{$_[1]};
    my $sample2taxa_abun_table = $_[2];
    mkdir ($output_dir);
    my %taxa2sample_abun_table_all;
    my %taxa2sample_abun_table_prokaryote;
    my %taxa2sample_abun_table_eukaryote;
    my %taxa2sample_abun_table_virus;
    
    foreach my $i (@samples){
	foreach my $j (keys(%{$sample2taxa_abun_table->{$i}})){   
	    my $abun = $sample2taxa_abun_table->{$i}->{$j};
	    $taxa2sample_abun_table_all{$j}{$i}=$abun;
	    if ($j =~/k__Bacteria;|k__Archaea;/) {
		$taxa2sample_abun_table_prokaryote{$j}{$i}=$abun;
	    }elsif($j =~ /k__Eukaryota;/){
		$taxa2sample_abun_table_eukaryote{$j}{$i}=$abun;
	    }elsif($j =~ /k__Viruses;/){
		$taxa2sample_abun_table_virus{$j}{$i}=$abun;
	    }
	    
	}
	
    }
    
    print_abun_table("$output_dir/all_abun_table.tsv", \%taxa2sample_abun_table_all, \@samples);
    print_abun_table("$output_dir/prokaryote_abun_table.tsv", \%taxa2sample_abun_table_prokaryote, \@samples);
    print_abun_table("$output_dir/eukaryote_abun_table.tsv", \%taxa2sample_abun_table_eukaryote, \@samples);
    print_abun_table("$output_dir/virus_abun_table.tsv", \%taxa2sample_abun_table_virus, \@samples);    
}

sub print_abun_table{
    my $file = shift;
    my $taxa2sample = shift;
    my $samples = shift;
    open(my $f1, ">$file") or die;
    my $title = join("\t", @{$samples});
    print ($f1 "#taxon_id\t$title\n");
    my @taxa = sort{$a cmp $b} keys(%{$taxa2sample});
    foreach my $i (@taxa){
	my @abun = ($i);
	foreach my $j (@{$samples}){
	    if (exists($taxa2sample->{$i}->{$j})) {
		push(@abun, $taxa2sample->{$i}->{$j});
	    }else{
		push(@abun, 0);
	    }
	    
	}
	my $line = join("\t", @abun);
	print($f1 "$line\n");
    }
    
    close $f1;
}

sub refine_centrifuge_results{
    my $output_dir=shift;
    my $hit_table=shift;
    my $nodes=shift;
    my $names=shift;
    my $phylogeny = shift;
    my $seq_type = shift;
    my $R1 = shift;
    my $R2 = shift;
    my $database = shift;
    my $threads = shift;
    my $sample_id = shift;
    my $script_dir = shift;
    my $reads_cutoff=10;
    my $acc2genome_file = read_acc2genome_table($database);
    system("echo step1 start: parse centrifuge results\n");
    my ($centrifuge_hit_info, $centrifuge_read_info, $lineage, $genome_info) = get_centrifuge_hits("$output_dir/$hit_table", $seq_type, $nodes, $names, $phylogeny, $acc2genome_file); # use genome_assembly file as the key
    my $centrifuge_hit_size = get_hit_size($centrifuge_hit_info);
    system("echo step2 start: hier centrifuge results\n");
    my ($centrifuge_hier_hit_info, $centrifuge_hit_fate) = centrifuge_hier_clustering($centrifuge_hit_info, $centrifuge_read_info, $reads_cutoff);
    %{$centrifuge_hit_info}=();
    system("mkdir $output_dir/$sample_id");
    $output_dir = "$output_dir/$sample_id";
    print_centrifuge_clustering_results($output_dir, $centrifuge_hit_fate);
    my %opt_hit_info=();
    # make a copy of R1 and R2 for subsequent analysis
    my $new_R1 = 'R1.fastq';
    my $new_R2 = 'R2.fastq';
    system("gunzip -c $R1 > $output_dir/$new_R1");
    system("gunzip -c $R2 > $output_dir/$new_R2");
    simplify_fastq_header($output_dir, $new_R1, $new_R2);
    
    my $centrifuge_hier_hit_size = get_hit_size($centrifuge_hier_hit_info);
    %{$centrifuge_hier_hit_info}=();#empty to release ram
    
    prepare_ref_genome($output_dir, $centrifuge_hier_hit_size, $genome_info);
    system("echo step3 start: bwa mapping\n");
    read_mapping_with_bwa($output_dir, $new_R1, $new_R2, $threads);
    system("echo step4 start: parse bwa results\n");
    my ($bwa_hit_info, $bwa_read_info, $reads_without_hits, $average_read_size)=parse_sam_file($output_dir); # use assemblyid as the key
    my $bwa_hit_size = get_hit_size($bwa_hit_info);
    
    system("echo step5 start: hier bwa results\n");
    my ($bwa_hier_hit_info, $bwa_hier_read_info) = bwa_hier_clustering($bwa_read_info, $bwa_hit_info, $reads_cutoff);
    my $bwa_hier_hit_size = get_hit_size($bwa_hier_hit_info);
    %{$bwa_hit_info} = ();#empty to release ram
    
    my $NR_hit_size = create_indiv_sam_files($output_dir, $bwa_hier_hit_size, $bwa_read_info);
    
    
    system("echo step6 start: genome coverage filtering\n");
    my ($opt_hit_info, $opt_read_info) = filter_hits_by_genome_coverage($bwa_hier_hit_info, $output_dir, $genome_info, $average_read_size, $bwa_hier_read_info, $lineage, $script_dir);
    %{$bwa_hier_hit_info} = ();
    my $opt_hit_size = get_hit_size($opt_hit_info);
    #summarize read fate into a table
    
    system("echo step7 start: print results\n");
    print_read_fate($bwa_read_info, $bwa_hier_read_info, $opt_read_info, $reads_without_hits, $output_dir, $lineage, $genome_info);
    print_hit_info($centrifuge_hit_size, $centrifuge_hier_hit_size, $bwa_hit_size, $bwa_hier_hit_size, $opt_hit_size, $output_dir, $lineage, $genome_info, $opt_hit_info);
    %{$opt_hit_info}=();
    my %taxa2abun;
    foreach my $i (keys%{$opt_hit_size}){
	my $taxid = $genome_info->{$i}->{'taxid'};
	my $taxa = taxid2tax($taxid, $phylogeny, $lineage);
	$taxa2abun{$taxa} += $NR_hit_size->{$i};
    }
    return (\%taxa2abun);
}


sub print_centrifuge_clustering_results{
    my $dir = shift;
    my $fate = shift;
    open(my $f1, ">$dir/centrifuge_fate.txt") or die;
    foreach my $i (keys%{$fate}){
	my $note = $fate->{$i};
	print($f1 "$i\t$note\n");
    }
    close $f1;
}

sub create_indiv_sam_files{
    my $dir = $_[0];
    my $hits_kept = $_[1];
    my $bwa_read_info = Storable::dclone($_[2]);
    foreach my $i (keys(%{$bwa_read_info})){
	my %hits = %{$bwa_read_info->{$i}};
	foreach my $j (keys(%hits)){
	    if (exists($hits_kept->{$j})) {
		
	    }else{
		delete($bwa_read_info->{$i}->{$j});
	    }
	    
	}
    }
    my %unique_hit_size;
    my @ambi_reads;
    foreach my $i (keys(%{$bwa_read_info})){
	my @hits = keys%{$bwa_read_info->{$i}};
	if (scalar(@hits)==1) {
	    $unique_hit_size{$hits[0]}++;
	}else{
	    push(@ambi_reads, $i);
	}
    }

    my %ambi_reads_reassigned;
    my %hit_size = %unique_hit_size;
    foreach my $i (@ambi_reads){
	my %hits = %{$bwa_read_info->{$i}};
	my $assigned_hit = assign_by_abun(\%unique_hit_size, \%hits);
	$ambi_reads_reassigned{$i}=$assigned_hit;
	$hit_size{$assigned_hit}++;
    }
    
    system("mkdir $dir/sam_files");
    my @hits = keys%hit_size;
    my $sam_fhs = create_hit_samfile(\@hits, $dir);
    open(my $f2, "<$dir/results.sam") or die;
    while (my $line = <$f2>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
	    my $read_id = $coln[0];
	    my $hit_id = $coln[2];
	    if ($hit_id eq '*') {
		next;
	    }else{
		if (exists($ambi_reads_reassigned{$read_id})) {
		    my $hit_reassigned = $ambi_reads_reassigned{$read_id};
		    if ($hit_reassigned ne $hit_id) {
			next;
		    }
		}
		my $fh = $sam_fhs->{$hit_id};
		print($fh "$line\n");
	    }
	    
	}else{
	    foreach my $i (keys%{$sam_fhs}){
		my $fh = $sam_fhs->{$i};
		print($fh "$line\n");
	    }
	}
	
    }
    foreach my $i (keys%{$sam_fhs}){
	my $fh = $sam_fhs->{$i};
	close $fh;
    }
    close $f2;
    my $bwa_hit_size = resolve_multi_loci_hits(\@hits, $dir);
    return(\%hit_size);
    
}

sub assign_by_abun{
    my $hit_abun = shift;
    my $hits = shift;
    my $sum=0;
    my @hits = keys%{$hits};

    foreach my $i (@hits){
	my $abun = $hit_abun->{$i};
	$sum+=$abun;
	$hits->{$i} = $sum;
    }
    
    my $random_pick = rand($sum);
    foreach my $i (@hits){
	if ($random_pick<=$hits->{$i}) {
	    return ($i);
	}
	
    }
    
}


sub resolve_multi_loci_hits{
    my $hits = shift;
    my $dir = shift;
    my @records;
    my $last_read='';
    my %hit_size;
    
    foreach my $i (@{$hits}){
	open(my $f1, "<$dir/sam_files/$i.sam") or die;
	open(my $f2, ">$dir/sam_files/$i.sam2") or die;
	my $counter=0;
	while (my $line = <$f1>) {
	    chomp $line;
	    if ($line =~ /^@/) {
		print($f2 "$line\n");
	    }else{
		my @coln = split(/\t/, $line);
		my $read_id = $coln[0];
		if (scalar(@records)>0) {
		    if ($read_id eq $last_read) {
		        push(@records, $line);
		    }else{
			my @reads = resolve_repeats(\@records);
			my $lines = join("\n", @reads);
			$hit_size{$i}+=scalar(@reads);
			print($f2 "$lines\n");
			@records =();
			push(@records, $line);
			$last_read = $read_id;
		    }
		}else{
		    push(@records, $line);
		    $last_read = $read_id;
		}
	    }
	    
	}
	close $f1;
	close $f2;
	system("mv $dir/sam_files/$i.sam2 $dir/sam_files/$i.sam");
    }
    
    
    return(\%hit_size);
}

sub resolve_repeats{ 
    # some reads can be mapped to multiple positions in the genome due to repeats
    # this function randomly assign the fwd and reverse read to one position
    my $records = shift;
    my @fwd;
    my @rev;
    foreach my $i (@{$records}){
	my @coln = split(/\t/, $i);
	my $flag = $coln[1];
	if ($flag & 0x40) { # binary operator to judge wether or not the read belong to the fwd or rev read.
	    push(@fwd, $i);
	}
	if ($flag & 0x80) {
	    push(@rev, $i);
	}
    }   
    my @reads;
    if (scalar(@fwd)>0) {
	my $index = rand(@fwd);
	my $fwd_choice = $fwd[$index];
	push(@reads, $fwd_choice);
   
    }
    if (scalar(@rev)>0) {
	my $index = rand(@rev);
	my $rev_choice = $rev[$index];
	push(@reads, $rev_choice);
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

sub simplify_fastq_header{
    my $dir=shift;
    my $R1=shift;
    my $R2=shift;
    open(my $f1, "<$dir/$R1") or die;
    open(my $f2, "<$dir/$R2") or die;
    open(my $f3, ">$dir/new.R1.fastq") or die;
    open(my $f4, ">$dir/new.R2.fastq") or die;
    open(my $list, ">$dir/read_list.txt") or die;		
    my $counter=0;
    while (my $fwd1 = <$f1>) {
	my $fwd2 = <$f1>; my $fwd3 = <$f1>; my $fwd4 = <$f1>;
	my $rev1 = <$f2>; my $rev2 = <$f2>; my $rev3 = <$f2>; my $rev4 = <$f2>;
	if ($fwd1 =~ /^@/) {
	    if ($rev1 =~ /^@/) {
		my @header_fwd = split(/\ /,$fwd1);
		my @header_rev = split(/\ /, $rev1);
		my $read_id = substr($header_fwd[0], 1, length($header_fwd[0])-1);
		$header_fwd[0] = '@'.$counter;
		$header_rev[0] = '@'.$counter;
		$counter++;
		$fwd1 = join(" ", @header_fwd);
		$rev1 = join(" ", @header_rev);
		print($list "$read_id\n");
		print($f3 $fwd1.$fwd2.$fwd3.$fwd4);
		print($f4 $rev1.$rev2.$rev3.$rev4);
		
	    }else{
		print("error in fastq file: $R2\n");
	    }
	    
	}else{
	    print("error in fastq file: $R1\n");
	}
	
	
    }
    close $f1;
    close $f2;
    close $f3;
    close $f4;
    system("mv $dir/new.R1.fastq $dir/$R1");
    system("mv $dir/new.R2.fastq $dir/$R2");
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

sub get_centrifuge_hits{
    my $hit_table = shift;
    my $seq_type = shift;
    my $nodes = shift;
    my $names = shift;
    my $phylogeny = shift;
    my $acc2genome_file = shift;
    open(my $f1, "<$hit_table") or die;
    my $title = <$f1>;
    my %read_info;# read info
    my %lineage;
    my %hit_info; # list of reads uniquely assigned to each hit
    my @reads_above_strain;
    my %genome_info;
    my %records;
    while (my $line = <$f1>) {
        chomp $line;
        my @coln = split(/\t/, $line);
	if (scalar@coln==0) {next;}
	my $taxid=$coln[2];
	my $acc_no = $coln[1];
	if ($taxid eq '703612') {
	    next;
	}
	if ($taxid==0) {next;}
	my $hit_length_cutoff=80;# filter reads that have short_hit_length
	if ($seq_type =~ /illumina\.pe/) {# illumina paired
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /illumina\.se/){#nanopore
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /nanopore/){#nanopore
	    $hit_length_cutoff = 80;
	}elsif($seq_type =~ /pacbio/){#nanopore
	    $hit_length_cutoff = 80;
	}
	if ($coln[5]<$hit_length_cutoff) {# filter hits that have a hit_length
	    next;
	}

	# enrol previous read infor if this line belongs to a new read.
	my $read_id = $coln[0];
	if ($acc_no !~ /\.|^zymo/) {
	    push(@reads_above_strain, $read_id);
	}
	push(@{$records{$read_id}}, $line);
    }
    
    foreach my $j (@reads_above_strain){
	delete($records{$j});
    }
    foreach my $read_id (keys%records){
	my @lines = @{$records{$read_id}};
	foreach my $j (@lines){
	    my @coln = split(/\t/, $j);
	    my $taxid = $coln[2];
	    my $phylo_lineage;
	    if (exists($lineage{$taxid})) {
		$phylo_lineage = $lineage{$taxid};
	    }else{
		$phylo_lineage = get_lineage($taxid, $nodes, $names, $phylogeny);
		if (scalar(@{$phylo_lineage})>0) { # skip the hits to return no lineage due to errors in the database
		    $lineage{$taxid}=$phylo_lineage;
		}else{
		    next;
		}
		if ($phylo_lineage->[5] =~/^g__Shigella$/) { # merge Shigella genus into Escherichia
		    $phylo_lineage->[5] = 'g__Escherichia';
		}
	    }
	    my $acc_no = $coln[1];
	    
	    if ($acc_no =~ /\.|^zymo/) {
		my $assembly_id;
		my $genome_file;
		if ($acc_no =~ /^zymo/) { # zymo genomes use taxid as assembly_id
		    my @name = split(/s[0-9]/, $acc_no);
		    $assembly_id = $name[0];
		    $genome_file = $name[0].'.fasta.gz';
		}else{
		    if (exists($acc2genome_file -> {$acc_no})) {
			$genome_file = $acc2genome_file -> {$acc_no};
			my @names = split(/\_/, $genome_file);
			$assembly_id = $names[0].'_'.$names[1];
		    }else{
			print("accession number error: $acc_no\n");
			next;
		    }
		    
		}
		$genome_info{$assembly_id}{'file'} = $genome_file;
		$genome_info{$assembly_id}{'taxid'} = $taxid;
		$hit_info{$assembly_id}{$read_id}='strain';
		$read_info{$read_id}{$assembly_id}='Centrifuge.unique';
	    }
	    
	}
    }
    
    
    
    	
    #use the array index as the read id
	
    #my %unique_hit_info;
    #foreach my $i (keys%read_info){
	#my %hits = %{$read_info{$i}};
	#my @hits_list = keys%hits;
	#if (scalar(@hits_list)<=3) {
	 #   foreach my $j (@hits_list){
	#	$unique_hit_info{$j}{'reads'}{$i}='strain';
	 #   }
	#}
	
    #}
    return (\%hit_info, \%read_info, \%lineage, \%genome_info);
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
    my $new_lineage = join(";", @new_ranks);
    return $new_lineage;
}

sub print_hit_info{
    my $centrifuge_hit_size = $_[0];
    my $centrifuge_hier_hit_size = $_[1];
    my $bwa_hit_size = $_[2];
    my $bwa_hier_hit_size = $_[3];
    my $opt_hit_size = $_[4];
    my $output_dir = $_[5];
    my $lineage = $_[6];
    my $genome_info=$_[7];
    my $opt_hit_info = $_[8];
    my @final_hits = keys(%{$opt_hit_size});
    my @others;
    foreach my $i (keys%{$centrifuge_hit_size}){
	if (exists($opt_hit_size->{$i})) {
	    
	}else{
	    push(@others, $i);
	}
	
    }
    push(@final_hits, @others);
    
    open(my $f1, ">$output_dir/hit_info.tsv") or die;
    print ($f1 "hit_assembly_id\treads\ttaxid\tlineage\n");
    open(my $f2, ">$output_dir/hits2reads.tsv") or die;
    open(my $f3, ">$output_dir/hit_fate.tsv") or die;
    print($f3 "assembly_id\ttaxid\tlineage\tcentrifuge\tcentrifuge_hier\tbwa\tbwa_hier\tgenome_cov\n");
    foreach my $i (@final_hits){
	my $taxid = $genome_info->{$i}->{'taxid'};
	my $taxa = join (";", @{$lineage->{$taxid}});
	if (exists($opt_hit_size->{$i})) {
	    my @reads = keys%{$opt_hit_info->{$i}};
	    my $size = scalar(@reads);
	    my $read_str = join(",",@reads);
	    print ($f1 "$i\t$size\ttaxid\t$taxa\n");
	    print ($f2 ">$i\n$read_str\n");
	}

	my $centrifuge_size=0;
	if (exists($centrifuge_hit_size->{$i})) {
	    $centrifuge_size = $centrifuge_hit_size->{$i};
	}
	my $centrifuge_hier_size=0;
	if (exists($centrifuge_hier_hit_size->{$i})) {
	    $centrifuge_hier_size = $centrifuge_hier_hit_size->{$i};
	}
	my $bwa_size=0;
	if (exists($bwa_hit_size->{$i})) {
	    $bwa_size = $bwa_hit_size->{$i};
	}
	my $bwa_hier_size=0;
	if (exists($bwa_hier_hit_size->{$i})) {
	    $bwa_hier_size = $bwa_hier_hit_size->{$i};
	}
	my $opt_size=0;
	if (exists($opt_hit_size->{$i})) {
	    $opt_size = $opt_hit_size->{$i};
	}
	
	print($f3 "$i\t$taxid\t$taxa\t$centrifuge_size\t$centrifuge_hier_size\t$bwa_size\t$bwa_hier_size\t$opt_size\n");
	
	
    }
    
    

    close $f3;
    close $f1;
    close $f2;
}

sub filter_hits_by_genome_coverage{
    my $hier_hit_info = $_[0];
    my $dir = $_[1];
    my $genome_info = $_[2];
    my $read_size = $_[3];
    my $hier_read_info =$_[4];
    my $lineage = $_[5];
    my $script_dir = $_[6];
    my $opt_hit_info;
    my $opt_read_info;
    foreach my $assemblyid (keys%{$hier_hit_info}){
	my $sam_file = "$dir/sam_files/$assemblyid.sam";
	my $ref_genome = "$dir/ref_genomes/$assemblyid.fasta";
	my ($tot_mapped_reads, $cov_table) = sam2cov($dir, $assemblyid, $sam_file, $ref_genome, $lineage, $script_dir, $read_size, $genome_info);
	if (removed_by_genome_cov($cov_table, $ref_genome, $tot_mapped_reads, $genome_info, $assemblyid, $read_size)) {
	}else{
	    $opt_hit_info->{$assemblyid} = $hier_hit_info->{$assemblyid};
	    foreach my $i (keys%{$opt_hit_info->{$assemblyid}}){
		$opt_read_info->{$i}=$assemblyid;
	    }
	}
	
    }
    return($opt_hit_info, $opt_read_info);
}

sub sam2cov{
    my $dir=shift;
    my $assemblyid = shift;
    my $sam_file = shift;
    my $ref_genome = shift;
    my $lineage = shift;
    my $script_dir = shift;
    my $read_size = shift;
    my $genome_info = shift;
    system("mkdir $dir/genom_cov");
    system("samtools view -S -b $sam_file >$dir/genom_cov/$assemblyid.bam");
    system("samtools rmdup $dir/genom_cov/$assemblyid.bam $dir/genom_cov/$assemblyid.rd.bam");# remove duplicatres
    system("samtools sort $dir/genom_cov/$assemblyid.rd.bam -o $dir/genom_cov/$assemblyid.sorted.bam");
    system("samtools index $dir/genom_cov/$assemblyid.sorted.bam");
    my $taxid = $genome_info->{$assemblyid}->{'taxid'};
    my $strain = $lineage->{$taxid}->[9];
    $strain = substr($strain,3, length($strain)-3);
    $strain =~ s/ /\\ /g;
    my $genome_size = $genome_info->{$assemblyid}->{'size'};
    system("Rscript $script_dir/bin/read_depth_plot.r -l $assemblyid -n $strain -r $read_size -b $dir/genom_cov/$assemblyid.sorted.bam -s $genome_size -o $dir/genom_cov/$assemblyid\_cov.pdf");
    #plot_genome_coverage("$dir/genom_cov/$taxid.sorted.bam", $ref_genome, $taxid);
    system("bedtools genomecov -bga -split -ibam $dir/genom_cov/$assemblyid.sorted.bam -g $ref_genome > $dir/genom_cov/$assemblyid.genom.cov.txt");
    my $tot_reads = `samtools view -c -F 260 $dir/genom_cov/$assemblyid.sorted.bam`;
    chomp $tot_reads;
    return($tot_reads, "$dir/genom_cov/$assemblyid.genom.cov.txt");
    
}

sub removed_by_genome_cov{
    my $cov_table = shift;
    my $ref_genome = shift;
    my $tot_mapped_reads = shift;
    my $genome_info = shift;
    my $assemblyid = shift;
    my $read_size = shift;
    my $no_fragments = 100;
    if ($tot_mapped_reads>=100) {
	$no_fragments = 100;
    }elsif($tot_mapped_reads>=10){
	$no_fragments = $tot_mapped_reads;
    }else{
	return(1); # remove the ones less than 10 reads;
    }   
    my $unit_size = int($genome_info->{$assemblyid}->{'size'}/$no_fragments);
    open(my $f1, "<$cov_table") or die;
    my $upper = $unit_size;
    my $bp = 0;
    my @read_dist;
    
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^$assemblyid/) {
	    next;
	}
	my @coln = split(/\t/, $line);
	my $start = $coln[1];
	my $end = $coln[2];
	my $depth = $coln[3];
	if ($end>=$upper &&$start<$upper) {
	    while ($upper<=$end) {
		$bp+=($upper-$start)*$depth;
		my $read_no = int($bp/$read_size)+1;
		push(@read_dist, $read_no);
		$start = $upper;
		$upper+=$unit_size;
		if ($upper>$end) {
		    $bp = ($end-$start)*$depth;
		}
	    }
	    
	}elsif($end<$upper){
	    $bp+=($end-$start)*$depth;
	    if ($end >= $genome_info->{$assemblyid}->{'size'}) {
		my $read_no = int($bp/$read_size)+1;
		push(@read_dist, $read_no);
	    }
	    
	}
	
	
    }
    # false postive: 90% reads located in <= 10% of the genome
    my $total_reads=0;
    foreach my $i (@read_dist){
	$total_reads += $i;
    }
    
    @read_dist = sort{$b<=>$a}@read_dist;
    my $accu_reads = 0;
    for (my $i=0; $i<=0.1*$no_fragments;$i++){ #10% regions of the genome
	$accu_reads+=$read_dist[$i];
    }
    if ($accu_reads/$total_reads>=0.9) { # 90% of the reads
	return (1);
    }else{
	return (0);
    }
    
}

sub print_read_fate{
    my $bwa_read_info = $_[0];
    my $bwa_hier_read_info = $_[1];
    my $opt_read_info = $_[2];
    my $reads_without_hits = $_[3];
    my $output_dir = $_[4];
    my $lineage = $_[5];
    my $genome_info=$_[6];
    open(my $f1, ">$output_dir/read_fate.tsv") or die;
    print($f1 "read_id\tbwa\tbwa_hier\tfinal_by_cov\tspecies\n");
    my @reads = keys%{$reads_without_hits};
    my @reads_with_hits = keys%{$bwa_read_info};
    push(@reads, @reads_with_hits);
    foreach my $i (@reads){
	my $bwa='';
	if (exists($bwa_read_info->{$i})) {
	    my @hits = keys%{$bwa_read_info->{$i}};
	    $bwa = join(",", @hits);
	}
	my $bwa_hier ='';
	if (exists($bwa_hier_read_info->{$i})) {
	    $bwa_hier = $bwa_hier_read_info->{$i};
	}
	my $cov_filtration='';
	my $species = '';
	if (exists($opt_read_info->{$i})) {
	    $cov_filtration = $opt_read_info->{$i};
	    my $taxid = $genome_info->{$cov_filtration}->{'taxid'};
	    $species = $lineage -> {$taxid}->[6];
	    $species = substr($species, 3, length($species)-3);
	}
	print($f1 "$i\t$bwa\t$bwa_hier\t$cov_filtration\t$species\n");
    }
    close $f1;
    
}

sub get_centrifuge_results{
    my $centrifuge_read_info = $_[0];
    my $centrifuge_note='';
    my $read_id = $_[1];
    if (exists($centrifuge_read_info->{$read_id})) {
	    my %hits = %{$centrifuge_read_info->{$read_id}};
	    my @assemblies = keys%hits;
	    $centrifuge_note = join(",", @assemblies);
    }
    return($centrifuge_note);
    
}


sub read_acc2genome_table {
    my $database = shift;
    my %acc2genome_file;
    open(my $f1, "<$database/acc2genome.tsv") or die;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	$acc2genome_file{$coln[0]}=$coln[1];
    }
    return \%acc2genome_file;
}



sub prepare_ref_genome{
    my $dir = $_[0];
    my $hit_info = $_[1];#assembly_id as the key
    my $genome_info = $_[2];
    system("mkdir $dir/ref_genomes");
    open(my $f1, ">$dir/ref_genomes/ref.fasta") or die;
    foreach my $assembly_id (keys%{$hit_info}){
	my $taxid = $genome_info->{$assembly_id}->{'taxid'};
	my $genome_file = $genome_info->{$assembly_id}->{'file'};
	my @taxid_to_exclude = (9606,10090,10116,9615,9685);
	my $exclude = 0;
	foreach my $m (@taxid_to_exclude){
	    if ($taxid == $m) {
		$exclude=1;
	    }
	    
	}
	if ($exclude == 1) {
	    next;
	}
	

	system("aws s3 cp s3://zymo-files/WGS_Pipeline/shotgun_database/20200327_genomes/$genome_file $dir/ref_genomes/$assembly_id.fasta.gz");
	#system("cp /home/stdell/Desktop/shotgun_pipeline/ubuntu/database/zymo_centrifuge/fasta/$genome_file $dir/ref_genomes/$assembly_id.fasta.gz");
	system("gunzip -c $dir/ref_genomes/$assembly_id.fasta.gz > $dir/ref_genomes/$assembly_id.fasta");
	system("rm $dir/ref_genomes/$assembly_id.fasta.gz");
	open(my $f2, "<$dir/ref_genomes/$assembly_id.fasta") or die;
	open(my $f3, ">$dir/ref_genomes/$assembly_id.cat.fasta") or die;
	print($f1 ">$assembly_id\n");
	print($f3 ">$assembly_id\n");
	my $polyN = '';
	my $size = 0;
	while (my $line = <$f2>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		if (length($polyN)==0){
		    $polyN = 'N' x 100;
		}else{
		    print($f1 "$polyN\n");
		    print($f3 "$polyN\n");
		    $size += length($polyN);
		}
			    
	    }else{
		print($f1 "$line\n");
		print($f3 "$line\n");
		$size += length($line);
			    
	    }
			
	}
	close $f2;
	close $f3;
	system("mv $dir/ref_genomes/$assembly_id.cat.fasta $dir/ref_genomes/$assembly_id.fasta");
	$genome_info->{$assembly_id}->{'size'}=$size;
	
    }
    
    close $f1;
}

sub centrifuge_hier_clustering{
    my $hit_info = $_[0];# hash reference to a copy of the hash
    my $read_info = $_[1];
    my $cutoff = $_[2];
    my $hit_fate;
    my $new_hit_info;
    #my $new_read_info;
    while (scalar(keys(%{$hit_info}))>0) {
	my $highest = get_highest_hit($hit_info);
	my %reads = %{$hit_info->{$highest}};
	$hit_fate -> {$highest} = 'centrifuge_best_hits';
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
		$hit_fate->{$x} = 'centrifuge_reassigned';
	    }elsif (scalar(keys%reads)<$cutoff) {
		delete($hit_info->{$x});
		$hit_fate->{$x} = 'centrifuge_below_cutoff';
	    }
	    
	}
	
    }
    
    return ($new_hit_info, $hit_fate);
}

sub bwa_hier_clustering{
    my $read_info = $_[0];
    my $hit_info = $_[1];
    my $reads_cutoff= $_[2];
    my %opt_hit_info;
    my %opt_read_info;
    my $rank = 1;
    while (scalar(keys(%{$hit_info}))>0) {
	my $highest = get_highest_hit($hit_info);
	$rank++;
	my %reads = %{$hit_info->{$highest}};
	$opt_hit_info{$highest} = $hit_info->{$highest};
	delete($hit_info->{$highest});
	foreach my $j (keys%reads){
	    $opt_read_info{$j}=$highest;
	    my %hits = %{$read_info->{$j}};
	    foreach my $i (keys%hits){
		if ($i ne $highest) {
		    delete($hit_info->{$i}->{$j});
		}
	    }
	}
	foreach my $x (keys%{$hit_info}){
	    my %reads = %{$hit_info->{$x}};	    
	    if (scalar(keys%reads)<$reads_cutoff) {
		delete($hit_info->{$x});
	    }
	    
	}
	
    }
    
    return (\%opt_hit_info, \%opt_read_info);
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
    open(my $f1, "<$dir/results.sam") or die;
    my %hit_info;
    my %read_info;
    my %reads_without_hits;
    my @headers;
    my $tot_bp = 0;
    my $reads = 0;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
	    if ($coln[2] ne '*') {
		    my $assembly_id = $coln[2];
		    push(@{$hit_info{$assembly_id}{$coln[0]}}, $coln[1]);
		    $read_info{$coln[0]}{$assembly_id}++;
		    if ($coln[9] ne '*') {
			$tot_bp+=length($coln[9]);
			$reads+=1;
		    }
	    }else{
		$reads_without_hits{$coln[0]}++;
	    }
	}else{
	    push(@headers, $line);
	}
	
    }
    close $f1;
    my $average_read_size = $tot_bp/$reads;
    
    return(\%hit_info, \%read_info, \%reads_without_hits, $average_read_size);
    
}



sub read_mapping_with_bwa{
    my $dir = shift;
    my $R1 = shift;
    my $R2 = shift;
    my $threads = shift;
    system("bwa index -a bwtsw $dir/ref_genomes/ref.fasta");
    system("bwa mem -a -P -t $threads $dir/ref_genomes/ref.fasta $dir/$R1 $dir/$R2 > $dir/results.sam");

}

sub get_lineage_nodes{
    my $node_id =shift;
    my $nodes=shift;
    my @lineage_nodes;
    my $num=0;
    while ($node_id != 1 and $num<50) {
	push(@lineage_nodes, $node_id);
	my $parent_id = $nodes -> {$node_id} -> {'parent'};
	$node_id = $parent_id;
	$num++;
    }
    return (\@lineage_nodes);
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


sub read_ncbi_taxa_database{
    my $ncbi_database_folder=shift;
    my %nodes;
    open(my $NCBI_nodes, "<$ncbi_database_folder/zymo.tree") or die;
    while (my $line = <$NCBI_nodes>) {
        chomp $line;
        my @coln = split(/\t/, $line);
        if (scalar@coln>1) {
            $nodes{$coln[0]}{'parent'}=$coln[2];
	    $nodes{$coln[0]}{'rank'}=$coln[4];
        }
    }
    close $NCBI_nodes;

    my %names;
    open(my $NCBI_names, "<$ncbi_database_folder/zymo.name") or die;
    while (my $line = <$NCBI_names>) {
        chomp $line;
        my @coln = split(/\t/, $line);
        my $name_class = $coln[6];
        if (exists($names{$coln[0]}{'name'})) {
	    if ($name_class eq 'scientific name') {
		$names{$coln[0]}{'name_class'}=$coln[6];
		$names{$coln[0]}{'note'}=$coln[4];
		$names{$coln[0]}{'name'}=$coln[2];
	    }
	}else{
	    $names{$coln[0]}{'name_class'}=$coln[6];
	    $names{$coln[0]}{'note'}=$coln[4];
	    $names{$coln[0]}{'name'}=$coln[2];
	}
    }
    close $NCBI_names;
    return(\%nodes, \%names);
}

1;