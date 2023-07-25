#!/usr/bin/perl
use strict;
use warnings;
#use Inputs;
use Storable 'dclone';
#package RefineCentrifuge;
#----------------------------
my $input_dir = 'Trimmomatic';
my $output_dir = 'Centrifuge';
#my $bin_path = $_[3];
my $threads = 4;
my $database_folder = '/home/stdell/Desktop/shotgun_pipeline/ubuntu/database/';
my $script_dir = '/home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/';
my $phylogeny = ranks_to_use($script_dir);
my $database = "$database_folder/zymo_centrifuge";
my $ncbi_tax_database = "$database_folder/NCBI_taxonomy";
my ($nodes, $names)=read_ncbi_taxa_database($ncbi_tax_database);
my $R1 = "$input_dir/in1717_1_R1.fastq";
my $R2 = "$input_dir/in1717_1_R2.fastq";
my $seq_type = "illumina.pe";
my $sample_id = 'in1717_1';
my $scirpt_dir = '/home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts';
my $cap_size = 2000000;
my ($taxa2abun, $hits_info, $hit_fate) = refine_centrifuge_results("$output_dir", "in1717_1_class.tsv", $nodes, $names, $phylogeny, $seq_type, $R1, $R2, $database, $threads, $sample_id, $script_dir, $cap_size);


#--------------------
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
    my $cap_size = shift;
    my $hit_length_cutoff= 90; # hits less than 80 bp match 
    my $acc2genome_file = read_acc2genome_table($database);
    my $centrifuge_read_cutoff=5; # discard hits with less than 10 reads
    my $centrifuge_subset_size=200000;
    my ($centrifuge_hit_info, $centrifuge_read_info, $centrifuge_hit_info_subset, $centrifuge_read_info_subset, $lineage, $genome_info, $host_reads) = get_centrifuge_hits("$output_dir/$hit_table", $seq_type, $nodes, $names, $phylogeny, $acc2genome_file, $hit_length_cutoff, $centrifuge_subset_size); # use genome_assembly file as the key
    my $centrifuge_hit_size = get_hit_size($centrifuge_hit_info);
    my ($centrifuge_hier_hit_info_subset, $centrifuge_hit_fate_subset) = centrifuge_hier_clustering($centrifuge_hit_info_subset, $centrifuge_read_info_subset, $centrifuge_read_cutoff);
    my $centrifuge_hier_hit_size_subset = get_hit_size($centrifuge_hier_hit_info_subset);
    my ($remain_read_info, $remain_hit_info, $centrifuge_hier_hit_info_subset_new) = substract_hits_and_reads($centrifuge_read_info, $centrifuge_hier_hit_size_subset);
    my ($centrifuge_hier_hit_info_remain, $centrifuge_hit_fate_remain) = centrifuge_hier_clustering($remain_hit_info, $remain_read_info, $centrifuge_read_cutoff);
    my %centrifuge_hits = (%{$centrifuge_hier_hit_info_subset_new}, %{$centrifuge_hier_hit_info_remain});
    my $centrifuge_hier_hit_info = \%centrifuge_hits;
    my %cent_hit_fate = (%{$centrifuge_hit_fate_subset}, %{$centrifuge_hit_fate_remain});
    my $centrifuge_hit_fate = \%cent_hit_fate;
    my $centrifuge_hier_hit_size = get_hit_size($centrifuge_hier_hit_info);
    
    %{$centrifuge_read_info}=();
    %{$centrifuge_hit_info}=();
    
    
    system("mkdir $output_dir/$sample_id");
    $output_dir = "$output_dir/$sample_id";
    my $new_R1 = 'R1.fastq';
    my $new_R2 = 'R2.fastq';
    
    my $ave_read_length = copy_and_filter_host($R1, $R2, $new_R1, $new_R2, $host_reads, $output_dir, $cap_size);
    my $bk_R1 = 'R1_bk.fastq';
    my $bk_R2 = 'R2_bk.fastq';
    system("cp $output_dir/$new_R1 $output_dir/$bk_R1");
    system("cp $output_dir/$new_R2 $output_dir/$bk_R2");
    my $subset_size = 200000;
    my $sampling_size = 100;
    my %bwa_hier_hit_size=();
    my %unmapped_reads=();
    my @centrifuge_hit_no=();
    my $read_cutoff=$cap_size*0.00001;
    if ($read_cutoff<10) {
	$read_cutoff=10;
    }
    
    while (scalar(keys%{$centrifuge_hier_hit_info})>0) {
	push(@centrifuge_hit_no, scalar(keys%{$centrifuge_hier_hit_info}));
	my $cent_hier_hit_size = get_hit_size($centrifuge_hier_hit_info);
	my $top_cent_hits = get_top_hits($cent_hier_hit_size, $sampling_size);
	download_genomes($output_dir, $top_cent_hits, $genome_info);# download and concatenate contigs
	format_ref_db($top_cent_hits, $output_dir);
	bwa_mapping("$output_dir/$new_R1", "$output_dir/$new_R2", $output_dir, $threads, $hit_length_cutoff);
	my ($hit_info, $read_info, $reads) = parse_sam_subset($output_dir, $subset_size, \%unmapped_reads); # do not enrol read information into memory -> out of memory
	my $survived_hit_size = bwa_hier_clustering($hit_info, $read_info, $read_cutoff); # parse the sam file again, reduce the abundance of a hit by 1 if a assigned read was reassigned
	my $left_hit_size=();
	my ($left_hit_info, $left_read_info);
	if ($reads == $subset_size) {
	    ($left_hit_info, $left_read_info, $survived_hit_size) = substract_sam($output_dir, $survived_hit_size, \%unmapped_reads); #update the survived hit size as well
	    $left_hit_size = bwa_hier_clustering($left_hit_info, $left_read_info, $read_cutoff);
	}
	$a = scalar(keys%unmapped_reads);
	save_hits(\%bwa_hier_hit_size, $survived_hit_size, $left_hit_size);
	update_reads($output_dir, $R1, $R2, $new_R1, $new_R2, \%unmapped_reads);
	update_hits($centrifuge_hier_hit_info, \%unmapped_reads, $centrifuge_read_cutoff, $top_cent_hits);
	%unmapped_reads=();
	
    }
    
    #redo mapping with curated hits to get accurate abundance and perform filtration by coverage.
    my $unmapped_read_no = scalar(keys%unmapped_reads);
    my @survival =  keys%bwa_hier_hit_size;
    format_ref_db(\@survival, $output_dir);
    bwa_mapping("$output_dir/$bk_R1", "$output_dir/$bk_R2", $output_dir, $threads, $hit_length_cutoff);
    my ($hit_size, $unique_hit_size) = parse_sam_file($output_dir, \%bwa_hier_hit_size);
    my ($NR_hit_size, $NR_hit_info) = divide_sam_file($output_dir, $unique_hit_size, $hit_size);
    system("echo step6 start: genome coverage filtering\n");
    my $opt_hit_size = filter_hits_by_genome_coverage($NR_hit_size, $output_dir, $genome_info, $ave_read_length, $lineage, $script_dir);
    
    #summarize read fate into a table
    system("echo step7 start: print results\n");
    my %taxa2abun;
    print_summary_results($opt_hit_size, $phylogeny, $lineage, $NR_hit_info, \%taxa2abun, $genome_info, $output_dir);
    print_hit_fate($centrifuge_hit_fate, $centrifuge_hier_hit_size, \%bwa_hier_hit_size, $opt_hit_size, $output_dir, $lineage, $genome_info, $phylogeny);
    return (\%taxa2abun);
    #delete the trimmomatic $R1 and $R2 to save disk space.
    system("rm $R1 $R2");
}

sub update_hits{
    my $hit2reads=$_[0];
    my $unmapped_reads = $_[1];
    my $read_cutoff=$_[2];
    
    my $hits_to_skip=$_[3];
    foreach my $i (@{$hits_to_skip}){
	delete($hit2reads->{$i});
    }
    foreach my $i (keys%{$hit2reads}){	
	my $reads = $hit2reads->{$i};
	foreach my $j (keys%{$reads}){
	    if (exists($unmapped_reads->{$j})) {
		
	    }else{
		delete($hit2reads->{$i}->{$j});
	    }
	    
	}
	if (scalar(keys%{$reads})<$read_cutoff) {
	    delete($hit2reads->{$i});
	}
	
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
    my $R2=$_[1];
    my $new_R1=$_[2];
    my $new_R2=$_[3];
    my $host_reads = $_[4];
    my $dir = $_[5];
    my $cap_size = $_[6];
    my @read_order = sort{$a <=> $b} keys(%{$host_reads});
    open(my $host, ">$dir/host_reads.txt") or die;
    my $host_read_count = scalar(@read_order);
    my $host_read_str = join(",", @read_order);
    print($host ">host_reads $host_read_count\n$host_read_str\n");
    open(my $f1, "<$R1") or die;
    open(my $f2, "<$R2") or die;
    open(my $f3, ">$dir/$new_R1") or die;
    open(my $f4, ">$dir/$new_R2") or die;
    my $counter=0;
    my $added=0;
    my $bp =0;
    while (my $r1 = <$f1>) {
	my $r2 = <$f1>; my $r3 = <$f1>; my $r4 = <$f1>;
	my $x1 = <$f2>; my $x2 = <$f2>; my $x3 = <$f2>; my $x4 = <$f2>;
	if (scalar(@read_order)>0) {
	    if ($counter==$read_order[0]) { #skip host reads
		shift(@read_order); 
	    }else{
		print($f3 "$r1$r2$r3$r4");
		print($f4 "$x1$x2$x3$x4");
		$added++;
	    }
	}else{
	    print($f3 "$r1$r2$r3$r4");
	    print($f4 "$x1$x2$x3$x4");
	    $added++;
	}
	    

	$counter++;
	$bp+=length($r2)-1;
	$bp+=length($x2)-1;
	
	if ($added == $cap_size) { # keep a cap size of the reads 
	    last;
	}
	
	
    }
    close $f1;
    close $f2;
    close $f3;
    close $f4;
    
    my $ave_read_length = int($bp/$counter/2);
    return($ave_read_length);
    
}

sub print_hit_fate{
    my $centrifuge_hit_fate=shift;
    my $centrifuge_hier_hit_size = shift;
    my $bwa_hier_hit_size = shift;
    my $opt_hit_size = shift;
    my $dir = shift;
    my $lineage = shift;
    my $genome_info=shift;
    my $phylogeny=shift;
    open(my $f1, ">$dir/hit_fate.tsv") or die;
    print($f1 "genome_id\ttaxid\ttaxa\tcentrifuge\tbwa\tcoverage\n");
    my @hits = sort {$centrifuge_hit_fate->{$a} cmp $centrifuge_hit_fate->{$b}} keys%{$centrifuge_hit_fate};
    foreach my $i (@hits){
	my $taxid = $genome_info->{$i}->{'taxid'};
	my $ranks = taxid2tax($taxid, $phylogeny, $lineage);
	my $taxa = join(";", @{$ranks});
	my $cent_fate = $centrifuge_hit_fate ->{$i};
	my $bwa_fate ='';
	my $opt_fate='';
	if (exists($centrifuge_hier_hit_size->{$i})) {
	    $cent_fate=$centrifuge_hier_hit_size->{$i};
	    if (exists($bwa_hier_hit_size->{$i})) {
	        $bwa_fate = "$bwa_hier_hit_size->{$i}";
		if (exists($opt_hit_size->{$i})) {
		    $opt_fate = "$opt_hit_size->{$i}";
		}else{
		    $opt_fate = 'removed';
		}
		
	    }else{
		$bwa_fate = 'removed'    
	    }
	}
	print ($f1 "$i\t$taxid\t$taxa\t$cent_fate\t$bwa_fate\t$opt_fate\n");
	
    }
    close $f1;
}


sub print_summary_results{
    my $opt_hit_size=$_[0];
    my $phylogeny = $_[1];
    my $lineage = $_[2];
    my $NR_hit_info = $_[3];
    my $taxa2abun = $_[4];
    my $genome_info = $_[5];
    my $dir = $_[6];
    open(my $f1, ">$dir/hit2read.txt") or die;
    open(my $f2, ">$dir/read2hit.tsv") or die;
    print($f2 "read_id\tgenome_id\ttaxid\tspecies\n");
    foreach my $i (keys%{$opt_hit_size}){
	my $taxid = $genome_info->{$i}->{'taxid'};
	my $ranks = taxid2tax($taxid, $phylogeny, $lineage);
	my $species = $ranks->[3];
	$species = substr($species, 3, length($species)-3);
	my $taxa = join(";", @{$ranks});
	my @reads = @{$NR_hit_info->{$i}};
	my $read_str = join(",", @reads);
	print($f1 ">$i\_$taxid\_$species\n$read_str\n");
	foreach my $j (@reads){
	    print($f2 "$j\t$i\t$taxid\t$species\n");
	}
	$taxa2abun->{$taxa} += $opt_hit_size->{$i};
    }
    close $f1;
    close $f2;
}

sub download_genomes{
    my $dir = $_[0];
    my $hit_list = $_[1];#assembly_id as the key
    my $genome_info = $_[2];
    system("mkdir $dir/ref_genomes");
    foreach my $assembly_id (@{$hit_list}){
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
	
	#system("aws s3 cp s3://zymo-files/WGS_Pipeline/shotgun_database/20200327_genomes/$genome_file $dir/ref_genomes/$assembly_id.fasta.gz");
	system("cp /home/stdell/Desktop/shotgun_pipeline/ubuntu/database/zymo_centrifuge/fasta/$assembly_id.fasta.gz $dir/ref_genomes/$assembly_id.fasta.gz");
	system("gunzip -c $dir/ref_genomes/$assembly_id.fasta.gz > $dir/ref_genomes/$assembly_id.fasta");
	system("rm $dir/ref_genomes/$assembly_id.fasta.gz");
	open(my $f2, "<$dir/ref_genomes/$assembly_id.fasta") or die;
	open(my $f3, ">$dir/ref_genomes/$assembly_id.cat.fasta") or die;
	print($f3 ">$assembly_id\n");
	my $polyN = '';
	my $size = 0;
	while (my $line = <$f2>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		if (length($polyN)==0){
		    $polyN = 'N' x 100;
		}else{
		    print($f3 "$polyN\n");
		    $size += length($polyN);
		}
			    
	    }else{
		print($f3 "$line\n");
		$size += length($line);
			    
	    }
			
	}
	close $f2;
	close $f3;
	if ($size < 100000000) {# remove big eukaryote genomes
	    $genome_info->{$assembly_id}->{'size'}=$size;
	    system("mv $dir/ref_genomes/$assembly_id.cat.fasta $dir/ref_genomes/$assembly_id.fasta");
	}else{
	    system("rm $dir/ref_genomes/$assembly_id.cat.fasta $dir/ref_genomes/$assembly_id.fasta");
	}
	
    }
    
}

sub format_ref_db {
    my $hit_list=$_[0];
    my $dir = $_[1];
    my @genomes='';
    foreach my $i (@{$hit_list}){
	my $file = "$dir/ref_genomes/$i.fasta";
	if (-e $file) {
	    push (@genomes, $file);
	}
    }
    my $str = join(" ", @genomes);
    system ("cat $str > $dir/ref_genomes/ref.fasta");
}

sub bwa_mapping{
    my $R1 = shift;
    my $R2 = shift;
    my $dir = shift;
    my $threads = shift;
    my $hit_length_cutoff = shift;
    system("bwa index -a bwtsw $dir/ref_genomes/ref.fasta");
    system("bwa mem -a -P -t $threads $dir/ref_genomes/ref.fasta $R1 $R2 > $dir/results1.sam");
    #remove the hits of short length from the sam file.
    open(my $f1, "<$dir/results1.sam") or die;
    open(my $f2, ">$dir/results.sam") or die;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^857/) {
	    my $x=0;
	}
	
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
			my $assigned_hit = assign_by_abun($unique_hit_size, \%hits);
			my $size = scalar(@{$hits{$assigned_hit}});
			my @reads;
			if ($size>2) {
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
	    foreach my $i (keys%sam_fhs){
		my $fh = $sam_fhs{$i};
		print($fh "$line\n");
	    }
	}
	
    }
    close $f1;
    
    foreach my $i (keys%sam_fhs){
	close $sam_fhs{$i};
    }
    
    return(\%NR_hit_size, \%NR_hit_info);
    
}


sub bwa_hier_clustering{
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
    my $bk_R2 = $_[2];
    my $R1 =$_[3];
    my $R2 =$_[4];
    my $unmapped_reads = $_[5];
    
    my @read_order = sort{$a <=> $b} keys(%{$unmapped_reads});
    open(my $f1, "<$bk_R1") or die;
    open(my $f2, "<$bk_R2") or die;
    open(my $f3, ">$dir/$R1") or die;
    open(my $f4, ">$dir/$R2") or die;
    my $counter=0;
    while (my $r1 = <$f1>) {
	my $r2 = <$f1>; my $r3 = <$f1>; my $r4 = <$f1>;
	my $x1 = <$f2>; my $x2 = <$f2>; my $x3 = <$f2>; my $x4 = <$f2>;
	if ($counter == $read_order[0]) {
	    print($f3 "$r1$r2$r3$r4");
	    print($f4 "$x1$x2$x3$x4");
	    shift(@read_order);
	    if (scalar(@read_order)==0) {
	        last;
	    }
	}
	$counter++;
	
    }
    close $f1;
    close $f2;
    close $f3;
    close $f4;
    
    
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
	if ($coln[2] !~ /\.|^zymo/) {
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
    my $hit_length_cutoff=shift;
    my $subset_size = shift;
    open(my $f1, "<$hit_table") or die;
    my $title = <$f1>;
    my %read_info;# read info
    my %lineage;
    my %hit_info; # list of reads uniquely assigned to each hit
    my @reads_above_strain;
    my %genome_info;
    my %records;
    my %host_reads;
    my $counter=0;
    my %hit_info_subset;
    my %read_info_subset;
    while (my $line = <$f1>) {
        chomp $line;
        my @coln = split(/\t/, $line);
	if (scalar@coln==0) {next;}
	my $taxid=$coln[2];
	my $acc_no = $coln[1];
	if ($taxid eq '703612') {
	    next;
	}
	my @taxid_to_exclude = (9606,10090,10116,9615,9685);
	my $exclude = 0;
	foreach my $m (@taxid_to_exclude){
	    if ($taxid == $m) {
		$exclude=1;
	    }
	    
	}
	if ($exclude == 1) {
	    $host_reads{$coln[0]}++;
	    next;
	}
	
	
	if ($taxid==0) {next;}
	# filter reads that have short_hit_length
	#if ($seq_type =~ /illumina\.pe/) {# illumina paired
	 #   $hit_length_cutoff = 80;
	#}elsif($seq_type =~ /illumina\.se/){#nanopore
	 #   $hit_length_cutoff = 80;
	#}elsif($seq_type =~ /nanopore/){#nanopore
	 #   $hit_length_cutoff = 80;
	#}elsif($seq_type =~ /pacbio/){#nanopore
	 #   $hit_length_cutoff = 80;
	#}
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
    my %assembly_to_exclude = (
	    'GCA_006364795.1' => 1, #610 P. aeruginosa
	    'GCA_006364815.1' => 1, #611 E. faecalis
	    'GCA_006370485.1' => 1, #612 S. aureus
	    'GCA_006370475.1' => 1, #608 E. coli
	    'GCA_006364035.1' => 1, #607 L. monocytogenes
	    'GCA_006364495.1' => 1, #606 B. subtilis
	    'GCA_012275245.1' => 1, #613 S. cerevisiae
	    'GCA_006370495.1' => 1, #609 S. enterica
	    'GCA_900148615.1' => 2, #contaminated genome
	    'GCA_000010425.1' => 2, #Bifidobacterium adolescentis ATCC 15703, the same strain as LMG 10502
	    'GCA_000016525.1' => 2, #Methanobrevibacter smithii ATCC 35061 is the same strain as DSM861
	    'GCA_000258145.1' => 2, #Ecoli strain W's genome closely related genome to E.coli B766
	    'GCA_000184185.1' => 2, #another genome of strain W
	    'GCA_000019385.1' => 2, #Ecoli strain ATCC 8739 genome closely related genome to E.coli B3008
	    'GCA_003591595.1' => 2, #another genome of strain 8739
    );
    foreach my $read_id (keys%records){
	$counter++;
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
		if (exists($assembly_to_exclude{$assembly_id})) {
		    next; #exclude FDAARGOS genomes that are duplicates of the genomes of the standard D6300
		}
		$genome_info{$assembly_id}{'file'} = $genome_file;
		$genome_info{$assembly_id}{'taxid'} = $taxid;
		$hit_info{$assembly_id}{$read_id}='strain';
		$read_info{$read_id}{$assembly_id}='Centrifuge.unique';
		if ($counter<=$subset_size) {
		    $hit_info_subset{$assembly_id}{$read_id}='strain';
		    $read_info_subset{$read_id}{$assembly_id}='Centrifuge.unique';
		}else{
		    my $a=0;
		}
		
	    }
	    
	}
    }
    
    return (\%hit_info, \%read_info, \%hit_info_subset, \%read_info_subset, \%lineage, \%genome_info, \%host_reads);
}

# my ($remain_read_info, $remain_hit_info, $centrifuge_hier_hit_info_subset) = substract_hits_and_reads($centrifuge_read_info, $centrifuge_hier_hit_info_subset);
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
    my $lineage = $_[4];
    my $script_dir = $_[5];
    my $opt_hit_size;
    foreach my $assemblyid (keys%{$NR_hit_size}){
	my $sam_file = "$dir/sam_files/$assemblyid.sam";
	my $ref_genome = "$dir/ref_genomes/$assemblyid.fasta";
	my ($tot_mapped_reads, $cov_table) = sam2cov($dir, $assemblyid, $sam_file, $ref_genome, $lineage, $script_dir, $read_size, $genome_info);
	if (removed_by_genome_cov($cov_table, $ref_genome, $tot_mapped_reads, $genome_info, $assemblyid, $read_size)) {
	}else{
	    $opt_hit_size->{$assemblyid} = $NR_hit_size->{$assemblyid};
	}
	
    }
    return($opt_hit_size);
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
    system("samtools sort $dir/genom_cov/$assemblyid.bam -o $dir/genom_cov/$assemblyid.sorted.bam");
    system("samtools rmdup $dir/genom_cov/$assemblyid.sorted.bam $dir/genom_cov/$assemblyid.rd.bam");# remove duplicatres
    system("samtools index $dir/genom_cov/$assemblyid.rd.bam");
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


#1;