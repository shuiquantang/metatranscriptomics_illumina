#!/usr/bin/perl
use strict;
use warnings;
#use Inputs;
use List::MoreUtils qw(any);
use Statistics::Basic qw(:all);
use POSIX qw/ceil/;
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
my $min_hit_length = 80;
my ($nodes, $names)=read_ncbi_taxa_database($ncbi_tax_database);
my $R1 = "$input_dir/in709_1_R1.paired.fastq.gz";
my $R2 = "$input_dir/in709_1_R2.paired.fastq.gz";
my $seq_type = "illumina.pe";
my $sample_id = 'in709_1';
my $scirpt_dir = '/home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts';

my ($taxa2abun, $hits_info, $hit_fate) = refine_centrifuge_results("$output_dir", "in709_1_class.tsv", $nodes, $names, $phylogeny, $seq_type, $R1, $R2, $database, $threads, $sample_id, $script_dir);

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
    my $acc2genome_file = read_acc2genome_table($database);
    my ($hit_info, $read_info, $lineage) = get_centrifuge_hits("$output_dir/$hit_table", $seq_type, $nodes, $names, $phylogeny); # use acc_no as the key
    system("mkdir $output_dir/$sample_id");
    $output_dir = "$output_dir/$sample_id";
    my %opt_hit_info=();
    # make a copy of R1 and R2 for subsequent analysis
    my $new_R1 = 'R1.fastq';
    my $new_R2 = 'R2.fastq';
    system("gunzip -c $R1 > $output_dir/$new_R1");
    system("gunzip -c $R2 > $output_dir/$new_R2");
    my $genome_size = prepare_ref_genome($output_dir, $hit_info, $acc2genome_file);
    read_mapping_with_bwa($output_dir, $new_R1, $new_R2, $threads);
    my ($bwa_hit_info, $bwa_read_info, $reads_without_hits, $average_read_size)=parse_sam_file($output_dir); # use taxid as the key
    my ($hier_hit_info, $hier_read_info, $hit_fate) = filter_reassign_reads($bwa_read_info, $bwa_hit_info);
    my ($opt_hit_info, $opt_read_info) = filter_hits_by_genome_coverage($hier_hit_info, $hit_fate, $output_dir, $genome_size, $average_read_size, $hier_read_info, $lineage, $script_dir);
    #summarize read fate into a table
    print_read_fate($opt_read_info, $hier_read_info, $bwa_read_info, $read_info, $reads_without_hits, $hit_info, $output_dir, $lineage);
    print_hit_info($output_dir, $opt_hit_info, $hit_fate, $lineage);
    my %taxa2abun;
    foreach my $i (keys%{$opt_hit_info}){
	my %reads = %{$opt_hit_info->{$i}};
	my $taxa = taxid2tax($i, $phylogeny, $lineage);
	$taxa2abun{$taxa} = scalar(keys%reads);
    }
    return (\%taxa2abun);
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
    my $output_dir = shift;
    my $hit_info = shift;
    my $hit_fate = shift;
    my $lineage = shift;
    my @hits = keys(%{$hit_info});
    @hits = sort{scalar(keys%{$hit_info->{$b}}) <=> scalar(keys%{$hit_info->{$a}})}@hits;
    open(my $f1, ">$output_dir/hit_info.tsv") or die;
    print ($f1 "hit_taxid\treads\tlineage\n");
    open(my $f2, ">$output_dir/hits2reads.tsv") or die;
    foreach my $i (@hits){
	my $taxa = join (";", @{$lineage->{$i}});
	my $size = scalar(keys%{$hit_info->{$i}});
	my @reads = keys%{$hit_info->{$i}};
	my $read_str = join(",",@reads);
	print ($f1 "$i\t$size\t$taxa\n");
	print ($f2 ">$i\n$read_str\n");
    }
    open(my $f3, ">$output_dir/hit_fate.tsv") or die;
    print($f3 "taxid\thit_fate\n");
    my @order = sort{$hit_fate->{$a} cmp $hit_fate->{$b}} keys%{$hit_fate};
    foreach my $i (@order){
	my $fate = $hit_fate->{$i};
	print($f3 "$i\t$fate\n");
    }
    close $f3;
    close $f1;
    close $f2;
}

sub filter_hits_by_genome_coverage{
    my $hier_hit_info = $_[0];
    my $hit_fate = $_[1];
    my $dir = $_[2];
    my $genome_size = $_[3];
    my $read_size = $_[4];
    my $hier_read_info =$_[5];
    my $lineage = $_[6];
    my $script_dir = $_[7];
    my $opt_hit_info;
    my $opt_read_info;
    foreach my $taxid (keys%{$hier_hit_info}){
	my $sam_file = "$dir/sam_files/$taxid.sam";
	my $ref_genome = "$dir/ref_genomes/$taxid.fasta";
	my ($tot_mapped_reads, $cov_table) = sam2cov($dir, $taxid, $sam_file, $ref_genome, $lineage, $script_dir, $read_size, $genome_size);
	if (removed_by_genome_cov($cov_table, $ref_genome, $tot_mapped_reads, $genome_size, $taxid, $read_size)) {
	    $hit_fate->{$taxid} = 'removed_by_genom_cov';
	}else{
	    $opt_hit_info->{$taxid} = $hier_hit_info->{$taxid};
	    foreach my $i (keys%{$opt_hit_info->{$taxid}}){
		$opt_read_info->{$i}=$taxid;
	    }
	}
	
    }
    return($opt_hit_info, $opt_read_info);
}

sub sam2cov{
    my $dir=shift;
    my $taxid = shift;
    my $sam_file = shift;
    my $ref_genome = shift;
    my $lineage = shift;
    my $script_dir = shift;
    my $read_size = shift;
    my $genome_size = shift;
    system("mkdir $dir/genom_cov");
    system("samtools view -S -b $sam_file >$dir/genom_cov/$taxid.bam");
    system("samtools rmdup $dir/genom_cov/$taxid.bam $dir/genom_cov/$taxid.rd.bam");# remove duplicatres
    system("samtools sort $dir/genom_cov/$taxid.rd.bam -o $dir/genom_cov/$taxid.sorted.bam");
    system("samtools index $dir/genom_cov/$taxid.sorted.bam");
    my $strain = $lineage->{$taxid}->[9];
    $strain = substr($strain,3, length($strain)-3);
    $strain =~ s/ /\\ /g;
    my $geno_size = $genome_size->{$taxid};
    system("Rscript $script_dir/bin/read_depth_plot.r -l $taxid -n $strain -r $read_size -b $dir/genom_cov/$taxid.sorted.bam -s $geno_size -o $dir/genom_cov/$taxid\_cov.pdf");
    #plot_genome_coverage("$dir/genom_cov/$taxid.sorted.bam", $ref_genome, $taxid);
    system("bedtools genomecov -bga -split -ibam $dir/genom_cov/$taxid.sorted.bam -g $ref_genome > $dir/genom_cov/$taxid.genom.cov.txt");
    my $tot_reads = `samtools view -c -F 260 $dir/genom_cov/$taxid.sorted.bam`;
    chomp $tot_reads;
    return($tot_reads, "$dir/genom_cov/$taxid.genom.cov.txt");
    
}

sub removed_by_genome_cov{
    my $cov_table = shift;
    my $ref_genome = shift;
    my $tot_mapped_reads = shift;
    my $genome_size = shift;
    my $taxid = shift;
    my $read_size = shift;
    my $no_fragments = 100;
    if ($tot_mapped_reads>=100) {
	$no_fragments = 100;
    }elsif($tot_mapped_reads>=10){
	$no_fragments = $tot_mapped_reads;
    }else{
	return(1); # remove the ones less than 10 reads;
    }   
    my $unit_size = int($genome_size->{$taxid}/$no_fragments);
    open(my $f1, "<$cov_table") or die;
    my $upper = $unit_size;
    my $bp = 0;
    my @read_dist;
    
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^$taxid/) {
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
	    if ($end >= $genome_size->{$taxid}) {
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
    my $opt_read_info= $_[0];
    my $hier_read_info = $_[1];
    my $bwa_fate = $_[2];
    my $centrifuge_fate = $_[3];
    my $reads_without_hits = $_[4];
    my $centrifuge_hit_info = $_[5];
    my $output_dir = $_[6];
    my $lineage = $_[7];
    open(my $f1, ">$output_dir/read_fate.tsv") or die;
    print($f1 "read_id\tfinal_by_genom_cov\thierarchical_read_filt\tbwa\tcentrifuge\tfate_codes\tspecies\n");
    foreach my $i (keys%{$reads_without_hits}){
	my $centrifuge_note = get_centrifuge_results($centrifuge_fate, $centrifuge_hit_info, $i);
	if (length($centrifuge_note)==0) {
	    print($f1 "$i\t\t\t\t$centrifuge_note\t1\t\n");
	}else{
	    print($f1 "$i\t\t\t\t$centrifuge_note\t2\t\n");
	}
	
    }
    foreach my $i (keys%{$hier_read_info}){
	my $hier = $hier_read_info->{$i};
	my @bwa_hits = keys%{$bwa_fate->{$i}};
	my $bwa_note = join(",", @bwa_hits);
	my $centrifuge_note = get_centrifuge_results($centrifuge_fate, $centrifuge_hit_info, $i);
	my $final = $opt_read_info ->{$i};
	my $species = '';
	if (exists($lineage->{$final})) {
	    $species = $lineage -> {$final}->[6];
	    $species = substr($species, 3, length($species)-3);
	}
	
	if ($hier ne $bwa_note) {
	    if ($final ne $hier) {
		print($f1 "$i\t$final\t$hier\t$bwa_note\t$centrifuge_note\t3\t$species\n");
	    }else{
		print($f1 "$i\t$final\t$hier\t$bwa_note\t$centrifuge_note\t4\t$species\n");
	    }
	}else{
	    if ($final ne $hier) {
		print($f1 "$i\t$final\t$hier\t$bwa_note\t$centrifuge_note\t5\t$species\n");
	    }else{
		print($f1 "$i\t$final\t$hier\t$bwa_note\t$centrifuge_note\t6\t$species\n");
	    }
	}
	
	
    }
    close $f1;
    
}

sub get_centrifuge_results{
    my $centrifuge_fate = $_[0];
    my $centrifuge_hit_info = $_[1];
    my $centrifuge_note='';
    my $read_id = $_[2];
    if (exists($centrifuge_fate->{$read_id})) {
	    my %hits = %{$centrifuge_fate->{$read_id}};
	    my @hits;
	    foreach my $j (keys%hits){
		if ($j =~ /\.|^zymo/) {
		    my $taxid = $centrifuge_hit_info->{$j}->{'taxid'};
		    push(@hits, $taxid);
		}else{
		    push(@hits, $j);
		}
		
	    }
	    $centrifuge_note = join(",", @hits);
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

sub get_centrifuge_hits{
    my $hit_table = shift;
    my $seq_type = shift;
    my $nodes = shift;
    my $names = shift;
    my $phylogeny = shift;
    open(my $f1, "<$hit_table") or die;
    my $title = <$f1>;
    my %read_info;# read info
    my %lineage;
    my %hit_info; # list of reads uniquely assigned to each hit
    my %reads_above_strain;
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
	# record the lineage information
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
	# enrol previous read infor if this line belongs to a new read.
	my $read_id = $coln[0];#use the array index as the read id
	
	if ($acc_no =~ /\.|^zymo/) {
	    $hit_info{$acc_no}{'reads'}{$read_id}='strain';
	    $hit_info{$acc_no}{'taxid'}=$taxid;
	    $read_info{$read_id}{$acc_no}='Centrifuge.unique';
	}else{
	    $reads_above_strain{$taxid}{$read_id}=$acc_no;
	}
	
    }
    
    #assign reads without strain-level resolution to existing strains.
    my @taxid_reassigned;
    foreach my $i (keys%hit_info){
	my $taxid = $hit_info{$i}{'taxid'};
	my $node_lineage = get_lineage_nodes($taxid, $nodes);
	foreach my $j (@{$node_lineage}){
	    if (exists($reads_above_strain{$j})) {
		push(@taxid_reassigned, $j);
		my %reads = %{$reads_above_strain{$j}};
		foreach my $m (keys%reads){
		    $hit_info{$i}{'reads'}{$m}='above_strain_reassigned';
		    $read_info{$m}{$i}='Centrifuge_above_strains_reassigned';
		}
		
	    }
	    
	}
    }
    # delete the ones that have been reassigned
    foreach my $i (@taxid_reassigned){
	delete($reads_above_strain{$i});
    }
    # for reads of taxid that cannot be assigned to a existing strain.
    foreach my $i (keys%reads_above_strain){
	my %reads = %{$reads_above_strain{$i}};
	foreach my $m(keys%reads){
	    $hit_info{$i}{'reads'}{$m}='cannot assigned to strain';
	    $hit_info{$i}{'taxid'}='above_strain';
	    $read_info{$m}{$i}='cannot assigned to strain';
	}
    }
       
    return (\%hit_info, \%read_info, \%lineage);
}


sub prepare_ref_genome{
    my $dir = shift;
    my $hit_info = shift;
    my $acc2genome_file=shift;
    my %genome_included;
    my %genome_size;
    my %genome_info;
    my $ref_genome_read_count_cutoff = 10;
    system("mkdir $dir/ref_genomes");
    open(my $f1, ">$dir/ref_genomes/ref.fasta") or die;
    foreach my $i (keys%{$hit_info}){
	if (length($i)>0) {
		my $genome_id;
		if ($i =~ /^zymo/) {
		    my @name = split(/s[0-9]/, $i);
		    $genome_id = $name[0].'.fasta.gz';
		}else{
		    if (exists($acc2genome_file->{$i})) {
			$genome_id = $acc2genome_file->{$i};
		    }else{
			print("error: the genome of $i was not found!\n");
			next;
		    }
		}
		my $taxid = $hit_info->{$i}->{'taxid'};
		$genome_info{$genome_id}{'taxid'}=$taxid;
		my $read_counts = keys(%{$hit_info->{$i}->{'reads'}});
		$genome_info{$genome_id}{'reads'}+=$read_counts;
		
		
	    
        }
    }
    
    foreach my $genome_id (keys%genome_info){
	if ($genome_info{$genome_id}{'reads'}<$ref_genome_read_count_cutoff) { # exclude genomes of small number of hits
	    next;
	}
	my $taxid = $genome_info{$genome_id}{'taxid'};
	#system("aws s3 cp s3://zymo-files/WGS_Pipeline/shotgun_database/20200327_genomes/$genome_id $dir/ref_genomes/$taxid.fasta.gz");
	system("cp /home/stdell/Desktop/shotgun_pipeline/ubuntu/database/zymo_centrifuge/fasta/$genome_id $dir/ref_genomes/$taxid.fasta.gz");
	system("gunzip -c $dir/ref_genomes/$taxid.fasta.gz > $dir/ref_genomes/$taxid.fasta");
	system("rm $dir/ref_genomes/$taxid.fasta.gz");
	open(my $f2, "<$dir/ref_genomes/$taxid.fasta") or die;
	open(my $f3, ">$dir/ref_genomes/$taxid.cat.fasta") or die;
	print($f1 ">$taxid\n");
	print($f3 ">$taxid\n");
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
	system("mv $dir/ref_genomes/$taxid.cat.fasta $dir/ref_genomes/$taxid.fasta");
	$genome_size{$taxid}=$size;
	
    }
    
    close $f1;
    return (\%genome_size);
}

sub filter_reassign_reads{
    my $read_info = $_[0];
    my $hit_info = $_[1];
    my $hit_fate;
    my %opt_hit_info;
    my %opt_read_info;
    my $rank = 1;
    while (scalar(keys(%{$hit_info}))>0) {
	my $highest = get_highest_hit($hit_info);
	$hit_fate->{$highest} = $rank;
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
	    if (scalar(keys%reads)==0) {
		delete($hit_info->{$x});
		$hit_fate->{$x} = 'reads_reassigned';
	    }
	    
	}
    }
    
    return (\%opt_hit_info, \%opt_read_info, $hit_fate);
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
    my %hit_records;
    my @headers;
    my $tot_bp = 0;
    my $reads = 0;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line !~ /^@/) {
	    my @coln = split(/\t/, $line);
	    if ($coln[2] ne '*') {
		    my $taxid = $coln[2];
		    push(@{$hit_info{$taxid}{$coln[0]}}, $coln[1]);
		    $read_info{$coln[0]}{$taxid}++;
		    push(@{$hit_records{$taxid}{$coln[0]}}, $line);
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
    
    my $ambi_reads_reassinged = resolve_ambi_reads(\%hit_info, \%read_info); #temperally assigned the ambiguous reads to potential hits by random based on hit abundance
    my $average_read_size = $tot_bp/$reads;
    system("mkdir $dir/sam_files");
    my $header = join("\n", @headers);
    foreach my $i (keys%hit_records){
	open(my $f1, ">$dir/sam_files/$i.sam") or die;
	print($f1 "$header\n");
	my %records = %{$hit_records{$i}};
	foreach my $j (keys%records){
	    my %hits = %{$read_info{$j}};
	    if (scalar(keys%hits)>1) {
		my $hit_reassigned = $ambi_reads_reassinged->{$j};
		if ($hit_reassigned ne $i) {
		    next;
		}
	    }
	    my @record = @{$hit_records{$i}{$j}};
	    
	    @record = resolve_repeats(\@record);
	    
	    foreach my $line (@record){
		
		print($f1 "$line\n");
	    }
	}
	close $f1;
    }
    
    return(\%hit_info, \%read_info, \%reads_without_hits, $average_read_size);
    
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
    if (scalar(@fwd)>1) {
	my $m=0;
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


sub resolve_ambi_reads{
    my $new_hit_info = $_[0];
    my $new_read_info = $_[1];
    my %ambi_reads;
    my %hit_abun;
    foreach my $i (keys%{$new_hit_info}){
	$hit_abun{$i} = scalar(keys%{$new_hit_info->{$i}});
	foreach my $j (keys%{$new_hit_info->{$i}}){
	    my %hits = %{$new_read_info->{$j}};
	    if (scalar(keys%hits)>1) {
		$ambi_reads{$j}++;
	    }
	}
    }
    
    foreach my $i (keys%ambi_reads){
	my %hits = %{$new_read_info->{$i}};
	my $hit_assigned = assign_by_abun(\%hit_abun, \%hits);
	$ambi_reads{$i} = $hit_assigned;
    }
    return (\%ambi_reads);
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


#1;