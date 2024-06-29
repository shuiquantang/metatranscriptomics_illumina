#!/usr/bin/perl -I /media/sf_VirtualBox/scripts/metaillu/scripts/bin
use strict;
use warnings;
use Inputs;
#use NCBITaxonomy;
#use SummarizeTaxa;
package IlluminaAnalysis;
use Parallel::Loops;

sub TrimAndFilter {
    #$rawdata, $trim_folder, $sample_info, $bin_path, $threads
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $bin_path = $_[3];
    my $reference_db = $_[4];
    my $threads_RNA = $_[5];
    my $QC_folder = $_[6];
    my $trim_quality = 15;
    my $min_size = 60;
    my $cpu = `nproc`;
    chomp $cpu;
    
    # run Trimmomatic
    my $max_memory = `free -m | awk 'FNR == 2 {print \$2}'`; # max memory available in the system in Mb
    chomp($max_memory);
    my @samples = keys%{$sample_info};
    my $pl = Parallel::Loops->new($cpu);
    my $processes = int($max_memory/1000/10);
    if ($processes>scalar(@samples)) {
	$processes = scalar(@samples);
    }elsif($processes==0){
	$processes = 1;
    }
    my $threads_per_process = int($cpu/$processes);
    if ($threads_per_process<1) {
	$threads_per_process = 1;
    }elsif($threads_per_process>5){
	$threads_per_process=5
    }
    
    
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (keys%{$sample_info}){
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1 = $sample_info->{$i}->{'R1'};
	my $adapter_seq_file = 'NexteraPE-PE.fa';
	my $adapter_file = "$QC_folder/$i\_adapter.txt";
	if (-e $adapter_file) {
	    open(my $f1, "<$adapter_file") or die;
	    my $line = <$f1>;
	    
	    if ($line =~ /truseq/) {
		$adapter_seq_file = 'TruSeq3-PE-2.fa';
	    }
	    close $f1;
	}
	if ($seq_type =~ /\.pe/) {
	    my $R2 = $sample_info->{$i}->{'R2'};
	    system("java -jar $bin_path/Trimmomatic/trimmomatic-0.33.jar PE -threads $threads_per_process -phred33 ".
		   "$input_dir/$R1 $input_dir/$R2 $output_dir/$i\_R1.paired.fastq.gz $output_dir/$i\_R1.unpaired.fastq.gz $output_dir/$i\_R2.paired.fastq.gz $output_dir/$i\_R2.unpaired.fastq.gz ".
		   "ILLUMINACLIP:$bin_path/Trimmomatic/$adapter_seq_file:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$trim_quality MINLEN:$min_size >/dev/null 2>$output_dir/$i.log");
	    system("cat $output_dir/$i\_R1.paired.fastq.gz $output_dir/$i\_R1.unpaired.fastq.gz $output_dir/$i\_R2.paired.fastq.gz $output_dir/$i\_R2.unpaired.fastq.gz > $output_dir/$i.TM.fastq.gz");
	    system("rm $output_dir/$i\_R1.paired.fastq.gz $output_dir/$i\_R1.unpaired.fastq.gz $output_dir/$i\_R2.paired.fastq.gz $output_dir/$i\_R2.unpaired.fastq.gz");
	}else{
	    system("java -jar $bin_path/Trimmomatic/trimmomatic-0.33.jar SE -threads $threads_per_process -phred33 ".
		   "$input_dir/$R1 $output_dir/$i.TM.fastq.gz ".
		   "ILLUMINACLIP:$bin_path/Trimmomatic/$adapter_seq_file:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:8:$trim_quality MINLEN:$min_size >/dev/null 2>$output_dir/$i.log");
	}
	
    #}
    });
    chdir($output_dir);
    
    
    # run ribodetector
    my $chunk_size = int($max_memory/1000/$cpu)*512;
    if ($chunk_size == 0 ) {
	$chunk_size = 256;
    }
    my $cpu_per_process = int($cpu/$threads_RNA);
    $pl = Parallel::Loops->new($threads_RNA);
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (keys%{$sample_info}){
	#system("ribodetector -t $cpu -l 80 -i $output_dir/$i.TM.fastq.gz -o $output_dir/$i.fastq.gz --chunk_size $chunk_size");
	my $cmd = "ribodetector_cpu -t $cpu_per_process -l 80 -i $output_dir/$i.TM.fastq.gz -o $output_dir/$i.RD.fastq.gz --chunk_size $chunk_size";
	my $log = "$output_dir/$i.log";
	Inputs::print_and_execute($cmd, $log);
    #}
    });
    $pl = Parallel::Loops->new($cpu);
    my %tot_reads;
    $pl->share(\%tot_reads);
    $pl -> foreach (\@samples, sub{
	my $i = $_;
	my $trimmomatic_survived_reads = `gunzip -c $output_dir/$i.TM.fastq.gz | wc -l`;
	my $non_RNA_reads = `gunzip -c $output_dir/$i.RD.fastq.gz | wc -l`;
	chomp $trimmomatic_survived_reads; chomp $non_RNA_reads;
	my $RNA_reads = $trimmomatic_survived_reads-$non_RNA_reads;
	$trimmomatic_survived_reads = $trimmomatic_survived_reads/4;
	$RNA_reads = $RNA_reads/4;
	$non_RNA_reads = $non_RNA_reads/4;
	$tot_reads{$i}="$trimmomatic_survived_reads\-$RNA_reads\-$non_RNA_reads";
        system("rm $output_dir/$i.TM.fastq.gz");
    #}
    });
    
    
    # use sourmash to quickly determine the host ID
    $processes = int($max_memory/1000/8);
    $pl = Parallel::Loops->new($processes);
    my %host_assembly_id;
    $pl->share(\%host_assembly_id);
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (@samples){
	my $input_seqs = "$i.RD.fastq.gz";
        my $host_id = predict_host_genome($input_seqs, $reference_db, $i);
	$host_assembly_id{$i}=$host_id;
    #}
    });
    
    #use Kraken2 to filter out host reads.
    
    $processes = int($max_memory/1000/9);
    if ($processes>scalar(@samples)) {
	$processes = scalar(@samples);
    }
    
    $pl = Parallel::Loops->new($processes);
    $threads_per_process = int($cpu/$processes);
    
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (keys%{$sample_info}){
	my $host_id = $host_assembly_id{$i};
	if (length($host_id)>0) {
	    host_dna_filtration($i, $reference_db, $host_id, $threads_per_process);
	}
	
    #}
    });
    
    chdir("..");
    my %read_length;
    $pl = Parallel::Loops->new($cpu);
    $pl -> share(\%read_length);
    $pl -> foreach (\@samples, sub{
	my $i = $_;
	my $ave_read_length=0;
    #foreach my $i (@samples){
	if (-e "$output_dir/$i.HDF.fastq") {
	    $ave_read_length = rename_fastq("$output_dir/$i.HDF.fastq", "$output_dir/$i.RN.fastq.gz", "$output_dir/$i.RN.fasta");
	    system("rm $output_dir/$i.RD.fastq.gz $output_dir/$i.HDF.fastq");
	}else{
	    $ave_read_length = rename_fastq("$output_dir/$i.RD.fastq.gz", "$output_dir/$i.RN.fastq.gz", "$output_dir/$i.RN.fasta");
	    system("rm $output_dir/$i.RD.fastq.gz");
	}
	my $final_fastq = "$output_dir/$i.fastq.gz";
	my $low_diversity_reads = remove_low_diversity_reads($output_dir, "$output_dir/$i.RN.fastq.gz", "$output_dir/$i.RN.fasta", $i, $bin_path, $final_fastq);
	system("rm $output_dir/$i.RN.fastq.gz $output_dir/$i.RN.fasta");
	$read_length{$i}="$ave_read_length\t$low_diversity_reads";
    #}
    });
    #remove low diversity reads
    
    
    
    open(my $f1, ">$output_dir/summary.tsv") or die;
    print($f1 "internal_id\traw_reads\tpaired_end\tboth_surviving(%)\tforward_only(%)\treverse_only(%)\tboth_dropped(%)\trRNA\thost_reads(%)\tread_length(bp)\tlow_diversity_reads\n");
    
    foreach my $i (@samples){
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1 = $sample_info->{$i}->{'R1'};
	my $survived_reads = $tot_reads{$i};
	my @reads = split(/\-/, $survived_reads);
	my $RNA_perc = int($reads[1]/$reads[0]*10000)/100;
	my ($TM_numbers, $host_reads_perc) = read_filtration_log("$output_dir/$i.log");
	my @info = split(/\t/, $read_length{$i});
	if ($seq_type =~ /\.pe/) {
	    print($f1 "$i\t$TM_numbers->[0]\tPE\t$TM_numbers->[2]\t$TM_numbers->[4]\t$TM_numbers->[6]\t$TM_numbers->[8]\t$RNA_perc\t$host_reads_perc\t$info[0]\t$info[1]\n");
	}else{
	    print($f1 "$i\t$TM_numbers->[0]\tSE\t0%\t$TM_numbers->[2]\t0%\t$TM_numbers->[4]\t$RNA_perc\t$host_reads_perc\t$info[0]\t$info[1]\n");
	}
	
    }
    close $f1;
}

sub rename_fastq{
    my $input_file = shift;
    my $output_file = shift;
    my $fasta = shift;
    my $f1;
    if ($input_file =~ /.gz$/) {
	open($f1, "gzip -dc $input_file |") or die;
    }else{
	open($f1, "<$input_file") or die;
    }
    open(my $f2, '|-', "gzip> $output_file") or die;
    open(my $f3, ">$fasta") or die;
    my $num=0;
    my $bp = 0;
    while (my $L1 = <$f1>) {
	chomp $L1;
	my $L2 = <$f1>;
	my $L3 = <$f1>;
	my $L4 = <$f1>;
	$bp+=(length($L2)-1);
	my $record = "\@$num\n".$L2.$L3.$L4;
	print ($f2 $record);
	print($f3 "\@$num\n".$L2);
	$num++;
    }
    close $f1;
    close $f2;
    close $f3;
    my $ave_length = 0;
    if ($num>0) {
	$ave_length = int($bp/$num);
    }
    return($ave_length);
}


sub read_filtration_log{
    my $log = shift;
    my @TM_results;
    my $host_reads_perc = 0;
    open(my $f1, "<$log") or die;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^Input Read/) {
	    $line =~ s/NaN/0/g;
	    my @words = split(/\ |\(|\)/, $line);
	    foreach my $j (@words){
		if ($j =~ /\d/) {
		    push(@TM_results, $j);
		}
	    }
	    if ($TM_results[0] == 0) {
		$TM_results[0]=1;
	    }
	}elsif($line =~ /^\s*[0-9]+ sequences classified/){
	    my @words = split(/\(|\)/, $line);
	    $host_reads_perc = $words[1];
	}
		
    }
	    
    close $f1;
    return (\@TM_results, $host_reads_perc);
}

#host_dna_filtration($input_seqs, $survived_seqs, $host_id);
sub host_dna_filtration{
    my $id = shift;
    my $reference_db = shift;
    my $host_id = shift;
    my $threads = shift;
    my $input_seqs = "$id.RD.fastq.gz";
    my $survived_seqs = "$id.HDF.fastq";
    my $kraken2_output = "$id.kraken2.output.txt";
    my $log = "$id.log";
    my $db = "$reference_db/kraken2/$host_id";
    my $cmd = "kraken2 --gzip-compressed --db $db --threads $threads --unclassified-out $survived_seqs --output $kraken2_output $input_seqs 2>>$log";
    Inputs::print_and_execute($cmd, $log);
    system("rm -r $kraken2_output");
}

#$host_id = predict_host_genome($input_seqs);
sub predict_host_genome{
    my $input_seqs=shift;
    my $reference_db = shift;
    my $sample_id = shift;
    my $log = "$sample_id.log";
    my $sourmash_database = "$reference_db/sourmash";
    #take a subset of the input fastq file
    my $subset_file = "$sample_id.subset.fastq.gz";
    my $subset_size = 1000000*4;
    my $cmd = "gunzip -c $input_seqs | head -$subset_size | gzip > $subset_file 2>>$log";
    Inputs::print_and_execute($cmd, $log);
    #Build match with reference databases
    my $host_genomes = "$sourmash_database/host_genomes.k51.zip";
    my $host_genome_tax = "$sourmash_database/lineage/host_genomes_tax.csv";
    my $threshold_bp = 50000; #50000 by default
    my $fastq_sig = "$sample_id.sig.gz";
    $cmd = "sourmash sketch dna -p k=51,scaled=1000,abund $subset_file -o $fastq_sig --name $sample_id 2>&1 > /dev/null";
    Inputs::print_and_execute($cmd, $log);
    system("rm $subset_file");
    my $sketch_output = "$sample_id.sketch.csv";
    $cmd = "sourmash gather $fastq_sig $host_genomes --dna --ksize 51 --threshold-bp $threshold_bp -o $sketch_output>>$log";
    Inputs::print_and_execute($cmd, $log);
    my $host_assmebly_id = '';
    if (-e "$sketch_output") {
	open(my $f1, "<$sketch_output") or die;
	my $header = <$f1>;
	my $best_hit = <$f1>;
	my @coln = split(/\,/, $best_hit);
	if (@coln>10) {
	    $host_assmebly_id = $coln[9];
	}
	#system("rm $sketch_output");
    }
    system("rm $fastq_sig");
    return($host_assmebly_id);
}


# subroutine to read input parameters, obtained on-line
sub SubsetRawdata {
    # $rawdata, $QC_folder, $sample_info, $bin_path
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $subset_size = $_[3];
    my $cpu = `nproc`;
    chomp $cpu;
    my @raw_seq_files;
    my $pl = Parallel::Loops->new($cpu);
    my @samples = sort(keys%{$sample_info});
    my %tot_reads;
    $pl -> share(\%tot_reads);
    $pl -> foreach (\@samples, sub{
    #foreach my $i (keys%{$sample_info}){
        my $i =$_;
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1 = $sample_info->{$i}->{'R1'};
	push (@raw_seq_files, "$input_dir/$R1");
	my $R2 = $sample_info->{$i}->{'R2'};
	$subset_size = $subset_size*4;
	system("gunzip -c $input_dir/$R1 | head -$subset_size | gzip > $output_dir/$R1");
	system("gunzip -c $input_dir/$R2 | head -$subset_size | gzip > $output_dir/$R2");
	my $tot_reads = `gunzip -c $input_dir/$R1 | wc -l`;
	chomp $tot_reads;
	$tot_reads{$i}=$tot_reads/4;
    #}
    });
    
    foreach my $i (@samples){
	$sample_info->{$i}->{'tot_reads'} = $tot_reads{$i};
    }
    
}

# subroutine to read input parameters, obtained on-line
sub FastQC {
    # $rawdata, $QC_folder, $sample_info, $bin_path
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $bin_path = $_[3];
    my $cpu = `nproc`;
    chomp $cpu;
    my @raw_seq_files;
    my $subset_size = 1000000;
    my @samples = sort(keys%{$sample_info});
    my $pl = Parallel::Loops->new(int($cpu/5));
    #foreach my $i (@samples){
    $pl -> foreach (\@samples, sub{
        my $i =$_;
	my $R1 = $sample_info->{$i}->{'R1'};
	system("gunzip -c $input_dir/$R1 | head -$subset_size | gzip > $output_dir/$R1");
	system("perl $bin_path/FastQC/fastqc -o $output_dir -t 5 --quiet $output_dir/$R1");
	my $adapter = get_adapter_from_fastqc($i, $R1, $output_dir);
	system("rm -r $output_dir/$R1");
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	if ($seq_type =~ /\.pe/) {
	    my $R2 = $sample_info->{$i}->{'R2'};
	    system("gunzip -c $input_dir/$R2 | head -$subset_size | gzip > $output_dir/$R2");
	    system("perl $bin_path/FastQC/fastqc -o $output_dir -t 5 --quiet $output_dir/$R2");
	    system("rm -r $output_dir/$R2");
	}
    });
    #}
    system("rm $output_dir/*fastqc.zip");
}

sub remove_low_diversity_reads{
    my $output_dir=shift;
    my $fastq = shift;
    my $fasta = shift;
    my $sample_id = shift;
    my $bin_path = shift;
    my $final_fastq = shift;
    my $log = "$output_dir/$sample_id.log";
    my $cpu = `nproc`;
    chomp $cpu;
    my $cmd = "$bin_path/sdust/sdust -w 64 -t $cpu $fasta > $output_dir/$sample_id\_low_complexity_reads.txt";
    Inputs::print_and_execute($cmd, $log);
    my %low_complex_reads;
    read_sdust_output("$output_dir/$sample_id\_low_complexity_reads.txt", \%low_complex_reads);
    my @low_complex_reads = sort {$a <=> $b} keys(%low_complex_reads);
    my $size = scalar(@low_complex_reads);
    open(my $f1, "gzip -dc $fastq |") or die;
    open(my $f2, '|-', "gzip> $final_fastq") or die;
    while (my $r1 = <$f1>) {
	my $r2 = <$f1>; my $r3 = <$f1>; my $r4 = <$f1>;
	my $id = $r1; chomp $id;
	if (scalar(@low_complex_reads)>0) {
	    my $query = $low_complex_reads[0];
	    if ($id eq "\@$query") { #skip low diversity reads
		shift(@low_complex_reads);
	    }else{
		print($f2 "$r1$r2$r3$r4");	
	    }
	}else{
	    print($f2 "$r1$r2$r3$r4");
	}
    }
    close $f1;
    close $f2;
    return($size);
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

sub get_adapter_from_fastqc{
    my $internal_id = shift;
    my $R1=shift;
    my $dir = shift;
    my @id = split(/\./, $R1);
    pop@id; pop@id;
    my $new_id = join(".", @id);
    my $fastqc_zip = "$dir/$new_id\_fastqc.zip";
    system("unzip -o -d $dir $fastqc_zip");
    open(my $f1, "<$dir/$new_id\_fastqc/fastqc_data.txt") or return();
    my @adapters;
    my @percentage;
    my $recording_controller=0;
    while (my $line = <$f1>) {
	chomp $line;
	if ($line =~ /^>>Adapter Content/) {
	    $recording_controller = 1;
	}elsif($line =~ /^>>/){
	    $recording_controller = 0;
	}else{
	    if ($recording_controller == 1) {
		my @coln = split(/\t/, $line);
		if ($line =~ /^#Position/) {
		    @adapters = @coln;
		}else{
		    if (scalar(@coln) eq scalar(@adapters)) {
			@percentage = @coln;
		    }
		}
		
	    }
	    
	}
	
    }
    system("rm -r $dir/$new_id\_fastqc/");
    my $nextera_perc = 0;
    my $truseq_perc = 0;
    for (my $i=1; $i<scalar(@adapters); $i++){
	my $adapter = $adapters[$i];
	my $perc = $percentage[$i];
	if ($adapter =~ /Illumina Universal Adapter/) {
	    $truseq_perc = $perc;
	}elsif($adapter =~ /Nextera/){
	    $nextera_perc = $perc;
	}
	
    }
    open(my $f2, ">$dir/$internal_id\_adapter.txt") or die;
    if ($truseq_perc > $nextera_perc ) {
	print($f2 "truseq");
    }else{
	print($f2 "nextera");
    }
    close $f2;
    
}
1;
