#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/bin
use strict;
use warnings;
#use NCBITaxonomy;
#use SummarizeTaxa;
package IlluminaAnalysis;
use Parallel::Loops;

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
    my @samples = keys%{$sample_info};
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
    my @samples = keys%{$sample_info};
    my $pl = Parallel::Loops->new(int($cpu/5));
    #foreach my $i (@samples){
    $pl -> foreach (\@samples, sub{
        my $i =$_;
	my $R1 = $sample_info->{$i}->{'R1'};
	system("gunzip -c $input_dir/$R1 | head -$subset_size | gzip > $output_dir/$R1");
	system("perl $bin_path/FastQC/fastqc -o $output_dir -t 5 --quiet $output_dir/$R1");
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


sub trim_and_filter {
    #$rawdata, $trim_folder, $sample_info, $bin_path, $threads
    my $input_dir = $_[0];
    my $output_dir = $_[1];
    my $sample_info = $_[2];
    my $bin_path = $_[3];
    my $threads_trimmomatic = $_[4];
    my $threads_SortMeRNA = $_[5];
    my $sortmerna_database = $_[6];
    my $trim_quality = 15;
    my $min_size = 60;
    my $cpu = `nproc`;
    chomp $cpu;
    my $summary_file = "trimming.summary.txt";
    open(my $f1, ">$output_dir/summary.tsv") or die;
    print($f1 "internal_id\traw_reads\tpaired_end\tboth_surviving(%)\tforward_only(%)\treverse_only(%)\tdropped(%)\trRNAs_removed(%)\tfinal_reads(%)\n");
    # run Trimmomatic
    my $pl = Parallel::Loops->new($threads_trimmomatic);
    my @samples = keys%{$sample_info};
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (keys%{$sample_info}){
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1 = $sample_info->{$i}->{'R1'};
	if ($seq_type =~ /\.pe/) {
	    my $R2 = $sample_info->{$i}->{'R2'};
	    system("java -jar $bin_path/Trimmomatic/trimmomatic-0.33.jar PE -threads $cpu -phred33 ".
		   "$input_dir/$R1 $input_dir/$R2 $output_dir/$i\_R1.paired.fastq.gz $output_dir/$i\_R1.unpaired.fastq.gz $output_dir/$i\_R2.paired.fastq.gz $output_dir/$i\_R2.unpaired.fastq.gz ".
		   "ILLUMINACLIP:$bin_path/Trimmomatic/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:$trim_quality MINLEN:$min_size >/dev/null 2>$output_dir/$i.log");
	    system("cat $output_dir/$i\_R1.paired.fastq.gz $output_dir/$i\_R1.unpaired.fastq.gz $output_dir/$i\_R2.paired.fastq.gz $output_dir/$i\_R2.unpaired.fastq.gz > $output_dir/$i.TM.fastq.gz");
	    system("rm $output_dir/$i\_R1.paired.fastq.gz $output_dir/$i\_R1.unpaired.fastq.gz $output_dir/$i\_R2.paired.fastq.gz $output_dir/$i\_R2.unpaired.fastq.gz");
	}else{
	    system("java -jar $bin_path/Trimmomatic/trimmomatic-0.33.jar SE -threads $cpu -phred33 ".
		   "$input_dir/$R1 $output_dir/$i.TM.fastq.gz ".
		   "ILLUMINACLIP:$bin_path/Trimmomatic/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:8:$trim_quality MINLEN:$min_size >/dev/null 2>$output_dir/$i.log");
	}
	
    #}
    });
    
    # run SortMeRNA
    
    $pl = Parallel::Loops->new($threads_SortMeRNA);
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    # every run of sortmerna will check if index are available in the workdir, if so, it will skip indexing 
    #foreach my $i (keys%{$sample_info}){
        system("mkdir $output_dir/$i");
	system("sortmerna --workdir $output_dir/$i --threads $cpu -ref /$sortmerna_database/smr_v4.3_default_db.fasta --idx-dir /$sortmerna_database/idx -reads $output_dir/$i.TM.fastq.gz -aligned $i.rRNA -other $i -fastx");
	system("rm -r $output_dir/$i");
	system("mv $i.rRNA.fq.gz $output_dir/$i.rRNA.fastq.gz");
	system("mv $i.fq.gz $output_dir/$i.fastq.gz");
	system("mv $i.rRNA.log $output_dir/");
    #}
    });
    
    
    
    $pl = Parallel::Loops->new($cpu);
    my %tot_reads;
    $pl->share(\%tot_reads);
    $pl -> foreach (\@samples, sub{
	my $i = $_;
	my $trimmomatic_survived_reads = `gunzip -c $output_dir/$i.TM.fastq.gz | wc -l`;
	my $RNA_reads = `gunzip -c $output_dir/$i.rRNA.fastq.gz | wc -l`;
	my $non_RNA_reads = `gunzip -c $output_dir/$i.fastq.gz | wc -l`;
	chomp $trimmomatic_survived_reads; chomp $RNA_reads; chomp $non_RNA_reads;
	$trimmomatic_survived_reads = $trimmomatic_survived_reads/4;
	$RNA_reads = $RNA_reads/4;
	$non_RNA_reads = $non_RNA_reads/4;
	$tot_reads{$i}="$trimmomatic_survived_reads\-$RNA_reads\-$non_RNA_reads";
        system("rm $output_dir/$i.TM.fastq.gz $output_dir/$i.rRNA.fastq.gz");
    #}
    });
    
    foreach my $i (@samples){
	my $seq_type = $sample_info -> {$i}->{'seq_type'};
	my $R1 = $sample_info->{$i}->{'R1'};
	my $survived_reads = $tot_reads{$i};
	my @reads = split(/\-/, $survived_reads);
	my $RNA_perc = int($reads[1]/$reads[0]*10000)/100;
	if ($seq_type =~ /\.pe/) {
	    open(my $f2, "<$output_dir/$i.log") or die;
	    while (my $line = <$f2>) {
		if ($line =~ /^Input Read/) {
		    $line =~ s/NaN/0/g;
		    my @words = split(/\ |\(|\)/, $line);
		    my @numbers;
		    foreach my $i (@words){
			if ($i =~ /\d/) {
			    push(@numbers, $i);
			}
		    }
		    if ($numbers[0] == 0) {
			$numbers[0]=1;
		    }
		    
		    my $final_survival = int($reads[2]/$numbers[0]/2*10000)/100;
		    print($f1 "$i\t$numbers[0]\tPE\t$numbers[2]\t$numbers[4]\t$numbers[6]\t$numbers[8]\t$RNA_perc%\t$final_survival%\n");
		}
		
	    }
	    close $f2;
	}else{
	    open(my $f2, "<$output_dir/$i.log") or die;
	    while (my $line = <$f2>) {
		if ($line =~ /^Input Read/) {
		    $line =~ s/NaN/0/g;
		    my @words = split(/\ |\(|\)/, $line);
		    my @numbers;
		    foreach my $i (@words){
			if ($i =~ /\d/) {
			    push(@numbers, $i);
			}
		    }
		    if ($numbers[0] == 0) {
			$numbers[0]=1;
		    }
		    
		    my $final_survival = int($reads[2]/$numbers[0]*10000)/100;
		    print($f1 "$i\t$numbers[0]\tSE\t\t$numbers[2]\t\t$numbers[4]\t$RNA_perc%\t$final_survival%\n");
		}
		
	    }
	    
	    close $f2;
	}
	
    }
    close $f1;
}

1;
