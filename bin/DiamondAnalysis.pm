#!/usr/bin/perl -I /home/stdell/Desktop/shotgun_pipeline/ubuntu/scripts/bin
use strict;
use warnings;
use NCBITaxonomy;
package DiamondAnalysis;
use Parallel::Loops;

sub AMR_Diamond {
    my $trim_dir = $_[0];
    my $input_dir = $_[1];
    my $output_dir = $_[2];
    my $sample_info = $_[3];
    my $species_abun = $_[4];
    my $threads = $_[5];
    my $cpu = `nproc`;
    chomp $cpu;
    mkdir("$output_dir/raw_results");
    mkdir("$output_dir/summary");
    my $pl = Parallel::Loops->new($threads);
    my @samples = sort(keys%{$sample_info});
    $pl -> foreach (\@samples, sub{
    #foreach my $i (keys%{$sample_info}){
	my $i = $_;
	my $R1 = "$trim_dir/$i.fastq.gz";
	my $command = "diamond blastx -d /diamond_db/NCBIATB.dmnd -q $R1 --outfmt 6 stitle qseqid slen -o $output_dir/raw_results/$i.results.txt ".
	       "--max-target-seqs 10 -e 0.000000000000001 --id 80 --query-cover 90 --threads $cpu --quiet >/dev/null 2>&1";
	system($command);
    #}
    });
    foreach my $i (@samples){
	my %abun;
	if (exists($species_abun->{$i})) {
	    %abun = %{$species_abun->{$i}};
	}
	summarize_results("$output_dir/raw_results/$i.results.txt",  "$input_dir/$i/read2hit.tsv", "atb", "$output_dir/summary/$i.summary.txt");
    }
    #system("rm -r $output_dir/raw_results/");
}

sub viru_Diamond {
    my $trim_dir = $_[0];
    my $input_dir = $_[1];
    my $output_dir = $_[2];
    my $sample_info = $_[3];
    my $species_abun = $_[4];
    my $threads = $_[5];
    my $cpu = `nproc`;
    chomp $cpu;
    mkdir("$output_dir/raw_results");
    mkdir("$output_dir/summary");
    my $pl = Parallel::Loops->new($threads);
    my @samples = sort(keys%{$sample_info});
    $pl -> foreach (\@samples, sub{
	my $i = $_;
    #foreach my $i (keys%{$sample_info}){
	my $R1 = "$trim_dir/$i.fastq.gz";
	my $command= "diamond blastx -d /diamond_db/VFDBS.dmnd -q $R1 --outfmt 6 stitle qseqid slen -o $output_dir/raw_results/$i.results.txt ".
	       "--max-target-seqs 10 -e 0.000000000000001 --id 80 --query-cover 90 --threads $cpu --quiet >/dev/null 2>&1";
	system($command);
    #}
    });
    foreach my $i (@samples){
	my %abun;
	if (exists($species_abun->{$i})) {
	    %abun = %{$species_abun->{$i}};
	}
	summarize_results("$output_dir/raw_results/$i.results.txt", "$input_dir/$i/read2hit.tsv", 'vf', "$output_dir/summary/$i.summary.txt");
    }
    #system("rm -r $output_dir/raw_results/");
    
}

sub summarize_results{
    my $diamond_file = shift;
    my $read2hit_file = shift;
    my $tag=shift;
    my $output_file=shift;
    my %read2taxa;
    if (-e $read2hit_file) {
	open(my $f1, "<$read2hit_file");
	while (my $line = <$f1>) {
	    chomp $line;
	    my @coln = split(/\t/, $line);
	    $read2taxa{$coln[0]}=$coln[3];
	}
	close $f1;
    }
    my %gene_abun;
    open(my $f2, "<$diamond_file") or return;
    my @records;
    my $last_read = '';
    my %taxa_abun;
    while (my $line = <$f2>) {
	chomp $line;
	my @coln = split(/\t/, $line);
	if (length($last_read)>0 && $last_read ne $coln[1]) {
	    my $gene=get_concensus(\@records, $tag);
	    my $taxa='NA';
	    if (exists($read2taxa{$last_read})) {
		$taxa = $read2taxa{$last_read};
	    }
	    $gene_abun{$taxa}{$gene}{'abun'}++;
	    $taxa_abun{$taxa}++;
	    @records=();
	    
	}
	push(@records, $line);
	$last_read = $coln[1];
    }
    close $f2;
    my $gene=get_concensus(\@records, $tag);
    my $taxa='NA';
    if (exists($read2taxa{$last_read})) {
	$taxa = $read2taxa{$last_read};
    }
    $gene_abun{$taxa}{$gene}{'abun'}++;
    $taxa_abun{$taxa}++;
    
    open(my $f3, ">$output_file") or die;
    print($f3 "species\tprotein\tread_counts\n");
    my @taxa = sort{$taxa_abun{$b}<=>$taxa_abun{$a}}keys%gene_abun;
    foreach my $i (@taxa){
	foreach my $j (keys%{$gene_abun{$i}}){
	    my $abun = $gene_abun{$i}{$j}{'abun'};
	    print($f3 "$i\t$j\t$abun\n");
	    
	}
    }
    close $f3;
    
    
}

sub get_concensus{
    my $records = $_[0];
    my $tag = $_[1];
    my %names;
    
    if ($tag eq 'vf') {
	foreach my $i (@{$records}){
	    my @arr = split(/\t/, $i);
	    $i = $arr[0];
            if ($i =~ /   /) {
                my @coln = split(/   /,$i);
                if (scalar(@coln)>1) {
                    $names{$coln[1]}++;
                }
            }else{
                $names{$i}++;
            }
            
	    
	}
	
    }elsif($tag eq 'atb'){
	foreach my $i (@{$records}){
	    my @arr = split(/\t/, $i);
	    $i = $arr[0];
            my $new='';
            if ($i =~ / \[/) {
                my @coln = split(/ \[/,$i);
                @coln = split(/\ /, $coln[0]);
                if ($i =~ /MULTISPECIES/) {
                    shift(@coln);
                    shift(@coln);
                }else{
                    shift(@coln);
                }
                $new = join(" ", @coln);
            }else{
                $new = $i;
            }
            
	    
	    $names{$new}++;
	}
    }
    my @name = sort{$names{$b}<=>$names{$a}}keys%names;
    if (scalar(@name)==0) {
        return('');
    }else{
        return($name[0]);
    }
}

1;
