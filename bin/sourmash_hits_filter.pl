#!/usr/bin/perl
#use strict;
#use warnings;
use vars qw($opt_i $opt_o $opt_u $opt_l $opt_L $opt_S $opt_s $opt_G $opt_g);
&Getopts('i:o:u:l:L:S:s:G:g:');
# -u unique_hits cutoff, -l min_coverage_length, -L good_coverage_length,
# -S the min abundance of strains within the same species, -s the top maximum strains allowed in the same species
# -G the min abundance of species within the same genus, -g the top maximum species allowed in the same genus
# -i the input sourmash hit table -o output_file
# -i in1717_1.k51.sketch.with-lineages.csv -o in1717_1.filtered.csv -u 5 -l 1000 -L 50000 -S 0.05 -s 5 -G 0.01 -g 15
my %lineage;
my %lineage_old;
my %strain_info;
my @title = read_csv_table($opt_i, \%strain_info);
phylogeny_summary(\%strain_info, \%lineage_old);
filter_by_coverage(\%strain_info, $opt_u, $opt_l, $opt_L);
phylogeny_summary(\%strain_info, \%lineage);
filter_by_phylogeny_cutoffs(\%lineage, \%strain_info, $opt_S, $opt_s, $opt_G, $opt_g);
output_filtered_table(\%lineage_old, \%strain_info, $opt_o, \@title);

sub output_filtered_table{
    my $lineage = shift;
    my $strain_info = shift;
    my $output_file = shift;
    my $titles = shift;
    open(my $f1, ">$output_file") or die;
    unshift(@{$titles},'status');
    my $header = join("\,", @{$titles});
    print($f1 "$header\n");
    foreach my $i (sort(keys(%{$lineage}))){
        foreach my $j (sort(keys(%{$lineage->{$i}}))){
            my %strain_abun = %{$lineage->{$i}->{$j}};
            my @strains = sort{$strain_abun{$b}<=>$strain_abun{$a}}keys%strain_abun;
            foreach my $z (@strains){
                my @coln;
                foreach my $x (@{$titles}){
                    my $item = $strain_info->{$z}->{$x};
		    if ($item =~ /\,/) {
			$item = "\"$item\"";
		    }
		    
                    push(@coln, $item);
                }
                my $line = join("\,", @coln);
                print($f1 "$line\n");
            }
        }
    }
    close $f1;
}

sub filter_by_phylogeny_cutoffs{
    my $lineage = shift;
    my $strain_info = shift;
    my $min_strain_perc_in_species = shift;
    my $max_strians_in_species = shift;
    my $min_species_perc_in_genus = shift;
    my $max_species_in_genus = shift;
    my %genus2sp;
    my %sp2strain;
    foreach my $i (keys(%{$lineage})){
        foreach my $j (keys(%{$lineage->{$i}})){
            my %strain_abun = %{$lineage->{$i}->{$j}};
            my %strain_perc_abun = abs2perc(\%strain_abun);
            my @stains = sort{$strain_perc_abun{$b}<=>$strain_perc_abun{$a}}keys(%strain_perc_abun);
            my $num = 1;
            $genus2sp{$i}{$j}=0;
            foreach my $z (@stains){
                my $abun = $strain_perc_abun{$z};
                my $abs = $strain_abun{$z};
                my $status = $strain_info->{$z}->{'status'};
                if ($status eq 'Y') {
                    if ($abun <$min_strain_perc_in_species) {
                        $strain_info->{$z}->{'status'}="less than $min_strain_perc_in_species within the same species";
                        next;
                    }
                    if ($num > $max_strians_in_species) {
                        $strain_info->{$z}->{'status'}="not top $max_strians_in_species strains within the same species";
                        $num++;
                        next;
                    }
                }
                $genus2sp{$i}{$j}+=$abs;
                push(@{$sp2strain{$j}}, $z);
                $num++;
            }
        }
    }
    foreach my $i (keys(%genus2sp)){
        my %species_abun = %{$genus2sp{$i}};
        my %species_perc_abun = abs2perc(\%species_abun);
        my @species = sort{$species_perc_abun{$b} <=> $species_perc_abun{$a}}keys(%species_perc_abun);
        my $num = 1;
        foreach my $j (@species){
            my $abun = $species_perc_abun{$j};
            if ($abun <$min_species_perc_in_genus) {
                    foreach my $z (@{$sp2strain{$j}}){
                        my $status = $strain_info->{$z}->{'status'};
                        if ($status eq 'Y') {
                            $strain_info->{$z}->{'status'}="less than $min_species_perc_in_genus within the same genus";
                        }
                        
                    }
                    next;
            }
            if ($num > $max_species_in_genus) {
                    foreach my $z (@{$sp2strain{$j}}){
                        my $status = $strain_info->{$z}->{'status'};
                        if ($status eq 'Y') {
                            $strain_info->{$z}->{'status'}="not top $max_strians_in_species species within the same genus";
                        }
                    }
                    $num++;
                    next;
            }
            $num++;
        }
    }
     
}

sub abs2perc{
    my $old = shift;
    my %new;
    my $sum = 0;
    foreach my $i (keys(%{$old})){
        my $value = $old->{$i};
        $sum+=$value;
    }
    if ($sum>0) {
        foreach my $i (keys(%{$old})){
            my $value = $old->{$i};
            $new{$i}= $value/$sum;
        }
    }
    return %new;
}

sub phylogeny_summary{
    my $strain_info = shift;
    my $lineage = shift;
    foreach my $i (keys(%{$strain_info})){
        my $status = $strain_info->{$i}->{'status'};
        if ($status !~ /^Y/) {
            next;
        }
        my $taxa = $strain_info->{$i}->{'lineage'};
        my $abun = $strain_info->{$i}->{'f_unique_weighted'};
        my @taxa = split(/\;/, $taxa);
        my $species='';
        my $genus='';
        my $strain = $i;
        if (scalar(@taxa)>=6) {
            if ($taxa[5] =~ /__/) {
                $genus = substr($taxa[5], 3, length($taxa[5])-3);
            }else{
                $genus = $taxa[5];
            }
        }
        if (scalar(@taxa)>=7) {
            if ($taxa[6] =~ /__/) {
                $species = substr($taxa[6], 3, length($taxa[6])-3);
            }else{
                $species = $taxa[6];
            }
        }
        $lineage->{$genus}->{$species}->{$strain}=$abun;
    }
}


sub filter_by_coverage{
    my $strain_info = shift;
    my $min_unique_seqs = shift;
    my $min_coverage_length = shift;
    my $good_coverage_length = shift;
    foreach my $i (keys(%{$strain_info})){
        my $unique_seqs = $strain_info->{$i}->{'n_unique_weighted_found'};
        my $coverage_length = $strain_info->{$i}->{'unique_intersect_bp'};
        my $sourmash_potential_false_negative= $strain_info->{$i}->{'potential_false_negative'};
        if ($coverage_length<$min_coverage_length) {
            $strain_info->{$i}->{'status'} = "coverage length less than $min_coverage_length bp";
            next;
        }
        if ($unique_seqs<$min_unique_seqs) {
            $strain_info->{$i}->{'status'} = "less than $min_unique_seqs unique seqs";
            next;
        }
        if ($coverage_length>=$good_coverage_length) {
            $strain_info->{$i}->{'status'} = "Y: coverage above 50kb";
            next;
        }
        if ($sourmash_potential_false_negative eq 'True') {
            $strain_info->{$i}->{'status'} = "Sourmash potential false positive";
        }
        
        
    }
}

sub read_csv_table{
    my $file = shift;
    my $strain_info = shift;
    my $title = shift;
    open(my $f1, "<$file") or die;
    $title = <$f1>;
    chomp $title;
    $title =~ s/\r//g;
    my @titles = split(/\,/, $title);
    my $name_index = 0;
    for (my $j =0; $j<scalar(@titles);$j++){
        if ($titles[$j] eq 'name') {
            $name_index = $j;
        }
    }
    while (my $line = <$f1>) {
        chomp $line;
        $line  =~ s/\r//g;
        my @coln = parse_csv_line($line);
        my $strain = $coln[$name_index];
        my @names = split(/\ /, $strain);
        $strain = $names[0];
        for (my $i=0; $i<scalar(@titles); $i++){
            $strain_info->{$strain}->{$titles[$i]}=$coln[$i];
            $strain_info->{$strain}->{'status'}='Y';
        }
    }
    return(@titles);
}

sub parse_csv_line {
    my ($line) = @_;
    my @fields;
    my $current_field = '';
    my $in_quotes = 0;
    foreach my $char (split //, $line) {
        if ($char eq '"') {
            $in_quotes = !$in_quotes;
        } elsif ($char eq ',' && !$in_quotes) {
            push @fields, $current_field;
            $current_field = '';
        } else {
            $current_field .= $char;
        }
    }
    push @fields, $current_field;

    return @fields;
}

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
