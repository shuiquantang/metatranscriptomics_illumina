use strict;
use warnings;

package LefseAnalysis;

sub plot_lefse_by_category{
    my $bin_path = $_[0];
    my $category_number = $_[1];
    my $cat_titles = $_[2];
    my $group_id = $_[3];
    my $group_info = $_[4];
    my $abun_table_folder = $_[5];
    my $phylogeny = $_[6];
    my $sample_no = scalar(@{$group_info -> {$group_id} -> {'sample_order'}});
    if ($category_number>0 and $sample_no>5) {
        for (my $j = 0; $j<$category_number; $j++){
            my $cat = $cat_titles->[$j];
            my $lefse_dir = "lefse_$cat";
            my $map_file = "$cat.mapping.file.txt";
            system("mkdir -p $lefse_dir");
	    my @ranks= @{$phylogeny->{'ranks'}};
	    my $rank = $ranks[10-1];
	    my $abun_file = "$abun_table_folder/$rank/$rank.tsv";
            lefse($bin_path, $lefse_dir, $abun_file, $map_file);
        }
    }else{
	print ("skip lefse analysis because group information is not given or there are less than 6 samples in total.\n");
    }
}


sub plot_humann_lefse_by_category{
    my $bin_path = $_[0];
    my $category_number = $_[1];
    my $cat_titles = $_[2];
    my $group_id = $_[3];
    my $group_info = $_[4];
    my $abun_file = $_[5];
    my $tag = $_[6];
    my $sample_no = scalar(@{$group_info -> {$group_id} -> {'sample_order'}});
    if ($category_number>0 and $sample_no>5) {
        for (my $j = 0; $j<$category_number; $j++){
            my $cat = $cat_titles->[$j];
            my $lefse_dir = "lefse_$tag\_$cat";
            my $map_file = "$cat.mapping.file.txt";
            system("mkdir -p $lefse_dir");
	    lefse($bin_path, $lefse_dir, $abun_file, $map_file);
        }
    }else{
	print ("skip lefse analysis because group information is not given or there are less than 6 samples in total.\n");
    }
}

sub lefse{
    my $bin_path = $_[0];
    my $dir = $_[1];
    my $old_abun_file =$_[2];
    my $old_map_file = $_[3];
    my $map_file = "map.txt";
    my $abun_file = "abun.txt";
    my $samples_in_abun_file =create_abun_file($old_abun_file, "$dir/$abun_file");
    create_map_file($old_map_file, "$dir/$map_file", $samples_in_abun_file);
    
    # run lefse
    my $command = "python $bin_path/lefse/transform_data_for_lefse.py -i $dir/$abun_file -m $dir/$map_file  -o $dir/lefse.input.txt >/dev/null 2>&1";
    system($command);
    $command = "python $bin_path/lefse/format_input.py $dir/lefse.input.txt $dir/lefse.in -c 2 -u 1 -o 1000000 >/dev/null 2>&1";
    system($command);
    $command = "python $bin_path/lefse/run_lefse.py $dir/lefse.in $dir/lefse.res.xls >/dev/null 2>&1";
    system($command);
    if (! -e "$dir/lefse.res.xls") {
	print ("lefse failed\n");
	return;
    }
    #system("python $ScriptDir/lefse/plot_res.py lefse/results/lefse.res.xls lefse/results/markers.pdf --format pdf --subclades -1 --max_feature_len 200 --width 40 --height 40 >/dev/null 2>&1");
    $command = "python $bin_path/lefse/plot_cladogram.py $dir/lefse.res.xls $dir/cladogram.pdf --format pdf --label_font_size 3 --class_legend_font_size 6 --labeled_stop_lev 5 --abrv_stop_lev 5 >/dev/null 2>&1";
    system($command);
    $command = "python $bin_path/draw_lefse_plot.py -i $dir/lefse.res.xls -o $dir/biomarkers.pdf >/dev/null 2>&1";
    system($command);
    $command = "python $bin_path/lefse/plot_features.py -f diff --archive zip $dir/lefse.in $dir/lefse.res.xls $dir/biomarkers.zip --format pdf >/dev/null 2>&1";
    system($command);
    
    system("mkdir $dir/biomarkers");
    
    system ("unzip -q $dir/biomarkers.zip -d $dir/biomarkers");
    
    system ("rm $dir/biomarkers.zip");
    rename_lefse_files("$dir");
}

sub rename_lefse_files{
    my $dir = $_[0];
    my $biomarker_folder = "$dir/biomarkers";
    open(my $f1, "<$dir/lefse.res.xls") or die;
    my @lines;
    my $n=1;
    while (my $line = <$f1>) {
	chomp $line;
	my @coln = split("\t", $line);
	my $new_file='';
	my $cat = $coln[2];
	$cat =~ s/\s+//g;
	if (length($cat)>0) {
	    my $taxa = $coln[0];
	    $taxa =~ s/\./\-/g;
	    my $old_file = "1_$taxa.pdf";
	    $new_file = "$n.pdf";
	    if (-e "$biomarker_folder/$old_file") {
		rename("$biomarker_folder/$old_file", "$biomarker_folder/$new_file");
	    }else{
		print("error in lefse analysis: file $biomarker_folder/$old_file is not found.\n")
	    }
	    $n++;
	}
	push (@coln, $new_file);
	push (@lines, \@coln);
    }
    close $f1;
    open(my $f2, ">$dir/lefse.res.xls") or die;
    foreach my $line (@lines){
	my $line = join("\t", @{$line});
	print ($f2 "$line\n");
    }
    close $f2;
}

sub create_map_file{
    my $input = shift;
    my $output = shift;
    my $samples_found =shift;
    open(my $f1, "<$input") or die ("failed to open $input");
    open(my $f2, ">$output") or die ("failed to open $output");
    print($f2 "SampleID\tGroup\n");
    while (my $line = <$f1>) {
        chomp$line;
	my @col = split('\t', $line);
	if ($col[0] =~ /^#/) {
	}else{
	    if (exists $samples_found->{$col[0]}) {
		print($f2 "$col[0]\t$col[3]\n");
	    }
	}
    }
    
    close $f1;
    close $f2;
}

sub create_abun_file{
    my $input = shift;
    my $output = shift;
    my %samples_found;
    open(my $f1, "<$input") or die ("failed to open $input");
    open(my $f2, ">$output") or die ("failed to open $output");
    my $line1 = <$f1>;
    $line1 = <$f1>;
    chomp $line1;
    my @col = split("\t", $line1);
    for (my $i=1;$i<scalar(@col);$i++){
	$samples_found{$col[$i]}++;
    }
    $col[0] = "SampleID";
    $line1 = join("\t", @col);
    print ($f2 "$line1\n");
    while (my $line = <$f1>){
	print ($f2 $line);
    }
    close $f1;
    close $f2;
    return \%samples_found;
}

1;
