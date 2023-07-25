sub trim_by_abundance{
    my $taxa_link = shift;
    my $abs_cutoff = shift;
    my $rel_cutoff = shift;
    my %taxa_abun = %{$taxa_link};
    my @taxa = keys%taxa_abun;
    my @taxa_to_delete;
    foreach my $i (@taxa){
	my $abun = $taxa_abun{$i};
	if ($abun<$abs_cutoff) {
	    delete($taxa_abun{$i});
	    push(@taxa_to_delete, $i);
	}elsif($i =~ /unknown/){# skip the unknown taxa that will be treated separately in the function of resolve_unknown()
	    delete($taxa_abun{$i});
	}
	
    }
    @taxa = sort{$taxa_abun{$b}<=>$taxa_abun{$a}} keys%taxa_abun;
    my $abun_of_last_record=0.01;
    
    my $delete_controller = 0;
    foreach my $i (@taxa){
	if ($i =~ 'Listeria') {
	    my $m=0;
	}
	
	my $abun =$taxa_abun{$i};
	if ($delete_controller == 0) {
	    if ($abun/$abun_of_last_record<$rel_cutoff) {
	        push(@taxa_to_delete, $i);
	        $delete_controller = 1;
	    }else{
		$abun_of_last_record=$abun;
	    }
	}else{
	    push(@taxa_to_delete, $i);
	}
	
    }
    return @taxa_to_delete;
}
