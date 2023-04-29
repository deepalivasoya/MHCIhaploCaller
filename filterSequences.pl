use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Bio::SeqIO;
use Bio::Root::Exception;
use Scalar::Util qw(looks_like_number);
use File::Path qw( make_path );

my $samplesheet = shift or die "Cannot find samplesheet";
my $work_dir = shift or die "Cannot find working directory";
my $prefix = shift or die "Cannot find prefix";
my $primer_prefix = shift or die "Cannot find primer prefix";
my $cutoff = shift or die "Cannot find cutoff";
my $ampliconeSize = shift or die "Cannot find amplicone size";
my $primers = shift;

my $in = Bio::SeqIO -> new( -file => "$primers", -format => 'fasta');
while ( $seq = $in->next_seq() ) {
	my $id = $seq->id;
	if ($id =~ /for/i){
		$for = length($seq->seq());
	} 
	elsif ($id =~ /rev/i){
		$rev = length($seq->seq());
	}
}
print "For: $for\nRev: $rev\n";

sub possible_chimaera {
	my $sequence1 = $_[0];
	my $sequence2 = $_[1];
	my $candidate = $_[2];
	$candidate = uc($candidate);
	$sequence1 = uc($sequence1);
	$sequence2 = uc($sequence2);
	my %base1 = ();
	my %base2 = ();
	for (my $i = 0; $i < length($sequence1); $i++) {
		$base1{$i} = substr($sequence1, $i, 1 );
	}
	for (my $i = 0; $i < length($sequence2); $i++) {
		$base2{$i} = substr($sequence2, $i, 1 );
	}
	$len_diff1 = length($candidate) - length($sequence1);
	$len_diff2 = length($candidate) - length($sequence2);

	$stop_base = 0;
	$stop_base1 = 0;
	LOOP1: for (my $i = 0; $i < length($candidate); $i++) {
		my $base = substr($candidate, $i, 1);
		if ($base1{$i} ne $base) {
			$stop_base = $i;
			last LOOP1;
		}
	}
	if ($stop_base > 0){
		LOOP2: for (my $i = (length($candidate) - 1); $i >= 0; $i--) {
			my $base = substr($candidate, $i, 1);
			if ($base2{$i - $len_diff2} ne $base) {
				$stop_base1 = $i;
				last LOOP2;
			}
		}
	}
	else {
		LOOP2: for (my $i = 0; $i < length($candidate); $i++) {
			my $base = substr($candidate, $i, 1);
			if ($base2{$i} ne $base) {
				$stop_base = $i;
				last LOOP2;
			}
		}
		if ($stop_base > 0){
			LOOP3: for (my $i = (length($candidate) - 1); $i >= 0; $i--) {
				my $base = substr($candidate, $i, 1);
				if ($base1{$i - $len_diff1} ne $base) {
					$stop_base1 = $i;
					last LOOP3;
				}
			}
		}
	}
	if ($stop_base1 > 0 and $stop_base > 0 and $stop_base1 < $stop_base){
		return 1;
	}
	else{
		return 0;
	}
		
}


sub filterSeq{
	my $in_fasta = $_[0];
	my $out_log = $_[1];
	my $out_selected_fasta = $_[2];
	my $out_1snp = $_[3];
	my $sa = $_[4];

	my %selected_variants = ();
	my %selected_variants_ids = ();
	my $countSingletons = 0;
	my $countLowCoverage = 0;
	my $countNs = 0;
	my $countChimaera = 0;
	my $countLenDiff = 0;
	my $countPCRerrors1 = 0;
	my $countPCRerrors2 = 0;
	my $count1snp = 0;
	my $readsLowCoverage = 0;
	my $readsNs = 0;
	my $readsChimaera = 0;
	my $readsLenDiff = 0;
	my $readsPCRerrors1 = 0;
	my $readsPCRerrors2 = 0;
	my $reads1snp = 0;
	my $total = 0;
	my $read_freqency = 0;
	my $final_selected = 0;
	my $len = 0;
	my $total_len = 0;
	my $total_selected_variants = 0;
	my $ave_len = 0;
	my %lengths = ();
	my $flag_gap = 0;
	my $total_amplicone_size = 0;

	open (LOG, ">$out_log") or die "Cannot write $out_log\n";
	open (OUT_SELECTED, ">$out_selected_fasta") or die "Cannot write $out_selected_fasta\n";
	open (OUT_1SNP, ">$out_1snp") or die "Cannot write $out_1snp\n";

	my $seqio = Bio::SeqIO->new(-file => "$in_fasta", -format => "fasta");
	LOOP8: while(my $seq = $seqio -> next_seq) {
		$id = $seq->id;
		@words = split(":", $id);
		if ($words[3] >= $cutoff){
			# print "$id\n";
			$total = $total + $words[2];
			$selected_variants{$seq->seq()} = $words[2];
			$selected_variants_ids{$seq->seq()} = $id;
			$len = length($seq->seq());
			$lengths{$len}++;
		}
		else{
			last LOOP8;
		}
	}

	$total_selected_variants = keys %selected_variants;
	# print "$total_selected_variants\n";
	if ($ampliconeSize eq 0){
		LOOP9: foreach my $len (sort {$lengths{$b} <=> $lengths{$a}} keys %lengths){
			$ave_len = $len;
			print LOG "Average amplicone length = $len\n";
			last LOOP9;
		}
	}
	else{
		$ave_len = $ampliconeSize;
	}
	$total_amplicone_size = $for + $ampliconeSize + $rev;

	foreach my $seq (sort {$selected_variants{$b} <=> $selected_variants{$a}} keys %selected_variants){
		# print "Variant: $selected_variants_ids{$seq}...\n";
		$flag_gap = 0;
		my @nucl = split //, $seq;
		foreach my $n (@nucl){
			if ($n !~ /[A|T|C|G]/i){
				$flag_gap = 1;
			}
		}

		if ($flag_gap == 1){
			print LOG "$selected_variants_ids{$seq}\tN\t-\t$seq\n";
			$countNs++;
			$readsNs = $readsNs + $selected_variants{$seq};
		}
		elsif ($selected_variants{$seq} == 1){
			print LOG "$selected_variants_ids{$seq}\tsingleton\t-\t$seq\n";
			$countSingletons++;
		}
		elsif ($selected_variants{$seq} <= 10){
			print LOG "$selected_variants_ids{$seq}\tLow coverage\t-\t$seq\n";
			$countLowCoverage++;
			$readsLowCoverage = $readsLowCoverage + $selected_variants{$seq};
		}
		elsif (length($seq) >= ($ave_len + 9) or length($seq) <= ($ave_len - 9)){
			my $lendiff = $ampliconeSize - length($seq);
			print LOG "$selected_variants_ids{$seq}\tLength Difference\t$ampliconeSize:$lendiff\t$seq\n";
			$countLenDiff++;
			$readsLenDiff = $readsLenDiff + $selected_variants{$seq};
		}
		else{
			my $flag_chimaera = 0;
			LOOP6: foreach $seq1 (sort {$selected_variants{$b} <=> $selected_variants{$a}} keys %selected_variants){
				# print "SEQ1: $selected_variants_ids{$seq1}\n";
				if ($seq ne $seq1){
					if ($selected_variants{$seq1} > $selected_variants{$seq}){
						LOOP8: foreach $seq2 (sort {$selected_variants{$b} <=> $selected_variants{$a}} keys %selected_variants){
						 	# print "SEQ2: $selected_variants_ids{$seq2}\n";
							if ($seq ne $seq2){
								if ($seq1 ne $seq2){
									if ($selected_variants{$seq2} > $selected_variants{$seq}){
										# print "Checking chimaera $selected_variants_ids{$seq} of $selected_variants_ids{$seq1} and $selected_variants_ids{$seq2}...\n";
										$check = possible_chimaera($seq1, $seq2, $seq);
										if ($check == 1){
											$flag_chimaera = 1;
											print LOG "$selected_variants_ids{$seq}\tChimaera\t$selected_variants_ids{$seq1},$selected_variants_ids{$seq2}\t$seq\n";
											last LOOP6;
										}
									}
								}
							}
							else{
								last LOOP8;
							}
						}
					}
					else{
						last LOOP6;
					}
				}
				else{
					last LOOP6;
				}
			}
			if ($flag_chimaera == 1){
				$countChimaera++;
				$readsChimaera = $readsChimaera + $selected_variants{$seq};
			}
			else{
				my $flag_pcr_error = 0;
				LOOP7: foreach $seq1 (sort {$selected_variants{$b} <=> $selected_variants{$a}} keys %selected_variants){
					if (length($seq1) eq length($seq) and $selected_variants{$seq1} > $selected_variants{$seq}){
						my $i = 0;
						my $match = 0;
						my $mismatch = "";
						while ($i < length($seq)){
							$sub1 = substr ($seq, $i, 1);
							$sub2 = substr ($seq1, $i, 1);
							if ( $sub1 eq $sub2 ){
								$match++;
							}
							else{
								if ($mismatch eq ""){
									$j = $i + 1;
									$mismatch = "$j/$sub2->$sub1";
								}else{
									$j = $i + 1;
									$mismatch = $mismatch.",$j/$sub2->$sub1";
								}
							}
							$i++;
						}
						my $fc = sprintf('%.2f', $selected_variants{$seq1} / $selected_variants{$seq});
						if (($selected_variants{$seq1} / $selected_variants{$seq}) >= 30 and $match eq length($seq) - 1){
							$flag_pcr_error = 1;
							open (SORT, "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.primers.hist") or die "Cannot open $work_dir/$sample/$primer_prefix/$sample.$primer_prefix.primers.hist\n";
							LOOP10: while(<SORT>){
								chomp $_;
								@info1 = split(/\t/, $_);
								@info2 = split(/\s/, $info1[0]);
								$for_stagger = $info2[1];
								$for_primer = $info2[2];
								$rev_primer = $info2[4];
								$rev_stagger = $info2[5];
								last LOOP10;
							}
							$ampliconeseq = $for_stagger.$for_primer.$seq.$rev_primer.$rev_stagger;
							$ampliconeseq1 = $for_stagger.$for_primer.$seq1.$rev_primer.$rev_stagger;
							$total_amplicone_size = length($for_stagger) + $for + $ampliconeSize + $rev + length($rev_stagger);
							$non_overlapping_len = $total_amplicone_size - 300;
							@mismatch_info = split("/", $mismatch);
							if ($mismatch_info[0] > 300){
								$read1_pos = "";
								$read2_pos = $total_amplicone_size - ($mismatch_info[0] + $for + length($for_stagger) + 1);
								$subseq = substr $ampliconeseq, (length($for_stagger) + $for + $mismatch_info[0])-5, 9;
								$subseq1 = substr $ampliconeseq1, (length($for_stagger) + $for + $mismatch_info[0])-5, 9;
							}
							elsif ($mismatch_info[0] <= $non_overlapping_len){
								$read1_pos = $mismatch_info[0] + length($for_stagger) + $for;
								$read2_pos = "";
								$subseq = substr $ampliconeseq, $read1_pos-5, 9;
								$subseq1 = substr $ampliconeseq1, $read1_pos-5, 9;
							}
							else{
								$read1_pos = $mismatch_info[0] + length($for_stagger) + $for;
								$read2_pos = $total_amplicone_size - ($mismatch_info[0] + $for + length($for_stagger) + 1);
								$subseq = substr $ampliconeseq, $read1_pos-5, 9;
								$subseq1 = substr $ampliconeseq1, $read1_pos-5, 9;
							}
							print LOG "$selected_variants_ids{$seq}\t1 PCR error\t$selected_variants_ids{$seq1}/$fc/$mismatch/$read1_pos/$read2_pos/$subseq/$subseq1/$for_stagger/$rev_stagger\t$seq\n";
							$countPCRerrors1++;
							$readsPCRerrors1 = $readsPCRerrors1 + $selected_variants{$seq};
							last LOOP7;
						}
						elsif (($selected_variants{$seq1} / $selected_variants{$seq}) >= 100 and $match eq length($seq) - 2){
							$flag_pcr_error = 1;

							print LOG "$selected_variants_ids{$seq}\t2 PCR error\t$selected_variants_ids{$seq1}/$fc/$mismatch\t$seq\n";
							$countPCRerrors2++;
							$readsPCRerrors2 = $readsPCRerrors2 + $selected_variants{$seq};
							last LOOP7;
						}
						elsif ($match eq length($seq) - 1){
							$flag_pcr_error = 1;
							open (SORT, "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.primers.hist") or die "Cannot open $work_dir/$sample/$primer_prefix/$sample.$primer_prefix.primers.hist\n";
							LOOP10: while(<SORT>){
								chomp $_;
								@info1 = split(/\t/, $_);
								@info2 = split(/\s/, $info1[0]);
								$for_stagger = $info2[1];
								$for_primer = $info2[2];
								$rev_primer = $info2[4];
								# print "$for_primer\t$rev_primer\n";
								$rev_stagger = $info2[5];
								print "$for_stagger\t$for_primer\t$rev_primer\t$rev_stagger\n";
								last LOOP10;
							}
							$ampliconeseq = $for_stagger.$for_primer.$seq.$rev_primer.$rev_stagger;
							$ampliconeseq1 = $for_stagger.$for_primer.$seq1.$rev_primer.$rev_stagger;
							# print length($ampliconeseq)."\t".length($ampliconeseq1)."\n";
							$total_amplicone_size = length($for_stagger) + $for + $ampliconeSize + $rev + length($rev_stagger);
							print length($ampliconeseq)."\t".length($ampliconeseq1)."\t".$total_amplicone_size."\n";
							$non_overlapping_len = $total_amplicone_size - 300;
							@mismatch_info = split("/", $mismatch);
							if ($mismatch_info[0] > 300){
								$read1_pos = "";
								$read2_pos = $total_amplicone_size - ($mismatch_info[0] + $for + length($for_stagger) + 1);
								$subseq = substr $ampliconeseq, (length($for_stagger) + $for + $mismatch_info[0])-5, 9;
								$subseq1 = substr $ampliconeseq1, (length($for_stagger) + $for + $mismatch_info[0])-5, 9;
							}
							elsif ($mismatch_info[0] <= $non_overlapping_len){
								$read1_pos = $mismatch_info[0] + length($for_stagger) + $for;
								$read2_pos = "";
								$subseq = substr $ampliconeseq, $read1_pos-5, 9;
								$subseq1 = substr $ampliconeseq1, $read1_pos-5, 9;
							}
							else{
								$read1_pos = $mismatch_info[0] + length($for_stagger) + $for;
								$read2_pos = $total_amplicone_size - ($mismatch_info[0] + $rev + length($rev_stagger) + 1);
								$subseq = substr $ampliconeseq, $read1_pos-5, 9;
								$subseq1 = substr $ampliconeseq1, $read1_pos-5, 9;
							}
							
							print LOG "$selected_variants_ids{$seq}\t1SNP\t$selected_variants_ids{$seq1}/$fc/$mismatch/$read1_pos/$read2_pos/$subseq/$subseq1/$for_stagger/$rev_stagger\t$seq\n";
							$count1snp++;
							$reads1snp = $reads1snp + $selected_variants{$seq};
							print OUT_SELECTED ">$selected_variants_ids{$seq}\n$seq\n";
							print OUT_1SNP "$sample\t$selected_variants_ids{$seq}\t$selected_variants_ids{$seq1}\t$mismatch/$read1_pos/$read2_pos/$subseq/$subseq1/$for_stagger/$rev_stagger\t$seq\t$seq1\n";
							last LOOP7;
						}
					}
				}
				if ($flag_pcr_error == 0){
					print LOG "$selected_variants_ids{$seq}\tUnknown\t-\t$seq\n";
					print OUT_SELECTED ">$selected_variants_ids{$seq}\n$seq\n";
					$final_selected++;
					$read_freqency = $read_freqency + $selected_variants{$seq};
				}
			}
		}
	}
	close (OUT_SELECTED);

	# print "Writting unique variants log $out_log...\n";
	print LOG "\n\n****** SUMMARY****** \n";
	print LOG "Animal ID\t$sample\n";
	print LOG "Primer\t$primer_prefix\n";
	print LOG "Total variants above cutoff\t$total_selected_variants\t$total\n";
	print LOG "Singletons\t$countSingletons\t$countSingletons\n";
	print LOG "Low coverage\t$countLowCoverage\t$readsLowCoverage\n";
	print LOG "N in sequences\t$countNs\t$readsNs\n";
	print LOG "Chimaeras\t$countChimaera\t$readsChimaera\n";
	print LOG "1 PCR errors (30 fold change)\t$countPCRerrors1\t$readsPCRerrors1\n";
	print LOG "2 PCR errors (100 fold change)\t$countPCRerrors2\t$readsPCRerrors2\n";
	print LOG "1 snp\t$count1snp\t$reads1snp\n";
	print LOG "Length difference\t$countLenDiff\t$readsLenDiff\n";
	print LOG "Selected\t$final_selected\t$read_freqency\n";
	close (LOG);
	print STAT "$total_selected_variants\t$countSingletons\t$countLowCoverage\t$countNs\t$countChimaera\t$countLenDiff\t$countPCRerrors1\t$countPCRerrors2\t$count1snp\t$final_selected\t\t$total\t$countSingletons\t$readsLowCoverage\t$readsNs\t$readsChimaera\t$readsLenDiff\t$readsPCRerrors1\t$readsPCRerrors2\t$reads1snp\t$read_freqency\n";
}

open (STAT, ">", "$work_dir/$prefix.$primer_prefix.filterSequences.stats.txt");
print STAT "Animal ID\tPrimer\tTotal variants above cutoff\tSingletons\tLow coverage\tN in sequence\tChimaeras\tLength difference\t1 PCR errors (30 fold change)\t2 PCR errors (100 fold change)\t1SNP\tSelected\t\tTotal: #reads\tSingletons\tLow coverage: #reads\tN in sequence: #reads\tChimaera: #reads\tLength difference: #reads\t1 PCR error: #reads\t2 PCR error: #reads\t1SNP: #reads\tSelected: #reads\n";
open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	if (! defined $words[0]){
		die "No sample loaded from samplesheet. Please check the samplesheet again and try again\n";
	}
	$sample = $words[0];
	print "$sample filtering $primer_prefix...\n";

	print STAT "$sample\t$primer_prefix\t";
	my $file_to_check = "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.unique.variants.fasta";
	if ( -e $file_to_check ){
		filterSeq($file_to_check, "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.filterSequences.log", "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.selected.variants.fasta", "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.selected.variants.1snp.info.txt", $sample);

	} else {
		die "Cannot find the overlapping reads: $file_to_check\n";
	}
	# close(STAT);
}










