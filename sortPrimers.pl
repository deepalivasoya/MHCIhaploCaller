use warnings;
no warnings ('uninitialized', 'substr');
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long;
use Bio::SeqIO;
use Bio::Root::Exception;
use Bio::PrimarySeq;
use Bio::Tools::IUPAC;
use File::Path qw( make_path );

my $sample = shift;
my $work_dir = shift;
my $primer_prefix_info = shift;
my $primers = shift;

# Bovine_for1	forward	For1Rev2
# Bovine_rev2	reverse	For1Rev2
# Bovine_for3	forward	For3Rev1
# Bovine_rev1	reverse	For3Rev1
# DQA_for3	forward	DQA
# DQA_rev2	reverse	DQA
# DQB_for1	forward	DQB
# DQB_rev2	reverse	DQB
# drb3_for	forward	DRB3
# drb3_rev	reverse	DRB3

my %type = ();
my %primers_id = ();

open (INFO, "$primer_prefix_info") or die "Cannot open Primer prefix information file: $primer_prefix_info\n";
while(<INFO>){
	chomp $_;
	@words = split("\t", $_);
	$type{$words[0]} = $words[1];
	$primers_id{$words[0]} = $words[2];
	if (! exists $primer_pairs{$words[2]}){
		$primer_pairs{$words[2]} = 1;
		open ($words[2], ">$work_dir/$sample/$sample.$words[2].fasta") or die "Cannot write $out_identified_primers\n";
	}
}

my $removePrimers = 1;

my %forward_primer = ();
my %reverse_primer_revcomp = ();

my $for = 0;
my $rev = 0;
my $primer_prefix = "";
my $in = Bio::SeqIO -> new( -file => "$primers", -format => 'fasta');
while ( $seq = $in->next_seq() ) {
	my $ambiseq = Bio::PrimarySeq->new(-seq => $seq->seq(), -alphabet => 'dna');

	# Create all possible IUPAC sequences
	my $iupac = Bio::Tools::IUPAC->new(-seq => $ambiseq);
	while ($uniqueseq = $iupac->next_seq()) {
		# process the unique Bio::Seq object.
	 	my $revcomp = $uniqueseq->revcom;
		my $id = $seq->id;
		$count++;
		if ($type{$id} =~ /for/i){
			$for++;
			my $text = "For$for";
			$primer_prefix = $primers_id{$id};
			$$primer_prefix{$text} = $uniqueseq->seq();
			# $forward_primer_revcomp{$for} = $revcomp->seq();
		} elsif ($id =~ /rev/i){
			$rev++;
			my $text = "Rev$rev";
			$primer_prefix = $primers_id{$id};
			$$primer_prefix{$text} = $revcomp->seq();
			# $reverse_primer{$rev} = $uniqueseq->seq();
		} else {
			die "$id in $primers is not valid primer id Read the user manual to see the correct format of file and Try again\n";
		}
	}
}
foreach $prime_pair (keys %primer_pairs){
	foreach $primer (keys %$prime_pair){
		print "$primer\t$$prime_pair{$prime_pair}\n";
	}
}

my $in_fastq = "$work_dir/$sample/$sample.extendedFrags.fastq";
if ( -e $in_fastq ){

	open (OUT_HISTO, ">$work_dir/$sample/$sample.primers.hist");

	my %unique_sequences = ();
	my %primer_group = ();
	my $total_primer_reads = 0;
	my $total_unique_seq = 0;
	my $len = 0;
	my $per = 0;
	my $flag = 0;
	my %stagger1 = ();
	my %stagger2 = ();
	my $stagger = "";

	
	my $in = Bio::SeqIO -> new( -file => "$in_fastq", -format => 'fastq');
	print "Identifiying primers and writing unique variants...\n";
	LOOP5: while ( $seq = $in->next_seq() ) {
		$stagger = "";
		my $sequence = $seq->seq();
		my $id = $seq->id;
		foreach $primer (keys %primer_pairs){
			foreach my $forward (keys %$primer){
				if ($forward =~ /For/ and $sequence =~ /^(AAA|AA|A|)$$primer{$forward}(\w+)/) {
					$forward_stagger = $1;
					foreach my $reverse (keys %$primer){
						if ($reverse =~ /Rev/ and $2 =~ /(\w+)$primer{$reverse}(TTT|TT|T|)$/) {
							if ($removePrimers == 1){
								print $primer ">$id\n$1\n";
							}
							$reverse_stagger = $2;
							$primer_group{"$primer\t$forward_stagger\t$$primer{$forward}\t$$primer{$reverse}\t$reverse_stagger"}++;
							next LOOP5;
						}
					}
				}
			}
		}
		$rev_comp = $seq->revcom();
		$sequence = $rev_comp->seq();
		foreach $primer (keys %primer_pairs){
			foreach my $forward (keys %$primer){
				if ($forward =~ /For/ and $sequence =~ /^(AAA|AA|A|)$$primer{$forward}(\w+)/) {
					$forward_stagger = $1;
					foreach my $reverse (keys %$primer){
						if ($reverse =~ /Rev/ and $2 =~ /(\w+)$primer{$reverse}(TTT|TT|T|)$/) {
							if ($removePrimers == 1){
								print $primer ">$id\n$1\n";
							}
							$reverse_stagger = $2;
							$primer_group{"$primer\t$forward_stagger\t$$primer{$forward}\t$$primer{$reverse}\t$reverse_stagger"}++;
							next LOOP5;
						}
					}
				}
			}
		}
	}

	# $total_unique_seq = keys %unique_sequences;
	# $count = 0;
	# print "Total unique sequences with valid primers: $total_unique_seq\n";
	# foreach $seq (sort {$unique_sequences{$b} <=> $unique_sequences{$a}} keys %unique_sequences){
	# 	$count++;
	# 	$len = length($seq);
	# 	# $frequency_cutoff = ($total_primer_reads * $cutoff) / 100; 
	# 	$per = sprintf("%.2f", ($unique_sequences{$seq} / $total_primer_reads) * 100);
	# 	print OUT_IDENT_PRIMERS ">$count:$len:$unique_sequences{$seq}:$per\n$seq\n";
	# }
	# close (OUT_IDENT_PRIMERS);

	# open (LOG, ">$out_log") or die "Cannot write $out_log\n";
	# print "Writting unique variants log $out_log...\n";
	# foreach $primer (sort {$primer_group{$b} <=> $primer_group{$a}} keys %primer_group){
	# 	print LOG "$primer\t$primer_group{$primer}\n";
	# }
	# close (LOG);
	# print STAT "\t$total_primer_reads\t$total_unique_seq";
	# foreach $stagger (sort {$stagger1{$b} <=> $stagger1{$a}} keys %stagger1){
	# 	print STAT "\tFor:$stagger:$stagger1{$stagger}";
	# 	print "For:$stagger:$stagger1{$stagger}\n";
	# }
	# foreach $stagger (sort {$stagger2{$b} <=> $stagger2{$a}} keys %stagger2){
	# 	print STAT "\tRev:$stagger:$stagger2{$stagger}";
	# 	print "Rev:$stagger:$stagger2{$stagger}\n";
	# }
	# print STAT "\n";

	# foreach $len (sort {$a <=> $b} keys %histo){
	# 	print OUT_HISTO "$len\t$histo{$len}\n";
	# }
	# close (OUT_HISTO);
}
else{
	die "Cannot find the overlapping reads: $work_dir/$sample/$sample.extendedFrags.fastq\n";
}

# open (STAT, ">", "$work_dir/$prefix.$primer_prefix.sortPrimers.stats.txt");
# print STAT "Animal ID\tPrimer\tTotal reads with primer\tTotal unique clustered variants\n";

# if ($isPaired == 1 and $isOverlap == 1){
# 	open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
# 	my $flag_samplesheet = 0;
# 	while(<SAMPLESHEET>){
# 		chomp $_;
# 		@words = split("\t", $_);
# 		if (! defined $words[0]){
# 			die "No sample loaded from samplesheet. Please check the samplesheet again and try again\n";
# 		}
# 		$sample = $words[0];
# 		print "$sample running...\n";

# 		$error = 0;
# 		my $prefix_dir = "$work_dir/$sample/$primer_prefix";
# 		if ( -d $prefix_dir ){
			
# 		} else {
# 			mkdir $prefix_dir or $error = 1;
# 		}
# 		if ($error == 1){
# 			die "Cannot create $prefix_dir directory.\nSomething is wrong while creating $sample directory.Try Again\n";
# 		}

# 		print STAT "$sample\t$primer_prefix";
# 		my $file_to_check = "$work_dir/$sample/$sample.extendedFrags.fastq";
# 		if ( -e $file_to_check ){
# 			if ($isPrimer == 1){
# 				findPrimers($file_to_check, "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.primerSorted.details.txt","$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.unique.variants.fasta", "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.primers.hist", "$work_dir/$sample/$primer_prefix/$sample.$primer_prefix.length.hist");
# 			}
# 			else{

# 			}
# 		} else {
# 			die "Cannot find the overlapping reads: $work_dir/$sample/$sample.extendedFrags.fastq\n";
# 		}
# 	}
# }









