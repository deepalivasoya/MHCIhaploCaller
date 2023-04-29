use warnings;
use Cwd;
use Bio::SeqIO;
use Bio::Root::Exception;

my $samplesheet = shift;
my $work_dir = shift;
my $prefix = shift;
my $min_overlap = shift;

my $isOverlap = 1;
my $isPaired = 1;

my $sample;
my $flash = "flash";
my $read1;
my $read2;

# open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
# while(<SAMPLESHEET>){
# 	chomp $_;
# 	@words = split("\t", $_);
# 	$sample = $words[0];
# 	$error = 0;
# 	my $sample_dir = "$work_dir/$sample";
# 	if (! -d $sample_dir ){
# 		mkdir $sample_dir or $error = 1;
# 	}
# 	if ($error == 1){
# 		die "Cannot create $sample_dir directory.\nSomething is wrong while creating $sample directory.Try Again\n";
# 	}
# 	if (defined $words[0]){
# 		if ($isPaired == 1){
# 			if (defined $words[1] and defined $words[2] and $words[1] ne "" and $words[2] ne ""){
# 				$read1{$words[0]} = $words[1];
# 				$read2{$words[0]} = $words[2];
# 				$r1 = "$fastqdir/$read1{$sample}";
# 				$r2 = "$fastqdir/$read2{$sample}";
# 				if (! -e $r1 or ! -e $r2){
# 					die "Something is wrong in either given samplesheet or seuqencing read files.\nCannot read $r2 or $r2 for sample $samples{$sample}. Try again.\n";
# 	    		}
# 			} else {
# 				die "No samples loaded from samplesheet. Please check the samplesheet again and try again\n";
# 			}
# 		}
# 	}
# 	else{
# 		die "No samples loaded from samplesheet. Please check the samplesheet again and try again\n";
# 	}
# }

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
open (READS_STAT, ">", "$work_dir/$prefix.ovelapping.stats.txt");
print READS_STAT "Sample\tAnimal_ID\tIN\tOUT\n";
if ($isPaired == 1 and $isOverlap == 1){
	while(<SAMPLESHEET>){
		chomp $_;
		@words = split("\t", $_);
		$sample = $words[0];
		print "Running FLASH for $sample...\n";
		$read1 = "$work_dir/$sample/$sample.trimmed.R1.fastq.gz";
		$read2 = "$work_dir/$sample/$sample.trimmed.R2.fastq.gz";
		# $read1 = $fastqdir/$read1{$sample};
		# $read2 = $fastqdir/$read2{$sample};
		print "$read1\t$read2\n";
		if (-e $read1 and -e $read2){
			system ("$flash -m $min_overlap -M 300 -O -o $sample -d $work_dir/$sample $read1 $read2 > $work_dir/$sample/$sample.flash.log");
			my $file_to_check = "$work_dir/$sample/$sample.extendedFrags.fastq";
			if ( -e $file_to_check ){
				$file_to_check1 = "$work_dir/$sample/$sample.flash.log";
				if ( -e $file_to_check1 ){
					open (IN, "$file_to_check1") or die "Cannot open $file_to_check1\n";
					print READS_STAT "$sample\t";
					while(<IN>){
						chomp $_;
						if (/Total pairs:\s+(\d+)/){
							print READS_STAT "\t$1";
						}
						elsif (/Combined pairs:\s+(\d+)/){
							print READS_STAT "\t$1";
						}
					}
					close (IN);
				}
				else{
					print "Cannot find $file_to_check1\n";
				}
			}
			else{
				print "Cannot find $file_to_check\n";
			}
		}
		else{
			print "Cannot find $read1 or/and $read2\n";
		}
		print READS_STAT "\n";
	}
}
