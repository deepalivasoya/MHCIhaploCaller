use warnings;
use Cwd;
use Bio::SeqIO;
use Bio::Root::Exception;

my $samplesheet = shift;
my $work_dir = shift;
my $prefix = shift;
my $fastqdir = shift;
my $minLength = shift;
my $minQuality = shift;

my $isOverlap = 1;
my $isPaired = 1;

my %read1 = ();
my %read2 = ();
my %samples = ();
my $sickle = "sickle";
my $error = 0;

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	$sample = $words[0];
	$error = 0;
	my $sample_dir = "$work_dir/$sample";
	if (! -d $sample_dir ){
		mkdir $sample_dir or $error = 1;
	}
	if ($error == 1){
		die "Cannot create $sample_dir directory.\nSomething is wrong while creating $sample directory.Try Again\n";
	}
	if (defined $words[0]){
		if ($isPaired == 1){
			if (defined $words[1] and defined $words[2] and $words[1] ne "" and $words[2] ne ""){
				$read1{$words[0]} = $words[1];
				$read2{$words[0]} = $words[2];
				$r1 = "$fastqdir/$read1{$sample}";
				$r2 = "$fastqdir/$read2{$sample}";
				if (! -e $r1 or ! -e $r2){
					die "Something is wrong in either given samplesheet or seuqencing read files.\nCannot read $r2 or $r2 for sample $samples{$sample}. Try again.\n";
	    		}
			} else {
				die "No samples loaded from samplesheet. Please check the samplesheet again and try again\n";
			}
		}
	}
	else{
		die "No samples loaded from samplesheet. Please check the samplesheet again and try again\n";
	}
}

my $paired_kept = 0;
my $paired_dis = 0;
my $singles_dis = 0;
my $single_kept = 0;
open (READS_STAT, ">", "$work_dir/$prefix.trimming.stats.txt");
print READS_STAT "Sample\tAnimal ID\tTotal sequencing read pairs\tTotal trimmed pairs\tTotal trimmed singles\n";
if ($isPaired == 1 and $isOverlap == 1){
	open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
	while(<SAMPLESHEET>){
		@words = split("\t", $_);
		$sample = $words[0];
		print READS_STAT "$sample\t";
		print "Running sickle for $sample...\n";
		system ("$sickle pe -f $fastqdir/$read1{$sample} -r $fastqdir/$read2{$sample} -o $work_dir/$sample/$sample.trimmed.R1.fastq.gz -p $work_dir/$sample/$sample.trimmed.R2.fastq.gz -s $work_dir/$sample/$sample.trimmed.singles.fastq.gz -t sanger -q $minQuality -l $minLength -g > $work_dir/$sample/$sample.trimming.log");
		$file_to_check1 = "$work_dir/$sample/$sample.trimmed.R1.fastq.gz";
		$file_to_check2 = "$work_dir/$sample/$sample.trimmed.R2.fastq.gz";
		if ( -e $file_to_check1 and -e $file_to_check2 ){
			$file_to_check = "$work_dir/$sample/$sample.trimming.log";
			if ( -e $file_to_check ){ 
				open (IN, "$file_to_check") or print "Cannot open $file_to_check\n";
				$paired_kept = 0;
				$paired_dis = 0;
				$singles_dis = 0;
				$single_kept = 0;
				$total = 0;
				while(<IN>){
					chomp $_;
					if (/FastQ paired records kept: (\d+)/){
						$paired_kept = $1/2;
					}
					elsif (/FastQ single records kept: (\d+)/){
						$single_kept = $1;
					}
					elsif (/FastQ paired records discarded: (\d+)/){
						$paired_dis = $1/2;
					}
					elsif (/FastQ single records discarded: (\d+)/){
						$singles_dis = $1;
					}
					$total = $paired_kept + $single_kept + $paired_dis;
				}
				close (IN);
				print READS_STAT "\t$total\t$paired_kept\t$single_kept";		
			}
		}
		print READS_STAT "\n";
	}
}
