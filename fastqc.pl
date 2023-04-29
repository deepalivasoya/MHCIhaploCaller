use warnings;
use Cwd;
use Bio::SeqIO;
use Bio::Root::Exception;

my $samplesheet = shift;
my $work_dir = shift;
my $prefix = shift;
my $fastqdir = shift;

my $isPaired = 1;

my %read1 = ();
my %read2 = ();	
my %samples = ();
my $fastqc = "fastqc";
my $error = 0;

open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
while(<SAMPLESHEET>){
	chomp $_;
	@words = split("\t", $_);
	$sample = $words[0];
	$error = 0;
	my $sample_dir = "$work_dir/$sample";
	my $fastqc_dir = "$work_dir/$sample/FastQC";
	if (! -d $sample_dir ){
		mkdir $sample_dir or $error = 1;
		mkdir $fastqc_dir or $error = 1;
	}
	if (! -d $fastqc_dir ){
		mkdir $fastqc_dir or $error = 1;
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

if ($isPaired == 1){
	open (SAMPLESHEET, "$samplesheet") or die "Cannot open uploaded $samplesheet. Try again.\n";
	while(<SAMPLESHEET>){
		@words = split("\t", $_);
		$sample = $words[0];
		print "Running Fastqc for $sample R1...\n";
		system ("$fastqc -f fastq -o $work_dir/$sample/FastQC $fastqdir/$read1{$sample}");
		system ("$fastqc -f fastq -o $work_dir/$sample/FastQC $work_dir/$sample/$sample.trimmed.R1.fastq.gz");
		print "Running Fastqc for $sample R2...\n";
		system ("$fastqc -f fastq -o $work_dir/$sample/FastQC $fastqdir/$read2{$sample}");
		system ("$fastqc -f fastq -o $work_dir/$sample/FastQC $work_dir/$sample/$sample.trimmed.R2.fastq.gz");
		print "Running Fastqc for $sample overlapped reads...\n";
		system ("$fastqc -f fastq -o $work_dir/$sample/FastQC $work_dir/$sample/$sample.extendedFrags.fastq");
	}
}
