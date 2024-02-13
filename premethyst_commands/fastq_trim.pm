package premethyst_commands::fastq_trim;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_trim");

sub fastq_trim {

# defaults
$min_RL = 30;
$threads = 1;
$r2_trim = 0;
$a1 = "CTATCTCTTATA";
$a2 = "AGATCGGAAGAGC";

getopts("O:1:2:m:t:e:ur:a:b:", \%opt);

$die = "

premethyst fastq-trim (options) -O OutputPrefix -1 read1.fq.gz -2 read2.fq.gz
       or  trim

Runs trim galore with settings specified for sciMETv2/3.
Uses sciMETv2/3 adapter sequences (hard coded) and many
of the trim galore defaults in paired mode.

Trimming will complete and produce a 'trim.complete'
file, at which point subsequent processing can happen.
The command will continue to run, producing additional
QC stats that are not required prior to alignment.

Options:
-O   [STR]   Output Prefix (req)
-1   [STR]   Read1 fastq (req)
-2   [STR]   Read2 fastq (req)
-m   [INT]   Min read length (def = $min_RL)
-t   [INT]   Threads to use (def = $threads)
-e   [INT]   Trim bases from the end of read 2 after adapter trim (def = $r2_trim)
               (set to 0 for SL preps since Hmer is not incorporated, 10 for LA)
-u           Retain unpaired reads (def = discard)
-r   [STR]   Report to specified slack channel when trimming is complete.

-a   [STR]   Adapter 1 sequence (def = $a1)
-b   [STR]   Adapter 2 sequence (def = $a2)

Executable Commands (from $DEFAULTS_FILE)
   trim_galore: $trim_galore
   slack:       $slack
   zcat:        $zcat

";

# universal options

if (!defined $opt{'O'}) {die "\nERROR: Specify output prefix as -O\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'e'}) {$r2_trim = $opt{'e'}};
if (defined $opt{'a'}) {$a1 = $opt{'a'}};
if (defined $opt{'b'}) {$a2 = $opt{'b'}};

if (!defined $opt{'1'} || !defined $opt{'2'}) {die "\nERROR: Reads 1 and 2 must be specified.\n$die"};

if (defined $opt{'u'}) {
	if ($r2_trim > 0) {
		$trim_command = "$trim_galore -a $a1 -a2 $a2 --three_prime_clip_R2 $r2_trim -j $threads --paired --retain_unpaired $opt{'1'} $opt{'2'} >> $opt{'O'}.trim.log 2>> $opt{'O'}.trim.log";
	} else {
		$trim_command = "$trim_galore -a $a1 -a2 $a2 -j $threads --paired --retain_unpaired $opt{'1'} $opt{'2'} >> $opt{'O'}.trim.log 2>> $opt{'O'}.trim.log";
	}
} else {
	if ($r2_trim > 0) {
		$trim_command = "$trim_galore -a $a1 -a2 $a2 --three_prime_clip_R2 $r2_trim -j $threads --paired $opt{'1'} $opt{'2'} >> $opt{'O'}.trim.log 2>> $opt{'O'}.trim.log";
	} else {
		$trim_command = "$trim_galore -a $a1 -a2 $a2 -j $threads --paired $opt{'1'} $opt{'2'} >> $opt{'O'}.trim.log 2>> $opt{'O'}.trim.log";
	}
}

system("echo 'RUNNING: $trim_command' >> $opt{'O'}.trim.log");
system($trim_command);

# rename to output prefix names
$out1 = $opt{'1'}; $out1 =~ s/\.fq\.gz$//; $out1 .= "_val_1.fq.gz";
$out2 = $opt{'2'}; $out2 =~ s/\.fq\.gz$//; $out2 .= "_val_2.fq.gz";
if  (defined $opt{'u'}) {
	$out1u = $opt{'1'}; $out1u =~ s/\.fq\.gz$//; $out1u .= "_unpaired_1.fq.gz";
	$out2u = $opt{'2'}; $out2u =~ s/\.fq\.gz$//; $out2u .= "_unpaired_2.fq.gz";
}

system("mv $out1 $opt{'O'}.trimmed.paired.R1.fq.gz");
system("mv $out2 $opt{'O'}.trimmed.paired.R2.fq.gz");
if (defined $opt{'u'}) {
	system("mv $out1u $opt{'O'}.trimmed.unpaired.R1.fq.gz");
	system("mv $out2u $opt{'O'}.trimmed.unpaired.R2.fq.gz");
}

# write trimming complete file to note further processing can happen
system("date > $opt{'O'}.trim.complete");

if (defined $opt{'r'}) {
	$message = "Trimming complete for $opt{'O'}! Generating per-cell stats, but alignment can proceed.";
	system("$slack -m \"$message\" $opt{'r'} >/dev/null 2>/dev/null");
}

# use single core to get stats on the trim in barcode-aware manner.
%BARC_IN_ct = (); %BARC_OUT_ct = (); %BARC_R1_up = (); %BARC_R2_up = ();
open IN, "$zcat $opt{'1'} |";
while ($tag = <IN>) {
	chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
	$BARC_IN_ct{$tag}++;
	$BARC_OUT_ct{$tag} = 0;
	$null = <IN>; $null = <IN>; $null = <IN>;
} close IN;

open IN, "$zcat $opt{'O'}.trimmed.paired.R1.fq.gz |";
while ($tag = <IN>) {
	chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
	$BARC_OUT_ct{$tag}++;
	$null = <IN>; $null = <IN>; $null = <IN>;
} close IN;

if (defined $opt{'u'}) {
	open IN, "$zcat $opt{'O'}.trimmed.unpaired.R1.fq.gz |";
	while ($tag = <IN>) {
		chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
		$BARC_R1_up{$tag}++;
		$null = <IN>; $null = <IN>; $null = <IN>;
	} close IN;

	open IN, "$zcat $opt{'O'}.trimmed.unpaired.R2.fq.gz |";
	while ($tag = <IN>) {
		chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
		$BARC_R2_up{$tag}++;
		$null = <IN>; $null = <IN>; $null = <IN>;
	} close IN;
}

open RPT, ">$opt{'O'}.trimmed.stats.txt";
foreach $barc (keys %BARC_IN_ct) {
	$pct = sprintf("%.2f", ($BARC_OUT_ct{$barc}/$BARC_IN_ct{$barc})*100);
	if (!defined $BARC_R1_up{$barc}) {$BARC_R1_up{$barc} = 0};
	if (!defined $BARC_R2_up{$barc}) {$BARC_R2_up{$barc} = 0};
	$pct1u = sprintf("%.2f", ($BARC_R1_up{$barc}/$BARC_IN_ct{$barc})*100);
	$pct2u = sprintf("%.2f", ($BARC_R2_up{$barc}/$BARC_IN_ct{$barc})*100);
	
	print RPT "$barc\t$BARC_IN_ct{$barc}\t$BARC_OUT_ct{$barc}\t$pct\t$BARC_R1_up{$barc}\t$pct1u\t$BARC_R2_up{$barc}\t$pct2u\n";
} close RPT;

exit;

}

1;