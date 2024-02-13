package premethyst_commands::fastq_align;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fastq_align");

sub fastq_align {

getopts("O:1:2:t:o:R:r:X", \%opt);

$threads = 1;
$o_threads = 1;

$die = "

premethyst fastq-align (options) -R [reference path] -O [output prefix] -1 [read1.trimmed.fq.gz] -2 [read2.trimmed.fq.gz]
      or   align

Wrapper for BSBOLT to run alignment of sciMETv2 reads.
reads can be a list that is comma-separated.

Will sort output bam by read name.

Options:

-R   [STR]   Reference path (required)
    Shortcuts:
$ref_shortcuts
			   
-O   [STR]   Output prefix (required)

-1   [STR]   Trimmed read 1 (paired, req)
-2   [STR]   Trimmed read 2 (paired, req)

-t   [INT]   Number of threads for alignment.
-o   [INT]   Number of threads for output.
               (also thread count for sorting by name)

-r   [STR]   Report alignment stats to slack channel
              Requires 'slack' as cli callable
-X           Retain coord sorted bam (def is only name sorted)

Executable Commands (from $DEFAULTS_FILE)
   bsbolt:   $bsbolt
   slack:    $slack
   samtools: $samtools
	
";

$start_time = localtime(time);

if (!defined $opt{'R'}) {
	die "\nERROR: Provide a reference as -R\n$die";
} else {
	if (defined $REF{$opt{'R'}}) {$ref = $REF{$opt{'R'}}}
	else {$ref = $opt{'R'}};
}

if (!defined $opt{'1'} || !defined $opt{'2'}) {die "\nERROR: Reads 1 and 2 MUST be specified!\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'o'}) {$o_threads = $opt{'o'}};

open LOG, ">$opt{'O'}.bsbolt.log";
$ts = localtime(time);
print LOG "$ts\tAlignment called.\n";

#bsbolt Align -F1 ../220727_MB_RPA.trimmed.paired.R2.fq.gz -F2 ../220727_MB_RPA.trimmed.paired.R1.fq.gz -t 44 -OT 2 -T 10 -O bsbolt/220727.bsboltAutoC2T -DB /home/groups/oroaklab/refs/mm10/bsbolt/
$align_call = "$bsbolt Align -F1 $opt{'2'} -F2 $opt{'1'} -t $threads -OT $o_threads -O $opt{'O'} -DB $ref >> $opt{'O'}.bsbolt.log 2>> $opt{'O'}.bsbolt.log";

print LOG "Command: $align_call\n";
system("$align_call");

$ts = localtime(time);
print LOG "$ts\tDone.\n";

if (defined $opt{'r'}) {
	$message = "Alignment complete for $opt{'O'}\nStart time: $start_time\nEnd time: $ts\nCall: $pe_align_call\n";
	system("$slack -F $opt{'O'}.align_report.txt -c \"$message\" $opt{'r'} >/dev/null 2>/dev/null");
}

$ts = localtime(time);
print LOG "$ts\tSorting by read name.\n";

$sort_call = "$samtools sort -@ $o_threads -n -m 2G $opt{'O'}.bam > $opt{'O'}.nsrt.bam";
print LOG "Command: $sort_call\n";
system("$sort_call");

if (!defined $opt{'X'}) {
	print LOG "Deleting $opt{'O'}.bam\n";
	system("rm -f $opt{'O'}.bam");
}

exit;

}

1;