package premethyst_commands::calls_filter;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("calls_filter");

sub calls_filter {
	
getopts("C:H:G:h:g:ps", \%opt);

$die = "

premethyst calls-filter (options) [cellInfo.txt file (req)] [folder with cellcall files(req)] [folder to deposit filtered cell files (req)]
      or   filter

Uses the cellInfo file to filter cellCell files and move those missing
certain criteria to a separate calls folder.

This is only to reduce file size or remove contexts etc... All filters
can also be applied later when performing analysis with amethyst.

Options:

   -C   [STR],[STR], ...  Contexts to retain. If not specified, all context
                  files will be kept. Canbe multiple comma separated.
   -H   [INT]   Minimum CH sites covered.
   -G   [INT]   Minimum CG sites covered.
   -h   [MIN,MAX]  Min & Max CH global methylation.
   -g   [MIN,MAX]  Min & Max CG global methylation.
   -p           Prefix is present in CellID names in cellcall files.
                  (will ignore everything before '-' when matching to cellinfo file)
   -s           Suffix is present in CellID names in cellcall files.
                  (will ignore everything after '-' when matching to cellinfo file)
				   
";

if (!defined $ARGV[2]) {die $die};

if (defined $opt{'h'}) {
	($minH,$maxH) = split(/,/, $opt{'h'});
}
if (defined $opt{'g'}) {
	($minG,$maxG) = split(/,/, $opt{'g'});
}

if (defined $opt{'C'}) {
	@CONTEXTS = split(/,/, %opt{'C'});
	foreach $context (@CONTEXTS) {
		$CONTEXT_keep{$context} = 1;
	}
}

open LOG, ">> $ARGV[2].filt.log";

$ts = localtime(time);

print LOG "

$ts Filtering $ARGV[1] using cellInfo $ARGV[0], output to $ARGV[2].
\tOptions:
";
foreach $option (keys %opt) {
	print LOG "\t\t$option\t$opt{$option}\n";
}

if (-d "$ARGV[2]") {
	print LOG "Fail directory already exists and will be used: $ARGV[2]\n";
} else {
	print LOG "Creating fail directory: $ARGV[2]\n";
}

open INFO, "$ARGV[0]";
while ($l = <INFO>) {
	chomp $l;
	($cellID,$AllCov,$CG_cov,$CGpct,$CH_cov,$CHpct) = split(/\t/, $l);
	$status = 'P';
	if (defined $opt{'H'} && $CH_cov<$opt{'H'}) {$status = 'F'};
	if (defined $opt{'G'} && $CG_cov<$opt{'G'}) {$status = 'F'};
	if (defined $opt{'h'} && ($CHpct > $maxH || $CHpct < $minH)) {$status = 'F'};
	if (defined $opt{'g'} && ($CGpct > $maxG || $CGpct < $minG)) {$status = 'F'};
	$CELLID_status{$cellID} = $status;
} close INFO;

opendir(CCD, "$ARGV[0]");
@CELLCALLS = readdir(CCD);
closedir(CCD);
$pass = 0; $fail = 0;
foreach $cellCallFile (@CELLCALLS) {
	if ($cellCallFile =~ /\.cov/) {
		($cellID,$context,$cov) = split(/\./, $cellCallFile);
		
		if ($pass == 0 && $fail == 0) {
			print LOG "First cellID pre prefix/suffix trimming:  $cellID\n";
		}
		
		if (defined $opt{'p'}) {$cellID =~ s/^.+-//};
		if (defined $opt{'s'}) {$cellID =~ s/-.+$//};
		
		if ($pass == 0 && $fail == 0) {
			print LOG "First cellID post prefix/suffix trimming: $cellID\n";
		}
		
		$status = 'P';
		if (defined $opt{'C'} && !defined $CONTEXT_keep{$context}) {
			$fail++;
			$status = 'F';
		} else {
			if (!defined $CELLID_status{$cellID}) {
				print LOG "WARNING: CellID $cellID not present in cellInfo file! Retaining by default.\n";
			} else {
				$status = $CELLID_status{$cellID};
			}
		}
		
		if ($status eq 'F') {
			system("mv $ARGV[1]/$cellCallFile $ARGV[2]/$cellCallFile");
		}
	}
}

$totalFiles = $pass+$fail;
$ts = localtime(time);
print LOG "$ts Filtering complete. $totalFiles cov files processed
\t$pass passing and remain in $ARGV[1]
\t$fail filtered and moved to $ARGV[2]\n\n";

exit;

}

1;