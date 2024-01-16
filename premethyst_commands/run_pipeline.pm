package premethyst_commands::run_pipeline;

use Cwd;
use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("run_pipeline");

sub run_pipeline {
	
getopts("C:O:S:E:r:d", \%opt);

$out = "PremethystOut";
$start = "1";
$end = "7";

%STEPS = ('0' => 'unidex', '1' => 'trim',
          '2' => 'align', '3' => 'rmdup',
          '4' => 'plot', '5' => 'extract',
          '6' => 'filter', '7' => 'export');

$die = "

premethyst run-pipeline (options) [run.cfg]
      or   run

Will run the full pipeline from demultiplexed fastq files through
generation of final h5 file for amethyst processing.

Options:
  -C  [run.cfg] Print out an example configuration file with
                  recommended base parameters.
  -O  [STR]     Output prefix to output for the example
                  cfg file (def = $out)
  -r  [STR]     Report progress to specified slack channel (opt)
                 (Requires slack as cli callable)
  -d            Debug mode: Generate fodlers & log file with list
                  of commands that would be run, just do not
                  execute the actual commands.

  -S  [0-6]     Step to start pipeline (def = $start)
  -E  [1-7]     Step to end after (def = $end)

Steps:
 (0 - unidex)    4 - plot
  1 - trim       5 - extract
  2 - align      6 - filter
  3 - rmdup      7 - export
  
Executable Commands (from $DEFAULTS_FILE)
   premethyst: $premethyst
   slack:      $slack

";

if ((defined $opt{'C'} || defined $opt{'O'}) && defined $ARGV[0]) {
	die "\nERROR: Cannot run with a cfg file and output a cfg file in the same run.\n$die";
} elsif (!defined $opt{'C'} && !defined $ARGV[0]) {die $die};

# Print sample config and exit.
if (defined $opt{'C'}) {
	if (defined $opt{'O'}) {$out = $opt{'O'}};
	print_cfg();
	exit;
}

# Pipeline mode

# Read config and make sure required O is present, set start and end
%OPTS = ();
read_cfg();
if (!defined $OPTS{'global'}{'O'}) {
	die "\nERROR: Config file must contain global output prefix as 'global_O=[pfx]'\n$die";
}
if (!defined $OPTS{'global'}{'D'}) {
	die "\nERROR: Config file must contain global output folder as 'global_D=[pfx]'\n$die";
}
if (-d "$OPTS{'global'}{'D'}") {
	print STDERR "\nWARNING: $OPTS{'global'}{'D'} directory exists! Will overwrite contents unless -O is new. Proceeding in 30 seconds.\n";
	sleep(30)
} else {
	system("mkdir $OPTS{'global'}{'D'}");
}


# set start and end
if (defined $opt{'S'}) {$start = $opt{'S'}};
if (defined $opt{'E'}) {$end = $opt{'E'}};
if ($start < 0 || $start > 6) {die "\nERROR: starting step (-S) must be 0-6.\n$die"};
if ($end < 1 || $end > 7) {die "\nERROR: ending step (-E) must be 1-7.\n$die"};

# start log file
$ts = localtime(time);
if (-e "$OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.premethyst.log") {
	open LOG, ">>$OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.premethyst.log";
	print LOG "
=======================================================
========== $ts premethyst run, adding to log.\n";
} else {
	open LOG, ">$OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.premethyst.log";
	print LOG "========== $ts premethyst run initiatied.\n";
}
print LOG "\tCONFIG_FILE=$ARGV[0]
\tSTART_STEP=$start
\tEND_STEP=$end
Parameters specified in cfg file:\n";

foreach $step (keys %OPTS) {
	foreach $opt (keys %{$OPTS{$step}}) {
		print LOG "\t$step -$opt $OPTS{$step}{$opt}\n";
	}
}

# STEP 0: UNIDEX
if ($start==0) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 0: UNIDEX\n";
	if (!defined $OPTS{'unidex'}{'R'} ||
		(!defined $OPTS{'unidex'}{'1'} &&
		 !defined $OPTS{'unidex'}{'2'} &&
		 !defined $OPTS{'unidex'}{'3'} &&
		 !defined $OPTS{'unidex'}{'4'})) {
		print LOG "ERROR: -R or all four read files (-1 to -4) must be specified!\n";
		exit;
	}
	if (!defined $OPTS{'unidex'}{'M'}) {
		print LOG "ERROR: -M [mode(s)] must be provided.\n";
		exit;
	}
	$options = build_options("unidex");
	$command = "unidex $options";
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 0: UNIDEX.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
}

# STEP 1: TRIM
if ($start <= 1) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 1: TRIM\n";
	if (!defined $OPTS{'trim'}{'1'} ||
		!defined $OPTS{'trim'}{'2'}) {
		print LOG "ERROR: Reads 1 and 2, -1 and -2 must be provided. This is true even if unidex is the first step.\n";
		exit;
	}
	$options = build_options("trim");
	$command = "premethyst trim -O $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'} $options";
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 1: TRIM.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
}

# STEP 2: ALIGN
if ($start <= 2 && $end >= 2) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 2: ALIGN\n";
	if (!defined $OPTS{'align'}{'R'}) {
		print LOG "ERROR: Reference (R) must be provided.\n";
		exit;
	}
	$options = build_options("align");
	$command = "$premethyst align $options -O $OPTS{'global'}{'O'} -1 $OPTS{'global'}{'O'}.trimmed.paired.R1.fq.gz -2 $OPTS{'global'}{'O'}.trimmed.paired.R2.fq.gz";
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 2: ALIGN.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
}

# STEP 3: RMDUP
if ($start <= 3 && $end >= 3) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 3: RMDUP\n";
	if (defined $OPTS{'rmdup'}) {
		$options = build_options("rmdup");
		$command = "$premethyst rmdup $options $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.nsrt.bam";
	} else {
		$command = "$premethyst rmdup $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.nsrt.bam";
	}
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 3: RMDUP.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
}

# STEP 4: PLOT COMPLEXITY
if ($start <= 5 && $end >= 5) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 4: PLOT COMPLEXITY\n";
	if (defined $OPTS{'cplx'}) {
		$options = build_options("cplx");
		$command = "$premethyst complexity $options $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.txt";
	} else {
		$command = "$premethyst complexity $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.txt";
	}
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 4: PLOT COMPLEXITY.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
	if (defined $opt{'r'}) {system("$slack -F $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.png $opt{'r'}")};
}

# STEP 5: EXTRACT
if ($start <= 5 && $end >= 5) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 5: EXTRACT\n";
	if (!defined $OPTS{'extract'}{'C'}) {
		if (defined $OPTS{'extract'}{'N'} ||
			defined $OPTS{'extract'}{'P'} ||
			defined $OPTS{'extract'}{'p'}) {
			if (-e "$OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.txt" || defined $opt{'d'}) {
				print LOG "Complexity file auto-detected as $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.txt.\n";
				$OPTS{'extract'}{'C'} = "$OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.txt";
			} elsif (!defined $opt{'d'}) {
				print LOG "ERROR: Cannot find complexity file ($OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.complexity.txt) to perform filtering, but filter parameters are specified.\n";
				exit;
			}
		}
	}
	if (defined $OPTS{'extract'}) {
		$options = build_options("extract");
		$command = "$premethyst extract $options $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.bbrd.q10.nsrt.bam";
	} else {
		$command = "$premethyst extract $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.bbrd.q10.nsrt.bam";
	}
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 5: EXTRACT.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
}

# STEP 6: FILTER ETC...
if ($start <= 6 && $end >= 6) {
	# rename
	if (defined $OPTS{'rename'}) {
		$ts = localtime(time);
		print LOG "\n========== $ts STEP 6a: RENAME\n";
		$options = build_options("rename");
		$command = "$premethyst rename $options $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}";
		print LOG "Command: $command\n";
		if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 6a: RENAME.\nCommand:\n\t$command\" $opt{'r'}")};
		if (!defined $opt{'d'}) {system("$command")};
	}
	# extract additional context
	if (defined $OPTS{'context'}) {
		$ts = localtime(time);
		print LOG "\n========== $ts STEP 6b: CONTEXT EXTRACT\n";
		if (!defined $OPTS{'context'}{'C'} ||
			!defined $OPTS{'context'}{'C'} ||
			!defined $OPTS{'context'}{'C'}) {
			print LOG "ERROR: If extracting an additional context, options C, N and B must be provided.\n";
			exit;
		}
		$options = build_options("context");
		$command = "$premethyst context $options $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}";
		print LOG "Command: $command\n";
		if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 6b: CONTEXT EXTRACT.\nCommand:\n\t$command\" $opt{'r'}")};
		if (!defined $opt{'d'}) {system("$command")};
	}
	# filter
	if (defined $OPTS{'filter'}) {
		$ts = localtime(time);
		print LOG "\n========== $ts STEP 6c: FILTER\n";
		if (-e "$OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.cellInfo.txt") {
			print LOG "Using $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.cellInfo.txt for filtering.\n";
		} else {
			print LOG "ERROR: cannot find $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.cellInfo.txt for filtering.\n";
		}
		if (defined $OPTS{'rename'}{'P'} && !defined $OPTS{'filter'}{'p'}) {
			print LOG "Prefix addition was specified in rename step, toggling -p for filtering.\n";
			$OPTS{'filter'}{'p'} = "true";
		}
		if (defined $OPTS{'rename'}{'S'} && !defined $OPTS{'filter'}{'s'}) {
			print LOG "Suffix addition was specified in rename step, toggling -s for filtering.\n";
			$OPTS{'filter'}{'s'} = "true";
		}
		$options = build_options("filter");
		$command = "$premethyst filter $options $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.cellInfo.txt $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'} $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}.filtered";
		print LOG "Command: $command\n";
		if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 6c: FILTER.\nCommand:\n\t$command\" $opt{'r'}")};
		if (!defined $opt{'d'}) {system("$command")};
	}
}

# STEP 7: EXPORT
if ($start <= 7 && $end >= 7) {
	$ts = localtime(time);
	print LOG "\n========== $ts STEP 7: EXPORT\n";
	$command = "$premethyst calls2h5 $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'} $OPTS{'global'}{'D'}/$OPTS{'global'}{'O'}";
	print LOG "Command: $command\n";
	if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], starting STEP 7: EXPORT.\nCommand:\n\t$command\" $opt{'r'}")};
	if (!defined $opt{'d'}) {system("$command")};
}

$ts = localtime(time);
print LOG "\n========== $ts COMPLETE!\n";
close LOG;
if (defined $opt{'r'}) {system("$slack -c \"$ts premethyst run $ARGV[0], COMPLETE!\" $opt{'r'}")};

#### SUBROUTINES

sub read_cfg {
	open IN, "$ARGV[0]";
	while ($l = <IN>) {
		if ($l !~ /^#/) {
			chomp $l;
			($stepOption,$val) = split(/=/, $l);
			($stepName,$option) = split(/_/, $stepOption);
			$OPTS{$stepName}{$option} = $val;
		}
	} close IN;
}

sub build_options {
	$build = "";
	$step = $_[0];
	foreach $opt (keys %{$OPTS{$step}}) {
		if ($OPTS{$step}{$opt} =~ /true/i) {
			$build .= " -$opt";
		} else {
			$build .= " -$opt $OPTS{$step}{$opt}";
		}
	}
	$build =~ s/^\s//;
	return $build;
}

sub print_cfg {
open OUT, ">$opt{'C'}";
$cwd = cwd();
print OUT "## PREMETHYST RUN CONFIG FILE
##
## Alter parameters to your own system and experiment parameters.
## If generated with 'premethyst run -C' this file will contain
## current working directory locations.
##
## Each step in the pipeline will use defaults unless specified.
## To alter any one of the additional options availble to the
## step, it can be added by including aline with [stepName]_[optID]
## For example, 'trim_m=50' would set the -m option in the trim
## command (which sets min read length) to 50.
## For paramters that do not have a value, set to true, e.g:
## 'trim_u=true', which would toggle -u (retain unpaired reads)
## for the trim step.
##
## GLOBAL PARAMETERS
## Output directory (D) and prefix (O) used for all outputs
## Do not specify O option for individual commands.
global_D=$cwd/$out
global_O=$out
#
## STEP 0: UNIDEX (optional; stepName=unidex)
## If using unidex to demultiplex, the parameters can be set in
## this block.
#unidex_R=[myRunID]
#unidex_r=[my/path/toR]
#unidex_M=sciMET
#unidex_o=$cwd
#
## STEP 1: TRIMMING (stepName=trim)
## read 1 and 2 fastq files (required)
trim_1=$out.sciMET.R1.fq.gz
trim_2=$out.sciMET.R2.fq.gz
#
## STEP 2: ALIGNMENT (stepName=align)
## Do not specify fastq files, they are inferred from previous step.
## Reference genome (required):
align_R=hg38
## Alignment threads
align_t=24
## Output threads (separate from alignment threads)
align_o=8
#
## STEP 3: RMDUP (stepName=rmdup)
## Use defaults.
#
## STEP 4: PLOT COMPLEXITY (stepName=cplx)
## Use defaults.
#
## STEP 5: EXTRACT CALLS (stepName=extract)
## Filter to min read count (recommended)
extract_N=10000
## Threads, recommended to use no more than 12.
extract_t=12
#
## STEP 6: MODIFY CALL FOLDER (optional steps)
## 6A) Rename cellIDs (stepName=rename)
# Prefix or suffix:
#rename_P=PREFIXX
#rename_S=SUFFIX
## 6B) Extract additional context (stepName=context)
## Searched context (req)
#context_C=CH
## New context name (req)
#context_N=CAC
## New context bed of sites (req)
#context_B=[bed_of_sites]
## 6C) Filter cellIDs (stepName=filter)
## Filter to max CH of 1%
#filter_h=0,1
#
## STEP 7: EXPORT (stepName=export)
## No options
#
";
}

}

1;