package premethyst_commands::bam_extract;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("bam_extract");

sub bam_extract {

getopts("O:t:sm:G:C:N:P:p:T:xs:L:M:BEH", \%opt);

# dfaults
$minSize = 10000000;
$minReads = 10000;
$maxPct = 100;
$minPct = 0;
$in_threads = 1;
$sleep = 3;
$read_mCHmax = 0.4;

$die = "

premethyst bam-extract (options) -O [out prefix (req)] [rmdup & filtered bam file, name-sorted]
      or   extract

Extracts methylation from a BSBolt bam and outputs a cellCall folder.

Strongly recommended to run with the -C [compelxity file] and -N [cutoff]
options, otherwise noise barcodes will be processed.

Alternatively, the bam file can be pre-filtered to exclude noise
barcodes; however -C/-N is much simpler.

Options:
   -O   [STR]   Out prefix
   -m   [INT]   Minimum chromosome size to retain (def = $minSize)
                  Used to exclude random and other small contigs
   -M   [FLT]   Max allowed fraction mCH sites methylated (def = $read_mCHmax)
                  For brain, rec: 0.4; for non mCH cell types, rec: 0.1
   -B           Methylation call field is in Bismark format (XM:Z:)
                  Eg. Bismark or UA-Meth as aligner.
                  Default assumes BSBolt (XB:Z:)
   -H           Exclude CH context

Threading:
   -t   [INT]   Max number of concurrent extract threads (def = 1)
   -T   [INT]   Number of threads for reading input bam (def = $in_threads)

Cell Filtering:
   -C   [STR]   Complexity.txt file from rmdup (recommened)
   -L   [STR]   List of cell barcodes to include (further filtered if -C)
   -N   [INT]   Minimum unique reads per cell to include (def = $minReads)
                  Note: can be inclusive, additional filtering can
                        be carried out on the calls folder.
   -P   [FLT]   Max percent unique reads to include (def = $maxPct)
   -p   [FLT]   Min percent unique reads to include (def = $minPct)

Other Options:
   -x           Only generate stats, not cellCall files.
   -s   [INT]   Seconds to wait between thread checks (def = $sleep)

Executable Commands (from $DEFAULTS_FILE)
   samtools:   $samtools
   premethyst: $premethyst
   
";

if (!defined $ARGV[0] || (!defined $opt{'O'} && !defined $opt{'s'})) {die $die};
if (defined $opt{'N'}) {$minReads = $opt{'N'}};
if (defined $opt{'P'}) {$maxPct = $opt{'P'}};
if (defined $opt{'p'}) {$minPct = $opt{'p'}};
if (defined $opt{'s'}) {$sleep = $opt{'s'}};
if (defined $opt{'T'}) {$in_threads = $opt{'T'}};
if (defined $opt{'M'}) {$read_mCHmax = $opt{'M'}};
if (defined $opt{'B'}) {$methField = "XM:Z:"} else {$methField = "XB:Z:"};

if (defined $opt{'x'} && defined $opt{'E'}) {
	die "\nERROR: Cannot run with no cov file output (-x) while also exporting the output (-E)!\n\n";
}

if (defined $opt{'L'}) {
	open IN, "$opt{'L'}";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$INCLUDE{$P[1]} = 1;
	}
}

if (defined $opt{'C'}) {
	open IN, "$opt{'C'}";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if ($P[3] >= $minReads && $P[4] <= $maxPct && $P[4] >= $minPct) {
			$INCLUDE{$P[1]} = 1;
		} elsif (defined $INCLUDE{$P[1]}) {
			$INCLUDE{$P[1]} = 0;
		}
	} close IN;
}

%CALL_CONV = ("x" => "z", "X" => "Z", "y" => "x", "Y" => "X", "z" => "h", "Z" => "H");

if (!defined $opt{'s'}) { # main thread
	#defaults
	if (!defined $opt{'t'}) {$opt{'t'} = 1};
	
	$thread_opts = "";
	if (defined $opt{'x'}) {$thread_opts .= "-x "};
	if (defined $opt{'B'}) {$thread_opts .= "-B "};
	if (defined $opt{'H'}) {$thread_opts .= "-H "};
	if (defined $opt{'E'}) {$thread_opts .= "-E "};
	$thread_opts .= "-s";
	
	# open log file
	open LOG, ">$opt{'O'}.log";
	$ts = localtime(time);
	print LOG "$ts\tProgram called.\n\n============== Phase 1: Parsing bam and processing per-cell calls ==============\n\n";

	# setup output directory
	system("mkdir $opt{'O'}");
	$ts = localtime(time);
	print LOG "$ts\tOutput directory created: $opt{'O'}\n";
	
	if (defined $opt{'E'}) {
		print LOG "$ts\tInitializing h5 export streaming (option -E)\n";
		system("mkdir $opt{'O'}.export");
		if (defined $opt{'H'}) {
			system("$premethyst stream-h5 $opt{'O'}.export $opt{'O'} CH &");
		} else {
			system("$premethyst stream-h5 $opt{'O'}.export $opt{'O'} &");
		}
	}
	
	if (defined $opt{'x'}) {
		print LOG "$ts\tWARNING: -x toggled, will NOT output cov files, only stats files.\n";
	}
	
	# setup threading stuff
	%QUEUE = ();
	
	# start parsing bam
	$ts = localtime(time);
	print LOG "$ts\tParsing bam file: $ARGV[0]\n";
	
	open IN, "$samtools view -@ $in_threads $ARGV[0] |";
	$currentBarc = "null";
	$methCol = 0;

	$prevID = "null";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$barc = $P[0]; $barc =~ s/:.+$//;
		$readID = $P[0]; $readID =~ s/#.+$//;
		
		# check if barc is included, ortherwise skip line
		if (!defined $opt{'C'} || $INCLUDE{$barc} > 0) {
			# save barc
			$ALL_CELLS{$barc} = 1;
			
			# check to close out barcode
			if ($currentBarc ne $barc) {
				# close out prev barc
				if ($currentBarc ne "null") {
					close OUT;
					$QUEUE{$currentBarc} = 0;
				}
				# check thread stats
				check_threads();
				
				# set next barcode & open outfile
				$currentBarc = $barc;
				$prevID = "null";
				open OUT, ">$opt{'O'}/$currentBarc.meth";
			}
			# parse read
			$chr = $P[2];
			$pos = $P[3];
			$meth = get_field($methCol,$methField);
			
			if ($prevID eq "null") {
				$prevChr = $chr;
				$prevPos = $pos;
				$prevMeth = $meth;
				$prevID = $readID;
			} elsif ($prevID eq $readID) { # print pair
				print OUT "$prevChr\t$prevPos\t$prevMeth\t$chr\t$pos\t$meth\n";
				$prevID = "null";
			} else { # solo - print previous as solo & save new
				print OUT "$prevChr\t$prevPos\t$prevMeth\n";
				$prevChr = $chr;
				$prevPos = $pos;
				$prevMeth = $meth;
				$prevID = $readID;
			}
			
			#print OUT "$chr\t$pos\t$meth\n";
			
		}
	} close IN;
	
	# finish last one
	close OUT;
	$QUEUE{$currentBarc} = 0;
	
	# process rest of queue
	$incomplete = 1;
	while ($incomplete > 0) {
		$incomplete = 0;
		sleep($sleep); # wait
		check_threads();
		foreach $queued_barc (keys %QUEUE) {
			if ($QUEUE{$queued_barc} < 2) { # active or waiting
				$incomplete++;
			}
		}
	}
	$ts = localtime(time);
	print LOG "$ts\tAll threads complete! Collating CellInfo.\n";

	foreach $barc (keys %ALL_CELLS) {
		if (defined $opt{'x'}) {
			system("cat $opt{'O'}/$barc.complete >> $opt{'O'}.cellInfo.txt");
		} else {
			system("cat $opt{'O'}/$barc.complete >> $opt{'O'}.cellInfo.txt && rm -f $opt{'O'}/$barc.complete $opt{'O'}/$barc.meth");
		}
	}

	exit;
	
} else { # subthread
	open IN, "$ARGV[0]/$ARGV[1].meth";
	if (!defined $opt{'x'}) {
		open CG, "| sort -k 1,1 -k 2,2n > $ARGV[0]/$ARGV[1].CG.cov";
		if (!defined $opt{'H'}) {open CH, "| sort -k 1,1 -k 2,2n > $ARGV[0]/$ARGV[1].CH.cov"};
	}
	$pair_ct = 0; $single_ct = 0; $frag_ct = 0; $excluded_mCH = 0;
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		
		%read_cov = (); %read_meth = ();
		$read_mCH = 0; $read_CH = 0;
		if (defined $P[3]) { # pair
			load_read($P[0],$P[1],$P[2]);
			load_read($P[3],$P[4],$P[5]);
			$pair_ct++;
		} else {
			load_read($P[0],$P[1],$P[2]);
			$single_ct++;
		}
		
		if ($read_CH>0) {
			if (($read_mCH/$read_CH)<=$read_mCHmax) {
				foreach $coord (keys %read_cov) {
					if (!defined $COORD_cov{$coord}) {
						if (defined $read_cov{$coord}{'x'}) {$COORD_cov{$coord}{'x'} = 1};
						if (defined $read_cov{$coord}{'h'}) {$COORD_cov{$coord}{'h'} = 1};
						if (defined $read_meth{$coord}{'x'}) {$COORD_meth{$coord}{'x'} = 1};
						if (defined $read_meth{$coord}{'h'}) {$COORD_meth{$coord}{'h'} = 1};
					} else {
						if (defined $read_cov{$coord}{'x'}) {$COORD_cov{$coord}{'x'}++};
						if (defined $read_cov{$coord}{'h'}) {$COORD_cov{$coord}{'h'}++};
						if (defined $read_meth{$coord}{'x'}) {$COORD_meth{$coord}{'x'}++};
						if (defined $read_meth{$coord}{'h'}) {$COORD_meth{$coord}{'h'}++};
					}
				}
			}
		} else {
			$excluded_mCH++;
		}
		
	} close IN;
	
	foreach $coord (keys %COORD_cov) {
		if (defined $COORD_cov{$coord}{'x'}) {
			if (!defined $COORD_meth{$coord}{'x'}) {$COORD_meth{$coord}{'x'}=0};
			$pct = sprintf("%.2f", ($COORD_meth{$coord}{'x'}/$COORD_cov{$coord}{'x'})*100);
			$unmeth = $COORD_cov{$coord}{'x'}-$COORD_meth{$coord}{'x'};
			if (!defined $opt{'x'}) {print CG "$coord\t$pct\t$unmeth\t$COORD_meth{$coord}{'x'}\n"};
			$CG_cov++; $CG_bases+=$COORD_cov{$coord}{'x'}; $CG_meth+=$COORD_meth{$coord}{'x'};
		}
		if (defined $COORD_cov{$coord}{'h'}) {
			if (!defined $COORD_meth{$coord}{'h'}) {$COORD_meth{$coord}{'h'}=0};
			$pct = sprintf("%.2f", ($COORD_meth{$coord}{'h'}/$COORD_cov{$coord}{'h'})*100);
			$unmeth = $COORD_cov{$coord}{'h'}-$COORD_meth{$coord}{'h'};
			if (!defined $opt{'x'} && !defined $opt{'H'}) {print CH "$coord\t$pct\t$unmeth\t$COORD_meth{$coord}{'h'}\n"};
			$CH_cov++; $CH_bases+=$COORD_cov{$coord}{'h'}; $CH_meth+=$COORD_meth{$coord}{'h'};
		}
	}
	
	close CG; close CH;

	$AllCov = $CG_bases+$CH_bases;
	$frag_ct = $single_ct+$pair_ct;
	if ($CG_bases>0) {$CGpct = sprintf("%.2f", ($CG_meth/$CG_bases)*100)} else {$CGpct = "0.00"};
	if ($CH_bases>0) {$CHpct = sprintf("%.2f", ($CH_meth/$CH_bases)*100)} else {$CHpct = "0.00"};
	
	$cellInfo = "$ARGV[1]\t$AllCov\t$CG_cov\t$CGpct\t$CH_cov\t$CHpct\t$frag_ct\t$pair_ct\t$single_ct\t$excluded_mCH";

	system("echo '$cellInfo' > $ARGV[0]/$ARGV[1].complete");
	
	if (defined $opt{'E'}) {
		system("mv $ARGV[0]/$ARGV[1].CG.cov $ARGV[0].export/$ARGV[1].CG.cov");
		if (!defined $opt{'H'}) {system("mv $ARGV[0]/$ARGV[1].CH.cov $ARGV[0].export/$ARGV[1].CH.cov")};
	}
	
	exit;
}

sub get_field {
	($init_field,$flag) = @_;
	if ($P[$init_field] =~  /^$flag/) {
		$field_out = $P[$init_field];
		$field_out =~ s/^.+://;
		return($init_field,$field_out);
	} else {
		for ($fieldID = 10; $fieldID < @P; $fieldID++) {
			if ($P[$fieldID] =~  /^$flag/) {
				$field_out = $P[$fieldID];
				$field_out =~ s/^.+://;
				return($fieldID,$field_out);
			}
		}
	}
}

sub load_read {
	$chr = $_[0]; $pos = $_[1]; $meth = $_[2];
	if ($pos > 0) { # aligned
		@M = split(//, $meth);
		if (defined $opt{'B'}) { # Bismark format
			for ($i = 0; $i < @M; $i++) {
				$coord = "$chr\t$pos";
				if ($M[$i] eq "z") { # CG, unmeth
					$read_cov{$coord}{'x'} = 1;
				} elsif ($M[$i] eq "Z") { # CG, meth
					$read_cov{$coord}{'x'} = 1;
					$read_meth{$coord}{'x'} = 1;
				} elsif ($M[$i] eq "x" || $M[$i] eq "h") { # CH , unmeth
					$read_cov{$coord}{'h'} = 1;
					$read_CH++;
				} elsif ($M[$i] eq "X" || $M[$i] eq "H") { # CH , meth
					$read_cov{$coord}{'h'} = 1;
					$read_meth{$coord}{'h'} = 1;
					$read_mCH++; $read_CH++;
				}
				$pos++;
			}
		} else { # BSBolt format
			for ($i = 0; $i < @M; $i++) {
				if ($M[$i] !~ /[0-9]/) {
					
					$coord = "$chr\t$pos";
					if ($M[$i] eq "x") { #CG, unmeth
						$read_cov{$coord}{'x'} = 1;
					} elsif ($M[$i] eq "X") { #CG, meth
						$read_cov{$coord}{'x'} = 1;
						$read_meth{$coord}{'x'} = 1;
					} elsif ($M[$i] eq "y" || $M[$i] eq "z") { # CH , unmeth
						$read_cov{$coord}{'h'} = 1;
						$read_CH++;
					} elsif ($M[$i] eq "Y" || $M[$i] eq "Z") { # CH , meth
						$read_cov{$coord}{'h'} = 1;
						$read_meth{$coord}{'h'} = 1;
						$read_mCH++; $read_CH++;
					}
					
					$pos++;
					
				} elsif ($M[$i] =~ /[0-9]/) {
					$add = $M[$i];
					while ($M[$i+1] !~ /xyz/i && $M[$i+1] =~ /[0-9]/) {
						$i++;
						$add .= $M[$i];
					}
					$pos += $add;
				}
			}
		}
	}
}

sub check_threads {
	$running = 0; @WAITING = ();
	foreach $queued_barc (keys %QUEUE) {
		if ($QUEUE{$queued_barc} == 1) { # active
			$running++;
			if (-e "$opt{'O'}/$queued_barc.complete") { # finished
				$running--;
				$QUEUE{$queued_barc} = 2;
				print LOG "\t\t\t$queued_barc COMPLETED!\n";
			}
		} elsif ($QUEUE{$queued_barc} == 0) { # store to kick off
			push @WAITING, $queued_barc;
		}
	}
	if ($running < $opt{'t'}) { # if not at full threading
		$ts = localtime(time);
#		print LOG "\t$ts\tThreads active: $running, adding...\n";
		for ($waitingID = 0; $waitingID < ($opt{'t'} - $running); $waitingID++) { # see how many short
			if (defined $QUEUE{$WAITING[$waitingID]}) { # if wiating is occupied
				system("$premethyst bam-extract $thread_opts $opt{'O'} $WAITING[$waitingID] &"); # start the thread
				$QUEUE{$WAITING[$waitingID]} = 1; # log it as active
				print LOG "\t\t$WAITING[$waitingID] now running.\n";
			}
		}
	}
}

}

1;