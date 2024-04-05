package premethyst_commands::ultima2fastq;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("ultima2fastq");

sub ultima2fastq {

getopts("C:O:a:b:A:B:t:T:vr", \%opt);

$threads = 4;
$o_threads = 2;
$BA_ED = 2;
$BB_ED = 2;

$die = "

premethyst ultima2fastq (options) -C [ultima cram file] -O [output prefix]

Takes in a sample-split unaligned UG cram file and edit-distance
matches to expected barcode sequences then outputs a fastq file
with the error-corrected barcodes in a 'sci' format.

Options:

-C   [STR]   Ultima cram file from UG pipeline (required)
-O   [STR]   Output prefix (required)

-r           Error-corrected index already provided as rx:Z: field
              Just reformat using the rx field
-a   [STR]   Index list expected for 'ba' barcode field
              (required if not -r)
-b   [STR]   Index list expected for 'bb' barcode field
              (required if not -r)

-A   [INT]   Max matching edit distance for ba index (def = $BA_ED)
-B   [INT]   Max matching edit distance for bb index (def = $BB_ED)

-t   [INT]   Threads for cram read input (def = $threads)
-T   [INT]   Threads for fastq output (def = $o_threads)

Executable Commands (from $DEFAULTS_FILE)
   samtools: $samtools
   pigz:     $pigz
   
";

if (!defined $opt{'C'}) {die "\nERROR: CRAM MUST be specified!\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (!defined $opt{'r'} && (!defined $opt{'a'} || !defined $opt{'b'})) {die "\nERROR: If edit distance matching mode (ie no rx field), -a and -b must be specified!\n$die"};

if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'T'}) {$o_threads = $opt{'T'}};
if (defined $opt{'A'}) {$BA_ED = $opt{'A'}};
if (defined $opt{'B'}) {$BB_ED = $opt{'B'}};
if (defined $opt{'a'}) {$BA_LIST = $opt{'a'}};
if (defined $opt{'b'}) {$BB_LIST = $opt{'b'}};

if (defined $opt{'r'}) { # error correction done - just reformat
	open IN, "samtools view -@ $threads $opt{'C'} |";
	open PASS, "| $pigz -p $o_threads > $opt{'O'}.fq.gz";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$readCT++;
		
		# pull indexes & RG
		$RG = "null"; $RX = "null";
		($rg_field_new,$RG) = get_field($rg_field,"RG:Z:"); $rg_field = $rg_field_new;
		($rx_field_new,$RX) = get_field($rx_field,"rx:Z:"); $rx_field = $rx_field_new;
		
		if ($RG ne "null" && $RX ne "null") {
			$RX =~ s/-/_/;
			print PASS "\@$RG\_$RX:$readCT\n$P[9]\n\+\n$P[10]\n";
		}
	} close IN; close PASS;
	exit;
} else {

$time = localtime(time);
open LOG, ">$opt{'O'}.ultima2fastq.log";
print LOG "premethyst ultima2fastq called at $time
Edit distance error correction mode toggled.
Input cram: $opt{'C'}
Output Prefix: $opt{'O'}
BA File: $BA_LIST, ED: $BA_ED
BB File: $BB_LIST, ED: $BB_ED
In Threads: $threads
Out Threads: $o_threads
Running ...
";


# read in index whitelists
open BA, "$BA_LIST";
while ($l = <BA>) {
	chomp $l;
	@P = split(/t/, $l);
	$seq = uc($P[-1]);
	if ($seq ne "") {
		$BA_COUNT{$seq} = 0;
		$BA_MATCH{$seq} = $seq;
		$ba_total++;
	}
} close BA;

open BB, "$BB_LIST";
while ($l = <BB>) {
	chomp $l;
	@P = split(/t/, $l);
	$seq = uc($P[-1]);
	if ($seq ne "") {
		$BB_COUNT{$seq} = 0;
		$BB_MATCH{$seq} = $seq;
		$bb_total++;
	}
} close BB;

if (defined $opt{'v'}) {print STDERR "$ba_total total BA indexes, $bb_total total BB indexes\n"};

$ba_field = 0; $bb_field = 0; $rg_field = 0;
$passCT = 0; $failCT = 0; $readCT = 0;
$ba_exact = 0; $bb_exact = 0;
%BA_FAIL = (); %BB_FAIL = ();
open IN, "samtools view -@ $threads $opt{'C'} |";
open PASS, "| $pigz -p $o_threads > $opt{'O'}.fq.gz";
open FAIL, "| $pigz -p $o_threads > $opt{'O'}.fail.fq.gz";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$readCT++;
	
	# pull indexes & RG
	$RG = "null"; $ba_seq = "null"; $bb_seq = "null";
	($rg_field_new,$RG) = get_field($rg_field,"RG:Z:"); $rg_field = $rg_field_new;
	($ba_field_new,$ba_seq) = get_field($ba_field,"ba:Z:"); $ba_field = $ba_field_new;
	($bb_field_new,$bb_seq) = get_field($bb_field,"bb:Z:"); $bb_field = $bb_field_new;
	
	$ba_pass = 0; $bb_pass = 0;
	if ($RG ne "null" && $ba_seq ne "null" && $bb_seq ne "null") {
	
		if (defined $BA_MATCH{$ba_seq}) {
			$ba_corrected = $BA_MATCH{$ba_seq};
			$BA_COUNT{$ba_corrected}++;
			$ba_pass++;
			if (defined $BA_COUNT{$ba_seq}) {$ba_exact++};
			if (defined $opt{'v'}) {print STDERR "BA: $ba_seq prevID match!\n"};
		} elsif (!defined $BA_FAIL{$ba_seq}) { # no match recorded - check ED
			%checked = (); $found_match = "null";
			$pos123 = substr($ba_seq,0,3);
			# run for perfect latch of 1st 3bp
			foreach $check_seq (grep {/^$pos123/} keys %BA_COUNT) {
				if ($found_match eq "null") {
					$checked{$check_seq} = edit_distance($check_seq,$ba_seq);
				}
				if ($checked{$check_seq} <= $BA_ED) {
					$found_match = $check_seq;
					if (defined $opt{'v'}) {print STDERR "BA:	123 $ba_seq matches $found_match ED = $checked{$check_seq}\n"};
					last;
				}
			}
			if ($found_match eq "null") { # for perfect 1st 2 bp
				$pos12 = substr($ba_seq,0,2);
				foreach $check_seq (grep {/^$pos12/} keys %BA_COUNT) {
					if ($found_match eq "null" && !defined $checked{$check_seq}) {
						$checked{$check_seq} = edit_distance($check_seq,$ba_seq);
					}
					if ($checked{$check_seq} <= $BA_ED) {
						$found_match = $check_seq;
						if (defined $opt{'v'}) {print STDERR "BA:		12 $ba_seq matches $found_match ED = $checked{$check_seq}\n"};
						last;
					}
				}
			}
			if ($found_match eq "null") { # for perfect 1st 1 bp
				$pos1 = substr($ba_seq,0,1);
				foreach $check_seq (grep {/^$pos1/} keys %BA_COUNT) {
					if ($found_match eq "null" && !defined $checked{$check_seq}) {
						$checked{$check_seq} = edit_distance($check_seq,$ba_seq);
					}
					if ($checked{$check_seq} <= $BA_ED) {
						$found_match = $check_seq;
						if (defined $opt{'v'}) {print STDERR "BA:			1 $ba_seq matches $found_match ED = $checked{$check_seq}\n"};
						last;
					}
				}
			}
			if ($found_match eq "null") { # for mismatched first
				foreach $check_seq (keys %BA_COUNT) {
					if ($found_match eq "null" && !defined $checked{$check_seq}) {
						$checked{$check_seq} = edit_distance($check_seq,$ba_seq);
					}
					if ($checked{$check_seq} <= $BA_ED) {
						$found_match = $check_seq;
						if (defined $opt{'v'}) {print STDERR "BA:				0 $ba_seq matches $found_match ED = $checked{$check_seq}\n"};
						last;
					}
				}
			}
			if ($found_match ne "null") {
				$BA_MATCH{$ba_seq} = $found_match;
				$ba_corrected = $found_match;
				$BA_COUNT{$ba_corrected}++;
				$ba_pass++;
			} else {
				$BA_FAIL{$ba_seq} = 1;
			}
		}
		
		if (defined $BB_MATCH{$bb_seq}) {
			$bb_corrected = $BB_MATCH{$bb_seq};
			$BB_COUNT{$bb_corrected}++;
			$bb_pass++;
			if (defined $BB_COUNT{$bb_seq}) {$bb_exact++};
			if (defined $opt{'v'}) {print STDERR "BB: $bb_seq prevID match!\n"};
		} elsif (!defined $BB_FAIL{$bb_seq}) { # no match recorded - check ED
			%checked = (); $found_match = "null";
			
			$pos123 = substr($bb_seq,0,3);
			# run for perfect latch of 1st 3bp
			foreach $check_seq (grep {/^$pos123/} keys %BB_COUNT) {
				if ($found_match eq "null") {
					$checked{$check_seq} = edit_distance($check_seq,$bb_seq);
				}
				if ($checked{$check_seq} <= $bb_corrected) {
					$found_match = $check_seq;
					if (defined $opt{'v'}) {print STDERR "BB:	123 $bb_seq matches $found_match ED = $checked{$check_seq}\n"};
					last;
				}
			}
			if ($found_match eq "null") { # for perfect 1st 2 bp
				$pos12 = substr($bb_seq,0,2);
				foreach $check_seq (grep {/^$pos12/} keys %BB_COUNT) {
					if ($found_match eq "null" && !defined $checked{$check_seq}) {
						$checked{$check_seq} = edit_distance($check_seq,$bb_seq);
					}
					if ($checked{$check_seq} <= $BB_ED) {
						$found_match = $check_seq;
						if (defined $opt{'v'}) {print STDERR "BB:		12 $bb_seq matches $found_match ED = $checked{$check_seq}\n"};
						last;
					}
				}
			}
			if ($found_match eq "null") { # for perfect 1st 1 bp
				$pos1 = substr($bb_seq,0,1);
				foreach $check_seq (grep {/^$pos1/} keys %BB_COUNT) {
					if ($found_match eq "null" && !defined $checked{$check_seq}) {
						$checked{$check_seq} = edit_distance($check_seq,$bb_seq);
					}
					if ($checked{$check_seq} <= $BB_ED) {
						$found_match = $check_seq;
						if (defined $opt{'v'}) {print STDERR "BB:			1 $bb_seq matches $found_match ED = $checked{$check_seq}\n"};
						last;
					}
				}
			}
			if ($found_match eq "null") { # for mismatched first
				foreach $check_seq (keys %BB_COUNT) {
					if ($found_match eq "null" && !defined $checked{$check_seq}) {
						$checked{$check_seq} = edit_distance($check_seq,$bb_seq);
					}
					if ($checked{$check_seq} <= $BB_ED) {
						$found_match = $check_seq;
						if (defined $opt{'v'}) {print STDERR "BB:				0 $bb_seq matches $found_match ED = $checked{$check_seq}\n"};
						last;
					}
				}
			}
			if ($found_match ne "null") {
				$BB_MATCH{$bb_seq} = $found_match;
				$bb_corrected = $found_match;
				$BB_COUNT{$bb_corrected}++;
				$bb_pass++;
			} else {
				$BB_FAIL{$bb_seq} = 1;
			}
		}
	}
	
	if ($ba_pass>0 && $bb_pass>0) {
		print PASS "\@$RG-$ba_corrected-$bb_corrected:$readCT\n$P[9]\n\+\n$P[10]\n";
		$passCT++;
	} else {
		if ($ba_pass<1) {$BA_FAIL{$ba_seq}++; if (defined $opt{'v'}) {print STDERR "BA: FAIL! $ba_seq\n"}};
		if ($bb_pass<1) {$BB_FAIL{$bb_seq}++; if (defined $opt{'v'}) {print STDERR "BB: FAIL! $bb_seq\n"}};
		print FAIL "\@$RG-$ba_seq-$bb_seq:$readCT\n$P[9]\n\+\n$P[10]\n";
		$failCT++;
	}
} close IN; close PASS; close FAIL;


$pass_pct = sprintf("%.2f", ($passCT/$readCT)*100);
$fail_pct = sprintf("%.2f", ($failCT/$readCT)*100);
$time = localtime(time);
print LOG "... Completed at $time
Total Reads: $readCT
Passing Reads: $passCT ($pass_pct %)
Fail Reads: $failCT ($fail_pct %)
BA Exact Match: $ba_exact
BB Exact Match: $bb_exact
";
close LOG;

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

sub edit_distance {
    ($str1, $str2) = @_;

    # Create a matrix to store distances
    @distance;
    $distance[$_][0] = $_ for (0 .. length($str1));
    $distance[0][$_] = $_ for (0 .. length($str2));

    for $i (1 .. length($str1)) {
        for $j (1 .. length($str2)) {
            $cost = (substr($str1, $i-1, 1) eq substr($str2, $j-1, 1)) ? 0 : 1;
            
            $distance[$i][$j] = min(
                $distance[$i-1][$j] + 1,                  # Deletion
                $distance[$i][$j-1] + 1,                  # Insertion
                $distance[$i-1][$j-1] + $cost             # Substitution
            );
        }
    }

    return $distance[length($str1)][length($str2)];
}

sub min {
    ($x, $y, $z) = @_;
    return $x < $y ? ($x < $z ? $x : $z) : ($y < $z ? $y : $z);
}

}

1;