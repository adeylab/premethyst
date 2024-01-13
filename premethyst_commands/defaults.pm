package premethyst_commands::defaults;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use FindBin '$RealBin';
use Exporter "import";
@EXPORT = ("defaults");

sub defaults {

$ts = localtime(time);
print "
premethyst defaults
Run at $ts
premethyst location: $RealBin
Config file: $ARGV[0]

Default command executable calls:
   premethyst:  $premethyst
   gzip:        $gzip
   zcat:        $zcat
   samtools:    $samtools
   bedtools:    $bedtools
   R scripts:   $Rscript
   Py scripts:  $Pscript
   slack:       $slack
   trim_galore: $trim_galore
   bsbolt:      $bsbolt

Reference genomes:\n";
foreach $refID (sort keys %REF) {
	print "   $refID = $REF{$refID}\n";
}
print "
Other defaults:\n";
foreach $var (sort keys %VAR) {
	print "   $var = $VAR{$var}\n";
}
print "\n";

}
1;