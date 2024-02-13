package premethyst_commands::general;

use Getopt::Std; %opt = ();
use Exporter "import";

# Export all of the subroutines and variables that should be passed to scitools
@EXPORT = (
	"load_defaults",
		qw($ref_shortcuts),qw(@BASES),qw(%REF),qw(%VAR),qw($gzip),qw($zcat),
		qw($pigz),qw($samtools),qw($bedtools),qw($Rscript),qw($Pscript),
		qw($slack),qw($trim_galore),qw($bsbolt),qw($premethyst),
		qw($DEFAULTS_FILE)
);

# UTILITY COMMANDS

sub load_defaults {
	# Some global ones that are not configured
	$DEFAULTS_FILE = $_[0];
	$ref_shortcuts = "";
	@BASES = ("A", "C", "G", "T", "N");
	%REF; %VAR; %LOG_DIR;
	if (-e "$_[0]") {
		open DEF, "$_[0]";
		while ($def = <DEF>) {
			if ($def !~ /^#/) {
				chomp $def;
				$def =~ s/\r//g;
				($var,$val) = split(/=/, $def);
				if ($var =~ /_ref$/) { # reference shortcuts
					$refID = $var; $refID =~ s/_ref$//;
					$REF{$refID} = $val;
					$ref_shortcuts .= "\t$refID = $val\n";
				} else { # software calls & others
					if ($var eq "gzip") {$gzip = $val}
					elsif ($var eq "zcat") {$zcat = $val}
					elsif ($var eq "pigz") {$pigz = $val}
					elsif ($var eq "samtools") {$samtools = $val}
					elsif ($var eq "bedtools") {$bedtools = $val}
					elsif ($var eq "Rscript") {$Rscript = $val}
					elsif ($var eq "Pscript") {$Pscript = $val}
					elsif ($var eq "slack") {$slack = $val}
					elsif ($var eq "trim_galore") {$trim_galore = $val}
					elsif ($var eq "bsbolt") {$bsbolt = $val}
					elsif ($var eq "premethyst") {$premethyst = $val}
					else {$VAR{$var} = $val};
				}
			}
		} close DEF;
	}
	# Load vars that need to be specified for functionality if they are not found in the config file
	if (!defined $gzip) {$gzip = "gzip"};
	if (!defined $zcat) {$zcat = "zcat"};
	if (!defined $pigz) {$zcat = "pigz"};
	if (!defined $samtools) {$samtools = "samtools"};
	if (!defined $bedtools) {$bedtools = "bedtools"};
	if (!defined $Rscript) {$Rscript = "Rscript"};
	if (!defined $Pscript) {$Pscript = "python"};
	if (!defined $slack) {$slack = "slack"};
	if (!defined $trim_galore) {$trim_galore = "trim_galore"};
	if (!defined $bsbolt) {$bsbolt = "bsbolt"};
	if (!defined $premethyst) {$premethyst = "premethyst"};
}

1;