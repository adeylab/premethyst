package premethyst_commands::calls_rename;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("calls_rename");

sub calls_rename {
	
getopts("P:S:", \%opt);

$die = "

premethyst calls-rename (options) [folder with cellcall files(req)]
      or   rename

Simply prepends or appends what is provided to the .cov files.
Automatically includes a '-' separator. Do not include a '-' in
any of the added fields, since this is used to split out the
original barcode/cellID for any filtering.

Cov files MUST be in the standard format of:
   [cellID].[context].cov
e.g.
   ACGTACGT.CG.cov

Options (at least 1 required):

   -P  [STR]  Prefix to add to cell IDs.
                Should not start with a number or special character.
   -S  [STR]  Suffix to add to cell IDs.

";

if (!defined $ARGV[0] || (!defined $opt{'P'} && !defined $opt{'S'})) {die $die};

opendir(CCD, "$ARGV[0]");
@CELLCALLS = readdir(CCD);
closedir(CCD);

foreach $cellCallFile (@CELLCALLS) {
	if ($cellCallFile =~ /\.cov/) {
		($cellID,$context,$cov) = split(/\./, $cellCallFile);
		if (defined $opt{'P'}) {
			$newID = $opt{'P'}."-".$cellID;
		} else {
			$newID = $cellID;
		}
		if (defined $opt{'S'}) {
			$newID .= "-".$opt{'S'};
		}
		system("mv $ARGV[0]/$cellCallFile $ARGV[0]/$newID.$context.cov");
	}
}

}

1;