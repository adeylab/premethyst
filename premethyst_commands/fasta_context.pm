package premethyst_commands::fasta_context;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("fasta_context");

sub fasta_context {

getopts("O:C:G:z", \%opt);

$die = "

premethyst context-extract (options) -C [Req: Context] -G [genome fasta] -O [output prefix]
      or   context

Options:

   -O  [STR]  Output prefix. (req) Adds .bed (or .bed.gz)
   -C  [STR]  Context that will be searched. (req)
                  Can be multiple, semicolon separated.
                Format: [context]:[position of mC in motif, can be multiple, comma separated]
                  Can include A C G T or N, max length = 10
                  If mC position is a C it will be + strand, if G then -
   -G  [STR]  Fasta file that will be searched. (req)
   -z         Gzip output (gzip call: $gzip)
   
Example: Extracting CAG, where I want all sites for the C and the G reported,
         and matching both strands. Gzip the output:

premethyst context-extract -C CAG:1,3 -G hg38.fa.gz -O hg38.CAG -z

";

#### DIE
die "\nERROR: This tools is still in development. Use the comtext-specific scripts in the sciMETv2 repository.$die";

if (!defined $opt{'C'} || !defined $opt{'G'} || !defined $opt{'O'}) {die $die};
$opt{'O'} =~ s/\.gz$//; $opt{'O'} =~ s/\.bed$//;

@CONTEXTS = split(/;/, $opt{'C'});

foreach $context (@CONTEXTS) {
	($ctx,$posns) = split(/:/, $context);
	@P = split(/,/, $posns);
	$CONTEXT{'length'} = length($ctx);
	@{$CONTEXT{'mC'}} = @P;
	%{$CONTEXT{'N'}} = ();
	set_N_pos($ctx);
}

if ($opt{'G'} =~ /gz$/) {
	open IN, "$zcat $opt{'G'} |";
} else {
	open IN, "$opt{'G'}";
}

###


sub set_N_pos {
    ($sequence) = @_;
    $index = -1;
    while (($index = index($sequence, 'N', $index + 1)) != -1) {
		$CONTEXT{'N'}{$index + 1} = 1; # Adding 1 to convert from 0-based to 1-based position
    }
}

}

1;