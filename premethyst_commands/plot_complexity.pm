package premethyst_commands::plot_complexity;

use premethyst_commands::general;
use Getopt::Std; %opt = ();
use Exporter "import";
@EXPORT = ("plot_complexity");

sub plot_complexity {

@ARGV = @_;
getopts("O:A:a:R:f:p:k:h:w:T:m:y:t:D:H", \%opt);

#defaults
$alpha = 0.3;
$ptSize = 1;
$height = 6;
$hexBins = 300;
$width = 7;
$max = 7;
$title = "Library Complexity";
$contourCT = 50;
$die2 = "
premethyst plot-complexity [options] [complexity file(s) can be comma separated]

Options:
   -O   [STR]   Output prefix (default is complexity file 1 prefix)
   -H           Plot as hexplot (will not plot contours)
   -D   [INT]   Number of hexbins (def=$hexBins)
   -A   [STR]   Annotation file (to color code points)
   -a   [STR]   Comma separated list of annoations to include in plot
                 (requires -A to be specified)
   -p   [FLT]   Point size (def = $ptSize)
   -f   [FLT]   Alpha for plotting points (def = $alpha)
   -k   [STR]   If defined will color density lines the specified color (def = same as points)
                 either #hexcolor, or colorName
   -t   [INT]   Number of contours for 2d density (def = $contourCT)
   -w   [FLT]   Plot width (def = $width)
   -h   [FLT]   Plot height (def = $height)
   -y   [INT]   Max scale for plot in log10 unique reads (def = $max)
   -T   [STR]   Title (def = $title)
   -R   [STR]   Rscript call (def = $Rscript)

Note: Requires ggplot2 R package

";

if (!defined $ARGV[0]) {die $die2};
if (defined $opt{'a'} && !defined $opt{'A'}) {die "\nMust provide an annotaiton file (-A) if specifying annotations to plot (-a)!\n$die2"};
if (defined $opt{'R'}) {$Rscript = $opt{'R'}};
if (defined $opt{'A'}) {read_annot($opt{'A'})};
if (defined $opt{'a'}) {
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}
if (defined $opt{'m'}) {$kmean_centers = $opt{'m'}};
if (defined $opt{'p'}) {$ptSize = $opt{'p'}};
if (defined $opt{'f'}) {$alpha = $opt{'f'}};
if (defined $opt{'k'}) {
	if ($opt{'k'} =~ /^#/) {$cont_col = $opt{'k'}}
	else {$cont_col = "\"$opt{'k'}\""};
}
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'T'}) {$title = $opt{'T'}; $title =~ s/"//g};
if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
if (defined $opt{'y'}) {$max = $opt{'y'}};
if (defined $opt{'t'}) {$contourCT = $opt{'t'}};
if (defined $opt{'D'}) {$hexBins = $opt{'D'}; $opt{'H'} = 1};

if (!defined $opt{'O'}) {$opt{'O'} = $ARGV[0]};
$opt{'O'} =~ s/\.txt$//;

read_complexity($ARGV[0]);

open OUT, ">$opt{'O'}.plot.txt";
foreach $cellID (keys %CELLID_complexity) {
	$annot_cellID = $cellID;
	if (defined $opt{'a'}) {
		$annot = $CELLID_annot{$annot_cellID};
		if (defined $ANNOT_include{$annot} && defined $CELLID_annot{$annot_cellID}) {
			print OUT "$cellID\t$CELLID_annot{$annot_cellID}\t$CELLID_uniq_reads{$annot_cellID}\t$CELLID_complexity{$annot_cellID}\n";
		}
	} elsif (defined $opt{'A'} && defined $CELLID_annot{$cellID}) {
		print OUT "$cellID\t$CELLID_annot{$annot_cellID}\t$CELLID_uniq_reads{$annot_cellID}\t$CELLID_complexity{$annot_cellID}\n";
	} else {
		print OUT "$cellID\tCell\t$CELLID_uniq_reads{$annot_cellID}\t$CELLID_complexity{$annot_cellID}\n";
	}
} close OUT;

open R, ">$opt{'O'}.plot.r";
print R "
library(ggplot2)
IN<-read.table(\"$opt{'O'}.plot.txt\")
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
";
if (!defined $opt{'H'}) {
	print R "   geom_point(aes(V4,log10(V3),color=V2),size=$ptSize,alpha=$alpha,shape=16) +\n";
	if (defined $opt{'k'}) {
		print R "   geom_density2d(aes(V4,log10(V3),color=$cont_col,bins=$contourCT),size=0.3) +\n";
		} else {
		print R "   geom_density2d(aes(V4,log10(V3),color=V2,bins=$contourCT),size=0.3) +\n";
	}

} else {
	if (!defined $opt{'A'}) {
		print R "   geom_hex(aes(V4,log10(V3)),fill=\"lightsteelblue4\",bins=$hexBins) +\n";
	} else {
		print R "	geom_hex(aes(V4,log10(V3),fill=annot),bins=$hexBins) +
	guides(colour = guide_legend(override.aes = list(size=4))) +\n";
	}
}

print R "
   scale_x_continuous(limits=c(0,100)) +
   scale_y_continuous(limits=c(0,$max)) +
   xlab(\"Percent Passing Reads\") +
   ylab(\"log10 Passing Reads\") +
   ggtitle(\"$title\") +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=$width,height=$height)
";

print R "
PLT<-ggplot(data=subset(IN,V4<100&V4>0)) + theme_bw() +
	geom_histogram(aes(log10(V3),fill=V2)) +
	xlab(\"log10 Passing Reads\") +
	ylab(\"Counts\") +
	ggtitle(\"$title\") +
	scale_x_continuous(limits=c(0,$max)) +";
if (defined $opt{'A'}) {
print R "
	theme(legend.background=element_blank(),legend.title=element_blank())";
} else {
print R "
	theme(legend.position=\"none\")";
}
print R "
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.hist.pdf\",width=$width,height=$height)
";

close R;

system("$Rscript $opt{'O'}.plot.r");

sub read_complexity {
	%CELLID_uniq_reads = ();
	%CELLID_raw_reads = ();
	%CELLID_complexity = ();
	%CELLID_complexity_rank = ();
	@COMPLEXITY_FILES = split(/,/, $_[0]);
	foreach $complexity_file (@COMPLEXITY_FILES) {
		open COMPL, "$complexity_file";
		while ($comp_line = <COMPL>) {
			chomp $comp_line;
			($num,$cellID,$raw,$uniq,$pct) = split(/\t/, $comp_line);
			$CELLID_complexity_rank{$cellID} = $num;
			$CELLID_uniq_reads{$cellID} = $uniq;
			$CELLID_raw_reads{$cellID} = $raw;
			$CELLID_complexity{$cellID} = $pct;
		} close COMPL;
	}
}

}
1;
