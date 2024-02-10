## WARNING: Premethyst is still in development.
We recommend using sciMETv2 scripts until development is completed.
https://github.com/adeylab/sciMETv2

# Premethyst
Premethyst is a pre-processing pipeline for single-cell DNA methylaiton data generated using sciMET workflows. It starts with demultiplexed fastq files and proceeds through alignment, methylation call extraction, and produces a compressed h5 file containing all base-level methylaiton calls for each cell with CG and CH contexts as defaults. The h5 file can the be directly analyzed using amethyst.

Alternatively, if data are generated using the ScaleBio commercial kit, the associated analysis workflow produces bam files that can be directly processed from the extract step on. This requires less prerequisites due to the container-based nature of the comemrcially developed workflow.

## External Software Requirements
Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

BSBolt (https://github.com/NuttyLogic/BSBolt)

Samtools (https://www.htslib.org/)

Recommended: unidex (https://github.com/adeylab/unidex) for demultiplexing

## Setup
Download the repository and add the premethyst executable to your path. Premethyst has a number of default parameters that you can set using the configuration file which should be in the same directory as the premethyst executable and can be edited as needed.
## Processing Workflow
Starting files: Pair of fastq files where the read name contains the error-corrected cell barcode. Ideally generated using unidex (https://github.com/adeylab/unidex)
```
@[CELL_BARCODE]:NNN#0/1
```
### Read Trimming
Paired fastq files are trimmed, with defaults set for sciMET-generated libraries. Leverages TrimGalore.
```
premethyst fastq-trim (options) -O OutputPrefix -1 read1.fq.gz -2 read2.fq.gz
```
### Read Alignment
Trimmed reads are aligned using BSBolt to a reference genome and sorted by read name.
```
premethyst fastq-align (options) -R [reference path] -O [output prefix] -1 [read1.trimmed.fq.gz] -2 [read2.trimmed.fq.gz]
```
### Read Deduplication
Aligned reads are quality filtered and PCR duplicates removed. Also produces a complexity file for QC purposes.
```
premethyst bam-rmdup (options) [namesort bam file]
```
### Plot Complexity
Plot the unique, passing reads by percent unique, passing reads per cell to assess sequencing saturation and read count cutoff for cells. Can also include an annotation file (tab-delimited: cellID (tab) Annotation) to color by metadata.
```
premethyst plot-complexity [options] [complexity file(s) can be comma separated]
```
### Extract Methylation Calls
Extract mC calls to generate a folder of cellCall files for each cell.
```
premethyst bam-extract (options) [rmdup & filtered bam file, name-sorted]
```
### Further Filter and Adjust CellCall Folder
Append additional naming criteria, extract specific contexts to be included. Remove alternative contexts. Further filer to remove certain cells. Additional filtering can be carried out later in amethyst.
```
# Add on a name for cell IDs for future multiplexing
premethyst rename-calls -P (name_prefix) -S (name_suffix) [cellCalls Folder]

# Filter to specific criteria
premethyst filter-calls (options) [cellInfo file] [input cellCall folder] [filtered cellCall folder]

# Extract an additional context, e.g. GCH
premethyst context-extract -O MyOut -C CH -N GCH -B myGCsites.bed -X HCH MyIn
```
### Convert to a Compressed h5 File
Produce a compressed h5 file that can be loaded into amethyst for subsequent analysis, dataset interaction and visualization in a personal compuer environment (ie RStudio)
```
premethyst calls2h5 [input calls folder] [output h5 file prefix]
```

