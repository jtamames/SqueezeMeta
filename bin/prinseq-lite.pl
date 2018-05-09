#!/usr/bin/perl

#===============================================================================
#   Author: Robert SCHMIEDER, Computational Science Research Center @ SDSU, CA
#
#   File: prinseq-lite
#   Date: 2013-11-10
#   Version: 0.20.4 lite
#
#   Usage:
#      prinseq-lite [options]
#
#      Try 'prinseq-lite -h' for more information.
#
#    Purpose: PRINSEQ will help you to preprocess your genomic or metagenomic
#             sequence data in FASTA or FASTQ format. The lite version does not
#             require any non-core perl modules for processing.
#
#    Bugs: Please use http://sourceforge.net/p/prinseq/bugs/ or
#          https://groups.google.com/d/forum/edwardslabtools
#
#===============================================================================

use strict;
use warnings;

#use Data::Dumper; ###
use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile); #for output files
use Fcntl qw(:flock SEEK_END); #for log file
use Digest::MD5 qw(md5_hex); #for dereplication
use Cwd;
use List::Util qw(sum min max);

$| = 1; # Do not buffer output

my $WINDOWSIZE = 64;
my $WINDOWSTEP = 32;
my $WORDSIZE = 3;
my @WINDOWSIZEARRAY = (0..61);
my $LOG62 = log(62);
my $ONEOVERLOG62 = 1/log(62);
my $POINTFIVE = 1/2;
my $LINE_WIDTH = 60;
my $GRAPH_DATA_SEQ_MAX_LENGTH = 1000;
my $TRIM_QUAL_WINDOW = 1;
my $TRIM_QUAL_STEP = 1;
my $TRIM_QUAL_TYPE = 'min';
my $TRIM_QUAL_RULE = 'lt';
my $TAG_LENGTH = 20;
my %MIDS = (ACGAGTGCGT => 0,
            ACGCTCGACA => 0,
            AGACGCACTC => 0,
            AGCACTGTAG => 0,
            ATCAGACACG => 0,
            ATATCGCGAG => 0,
            CGTGTCTCTA => 0,
            CTCGCGTGTC => 0,
            TAGTATCAGC => 0,
            TCTCTATGCG => 0,
            TGATACGTCT => 0,
            TACTGAGCTA => 0,
            CATAGTAGTG => 0,
            CGAGAGATAC => 0,
            ACACGACGACT => 0,
            ACACGTAGTAT => 0,
            ACACTACTCGT => 0,
            ACGACACGTAT => 0,
            ACGAGTAGACT => 0,
            ACGCGTCTAGT => 0,
            ACGTACACACT => 0,
            ACGTACTGTGT => 0,
            ACGTAGATCGT => 0,
            ACTACGTCTCT => 0,
            ACTATACGAGT => 0,
            ACTCGCGTCGT => 0);
my $MIDCHECKLENGTH = 15; #maximum MID length plus possible key length (by default 4 bp for 454)
my %DN_DI = ('AA' => 0, 'AC' => 0, 'AG' => 0, 'AT' => 0, 'CA' => 0, 'CC' => 0, 'CG' => 0, 'CT' => 0, 'GA' => 0, 'GC' => 0, 'GG' => 0, 'GT' => 0, 'TA' => 0, 'TC' => 0, 'TG' => 0, 'TT' => 0);
my %GRAPH_OPTIONS = map {$_ => 1} qw(ld gc qd ns pt ts aq de da sc dn);
my $VERSION = '0.20.4';
my $WHAT = 'lite';

my $man = 0;
my $help = 0;
my %params = ('help' => \$help, 'h' => \$help, 'man' => \$man);
GetOptions( \%params,
            'help|h',
            'man',
            'verbose',
            'version' => sub { print "PRINSEQ-$WHAT $VERSION\n"; exit; },
            'fastq=s',
            'fasta=s',
            'fastq2=s',
            'fasta2=s',
            'qual=s',
            'min_len=i',
            'max_len=i',
            'range_len=s',
            'min_gc=i',
            'max_gc=i',
            'range_gc=s',
            'min_qual_score=i',
            'max_qual_score=i',
            'min_qual_mean=i',
            'max_qual_mean=i',
            'ns_max_p=i',
            'ns_max_n=i',
            'noniupac',
            'seq_num=i',
            'derep=i',
            'derep_min=i',
            'lc_method=s',
            'lc_threshold=i',
            'trim_to_len=i',
            'trim_left=i',
            'trim_right=i',
            'trim_left_p=i',
            'trim_right_p=i',
            'trim_tail_left=i',
            'trim_tail_right=i',
            'trim_ns_left=i',
            'trim_ns_right=i',
            'trim_qual_left=i',
            'trim_qual_right=i',
            'trim_qual_type=s',
            'trim_qual_rule=s',
            'trim_qual_window=i',
            'trim_qual_step=i',
            'seq_case=s',
            'dna_rna=s',
            'line_width=i',
            'rm_header',
            'seq_id=s',
            'seq_id_mappings:s',
            'out_format=i',
            'out_good=s',
            'out_bad=s',
            'stats_len',
            'stats_dinuc',
            'stats_info',
            'stats_tag',
            'stats_dupl',
            'stats_ns',
            'stats_assembly',
            'stats_all',
            'aa',
            'log:s',
            'graph_data:s',
            'graph_stats=s',
            'phred64',
            'qual_noscale',
            'no_qual_header',
            'exact_only',
            'web:s',
            'filename1=s',
            'filename2=s',
            'custom_params=s',
            'params=s'
            ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

=head1 NAME

PRINSEQ - PReprocessing and INformation of SEQuence data

=head1 VERSION

PRINSEQ-lite 0.20.4

=head1 SYNOPSIS

perl prinseq-lite.pl [-h] [-help] [-version] [-man] [-verbose] [-fastq input_fastq_file] [-fasta input_fasta_file] [-fastq2 input_fastq_file_2] [-fasta2 input_fasta_file_2] [-qual input_quality_file] [-min_len int_value] [-max_len int_value] [-range_len ranges] [-min_gc int_value] [-max_gc int_value] [-range_gc ranges] [-min_qual_score int_value] [-max_qual_score int_value] [-min_qual_mean int_value] [-max_qual_mean int_value] [-ns_max_p int_value] [-ns_max_n int_value] [-noniupac] [-seq_num int_value] [-derep int_value] [-derep_min int_value] [-lc_method method_name] [-lc_threshold int_value] [-trim_to_len int_value] [-trim_left int_value] [-trim_right int_value] [-trim_left_p int_value] [-trim_right_p int_value] [-trim_ns_left int_value] [-trim_ns_right int_value] [-trim_tail_left int_value] [-trim_tail_right int_value] [-trim_qual_left int_value] [-trim_qual_right int_value] [-trim_qual_type type] [-trim_qual_rule rule] [-trim_qual_window int_value] [-trim_qual_step int_value] [-seq_case case] [-dna_rna type] [-line_width int_value] [-rm_header] [-seq_id id_string] [-out_format int_value] [-out_good filename_prefix] [-out_bad filename_prefix] [-phred64] [-stats_info] [-stats_len] [-stats_dinuc] [-stats_tag] [-stats_dupl] [-stats_ns] [-stats_assembly] [-stats_all] [-aa] [-graph_data file] [-graph_stats string] [-qual_noscale] [-no_qual_header] [-exact_only] [-log file] [-custom_params string] [-params file] [-seq_id_mappings file]

=head1 DESCRIPTION

PRINSEQ will help you to preprocess your genomic or metagenomic sequence data in FASTA (and QUAL) or FASTQ format. The lite version does not require any non-core perl modules for processing.

=head1 OPTIONS

=over 8

=item B<-help> | B<-h>

Print the help message; ignore other arguments.

=item B<-man>

Print the full documentation; ignore other arguments.

=item B<-version>

Print program version; ignore other arguments.

=item B<-verbose>

Prints status and info messages during processing.

=item B<***** INPUT OPTIONS *****>

=item B<-fastq> <file>

Input file in FASTQ format that contains the sequence and quality data. Use stdin instead of a file name to read from STDIN (-fasta stdin). This can be useful to process compressed files using Unix pipes.

=item B<-fasta> <file>

Input file in FASTA format that contains the sequence data. Use stdin instead of a file name to read from STDIN (-fastq stdin). This can be useful to process compressed files using Unix pipes.

=item B<-qual> <file>

Input file in QUAL format that contains the quality data.

=item B<-fastq2> <file>

For paired-end data only. Input file in FASTQ format that contains the sequence and quality data. The sequence identifiers for two matching paired-end sequences in separate files can be marked by /1 and /2, or _L and _R, or _left and _right, or must have the exact same identifier in both input files. The input sequences must be sorted by their sequence identifiers. Singletons are allowed in the input files.

=item B<-fasta2> <file>

For paired-end data only. Input file in FASTA format that contains the sequence data. The sequence identifiers for two matching paired-end sequences in separate files can be marked by /1 and /2, or _L and _R, or _left and _right, or must have the exact same identifier in both input files. The input sequences must be sorted by their sequence identifiers. Singletons are allowed in the input files.

=item B<-params> <file>

Input file in text format that contains PRINSEQ parameters. Each parameter should be specified on a new line and arguments should be separated by spaces or tabs. Comments can be specified on lines starting with the # sign. Can be combined with command line parameters. Parameters specified on the command line will overwrite the arguments in the file (if any).

=item B<-si13>

This option was replaced by option -phred64.

=item B<-phred64>

Quality data in FASTQ file is in Phred+64 format (http://en.wikipedia.org/wiki/FASTQ_format#Encoding). Not required for Illumina 1.8+, Sanger, Roche/454, Ion Torrent, PacBio data.

=item B<-aa>

Input is amino acid (protein) sequences instead of nucleic acid (DNA or RNA) sequences. Allowed amino acid characters: ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*- and allowed nucleic acid characters: ACGTURYKMSWBDHVNXacgturykmswbdhvnx-

The following options are ignored for -aa: stats_dinuc,stats_tag,stats_ns,dna_rna

=item B<***** OUTPUT OPTIONS *****>

=item B<-out_format> <integer>

To change the output format, use one of the following options. If not defined, the output format will be the same as the input format.

1 (FASTA only), 2 (FASTA and QUAL), 3 (FASTQ), 4 (FASTQ and FASTA), or 5 (FASTQ, FASTA and QUAL)

=item B<-out_good> <string>

By default, the output files are created in the same directory as the input file containing the sequence data with an additional "_prinseq_good_XXXX" in their name (where XXXX is replaced by random characters to prevent overwriting previous files). To change the output filename and location, specify the filename using this option. The file extension will be added automatically (either .fasta, .qual, or .fastq). For paired-end data, filenames contain additionally "_1", "_1_singletons", "_2", and "_2_singletons" before the file extension. Use "-out_good null" to prevent the program from generating the output file(s) for data passing all filters. Use "-out_good stdout" to write data passing all filters to STDOUT (only for FASTA or FASTQ output files).

Example: use "file_passed" to generate the output file file_passed.fasta in the current directory

=item B<-out_bad> <string>

By default, the output files are created in the same directory as the input file containing the sequence data with an additional "_prinseq_bad_XXXX" in their name (where XXXX is replaced by random characters to prevent overwriting previous files). To change the output filename and location, specify the filename using this option. The file extension will be added automatically (either .fasta, .qual, or .fastq). For paired-end data, filenames contain additionally "_1" and "_2" before the file extension. Use "-out_bad null" to prevent the program from generating the output file(s) for data not passing any filter. Use "-out_bad stdout" to write data not passing any filter to STDOUT (only for FASTA or FASTQ output files).

Example: use "file_filtered" to generate the output file file_filtered.fasta in the current directory

Example: "-out_good stdout -out_bad null" will write data passing filters to STDOUT and data not passing any filter will be ignored

=item B<-log> <file>

Log file to keep track of parameters, errors, etc. The log file name is optional. If no file name is given, the log file name will be "inputname.log". If the log file already exists, new content will be added to the file.

=item B<-graph_data> <file>

File that contains the necessary information to generate the graphs similar to the ones in the web version. The file name is optional. If no file name is given, the file name will be "inputname.gd". If the file already exists, new content will overwrite the file. Use "-out_good null -out_bad null" to prevent generating any additional outputs. (See below for more options related to the graph data.)

The graph data can be used as input for the prinseq-graphs.pl file to generate the PNG graph files or an HTML report file. If you have trouble installing the required prinseq-graphs.pl modules or want to see an output
example report, upload the graph data file at: http://edwards.sdsu.edu/prinseq/ -> Choose "Get Report"

=item B<-graph_stats> <string>

Use this option to select what statistics should be calculated and included in the graph_data file. This is useful if you e.g. do not need sequence complexity information, which requires a lot of computation. Requires to have graph_data specified. Default is all selected.

Allowed option are (separate multiple by comma with no spaces): ld (Length distribution), gc (GC content distribution), qd (Base quality distribution), ns (Occurence of N), pt (Poly-A/T tails), ts (Tag sequence check), aq (Assembly quality measure), de (Sequence duplication - exact only), da (Sequence duplication - exact + 5'/3'), sc (Sequence complexity), dn (Dinucleotide odds ratios, includes the PCA plots)

Example use: -graph_stats ld,gc,qd,de

=item B<-qual_noscale>

Use this option if all your sequences are shorter than 100bp as they do not require to scale quality data to 100 data points in the graph. By default, quality scores of sequences shorter than 100bp or longer than 100bp are fit to 100 data points. (To retrieve this information and calculate the graph data would otherwise require to parse the data two times or store all the quality data in memory.)

=item B<-no_qual_header>

In order to reduce the file size, this option will generate an empty header line for the quality data in FASTQ files. Instead of +header, only the + sign will be output. The header of the sequence data will be left unchanged. This option applies to FASTQ output files only.

=item B<-exact_only>

Use this option to check for exact (forward and reverse) duplicates only when generating the graph data. This allows to keep the memory requirements low for large input files and is faster. This option will automatically be applied when using -derep options 1 and/or 4 only. Specify option -derep 1 or -derep 4 if you do not want to apply both at the same time.

=item B<-seq_id_mappings> <file>

Text file containing the old and new (specified with -seq_id) identifiers for later reference. This option is useful if e.g. a renamed sequence has to be identified based on the new sequence identifier. The file name is optional. If no file name is given, the file name will be "inputname_prinseq_good.ids" (only good sequences are renamed). If a file with the same name already exists, new content will overwrite the old file. The text file contains one sequence identifier pair per line, separated by tabs (old-tab-new). Requires option -seq_id.


=item B<***** FILTER OPTIONS *****>

=item B<-min_len> <integer>

Filter sequence shorter than min_len.

=item B<-max_len> <integer>

Filter sequence longer than max_len.

=item B<-range_len> <string>

Filter sequence by length range. Multiple range values should be separated by comma without spaces.

Example: -range_len 50-100,250-300

=item B<-min_gc> <integer>

Filter sequence with GC content below min_gc.

=item B<-max_gc> <integer>

Filter sequence with GC content above max_gc.

=item B<-range_gc> <string>

Filter sequence by GC content range. Multiple range values should be separated by comma without spaces.

Example: -range_gc 50-60,75-90

=item B<-min_qual_score> <integer>

Filter sequence with at least one quality score below min_qual_score.

=item B<-max_qual_score> <integer>

Filter sequence with at least one quality score above max_qual_score.

=item B<-min_qual_mean> <integer>

Filter sequence with quality score mean below min_qual_mean.

=item B<-max_qual_mean> <integer>

Filter sequence with quality score mean above max_qual_mean.

=item B<-ns_max_p> <integer>

Filter sequence with more than ns_max_p percentage of Ns.

=item B<-ns_max_n> <integer>

Filter sequence with more than ns_max_n Ns.

=item B<-noniupac>

Filter sequence with characters other than A, C, G, T or N.

=item B<-seq_num> <integer>

Only keep the first seq_num number of sequences (that pass all other filters).

=item B<-derep> <integer>

Type of duplicates to filter. Allowed values are 1, 2, 3, 4 and 5. Use integers for multiple selections (e.g. 124 to use type 1, 2 and 4). The order does not matter. Option 2 and 3 will set 1 and option 5 will set 4 as these are subsets of the other option.

1 (exact duplicate), 2 (5' duplicate), 3 (3' duplicate), 4 (reverse complement exact duplicate), 5 (reverse complement 5'/3' duplicate)

=item B<-derep_min> <integer>

This option specifies the number of allowed duplicates. If you want to remove sequence duplicates that occur more than x times, then you would specify x+1 as the -derep_min values. For examples, to remove sequences that occur more than 5 times, you would specify -derep_min 6. This option can only be used in combination with -derep 1 and/or 4 (forward and/or reverse exact duplicates). [default : 2]

=item B<-lc_method> <string>

Method to filter low complexity sequences. The current options are "dust" and "entropy". Use "-lc_method dust" to calculate the complexity using the dust method.

=item B<-lc_threshold> <integer>

The threshold value (between 0 and 100) used to filter sequences by sequence complexity. The dust method uses this as maximum allowed score and the entropy method as minimum allowed value.

=item B<-custom_params> <string>

Can be used to specify additional filters. The current set of possible rules is limited and has to follow the specifications below. The custom parameters have to be specified within quotes (either ' or ").

Please separate parameter values with a space and separate new parameter sets with semicolon (;). Parameters are defined by two values:
  (1) the pattern (any combination of the letters "ACGTN"),
  (2) the number of repeats or percentage of occurence
Percentage values are defined by a number followed by the %-sign (without space).
If no %-sign is given, it is assumed that the given number specifies the number of repeats of the pattern.

Examples: "AAT 10" (filters out sequences containing AATAATAATAATAATAATAATAATAATAAT anywhere in the sequence), "T 70%" (filters out sequences with more than 70% Ts in the sequence), "A 15" (filters out sequences containing AAAAAAAAAAAAAAA anywhere in the sequence), "AAT 10;T 70%;A 15" (apply all three filters)

=item B<***** TRIM OPTIONS *****>

=item B<-trim_to_len> <integer>

Trim all sequence from the 3'-end to result in sequence with this length.

=item B<-trim_left> <integer>

Trim sequence at the 5'-end by trim_left positions.

=item B<-trim_right> <integer>

Trim sequence at the 3'-end by trim_right positions.

=item B<-trim_left_p> <integer>

Trim sequence at the 5'-end by trim_left_p percentage of read length. The trim length is rounded towards the lower integer (e.g. 143.6 is rounded to 143 positions). Use an integer between 1 and 100 for the percentage value.

=item B<-trim_right_p> <integer>

Trim sequence at the 3'-end by trim_right_p percentage of read length. The trim length is rounded towards the lower integer (e.g. 143.6 is rounded to 143 positions). Use an integer between 1 and 100 for the percentage value.

=item B<-trim_tail_left> <integer>

Trim poly-A/T tail with a minimum length of trim_tail_left at the 5'-end.

=item B<-trim_tail_right> <integer>

Trim poly-A/T tail with a minimum length of trim_tail_right at the 3'-end.

=item B<-trim_ns_left> <integer>

Trim poly-N tail with a minimum length of trim_ns_left at the 5'-end.

=item B<-trim_ns_right> <integer>

Trim poly-N tail with a minimum length of trim_ns_right at the 3'-end.

=item B<-trim_qual_left> <integer>

Trim sequence by quality score from the 5'-end with this threshold score.

=item B<-trim_qual_right> <integer>

Trim sequence by quality score from the 3'-end with this threshold score.

=item B<-trim_qual_type> <string>

Type of quality score calculation to use. Allowed options are min, mean, max and sum. [default: min]

=item B<-trim_qual_rule> <string>

Rule to use to compare quality score to calculated value. Allowed options are lt (less than), gt (greater than) and et (equal to). [default: lt]

=item B<-trim_qual_window> <integer>

The sliding window size used to calculate quality score by type. To stop at the first base that fails the rule defined, use a window size of 1. [default: 1]

=item B<-trim_qual_step> <integer>

Step size used to move the sliding window. To move the window over all quality scores without missing any, the step size should be less or equal to the window size. [default: 1]

=item B<***** REFORMAT OPTIONS *****>

=item B<-seq_case> <string>

Changes sequence character case to upper or lower case. Allowed options are "upper" and "lower". Use this option to remove soft-masking from your sequences.

=item B<-dna_rna> <string>

Convert sequence between DNA and RNA. Allowed options are "dna" (convert from RNA to DNA) and "rna" (convert from DNA to RNA).

=item B<-line_width> <integer>

Sequence characters per line. Use 0 if you want each sequence in a single line. Use 80 for line breaks every 80 characters. Note that this option only applies to FASTA output files, since FASTQ files store sequences without additional line breaks. [default: 60]

=item B<-rm_header>

Remove the sequence header. This includes everything after the sequence identifier (which is kept unchanged).

=item B<-seq_id> <string>

Rename the sequence identifier. A counter is added to each identifier to assure its uniqueness. Use option -seq_id_mappings to generate a file containing the old and new identifiers for later reference.

Example: "mySeq_10" will generate the IDs (in FASTA format) >mySeq_101, >mySeq_102, >mySeq_103, ...

=item B<***** SUMMARY STATISTIC OPTIONS *****>

The summary statistic values are written to STDOUT in the form: "parameter_name statistic_name value" (without the quotes). For example, "stats_info reads 10000" or "stats_len max 500". Only one statistic is written per line and values are separated by tabs.

If you specify any statistic option, no other ouput will be generated. To preprocess data, do not specify a statistics option.

=item B<-stats_info>

Outputs basic information such as number of reads (reads) and total bases (bases).

=item B<-stats_len>

Outputs minimum (min), maximum (max), range (range), mean (mean), standard deviation (stddev), mode (mode) and mode value (modeval), and median (median) for read length.

=item B<-stats_dinuc>

Outputs the dinucleotide odds ratio for AA/TT (aatt), AC/GT (acgt), AG/CT (agct), AT (at), CA/TG (catg), CC/GG (ccgg), CG (cg), GA/TC (gatc), GC (gc) and TA (ta).

=item B<-stats_tag>

Outputs the probability of a tag sequence at the 5'-end (prob5) and 3'-end (prob3) in percentage (0..100). Provides the number of predefined MIDs (midnum) and the MID sequences (midseq, separated by comma, only provided if midnum > 0) that occur in more than 34/100 (approx. 3%) of the reads.

=item B<-stats_dupl>

Outputs the number of exact duplicates (exact), 5' duplicates (5), 3' duplicates (3), exact duplicates with reverse complements (exactrevcom) and 5'/3' duplicates with reverse complements (revcomp), and total number of duplicates (total). The maximum number of duplicates is given under the value name with an additional "maxd" (e.g. exactmaxd or 5maxd).

=item B<-stats_ns>

Outputs the number of reads with ambiguous base N (seqswithn), the maximum number of Ns per read (maxn) and the maximum percentage of Ns per read (maxp). The maxn and maxp value are not necessary from the same sequence.

=item B<-stats_assembly>

Outputs the N50, N90, etc contig sizes. The Nxx contig size is a weighted median that is defined as the length of the smallest contig C in the sorted list of all contigs where the cumulative length from the largest contig to contig C is at least xx% of the total length (sum of contig lengths).

=item B<-stats_all>

Outputs all available summary statistics.

=item B<***** ORDER OF PROCESSING *****>

The available options are processed in the following order:

seq_num, trim_left, trim_right, trim_left_p, trim_right_p, trim_qual_left, trim_qual_right, trim_tail_left, trim_tail_right, trim_ns_left, trim_ns_right, trim_to_len, min_len, max_len, range_len, min_qual_score, max_qual_score, min_qual_mean, max_qual_mean, min_gc, max_gc, range_gc, ns_max_p, ns_max_n, noniupac, lc_method, derep, seq_id, seq_case, dna_rna, out_format

=back

=head1 AUTHOR

Robert SCHMIEDER, C<< <rschmieder_at_gmail_dot_com> >>

=head1 BUGS

If you find a bug please email me at C<< <rschmieder_at_gmail_dot_com> >> or use http://sourceforge.net/tracker/?group_id=315449 so that I can make PRINSEQ better.

=head1 COPYRIGHT

Copyright (C) 2010-2012  Robert SCHMIEDER

=head1 LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

#
################################################################################
## DATA AND PARAMETER CHECKING
################################################################################
#

my ($file1,$file2,$command,@dataread,$aa,%filtercount,$slashnum,$trimnum1,$trimnum2);

#check if params file
if(exists $params{params}) {
    my $ps = &readParamsFile($params{params});
    foreach my $p (keys %$ps) {
        next if(exists $params{$p});
        $params{$p} = $ps->{$p};
    }
}

#check if amino acid or nucleic acid input
if(exists $params{aa}) {
    $aa = 1;
    $command .= ' -aa';
} else {
    $aa = 0;
}

#Check if input file exists and check if file format is correct
if(exists $params{fasta} && exists $params{fastq}) {
    &printError('fasta and fastq cannot be used together');
} elsif(exists $params{fasta}) {
    $command .= ' -fasta '.$params{fasta};
    $file1 = $params{fasta};
    if($params{fasta} eq 'stdin') {
        if(exists $params{qual} && $params{qual} eq 'stdin') {
            &printError('input from STDIN is only allowed for either of the input files');
        } else {
            my $format = &checkInputFormat();
            unless($format eq 'fasta') {
                &printError('input data for -fasta is in '.uc($format).' format not in FASTA format');
            }
        }
    } elsif(-e $params{fasta}) {
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'fasta') {
            &printError('input file for -fasta is in '.uc($format).' format not in FASTA format');
        }
    } else {
        &printError("could not find input file \"".$params{fasta}."\"");
    }
} elsif(exists $params{fastq}) {
    $command .= ' -fastq '.$params{fastq};
    $file1 = $params{fastq};
    if($params{fastq} eq 'stdin') {
        my $format = &checkInputFormat();
        unless($format eq 'fastq') {
            &printError('input data for -fastq is in '.uc($format).' format not in FASTQ format');
        }
    } elsif(-e $params{fastq}) {
        #check for file format
        my $format = &checkFileFormat($file1);
        unless($format eq 'fastq') {
            &printError('input file for -fastq is in '.uc($format).' format not in FASTQ format');
        }
    } else {
        &printError("could not find input file \"".$params{fastq}."\"");
    }
} else {
    &printError("you did not specify an input file containing the query sequences");
}
if(exists $params{fastq} && exists $params{qual}) {
    &printError('fastq and qual cannot be used together');
} elsif(exists $params{qual}) {
    $command .= ' -qual '.$params{qual};
    if($params{qual} eq 'stdin') {
        &printError('QUAL data cannot be read from STDIN');
    } elsif(-e $params{qual}) {
        #check for file format
        my $format = &checkFileFormat($params{qual});
        unless($format eq 'qual') {
            &printError('input file for -qual is in '.uc($format).' format not in QUAL format');
        }
    } else {
        &printError("could not find input file \"".$params{qual}."\"");
    }
}
if(exists $params{fasta2} && exists $params{fastq2}) {
    &printError('fasta2 and fastq2 cannot be used together');
} elsif(exists $params{fasta2}) {
    if(!exists $params{fasta}) {
        &printError('option fasta2 requires option fasta');
    } elsif($params{fasta} eq $params{fasta2}) {
        &printError('option fasta and fasta2 cannot be the same input file');
    } else {
        $command .= ' -fasta2 '.$params{fasta2};
        $file2 = $params{fasta2};
        if($params{fasta} eq 'stdin' || $params{fasta2} eq 'stdin') {
            &printError('paired-end data cannot be processed from STDIN');
        } elsif(-e $params{fasta2}) {
            #check for file format
            my $format = &checkFileFormat($file2);
            unless($format eq 'fasta') {
                &printError('input file for -fasta2 is in '.uc($format).' format not in FASTA format');
            }
        } else {
            &printError("could not find input file \"".$params{fasta2}."\"");
        }
    }
    ($slashnum,$trimnum1,$trimnum2) = &checkSlashnum($file2);
} elsif(exists $params{fastq2}) {
    if(!exists $params{fastq}) {
        &printError('option fastq2 requires option fastq');
    } elsif($params{fastq} eq $params{fastq2}) {
        &printError('option fastq and fastq2 cannot be the same input file');
    } else {
        $command .= ' -fastq2 '.$params{fastq2};
        $file2 = $params{fastq2};
        if($params{fastq} eq 'stdin' || $params{fastq2} eq 'stdin') {
            &printError('paired-end data cannot be processed from STDIN');
        } elsif(-e $params{fastq2}) {
            #check for file format
            my $format = &checkFileFormat($file2);
            unless($format eq 'fastq') {
                &printError('input file for -fastq2 is in '.uc($format).' format not in FASTQ format');
            }
        } else {
            &printError("could not find input file \"".$params{fastq2}."\"");
        }
    }
    ($slashnum,$trimnum1,$trimnum2) = &checkSlashnum($file2);
}


#check if stats_all
if(exists $params{stats_all}) {
    $params{stats_info} = 1;
    $params{stats_len} = 1;
    $params{stats_dupl} = 1 unless($file2);
    $params{stats_dinuc} = 1;
    $params{stats_tag} = 1;
    $params{stats_ns} = 1;
    $params{stats_assembly} = 1;
    delete($params{stats_all});
}
if($aa) {
    delete($params{stats_dinuc});
    delete($params{stats_tag});
    delete($params{stats_ns});
}
if($file2) {
    delete($params{stats_dupl});
    delete($params{stats_assembly});
}

#check if anything todo
unless( exists $params{min_len} ||
        exists $params{max_len} ||
        exists $params{range_len} ||
        exists $params{min_gc} ||
        exists $params{max_gc} ||
        exists $params{range_gc} ||
        exists $params{min_qual_score} ||
        exists $params{max_qual_score} ||
        exists $params{min_qual_mean} ||
        exists $params{max_qual_mean} ||
        exists $params{ns_max_p} ||
        exists $params{ns_max_n} ||
        exists $params{noniupac} ||
        exists $params{seq_num} ||
        exists $params{derep} ||
        exists $params{lc_method} ||
        exists $params{lc_threshold} ||
        exists $params{trim_to_len} ||
        exists $params{trim_left} ||
        exists $params{trim_right} ||
        exists $params{trim_left_p} ||
        exists $params{trim_right_p} ||
        exists $params{trim_tail_left} ||
        exists $params{trim_tail_right} ||
        exists $params{trim_ns_left} ||
        exists $params{trim_ns_right} ||
        exists $params{trim_qual_left} ||
        exists $params{trim_qual_right} ||
        exists $params{trim_qual_type} ||
        exists $params{trim_qual_rule} ||
        exists $params{trim_qual_window} ||
        exists $params{trim_qual_step} ||
        exists $params{seq_case} ||
        exists $params{dna_rna} ||
        exists $params{exact_only} ||
        exists $params{line_width} ||
        exists $params{rm_header} ||
        exists $params{seq_id} ||
        exists $params{out_format} ||
        exists $params{stats_info} ||
        exists $params{stats_len} ||
        exists $params{stats_dinuc} ||
        exists $params{stats_tag} ||
        exists $params{stats_dupl} ||
        exists $params{stats_ns} ||
        exists $params{stats_assembly} ||
        exists $params{phred64} ||
        exists $params{no_qual_header} ||
        exists $params{graph_data} ||
        exists $params{custom_params}
        ) {
    &printError('nothing to do with input data');
}
#prevent out of files for stats
if(exists $params{stats_info} || exists $params{stats_len} || exists $params{stats_dinuc} || exists $params{stats_tag} || exists $params{stats_dupl} || exists $params{stats_ns} || exists $params{stats_assembly}) {
    $params{out_good} = 'null';
    $params{out_bad} = 'null';
    $params{stats} = 1;
} elsif(exists $params{out_good} && $params{out_good} eq 'null' && exists $params{out_bad} && $params{out_bad} eq 'null' && !exists $params{graph_data}) {
    &printError('no output selected (both set to null)');
}

#check if FASTQ file is given for option phred64
if(exists $params{phred64}) {
    $command .= ' -phred64';
    unless(exists $params{fastq}) {
        &printError('option -phred64 can only be used for FASTQ input files');
    }
}

#check if output format is possible
if(exists $params{out_format}) {
    $command .= ' -out_format '.$params{out_format};
    if($params{out_format} =~ /\D/) {
        &printError('output format option has to be an integer value');
    } elsif($params{out_format} == 2 || $params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5) {
        unless(exists $params{fastq} || exists $params{qual}) {
            &printError('cannot use this output format option without providing quality data as input');
        }
        if(exists $params{fastq2} && ($params{out_format} == 2 || $params{out_format} == 4 || $params{out_format} == 5)) {
            &printError('cannot use this output format option for paired-end input data. Only values 1 and 3 are allowed');
        }
    } elsif($params{out_format} != 1) {
        &printError('output format option not available');
    }
} else {
    if(exists $params{fastq}) {
        $params{out_format} = 3;
    } elsif(exists $params{fasta} && exists $params{qual}) {
        $params{out_format} = 2;
    } else {
        $params{out_format} = 1;
    }
}
if(exists $params{no_qual_header} && $params{out_format} != 3) {
    &printError('the option -no_qual_header can only be used for FASTQ outputs');
}

#check if output names are different
if(exists $params{out_good} && exists $params{out_bad} && $params{out_good} eq $params{out_bad} && $params{out_good} ne 'null' && $params{out_good} ne 'stdout') {
    &printError('the output names for -out_good and -out_bad have to be different');
}
#check if output can be written to standard output
if(($params{out_format} == 2 || $params{out_format} == 4 || $params{out_format} == 5) && ((exists $params{out_good} && $params{out_good} eq 'stdout') || (exists $params{out_bad} && $params{out_bad} eq 'stdout'))) {
    &printError('the output cannot be written to STDOUT for multiple output files. This option can only be used for FASTA only (-out_format 1) or FASTQ output (-out_format 3)');
}

#check dereplication option
#1 - exact dub, 2 - prefix, 3 - suffix, 4 - revcomp exact, 5 - revcomp prefix/suffix
my $derep = 0;
my %dereptypes;
my $derepmin = 2;
if(exists $params{derep}) {
    $command .= ' -derep '.$params{derep};
    if($params{derep} < 0 || $params{derep} > 54321) {
        &printError('invalid option for dereplication');
    } else {
        my @tmp = split('',$params{derep});
        foreach(@tmp) {
            if($_ < 1 || $_ > 5) {
                &printError('invalid option '.$_.'for dereplication');
            } else {
                $derep = 1;
                $dereptypes{($_-1)} = 0;
            }
        }
    }
}
if(!exists $dereptypes{0} && (exists $dereptypes{1} || exists $dereptypes{2})) {
    $dereptypes{0} = 0;
}
if(!exists $dereptypes{3} && exists $dereptypes{4}) {
    $dereptypes{3} = 0;
}
my $exactonly = 0;
if(exists $params{exact_only}) {
    $command .= ' -exact_only';
    if(exists $dereptypes{1} || exists $dereptypes{2} || exists $dereptypes{4}) {
        &printError('option -exact_only can only be used with -derep options 1 and/or 4');
    }
    if(!exists $params{graph_data}) {
        &printError('option -exact_only requires option -graph_data');
    }
    if(!exists $params{derep}) {
        $dereptypes{0} = 1;
        $dereptypes{3} = 1;
    }
    $exactonly = 1;
} elsif((exists $dereptypes{0} || exists $dereptypes{3}) && !exists $dereptypes{1} && !exists $dereptypes{2} && !exists $dereptypes{4}) {
    $command .= ' -exact_only';
    $exactonly = 1;
    $params{exact_only} = 1;
}
if(exists $params{derep_min}) {
    $command .= ' -derep_min '.$params{derep_min};
    if($params{derep_min} < 2) {
        &printError('invalid option '.$params{derep_min}.'for derep_min. The values has to be greater than 1');
    } elsif(exists $dereptypes{1} || exists $dereptypes{2} || exists $dereptypes{4}) {
        &printError('option -derep_min can only be used with -derep options 1 and/or 4');
    } else {
        $derepmin = $params{derep_min};
    }
}

#check for low complexity method
my $complval;
if(exists $params{lc_method}) {
    $command .= ' -lc_method '.$params{lc_method};
    unless($params{lc_method} eq 'dust' || $params{lc_method} eq 'entropy') {
        &printError('invalid low complexity method');
    }
    unless(exists $params{lc_threshold}) {
        &printError('the low complexity method requires a threshold value specified by -lc_threshold');
    }
    $command .= ' -lc_threshold '.$params{lc_threshold};
    $complval = $params{lc_threshold};
}
if(exists $params{lc_threshold} && !exists $params{lc_method}) {
    &printError('the low complexity threshold requires a method specified by -lc_method');
}

#check for quality trimming
my $trimscore;
if(exists $params{trim_qual_left} || exists $params{trim_qual_right}) {
    $command .= ' -trim_qual_right '.$params{trim_qual_right} if(exists $params{trim_qual_right});
    $command .= ' -trim_qual_left '.$params{trim_qual_left} if(exists $params{trim_qual_left});
    if(exists $params{trim_qual_type}) {
        unless($params{trim_qual_type} eq 'min' || $params{trim_qual_type} eq 'mean' || $params{trim_qual_type} eq 'max' || $params{trim_qual_type} eq 'sum') {
            &printError('invalid value for trim_qual_type');
        }
    } else {
        $params{trim_qual_type} = $TRIM_QUAL_TYPE;
    }
    $command .= ' -trim_qual_type '.$params{trim_qual_type};
    if(exists $params{trim_qual_rule}) {
        unless($params{trim_qual_rule} eq 'lt' || $params{trim_qual_rule} eq 'gt' || $params{trim_qual_rule} eq 'et') {
            &printError('invalid value for trim_qual_rule');
        }
    } else {
        $params{trim_qual_rule} = $TRIM_QUAL_RULE;
    }
    $command .= ' -trim_qual_rule '.$params{trim_qual_rule};
    unless(exists $params{trim_qual_window}) {
        $params{trim_qual_window} = $TRIM_QUAL_WINDOW;
    }
    $command .= ' -trim_qual_window '.$params{trim_qual_window};
    unless(exists $params{trim_qual_step}) {
        $params{trim_qual_step} = $TRIM_QUAL_STEP;
    }
    $command .= ' -trim_qual_step '.$params{trim_qual_step};
    $trimscore = 1;
}

#check sequence case
if(exists $params{seq_case}) {
    $command .= ' -seq_case '.$params{seq_case};
    unless($params{seq_case} eq 'upper' || $params{seq_case} eq 'lower') {
        &printError('invalid sequence case option');
    }
}

#check for dna/rna
if(exists $params{dna_rna}) {
    $command .= ' -dna_rna '.$params{dna_rna};
    unless($params{dna_rna} eq 'dna' || $params{dna_rna} eq 'rna') {
        &printError('invalid option for -dna_rna');
    }
    if($aa) {
        &printError('option -dna_rna cannot be used with option -aa');
    }
}

#set remaining parameters
my $linelen;
if($params{out_format} == 3) {
    $linelen = 0;
} elsif(exists $params{line_width}) {
    $linelen = $params{line_width};
    $command .= ' -line_width '.$params{line_width};
} else {
    $linelen = $LINE_WIDTH;
}

if(exists $params{seq_id}) {
    $command .= ' -seq_id '.$params{seq_id};
    #remove spaces, ">" and quotes from sequence ids
    $params{seq_id} =~ s/[\s\>\"\'\`]//g;
} elsif(exists $params{seq_id_mappings}) {
    &printError('option -seq_id_mappings requires option -seq_id');
}
if(exists $params{seq_id_mappings}) {
    $command .= ' -seq_id_mappings'.($params{seq_id_mappings} ? ' '.$params{seq_id_mappings} : '');
}

my ($repAleft,$repTleft,$repAright,$repTright,$repNleft,$repNright);
if(exists $params{trim_tail_left}) {
    $command .= ' -trim_tail_left '.$params{trim_tail_left};
    $repAleft = 'A'x$params{trim_tail_left};
    $repAleft = qr/^$repAleft/;
    $repTleft = 'T'x$params{trim_tail_left};
    $repTleft = qr/^$repTleft/;
}
if(exists $params{trim_tail_right}) {
    $command .= ' -trim_tail_right '.$params{trim_tail_right};
    $repAright = 'A'x$params{trim_tail_right};
    $repAright = qr/$repAright$/;
    $repTright = 'T'x$params{trim_tail_right};
    $repTright = qr/$repTright$/;
}
if(exists $params{trim_ns_left}) {
    $command .= ' -trim_ns_left '.$params{trim_ns_left};
    $repNleft = 'N'x$params{trim_ns_left};
    $repNleft = qr/^$repNleft/;
}
if(exists $params{trim_ns_right}) {
    $command .= ' -trim_ns_right '.$params{trim_ns_right};
    $repNright = 'N'x$params{trim_ns_right};
    $repNright = qr/$repNright$/;
}

#graph data file
if(exists $params{graph_data}) {
    if(exists $params{stats}) {
        &printError("The graph data cannot be generated at the same time as the statistics");
    }
    $command .= ' -graph_data'.($params{graph_data} ? ' '.$params{graph_data} : '');
    unless($params{graph_data}) {
        $params{graph_data} = join("__",$file1||'nonamegiven').'.gd';
    }
    $params{graph_data} = cwd().'/'.$params{graph_data} unless($params{graph_data} =~ /^\//);
}
my $scale = 1;
if(exists $params{qual_noscale}) {
    $command .= ' -qual_noscale';
    $scale = 0;
}
#graph data selection
if(exists $params{graph_stats} && !exists $params{graph_data}) {
    &printError('option -graph_stats requires option -graph_data');
}
my %webstats = %GRAPH_OPTIONS;
my %graphstats = %GRAPH_OPTIONS;
if(exists $params{graph_stats}) {
    $command .= ' -graph_stats'.($params{graph_stats} ? ' '.$params{graph_stats} : '');
    if($params{graph_stats}) {
        #set all zeroto reset default selection
        foreach my $s (keys %graphstats) {
            $graphstats{$s} = 0;
            $webstats{$s} = 0;
        }
        my @tmp = split(',',$params{graph_stats});
        foreach my $s (@tmp) {
            if(exists $graphstats{$s}) {
                $graphstats{$s} = 1;
            } else {
                &printError('unknown option "'.$s.'" for -graph_stats');
            }
        }
    } else {
        &printError('please specify at least one option for -graph_stats');
    }
}
if(exists $params{graph_stats} && exists $params{web}) {
    &printError('option -graph_stats cannot be used in combination with -web');
}
#web output
my $webnoprocess = 0;
if(exists $params{web}) {
    $command .= ' -web'.($params{web} ? ' '.$params{web} : '');
    if($params{web}) {
        unless($params{web} eq 'process') {
            $webnoprocess = 1;
            foreach my $s (keys %webstats) {
                $webstats{$s} = 0;
                $graphstats{$s} = 0;
            }
            my @tmp = split(',',$params{web});
            foreach my $s (@tmp) {
                $webstats{$s} = 1;
            }
        }
    }
}
if(exists $params{graph_stats} && $graphstats{da}) {
    if($exactonly) {
        &printError('"-exact_only" and "-graph_stats da" cannot be specified at the same time');
    } else {
        $graphstats{de} = 0;
    }
}
#do not calculate all duplicates for paired-end data (at least for now)
if($file2) {
    if(exists $webstats{da}) {
        $webstats{da} = 0;
    }
    if(exists $graphstats{da}) {
        $graphstats{da} = 0;
    }
}
if($webstats{da}) {
    $webstats{de} = 0;
}
if($graphstats{da}) {
    $graphstats{de} = 0;
}
if((exists $params{graph_data} && $graphstats{de}) || (exists $params{web} && $webstats{de})) {
    $exactonly = 1;
}

#custom params
my @cps = ();
if(exists $params{custom_params}) {
    $command .= ' -custom_params "'.$params{custom_params}.'"';
    my ($repeats,@tmp,$bases);
    foreach my $rule (split(/\s*\;\s*/,$params{custom_params})) {
	$repeats = 1;
	@tmp = split(/\s+/,$rule);
	next unless(scalar(@tmp) == 2);
	$bases = ($tmp[0] =~ tr/ACGTN//);
	next if($bases < length($tmp[0]));
	if(index($tmp[1],'%') != -1) {
	    $tmp[1] =~ s/\%//g;
	    $repeats = 0;
	}
	next unless($tmp[1] =~ m/^\d+$/o);
	push(@cps,[$repeats,$tmp[0],$tmp[1]]);
    }
}

#add remaining to log command
if(exists $params{log} || exists $params{graph_data}) {
    if(exists $params{log}) {
        $command .= ' -log'.($params{log} ? ' '.$params{log} : '');
    }
    if(exists $params{min_len}) {
        $command .= ' -min_len '.$params{min_len};
    }
    if(exists $params{max_len}) {
        $command .= ' -max_len '.$params{max_len};
    }
    if(exists $params{range_len}) {
        $command .= ' -range_len '.$params{range_len};
    }
    if(exists $params{min_gc}) {
        $command .= ' -min_gc '.$params{min_gc};
    }
    if(exists $params{max_gc}) {
        $command .= ' -max_gc '.$params{max_gc};
    }
    if(exists $params{range_gc}) {
        $command .= ' -range_gc '.$params{range_gc};
    }
    if(exists $params{min_qual_score}) {
        $command .= ' -min_qual_score '.$params{min_qual_score};
    }
    if(exists $params{max_qual_score}) {
        $command .= ' -max_qual_score '.$params{max_qual_score};
    }
    if(exists $params{min_qual_mean}) {
        $command .= ' -min_qual_mean '.$params{min_qual_mean};
    }
    if(exists $params{max_qual_mean}) {
        $command .= ' -max_qual_mean '.$params{max_qual_mean};
    }
    if(exists $params{ns_max_p}) {
        $command .= ' -ns_max_p '.$params{ns_max_p};
    }
    if(exists $params{ns_max_n}) {
        $command .= ' -ns_max_n '.$params{ns_max_n};
    }
    if(exists $params{noniupac}) {
        $command .= ' -noniupac';
    }
    if(exists $params{seq_num}) {
        $command .= ' -seq_num '.$params{seq_num};
    }
    if(exists $params{trim_to_len}) {
        $command .= ' -trim_to_len '.$params{trim_to_len};
    }
    if(exists $params{trim_left}) {
        $command .= ' -trim_left '.$params{trim_left};
    }
    if(exists $params{trim_right}) {
        $command .= ' -trim_right '.$params{trim_right};
    }
    if(exists $params{trim_left_p}) {
        $command .= ' -trim_left_p '.$params{trim_left_p};
    }
    if(exists $params{trim_right_p}) {
        $command .= ' -trim_right_p '.$params{trim_right_p};
    }
    if(exists $params{rm_header}) {
        $command .= ' -rm_header';
    }
    if(exists $params{stats_len}) {
        $command .= ' -stats_len';
    }
    if(exists $params{stats_dinuc}) {
        $command .= ' -stats_dinuc';
    }
    if(exists $params{stats_info}) {
        $command .= ' -stats_info';
    }
    if(exists $params{stats_tag}) {
        $command .= ' -stats_tag';
    }
    if(exists $params{stats_dupl}) {
        $command .= ' -stats_dupl';
    }
    if(exists $params{stats_ns}) {
        $command .= ' -stats_ns';
    }
    if(exists $params{stats_assembly}) {
        $command .= ' -stats_assembly';
    }
    if(exists $params{verbose}) {
        $command .= ' -verbose';
    }
    if(exists $params{out_good}) {
        $command .= ' -out_good '.$params{out_good};
    }
    if(exists $params{out_bad}) {
        $command .= ' -out_bad '.$params{out_bad};
    }
    if(exists $params{no_qual_header}) {
        $command .= ' -no_qual_header';
    }

    if(exists $params{log}) {
        unless($params{log}) {
            $params{log} = join("__",$file1||'nonamegiven').'.log';
        }
        $params{log} = cwd().'/'.$params{log} unless($params{log} =~ /^\//);
        if(exists $params{web}) {
            &printLog("Executing PRINSEQ using params file");
        } else {
            &printLog("Executing PRINSEQ with command: \"perl prinseq-".$WHAT.".pl".$command."\"");
        }
    }
}

#
################################################################################
## DATA PROCESSING
################################################################################
#

#order of processing:
#seq_num, trim_left, trim_right, trim_left_p, trim_right_p, trim_qual_left, trim_qual_right, trim_tail_left, trim_tail_right, trim_ns_left, trim_ns_right, trim_to_len, min_len, max_len, range_len, min_qual_score, max_qual_score, min_qual_mean, max_qual_mean, min_gc, max_gc, range_gc, ns_max_p, ns_max_n, noniupac, lc_method, derep, seq_id, seq_case, dna_rna, out_format

my $filename = $file1;
while($filename =~ /[\w\d]+\.[\w\d]+$/) {
    $filename =~ s/\.[\w\d]+$//;
    last if($filename =~ /\/[^\.]+$/);
}
my $filename2;
if($file2) {
    $filename2 = $file2;
    while($filename2 =~ /[\w\d]+\.[\w\d]+$/) {
        $filename2 =~ s/\.[\w\d]+$//;
        last if($filename2 =~ /\/[^\.]+$/);
    }
}

#create filehandles for the output data
my ($fhgood,$fhgood2,$fhgood3,$fh2good,$fh2good2,$fhbad,$fhbad2,$fhbad3,$fh2bad,$fhmappings);
my ($filenamegood,$filenamegood2,$filenamegood3,$filename2good,$filename2good2,$filenamebad,$filenamebad2,$filenamebad3,$filename2bad,$filenamemappings);
my ($nogood,$nobad,$stdoutgood,$stdoutbad,$mappings);
$nogood = $nobad = $stdoutgood = $stdoutbad = $mappings = 0;
if(exists $params{out_good}) {
    if($params{out_good} eq 'null') {
        $nogood = 1;
    } elsif($params{out_good} eq 'stdout') {
        $stdoutgood = 1;
    } else {
        if($filename2) {
            #first input file outputs
            open($fhgood,">".$params{out_good}.'_1.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
            $filenamegood = $params{out_good}.'_1.fast'.($params{out_format} == 3 ? 'q' : 'a');
            open($fhgood2,">".$params{out_good}.'_1_singletons.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
            $filenamegood2 = $params{out_good}.'_1_singletons.fast'.($params{out_format} == 3 ? 'q' : 'a');
            #second input file outputs
            open($fh2good,">".$params{out_good}.'_2.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
            $filename2good = $params{out_good}.'_2.fast'.($params{out_format} == 3 ? 'q' : 'a');
            open($fh2good2,">".$params{out_good}.'_2_singletons.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
            $filename2good2 = $params{out_good}.'_2_singletons.fast'.($params{out_format} == 3 ? 'q' : 'a');
        } else {
            open($fhgood,">".$params{out_good}.'.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a')) or &printError('cannot open output file');
            $filenamegood = $params{out_good}.'.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a');
        }
    }
} else {
    if($filename2) {
        $fhgood = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_XXXX',
                                   SUFFIX => '.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a'),
                                   UNLINK => 0);
        $filenamegood = $fhgood->filename;
        $fhgood2 = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_singletons_XXXX',
                                    SUFFIX => '.fast'.($params{out_format} == 3 ? 'q' : 'a'),
                                    UNLINK => 0);
        $filenamegood2 = $fhgood2->filename;
        $fh2good = File::Temp->new( TEMPLATE => $filename2.'_prinseq_good_XXXX',
                                   SUFFIX => '.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a'),
                                   UNLINK => 0);
        $filename2good = $fh2good->filename;
        $fh2good2 = File::Temp->new( TEMPLATE => $filename2.'_prinseq_good_singletons_XXXX',
                                    SUFFIX => '.fast'.($params{out_format} == 3 ? 'q' : 'a'),
                                    UNLINK => 0);
        $filename2good2 = $fh2good2->filename;
    } else {
        $fhgood = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_XXXX',
                                   SUFFIX => '.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a'),
                                   UNLINK => 0);
        $filenamegood = $fhgood->filename;
    }
}
if(exists $params{out_bad}) {
    if($params{out_bad} eq 'null') {
        $nobad = 1;
    } elsif($params{out_bad} eq 'stdout') {
        $stdoutbad = 1;
    } else {
        if($filename2) {
            open($fhbad,">".$params{out_bad}.'_1.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
            $filenamebad = $params{out_bad}.'_1.fast'.($params{out_format} == 3 ? 'q' : 'a');
            open($fh2bad,">".$params{out_bad}.'_2.fast'.($params{out_format} == 3 ? 'q' : 'a')) or &printError('cannot open output file');
            $filename2bad = $params{out_bad}.'_2.fast'.($params{out_format} == 3 ? 'q' : 'a');
        } else {
            open($fhbad,">".$params{out_bad}.'.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a')) or &printError('cannot open output file');
            $filenamebad = $params{out_bad}.'.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a');
        }
    }
} else {
    $fhbad = File::Temp->new( TEMPLATE => $filename.'_prinseq_bad_XXXX',
                              SUFFIX => '.fast'.(($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5 ) ? 'q' : 'a'),
                              UNLINK => 0);
    $filenamebad = $fhbad->filename;
    if($filename2) {
        $fh2bad = File::Temp->new( TEMPLATE => $filename2.'_prinseq_bad_XXXX',
                                   SUFFIX => '.fast'.($params{out_format} ? 'q' : 'a'),
                                   UNLINK => 0);
        $filename2bad = $fh2bad->filename;
    }
}
if($params{out_format} == 2 || $params{out_format} == 5) {
    if(exists $params{out_good}) {
        unless($nogood) {
            open($fhgood2,">".$params{out_good}.'.qual') or &printError('cannot open output file');
            $filenamegood2 = $params{out_good}.'.qual';
        }
    } else {
        $fhgood2 = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_XXXX',
                                    SUFFIX => '.qual',
                                    UNLINK => 0);
        $filenamegood2 = $fhgood2->filename;
    }
    if(exists $params{out_bad}) {
        unless($nobad) {
            open($fhbad2,">".$params{out_bad}.'.qual') or &printError('cannot open output file');
            $filenamebad2 = $params{out_bad}.'.qual';
        }
    } else {
        $fhbad2 = File::Temp->new( TEMPLATE => $filename.'_prinseq_bad_XXXX',
                                   SUFFIX => '.qual',
                                   UNLINK => 0);
        $filenamebad2 = $fhbad2->filename;
    }
}
if($params{out_format} == 4 || $params{out_format} == 5) {
    if(exists $params{out_good}) {
        unless($nogood) {
            open($fhgood3,">".$params{out_good}.'.fasta') or &printError('cannot open output file');
            $filenamegood3 = $params{out_good}.'.fasta';
        }
    } else {
        $fhgood3 = File::Temp->new( TEMPLATE => $filename.'_prinseq_good_XXXX',
                                    SUFFIX => '.fasta',
                                    UNLINK => 0);
        $filenamegood3 = $fhgood3->filename;
    }
    if(exists $params{out_bad}) {
        unless($nobad) {
            open($fhbad3,">".$params{out_bad}.'.fasta') or &printError('cannot open output file');
            $filenamebad3 = $params{out_bad}.'.fasta';
        }
    } else {
        $fhbad3 = File::Temp->new( TEMPLATE => $filename.'_prinseq_bad_XXXX',
                                   SUFFIX => '.fasta',
                                   UNLINK => 0);
        $filenamebad3 = $fhbad3->filename;
    }
}
if(exists $params{seq_id_mappings}) {
    $mappings = 1;
    if($params{seq_id_mappings}) {
        open($fhmappings,">".$params{seq_id_mappings}) or &printError('cannot open output file');
        $filenamemappings = $params{seq_id_mappings};
    } else {
        open($fhmappings,">".$filename.'_prinseq_good.ids') or &printError('cannot open output file');
        $filenamemappings = $filename.'_prinseq_good.ids';
    }
}

my $numlines = 0;
$webnoprocess = 1 if(!$webnoprocess && exists $params{graph_data} && !exists $params{web} && $nogood && $nobad);
my ($progress,$counter,$part);
$progress = 0;
$counter = $part = 1;
if(exists $params{verbose}) {
    print STDERR "Estimate size of input data for status report (this might take a while for large files)\n";
    $numlines = ($file1 eq 'stdin' ? 1 : &getLineNumber($file1));
    $numlines += &getLineNumber($file2) if($file2);
    print STDERR "\tdone\n";
}
if(exists $params{web}) {
    &printWeb("STATUS: Estimate size of input data for status report (this might take a while for large files)");
    $numlines = &getLineNumber($file1) unless($numlines);
    $numlines += &getLineNumber($file2) if($file2)
}
#for progress bar
if($numlines) {
    $part = int($numlines/100);
}

#parse input data
print STDERR "Parse and process input data\n" if(exists $params{verbose});
print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
if(exists $params{web}) {
    &printWeb("STATUS: Parsing and processing input data");
    &printWeb("STATUS: process-status $progress");
    &printLog("Parsing and processing input data");
} else {
    &printLog("Parsing and processing input data: \"".$file1."\"".(exists $params{qual} ? " and \"".$params{qual}."\"" : "").($file2 ? " and \"".$file2."\"" : ""));
}
my $numseqs = 0;
my $numseqs2 = 0;
my $pairs = 0;
my $goodcount = 0;
my $badcount = 0;
my $badcount2 = 0;
my ($tsvfile,$seqid,$header,$seq,$qual,$count,$numbases,$numbases2,$length,$seqgd);
$count = 0;
$qual = '';
$seq = '';
#stats data
my (%stats,%kmers,%odds,%counts,%graphdata,%kmers2,%counts2,%graphdata2,$md5,$md5r,%md5s,%md5sg,$ucseq,$maxlength);
$maxlength = 0;
#parse data
my $seqcount = 0;
my $seqcount1 = 0; #singleton file 1
my $seqcount2 = 0; #singleton file 2
my $seqbases = 0;
my $seqbases1 = 0; #singleton file 1
my $seqbases2 = 0; #singleton file 2
my $badbases = 0;
my $badbases2 = 0;
my (@seqs,@seqsP,@printtmp,$line); #seqs for stats and graph data, seqsP for preprocessing sequences - all used for duplicate checking/removal

if($file2) {
    open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file1 |") or &printError("Could not open file $file1: $!");
    open(FILE2,"perl -pe 's/\r\n|\r/\n/g' < $file2 |") or &printError("Could not open file $file2: $!");
} else {
    if($file1 eq 'stdin') {
        *FILE = *STDIN;
    } else {
        open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file1 |") or &printError("Could not open file $file1: $!");
    }
    if(exists $params{qual}) {
        if($params{qual} eq 'stdin') {
            &printError('QUAL data cannot be read from STDIN');
        } else {
            open(FILE2,"perl -pe 's/\r\n|\r/\n/g' < ".$params{qual}." |") or &printError("Could not open file ".$params{qual}.": $!");
        }
    }
}

my $exists_stats = (exists $params{stats} ? 1 : 0);
my $exists_graphdata = (exists $params{graph_data} ? 1 : 0);

if($file2) {
    my ($seqid2,$seq2,$header2,$qual2,$length2,%tmpids,@tmpdata1,@tmpdata2,$tmpindex,$tmpentry,$skip1,$skip2,$tmpid,$tmpid2);
    $skip1 = $skip2 = 0;
    if(exists $params{fastq}) {
        while(1) {
            ($seq,$seqid,$qual,$header,$tmpid)      = &readEntryFastq(*FILE,$skip1--,$trimnum1);
            ($seq2,$seqid2,$qual2,$header2,$tmpid2) = &readEntryFastq(*FILE2,$skip2--,$trimnum2);
            last unless($seq || $seq2);
            #check if both have same id; if not, store in tmpdata
            if(defined $tmpid && defined $tmpid2 && $tmpid eq $tmpid2) { #same ids
                &processEntryPairedEnd(length($seq),$seq,$seqid,$header,$qual,length($seq2),$seq2,$seqid2,$header2,$qual2);
                if(keys %tmpids) { #empty tmpdata
                    %tmpids = ();
                    while(@tmpdata1) {
                        &processEntryPairedEnd(@{(shift(@tmpdata1))}[0..4]);
                    }
                    while(@tmpdata2) {
                        &processEntryPairedEnd(undef,undef,undef,undef,undef,@{(shift(@tmpdata2))}[0..4]);
                    }
                }
                $skip1 = $skip2 = 0;
            } elsif(defined $tmpid && exists $tmpids{$tmpid}) {
                while(@tmpdata1) {
                    $tmpentry = shift(@tmpdata1);
                    &processEntryPairedEnd(@$tmpentry[0..4]);
                    delete($tmpids{$tmpentry->[5]});
                }
                $tmpindex = $tmpids{$tmpid};
                while($tmpindex-- > 0) {
                    $tmpentry = shift(@tmpdata2);
                    &processEntryPairedEnd(undef,undef,undef,undef,undef,@$tmpentry[0..4]);
                    delete($tmpids{$tmpentry->[5]});
                }
                #correct indices
                $tmpindex = $tmpids{$tmpid}+1;
                foreach my $k (keys %tmpids) {
                    $tmpids{$k} -= $tmpindex;
                }
                &processEntryPairedEnd(length($seq),$seq,$seqid,$header,$qual,@{(shift(@tmpdata2))}[0..4]);
                delete($tmpids{$tmpid});
                if(defined $seqid2) {
                    $tmpids{$tmpid2} = scalar(@tmpdata2);
                    push(@tmpdata2,[length($seq2),$seq2,$seqid2,$header2,$qual2,$tmpid2]);
                }
                #equal out tmpdata1 and tmpdata2
                unless(scalar(@tmpdata1) == scalar(@tmpdata2)) {
                    $skip1 = 0;
                    $skip2 = scalar(@tmpdata2);
                }
            } elsif(defined $seqid2 && exists $tmpids{$tmpid2}) {
                while(@tmpdata2) {
                    $tmpentry = shift(@tmpdata2);
                    &processEntryPairedEnd(undef,undef,undef,undef,undef,@$tmpentry[0..4]);
                    delete($tmpids{$tmpentry->[5]});
                }
                $tmpindex = $tmpids{$tmpid2};
                while($tmpindex-- > 0) {
                    $tmpentry = shift(@tmpdata1);
                    &processEntryPairedEnd(@$tmpentry[0..4]);
                    delete($tmpids{$tmpentry->[5]});
                }
                #correct indices
                $tmpindex = $tmpids{$tmpid2}+1;
                foreach my $k (keys %tmpids) {
                    $tmpids{$k} -= $tmpindex;
                }
                &processEntryPairedEnd(@{(shift(@tmpdata1))}[0..4],length($seq2),$seq2,$seqid2,$header2,$qual2);
                delete($tmpids{$tmpid2});
                if(defined $seqid) {
                    $tmpids{$tmpid} = scalar(@tmpdata1);
                    push(@tmpdata1,[length($seq),$seq,$seqid,$header,$qual,$tmpid]);
                }
                #equal out tmpdata1 and tmpdata2
                unless(scalar(@tmpdata1) == scalar(@tmpdata2)) {
                    $skip1 = scalar(@tmpdata1);
                    $skip2 = 0;
                }
            } else { #store in tmpdata
                if(defined $seqid) {
                    $tmpids{$tmpid} = scalar(@tmpdata1);
                    push(@tmpdata1,[length($seq),$seq,$seqid,$header,$qual,$tmpid]);
                }
                if(defined $seqid2) {
                    $tmpids{$tmpid2} = scalar(@tmpdata2);
                    push(@tmpdata2,[length($seq2),$seq2,$seqid2,$header2,$qual2,$tmpid2]);
                }
            }
            #progress bar stuff
            if($counter > $part) {
                $counter = 1;
                $progress++;
                $progress = 99 if($progress > 99);
                print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                &printWeb("STATUS: process-status $progress");
            }
        }
        #process remaining in tmpdata
        while(@tmpdata1) {
            $tmpentry = shift(@tmpdata1);
            &processEntryPairedEnd(@$tmpentry[0..4]);
            delete($tmpids{$tmpentry->[5]});
        }
        while(@tmpdata2) {
            $tmpentry = shift(@tmpdata2);
            &processEntryPairedEnd(undef,undef,undef,undef,undef,@$tmpentry[0..4]);
            delete($tmpids{$tmpentry->[5]});
        }
    } else { #format is FASTA
        my ($nextid,$nextheader,$nextid2,$nextheader2);
        while(1) {
            ($seq,$seqid,$header,$tmpid,$nextid,$nextheader)        = &readEntryFasta(*FILE,$skip1--,$trimnum1,$nextid,$nextheader);
            ($seq2,$seqid2,$header2,$tmpid2,$nextid2,$nextheader2)  = &readEntryFasta(*FILE2,$skip2--,$trimnum2,$nextid2,$nextheader2);
            last unless($seq || $seq2);
            #check if both have same id; if not, store in tmpdata
            if(defined $tmpid && defined $tmpid2 && $tmpid eq $tmpid2) { #same ids
                &processEntryPairedEnd(length($seq),$seq,$seqid,$header,undef,length($seq2),$seq2,$seqid2,$header2,undef);
                if(keys %tmpids) { #empty tmpdata
                    %tmpids = ();
                    while(@tmpdata1) {
                        &processEntryPairedEnd(@{(shift(@tmpdata1))}[0..3]);
                    }
                    while(@tmpdata2) {
                        &processEntryPairedEnd(undef,undef,undef,undef,undef,@{(shift(@tmpdata2))}[0..3]);
                    }
                }
                $skip1 = $skip2 = 0;
            } elsif(defined $tmpid && exists $tmpids{$tmpid}) {
                while(@tmpdata1) {
                    $tmpentry = shift(@tmpdata1);
                    &processEntryPairedEnd(@$tmpentry[0..3]);
                    delete($tmpids{$tmpentry->[4]});
                }
                $tmpindex = $tmpids{$tmpid};
                while($tmpindex-- > 0) {
                    $tmpentry = shift(@tmpdata2);
                    &processEntryPairedEnd(undef,undef,undef,undef,undef,@$tmpentry[0..3]);
                    delete($tmpids{$tmpentry->[4]});
                }
                #correct indices
                $tmpindex = $tmpids{$tmpid}+1;
                foreach my $k (keys %tmpids) {
                    $tmpids{$k} -= $tmpindex;
                }
                &processEntryPairedEnd(length($seq),$seq,$seqid,$header,undef,@{(shift(@tmpdata2))}[0..3]);
                delete($tmpids{$tmpid});
                if(defined $seqid2) {
                    $tmpids{$tmpid2} = scalar(@tmpdata2);
                    push(@tmpdata2,[length($seq2),$seq2,$seqid2,$header2,$tmpid2]);
                }
                #equal out tmpdata1 and tmpdata2
                unless(scalar(@tmpdata1) == scalar(@tmpdata2)) {
                    $skip1 = 0;
                    $skip2 = scalar(@tmpdata2);
                }
            } elsif(defined $seqid2 && exists $tmpids{$tmpid2}) {
                while(@tmpdata2) {
                    $tmpentry = shift(@tmpdata2);
                    &processEntryPairedEnd(undef,undef,undef,undef,undef,@$tmpentry[0..3]);
                    delete($tmpids{$tmpentry->[4]});
                }
                $tmpindex = $tmpids{$tmpid2};
                while($tmpindex-- > 0) {
                    $tmpentry = shift(@tmpdata1);
                    &processEntryPairedEnd(@$tmpentry[0..3]);
                    delete($tmpids{$tmpentry->[4]});
                }
                #correct indices
                $tmpindex = $tmpids{$tmpid2}+1;
                foreach my $k (keys %tmpids) {
                    $tmpids{$k} -= $tmpindex;
                }
                &processEntryPairedEnd(@{(shift(@tmpdata1))}[0..3],undef,length($seq2),$seq2,$seqid2,$header2);
                delete($tmpids{$tmpid2});
                if(defined $seqid) {
                    $tmpids{$tmpid} = scalar(@tmpdata1);
                    push(@tmpdata1,[$length,$seq,$seqid,$header,$tmpid]);
                }
                #equal out tmpdata1 and tmpdata2
                unless(scalar(@tmpdata1) == scalar(@tmpdata2)) {
                    $skip1 = scalar(@tmpdata1);
                    $skip2 = 0;
                }
            } else { #store in tmpdata
                if(defined $seqid) {
                    $tmpids{$tmpid} = scalar(@tmpdata1);
                    push(@tmpdata1,[$length,$seq,$seqid,$header,$tmpid]);
                }
                if(defined $seqid2) {
                    $tmpids{$tmpid2} = scalar(@tmpdata2);
                    push(@tmpdata2,[length($seq2),$seq2,$seqid2,$header2,$tmpid2]);
                }
            }
            #progress bar stuff
            if($counter > $part) {
                $counter = 1;
                $progress++;
                $progress = 99 if($progress > 99);
                print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                &printWeb("STATUS: process-status $progress");
            }
        }
        #process remaining in tmpdata
        while(@tmpdata1) {
            $tmpentry = shift(@tmpdata1);
            &processEntryPairedEnd(@$tmpentry[0..3]);
            delete($tmpids{$tmpentry->[4]});
        }
        while(@tmpdata2) {
            $tmpentry = shift(@tmpdata2);
            &processEntryPairedEnd(undef,undef,undef,undef,undef,@$tmpentry[0..3]);
            delete($tmpids{$tmpentry->[4]});
        }
    }
} else {
    if(exists $params{fastq}) {
        foreach(@dataread) { #only used for stdin inputs
            chomp();
            if($count == 0 && /^\@(\S+)\s*(.*)$/o) {
                $length = length($seq);
                if($length) {
                    &processEntry($length,$seq,$seqid,$qual,$header);
                }
                $seqid = $1;
                $header = $2 || '';
                $seq = '';
                $qual = '';
            } elsif($count == 1) {
                $seq = $_;
            } elsif($count == 3) {
                $qual = $_;
                $count = -1;
            }
            $count++;
        }
        while(<FILE>) {
            chomp();
            if($count == 0 && /^\@(\S+)\s*(.*)$/o) {
                $length = length($seq);
                if($length) {
                    &processEntry($length,$seq,$seqid,$qual,$header);
                }
                $seqid = $1;
                $header = $2 || '';
                $seq = '';
                $qual = '';
            } elsif($count == 1) {
                $seq = $_;
            } elsif($count == 3) {
                $qual = $_;
                $count = -1;
            }
            $count++;
            #progress bar stuff
            $counter++;
            if($counter > $part) {
                $counter = 1;
                $progress++;
                $progress = 99 if($progress > 99);
                print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                &printWeb("STATUS: process-status $progress");
            }
        }
        #add last one
        $length = length($seq);
        if($length) {
            &processEntry($length,$seq,$seqid,$qual,$header);
        }
    } elsif(exists $params{fasta}) {
        foreach(@dataread) { #only used for stdin inputs
            chomp();
            if(/^>(\S+)\s*(.*)$/o) {
                $length = length($seq);
                if($length) {
                    #get qual data if provided
                    if(exists $params{qual}) {
                        while(<FILE2>) {
                            chomp();
                            last if(/^>/ && $qual);
                            next if(/^>/);
                            $qual .= $_.' ';
                        }
                        $qual = &convertQualNumsToAsciiString($qual);
                    }
                    &processEntry($length,$seq,$seqid,$qual,$header);
                    $qual = '';
                }
                $seqid = $1;
                $header = $2 || '';
                $seq = '';
            } else {
                $seq .= $_;
            }
        }
        while(<FILE>) {
            chomp();
            if(/^>(\S+)\s*(.*)$/o) {
                $length = length($seq);
                if($length) {
                    #get qual data if provided
                    if(exists $params{qual}) {
                        while(<FILE2>) {
                            chomp();
                            last if(/^>/ && $qual);
                            next if(/^>/);
                            $qual .= $_.' ';
                        }
                        $qual = &convertQualNumsToAsciiString($qual);
                    }
                    &processEntry($length,$seq,$seqid,$qual,$header);
                    $qual = '';
                }
                $seqid = $1;
                $header = $2 || '';
                $seq = '';
            } else {
                $seq .= $_;
            }
            #progress bar stuff
            $counter++;
            if($counter > $part) {
                $counter = 1;
                $progress++;
                $progress = 99 if($progress > 99);
                print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                &printWeb("STATUS: process-status $progress");
            }
        }
        #add last one
        $length = length($seq);
        if($length) {
            #get qual data if provided
            if(exists $params{qual}) {
                while(<FILE2>) {
                    chomp();
                    last if(/^>/ && $qual);
                    next if(/^>/);
                    $qual .= $_.' ';
                }
                $qual = &convertQualNumsToAsciiString($qual);
            }
            &processEntry($length,$seq,$seqid,$qual,$header);
            $qual = '';
        }
    }
}

print STDERR "\r\tdone          \n" if(exists $params{verbose});
&printWeb("STATUS: process-status 100");
&printWeb("Done parsing and processing input data");
&printLog("Done parsing and processing input data");

if($derep && !exists $params{stats} && !$exactonly && !$webnoprocess) {
    &printLog("Remove duplicates");
    &derepSeqs();
    &printLog("Done removing duplicates");
}

close(FILE);
close(FILE2) if(exists $params{qual} || $file2);

print STDERR "Clean up empty files\n" if(exists $params{verbose} && !exists $params{stats});
#close filehandles
unless($nogood || $stdoutgood) {
    close($fhgood);
    if($file2) {
        close($fhgood2);
        close($fh2good);
        close($fh2good2);
    }
}
unless($nobad || $stdoutbad) {
    close($fhbad);
    if($file2) {
        close($fh2bad);
    }
}
if($params{out_format} == 2 || $params{out_format} == 5) {
    close($fhgood2) unless($nogood || $stdoutgood);
    close($fhbad2) unless($nobad || $stdoutbad);
}
if($params{out_format} == 4 || $params{out_format} == 5) {
    close($fhgood3) unless($nogood || $stdoutgood);
    close($fhbad3) unless($nobad || $stdoutbad);
}
if(exists $params{seq_id_mappings}) {
    close($fhmappings);
}
#remove empty files
if($seqcount == 0 && !$nogood && !$stdoutgood) {
    unlink($filenamegood);
    if($params{out_format} == 2 || $params{out_format} == 5) {
        unlink($filenamegood2);
    }
    if($params{out_format} == 4 || $params{out_format} == 5) {
        unlink($filenamegood3);
    }
    unlink($fhmappings);
}
if($badcount == 0 && !$nobad && !$stdoutbad) {
    unlink($filenamebad);
    if($params{out_format} == 2 || $params{out_format} == 5) {
        unlink($filenamebad2);
    }
    if($params{out_format} == 4 || $params{out_format} == 5) {
        unlink($filenamebad3);
    }
}
if($file2) {
    foreach my $f ($filenamegood2,$filename2good,$filename2good2,$filename2bad) {
        if(defined $f && -e $f && &getLineNumber($f) == 0) {
            unlink($f);
        }
    }
}
print STDERR "\tdone\n" if(exists $params{verbose} && !exists $params{stats});

#print processing stats
if((exists $params{verbose} || !$webnoprocess) && !exists $params{stats}) {
    print STDERR "Input and filter stats:\n";
    if($file2) {
        print STDERR "\tInput sequences (file 1): ".&addCommas($numseqs)."\n";
        print STDERR "\tInput bases (file 1): ".&addCommas($numbases)."\n";
        print STDERR "\tInput mean length (file 1): ".sprintf("%.2f",$numbases/$numseqs)."\n" if($numseqs);
        print STDERR "\tInput sequences (file 2): ".&addCommas($numseqs2)."\n";
        print STDERR "\tInput bases (file 2): ".&addCommas($numbases2)."\n";
        print STDERR "\tInput mean length (file 2): ".sprintf("%.2f",$numbases2/$numseqs2)."\n" if($numseqs2);
        print STDERR "\tGood sequences (pairs): ".&addCommas($seqcount)."\n" if($numseqs);
        print STDERR "\tGood bases (pairs): ".&addCommas($seqbases)."\n" if($seqcount);
        print STDERR "\tGood mean length (pairs): ".sprintf("%.2f",$seqbases/$seqcount)."\n" if($seqcount);
        print STDERR "\tGood sequences (singletons file 1): ".&addCommas($seqcount1)." (".sprintf("%.2f",(100*$seqcount1/$numseqs))."%)\n" if($numseqs);
        print STDERR "\tGood bases (singletons file 1): ".&addCommas($seqbases1)."\n" if($seqcount1);
        print STDERR "\tGood mean length (singletons file 1): ".sprintf("%.2f",$seqbases1/$seqcount1)."\n" if($seqcount1);
        print STDERR "\tGood sequences (singletons file 2): ".&addCommas($seqcount2)." (".sprintf("%.2f",(100*$seqcount2/$numseqs2))."%)\n" if($numseqs2);
        print STDERR "\tGood bases (singletons file 2): ".&addCommas($seqbases2)."\n" if($seqcount2);
        print STDERR "\tGood mean length (singletons file 2): ".sprintf("%.2f",$seqbases2/$seqcount2)."\n" if($seqcount2);
        print STDERR "\tBad sequences (file 1): ".&addCommas($badcount)." (".sprintf("%.2f",(100*$badcount/$numseqs))."%)\n" if($numseqs);
        print STDERR "\tBad bases (file 1): ".&addCommas($badbases)."\n" if($badcount);
        print STDERR "\tBad mean length (file 1): ".sprintf("%.2f",$badbases/$badcount)."\n" if($badcount);
        print STDERR "\tBad sequences (file 2): ".&addCommas($badcount2)." (".sprintf("%.2f",(100*$badcount2/$numseqs2))."%)\n" if($numseqs2);
        print STDERR "\tBad bases (file 2): ".&addCommas($badbases2)."\n" if($badcount2);
        print STDERR "\tBad mean length (file 2): ".sprintf("%.2f",$badbases2/$badcount2)."\n" if($badcount2);
    } else {
        print STDERR "\tInput sequences: ".&addCommas($numseqs)."\n";
        print STDERR "\tInput bases: ".&addCommas($numbases)."\n";
        print STDERR "\tInput mean length: ".sprintf("%.2f",$numbases/$numseqs)."\n" if($numseqs);
        print STDERR "\tGood sequences: ".&addCommas($seqcount)." (".sprintf("%.2f",(100*$seqcount/$numseqs))."%)\n" if($numseqs);
        print STDERR "\tGood bases: ".&addCommas($seqbases)."\n" if($seqcount);
        print STDERR "\tGood mean length: ".sprintf("%.2f",$seqbases/$seqcount)."\n" if($seqcount);
        print STDERR "\tBad sequences: ".&addCommas($badcount)." (".sprintf("%.2f",(100*$badcount/$numseqs))."%)\n" if($numseqs);
        print STDERR "\tBad bases: ".&addCommas($badbases)."\n" if($badcount);
        print STDERR "\tBad mean length: ".sprintf("%.2f",$badbases/$badcount)."\n" if($badcount);
    }
    my $tmp = &getFiltercounts();
    if($tmp) {
        print STDERR "\tSequences filtered by specified parameters:\n";
        if(scalar(@$tmp)) {
            print STDERR "\t$_\n" foreach(@$tmp);
        } else {
            print STDERR "\tnone\n";
        }
    }
    if(exists $params{log}) {
        if($file2) {
            &printLog("Input sequences (file 1): ".&addCommas($numseqs));
            &printLog("Input bases (file 1): ".&addCommas($numbases));
            &printLog("Input mean length (file 1): ".sprintf("%.2f",$numbases/$numseqs)) if($numseqs);
            &printLog("Input sequences (file 2): ".&addCommas($numseqs2));
            &printLog("Input bases (file 2): ".&addCommas($numbases2));
            &printLog("Input mean length (file 2): ".sprintf("%.2f",$numbases2/$numseqs2)) if($numseqs2);
            &printLog("Good sequences (pairs): ".&addCommas($seqcount)) if($numseqs);
            &printLog("Good bases (pairs): ".&addCommas($seqbases)) if($seqcount);
            &printLog("Good mean length (pairs): ".sprintf("%.2f",$seqbases/$seqcount)) if($seqcount);
            &printLog("Good sequences (singletons file 1): ".&addCommas($seqcount1)." (".sprintf("%.2f",(100*$seqcount1/$numseqs))."%)") if($numseqs);
            &printLog("Good bases (singletons file 1): ".&addCommas($seqbases1)) if($seqcount1);
            &printLog("Good mean length (singletons file 1): ".sprintf("%.2f",$seqbases1/$seqcount1)) if($seqcount1);
            &printLog("Good sequences (singletons file 2): ".&addCommas($seqcount2)." (".sprintf("%.2f",(100*$seqcount2/$numseqs2))."%)") if($numseqs2);
            &printLog("Good bases (singletons file 2): ".&addCommas($seqbases2)) if($seqcount2);
            &printLog("Good mean length (singletons file 2): ".sprintf("%.2f",$seqbases2/$seqcount2)) if($seqcount2);
            &printLog("Bad sequences (file 1): ".&addCommas($badcount)." (".sprintf("%.2f",(100*$badcount/$numseqs))."%)") if($numseqs);
            &printLog("Bad bases (file 1): ".&addCommas($badbases)) if($badcount);
            &printLog("Bad mean length (file 1): ".sprintf("%.2f",$badbases/$badcount)) if($badcount);
            &printLog("Bad sequences (file 2): ".&addCommas($badcount2)." (".sprintf("%.2f",(100*$badcount2/$numseqs2))."%)") if($numseqs2);
            &printLog("Bad bases (file 2): ".&addCommas($badbases2)) if($badcount2);
            &printLog("Bad mean length (file 2): ".sprintf("%.2f",$badbases2/$badcount2)) if($badcount2);
        } else {
            &printLog("Input sequences: ".&addCommas($numseqs));
            &printLog("Input bases: ".&addCommas($numbases));
            &printLog("Input mean length: ".sprintf("%.2f",$numbases/$numseqs)) if($numseqs);
            &printLog("Good sequences: ".&addCommas($seqcount)." (".sprintf("%.2f",(100*$seqcount/$numseqs))."%)") if($numseqs);
            &printLog("Good bases: ".&addCommas($seqbases)) if($seqcount);
            &printLog("Good mean length: ".sprintf("%.2f",$seqbases/$seqcount)) if($seqcount);
            &printLog("Bad sequences: ".&addCommas($badcount)." (".sprintf("%.2f",(100*$badcount/$numseqs))."%)") if($numseqs);
            &printLog("Bad bases: ".&addCommas($badbases)) if($badcount);
            &printLog("Bad mean length: ".sprintf("%.2f",$badbases/$badcount)) if($badcount);
        }
        if($tmp) {
            &printLog("Sequences filtered by specified parameters:");
            if(scalar(@$tmp)) {
                &printLog($_) foreach(@$tmp);
            } else {
                &printLog("none");
            }

        }
    }
}

#print summary stats
if(exists $params{stats}) {
    if(exists $params{stats_info}) {
        $stats{stats_info}->{reads} = $numseqs;
        $stats{stats_info}->{bases} = $numbases;
        if($file2) {
            $stats{stats_info2}->{reads} = $numseqs2;
            $stats{stats_info2}->{bases} = $numbases2;
        }
    }
    if(exists $params{stats_len}) {
        $stats{stats_len} = &generateStats($counts{length});
        $stats{stats_len2} = &generateStats($counts2{length}) if($file2);
    }
    if(exists $params{stats_dinuc}) {
        foreach my $i (keys %odds) {
	    $stats{stats_dinuc}->{lc($i)} = sprintf("%.9f",$odds{$i}/($numseqs+$numseqs2));
	}
    }
    if(exists $params{stats_tag}) {
        #calculate frequency of 5-mers
        my $kmersum = &getTagFrequency(\%kmers);
        #check for frequency of MID tags
        my $midsum = 0;
        my $midcount = 0;
        my @midseqs;
        foreach my $mid (keys %MIDS) {
            $midsum += $MIDS{$mid};
            if($MIDS{$mid} > $numseqs/34) { #in more than 34/100 (approx. 3%) as this is estimated average error for MIDs
                $midcount++;
                push(@midseqs,$mid);
            }
        }
        $stats{stats_tag}->{midnum} = $midcount;
        if($midcount) {
            $stats{stats_tag}->{midseq} = join(",",@midseqs);
        }
        if($midsum > $kmersum->{5}) {
            $kmersum->{5} = $midsum;
        }
        foreach my $kmer (keys %$kmersum) {
            $stats{stats_tag}->{'prob'.$kmer} = sprintf("%d",(100/$numseqs*$kmersum->{$kmer}));
	}
        if($file2) {
            my $kmersum2 = &getTagFrequency(\%kmers2);
            foreach my $kmer (keys %$kmersum2) {
                $stats{stats_tag2}->{'prob'.$kmer} = sprintf("%d",(100/$numseqs2*$kmersum2->{$kmer}));
            }
        }
    }
    if(exists $params{stats_assembly}) {
        #calculate N50, N90, etc
        #sort in decreasing order
	my @sortvals = sort {$b <=> $a} keys %{$counts{length}};
	#calculate nx values
	my $n50 = $numbases*0.5;
	my $n75 = $numbases*0.75;
	my $n90 = $numbases*0.9;
	my $n95 = $numbases*0.95;
	my $curlen = 0;
	foreach my $len (@sortvals) {
	    foreach my $i (1..$counts{length}->{$len}) {
		$curlen += $len;
		if($curlen >= $n50 && !exists $stats{stats_assembly}->{N50}) {
		    $stats{stats_assembly}->{N50} = $len;
		} elsif($curlen >= $n75 && !exists $stats{stats_assembly}->{N75}) {
		    $stats{stats_assembly}->{N75} = $len;
		} elsif($curlen >= $n90 && !exists $stats{stats_assembly}->{N90}) {
		    $stats{stats_assembly}->{N90} = $len;
		} elsif($curlen >= $n95 && !exists $stats{stats_assembly}->{N95}) {
		    $stats{stats_assembly}->{N95} = $len;
		}
	    }
	}
        foreach my $i (50,75,90,95) {
            unless(exists $stats{stats_assembly}->{'N'.$i}) {
                $stats{stats_assembly}->{'N'.$i} = '-';
            }
        }
    }
    if(exists $params{stats_dupl}) {
        #empty vars before n-plicate check
        %counts = %kmers = %odds = ();
        #0 - exact dub, 1 - prefix, 2 - suffix, 3 - revcomp exact, 4 - revcomp prefix/suffix
        my %types = (0 => 'exact', 1 => '5', 2 => '3', 3 => 'exactrevcomp', 4 => 'revcomp');
        my ($dupls,undef,undef) = &checkForDupl(\@seqs,\%types,$numseqs);
        #set zero counts
        foreach my $s (keys %types) {
            $stats{stats_dupl}->{$types{$s}} = 0;
            $stats{stats_dupl}->{$types{$s}.'maxd'} = 0;
        }
        foreach my $n (keys %$dupls) {
	    foreach my $s (keys %{$dupls->{$n}}) {
		$stats{stats_dupl}->{$types{$s}} += $dupls->{$n}->{$s} * $n;
		$stats{stats_dupl}->{$types{$s}.'maxd'} = $n unless($stats{stats_dupl}->{$types{$s}.'maxd'} > $n);
		$stats{stats_dupl}->{total} += $dupls->{$n}->{$s} * $n;
	    }
	}
    }
    foreach my $type (sort keys %stats) {
        foreach my $value (sort keys %{$stats{$type}}) {
            print STDOUT join("\t",$type,$value,(defined $stats{$type}->{$value} ? $stats{$type}->{$value} : '-'))."\n";
        }
    }
}

if(exists $params{graph_data}) {
    &printLog("Generate graph data");
    print STDERR "Generate graph data\n" if(exists $params{verbose});
    &printWeb("Start generating statistics from data");
    #get qual stats
    my $binval = &getBinVal($maxlength);
    if($graphdata{quals} || $graphdata{quala}) {
        &printWeb("Generating statistics from quality data");
        if($graphdata{quals}) {
            $graphdata{quals} = &generateStatsType($graphdata{quals});
        }
        if($graphdata{quala}) {
            #calculate bin values
            my $tmppos;
            foreach my $pos (keys %{$graphdata{quala}}) {
                $tmppos = int(($pos-1)/$binval);
                foreach my $val (keys %{$graphdata{quala}->{$pos}}) {
                    $graphdata{qualsbin}->{$tmppos}->{$val} += $graphdata{quala}->{$pos}->{$val};
                }
            }
            $graphdata{qualsbin} = &generateStatsType($graphdata{qualsbin});
        }
    }
    if($file2 && ($graphdata2{quals} || $graphdata2{quala})) {
        $graphdata{qualsmean2} = $graphdata2{qualsmean};
        if($graphdata2{quals}) {
            $graphdata{quals2} = &generateStatsType($graphdata2{quals});
        }
        if($graphdata2{quala}) {
            #calculate bin values
            my $tmppos;
            foreach my $pos (keys %{$graphdata2{quala}}) {
                $tmppos = int(($pos-1)/$binval);
                foreach my $val (keys %{$graphdata2{quala}->{$pos}}) {
                    $graphdata{qualsbin2}->{$tmppos}->{$val} += $graphdata2{quala}->{$pos}->{$val};
                }
            }
            $graphdata{qualsbin2} = &generateStatsType($graphdata{qualsbin2});
        }
    }
    #get length stats
    if($graphdata{counts}) {
        &printWeb("Generating statistics from basic counts");
        $graphdata{stats} = &generateStatsType($graphdata{counts});
    }
    if($file2 && $graphdata2{counts}) {
        $graphdata{counts2} = $graphdata2{counts};
        $graphdata{stats2} = &generateStatsType($graphdata2{counts});
    }
    #check for ns
    if(($webstats{ns} || $graphstats{ns}) && scalar(keys %{$graphdata{counts}->{ns}}) == 0) {
        $graphdata{counts}->{ns}->{0} = 0;
    }
    if($file2 && ($webstats{ns} || $graphstats{ns}) && scalar(keys %{$graphdata{counts2}->{ns}}) == 0) {
        $graphdata{counts2}->{ns}->{0} = 0;
    }
    #add base frequencies
    foreach my $site (keys %{$graphdata{freqs}}) {
	foreach my $i (0..$TAG_LENGTH-1) {
	    foreach my $base ('A','C','G','T','N') {
		if(exists $graphdata{freqs}->{$site}->{$i}->{$base}) {
		    $graphdata{freqs}->{$site}->{$i}->{$base} = int($graphdata{freqs}->{$site}->{$i}->{$base}*100/$numseqs);
		} else {
		    $graphdata{freqs}->{$site}->{$i}->{$base} = 0;
		}
	    }
	}
    }
    if($file2) {
        foreach my $site (keys %{$graphdata2{freqs}}) {
            foreach my $i (0..$TAG_LENGTH-1) {
                foreach my $base ('A','C','G','T','N') {
                    if(exists $graphdata2{freqs}->{$site}->{$i}->{$base}) {
                        $graphdata{freqs2}->{$site}->{$i}->{$base} = int($graphdata2{freqs}->{$site}->{$i}->{$base}*100/$numseqs2);
                    } else {
                        $graphdata{freqs2}->{$site}->{$i}->{$base} = 0;
                    }
                }
            }
        }
    }
    #calculate possibility for tag sequences
    if(scalar(keys %{$graphdata{kmers}})) {
        &printWeb("Generating statistics for tag sequences");
	my %prob;
	#calculate frequency of 5-mers
        my $kmersum = &getTagFrequency($graphdata{kmers},$numseqs);
        #check for frequency of MID tags
        my $midsum = 0;
        my $midcount = 0;
        my @midseqs;
        foreach my $mid (keys %{$graphdata{mids}}) {
            $midsum += $graphdata{mids}->{$mid};
            if($graphdata{mids}->{$mid} > $numseqs/34) { #in more than 34/100 (approx. 3%) as this is estimated average error for MIDs
                $midcount++;
                push(@midseqs,$mid);
            }
        }
        if($midcount) {
            $graphdata{tagmidseq} = join(",",@midseqs);
        }
	$graphdata{tagmidnum} = $midcount;
        if($midsum > $kmersum->{5}) {
            $kmersum->{5} = $midsum;
        }
        foreach my $kmer (keys %$kmersum) {
            $prob{$kmer} = sprintf("%d",(100/$numseqs*$kmersum->{$kmer}));
	}
        $graphdata{tagprob} = \%prob;
        delete($graphdata{kmers});
    }
    if($file2 && scalar(keys %{$graphdata2{kmers}})) {
	my %prob;
	#calculate frequency of 5-mers
        my $kmersum = &getTagFrequency($graphdata2{kmers},$numseqs2);
        foreach my $kmer (keys %$kmersum) {
            $prob{$kmer} = sprintf("%d",(100/$numseqs2*$kmersum->{$kmer}));
	}
        $graphdata{tagprob2} = \%prob;
        delete($graphdata2{kmers});
    }
    #add dinucleotide odd ratios
    if(($webstats{dn} || $graphstats{dn}) && scalar(keys %{$graphdata{dinucodds}})) {
        &printWeb("Generating statistics from dinucleotide counts");
	foreach my $i (keys %{$graphdata{dinucodds}}) {
	    $graphdata{dinucodds}->{$i} = sprintf("%.9f",$graphdata{dinucodds}->{$i}/($numseqs+$numseqs2));
	}
    }
    #check for n-plicates (for paired-end data, not separated by input file)
    if($webstats{de} || $webstats{da} || $graphstats{de} || $graphstats{da}) {
        &printWeb("Generating statistics from duplicate counts");
        if($webstats{de} || $graphstats{de}) {
            foreach my $m (keys %md5sg) {
                if(exists $md5sg{$m}->{0} && $md5sg{$m}->{0} > 0) {
                    $graphdata{dubscounts}->{$md5sg{$m}->{0}}->{0}++;
                }
                if(exists $md5sg{$m}->{3} && $md5sg{$m}->{3} > 0) {
                    $graphdata{dubscounts}->{$md5sg{$m}->{3}}->{3}++;
                }
            }
        } else {
            my %types = (0 => 0, 1 => 0, 2 => 0, 3 => 0, 4 => 0);
            ($graphdata{dubscounts},$graphdata{dubslength},undef) = &checkForDupl(\@seqs,\%types,$numseqs); #0 - exact dub, 1 - prefix, 2 - suffix, 3 - revcomp exact, 4 - revcomp prefix/suffix
        }
    }

    #generate JSON string without the need of the JSON module
    my $str = '';
    $str .= '{"numseqs":'.$numseqs.',"numbases":'.$numbases.',"pairedend":'.($file2 ? 1 : 0).($file2 ? ',"numseqs2":'.$numseqs2.',"numbases2":'.$numbases2.',"pairs":'.$pairs : '').',"maxlength":'.$maxlength.',"binval":'.$binval.',"exactonly":'.$exactonly.',"tagmidnum":'.($graphdata{tagmidnum}||0).',"scale":'.$scale.',"filename1":"'.(exists $params{filename1} ? $params{filename1} : &convertStringToInt(&getFileName($file1))).'","format1":"'.(exists $params{fasta} ? 'fasta' : 'fastq').'"'.(exists $params{qual} ? ',"filename2":"'.(exists $params{filename2} ? $params{filename2} : &convertStringToInt(&getFileName($params{qual}))).'","format2":"qual"' : '').($file2 ? ',"filename2":"'.(exists $params{filename2} ? $params{filename2} : &convertStringToInt(&getFileName($file2))).'"' : '');
    foreach my $s (qw(counts counts2 stats stats2 quals quals2 qualsbin qualsbin2 complvals dubscounts dubslength)) {
        next unless exists($graphdata{$s});
        $str .= ',"'.$s.'":{';
        foreach my $t (sort keys %{$graphdata{$s}}) {
            $str .= '"'.$t.'":{';
            while (my ($k,$v) = each(%{$graphdata{$s}->{$t}})) {
                if($v =~ /^\d+$/) {
                    $str .= '"'.$k.'":'.$v.',';
                } else {
                    $str .= '"'.$k.'":"'.$v.'",';
                }
            }
            $str =~ s/\,$//;
            $str .= '},';
        }
        $str =~ s/\,$//;
        $str .= '}';
    }
    foreach my $s (qw(qualsmean qualsmean2 tagprob tagprob2 compldust complentropy dinucodds)) {
        next unless exists($graphdata{$s});
        $str .= ',"'.$s.'":{';
        while (my ($k,$v) = each(%{$graphdata{$s}})) {
            if($v =~ /^\d+$/) {
                $str .= '"'.$k.'":'.$v.',';
            } else {
                $str .= '"'.$k.'":"'.$v.'",';
            }
        }
        $str =~ s/\,$//;
        $str .= '}';
    }
    $str .= ',"tail":'.(exists $graphdata{counts}->{tail5} || exists $graphdata{counts}->{tail3} ? 1 : 0);
    $str .= ',"tail2":'.(exists $graphdata2{counts}->{tail5} || exists $graphdata2{counts}->{tail3} ? 1 : 0) if($file2);
    if($webstats{ts} || $graphstats{ts}) {
        $str .= ',"freqs":{';
        foreach my $i (keys %{$graphdata{freqs}}) { # 5, 3
            $str .= '"'.$i.'":{';
            foreach my $pos (keys %{$graphdata{freqs}->{$i}}) {
                $str .= '"'.$pos.'":{';
                while (my ($base,$v) = each(%{$graphdata{freqs}->{$i}->{$pos}})) {
                    $str .= '"'.$base.'":'.$v.',';
                }
                $str =~ s/\,$//;
                $str .= '},';
            }
            $str =~ s/\,$//;
            $str .= '},';
        }
        $str =~ s/\,$//;
        $str .= '}';
        if($file2) {
            $str .= ',"freqs2":{';
            foreach my $i (keys %{$graphdata{freqs2}}) { # 5, 3
                $str .= '"'.$i.'":{';
                foreach my $pos (keys %{$graphdata{freqs2}->{$i}}) {
                    $str .= '"'.$pos.'":{';
                    while (my ($base,$v) = each(%{$graphdata{freqs2}->{$i}->{$pos}})) {
                        $str .= '"'.$base.'":'.$v.',';
                    }
                    $str =~ s/\,$//;
                    $str .= '},';
                }
                $str =~ s/\,$//;
                $str .= '},';
            }
            $str =~ s/\,$//;
            $str .= '}';
        }
    }
    foreach my $s (qw(tagmidseq)) {
        next unless exists($graphdata{$s});
        $str .= ',"'.$s.'":"'.$graphdata{$s}.'"';
    }

    $str =~ s/\,$//;
    $str .= '}';

    #write data to file
    my $time = sprintf("%02d/%02d/%04d %02d:%02d:%02d",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
    open(FH, ">", $params{graph_data}) or &printError("Cannot open file ".$params{graph_data}.": $!");
    flock(FH, LOCK_EX) or &printError("Cannot lock file ".$params{graph_data}.": $!");
    print FH "#Graph data\n";
    print FH "#[prinseq-".$WHAT."-$VERSION] [$time] Command: \"perl prinseq-".$WHAT.".pl".$command."\"\n" unless(exists $params{web});
    print FH $str;
    flock(FH, LOCK_UN) or &printError("Cannot unlock ".$params{graph_data}.": $!");
    close(FH);
    &printLog("Done with graph data");
    print STDERR "\tdone\n" if(exists $params{verbose});
}

&printWeb("STATUS: done");

##
#################################################################################
### MISC FUNCTIONS
#################################################################################
##

sub printError {
    my $msg = shift;
    print STDERR "\nERROR: ".$msg.".\n\nTry \'perl prinseq-".$WHAT.".pl -h\' for more information.\nExit program.\n";
    &printLog("ERROR: ".$msg.". Exit program.\n");
    exit(0);
}

sub printWarning {
    my $msg = shift;
    print STDERR "WARNING: ".$msg.".\n";
    &printLog("WARNING: ".$msg.".\n");
}

sub printWeb {
    my $msg = shift;
    if(exists $params{web}) {
        print STDERR &getTime()."$msg\n";
    }
}

sub getTime {
    return sprintf("[%02d/%02d/%04d %02d:%02d:%02d] ",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
}

sub printLog {
    my $msg = shift;
    if(exists $params{log}) {
        my $time = sprintf("%02d/%02d/%04d %02d:%02d:%02d",sub {($_[4]+1,$_[3],$_[5]+1900,$_[2],$_[1],$_[0])}->(localtime));
        open(FH, ">>", $params{log}) or &printError("Cannot open file ".$params{log}.": $!");
        flock(FH, LOCK_EX) or &printError("Cannot lock file ".$params{log}.": $!");
        print FH "[prinseq-".$WHAT."-$VERSION] [$time] $msg\n";
        flock(FH, LOCK_UN) or &printError("Cannot unlock ".$params{log}.": $!");
        close(FH);
    }
}

sub addCommas {
    my $num = shift;
    return unless(defined $num);
    return $num if($num < 1000);
    $num = scalar reverse $num;
    $num =~ s/(\d{3})/$1\,/g;
    $num =~ s/\,$//;
    $num = scalar reverse $num;
    return $num;
}

sub getLineNumber {
    my $file = shift;
    my $lines = 0;
    open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file |") or &printError("Could not open file $file: $!");
    $lines += tr/\n/\n/ while sysread(FILE, $_, 2 ** 16);
    close(FILE);
    return $lines;
}

sub readParamsFile {
    my $file = shift;
    my @args;
    my %parameters = ();
    open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file |") or &printError("Could not open file $file: $!");
    while(<FILE>) {
        next if(/^\#/);
        chomp();
        @args = split(/\s+/);
        if(@args) {
            $args[0] =~ s/^\-//;
            $parameters{$args[0]} = (defined $args[1] ? join(" ",@args[1..scalar(@args)-1]) : '');
        }
    }
    close(FILE);
    return \%parameters;
}

sub checkFileFormat {
    my $file = shift;

    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file |") or &printError("Could not open file $file: $!");
    while (<FILE>) {
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/o) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/o) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o))) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/o) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/o) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/o) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o))) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/o) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/o);
        }
    }
    close(FILE);
    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}

sub checkSlashnum {
    my $file = shift;
    open(FILE,"perl -pe 's/\r\n|\r/\n/g' < $file |") or &printError("Could not open file $file: $!");
    while(<FILE>) {
        chomp();
        next unless(length($_));
        if(/^\S+\/[12]\s*/ || /^\S+\_[LR]\s*/) {
            close(FILE);
            return (1,2,2);
        } elsif(/^\S+\_left\s*/) {
            close(FILE);
            return (1,6,5);
        } elsif(/^\S+\_right\s*/) {
            close(FILE);
            return (1,5,6);
        } else {
            close(FILE);
            return 0;
        }
    }
    return 0;
}

sub checkInputFormat {
    my ($format,$count,$id,$fasta,$fastq,$qual);
    $count = 3;
    $fasta = $fastq = $qual = 0;
    $format = 'unknown';

    while (<STDIN>) {
        push(@dataread,$_);
#        chomp();
 #       next unless(length($_));
        if($count-- == 0) {
            last;
        } elsif(!$fasta && /^\>\S+\s*/o) {
            $fasta = 1;
            $qual = 1;
        } elsif($fasta == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/o) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o))) {
            $fasta = 2;
        } elsif($qual == 1 && /^\s*\d+/) {
            $qual = 2;
        } elsif(!$fastq && /^\@(\S+)\s*/) {
            $id = $1;
            $fastq = 1;
        } elsif($fastq == 1 && (($aa && /^[ABCDEFGHIKLMNOPQRSTUVWYZXabcdefghiklmmopqrstuvwyzx*-]+/o) || (!$aa && /^[ACGTURYKMSWBDHVNXacgturykmswbdhvnx-]+/o))) {
            $fastq = 2;
        } elsif($fastq == 2 && /^\+(\S*)\s*/o) {
            $fastq = 3 if($id eq $1 || /^\+\s*$/o);
        }
    }

    if($fasta == 2) {
        $format = 'fasta';
    } elsif($qual == 2) {
        $format = 'qual';
    } elsif($fastq == 3) {
        $format = 'fastq';
    }

    return $format;
}

sub getArrayMean {
    return @_ ? sum(@_) / @_ : 0;
}

sub convertQualNumsToAscii {
    my $qual = shift;
    my @ascii;
    $qual =~ s/^\s+//;
    $qual =~ s/\s+$//;
    my @nums = split(/\s+/,$qual);
    foreach(@nums) {
	push(@ascii,chr(($_ <= 93 ? $_ : 93) + 33));
    }
    return \@ascii;
}

sub convertQualNumsToAsciiString {
    my $qual = shift;
    my $ascii;
    $qual =~ s/^\s+//;
    $qual =~ s/\s+$//;
    my @nums = split(/\s+/,$qual);
    foreach(@nums) {
	$ascii .= chr(($_ <= 93 ? $_ : 93) + 33);
    }
    return $ascii;
}

sub convertQualAsciiToNums {
    my $qual = shift;
    my @nums;
    my @ascii = split(//,$qual);
    foreach(@ascii) {
	push(@nums,(ord($_) - 33));
    }
    return \@nums;
}

sub convertQualAsciiToNumsPhred64 {
    my $qual = shift;
    my @nums;
    my $tmp;
    my $err = 0;
    my @ascii = split('',$qual);
    foreach(@ascii) {
        $tmp = (ord($_) - 64);
        if($tmp < 0) {
	    $err = 1;
	    last;
	}
	push(@nums,$tmp);
    }
    return (\@nums,$err);
}

sub convertQualArrayToString {
    my ($nums,$linelen) = @_;
    $linelen = 80 unless($linelen);
    my $str;
    my $count = 1;
    foreach my $n (@$nums) {
        $str .= ($n < 10 ? ' '.$n : $n).' ';
        if(++$count > $linelen) {
            $count = 1;
            $str =~ s/\s$//;
            $str .= "\n";
        }
    }
    $str =~ s/[\s\n]$//;
    return $str;
}

sub checkRange {
    my ($range,$val) = @_;
    my @ranges = split(/\,/,$range);
    foreach my $r (@ranges) {
	my @tmp = split(/\-/,$r);
	return 0 if($val < $tmp[0] || $val > $tmp[1]);
    }
    return 1;
}

sub getFiltercounts {
    my @order = qw(seq_num trim_left trim_right trim_left_p trim_right_p trim_qual_left trim_qual_right trim_tail_left trim_tail_right trim_ns_left trim_ns_right zero_length min_len max_len range_len min_qual_score max_qual_score min_qual_mean max_qual_mean min_gc max_gc range_gc ns_max_p ns_max_n noniupac custom_params lc_method derep);
    my @counts;
    foreach my $p (@order) {
        if(exists $filtercount{$p}) {
            push(@counts,$p.': '.$filtercount{$p});
        }
    }
    return \@counts;
}

sub readEntryFastq {
    my ($fh,$skip,$trimnum) = @_;
    if($skip > 0) {
        return 1;
    }
    my ($seq,$seqid,$qual,$header,$line,$tmpid);
    if(defined($line = <$fh>)) {
        $line =~ /^\@(\S+)\s*(.*)$/o;
        $seqid = $1;
        $header = $2 || '';
        $seq = readline($fh);
        chomp($seq);
        readline($fh);
        $qual = readline($fh);
        chomp($qual);
        #progress bar stuff
        $counter += 4;
    }
    if($slashnum && $seqid) {
        $tmpid = substr($seqid, 0, -$trimnum)
    } else {
        $tmpid = $seqid;
    }
    return ($seq,$seqid,$qual,$header,$tmpid);
}

sub readEntryFasta {
    my ($fh,$skip,$trimnum,$nextid,$nextheader) = @_;
    if($skip > 0) {
        return (1,undef,undef,undef,$nextid,$nextheader);
    }
    my ($seq,$seqid,$header,$line,$tmpid);
    $seq = '';
    if($nextid) {
        $seqid = $nextid;
        $header = $nextheader;
    }
    while(<$fh>) {
        chomp();
        if(/^>(\S+)\s*(.*)$/o) {
            if($seq) {
                $nextid = $1;
                $nextheader = $2 || '';
                last;
            }
            $seqid = $1;
            $header = $2 || '';
        } else {
            $seq .= $_;
        }
    }
    if($slashnum) {
        $tmpid = substr($seqid, 0, -$trimnum)
    } else {
        $tmpid = $seqid;
    }
    return ($seq,$seqid,$header,$tmpid,$nextid,$nextheader);
}

sub processEntry {
    my ($length,$seq,$seqid,$qual,$header) = @_;
    #check that sequence and quality are same length
    if(defined $qual && length($qual) && $length != length($qual)) {
        &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"");
    }
    #remove anything non-alphabetic from sequences
    $seq =~ tr/a-zA-Z//cd;
    $numseqs++;
    $numbases += $length;
    #process entry
    if($exists_stats) { #calc summary stats
        $seq = uc($seq);
        &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
        if(exists $params{stats_dupl}) {
            push(@seqs,[$seq,$numseqs,$length]);
        }
    } else { #process data
        &processData($seqid,$seq,$qual,$header) unless($webnoprocess);
        #get graph data
        if($exists_graphdata) {
            $ucseq = uc($seq);
            &getSeqStats(\%graphdata,$ucseq,$length);
            if($qual) {
                &getQualStats(\%graphdata,$qual,$length);
            }
            if($webstats{de} || $graphstats{de}) {
                $md5 = md5_hex($ucseq);
                if(exists $md5sg{$md5}) { #forward duplicate
                    $graphdata{dubslength}->{$length}->{0}++;
                    $md5sg{$md5}->{0}++;
                } else {
                    $md5r = md5_hex(&revcompuc($ucseq));
                    if(exists $md5sg{$md5r}) { #reverse duplicate
                        $graphdata{dubslength}->{$length}->{3}++;
                        $md5sg{$md5}->{3}++;
                    }
                    unless(exists $md5sg{$md5}) {
                        $md5sg{$md5} = {0 => 0, 3 => 0};
                    }
                }
            } elsif($webstats{da} || $graphstats{da}) {
                push(@seqs,[$ucseq,$numseqs,$length]);
            }
        }
    }
}

sub processEntryPairedEnd {
    my ($length,$seq,$seqid,$header,$qual,$length2,$seq2,$seqid2,$header2,$qual2) = @_;
    my ($tmpseq1,$tmpseq2,$tmpgood1,$tmpgood2,$tmpbegin1,$tmpend1,$tmpbegin2,$tmpend2);
    $tmpgood1 = $tmpgood2 = 0;
    if($length) {
        #check that sequence and quality are same length
        if(defined $qual && $length != length($qual)) {
            &printError("The number of bases and quality scores are not the same for sequence \"$seqid\"");
        }
        #remove anything non-alphabetic from sequences
        $seq =~ tr/a-zA-Z//cd;
        $numseqs++;
        $numbases += $length;
        #process entry
        if($exists_stats) { #calc summary stats
            $seq = uc($seq);
            &calcSeqStats($seq,$length,\%stats,\%kmers,\%odds,\%counts);
        } else { #process data
            ($tmpseq1,$tmpgood1,$tmpbegin1,$tmpend1) = &processData($seqid,$seq,$qual,$header) unless($webnoprocess);
            #get graph data
            if($exists_graphdata) {
                $ucseq = uc($seq);
                &getSeqStats(\%graphdata,$ucseq,$length);
                if($qual) {
                    &getQualStats(\%graphdata,$qual,$length);
                }
                if($length2) {
                    $pairs++;
                }
                if($webstats{de} || $graphstats{de}) {
                    if($length2) { #both
                        $md5 = md5_hex($ucseq.'0'.uc($seq2));
                        if(exists $md5sg{$md5}) { #forward duplicate
                            $md5sg{$md5}->{0}++;
                            $graphdata{dubslength}->{($length+$length2)}->{0}++;
                        } else {
                            $md5sg{$md5}->{0} = 0;
                        }
                        $md5 = md5_hex($ucseq);
                        unless(exists $md5sg{$md5}) {
                            $md5sg{$md5}->{0} = 0;
                        }
                    } else { #only seq, no seq2
                        $md5 = md5_hex($ucseq);
                        if(exists $md5sg{$md5}) { #forward duplicate
                            $md5sg{$md5}->{0}++;
                            $graphdata{dubslength}->{$length}->{0}++;
                        } else {
                            $md5sg{$md5}->{0} = 0;
                        }
                    }
                }
            }
        }
    }
    if($length2) {
        if(defined $qual2 && $length2 != length($qual2)) {
            &printError("The number of bases and quality scores are not the same for sequence \"$seqid2\"");
        }
        #remove anything non-alphabetic from sequences
        $seq2 =~ tr/a-zA-Z//cd;
        $numseqs2++;
        $numbases2 += $length2;
        #process entry
        if($exists_stats) { #calc summary stats
            $seq2 = uc($seq2);
            &calcSeqStats($seq2,$length2,\%stats,\%kmers2,\%odds,\%counts2,1); #last one used to calculate stats for second input file
        } else { #process data
            ($tmpseq2,$tmpgood2,$tmpbegin2,$tmpend2) = &processData($seqid2,$seq2,$qual2,$header2) unless($webnoprocess);
            #get graph data
            if($exists_graphdata) {
                $ucseq = uc($seq2);
                &getSeqStats(\%graphdata2,$ucseq,$length2);
                if($qual2) {
                    &getQualStats(\%graphdata2,$qual2,$length2);
                }
                if($webstats{de} || $graphstats{de}) {
                    #only seq2, no seq
                    $md5 = md5_hex($ucseq);
                    unless($length) {
                        if(exists $md5sg{$md5}) { #forward duplicate
                            $md5sg{$md5}->{0}++;
                            $graphdata{dubslength}->{$length2}->{0}++;
                        } else {
                            $md5sg{$md5}->{0} = 0;
                        }
                    }
                }
            }
        }
    }

    return if($exists_stats || $webnoprocess);

    #check for duplicates
    if($derep) {
        if($tmpgood1 && $tmpgood2) {
            $md5 = md5_hex($tmpseq1.'0'.$tmpseq2);
            if(exists $md5s{$md5}) { #forward duplicate
                $md5s{$md5}++;
                if($derepmin <= $md5s{$md5}+1) {
                    $tmpgood1 = $tmpgood2 = 0;
                    $filtercount{derep}++;
                }
            } else {
                $md5s{$md5} = 0;
            }
            $md5 = md5_hex($tmpseq1);
            unless(exists $md5s{$md5}) {
                $md5s{$md5} = 0;
            }
            $md5 = md5_hex($tmpseq2);
            unless(exists $md5s{$md5}) {
                $md5s{$md5} = 0;
            }
        } elsif($tmpgood1) {
            $md5 = md5_hex($tmpseq1);
            if(exists $md5s{$md5}) { #forward duplicate
                $md5s{$md5}++;
                if($derepmin <= $md5s{$md5}+1) {
                    $tmpgood1 = 0;
                    $filtercount{derep}++;
                }
            } else {
                $md5s{$md5} = 0;
            }
        } elsif($tmpgood2) {
            $md5 = md5_hex($tmpseq2);
            if(exists $md5s{$md5}) { #forward duplicate
                $md5s{$md5}++;
                if($derepmin <= $md5s{$md5}+1) {
                    $tmpgood2 = 0;
                    $filtercount{derep}++;
                }
            } else {
                $md5s{$md5} = 0;
            }
        }
    }

    #write to outputs files
    if($tmpgood1) {
        if(exists $params{seq_id}) {
            if($mappings) {
                print $fhmappings join("\t",$seqid,$params{seq_id}.$seqcount)."\n";
            }
            $seqid = $params{seq_id}.$seqcount;
        }
        if(exists $params{rm_header}) {
            $header = undef;
        }
        #trim if necessary
        if($tmpbegin1) {
            $seq = substr($seq,$tmpbegin1);
            $qual = substr($qual,$tmpbegin1) if(defined $qual && length($qual));
        }
        if($tmpend1) {
            $length = length($seq);
            $seq = substr($seq,0,$length-$tmpend1);
            $qual = substr($qual,0,$length-$tmpend1) if(defined $qual && length($qual));
        }
        #change case
        if(exists $params{seq_case}) {
            if($params{seq_case} eq 'lower') { #lower case
                $seq = lc($seq);
            } elsif($params{seq_case} eq 'upper') { #upper case
                $seq = uc($seq);
            }
        }
        #convert between DNA and RNA
        if(exists $params{dna_rna}) {
            if($params{dna_rna} eq 'dna') { #RNA to DNA
                $seq =~ tr/Uu/Tt/;
            } elsif($params{dna_rna} eq 'rna') { #DNA to RNA
                $seq =~ tr/Tt/Uu/;
            }
        }
    }
    if($tmpgood2) {
        if(exists $params{seq_id}) {
            if($mappings) {
                print $fhmappings join("\t",$seqid2,$params{seq_id}.$seqcount)."\n";
            }
            $seqid2 = $params{seq_id}.$seqcount;
        }
        if(exists $params{rm_header}) {
            $header2 = undef;
        }
        #trim if necessary
        if($tmpbegin2) {
            $seq2 = substr($seq2,$tmpbegin2);
            $qual2 = substr($qual2,$tmpbegin2) if(defined $qual2 && length($qual2));
        }
        if($tmpend2) {
            $length2 = length($seq2);
            $seq2 = substr($seq2,0,$length2-$tmpend2);
            $qual2 = substr($qual2,0,$length2-$tmpend2) if(defined $qual2 && length($qual2));
        }
        #change case
        if(exists $params{seq_case}) {
            if($params{seq_case} eq 'lower') { #lower case
                $seq2 = lc($seq2);
            } elsif($params{seq_case} eq 'upper') { #upper case
                $seq2 = uc($seq2);
            }
        }
        #convert between DNA and RNA
        if(exists $params{dna_rna}) {
            if($params{dna_rna} eq 'dna') { #RNA to DNA
                $seq2 =~ tr/Uu/Tt/;
            } elsif($params{dna_rna} eq 'rna') { #DNA to RNA
                $seq2 =~ tr/Tt/Uu/;
            }
        }
    }
    if($tmpgood1 && $tmpgood2) { #pair
        $seqcount++;
        $seqbases += $length+$length2;
        return if($nogood);
        if($params{out_format} == 3) { # FASTQ
            &printError("missing quality data for sequence \"$seqid\" or greater number of sequences than available quality scores") unless(defined $qual);
            &printError("missing quality data for sequence \"$seqid2\" or greater number of sequences than available quality scores") unless(defined $qual2);
            if($stdoutgood) {
                print STDOUT '@'.$seqid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                print STDOUT $qual."\n";
                print STDOUT '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print STDOUT $seq2."\n";
                print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                print STDOUT $qual2."\n";
            } else {
                print $fhgood '@'.$seqid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
                print $fhgood '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                print $fhgood $qual."\n";
                print $fh2good '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print $fh2good $seq2."\n";
                print $fh2good '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                print $fh2good $qual2."\n";
            }
        } else { #FASTA
            #set line length
            if($linelen) {
                $seq =~ s/(.{$linelen})/$1\n/g;
                $seq =~ s/\n$//;
                $seq2 =~ s/(.{$linelen})/$1\n/g;
                $seq2 =~ s/\n$//;
            }
            if($stdoutgood) {
                print STDOUT '>'.$seqid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print STDOUT $seq2."\n";
            } else {
                print $fhgood '>'.$seqid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
                print $fh2good '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print $fh2good $seq2."\n";
            }
        }
    } elsif($tmpgood1) { #singleton
        $seqcount1++;
        $seqbases1 += $length;
        return if($nogood);
        if($params{out_format} == 3) { # FASTQ
            &printError("missing quality data for sequence \"$seqid\" or greater number of sequences than available quality scores") unless(defined $qual);
            if($stdoutgood) {
                print STDOUT '@'.$seqid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                print STDOUT $qual."\n";
            } else {
                print $fhgood2 '@'.$seqid.($header ? ' '.$header : '')."\n";
                print $fhgood2 $seq."\n";
                print $fhgood2 '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                print $fhgood2 $qual."\n";
            }
        } else { #FASTA
            #set line length
            if($linelen) {
                $seq =~ s/(.{$linelen})/$1\n/g;
                $seq =~ s/\n$//;
            }
            if($stdoutgood) {
                print STDOUT '>'.$seqid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
            } else {
                print $fhgood2 '>'.$seqid.($header ? ' '.$header : '')."\n";
                print $fhgood2 $seq."\n";
            }
        }
        if($length2) {
            $badcount2++;
            $badbases2 += length($seq2);
            return if($nobad);
            #write data
            if($params{out_format} == 3) { # FASTQ
                &printError("missing quality data for sequence \"$seqid2\" or greater number of sequences than available quality scores") unless(defined $qual2);
                if($stdoutbad) {
                    print STDOUT '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print STDOUT $seq2."\n";
                    print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                    print STDOUT $qual2."\n";
                } else {
                    print $fh2bad '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print $fh2bad $seq2."\n";
                    print $fh2bad '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                    print $fh2bad $qual2."\n";
                }
            } else { #FASTA
                #set line length
                if($linelen) {
                    $seq2 =~ s/(.{$linelen})/$1\n/g;
                    $seq2 =~ s/\n$//;
                }
                if($stdoutbad) {
                    print STDOUT '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print STDOUT $seq2."\n";
                } else {
                    print $fh2bad '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print $fh2bad $seq2."\n";
                }
            }
        }
    } elsif($tmpgood2) { #singleton
        $seqcount2++;
        $seqbases2 += $length2;
        return if($nogood);
        if($params{out_format} == 3) { # FASTQ
            &printError("missing quality data for sequence \"$seqid2\" or greater number of sequences than available quality scores") unless(defined $qual2);
            if($stdoutgood) {
                print STDOUT '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print STDOUT $seq2."\n";
                print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                print STDOUT $qual2."\n";
            } else {
                print $fh2good2 '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print $fh2good2 $seq2."\n";
                print $fh2good2 '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                print $fh2good2 $qual2."\n";
            }
        } else { #FASTA
            #set line length
            if($linelen) {
                $seq2 =~ s/(.{$linelen})/$1\n/g;
                $seq2 =~ s/\n$//;
            }
            if($stdoutgood) {
                print STDOUT '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print STDOUT $seq2."\n";
            } else {
                print $fh2good2 '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                print $fh2good2 $seq2."\n";
            }
        }
        if($length) {
            $badcount++;
            $badbases += length($seq);
            return if($nobad);
            #write data
            if($params{out_format} == 3) { # FASTQ
                &printError("missing quality data for sequence \"$seqid\" or greater number of sequences than available quality scores") unless(defined $qual);
                if($stdoutbad) {
                    print STDOUT '@'.$seqid.($header ? ' '.$header : '')."\n";
                    print STDOUT $seq."\n";
                    print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                    print STDOUT $qual."\n";
                } else {
                    print $fhbad '@'.$seqid.($header ? ' '.$header : '')."\n";
                    print $fhbad $seq."\n";
                    print $fhbad '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                    print $fhbad $qual."\n";
                }
            } else { #FASTA
                #set line length
                if($linelen) {
                    $seq =~ s/(.{$linelen})/$1\n/g;
                    $seq =~ s/\n$//;
                }
                if($stdoutbad) {
                    print STDOUT '>'.$seqid.($header ? ' '.$header : '')."\n";
                    print STDOUT $seq."\n";
                } else {
                    print $fhbad '>'.$seqid.($header ? ' '.$header : '')."\n";
                    print $fhbad $seq."\n";
                }
            }
        }
    } else {
        if($length) {
            $badcount++;
            $badbases += length($seq);
            return if($nobad);
            #write data
            if($params{out_format} == 3) { # FASTQ
                &printError("missing quality data for sequence \"$seqid\" or greater number of sequences than available quality scores") unless(defined $qual);
                if($stdoutbad) {
                    print STDOUT '@'.$seqid.($header ? ' '.$header : '')."\n";
                    print STDOUT $seq."\n";
                    print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                    print STDOUT $qual."\n";
                } else {
                    print $fhbad '@'.$seqid.($header ? ' '.$header : '')."\n";
                    print $fhbad $seq."\n";
                    print $fhbad '+'.(exists $params{no_qual_header} ? '' : $seqid.($header ? ' '.$header : ''))."\n";
                    print $fhbad $qual."\n";
                }
            } else { #FASTA
                #set line length
                if($linelen) {
                    $seq =~ s/(.{$linelen})/$1\n/g;
                    $seq =~ s/\n$//;
                }
                if($stdoutbad) {
                    print STDOUT '>'.$seqid.($header ? ' '.$header : '')."\n";
                    print STDOUT $seq."\n";
                } else {
                    print $fhbad '>'.$seqid.($header ? ' '.$header : '')."\n";
                    print $fhbad $seq."\n";
                }
            }
        }
        if($length2) {
            $badcount2++;
            $badbases2 += length($seq2);
            return if($nobad);
            #write data
            if($params{out_format} == 3) { # FASTQ
                &printError("missing quality data for sequence \"$seqid2\" or greater number of sequences than available quality scores") unless(defined $qual2);
                if($stdoutbad) {
                    print STDOUT '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print STDOUT $seq2."\n";
                    print STDOUT '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                    print STDOUT $qual2."\n";
                } else {
                    print $fh2bad '@'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print $fh2bad $seq2."\n";
                    print $fh2bad '+'.(exists $params{no_qual_header} ? '' : $seqid2.($header2 ? ' '.$header2 : ''))."\n";
                    print $fh2bad $qual2."\n";
                }
            } else { #FASTA
                #set line length
                if($linelen) {
                    $seq2 =~ s/(.{$linelen})/$1\n/g;
                    $seq2 =~ s/\n$//;
                }
                if($stdoutbad) {
                    print STDOUT '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print STDOUT $seq2."\n";
                } else {
                    print $fh2bad '>'.$seqid2.($header2 ? ' '.$header2 : '')."\n";
                    print $fh2bad $seq2."\n";
                }
            }
        }
    }
}

#process sequence (and qual) data
sub processData {
    my ($sid,$seq,$qual,$header) = @_;
    #assume sequence is good ;-)
    my $good = 1;
    my $seqn = uc($seq);
    my $qualn = $qual;
    my $begin = 0;
    my $end = 0;
    my ($length,$bylength,$qualsnums);

    #check for maximum number sequences requested
    if(exists $params{seq_num} && $params{seq_num} <= $seqcount) {
        $good = 0;
        $filtercount{seq_num}++;
    }

    #trim sequence ends
    if($good && exists $params{trim_left}) {
        $begin += $params{trim_left};
        if($begin >= length($seqn)) {
            $good = 0;
            $filtercount{trim_left}++;
        } else {
            $seqn = substr($seqn,$begin);
            $qualn = substr($qualn,$begin) if(defined $qualn && length($qualn));
        }
    }
    if($good && exists $params{trim_right}) {
        $end += $params{trim_right};
        $length = length($seqn);
        if($end >= $length) {
            $good = 0;
            $filtercount{trim_right}++;
        } else {
            $seqn = substr($seqn,0,$length-$end);
            $qualn = substr($qualn,0,$length-$end) if(defined $qualn && length($qualn));
        }
    }
    if($good && exists $params{trim_left_p}) {
        $length = length($seqn);
        my $begintmp = int($params{trim_left_p}/100*$length);
        if($begintmp >= $length) {
            $good = 0;
            $filtercount{trim_left_p}++;
        } else {
            $seqn = substr($seqn,$begintmp);
            $qualn = substr($qualn,$begintmp) if(defined $qualn && length($qualn));
            $begin += $begintmp;
        }
    }
    if($good && exists $params{trim_right_p}) {
        $length = length($seqn);
        my $endtmp = int($params{trim_right_p}/100*$length);
        if($endtmp >= $length) {
            $good = 0;
            $filtercount{trim_right_p}++;
        } else {
            $seqn = substr($seqn,0,$length-$endtmp);
            $qualn = substr($qualn,0,$length-$endtmp) if(defined $qualn && length($qualn));
            $end += $endtmp;
        }
    }

    #check for quality scores
    if($good && defined $qualn && $trimscore) {
        my ($err);
        if(exists $params{phred64}) { #scale data to Phred scale if necessary
            ($qualsnums,$err) = &convertQualAsciiToNumsPhred64($qualn);
            if($err) {
                &printError("The sequence quality scores are not in Phred+64 format");
            }
        } else {
            $qualsnums = &convertQualAsciiToNums($qualn);
        }
        $length = length($qualn);
        my $i = 0;
        my $begintmp = 0;
        my $endtmp = 0;
        my ($window,$val);
        #left
        if(exists $params{trim_qual_left}) {
            while($i < $length) {
                #calculate maximum window
                $window = ($i+$params{trim_qual_window} <= $length ? $params{trim_qual_window} : ($length-$i));
                #calculate value used to compare with given value
                if($window == 1) {
                    $val = $qualsnums->[$i];
                } elsif($params{trim_qual_type} eq 'min') {
                    $val = min(@$qualsnums[$i..($i+$window-1)]);
                } elsif($params{trim_qual_type} eq 'max') {
                    $val = max(@$qualsnums[$i..($i+$window-1)]);
                } elsif($params{trim_qual_type} eq 'mean') {
                    $val = &getArrayMean(@$qualsnums[$i..($i+$window-1)]);
                } elsif($params{trim_qual_type} eq 'sum') {
                    last if($window < $params{trim_qual_window});
                    $val = sum(@$qualsnums[$i..($i+$window-1)]);
                } else {
                    last;
                }
                #compare values
                if(($params{trim_qual_rule} eq 'lt' && $val < $params{trim_qual_left}) || ($params{trim_qual_rule} eq 'gt' && $val > $params{trim_qual_left}) || ($params{trim_qual_rule} eq 'et' && $val == $params{trim_qual_left})) {
                    $begintmp += $params{trim_qual_step};
                    $i += $params{trim_qual_step};
                } else {
                    last;
                }
            }
            if($begintmp >= $length) {
                $good = 0;
                $filtercount{trim_qual_left}++;
            } elsif($begintmp > 0) {
                $seqn = substr($seqn,$begintmp);
                $qualn = substr($qualn,$begintmp);
                $begin += $begintmp;
            }
        }
        #right
        if($good && exists $params{trim_qual_right}) {
            $length -= $begintmp;
            my @quals = reverse(@$qualsnums);
            $i = 0;
            while($i < $length) {
                #calculate maximum window
                $window = ($i+$params{trim_qual_window} <= $length ? $params{trim_qual_window} : ($length-$i));
                #calculate value used to compare with given value
                if($window == 1) {
                    $val = $quals[$i];
                } elsif($params{trim_qual_type} eq 'min') {
                    $val = min(@quals[$i..($i+$window-1)]);
                } elsif($params{trim_qual_type} eq 'max') {
                    $val = max(@quals[$i..($i+$window-1)]);
                } elsif($params{trim_qual_type} eq 'mean') {
                    $val = &getArrayMean(@quals[$i..($i+$window-1)]);
                } elsif($params{trim_qual_type} eq 'sum') {
                    last if($window < $params{trim_qual_window});
                    $val = sum(@quals[$i..($i+$window-1)]);
                } else {
                    last;
                }
                #compare values
                if(($params{trim_qual_rule} eq 'lt' && $val < $params{trim_qual_right}) || ($params{trim_qual_rule} eq 'gt' && $val > $params{trim_qual_right}) || ($params{trim_qual_rule} eq 'et' && $val == $params{trim_qual_right})) {
                    $endtmp += $params{trim_qual_step};
                    $i += $params{trim_qual_step};
                } else {
                    last;
                }
            }
            if($endtmp >= $length) {
                $good = 0;
                $filtercount{trim_qual_right}++;
            } elsif($endtmp > 0) {
                $seqn = substr($seqn,0,$length-$endtmp);
                $qualn = substr($qualn,0,$length-$endtmp);
                $end += $endtmp;
            }
        }
    }

    #check for tails with min trimtails char repeats
    if($good && exists $params{trim_tail_left}) {
        $length = length($seqn);
        my $begintmp = 0;
        if($seqn =~ $repAleft || $seqn =~ $repTleft) {
            my @tmp = split('',$seqn);
            my $tmpchar = $tmp[0]; #A or T
            $begintmp += $params{trim_tail_left};
            foreach ($params{trim_tail_left}..$length-1) {
                last unless($tmp[$_] eq $tmpchar || $tmp[$_] eq 'N');
                $begintmp++;
            }
            if($begintmp >= $length) {
                $good = 0;
                $filtercount{trim_tail_left}++;
            } else {
                $seqn = substr($seqn,$begintmp);
                $qualn = substr($qualn,$begintmp) if(defined $qualn && length($qualn));
                $length = length($seqn);
                $begin += $begintmp;
            }
        }
    }
    if($good && exists $params{trim_tail_right}) {
        $length = length($seqn);
        my $endtmp = 0;
        if($seqn =~ $repAright || $seqn =~ $repTright) {
            my @tmp = split('',$seqn);
            my $tmpchar = $tmp[$length-1]; #A or T
            $endtmp += $params{trim_tail_right};
            foreach (reverse 0..$length-$params{trim_tail_right}-1) {
                last unless($tmp[$_] eq $tmpchar || $tmp[$_] eq 'N');
                $endtmp++;
            }
            if($endtmp >= $length) {
                $good = 0;
                $filtercount{trim_tail_right}++;
            } else {
                $seqn = substr($seqn,0,$length-$endtmp);
                $qualn = substr($qualn,0,$length-$endtmp) if(defined $qualn && length($qualn));
                $end += $endtmp;
            }
        }
    }
    if($good && exists $params{trim_ns_left}) {
        $length = length($seqn);
        my $begintmp = 0;
        if($seqn =~ $repNleft) {
            my @tmp = split('',$seqn);
            $begintmp += $params{trim_ns_left};
            foreach ($params{trim_ns_left}..$length-1) {
                last unless($tmp[$_] eq 'N');
                $begintmp++;
            }
            if($begintmp >= $length) {
                $good = 0;
                $filtercount{trim_ns_left}++;
            } else {
                $seqn = substr($seqn,$begintmp);
                $qualn = substr($qualn,$begintmp) if(defined $qualn && length($qualn));
                $length = length($seqn);
                $begin += $begintmp;
            }
        }
    }
    if($good && exists $params{trim_ns_right}) {
        $length = length($seqn);
        my $endtmp = 0;
        if($seqn =~ $repNright) {
            my @tmp = split('',$seqn);
            $endtmp += $params{trim_ns_right};
            foreach (reverse 0..$length-$params{trim_ns_right}-1) {
                last unless($tmp[$_] eq 'N');
                $endtmp++;
            }
            if($endtmp >= $length) {
                $good = 0;
                $filtercount{trim_ns_right}++;
            } else {
                $seqn = substr($seqn,0,$length-$endtmp);
                $qualn = substr($qualn,0,$length-$endtmp) if(defined $qualn && length($qualn));
                $end += $endtmp;
            }
        }
    }

    #check if trim to certain length
    $length = length($seqn);
    if($good && exists $params{trim_to_len} && $length > $params{trim_to_len}) {
        $seqn = substr($seqn,0,$params{trim_to_len});
        $qualn = substr($qualn,0,$params{trim_to_len}) if(defined $qualn && length($qualn));
        $end += ($length-$params{trim_to_len});
    }

    #check for sequence length
    $length = length($seqn);
    $bylength = ($length ? 100/$length : 0);
    if($bylength == 0) {
        $good = 0;
        $filtercount{zero_length}++;
    }
    if($good && exists $params{min_len} && $length < $params{min_len}) {
        $good = 0;
        $filtercount{min_len}++;
    }
    if($good && exists $params{max_len} && $length > $params{max_len}) {
        $good = 0;
        $filtercount{max_len}++;
    }
    if($good && exists $params{range_len} && !&checkRange($params{range_len},$length)) {
        $good = 0;
        $filtercount{range_len}++;
    }

    #check for quality scores
    if($good && defined $qualn && (exists $params{min_qual_score} || exists $params{max_qual_score} || exists $params{min_qual_mean} || exists $params{max_qual_mean})) {
        my ($err);
        if($qualsnums) {
            if($begin > 0) {
                shift(@$qualsnums) foreach(1..$begin);
            }
            if($end > 0) {
                pop(@$qualsnums) foreach(1..$end);
            }
        } else {
            if(exists $params{phred64}) { #scale data to Phred scale if necessary
                ($qualsnums,$err) = &convertQualAsciiToNumsPhred64($qualn);
                if($err) {
                    &printError("The sequence quality scores are not in Phred+64 format");
                }
            } else {
                $qualsnums = &convertQualAsciiToNums($qualn);
            }
        }
        if($good && exists $params{min_qual_score} && min(@$qualsnums) < $params{min_qual_score}) {
            $good = 0;
            $filtercount{min_qual_score}++;
        }
        if($good && exists $params{max_qual_score} && max(@$qualsnums) < $params{max_qual_score}) {
            $good = 0;
            $filtercount{max_qual_score}++;
        }
        if($good && exists $params{min_qual_mean} && &getArrayMean(@$qualsnums) < $params{min_qual_mean}) {
            $good = 0;
            $filtercount{min_qual_mean}++;
        }
        if($good && exists $params{max_qual_mean} && &getArrayMean(@$qualsnums) < $params{max_qual_mean}) {
            $good = 0;
            $filtercount{max_qual_mean}++;
        }
    }

    #check for GC content
    if($good && (exists $params{min_gc} || exists $params{max_gc} || exists $params{range_gc})) {
        my $gc = ($seqn =~ tr/GC//);
        $gc = sprintf("%d",$gc*$bylength);
        if(exists $params{min_gc} && $gc < $params{min_gc}) {
            $good = 0;
            $filtercount{min_gc}++;
        }
        if($good && exists $params{max_gc} && $gc > $params{max_gc}) {
            $good = 0;
            $filtercount{max_gc}++;
        }
        if($good && exists $params{range_gc} && !&checkRange($params{range_gc},$gc)) {
            $good = 0;
            $filtercount{range_gc}++;
        }
    }

    #check for N's in sequence
    if($good && (exists $params{ns_max_p} || exists $params{ns_max_n})) {
        my $ns = ($seqn =~ tr/N//);
        if(exists $params{ns_max_p} && ($ns*$bylength) > $params{ns_max_p}) {
            $good = 0;
            $filtercount{ns_max_p}++;
        }
        if($good && exists $params{ns_max_n} && $ns > $params{ns_max_n}) {
            $good = 0;
            $filtercount{ns_max_n}++;
        }
    }

    #check for non IUPAC chars in sequence
    if($good && exists $params{noniupac} && $seqn =~ /[^ACGTN]/o) {
        $good = 0;
        $filtercount{noniupac}++;
    }

    #check for additional filter parameters
    if($good && @cps) {
        foreach my $p (@cps) {
            if($p->[0]) { #repeats
                if(index($seqn,$p->[1]x$p->[2]) != -1) {
                    $good = 0;
                    $filtercount{custom_params}++;
                    last;
                }
            } else { #percentage
                my $ns = 0;
                my $v = $p->[1];
                $ns++ while($seqn =~ /$v/g);
                if((100*$ns/$length) > $p->[2]) {
                    $good = 0;
                    $filtercount{custom_params}++;
                    last;
                }
            }
        }
    }

    #check for sequence complexity
    if($good && defined $complval) {
        my ($rest,$steps,@vals,$str,$num,$bynum);
        if($length <= $WINDOWSIZE) {
            $rest = $length;
            $steps = 0;
        } else {
            $steps = int(($length - $WINDOWSIZE) / $WINDOWSTEP) + 1;
            $rest = $length - $steps * $WINDOWSTEP;
            unless($rest > $WINDOWSTEP) {
                $rest += $WINDOWSTEP;
                $steps--;
            }
        }
        $num = $WINDOWSIZE-2;
        $bynum = 1/$num;
        $num--;
        my $mean = 0;
        if($params{lc_method} eq 'dust') {
            my $dustscore;
            foreach my $i (0..$steps-1) {
                $str = substr($seqn,($i * $WINDOWSTEP),$WINDOWSIZE);
                %counts = ();
                foreach my $i (@WINDOWSIZEARRAY) {
                    $counts{substr($str,$i,3)}++;
                }
                $dustscore = 0;
                foreach(values %counts) {
                    $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
                }
                push(@vals,($dustscore * $bynum));
            }
            #last step
            if($rest > 5) {
                $str = substr($seqn,($steps * $WINDOWSTEP),$rest);
                %counts = ();
                $num = $rest-2;
                foreach my $i (0..($num - 1)) {
                    $counts{substr($str,$i,3)}++;
                }
                $dustscore = 0;
                foreach(values %counts) {
                    $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
                }
                push(@vals,(($dustscore / ($num-1)) * (($WINDOWSIZE - 2) / $num)));
            } else {
                push(@vals,31); #to assign a maximum score based on the scaling factor 100/31
            }
            $mean = &getArrayMean(@vals);
            if(int($mean * 100 / 31) > $complval) {
                $good = 0;
                $filtercount{lc_method}++;
            }
        } elsif($params{lc_method} eq 'entropy') {
            my $entropyval;
            foreach my $i (0..$steps-1) {
                $str = substr($seqn,($i * $WINDOWSTEP),$WINDOWSIZE);
                %counts = ();
                foreach my $i (@WINDOWSIZEARRAY) {
                    $counts{substr($str,$i,3)}++;
                }
                $entropyval = 0;
                foreach(values %counts) {
                    $entropyval -= ($_ * $bynum) * log($_ * $bynum);
                }
                push(@vals,($entropyval * $ONEOVERLOG62));
            }
            #last step
            if($rest > 5) {
                $str = substr($seqn,($steps * $WINDOWSTEP),$rest);
                %counts = ();
                $num = $rest-2;
                foreach my $i (0..($num - 1)) {
                    $counts{substr($str,$i,3)}++;
                }
                $entropyval = 0;
                $bynum = 1/$num;
                foreach(values %counts) {
                    $entropyval -= ($_ * $bynum) * log($_ * $bynum);
                }
                push(@vals,($entropyval / log($num)));
            } else {
                push(@vals,0);
            }
            $mean = &getArrayMean(@vals);
            if(int($mean * 100) < $complval) {
                $good = 0;
                $filtercount{lc_method}++;
            }
        }
    }

    #stop here for paired-end reads
    if($file2) {
        return ($seqn,$good,$begin,$end);
    }

    #check for read duplicates
    if($good && $derep) {
        if($exactonly) {
            $md5 = md5_hex($seqn);
            if(exists $dereptypes{0}) {
                if(exists $md5s{$md5}) { #forward duplicate
                    $md5s{$md5}->{0}++;
                    if($derepmin <= $md5s{$md5}->{0}+1) {
                        $good = 0;
                        $filtercount{derep}++;
                    }
                }
            }
            if($good && exists $dereptypes{3}) {
                $md5r = md5_hex(&revcompuc($seqn));
                if(exists $md5s{$md5r}) { #reverse duplicate
                    $md5s{$md5}->{3}++;
                    if($derepmin <= $md5s{$md5}->{3}+1) {
                        $good = 0;
                        $filtercount{derep}++;
                    }
                }
            }
            unless(exists $md5s{$md5}) {
                $md5s{$md5} = {0 => 0, 3 => 0};
            }
        } else {
            push(@seqsP,[$seqn,$goodcount++,$length]);
            #keep write data for possible duplicates
            if($params{out_format} == 1) { #FASTA
                push(@printtmp,[$sid,$header,$seq,$begin,$end,'']);
            } else { # FASTQ or FASTA+QUAL or FASTQ+FASTA or FASTQ+FASTA+QUAL
                push(@printtmp,[$sid,$header,$seq,$begin,$end,$qual]);
            }
        }
    }

    if($good && (($derep && $exactonly) || !$derep)) { #passed filters
        $seqcount++;
        $seqbases += $length;
        return if($nogood);
        #check if change of sequence ID
        if(exists $params{seq_id}) {
            if($mappings) {
                print $fhmappings join("\t",$sid,$params{seq_id}.$seqcount)."\n";
            }
            $sid = $params{seq_id}.$seqcount;
        }
        if(exists $params{rm_header}) {
            $header = undef;
        }
        #trim if necessary
        if($begin) {
            $seq = substr($seq,$begin);
            $qual = substr($qual,$begin) if(defined $qual && length($qual));
        }
        if($end) {
            $length = length($seq);
            $seq = substr($seq,0,$length-$end);
            $qual = substr($qual,0,$length-$end) if(defined $qual && length($qual));
        }
        #change case
        if(exists $params{seq_case}) {
            if($params{seq_case} eq 'lower') { #lower case
                $seq = lc($seq);
            } elsif($params{seq_case} eq 'upper') { #upper case
                $seq = uc($seq);
            }
        }
        #convert between DNA and RNA
        if(exists $params{dna_rna}) {
            if($params{dna_rna} eq 'dna') { #RNA to DNA
                $seq =~ tr/Uu/Tt/;
            } elsif($params{dna_rna} eq 'rna') { #DNA to RNA
                $seq =~ tr/Tt/Uu/;
            }
        }
        #write data
        if($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5) { # FASTQ
            &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
            if($stdoutgood) {
                print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                print STDOUT $qual."\n";
            } else {
                print $fhgood '@'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
                print $fhgood '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                print $fhgood $qual."\n";
            }
        }
        if($params{out_format} == 1 || $params{out_format} == 2 || $params{out_format} == 4 || $params{out_format} == 5) { #FASTA
            #set line length
            if($linelen) {
                $seq =~ s/(.{$linelen})/$1\n/g;
                $seq =~ s/\n$//;
            }
            if($stdoutgood) {
                print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
            } elsif($params{out_format} == 1 || $params{out_format} == 2) {
                print $fhgood '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood $seq."\n";
            } else {
                print $fhgood3 '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhgood3 $seq."\n";
            }
        }
        if($params{out_format} == 2 || $params{out_format} == 5) { #QUAL
            &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
            print $fhgood2 '>'.$sid.($header ? ' '.$header : '')."\n";
            print $fhgood2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
        }
    } elsif($good && $derep && !$exactonly) {
        #do nothing as sequences will be used for duplicate check
    } else { #filtered out
        $badcount++;
        $badbases += length($seq);
        return if($nobad);
        #write data
        if($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5) { # FASTQ
            &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
            if($stdoutbad) {
                print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
                print STDOUT '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                print STDOUT $qual."\n";
            } else {
                print $fhbad '@'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad $seq."\n";
                print $fhbad '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                print $fhbad $qual."\n";
            }
        }
        if($params{out_format} == 1 || $params{out_format} == 2 || $params{out_format} == 4 || $params{out_format} == 5) { #FASTA
            #set line length
            if($linelen) {
                $seq =~ s/(.{$linelen})/$1\n/g;
                $seq =~ s/\n$//;
            }
            if($stdoutbad) {
                print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                print STDOUT $seq."\n";
            } elsif($params{out_format} == 1 || $params{out_format} == 2) {
                print $fhbad '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad $seq."\n";
            } else {
                print $fhbad3 '>'.$sid.($header ? ' '.$header : '')."\n";
                print $fhbad3 $seq."\n";
            }
        }
        if($params{out_format} == 2 || $params{out_format} == 5) { #QUAL
            &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
            print $fhbad2 '>'.$sid.($header ? ' '.$header : '')."\n";
            print $fhbad2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
        }
    }
}

#dereplicate sequences
sub derepSeqs {
    my $numseqs = scalar(@seqsP);
    if($derep && $numseqs) {
        my ($sid,$seq,$qual,$header,$begin,$end);
        my ($dcounts,undef,$dupls) = &checkForDupl(\@seqsP,\%dereptypes,$numseqs);

        print STDERR "Write results to output file(s)\n" if(exists $params{verbose});
        #for progress bar
        my $progress = 0;
        my $counter = 1;
        my $part = int($numseqs/100);
        print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});

        foreach my $i (0..$numseqs-1) {
            #progress bar stuff
	    $counter++;
	    if($counter > $part) {
		$counter = 1;
		$progress++;
		$progress = 99 if($progress > 99);
		print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
	    }
            #get data
            $sid = $printtmp[$i]->[0];
            $seq = $printtmp[$i]->[2];
            $qual = $printtmp[$i]->[5];
            $header = $printtmp[$i]->[1];
            $begin = $printtmp[$i]->[3];
            $end = $printtmp[$i]->[4];
            #write data
            if(exists $dupls->{$i} || (exists $params{seq_num} && $params{seq_num} <= $seqcount)) { #bad
                $filtercount{derep}++;
                $badcount++;
                $badbases += length($seq);
                next if($nobad);
                #write data
                if($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5) { # FASTQ
                    &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                    if($stdoutbad) {
                        print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                        print STDOUT '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                        print STDOUT $qual."\n";
                    } else {
                        print $fhbad '@'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad $seq."\n";
                        print $fhbad '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                        print $fhbad $qual."\n";
                    }
                }
                if($params{out_format} == 1 || $params{out_format} == 2 || $params{out_format} == 4 || $params{out_format} == 5) { #FASTA
                    #set line length
                    if($linelen) {
                        $seq =~ s/(.{$linelen})/$1\n/g;
                        $seq =~ s/\n$//;
                    }
                    if($stdoutbad) {
                        print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                    } elsif($params{out_format} == 1 || $params{out_format} == 2) {
                        print $fhbad '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad $seq."\n";
                    } else {
                        print $fhbad3 '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhbad3 $seq."\n";
                    }
                }
                if($params{out_format} == 2 || $params{out_format} == 5) { #QUAL
                    &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                    print $fhbad2 '>'.$sid.($header ? ' '.$header : '')."\n";
                    print $fhbad2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
                }
            } else { #good
                $seqcount++;
                $seqbases += (length($seq)-$begin-$end);
                next if($nogood);
                #check if change of sequence ID
                if(exists $params{seq_id}) {
                    if($mappings) {
                        print $fhmappings join("\t",$sid,$params{seq_id}.$seqcount)."\n";
                    }
                    $sid = $params{seq_id}.$seqcount;
                }
                if(exists $params{rm_header}) {
                    $header = undef;
                }
                #trim if necessary
                if($begin) {
                    $seq = substr($seq,$begin);
                    $qual = substr($qual,$begin) if(defined $qual && length($qual));
                }
                if($end) {
                    $length = length($seq);
                    $seq = substr($seq,0,$length-$end);
                    $qual = substr($qual,0,$length-$end) if(defined $qual && length($qual));
                }
                #change case
                if(exists $params{seq_case}) {
                    if($params{seq_case} eq 'lower') { #lower case
                        $seq = lc($seq);
                    } elsif($params{seq_case} eq 'upper') { #upper case
                        $seq = uc($seq);
                    }
                }
                #convert between DNA and RNA
                if(exists $params{dna_rna}) {
                    if($params{dna_rna} eq 'dna') { #RNA to DNA
                        $seq =~ tr/Uu/Tt/;
                    } elsif($params{dna_rna} eq 'rna') { #DNA to RNA
                        $seq =~ tr/Tt/Uu/;
                    }
                }
                #write data
                if($params{out_format} == 3 || $params{out_format} == 4 || $params{out_format} == 5) { # FASTQ
                    &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                    if($stdoutgood) {
                        print STDOUT '@'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                        print STDOUT '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                        print STDOUT $qual."\n";
                    } else {
                        print $fhgood '@'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood $seq."\n";
                        print $fhgood '+'.(exists $params{no_qual_header} ? '' : $sid.($header ? ' '.$header : ''))."\n";
                        print $fhgood $qual."\n";
                    }
                }
                if($params{out_format} == 1 || $params{out_format} == 2 || $params{out_format} == 4 || $params{out_format} == 5) { #FASTA
                    #set line length
                    if($linelen) {
                        $seq =~ s/(.{$linelen})/$1\n/g;
                        $seq =~ s/\n$//;
                    }
                    if($stdoutgood) {
                        print STDOUT '>'.$sid.($header ? ' '.$header : '')."\n";
                        print STDOUT $seq."\n";
                    } elsif($params{out_format} == 1 || $params{out_format} == 2) {
                        print $fhgood '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood $seq."\n";
                    } else {
                        print $fhgood3 '>'.$sid.($header ? ' '.$header : '')."\n";
                        print $fhgood3 $seq."\n";
                    }
                }
                if($params{out_format} == 2 || $params{out_format} == 5) { #QUAL
                    &printError("missing quality data for sequence \"$sid\" or greater number of sequences than available quality scores") unless(defined $qual);
                    print $fhgood2 '>'.$sid.($header ? ' '.$header : '')."\n";
                    print $fhgood2 &convertQualArrayToString(&convertQualAsciiToNums($qual),$linelen)."\n";
                }
            }
        }
        print STDERR "\r\tdone               \n" if(exists $params{verbose});
    }
}

#calculate summary statistics from sequences
sub calcSeqStats {
    my ($seq,$length,$stats,$kmers,$odds,$counts,$pair) = @_;

    #length related: min, max, range, mean, stddev, mode
    if(exists $params{stats_len} || exists $params{stats_assembly}) {
        $counts->{length}->{$length}++;
    }

    #dinucleotide odds ratio related: aatt, acgt, agct, at, catg, ccgg, cg, gatc, gc, ta
    if(exists $params{stats_dinuc}) {
        &dinucOdds($seq,$length,$odds);
    }

    #tag related: probability of 5' and 3' tag sequence based on kmer counts
    if(exists $params{stats_tag}) {
        #get kmers
        if($length >= 5) {
            #get 5' and 3' ends
            my $str5 = substr($seq,0,5);
            my $str3 = substr($seq,$length-5);
            unless($str5 eq 'AAAAA' || $str5 eq 'TTTTT' || $str5 eq 'CCCCC' || $str5 eq 'GGGGG' || $str5 eq 'NNNNN') {
                $kmers->{5}->{$str5}++;
            }
            unless($str3 eq 'AAAAA' || $str3 eq 'TTTTT' || $str3 eq 'CCCCC' || $str3 eq 'GGGGG' || $str3 eq 'NNNNN') {
                $kmers->{3}->{$str3}++;
            }
        }
        #check for MID tags
        if($length >= $MIDCHECKLENGTH) {
            my $str5 = substr($seq,0,$MIDCHECKLENGTH);
            foreach my $mid (keys %MIDS) {
                if(index($str5,$mid) != -1) {
                    $MIDS{$mid}++;
                    last;
                }
            }
        }
    }

    #ambiguous base N related: seqswithn, maxp
    if(exists $params{stats_ns}) {
        my $bylength = 100/$length;
        my $ns = ($seq =~ tr/N//);
        if($pair) {
            $stats->{stats_ns2}->{seqswithn}++ if($ns > 0);
            $stats->{stats_ns2}->{maxn} = $ns if($ns > ($stats->{stats_ns2}->{maxn}||0));
            $ns = ($ns > 0 && $ns*$bylength < 1 ? 1 : sprintf("%d",$ns*$bylength));
            $stats->{stats_ns2}->{maxp} = $ns if($ns > ($stats->{stats_ns2}->{maxp}||0));
        } else {
            $stats->{stats_ns}->{seqswithn}++ if($ns > 0);
            $stats->{stats_ns}->{maxn} = $ns if($ns > ($stats->{stats_ns}->{maxn}||0));
            $ns = ($ns > 0 && $ns*$bylength < 1 ? 1 : sprintf("%d",$ns*$bylength));
            $stats->{stats_ns}->{maxp} = $ns if($ns > ($stats->{stats_ns}->{maxp}||0));
        }
    }
}

#dinucleotide odds ratio calculation
sub dinucOdds {
    my ($seq,$length,$odds) = @_;
    my ($mononum,$dinum,$i,$x,$y);
    my %di = %DN_DI;
    my (%mono,$lengthtmp);
    my @tmp = split(/N+/,$seq);
    foreach(@tmp) {
        $lengthtmp = length($_)-1;
	next unless($lengthtmp > 0);
	$mono{AT} += ($_ =~ tr/AT//);
	$mono{GC} += ($_ =~ tr/GC//);
        $i = 0;
        while($i < $lengthtmp) {
            $di{substr($_,$i++,2)}++;
        }
    }
    $dinum = sum(values %di);

    if($dinum) {
        $mononum = sum(values %mono);
	my $factor = 2 * $mononum * $mononum / $dinum;
	my $AT = $mono{AT};
	my $GC = $mono{GC};
	if($AT) {
	    my $AT2 = $factor / ($AT * $AT);
	    $odds->{'AATT'} += ($di{'AA'} + $di{'TT'}) * $AT2;
	    $odds->{'AT'}   +=          2 * $di{'AT'}  * $AT2;
	    $odds->{'TA'}   +=          2 * $di{'TA'}  * $AT2;
	    if($GC) {
		my $ATGC = $factor / ($AT * $GC);
		$odds->{'ACGT'} += ($di{'AC'} + $di{'GT'}) * $ATGC;
		$odds->{'AGCT'} += ($di{'AG'} + $di{'CT'}) * $ATGC;
		$odds->{'CATG'} += ($di{'CA'} + $di{'TG'}) * $ATGC;
		$odds->{'GATC'} += ($di{'GA'} + $di{'TC'}) * $ATGC;
		my $GC2 = $factor / ($GC * $GC);
		$odds->{'CCGG'} += ($di{'CC'} + $di{'GG'}) * $GC2;
		$odds->{'CG'}   +=          2 * $di{'CG'}  * $GC2;
		$odds->{'GC'}   +=          2 * $di{'GC'}  * $GC2;
	    }
	} elsif($GC) {
	    my $GC2 = $factor / ($GC * $GC);
	    $odds->{'CCGG'} += ($di{'CC'} + $di{'GG'}) * $GC2;
	    $odds->{'CG'}   +=          2 * $di{'CG'}  * $GC2;
	    $odds->{'GC'}   +=          2 * $di{'GC'}  * $GC2;
	}
    }
}

#calculate basic stats from an hash of number->count values
sub generateStats {
    my $counts = shift;
    my ($min,$max,$modeval,$mode,$mean,$count,$std,$x,$c,@vals,$num,$median,%stats);

    #min, max, mode and modeval
    $min = -1;
    $max = $modeval = $mean = $count = $std = $num = 0;
    while (($x, $c) = each(%$counts)) {
        if($min == -1) {
            $min = $x;
        } elsif($min > $x) {
            $min = $x;
        }
        if($max < $x) {
            $max = $x;
        }
        if($modeval < $c) {
            $modeval = $c;
            $mode = $x;
        }
        $mean += $x*$c;
        $count += $c;
        foreach(1..$c) {
            push(@vals,$x);
            $num++;
        }
    }

    #mean and stddev
    $mean /= $count;
    while (($x, $c) = each(%$counts)) {
        $std += $c*(($x-$mean)**2);
    }

    #median
    if($num == 1) {
        $median = $vals[0];
    } elsif($num == 2) {
        $median = ($vals[0]+$vals[1])/2;
    } else {
        @vals = sort {$a <=> $b} @vals;
        if($num % 2) {
            $median = $vals[($num-1)/2];
        } else {
            $median = ($vals[$num/2]+$vals[$num/2-1])/2;
        }
    }

    #save stats
    $stats{min} = $min;
    $stats{max} = $max;
    $stats{range} = $max-$min+1;
    $stats{modeval} = $modeval;
    $stats{mode} = $mode;
    $stats{mean} = sprintf("%.2f",$mean);
    $stats{stddev} = sprintf("%.2f",($std/$count)**(1/2));
    $stats{median} = $median;

    return \%stats;
}

sub generateStatsType {
    my $counts = shift;
    my (%stats,$min,$max,$modeval,$mode,$mean,$std,$x,$c,@vals,$num,$median,$p25,$p75,$numq,$i,$j,$median1,$median2,$p251,$p252,$p751,$p752);

    foreach my $kind (keys %$counts) {
	@vals = ();
	$min = -1;
	$max = $modeval = $mean = $std = $num = 0;
        foreach my $x1 (sort {$a <=> $b} keys %{$counts->{$kind}}) {
            $c = $counts->{$kind}->{$x1};
	    if($min == -1) {
		$min = $x1;
	    }
	    if($max < $x1) {
		$max = $x1;
	    }
	    if($modeval < $c) {
		$modeval = $c;
		$mode = $x1;
	    }
	    $mean += $x1*$c;
	    $num += $c;
            push(@vals,[$c,$x1]); #count, values
	}

	$mean /= $num;
	while (($x, $c) = each(%{$counts->{$kind}})) {
	    $std += $c*(($x-$mean)**2);
	}

	if($num == 1) {
            $median = $p25 = $p75 = $vals[0]->[1];
	} elsif($num == 2) {
            if($vals[0]->[0] == 1) { #two different numbers
                $p25 = $vals[0]->[1];
                $p75 = $vals[1]->[1];
                $median = ($vals[0]->[1]+$vals[1]->[1])/2;
            } else {
                $p25 = $p75 = $median = $vals[0]->[1]; #both same
            }
	} elsif($num > 2) {
	    if($num % 2) {
                $i = 0;
                $j = 0;
                while($i <= ($num-1)/2) {
                    $median = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
		$numq = ($num+1)/2;
	    } else {
                $i = 0;
                $j = 0;
                while($i <= ($num/2-1)) {
                    $median1 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $median2 = $median1;
                while($i <= ($num/2)) {
                    $median2 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $median = ($median1 + $median2)/2;
		$numq = $num/2;
	    }
	    if($numq % 2) {
                $i = 0;
                $j = 0;
                while($i <= (($numq-1)/2)) {
                    $p25 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $p75 = $p25;
                while($i <= ($num-($numq-1)/2-1)) {
                    $p75 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
	    } else {
                $i = 0;
                $j = 0;
                while($i <= ($numq/2-1)) {
                    $p251 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $p252 = $p251;
                while($i <= ($numq/2)) {
                    $p252 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $p751 = $p252;
                while($i <= ($num-$numq/2-1)) {
                    $p751 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $p752 = $p751;
                while($i <= ($num-$numq/2)) {
                    $p752 = $vals[$j]->[1];
                    $i += $vals[$j]->[0];
                    $j++;
                }
                $p25 = ($p251 + $p252) / 2;
                $p75 = ($p751 + $p752) / 2;
	    }
	} else {
            $median = $p25 = $p75 = 0;
        }

	$stats{$kind}->{min} = $min;
	$stats{$kind}->{max} = $max;
	$stats{$kind}->{range} = $max-$min+1;
	$stats{$kind}->{modeval} = $modeval;
	$stats{$kind}->{mode} = $mode;
	$stats{$kind}->{mean} = sprintf("%.2f",$mean);
	$stats{$kind}->{std} = sprintf("%.2f",($std/$num)**(1/2));
	$stats{$kind}->{median} = int($median);
	$stats{$kind}->{p25} = int($p25);
	$stats{$kind}->{p75} = int($p75);
    }

    return \%stats;
}

#requires seqs array with [upper-case seq, array index, length] for each entry
sub checkForDupl {
    #requires seqs array with [upper-case seq, array index, length] for each entry
    my ($seqs,$types,$numseqs) = @_;
    my (@sort,$num,%dupls,$pretype,$precount,%counts,%lens);
    #precount = number duplicates for the same sequence

    print STDERR "Check for duplicates\n" if(exists $params{verbose});

    #for progress bar
    my $progress = 1;
    my $counter = 1;
    my $part = int($numseqs*4/100);
    print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
    &printWeb("STATUS: duplicate-status $progress");

    #exact duplicates and prefix duplicates
    if(exists $types->{0} || exists $types->{1} || exists $types->{2}) {
	$precount = 0;
	$pretype = -1;
        @sort = sort {$a->[0] cmp $b->[0]} @$seqs;
        foreach my $i (0..$numseqs-2) {
            if(exists $types->{0} && $sort[$i]->[2] == $sort[$i+1]->[2] && $sort[$i]->[0] eq $sort[$i+1]->[0]) {
                $dupls{$sort[$i]->[1]} = 0;
                $lens{$sort[$i]->[2]}->{0}++;
		if($pretype == 0) {
		    $precount++;
		} else {
		    if($pretype == 1 && $precount) {
			$counts{$precount}->{$pretype}++;
		    }
		    $pretype = 0;
		    $precount = 1;
		}
            } elsif(exists $types->{1} && $sort[$i]->[2] < $sort[$i+1]->[2] && $sort[$i]->[0] eq substr($sort[$i+1]->[0],0,$sort[$i]->[2])) {
                $dupls{$sort[$i]->[1]} = 1;
                $lens{$sort[$i]->[2]}->{1}++;
		if($pretype == 1) {
		    $precount++;
		} else {
		    if($pretype == 0 && $precount) {
			$counts{$precount}->{$pretype}++;
		    }
		    $pretype = 1;
		    $precount = 1;
		}
            } else {
		if($precount) {
		    $counts{$precount}->{$pretype}++;
		    $precount = 0;
		}
		$pretype = -1;
	    }
            $sort[$i] = undef;
	    #progress bar stuff
	    $counter++;
	    if($counter > $part) {
		$counter = 1;
		$progress++;
		$progress = 99 if($progress > 99);
		print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                &printWeb("STATUS: duplicate-status $progress");
	    }
        }
	if($precount) {
	    $counts{$precount}->{$pretype}++;
	}
    }
    #suffix duplicates
    if(exists $types->{2}) {
        $num = 0;
        @sort = ();
        foreach(@$seqs) {
            next if(exists $dupls{$_->[1]});
            push(@sort,[(scalar reverse $_->[0]),$_->[1],$_->[2]]);
            $num++;
        }
	if($num > 1) {
	    $precount = 0;
	    $pretype = -1;
	    @sort = sort {$a->[0] cmp $b->[0]} @sort;
	    foreach my $i (0..$num-2) {
		if($sort[$i]->[2] < $sort[$i+1]->[2] && $sort[$i]->[0] eq substr($sort[$i+1]->[0],0,$sort[$i]->[2])) {
		    $dupls{$sort[$i]->[1]} = 2;
                    $lens{$sort[$i]->[2]}->{2}++;
		    if($pretype == 2) {
			$precount++;
		    } else {
			$pretype = 2;
			$precount = 1;
		    }
		} else {
		    if($precount) {
			$counts{$precount}->{$pretype}++;
			$precount = 0;
		    }
		    $pretype = -1;
		}
		$sort[$i] = undef;
		#progress bar stuff
		$counter++;
		if($counter > $part) {
		    $counter = 1;
		    $progress++;
		    $progress = 99 if($progress > 99);
		    print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                    &printWeb("STATUS: duplicate-status $progress");
		}
	    }
	    if($precount) {
		$counts{$precount}->{$pretype}++;
	    }
	}
    }
    #reverse complement exact and prefix/suffix duplicates
    if(exists $types->{3} || exists $types->{4}) {
        $num = 0;
        @sort = ();
        foreach(@$seqs) {
            if(exists $dupls{$_->[1]}) {
		$counter++;
		next;
	    }
            push(@sort,[$_->[0],$_->[1],$_->[2],0]);
            push(@sort,[&revcompuc($_->[0]),$_->[1],$_->[2],1]);
            $num += 2;
        }
	if($num > 1) {
	    $precount = 0;
	    $pretype = -1;
	    @sort = sort {$a->[0] cmp $b->[0]} @sort;
	    foreach my $i (0..$num-2) {
		unless($sort[$i]->[3] == $sort[$i+1]->[3] || $sort[$i]->[1] eq $sort[$i+1]->[1] || exists $dupls{$sort[$i]->[1]}) { #don't check if both same (original or revcomp) or already counted as dubs
		    if(exists $types->{3} && $sort[$i]->[2] == $sort[$i+1]->[2] && $sort[$i]->[0] eq $sort[$i+1]->[0]) {
			$dupls{$sort[$i]->[1]} = 3;
                        $lens{$sort[$i]->[2]}->{3}++;
			if($pretype == 3) {
			    $precount++;
			} else {
			    if($pretype == 4 && $precount) {
				$counts{$precount}->{$pretype}++;
			    }
			    $pretype = 3;
			    $precount = 1;
			}
		    } elsif(exists $types->{4} && $sort[$i]->[2] < $sort[$i+1]->[2] && $sort[$i]->[0] eq substr($sort[$i+1]->[0],0,$sort[$i]->[2])) {
			$dupls{$sort[$i]->[1]} = 4;
                        $lens{$sort[$i]->[2]}->{4}++;
			if($pretype == 4) {
			    $precount++;
			} else {
			    if($pretype == 3 && $precount) {
				$counts{$precount}->{$pretype}++;
			    }
			    $pretype = 4;
			    $precount = 1;
			}
		    } else {
			if($precount) {
			    $counts{$precount}->{$pretype}++;
			    $precount = 0;
			}
			$pretype = -1;
		    }
		}
		$sort[$i] = undef;
		#progress bar stuff
		$counter++;
		if($counter > $part) {
		    $counter = 1;
		    $progress++;
		    $progress = 99 if($progress > 99);
		    print STDERR "\r\tstatus: ".int($progress)." \%" if(exists $params{verbose});
                    &printWeb("STATUS: duplicate-status $progress");
		}
	    }
	    if($precount) {
		$counts{$precount}->{$pretype}++;
	    }
	}
    }
    print STDERR "\r\tdone               \n" if(exists $params{verbose});
    &printWeb("STATUS: duplicate-status 100");
    return (\%counts,\%lens,\%dupls);
}

#get the frequency of possible tags by shifting kmers by max 2 positions when aligned
sub getTagFrequency {
    my ($kmers) = @_;

    #find most abundant kmer counts
    my $percentone = $numseqs/100;
    my $percentten = $numseqs/10;
    my %most;
    foreach my $sp (keys %$kmers) {
	$most{$sp}->{max} = 0;
	foreach(keys %{$kmers->{$sp}}) {
#	    next if($_ eq 'A'x5 || $_ eq 'T'x5 || $_ eq 'C'x5 || $_ eq 'G'x5 || $_ eq 'N'x5);
	    if($kmers->{$sp}->{$_} >= $percentten) {
		$most{$sp}->{ten}++;
	    } elsif($kmers->{$sp}->{$_} >= $percentone) {
		$most{$sp}->{one}++;
	    }
	    #get max count
	    $most{$sp}->{max} = $kmers->{$sp}->{$_} if($most{$sp}->{max} < $kmers->{$sp}->{$_});
	}
    }

    #filter kmers by frequency - threshold of >10% occurrence -> max of 9 different kmers or more if there is non with >10% occurrence
    my $numseqssub = $numseqs/10;
    my $onecount = 2;
    foreach my $sp (keys %$kmers) {
	foreach(keys %{$kmers->{$sp}}) {
	    if(exists $most{$sp}->{ten} && $most{$sp}->{ten} > 0) {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} < $percentten);
	    } elsif(exists $most{$sp}->{one} && $most{$sp}->{one} > 0) {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} < $percentone);
	    } else {
		delete $kmers->{$sp}->{$_} if($kmers->{$sp}->{$_} != $most{$sp}->{max});
	    }
	}
    }

    my (%kmersum,%kmershift);
    foreach my $sp (sort {$b <=> $a} keys %$kmers) { #5' before 3'
	#if more than one kmer in array, test if shifted by max 2 positions
	my $numkmer = scalar(keys %{$kmers->{$sp}});

	if($numkmer > 1) {
	    my @matrix;
	    my @kmersort = sort {$kmers->{$sp}->{$b} <=> $kmers->{$sp}->{$a}} keys %{$kmers->{$sp}};
	    foreach my $i (0..($numkmer-2)) {
		foreach my $j (($i+1)..($numkmer-1)) {
		    $matrix[$i]->[$j-($i+1)] = &align2seqs($kmersort[$j],$kmersort[$i]);
		}
	    }
	    my $countgood = 0;
	    foreach my $i (0..($numkmer-2)) {
		unless(@{$matrix[0]->[$i]}) { #not matching
		    my $count = 0;
		    foreach my $j (1..($numkmer-2)) {
			$count++;
			last if(defined $matrix[$j]->[$i-$j]); #found shift using other kmers
		    }
		    if($count < ($numkmer-1) && $i > 0) {
			my $sum = 0;
			my $sign;
			foreach my $j (0..$count) {
			    next unless(defined $matrix[$j] && defined $matrix[$j]->[$i-1]); #fix: 08/2010
			    if(defined $sign) {
				if(($sign < 0 && (defined $matrix[$j]->[$i-1]->[0] && $matrix[$j]->[$i-1]->[0] < 0)) || ($sign > 0 && (defined $matrix[$j]->[$i-1]->[0] && $matrix[$j]->[$i-1]->[0] > 0))) {
				    $sum += $matrix[$j]->[$i-1]->[0];
				} elsif(($sign < 0 && (defined $matrix[$j]->[$i-1]->[1] && $matrix[$j]->[$i-1]->[1] < 0)) || $sign > 0 && (defined $matrix[$j]->[$i-1]->[1] && $matrix[$j]->[$i-1]->[1] > 0)) {
				    $sum += $matrix[$j]->[$i-1]->[1];
				}
			    } elsif(defined $matrix[$j]->[$i-1]->[0]) {
				$sum += $matrix[$j]->[$i-1]->[0];
			    }
			    $sign = ((defined $matrix[$j]->[$i-1]->[0] && $matrix[$j]->[$i-1]->[0] < 0) ? -1 : 1);
			}
			$matrix[0]->[$i] = [$sum] if(defined $sign); #fix: 08/2010
		    }
		}
		unless(@{$matrix[0]->[$i]}) {
		    last;
		} else {
		    $countgood++;
		}
	    }
	    if($countgood) {
		my $min;
		if($sp == 3) { #3' prime end, 5 for 5' end
		    #find maximum shift to right (pos value)
		    $min = -100;
		    foreach my $i (0..($countgood-1)) {
			$min = ((defined $matrix[0]->[$i]->[0] && $min > $matrix[0]->[$i]->[0]) ? $min : $matrix[0]->[$i]->[0]);
		    }
		    if($min > 0) {
			$min = -$min;
		    } else {
			$min = 0;
		    }
		} else {
		    #find maximum shift to left (neg value)
		    $min = 100;
		    foreach my $i (0..($countgood-1)) {
			$min = ((defined $matrix[0]->[$i]->[0] && $min < $matrix[0]->[$i]->[0]) ? $min : $matrix[0]->[$i]->[0]);
		    }
		    if($min < 0) {
			$min = abs($min);
		    } else {
			$min = 0;
		    }
		}
#		$kmershift{$sp}->{$kmersort[0]} = $min;
		$kmersum{$sp} += $kmers->{$sp}->{$kmersort[0]};
		foreach my $i (0..($countgood-1)) {
#		    $kmershift{$sp}->{$kmersort[$i+1]} = $matrix[0]->[$i]->[0]+$min;
		    $kmersum{$sp} += $kmers->{$sp}->{$kmersort[$i+1]};
		}
	    } else {
		my $tmp = (sort {$kmers->{$sp}->{$b} <=> $kmers->{$sp}->{$a}} keys %{$kmers->{$sp}})[0];
#		$kmershift{$sp}->{$tmp} = 0;
		$kmersum{$sp} += $kmers->{$sp}->{$tmp};
	    }
	} elsif($numkmer == 1) {
	    my $tmp = (keys %{$kmers->{$sp}})[0];
#	    $kmershift{$sp}->{$tmp} = 0;
	    $kmersum{$sp} += $kmers->{$sp}->{$tmp};
	}
    }

   return \%kmersum;
}

sub align2seqs {
    my ($seq1,$seq2) = @_;
    my @shift;

    #get number of shifted positions
    if(substr($seq1,0,4) eq substr($seq2,1,4)) { #shift right by 1
	push(@shift,1);
    } elsif(substr($seq1,0,3) eq substr($seq2,2,3)) { #shift right by 2
	push(@shift,2);
    }
    if(substr($seq1,1,4) eq substr($seq2,0,4)) { #shift left by 1
	push(@shift,-1);
    } elsif(substr($seq1,2,3) eq substr($seq2,0,3)) { #shift left by 2
	push(@shift,-2);
    }

    return \@shift;
}

sub revcompuc {
    my $seq = shift;
    $seq = scalar reverse $seq;
    $seq =~ tr/GATC/CTAG/;
    return $seq;
}

sub compuc {
    my $seq = shift;
    $seq =~ tr/GATC/CTAG/;
    return $seq;
}

#get data for graphs
sub getSeqStats {
    my ($graphdata,$seqgd,$length) = @_;
    if($length > $maxlength) {
        $maxlength = $length;
    }
    my ($gc,$ns,$begin,$end,$str5,$str3,$bylength,$tmp);
    $begin = $end = $gc = $ns = 0;
    #get length
    $bylength = 100/$length;
    #get 5' and 3' ends
    if($webstats{pt} || $webstats{ts} || $graphstats{pt} || $graphstats{ts}) {
        $str5 = substr($seqgd,0,5);
        $str3 = substr($seqgd,$length-5);
    }
    #GC content
    if($webstats{gc} || $graphstats{gc}) {
        $gc = ($seqgd =~ tr/GC//);
        $gc = sprintf("%d",$gc*$bylength);
    }

    #N's
    if($webstats{ns} || $graphstats{ns}) {
        $ns = ($seqgd =~ tr/N//);
        $ns = ($ns > 0 && $ns*$bylength < 1 ? 1 : sprintf("%d",$ns*$bylength));
    }

    #tail stuff with min 5 char repeats
    if($webstats{pt} || $graphstats{pt}) {
        #at sequence 5'-end
        if($str5 eq 'AAAAA' || $str5 eq 'TTTTT') {
            my $tmpchar = substr($str5,0,1); #A or T
            $begin = 5;
            foreach(5..$length-1) {
                $tmp = substr($seqgd,$_,1);
                last unless($tmp eq $tmpchar || $tmp eq 'N');
                $begin++;
            }
        }
        #at sequence 3'-end
        if($str3 eq 'AAAAA' || $str3 eq 'TTTTT') {
            my $tmpchar = substr($str3,0,1); #A or T
            $end = 5;
            foreach (reverse 0..$length-6) {
                $tmp = substr($seqgd,$_,1);
                last unless($tmp eq $tmpchar || $tmp eq 'N');
                $end++;
            }
        }
    }

    #get base frequencies
    if($webstats{ts} || $graphstats{ts}) {
        if($length >= $TAG_LENGTH) {
            foreach my $i (0..$TAG_LENGTH-1) {
                $graphdata->{freqs}->{5}->{$i}->{substr($seqgd,$i,1)}++;
                $graphdata->{freqs}->{3}->{$i}->{substr($seqgd,$length-$TAG_LENGTH+$i,1)}++;
            }
        }
        #get kmers
        if($length >= 5) {
            unless($begin > 0 || $str5 eq 'CCCCC' || $str5 eq 'GGGGG' || $str5 eq 'NNNNN') {
                $graphdata->{kmers}->{5}->{$str5}++;
            }
            unless($end > 0 || $str3 eq 'CCCCC' || $str3 eq 'GGGGG' || $str3 eq 'NNNNN') {
                $graphdata->{kmers}->{3}->{$str3}++;
            }
        }
        #check for MID tags
        if($length >= $MIDCHECKLENGTH) {
            $str5 = substr($seqgd,0,$MIDCHECKLENGTH);
            foreach my $mid (keys %MIDS) {
                if(index($str5,$mid) != -1) {
                    $graphdata->{mids}->{$mid}++;
                    last;
                }
            }
        }
    }

    #calculate sequence complexity
    if($webstats{sc} || $graphstats{sc}) {
        my ($rest,$steps,@dust,@entropy,$mean,$str,%counts,$num,$dustscore,$entropyval,$bynum);
        if($length <= $WINDOWSIZE) {
            $rest = $length;
            $steps = 0;
        } else {
            $steps = int(($length - $WINDOWSIZE) / $WINDOWSTEP) + 1;
            $rest = $length - $steps * $WINDOWSTEP;
            unless($rest > $WINDOWSTEP) {
                $rest += $WINDOWSTEP;
                $steps--;
            }
        }
        #dust and entropy
        $num = $WINDOWSIZE-2;
        $bynum = 1/$num;
        $num--;
        foreach my $i (0..$steps-1) {
            $str = substr($seqgd,($i * $WINDOWSTEP),$WINDOWSIZE);
            %counts = ();
            foreach my $i (@WINDOWSIZEARRAY) {
                $counts{substr($str,$i,3)}++;
            }
            #dust and entropy
            $dustscore = $entropyval = 0;
            foreach(values %counts) {
                $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
                $entropyval -= ($_ * $bynum) * log($_ * $bynum);
            }
            push(@dust,($dustscore * $bynum));
            push(@entropy,($entropyval * $ONEOVERLOG62));
        }
        #last step
        if($rest > 5) {
            $str = substr($seqgd,($steps * $WINDOWSTEP),$rest);
            %counts = ();
            $num = $rest-2;
            foreach my $i (0..($num - 1)) {
                $counts{substr($str,$i,3)}++;
            }
            $dustscore = $entropyval = 0;
            $bynum = 1/$num;
            foreach(values %counts) {
                $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
                $entropyval -= ($_ * $bynum) * log($_ * $bynum);
            }
            push(@dust,(($dustscore / ($num-1)) * (($WINDOWSIZE - 2) / $num)));
            push(@entropy,($entropyval / log($num)));
        } else {
            push(@dust,31); #to assign a maximum score based on the scaling factor 100/31
            push(@entropy,0);
        }

        $mean = &getArrayMean(@dust);
        $mean = int($mean * 100 / 31); #scale to 100
        $graphdata{compldust}->{$mean}++;
        if(!exists $graphdata{complvals}->{dust}->{minval} || $graphdata{complvals}->{dust}->{minval} > $mean) {
            $graphdata{complvals}->{dust}->{minval} = $mean;
            $graphdata{complvals}->{dust}->{minseq} = &getGraphDataSequence($seqgd);
        }
        if(!exists $graphdata{complvals}->{dust}->{maxval} || $graphdata{complvals}->{dust}->{maxval} < $mean) {
            $graphdata{complvals}->{dust}->{maxval} = $mean;
            $graphdata{complvals}->{dust}->{maxseq} = &getGraphDataSequence($seqgd);
        }
        $mean = &getArrayMean(@entropy);
        $mean = int($mean * 100); #scale to 100
        $graphdata{complentropy}->{$mean}++;
        if(!exists $graphdata{complvals}->{entropy}->{minval} || $graphdata{complvals}->{entropy}->{minval} > $mean) {
            $graphdata{complvals}->{entropy}->{minval} = $mean;
            $graphdata{complvals}->{entropy}->{minseq} = &getGraphDataSequence($seqgd);
        }
        if(!exists $graphdata{complvals}->{entropy}->{maxval} || $graphdata{complvals}->{entropy}->{maxval} < $mean) {
            $graphdata{complvals}->{entropy}->{maxval} = $mean;
            $graphdata{complvals}->{entropy}->{maxseq} = &getGraphDataSequence($seqgd);
        }
    }

    #calculate dinucleotide odd ratios
    if($webstats{dn} || $graphstats{dn}) {
        &dinucOdds($seqgd,$length,\%{$graphdata{dinucodds}});
    }

    #store counts
    if($webstats{ld} || $webstats{ld} || $graphstats{ld} || $graphstats{ld}) {
        $graphdata->{counts}->{length}->{$length}++;
    }
    if($begin) {
	$graphdata->{counts}->{tail5}->{$begin}++;
    }
    if($end) {
	$graphdata->{counts}->{tail3}->{$end}++;
    }
    if($webstats{gc} || $graphstats{gc}) {
        $graphdata->{counts}->{gc}->{$gc}++;
    }
    if($ns) {
	$graphdata->{counts}->{ns}->{$ns}++;
    }

    return 1;
}

sub getGraphDataSequence {
    my $seq = shift;
    if(length($seq) > $GRAPH_DATA_SEQ_MAX_LENGTH) {
        return substr($seq,0,$GRAPH_DATA_SEQ_MAX_LENGTH)."...";
    } else {
        return $seq;
    }
}

sub getQualStats {
    my ($graphdata,$qual,$length) = @_;

    #check if quality values are available
    return 0 unless($qual && ($webstats{qd} || $graphstats{qd}));

    #calculate decimal values of quals
    my ($vals,$err);
    if(exists $params{phred64}) { #scale data to Phred scale if necessary
        ($vals,$err) = &convertQualAsciiToNumsPhred64($qual);
        if($err) {
            &printError("The sequence quality scores are not in Phred+64 format");
        }
    } else {
        $vals = &convertQualAsciiToNums($qual);
    }

    #mean quality score
    $graphdata->{qualsmean}->{int(&getArrayMean(@$vals))}++;

    my ($factor,$tmp,$count,$xmax,$bin,$tmpbin,$step);

    if($scale == 1) { #relative
	#qual
	if($length == 100) {
	    foreach my $i (0..99) {
		$graphdata->{quals}->{$i}->{$vals->[$i]}++;
	    }
	} elsif($length < 100) { #stretch
	    $factor = 100/$length;
	    foreach my $i (0..$length-1) {
		$tmp = $vals->[$i];
		foreach my $j (int($i*$factor)..int(($i+1)*$factor)-1) {
		    $graphdata->{quals}->{$j}->{$tmp}++;
		}
	    }
	} elsif($length > 100) { #shrink
	    $factor = $length/100;
	    foreach my $i (0..99) {
		$tmp = $count = 0;
		foreach my $j (int($i*$factor)..int(($i+1)*$factor)-1) {
		    $tmp += $vals->[$j];
		    $count++;
		}
		$graphdata->{quals}->{$i}->{int($tmp/$count)}++;
	    }
            #my $piece = int($length/100);
            #my $bypiece = 1/($piece+1);
            #my $start = 0;
            #my $end = 0;
            #my $rest = ($length % 100) - 1;
            #foreach my $i (0..$rest) {
            #    $end += $piece;
            #    $tmp = 0;
            #    foreach my $j ($start..$end) {
            #        $tmp += $vals->[$j];
            #    }
            #    $graphdata{quals}->{$i}->{int($tmp * $bypiece)}++;
            #    $start = ++$end;
            #}
            #$rest++;
            #$piece--;
            #$bypiece = 1/($piece+1);
            #foreach my $i ($rest..99) {
            #    $end += $piece;
            #    $tmp = 0;
            #    foreach my $j ($start..$end) {
            #        $tmp += $vals->[$j];
            #    }
            #    $graphdata{quals}->{$i}->{int($tmp * $bypiece)}++;
            #    $start = ++$end;
            #}
	}
    }
    #absolute
    foreach my $i (0..$length-1) {
        $graphdata->{quala}->{$i}->{$vals->[$i]}++;
    }
}

sub getBinVal {
    my $val = shift;
    my $step;
    if(!$val || $val <= 100) {
	return 1;
    } elsif($val < 10000) {
	return int($val/100)+($val % 100 ? 1 : 0);
    } elsif($val < 100000) {
	return 1000;
    } else {
	$step = 1000000;
	my $xmax = ($val % $step ? sprintf("%d",($val/$step+1))*$step : $val);
	return ($xmax/100);
    }
}

sub convertStringToInt {
    my $string = shift;
    $string =~ s/(.)/sprintf("%x",ord($1))/eg;
    return $string;
}

sub getFileName {
    my $str = shift;
    $str =~ s/^.*\/([^\/]+)$/$1/;
    return $str;
}
