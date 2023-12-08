#!/usr/bin/env perl

=head1 NAME

pogenom.pl - Calculates population genomic parameters from a VCF file

=head1 USAGE (minimum input)

  perl pogenom.pl --vcf_file <VCF_FILE> --out <OUTPUT_FILES_PREFIX> --genome_size <GENOME_SIZE>

or:
 
  perl pogenom.pl --vcf_file <VCF_FILE> --out <OUTPUT_FILES_PREFIX> --gff_file <GFF_FILE>
 
or:

  perl pogenom.pl --vcf_file <VCF_FILE> --out <OUTPUT_FILES_PREFIX> --fasta_file <FASTA_FILE>


=head1 REQUIRED ARGUMENTS

  --vcf_file <VCF_FILE          Specify vcf file with data from a single or multiple samples

  --out <OUTPUT_FILES_PREFIX    Specify the prefix of the output file name(s) (overwrites existing files with same names)

  --genome_size <GENOME_SIZE    Specify genome size (in bp; integer). Not required if --gff_file or --fasta_file with genome sequence is given
 


=head1 OPIONAL ARGUMENTS

  --gff_file <GFF_FILE>         Specify gff file. Either this, --genome_size or --fasta_file must be given
 
  --fasta_file <FASTA_FILE>     Specify fasta file. Either this, --genome_size or --gff_file must be given
 
  --genetic_code_file <GENETIC_CODE_FILE>   Specify genetic code file. E.g. standard_genetic_code.txt in the POGENOM distribution
 
  --loci_file <LOCI_FILE>       Specify file with ids of loci to be included

  --sample_file <SAMPLE_FILE>   Specify file with ids of samples to be included
 
  --min_count <MIN_COUNT>       Specify minimum coverage for a locus to be included for the sample
 
  --min_found <MIN_FOUND_IN>    Specify minimum number of samples that a locus needs to be present in to be included
 
  --subsample <SUBSAMPLE>       Specify coverage level at which to subsample
 
  --keep_haplotypes             If used, POGENOM will not split haplotypes into single-nucleotide variants, which is otherwise the default
 
  --vcf_format <VCF_FORMAT>     Specify VCF file format version. Can be set to freebayes (default) or GATK

  --fst_perm <FST_PERM>         Specify number of permutations (integer) for making randomised gene-wise Fst. Use with caution, output can be huge

  --pi_only                     If used, POGENOM will only calculate and output genome-wide pi
 
  --split_fasta_header          If used, if a fasta file header includes space character(s), POGENOM will only use the part preceeding the first space
 
  --help                        Prints this help message

[Press q to close this help message]

=cut

use Getopt::Long;
use List::MoreUtils qw(uniq);

$genome_size = undef;
$min_found_in = 1;
$min_count = 2;
$vcf_file = undef;
$gff_file = undef;
$fasta_file = undef;
$genetic_code_file = undef;
$loci_file = undef;
$sample_file = undef;
$outprefix = undef;
$reference = undef;
$keep_haplotypes = undef;
$pi_only = undef;
$split_fasta_header = undef;
$use_pseudocounts = undef;
$subsample = undef;
$vcf_format = "freebayes";
$na_if_missing_loci = 1;
$n_fst_permutations = undef;


&GetOptions('vcf_file=s' => \$vcf_file, 'vcf_format=s' => \$vcf_format, 'gff_file=s' => \$gff_file, 'fasta_file=s' => \$fasta_file, 'genetic_code_file=s' => \$genetic_code_file, 'output=s' => \$outprefix, 'min_count=i' => \$min_count, 'min_found=i' => \$min_found_in, 'ref=s' => \$reference, 'genome_size=i' => \$genome_size, 'keep_haplotypes!' => \$keep_haplotypes, 'loci_file=s' => \$loci_file, 'sample_file=s' => \$sample_file, 'subsample=s' => \$subsample, 'use_pseudocounts' => \$use_pseudocounts, 'fst_perm=i' => \$n_fst_permutations, 'pi_only!' => \$pi_only, 'split_fasta_header!' => \$split_fasta_header, 'h!' => \$help);

if (!$outprefix) {
    system ('perldoc', $0);
    exit;
}
if ($help) {
    system ('perldoc', $0);
    exit;
}
if (!$vcf_file) {
    system ('perldoc', $0);
    exit;
}
if (!$genome_size and !$gff_file and !$fasta_file) {
    system ('perldoc', $0);
    exit;
}
if ($min_count < 2) {
    print"Error: min_count cannot be set to <2\n";
    exit;
}
if ($vcf_format ne "freebayes") {
    if ($vcf_format ne "GATK") {
        print"\nError: Unrecognized vcf_format. Should be either freebayes (default) or GATK\n\n";
        exit;
    }
}

####################

$logtext = "vcf_file: $vcf_file\n";
if ($gff_file) {
    $logtext = $logtext."gff_file: $gff_file\n";
}
if ($fasta_file) {
    $logtext = $logtext."fasta_file: $fasta_file\n";
}
if ($genetic_code_file) {
    $logtext = $logtext."genetic_code_file: $genetic_code_file\n";
}
print"\n### Running pogenom ###\n";
if ($genome_size) {
    $logtext = $logtext."genome_size set to: $genome_size\n";
} elsif ($fasta_file) {
    $logtext = $logtext."genome_size calculated from fasta file\n";
} elsif ($gff_file) {
    $logtext = $logtext."genome_size calculated from GFF file\n";
}
$logtext = $logtext."min_count set to: $min_count\n";
if ($min_found_in == 0) {
    $logtext = $logtext."min_found set to: number of samples\n";
}
if ($min_found_in > 0) {
    $logtext = $logtext."min_found set to: $min_found_in\n";
}
if ($subsample) {
    $logtext = $logtext."Subsampling set to: $subsample reads per locus\n";
}
if ($pi_only) {
    $logtext = $logtext."Running in pi_only mode\n";
}
if ($loci_file) {
    $logtext = $logtext."Analysis restricted to loci in file: $loci_file\n";
    &read_loci_to_include;
}
if ($sample_file) {
    $logtext = $logtext."Analysis restricted to samples in file: $sample_file\n";
    &read_samples_to_include;
}
if ($n_fst_permutations) {
    $logtext = $logtext."Permuted gene-wise fst calculated with: $n_fst_permutations permutations\n";
}
print"$logtext\n";
print"\n### Read variant data ###\n";
if ($keep_haplotypes) {
    &get_snp_data_combined_vcf;
} else {
    &get_snp_data_combined_vcf_split_haplotypes;
}
if ($subsample) {
    &subsample_allele_counts;
    #&subsample_allele_counts_dupl;
}
if ($gff_file) {
    print"\n### Reading GFF file ###\n";
    &read_gff;
}
if ($fasta_file) {
    print"\n### Reading fasta sequence file ###\n";
    &read_fasta;
}
if ($genetic_code_file) {
    print"\n### Reading Genetic Code file ###\n";
    &read_genetic_code;
}
print"\n### Calculating Nucleotide Diversity (pi) ###\n";
&calc_pi;
if (@samples > 1) {
    &estimate_genome_coverage;
}
if ($gff_file and !$pi_only) {
    print"\n### Calculating Gene-wise Nucleotide Diversity (pi) ###\n";
    &calc_per_gene_pi;
    if ($genetic_code_file) {
        print"\n### Calculating Gene-wise Aminoacid Diversity (aa-pi) ###\n";
        &calc_per_gene_aminoacid_pi;
        print"\n### Calculating Aminoacid Frequencies ###\n";
        &calc_aminoacid_frequencies;
        print"\n### Calculating Gene-wise pN/pS ###\n";
        &calc_pN_pS; # comment this out to speed up
    }
}
if ((@samples > 1) and !$pi_only) {
    print"\n### Calculating Fixation Index (FST) ###\n";
    &calc_fst;
    if ($gff_file) {
        print"\n### Calculating Gene-wise Fixation Index (FST) ###\n";
        &calc_per_gene_fst;
        if ($n_fst_permutations) {
            print"\n### Calculating Permuted Gene-wise Fixation Index (FST) ###\n";
            &calc_per_gene_fst_permuted;
        }
        if ($genetic_code_file) {
            print"\n### Calculating Gene-wise Aminoacid Fixation Index (aa-FST) ###\n";
            &calc_per_gene_aminoacid_fst;
            #print"\n### Calculating Gene-wise Neutrality Index (NI) ###\n";
            #&calc_NI;
        }
    }
}
print"\n### Printing results to files ###\n";
&print_output_to_file;
print"\n### Finished pogenom succesfully ###\n\n";

####################

sub read_loci_to_include {
    open (INFILE, "$loci_file") || die ("Error: can't open $loci_file");
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        @fields = split(/\t/, $row);
        $locus = $fields[0]."|".$fields[1];
	print("$locus\n");
        $include_locus{$locus} = 1;
    }
    $temp = keys %include_locus;
    $logtext = $logtext."Number of loci specified in file: $temp\n";
}

sub read_samples_to_include {
    open (INFILE, "$sample_file") || die ("Error: can't open $sample_file");
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        $include_sample{$row} = 1;
    }
    $temp = keys %include_sample;
    $logtext = $logtext."Number of samples specified in file: $temp\n";
}

sub read_gff {
    $fasta_started = 0;
    $seq = "";
    open (INFILE, "$gff_file") || die ("Error: can't open $gff_file");
    print"Reading $gff_file\n";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if ($fasta_started == 1) {
            if (substr($row, 0, 1) eq ">") {
                if ($seq ne "") {
                    $contig_seq{$contig} = $seq;
                }
                $contig = $row;
                substr($contig, 0, 1) = "";
                #die;
                $seq = "";
            } else {
                $seq = $seq.$_;
                $contig_seq{$contig} = uc $seq;
            }
        }
        elsif (substr($row, 0, 7) eq "##FASTA") {
            $fasta_started = 1;
        } else {
            next if (substr($row, 0, 1) eq "#");
            @fields = split(/\t/, $row);
            next if (@fields == 1);
            next if ($fields[2] ne "CDS");
            $contig = $fields[0];
            $start = $fields[3];
            $end = $fields[4];
            $strand = $fields[6];
            $end = $end - 3 if ($strand eq "+"); # to exclude the stop codon
            $start = $start + 3 if ($strand eq "-"); # to exclude the stop codon
            $annotation = $fields[8];
            @subfields = split(/;/, $annotation);
            $gene = $subfields[0];
            if (defined $gene_start{$gene}) {
                print"Error: Non-unique gene identifiers in GFF file. The program will exit without finishing.\n\n"; exit;
            }
            $gene_start{$gene} = $start;
            $gene_end{$gene} = $end;
            $gene_strand{$gene} = $strand;
            $gene_length{$gene} = $gene_end{$gene} - $gene_start{$gene} + 1;
            $gene_contig{$gene} = $contig;
            if (!defined $contig_genes{$contig}) {
                push(@contigs, $contig);
            }
            push(@{ $contig_genes{$contig} }, $gene);
            for ($pos = $start; $pos < $end; $pos++) {
                $locus = $contig."|".$pos;
                if (defined $locus_found{$locus}) {
                    if ($locus_found{$locus} >= $min_found_in) {
                        $gene_locus{$gene}{$locus} = 1;
                    }
                }
            }
        }
    }
    close (INFILE);
    if ($fasta_started == 1) {
        $temp_genome_size = 0;
        foreach $contig (@contigs) {
            if (!defined $contig_seq{$contig}) {
                print"\nError: Missmatch between contig id ($contig) in upper and lower part of gff file\n\n"; exit;
            }
            $temp_genome_size = $temp_genome_size + length($contig_seq{$contig});
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                if ($gene_strand{$gene} eq "+") {
                    $gene_seq{$gene} = substr($contig_seq{$contig}, ($gene_start{$gene} - 1), $gene_length{$gene});
                    #$peptide = "";
                    #for ($j = 0; $j < $gene_length{$gene}; $j=$j+3) {
                    #    $codon = substr($gene_seq{$gene}, $j, 3);
                    #    $peptide = $peptide.$codon_aminoacid{$codon};
                    #}
                    #print">$gene\n$gene_seq{$gene}\n";
                    #print">$gene\n$peptide\n";
                }
                if ($gene_strand{$gene} eq "-") {
                    $gene_seq{$gene} = substr($contig_seq{$contig}, ($gene_start{$gene} - 1), $gene_length{$gene});
                    $gene_seq{$gene} = &make_revcomp($gene_seq{$gene});
                    #$peptide = "";
                    #for ($j = 0; $j < $gene_length{$gene}; $j=$j+3) {
                    #    $codon = substr($gene_seq{$gene}, $j, 3);
                    #    $peptide = $peptide.$codon_aminoacid{$codon};
                    #}
                    #print">$gene\n$gene_seq{$gene}\n";
                    #print">$gene\n$peptide\n";
                }
            }
        }
    }
    if (!$genome_size and !$fasta_file) {
        if ($fasta_started == 1) {
            $genome_size = $temp_genome_size;
            print"Genome size calculated from GFF to $genome_size bp\n";
            $logtext = $logtext."Genome size calculated from GFF to $genome_size bp\n";
        } else {
            print"Error: Genome size could not be calculated from GFF file\n\n"; exit;
        }
    }
}

sub read_fasta {
    %contig_seq = ();
    $seq = "";
    open (INFILE, "$fasta_file") || die ("Error: can't open $fasta_file");
    print"Reading $fasta_file\n";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        if (substr($row, 0, 1) eq ">") {
            if ($seq ne "") {
                $contig_seq{$contig} = $seq;
            }
            $contig = $row;
            if ($split_fasta_header) {
                @subfields = split(/\s+/, $row);
                $contig = $subfields[0];
            }
            substr($contig, 0, 1) = "";
            $seq = "";
        } else {
            $seq = $seq.$_;
            $contig_seq{$contig} = uc $seq;
        }
    }
    if (!$genome_size) {
        $genome_size = 0;
        foreach $contig (keys %contig_seq) {
            $genome_size = $genome_size + length($contig_seq{$contig});
        }
    }
    print"Genome size calculated from fasta file to $genome_size bp\n";
    $logtext = $logtext."Genome size calculated from fasta file to $genome_size bp\n";
    if (defined $gff_file) {
        foreach $contig (@contigs) {
            if (!defined $contig_seq{$contig}) {
                print"\nError: Missmatch between contig id in gff ($contig) and fasta file\n\n"; exit;
            }
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                if ($gene_strand{$gene} eq "+") {
                    $gene_seq{$gene} = substr($contig_seq{$contig}, ($gene_start{$gene} - 1), $gene_length{$gene});
                }
                if ($gene_strand{$gene} eq "-") {
                    $gene_seq{$gene} = substr($contig_seq{$contig}, ($gene_start{$gene} - 1), $gene_length{$gene});
                    $gene_seq{$gene} = &make_revcomp($gene_seq{$gene});
                }
            }
        }
    }
}

sub read_genetic_code {
    open(INFILE, "$genetic_code_file") || die ("Error: can't open $genetic_code_file");
    print"Reading $genetic_code_file\n";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        @fields = split(/\t/);
        $codon_aminoacid{$fields[0]} = $fields[1];
    }
    close (INFILE);
}

sub get_snp_data_combined_vcf {
    local($ok_sample);
    %nt = ('A', 1, 'T', 1, 'C', 1, 'G', 1);
    @samples = ();
    @samples_plus = ();
    open (INFILE, "$vcf_file") || die ("Error: can't open $vcf_file");
    print"Reading $vcf_file\n";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        @fields = split(/\t/, $row);
        if (substr($row, 0, 6) eq "#CHROM") {
            for ($i = 9; $i < @fields; $i++) {
                $sample = $fields[$i];
                $samples[$i - 9] = $sample;
                $include_sample{$sample} = 1 if (!$sample_file);
                $ok_sample++ if (defined $include_sample{$sample});
            }
            @samples_plus = (@samples, 'All_samples_combined');
            $min_found_in = $ok_sample if ($min_found_in == 0);
        }
        next if (substr($row, 0, 1) eq "#");
        #next if (@fields == 1);
        $contig = $fields[0];
        $pos = $fields[1];
        $locus = $contig."|".$pos;
        #print"$locus $locus_found{$locus} $min_found_in\n";
        if ($loci_file) {
            #print"$locus\n";
            next if (!defined $include_locus{$locus});
        }
        @subfields = split(/:/, $fields[8]);
        $nformat_fields = @subfields;
        if ($vcf_format eq "GATK") {
            $count_ix = undef;
            for ($i = 0; $i < @subfields; $i++) {
                $count_ix = $i if ($subfields[$i] eq "AD");
            }
            if (!$count_ix) {
                print"\nError: Field AD is lacking in the Format column of the VCF file at CHROM: $contig POS: $pos\n\n"; exit;
            }
        } elsif ($vcf_format eq "freebayes") {
            $tot_count_ix = $ref_count_ix = $alt_count_ix = undef;
            for ($i = 0; $i < @subfields; $i++) {
                #$tot_count_ix = $i if ($subfields[$i] eq "DP");
                $ref_count_ix = $i if ($subfields[$i] eq "RO");
                $alt_count_ix = $i if ($subfields[$i] eq "AO");
            }
            if (!$ref_count_ix or !$alt_count_ix) {
                print"\nError: Either field RO or AO is lacking in the Format column of the VCF file at CHROM: $contig POS: $pos\n\n"; exit;
            }
        }
        $locus_found{$locus} = 0;
        for ($i = 9; $i < @fields; $i++) {
            @subfields = split(/:/, $fields[$i]);
            $sample = $samples[$i - 9];
            next if (!defined $include_sample{$sample});
            if (@subfields == $nformat_fields) {
                if ($vcf_format eq "GATK") {
                    @allele_count = split(/,/, $subfields[$count_ix]);
                } elsif ($vcf_format eq "freebayes") {
                    @allele_count = split(/,/, $subfields[$alt_count_ix]);
                    unshift(@allele_count, $subfields[$ref_count_ix]);
                }
                if ($allele_count[0] ne ".") {
                    $tot_count = 0;
                    for ($j = 0; $j < @allele_count; $j++) {
                        $tot_count = $tot_count + $allele_count[$j];
                    }
                    if ($tot_count >= $min_count) {
                        $locus_found{$locus}++;
                        #print"$locus_found{$locus}\n";
                        $sample_foundlocus{$sample}{$locus} = 1; # this is for estimating genome coverage
                    }
                }
            }
        }
        next if ($locus_found{$locus} < $min_found_in);
        $contig_pos{$contig}{$pos} = 1;
        $ref = $fields[3];
        @alleles = split(/,/, $fields[4]);
        unshift(@alleles, $ref);
        $sample_locus_totcount{'All_samples_combined'}{$locus} = 0;
        for ($i = 9; $i < @fields; $i++) {
            $sample = $samples[$i - 9];
            next if (!defined $include_sample{$sample});
            @subfields = split(/:/, $fields[$i]);
            next if (@subfields != $nformat_fields);
            if ($vcf_format eq "GATK") {
                @allele_count = split(/,/, $subfields[$count_ix]);
            } elsif ($vcf_format eq "freebayes") {
                @allele_count = split(/,/, $subfields[$alt_count_ix]);
                unshift(@allele_count, $subfields[$ref_count_ix])
            }
            $tot_count = 0;
            if ($allele_count[0] ne ".") {
                for ($j = 0; $j < @allele_count; $j++) {
                    $tot_count = $tot_count + $allele_count[$j];
                }
            }
            next if ($allele_count[0] eq ".");
            next if ($tot_count < $min_count);
            $sample_locus_ref{$sample}{$locus} = $ref;
            $sample_locus_totcount{$sample}{$locus} = $tot_count;
            $sample_locus_totcount{'All_samples_combined'}{$locus} = $sample_locus_totcount{'All_samples_combined'}{$locus} + $tot_count;
            #print"Ref: #$ref# $ref_count\n";
            for ($j = 0; $j < @alleles; $j++) {
                $sample_locus_allel_counts{$sample}{$locus}{$alleles[$j]} = $allele_count[$j];
                $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alleles[$j]} = 0 if (!defined $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alleles[$j]});
                $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alleles[$j]} = $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$alleles[$j]} + $allele_count[$j];
            }
        }
        foreach $string (@alleles) {
            $contig_pos{$contig}{$pos} = 0 if (!defined $nt{$string});
        }
    }
    foreach $locus (keys %locus_found) {
        $found[$locus_found{$locus}]++;
    }
    $i = (@found - 1);
    $cumulative_found = $found[$i];
    print"Number of loci found $i times: $cumulative_found\n";
    for ($i = (@found - 2); $i > 0; $i--) {
        $cumulative_found = $cumulative_found + $found[$i];
        print"Number of loci found >= $i times: $cumulative_found\n";
    }
    if ((keys %sample_locus_allel_counts) == 0) {
        print"Zero loci found fulfilling criteria\n";
        print"No files will be generated\n";
        exit;
    }
    # removing unwanted samples from @samples
    local(@include);
    for ($i = 0; $i < @samples; $i++) {
        push(@include, $i) if (defined $include_sample{$samples[$i]});
    }
    @samples = @samples[@include];
    @samples_plus = (@samples, 'All_samples_combined');
}

sub get_snp_data_combined_vcf_split_haplotypes {
    local($ok_sample);
    %nt = ('A', 1, 'T', 1, 'C', 1, 'G', 1);
    @samples = ();
    @samples_plus = ();
    open (INFILE, "$vcf_file") || die ("Error: can't open $vcf_file");
    print"Reading $vcf_file\n";
    print"Haplotypes will be split into individual bases\n";
    while (<INFILE>) {
        $_ =~ s/\R//g;
        $row = $_;
        @fields = split(/\t/, $row);
        if (substr($row, 0, 6) eq "#CHROM") {
            for ($i = 9; $i < @fields; $i++) {
                $sample = $fields[$i];
                $samples[$i - 9] = $sample;
                $include_sample{$sample} = 1 if (!$sample_file);
                $ok_sample++ if (defined $include_sample{$sample});
            }
            @samples_plus = (@samples, 'All_samples_combined');
            $min_found_in = $ok_sample if ($min_found_in == 0);
        }
        next if (substr($row, 0, 1) eq "#");
        $contig = $fields[0];
        $pos = $fields[1];
        $ref = $fields[3];
        @alleles = split(/,/, $fields[4]);
        unshift(@alleles, $ref);
        $indel = undef;
        for ($i = 1; $i < @alleles; $i++) {
            if (length($alleles[$i]) != length($ref)) {
                $indel = 1;
            }
        }
        next if $indel;
        #print"Ref: $ref Alt: @alt\n";
        @subfields = split(/:/, $fields[8]);
        $nformat_fields = @subfields;
        if ($vcf_format eq "GATK") {
            $count_ix = undef;
            for ($i = 0; $j < @subfields; $i++) {
                $count_ix = $i if ($subfields[$i] eq "AD");
            }
            if (!$count_ix) {
                print"\nError: Field AD is lacking in the Format column of the VCF file at CHROM: $contig POS: $pos\n\n"; exit;
            }
        } elsif ($vcf_format eq "freebayes") {
            $tot_count_ix = $ref_count_ix = $alt_count_ix = undef;
            for ($i = 0; $i < @subfields; $i++) {
                #$tot_count_ix = $i if ($subfields[$i] eq "DP");
                $ref_count_ix = $i if ($subfields[$i] eq "RO");
                $alt_count_ix = $i if ($subfields[$i] eq "AO");
            }
            if (!$ref_count_ix or !$alt_count_ix) {
                print"\nError: Either field RO or AO is lacking in the Format column of the VCF file at CHROM: $contig POS: $pos\n\n"; exit;
            }
        }
        
        for ($i = 0; $i < length($ref); $i++) {
            $modpos = $pos + $i;
            $locus = $contig."|".$modpos;
            if ($loci_file) {
                next if (!defined $include_locus{$locus});
            }
            $variant = 0;
            for ($j = 0; $j < @alleles; $j++) {
                #print"$i $j $alleles[$j] $ref\n";
                if (substr($alleles[$j], $i, 1) ne substr($ref, $i, 1)) {
                    $variant = 1;
                }
            }
            #print"variant: $variant\n\n";
            next if ($variant == 0); # for excluding invariant loci (common within haplotypes)
            $locus_found{$locus} = 0;
            for ($j = 9; $j < @fields; $j++) {
                $sample = $samples[$j - 9];
                next if (!defined $include_sample{$sample});
                @subfields = split(/:/, $fields[$j]);
                if (@subfields == $nformat_fields) {
                    if ($vcf_format eq "GATK") {
                        @allele_count = split(/,/, $subfields[$count_ix]);
                    } elsif ($vcf_format eq "freebayes") {
                        @allele_count = split(/,/, $subfields[$alt_count_ix]);
                        unshift(@allele_count, $subfields[$ref_count_ix]);
                    }
                    if ($allele_count[0] ne ".") {
                        $tot_count = 0;
                        for ($k = 0; $k < @allele_count; $k++) {
                            $tot_count = $tot_count + $allele_count[$k];
                        }
                        if ($tot_count >= $min_count) {
                            $locus_found{$locus}++;
                            $sample_foundlocus{$sample}{$locus} = 1; # this is for estimating genome coverage
                        }
                    }
                }
            }
            next if ($locus_found{$locus} < $min_found_in);
            $contig_pos{$contig}{$modpos} = 1;
            $subref = substr($ref, $i, 1);
            $contig_pos{$contig}{$modpos} = 0 if (!defined $nt{$subref});
            $sample_locus_totcount{'All_samples_combined'}{$locus} = 0;
            for ($j = 9; $j < @fields; $j++) {
                $sample = $samples[$j - 9];
                @subfields = split(/:/, $fields[$j]);
                next if (@subfields != $nformat_fields);
                if ($vcf_format eq "GATK") {
                    @allele_count = split(/,/, $subfields[$count_ix]);
                } elsif ($vcf_format eq "freebayes") {
                    @allele_count = split(/,/, $subfields[$alt_count_ix]);
                    unshift(@allele_count, $subfields[$ref_count_ix])
                }
                $tot_count = 0;
                if ($allele_count[0] ne ".") {
                    for ($k = 0; $k < @allele_count; $k++) {
                        $tot_count = $tot_count + $allele_count[$k];
                    }
                }
                next if ($allele_count[0] eq ".");
                next if ($tot_count < $min_count);
                $sample_locus_ref{$sample}{$locus} = $subref;
                $sample_locus_totcount{$sample}{$locus} = $tot_count;
                $sample_locus_totcount{'All_samples_combined'}{$locus} = $sample_locus_totcount{'All_samples_combined'}{$locus} + $tot_count;
                #print"    Sample: $sample Subpos: $i Modpos: $modpos\n";
                #print"        Subref: $subref $ref_count\n";
                for ($k = 0; $k < @alleles; $k++) {
                    $subal = substr($alleles[$k], $i, 1);
                    if (defined $sample_locus_allel_counts{$sample}{$locus}{$subal}) {
                        $sample_locus_allel_counts{$sample}{$locus}{$subal} = $sample_locus_allel_counts{$sample}{$locus}{$subal} + $allele_count[$k];
                    } else {
                        $sample_locus_allel_counts{$sample}{$locus}{$subal} = $allele_count[$k];
                    }
                    $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subal} = 0 if (!defined $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subal});
                    $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subal} = $sample_locus_allel_counts{'All_samples_combined'}{$locus}{$subal} + $allele_count[$k];
                    $contig_pos{$contig}{$modpos} = 0 if (!defined $nt{$subal});
                    #print"        Subal: $subal $allele_count[$k] $sample_locus_allel_counts{$sample}{$locus}{$subal}\n";
                }
                
            }
        }
    }
    foreach $locus (keys %locus_found) {
        $found[$locus_found{$locus}]++;
    }
    $i = (@found - 1);
    $cumulative_found = $found[$i];
    print"Number of loci found $i times: $cumulative_found\n";
    for ($i = (@found - 2); $i > 0; $i--) {
        $cumulative_found = $cumulative_found + $found[$i] if (defined $found[$i]);
        print"Number of loci found >= $i times: $cumulative_found\n";
    }
    if ((keys %sample_locus_allel_counts) == 0) {
        print"Zero loci found fulfilling criteria\n";
        print"No files will be generated\n";
        exit;
    }
    # removing unwanted samples from @samples
    local(@include);
    for ($i = 0; $i < @samples; $i++) {
        push(@include, $i) if (defined $include_sample{$samples[$i]});
    }
    @samples = @samples[@include];
    @samples_plus = (@samples, 'All_samples_combined');
}


sub subsample_allele_counts {
    foreach $sample (@samples_plus) {
        foreach $locus (keys %{$sample_locus_allel_counts{$sample}}) {
            if ($sample_locus_totcount{$sample}{$locus} > $subsample) {
                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                @temp = ();
                foreach $allele (@alleles) {
                    $count = $sample_locus_allel_counts{$sample}{$locus}{$allele};
                    for ($i = 0; $i < $count; $i++) {
                        push(@temp, $allele);
                    }
                    $sample_locus_allel_counts{$sample}{$locus}{$allele} = 0;
                }
                for ($i = 0; $i < $subsample; $i++) {
                    $ix = int(rand(@temp));
                    $sample_locus_allel_counts{$sample}{$locus}{$temp[$ix]}++;
                }
                $sample_locus_totcount{$sample}{$locus} = $subsample;
            }
        }
    }
}

sub subsample_allele_counts_dupl {
    foreach $sample (@samples_plus) {
        $sample_dupl = $sample."_duplicate";
        push(@samples_plus_dupl, $sample);
        push(@samples_plus_dupl, $sample_dupl);
        foreach $locus (keys %{$sample_locus_allel_counts{$sample}}) {
            if ($sample_locus_totcount{$sample}{$locus} > $subsample) {
                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                @temp = ();
                foreach $allele (@alleles) {
                    $count = $sample_locus_allel_counts{$sample}{$locus}{$allele};
                    for ($i = 0; $i < $count; $i++) {
                        push(@temp, $allele);
                    }
                    $sample_locus_allel_counts{$sample}{$locus}{$allele} = 0;
                    $sample_locus_allel_counts{$sample_dupl}{$locus}{$allele} = 0;
                }
                
                for ($i = 0; $i < $subsample; $i++) {
                    $ix = int(rand(@temp));
                    $sample_locus_allel_counts{$sample}{$locus}{$temp[$ix]}++;
                }
                for ($i = 0; $i < $subsample; $i++) {
                    $ix = int(rand(@temp));
                    $sample_locus_allel_counts{$sample_dupl}{$locus}{$temp[$ix]}++;
                }
                $sample_locus_totcount{$sample}{$locus} = $subsample;
                $sample_locus_totcount{$sample_dupl}{$locus} = $subsample;
            }
        }
    }
    @samples = @samples_plus = @samples_plus_dupl;
}

sub estimate_genome_coverage {
    foreach $sample (@samples) {
        @loci = (keys %{$sample_locus_allel_counts{$sample}});
        if (@loci == 0) {
            return; # since estimating genome coverage does not work if any of the samples has 0 loci
        }
    }
    local($shared);
    $est_num_loci = 0;
    foreach $sample (@samples) {
        $sample_estcov{$sample} = 0;
        foreach $sample2 (@samples) {
            next if ($sample2 eq $sample);
            $shared = 0;
            foreach $locus (keys %{$sample_foundlocus{$sample2}}) {
                if (defined ($sample_foundlocus{$sample}{$locus})) {
                    $shared++;
                }
            }
            $sample_estcov{$sample} = $sample_estcov{$sample} + $shared/(keys %{$sample_foundlocus{$sample2}});
            #local($temp) = $shared/(keys %{$sample_foundlocus{$sample2}});
            #print " $sample $temp\n";
        }
        $sample_estcov{$sample} = $sample_estcov{$sample}/(@samples - 1); # this is the estimated proportion of the genome with coverage fulfilling the -min_count criteria
        $est_num_loci = $est_num_loci + (keys %{$sample_foundlocus{$sample}})/$sample_estcov{$sample};
        #local($temp) = (keys %{$sample_foundlocus{$sample}})/$sample_estcov{$sample};
        #print "##$sample $num_loci{$sample} $sample_estcov{$sample} $temp\n";
    }
    $est_num_loci = $est_num_loci/@samples; # this is the estimated total number of loci with variation for the pool of samples, after splitting up haplotypes
}

sub estimate_genome_coverage_old {
    local($shared);
    foreach $sample (@samples) {
        $sample_estcov{$sample} = 0;
        foreach $sample2 (@samples) {
            next if ($sample2 eq $sample);
            $shared = 0;
            foreach $locus (keys %{$sample_locus_totcount{$sample2}}) {
                if (defined ($sample_locus_totcount{$sample}{$locus})) {
                    $shared++;
                }
            }
            $sample_estcov{$sample} = $sample_estcov{$sample} + $shared/(keys %{$sample_locus_totcount{$sample2}});
            #local($temp) = $shared/(keys %{$sample_locus_totcount{$sample2}});
            #print " $sample $temp\n";
        }
        $sample_estcov{$sample} = $sample_estcov{$sample}/(@samples - 1);
        #print "#$sample $sample_estcov{$sample}\n\n";
    }
}

sub calc_pi {
    print "Sample\tpi\tTotal_num_loci\tTotal_num_alleles\tAverage_depth\n";
    foreach $sample (@samples_plus) {
        $intra_pi = 0;
        $tot_alleles = 0;
        $av_count = 0;
        @loci = (keys %{$sample_locus_allel_counts{$sample}});
        foreach $locus (@loci) {
            $tot_count = $sample_locus_totcount{$sample}{$locus};
            #print"Locus totcount $locus $tot_count\n";
            @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
            #print"Alleles: @alleles\n";
            $av_count = $av_count + $tot_count;
            $locus_intra_pi = 0;
            for ($i = 0; $i < @alleles; $i++) {
                #print"Allele 1: $alleles[$i]\n";
                $counts_1 = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$i]};
                #print"Counts1: $counts_1\n";
                $tot_alleles++ if ($counts_1 > 0);
                for ($j = 0; $j < @alleles; $j++) {
                    next if ($alleles[$i] eq $alleles[$j]);
                    #print"Allele 2: $alleles[$j]\n";
                    $counts_2 = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$j]};
                    #print"Counts2: $counts_2\n";
                    $locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count - 1)); # According to Schloising
                    #$locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count)); # According to our initial version
                }
            }
            #print"$locus $locus_intra_pi\n\n";
            $sample_locus_pi{$sample}{$locus} = $locus_intra_pi;
            $intra_pi = $intra_pi + $locus_intra_pi;
        }
        if (@loci > 0) {
            $intra_pi = $intra_pi/$genome_size;
            $av_count = $av_count/@loci;
        } else {
            $av_count = "NA";
            $intra_pi = "NA";
            $tot_alleles = "NA";
        }
        $sample_pi{$sample} = $intra_pi;
        $sample_avcount{$sample} = $av_count;
        $sample_totalleles{$sample} = $tot_alleles;
        $num_loci{$sample} = @loci;
        print "$sample\t$intra_pi\t$num_loci{$sample}\t$tot_alleles\t$av_count\n";
    }
}

sub calc_per_gene_pi {
    foreach $contig (@contigs) {
        #print "$contig\n";
        @genes = @{ $contig_genes{$contig} };
        #print "@genes\n"; die;
        foreach $gene (@genes) {
            foreach $sample (@samples_plus) {
                $intra_pi = 0;
                $missing_loci = 0;
                foreach $locus (keys %{$gene_locus{$gene}}) {
                    #print "$locus/n";
                    $missing_loci = 1 if (!defined $sample_locus_pi{$sample}{$locus});
                    next if (!defined $sample_locus_pi{$sample}{$locus});
                    $intra_pi = $intra_pi + $sample_locus_pi{$sample}{$locus};
                }
                $intra_pi = $intra_pi/$gene_length{$gene};
                $sample_gene_pi{$sample}{$gene} = $intra_pi;
                if ($na_if_missing_loci == 1 && $missing_loci > 0) {
                    $sample_gene_pi{$sample}{$gene} = "NA";
                }
            }
        }
    }
}

sub calc_per_gene_aminoacid_pi {
    foreach $sample (@samples_plus) {
        $tot_peptides = 0;
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $intra_pi = 0;
                $missing_loci = 0;
                #print"\n$gene\n";
                foreach $locus (keys %{$gene_locus{$gene}}) {
                    $missing_loci = 1 if (!defined $sample_locus_pi{$sample}{$locus});
                    next if (!defined $sample_locus_pi{$sample}{$locus});
                    @fields = split(/\|/, $locus);
                    $contig = $fields[0];
                    $pos = $fields[1];
                    $tot_count = $sample_locus_totcount{$sample}{$locus};
                    @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                    #print"Alleles: @alleles\n";
                    for ($i = 0; $i < @alleles; $i++) {
                        $mod_contig_seq = $contig_seq{$contig};
                        substr($mod_contig_seq, ($pos-1), length($alleles[$i])) = $alleles[$i];
                        $mod_gene_seq = substr($mod_contig_seq, ($gene_start{$gene} - 1), $gene_length{$gene}); # OBS!!!
                        $mod_gene_seq = &make_revcomp($mod_gene_seq) if ($gene_strand{$gene} eq "-");
                        $peptide = "";
                        for ($j = 0; $j < $gene_length{$gene}; $j=$j+3) {
                            $codon = substr($mod_gene_seq, $j, 3);
                            $peptide = $peptide.$codon_aminoacid{$codon};
                        }
                        if (defined $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide}) {
                            $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide} = $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide} + $sample_locus_allel_counts{$sample}{$locus}{$alleles[$i]};
                        } else {
                            $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptide} = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$i]};
                        }
                    }
                    $locus_intra_pi = 0;
                    @peptides = (keys %{$sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}});
                    for ($i = 0; $i < @peptides; $i++) {
                        #print"\nPeptide 1: $peptides[$i]\n";
                        $counts_1 = $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptides[$i]};
                        #print"Counts1: $counts_1\n";
                        $tot_peptides++ if ($counts_1 > 0);
                        for ($j = 0; $j < @peptides; $j++) {
                            next if ($peptides[$i] eq $peptides[$j]);
                            #print"Peptide 2: $peptides[$j]\n";
                            $counts_2 = $sample_locus_gene_peptide_counts{$sample}{$locus}{$gene}{$peptides[$j]};
                            #print"Counts2: $counts_2\n";
                            $locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count - 1)); # According to Schloising
                            #$locus_intra_pi = $locus_intra_pi + ($counts_1/$tot_count)*($counts_2/($tot_count)); # According to our initial version
                        }
                    }
                    #print"Sample: $sample $locus locus_intra_pi: $locus_intra_pi\n";
                    $sample_locus_aminoacid_pi{$sample}{$locus} = $locus_intra_pi;
                    $intra_pi = $intra_pi + $locus_intra_pi;
                    #print"Sample: $sample $locus intra_pi: $intra_pi\n\n";
                }
                $intra_pi = $intra_pi/$gene_length{$gene};
                $sample_gene_aminoacid_pi{$sample}{$gene} = $intra_pi;
                #$sample_totpeptides{$sample} = $tot_peptides;
                #print"Sample: $sample $locus intra_pi: $intra_pi\n\n";
                if ($na_if_missing_loci == 1 && $missing_loci == 1) {
                    $sample_gene_aminoacid_pi{$sample}{$gene} = "NA";
                }
            }
        }
    }
}

sub calc_aminoacid_frequencies {
    # calculate mean count (avrunda till närmaste heltal?) för varje aa, baserat på alla locus som överlappar en codon/aminosyra-position.
    # först ta reda på vilka loci som överlappar varje position.
    # skippa loci som har alleler av varierande längd (kan bli frameshifts)
    foreach $contig (@contigs) {
        @genes = @{ $contig_genes{$contig} };
        foreach $gene (@genes) {
            #next if ($gene ne "ID=PROKKA_MOD_00014");
            #print"$gene\n";
            $contig_seq = $contig_seq{$contig};
            if ($gene_strand{$gene} eq "+") {
                $gene_seq = substr($contig_seq, ($gene_start{$gene} - 1), $gene_length{$gene});
                #print"$gene_seq\n";
                @loci = (keys %{$gene_locus{$gene}});
                next if (@loci == 0);
                %codonposition_n_locus = ();
                #print"loci @loci\n";
                foreach $locus (@loci) {
                    next if ($sample_locus_totcount{'All_samples_combined'}{$locus} == 0);
                    @alleles = (keys %{$sample_locus_allel_counts{'All_samples_combined'}{$locus}});
                    @fields = split(/\|/, $locus);
                    #print"fields1 @fields\n";
                    $pos = $fields[1];
                    %lengths = ();
                    for ($k = 0; $k < @alleles; $k++) {
                        $length = length($alleles[$k]);
                        $lengths{$length} = 1;
                    }
                    next if ((keys %lengths) > 1); # skip loci that have alleles that differ in length
                    $gene_pos = $pos - $gene_start{$gene}; # 0 equals first nucleotide of first (i.e. start) codon
                    $first_codon_pos = int($gene_pos/3); # 0 equals first codon (start codon)
                    $last_codon_pos = int(($gene_pos + $length - 1)/3); # 0 equals first codon (start codon)
                    for ($k = $first_codon_pos; $k <= $last_codon_pos; $k++) {
                        $codonposition_n_locus{$k}{$locus} = 1;
                        #print"pos locus $k $locus\n";
                    }
                }
                foreach $codon_pos (keys %codonposition_n_locus) {
                    %sample_n_includedlocus = ();
                    foreach $locus (keys %{$codonposition_n_locus{$codon_pos}}) {
                        @fields = split(/\|/, $locus);
                        $pos = $fields[1];
                        #print"Pos $pos\n";
                        $gene_pos = $pos - $gene_start{$gene};
                        #print"Gene_pos $gene_pos\n";
                        @alleles = (keys %{$sample_locus_allel_counts{'All_samples_combined'}{$locus}});
                        for ($k = 0; $k < @alleles; $k++) {
                            $mod_gene_seq = $gene_seq;
                            substr($mod_gene_seq, $gene_pos, length($alleles[$k])) = $alleles[$k];
                            $mod_codon = substr($mod_gene_seq, $codon_pos*3, 3);
                            $aminoacid = $codon_aminoacid{$mod_codon};
                            foreach $sample (@samples_plus) {
                                if (defined $sample_locus_totcount{$sample}{$locus}) {
                                    $sample_n_includedlocus{$sample}{$locus} = 1;
                                    $count = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                    if (defined $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid}) {
                                        $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid} = $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid} + $count;
                                    } else {
                                        $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid} = $count;
                                    }
                                }
                            }
                        }
                    }
                    foreach $sample (@samples_plus) {
                        if (defined $sample_n_includedlocus{$sample}) {
                            $num = (keys %{$sample_n_includedlocus{$sample}});
                            foreach $aminoacid (keys %{$sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}}) {
                                $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid} = $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid}/$num;
                            }
                        }
                    }
                }
            }
        }
    }
}

sub calc_fst {
    print "Sample1\tSample2\tpi_1\tpi_2\tpi_1-2\tfst\n";
    for ($ix1 = 0; $ix1 < @samples; $ix1++) {
        $sample1 = $samples[$ix1];
        #print"Sample1 $sample1\n";
        for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
        #for ($ix2 = $ix1 + 0; $ix2 < @samples; $ix2++) {
            $sample2 = $samples[$ix2];
            #print"Sample2 $sample2\n";
            $inter_pi = 0;
            $sample1_pi = $sample2_pi = 0; # intra-pi based on only loci shared by both samples
            @loci = (keys %{$sample_locus_allel_counts{$sample1}});
            #$temp = @loci;
            #print"sample1: $sample1 num_loci: $temp\n"; die;
            #print"sample1: $sample1 loci: @loci\n"; die;
            foreach $locus (@loci) {
                next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                #print"Locus $locus\n";
                @alleles1 = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                @alleles2 = (keys %{$sample_locus_allel_counts{$sample2}{$locus}});
                $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                #print"Alleles: @alleles1 @alleles2\n";
                for ($i = 0; $i < @alleles1; $i++) {
                    #print"Allele 1: $alleles1[$i]\n";
                    $counts_1 = $sample_locus_allel_counts{$sample1}{$locus}{$alleles1[$i]};
                    #print"Counts1: $counts_1\n";
                    for ($j = 0; $j < @alleles2; $j++) {
                        next if ($alleles1[$i] eq $alleles2[$j]);
                        #print"Allele 2: $alleles2[$j]\n";
                        $counts_2 = $sample_locus_allel_counts{$sample2}{$locus}{$alleles2[$j]};
                        #print"Counts2: $counts_2\n";
                        $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                    }
                }
                $sample1_pi = $sample1_pi + $sample_locus_pi{$sample1}{$locus};
                $sample2_pi = $sample2_pi + $sample_locus_pi{$sample2}{$locus};
                
                #$sample1_pi = $sample1_pi;
                #$sample2_pi = $sample2_pi;
            }
            $inter_pi = $inter_pi/$genome_size;
            $sample_sample_pi{$sample1}{$sample2} = $inter_pi;
            $sample_sample_pi{$sample2}{$sample1} = $inter_pi;
            $sample1_pi = $sample1_pi/$genome_size;
            $sample2_pi = $sample2_pi/$genome_size;
            if ($inter_pi > 0) {
                #$fst = 1 - 0.5*($sample_pi{$sample1} + $sample_pi{$sample2})/$inter_pi; # intra pi values based on all loci
                $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi; # intra pi values only based on shared loci
                $fst = sprintf("%.4f", $fst);
            } else {
                $fst = "NA";
            }
            $sample_sample_fst{$sample1}{$sample2} = $fst;
            $sample_sample_fst{$sample2}{$sample1} = $fst;
            print"$sample1\t$sample2\t$sample_pi{$sample1}\t$sample_pi{$sample2}\t$inter_pi\t$fst\n";
        }
    }
}

sub calc_per_gene_fst {
    foreach $contig (@contigs) {
        #print "Contig $contig\n";
        @genes = @{ $contig_genes{$contig} };
        foreach $gene (@genes) {
            #print "Gene $gene\n";
            @loci = (keys %{$gene_locus{$gene}});
            next if (@loci == 0);
            for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                $sample1 = $samples[$ix1];
                #print"Sample1 $sample1\n";
                for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
                    $sample2 = $samples[$ix2];
                    #print"Sample2 $sample2\n";
                    $inter_pi = 0;
                    $sample1_pi = $sample2_pi = 0;
                    foreach $locus (@loci) {
                        #print"Locus $locus\n";
                        next if (!defined $sample_locus_allel_counts{$sample1}{$locus});
                        next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                        @alleles1 = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                        @alleles2 = (keys %{$sample_locus_allel_counts{$sample2}{$locus}});
                        $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                        $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                        #print"Alleles1: @alleles1 Alleles2: @alleles2\n";
                        for ($i = 0; $i < @alleles1; $i++) {
                            #print"Allele 1: $alleles1[$i]\n";
                            $counts_1 = $sample_locus_allel_counts{$sample1}{$locus}{$alleles1[$i]};
                            #print"Counts1: $counts_1\n";
                            for ($j = 0; $j < @alleles2; $j++) {
                                next if ($alleles1[$i] eq $alleles2[$j]);
                                #print"Allele 2: $alleles2[$j]\n";
                                $counts_2 = $sample_locus_allel_counts{$sample2}{$locus}{$alleles2[$j]};
                                #print"Counts2: $counts_2\n";
                                $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                                #print"inter_pi $inter_pi\n";
                            }
                        }
                        $sample1_pi = $sample1_pi + $sample_locus_pi{$sample1}{$locus};
                        $sample2_pi = $sample2_pi + $sample_locus_pi{$sample2}{$locus};
                    }
                    $inter_pi = $inter_pi/$gene_length{$gene};
                    $sample_sample_gene_pi{$sample1}{$sample2}{$gene} = $inter_pi;
                    $sample_sample_gene_pi{$sample2}{$sample1}{$gene} = $inter_pi;
                    $sample1_pi = $sample1_pi/$gene_length{$gene};
                    $sample2_pi = $sample2_pi/$gene_length{$gene};
                    if ($inter_pi > 0) {
                        #$fst = 1 - 0.5*($sample_gene_pi{$sample1}{$gene} + $sample_gene_pi{$sample2}{$gene})/$inter_pi; # intra pi based on all loci
                        $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi; # based on only shared loci
                        $fst = sprintf("%.4f", $fst);
                    } else {
                        $fst = "NA"; # i.e. no intra-pi in any of the two samples and consequently no inter-pi.
                    }
                    $sample_sample_gene_fst{$sample1}{$sample2}{$gene} = $fst;
                    $sample_sample_gene_fst{$sample2}{$sample1}{$gene} = $fst;
                    #print"$gene $sample1 $sample2 $sample_gene_pi{$sample1}{$gene} $sample_gene_pi{$sample2}{$gene} $inter_pi $fst\n";
                }
            }
        }
    }
}

sub calc_per_gene_fst_permuted {
    local $perm;
    for ($ix1 = 0; $ix1 < @samples; $ix1++) {
        $sample1 = $samples[$ix1];
        #print"Sample1 $sample1\n";
        for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
            $sample2 = $samples[$ix2];
            #print"Sample2 $sample2\n";
            local @locus_sample1_pi = ();
            local @locus_sample2_pi = ();
            local @locus_inter_pi = ();
            local %gene_npermloci = ();
            #$apa_test = -1; # for testing
            #%gene_apa_test = (); # for testing
            foreach $contig (@contigs) {
                #print "Contig $contig\n";
                @genes = @{ $contig_genes{$contig} };
                foreach $gene (@genes) {
                    #print "Gene $gene\n";
                    @loci = (keys %{$gene_locus{$gene}});
                    next if (@loci == 0);
                    foreach $locus (@loci) { # add locus inter and intra pi values to array
                        $inter_pi = 0;
                        #print"Locus $locus\n";
                        next if (!defined $sample_locus_allel_counts{$sample1}{$locus});
                        next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                        @alleles1 = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                        @alleles2 = (keys %{$sample_locus_allel_counts{$sample2}{$locus}});
                        $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                        $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                        #print"Alleles1: @alleles1 Alleles2: @alleles2\n";
                        for ($i = 0; $i < @alleles1; $i++) {
                            #print"Allele 1: $alleles1[$i]\n";
                            $counts_1 = $sample_locus_allel_counts{$sample1}{$locus}{$alleles1[$i]};
                            #print"Counts1: $counts_1\n";
                            for ($j = 0; $j < @alleles2; $j++) {
                                next if ($alleles1[$i] eq $alleles2[$j]);
                                #print"Allele 2: $alleles2[$j]\n";
                                $counts_2 = $sample_locus_allel_counts{$sample2}{$locus}{$alleles2[$j]};
                                #print"Counts2: $counts_2\n";
                                $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                                #print"inter_pi $inter_pi\n";
                            }
                        }
                        if ($inter_pi > 0) { # note, only loci with inter_pi for the specific sample pair will be included
                            #$apa_test++; # for testing
                            push(@locus_sample1_pi, $sample_locus_pi{$sample1}{$locus});
                            push(@locus_sample2_pi, $sample_locus_pi{$sample2}{$locus});
                            push(@locus_inter_pi, $inter_pi);
                            if (defined $gene_npermloci{$gene}) {
                                $gene_npermloci{$gene}++;
                            } else {
                                $gene_npermloci{$gene} = 1;
                            }
                            #$gene_apa_test{$gene}{$apa_test} = 1; # for testing
                        }
                    }
                }
            }
            foreach $contig (@contigs) {
                #print "Contig $contig\n";
                @genes = @{ $contig_genes{$contig} };
                foreach $gene (@genes) {
                    if (defined $gene_npermloci{$gene}) {
                        for ($perm = 0; $perm < $n_fst_permutations; $perm++) {
                            $inter_pi = 0;
                            $sample1_pi = $sample2_pi = 0;
                            for ($i = 0; $i < $gene_npermloci{$gene}; $i++) {
                                $ix = int(rand(@locus_sample1_pi));
                                $sample1_pi = $sample1_pi + $locus_sample1_pi[$ix];
                                #splice(@locus_sample1_pi, $ix, 1); # adding this gives sampling without replacement
                                $sample2_pi = $sample2_pi + $locus_sample2_pi[$ix];
                                #splice(@locus_sample2_pi, $ix, 1); # adding this gives sampling without replacement
                                $inter_pi = $inter_pi + $locus_inter_pi[$ix];
                                #splice(@locus_inter_pi, $ix, 1); # adding this gives sampling without replacement
                            }
                            $inter_pi = $inter_pi/$gene_length{$gene};
                            $sample1_pi = $sample1_pi/$gene_length{$gene};
                            $sample2_pi = $sample2_pi/$gene_length{$gene};
                            if ($inter_pi > 0) {
                                $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi;
                                $fst = sprintf("%.4f", $fst);
                            } else {
                                $fst = "NA"; # i.e. no intra-pi in any of the two samples and consequently no inter-pi.
                            }
                            $sample_sample_gene_fst_perm{$sample1}{$sample2}{$gene}{$perm} = $fst;
                            $sample_sample_gene_fst_perm{$sample2}{$sample1}{$gene}{$perm} = $fst;
                            #print"$gene $sample1 $sample2 $perm $sample1_pi $sample2_pi $inter_pi $fst\n";
                        }
                        
                        ## temporary testing code
                        #print"temporary testing code gene: $gene\n";
                        #$perm = 0;
                        #$inter_pi = 0;
                        #$sample1_pi = $sample2_pi = 0;
                        #foreach $ix (keys %{$gene_apa_test{$gene}}) {
                        #    $sample1_pi = $sample1_pi + $locus_sample1_pi[$ix];
                        #    #splice(@locus_sample1_pi, $ix, 1); # adding this gives sampling without replacement
                        #    $sample2_pi = $sample2_pi + $locus_sample2_pi[$ix];
                        #    #splice(@locus_sample2_pi, $ix, 1); # adding this gives sampling without replacement
                        #    $inter_pi = $inter_pi + $locus_inter_pi[$ix];
                        #    #splice(@locus_inter_pi, $ix, 1); # adding this gives sampling without replacement
                        #}
                        #$inter_pi = $inter_pi/$gene_length{$gene};
                        #$sample1_pi = $sample1_pi/$gene_length{$gene};
                        #$sample2_pi = $sample2_pi/$gene_length{$gene};
                        #if ($inter_pi > 0) {
                        #    $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi;
                        #    $fst = sprintf("%.4f", $fst);
                        #} else {
                        #    $fst = "NA"; # i.e. no intra-pi in any of the two samples and consequently no inter-pi.
                        #}
                        #$sample_sample_gene_fst_perm{$sample1}{$sample2}{$gene}{$perm} = $fst;
                        #$sample_sample_gene_fst_perm{$sample2}{$sample1}{$gene}{$perm} = $fst;
                        ##print"$gene $sample1 $sample2 $perm $sample1_pi $sample2_pi $inter_pi $fst\n";
                        ## end temporary testing code

                    }
                }
            }
        }
    }
}

sub calc_per_gene_aminoacid_fst {
    foreach $contig (@contigs) {
        #print "Contig $contig\n";
        @genes = @{ $contig_genes{$contig} };
        foreach $gene (@genes) {
            #print "Gene $gene\n";
            @loci = (keys %{$gene_locus{$gene}});
            next if (@loci == 0);
            for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                $sample1 = $samples[$ix1];
                #print"Sample1 $sample1\n";
                for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
                    $sample2 = $samples[$ix2];
                    #print"Sample2 $sample2\n";
                    $inter_pi = 0;
                    $sample1_pi = $sample2_pi = 0;
                    foreach $locus (@loci) {
                        #print"\nLocus $locus\n";
                        next if (!defined $sample_locus_allel_counts{$sample1}{$locus});
                        next if (!defined $sample_locus_allel_counts{$sample2}{$locus});
                        @peptides1 = (keys %{$sample_locus_gene_peptide_counts{$sample1}{$locus}{$gene}});
                        @peptides2 = (keys %{$sample_locus_gene_peptide_counts{$sample2}{$locus}{$gene}});
                        $tot_count1 = $sample_locus_totcount{$sample1}{$locus};
                        $tot_count2 = $sample_locus_totcount{$sample2}{$locus};
                        #print"Peptides1: @peptides1 Peptides2: @peptides2\n";
                        for ($i = 0; $i < @peptides1; $i++) {
                            #print"Peptide 1: $peptides1[$i]\n";
                            $counts_1 = $sample_locus_gene_peptide_counts{$sample1}{$locus}{$gene}{$peptides1[$i]};
                            #print"Counts1: $counts_1\n";
                            for ($j = 0; $j < @peptides2; $j++) {
                                next if ($peptides1[$i] eq $peptides2[$j]);
                                #print"Peptide 2: $peptides2[$j]\n";
                                $counts_2 = $sample_locus_gene_peptide_counts{$sample2}{$locus}{$gene}{$peptides2[$j]};
                                #print"Counts2: $counts_2\n";
                                $inter_pi = $inter_pi + ($counts_1/$tot_count1)*($counts_2/$tot_count2);
                            }
                        }
                        $sample1_pi = $sample1_pi + $sample_locus_aminoacid_pi{$sample1}{$locus};
                        $sample2_pi = $sample2_pi + $sample_locus_aminoacid_pi{$sample2}{$locus};
                    }
                    $inter_pi = $inter_pi/$gene_length{$gene};
                    #print"$gene $sample1 $sample2 $inter_pi\n";
                    $sample_sample_gene_aminoacid_pi{$sample1}{$sample2}{$gene} = $inter_pi;
                    $sample_sample_gene_aminoacid_pi{$sample2}{$sample1}{$gene} = $inter_pi;
                    $sample1_pi = $sample1_pi/$gene_length{$gene};
                    $sample2_pi = $sample2_pi/$gene_length{$gene};
                    if ($inter_pi > 0) {
                        #$fst = 1 - 0.5*($sample_gene_aminoacid_pi{$sample1}{$gene} + $sample_gene_aminoacid_pi{$sample2}{$gene})/$inter_pi;
                        $fst = 1 - 0.5*($sample1_pi + $sample2_pi)/$inter_pi; # based on only shared loci
                        $fst = sprintf("%.4f", $fst);
                    } else {
                        $fst = "NA"; # i.e. no intra-pi in any of the two samples and consequently no inter-pi.
                    }
                    $sample_sample_gene_aminoacid_fst{$sample1}{$sample2}{$gene} = $fst;
                    $sample_sample_gene_aminoacid_fst{$sample2}{$sample1}{$gene} = $fst;
                    #print"$gene $sample1 $sample2 $sample1_pi $sample2_pi $inter_pi $fst\n"; # based on only shared loci
                }
            }
        }
    }
}

sub calc_pN_pS {
    foreach $sample (@samples_plus) {
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $PN = 0; # non-synonomous polymorphisms
                $PS = 0; # synonomous polymorphisms
                $TN = 0; # possible non-synonomous polymorphisms
                $TS = 0; # possible synonomous polymorphisms
                # $pNpS = ($PN/$TN) / ($PS/$TS) # NA if $PS = 0;
                # pN equals the fraction of polymorphic nonsynonymous sites, pS equals the fraction of polymorphic synonymous sites
                $contig = $gene_contig{$gene};
                $gene_seq = $gene_seq{$gene};
                # Gene on "+" strand
                if ($gene_strand{$gene} eq "+") {
                    #print"\n$contig $gene $sample\n$gene_seq\n";
                    for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                        $codon = substr($gene_seq, $i, 3);
                        # modify codon to represent majority sequence
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_start{$gene} + $i + $j;
                            $locus = $contig."|".$pos;
                            $base = substr($gene_seq, $i + $j, 1);
                            if (defined $sample_locus_totcount{$sample}{$locus}) {
                                #print"$locus $gene_start{$gene} $i  $j\n";
                                #print"codon before: $codon aa before: $codon_aminoacid{$codon}\n";
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                #print"Alleles: @alleles\n";
                                #print"Base: $base\n";
                                $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$base};
                                for ($k = 0; $k < @alleles; $k++) {
                                    if (length($alleles[$k]) != 1) {
                                        print "\nError: Variant > 1 bp. Violates pN/pS calculation. Consider not running in --keep_haplotypes mode.\n\n"; die;
                                    }
                                    if ($sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]} > $base_counts) {
                                        $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                        $base = $alleles[$k];
                                    }
                                }
                                substr($codon, $j, 1) = $base;
                            }
                        }
                        next if ($codon_aminoacid{$codon} eq "*"); # CHECK!!!
                        # update $PN $TN $PS $TS
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_start{$gene} + $i + $j;
                            $locus = $contig."|".$pos;
                            $base = substr($codon, $j, 1);
                            foreach $nt (keys %nt) {
                                next if ($base eq $nt);
                                $mod_codon = $codon;
                                substr($mod_codon, $j, 1) = $nt;
                                next if ($codon_aminoacid{$mod_codon} eq "*"); # CHECK!!!
                                if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                    $TS++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$nt}) {
                                            #print"PS $j $nt codon $codon aa $codon_aminoacid{$codon} mod_codon $mod_codon mod_aa $codon_aminoacid{$mod_codon}\n";
                                            $PS++ if ($sample_locus_allel_counts{$sample}{$locus}{$nt} > 0);
                                        }
                                    }
                                } else {
                                    $TN++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$nt}) {
                                            #print"PN $j $nt codon $codon aa $codon_aminoacid{$codon} mod_codon $mod_codon mod_aa $codon_aminoacid{$mod_codon}\n";
                                            $PN++ if ($sample_locus_allel_counts{$sample}{$locus}{$nt} > 0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                
                # Gene on "-" strand
                if ($gene_strand{$gene} eq "-") {
                    $gene_seq = &make_revcomp($gene_seq);
                    #print"$gene_strand{$gene}\n$gene_seq\n\n";
                    #print"\n$contig $gene $sample\n";
                    for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                        $codon = substr($gene_seq, $i, 3);
                        # modify codon to represent majority sequence
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_end{$gene} - ($i + $j); # ok?
                            $locus = $contig."|".$pos;
                            $base = substr($gene_seq, $i + $j, 1);
                            if (defined $sample_locus_totcount{$sample}{$locus}) {
                                #print"$codon\n";
                                $rc_base = &make_revcomp($base);
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                #print"Alleles: @alleles\n";
                                #print"RC_base: $rc_base\n\n";
                                if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_base}) {
                                    $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$rc_base};
                                } else {
                                    $base_counts = 0;
                                }
                                for ($k = 0; $k < @alleles; $k++) {
                                    if (length($alleles[$k]) != 1) {
                                        print "\nError: Variant > 1 bp. Violates pN/pS calculation.\n\n"; die;
                                    }
                                    if ($sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]} > $base_counts) {
                                        $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                        $rc_base = $alleles[$k];
                                    }
                                }
                                $base = &make_revcomp($rc_base);
                                substr($codon, $j, 1) = $base;
                            }
                        }
                        #print"$codon\n";
                        next if ($codon_aminoacid{$codon} eq "*"); # CHECK!!!
                        # upgrade $PN $TN $PS $TS
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_end{$gene} - ($i + $j); # ok?
                            $locus = $contig."|".$pos;
                            $base = substr($codon, $j, 1);
                            foreach $nt (keys (%nt)) {
                                next if ($base eq $nt);
                                $mod_codon = $codon;
                                substr($mod_codon, $j, 1) = $nt;
                                next if ($codon_aminoacid{$mod_codon} eq "*"); # CHECK!!!
                                $rc_nt = &make_revcomp($nt);
                                if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                    $TS++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_nt}) {
                                            $PS++ if ($sample_locus_allel_counts{$sample}{$locus}{$rc_nt} > 0);
                                        }
                                    }
                                } else {
                                    $TN++;
                                    if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_nt}) {
                                            $PN++ if ($sample_locus_allel_counts{$sample}{$locus}{$rc_nt} > 0);
                                        }
                                    }
                                }
                                
                            }
                        }
                    }
                }

                if ($use_pseudocounts) {
                    if ($TS > 0) {
                        $pNpS = (($PN + 1)/$TN) / (($PS + 1)/$TS);
                    } else {
                        $pNpS = "NA";
                    }
                } elsif ($PS > 0) {
                    $pNpS = ($PN/$TN) / ($PS/$TS);
                } else {
                    $pNpS = "NA";
                }
                $sample_gene_pNpS{$sample}{$gene} = $pNpS;
                #print"$gene_strand{$gene} PN: $PN PS: $PS TN: $TN TS: $TS pNpS: $pNpS\n";
            }
        }
    }
}

sub calc_NI { # calculates Neutrality Index (NI)
    # 1. spara positioner med PN och PS för varje prov
    # 2. för varje par av prov, använd bara positioner där båda prov har data, räkna ut DN och DS för dessa, räkna ut PN och PS för varje prov baserat på 1.
    foreach $contig (@contigs) {
        #$contig = "P1994_127_bin61_k141_448048";
        @genes = @{ $contig_genes{$contig} };
        foreach $gene (@genes) {
            #print"$gene\n";
            $contig = $gene_contig{$gene};
            $gene_seq = $gene_seq{$gene};
            %sample_codonix_codon = ();
            # record synonomous and non-synonomous polymorhisms iin every sample
            foreach $sample (@samples) {
                # Gene on "+" strand
                if ($gene_strand{$gene} eq "+") {
                    for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                        $codon = substr($gene_seq, $i, 3);
                        # modify codon to represent majority sequence
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_start{$gene} + $i + $j;
                            $locus = $contig."|".$pos;
                            $base = substr($gene_seq, $i + $j, 1);
                            if (defined $sample_locus_totcount{$sample}{$locus}) {
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                if (defined $sample_locus_allel_counts{$sample}{$locus}{$base}) {
                                    $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$base};
                                } else {
                                    $base_counts = 0;
                                }
                                for ($k = 0; $k < @alleles; $k++) {
                                    if (length($alleles[$k]) != 1) {
                                        print "\nError: Variant > 1 bp. Violates pN/pS calculation. Consider not running in --keep_haplotypes mode.\n\n"; die;
                                    }
                                    if ($sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]} > $base_counts) {
                                        $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                        $base = $alleles[$k];
                                    }
                                }
                                substr($codon, $j, 1) = $base;
                            }
                        }
                        $sample_codonix_codon{$sample}{$i} = $codon;
                        next if ($codon_aminoacid{$codon} eq "*"); # CHECK!!!
                        # record synonomous and non-synonomous polymorhisms, only considering bi-allelic loci
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_start{$gene} + $i + $j;
                            $locus = $contig."|".$pos;
                            if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                if (@alleles < 3) {
                                    $sample_locus_syn{$sample}{$locus} = 0;
                                    $sample_locus_nonsyn{$sample}{$locus} = 0;
                                    $base = substr($codon, $j, 1);
                                    foreach $nt (@alleles) {
                                        next if ($base eq $nt);
                                        next if ($sample_locus_allel_counts{$sample}{$locus}{$nt} == 0);
                                        $mod_codon = $codon;
                                        substr($mod_codon, $j, 1) = $nt;
                                        next if ($codon_aminoacid{$mod_codon} eq "*"); # CHECK!!!
                                        if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                            $sample_locus_syn{$sample}{$locus} = 1;
                                        } else {
                                            $sample_locus_nonsyn{$sample}{$locus} = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                # Gene on "-" strand
                if ($gene_strand{$gene} eq "-") {
                    for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                        $codon = substr($gene_seq, $i, 3);
                        # modify codon to represent majority sequence
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_end{$gene} - ($i + $j); # ok?
                            $locus = $contig."|".$pos;
                            $base = substr($gene_seq, $i + $j, 1);
                            if (defined $sample_locus_totcount{$sample}{$locus}) {
                                $rc_base = &make_revcomp($base);
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                if (defined $sample_locus_allel_counts{$sample}{$locus}{$rc_base}) {
                                    $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$rc_base};
                                } else {
                                    $base_counts = 0;
                                }
                                for ($k = 0; $k < @alleles; $k++) {
                                    if (length($alleles[$k]) != 1) {
                                        print "\nError: Variant > 1 bp. Violates pN/pS calculation.\n\n"; die;
                                    }
                                    if ($sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]} > $base_counts) {
                                        $base_counts = $sample_locus_allel_counts{$sample}{$locus}{$alleles[$k]};
                                        $rc_base = $alleles[$k];
                                    }
                                }
                                $base = &make_revcomp($rc_base);
                                substr($codon, $j, 1) = $base;
                            }
                        }
                        $sample_codonix_codon{$sample}{$i} = $codon;
                        next if ($codon_aminoacid{$codon} eq "*"); # CHECK!!!
                        # record synonomous and non-synonomous polymorhisms, only considering bi-allelic loci
                        for ($j = 0; $j < 3; $j++) {
                            $pos = $gene_end{$gene} - ($i + $j); # ok?
                            $locus = $contig."|".$pos;
                            if (defined $sample_locus_allel_counts{$sample}{$locus}) {
                                @alleles = (keys %{$sample_locus_allel_counts{$sample}{$locus}});
                                if (@alleles < 3) {
                                    $sample_locus_syn{$sample}{$locus} = 0;
                                    $sample_locus_nonsyn{$sample}{$locus} = 0;
                                    $base = substr($codon, $j, 1);
                                    foreach $nt (@alleles) {
                                        $rc_nt = &make_revcomp($nt);
                                        next if ($base eq $rc_nt);
                                        next if ($sample_locus_allel_counts{$sample}{$locus}{$nt} == 0);
                                        $mod_codon = $codon;
                                        substr($mod_codon, $j, 1) = $rc_nt;
                                        next if ($codon_aminoacid{$mod_codon} eq "*"); # CHECK!!!
                                        #die if ($locus eq "P1994_127_bin61_k141_448048|1237");
                                        if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                            $sample_locus_syn{$sample}{$locus} = 1;
                                        } else {
                                            $sample_locus_nonsyn{$sample}{$locus} = 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            # record synonomous and non-synonomous "divergence" in every pair of samples, and calculate NI (neutrality index)
            for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                $sample1 = $samples[$ix1];
                for ($ix2 = $ix1 + 1; $ix2 < @samples; $ix2++) {
                    $sample2 = $samples[$ix2];
                    $ps1 = $ps2 = $pn1 = $pn2 = $ds = $dn = 0;
                    if ($gene_strand{$gene} eq "+") {
                        for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                            for ($j = 0; $j < 3; $j++) {
                                next if ($codon_aminoacid{$sample_codonix_codon{$sample1}{$i}} eq "*");
                                next if ($codon_aminoacid{$sample_codonix_codon{$sample2}{$i}} eq "*");
                                $pos = $gene_start{$gene} + $i + $j;
                                $locus = $contig."|".$pos;
                                if (defined $sample_locus_totcount{$sample1}{$locus}) {
                                    if (defined $sample_locus_totcount{$sample2}{$locus}) {
                                        @alleles = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                                        if (@alleles == 2) {
                                            if (($sample_locus_allel_counts{$sample1}{$locus}{$alleles[0]} == 0) and ($sample_locus_allel_counts{$sample2}{$locus}{$alleles[1]} == 0)) {
                                                $codon = $sample_codonix_codon{$sample1}{$i};
                                                $mod_codon = $codon;
                                                if (substr($codon, $j, 1) eq $alleles[0]) {
                                                    substr($mod_codon, $j, 1) = $alleles[1];
                                                } else {
                                                    substr($mod_codon, $j, 1) = $alleles[0];
                                                }
                                            
                                                if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                                    $ds++;
                                                    print"+ $gene ds\n";
                                                } else {
                                                    $dn++;
                                                    print"+ $gene dn\n";
                                                }
                                            } else {
                                                $ps1 = $ps1 + $sample_locus_syn{$sample1}{$locus};
                                                $ps2 = $ps2 + $sample_locus_syn{$sample2}{$locus};
                                                $pn1 = $pn1 + $sample_locus_nonsyn{$sample1}{$locus};
                                                $pn2 = $pn2 + $sample_locus_nonsyn{$sample2}{$locus};
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if ($gene_strand{$gene} eq "-") {
                        $gene_seq = &make_revcomp($gene_seq);
                        for ($i = 3; $i < $gene_length{$gene}; $i = $i + 3) { # skipping the start codon
                            for ($j = 0; $j < 3; $j++) {
                                next if ($codon_aminoacid{$sample_codonix_codon{$sample1}{$i}} eq "*");
                                next if ($codon_aminoacid{$sample_codonix_codon{$sample2}{$i}} eq "*");
                                $pos = $gene_end{$gene} - ($i + $j); # ok?
                                $locus = $contig."|".$pos;
                                if (defined $sample_locus_totcount{$sample1}{$locus}) {
                                    if (defined $sample_locus_totcount{$sample2}{$locus}) {
                                        @alleles = (keys %{$sample_locus_allel_counts{$sample1}{$locus}});
                                        if (@alleles == 2) {
                                            if (($sample_locus_allel_counts{$sample1}{$locus}{$alleles[0]} == 0) and ($sample_locus_allel_counts{$sample2}{$locus}{$alleles[1]} == 0)) {
                                                $codon = $sample_codonix_codon{$sample1}{$i};
                                                $mod_codon = $codon;
                                                $rc_nt = &make_revcomp($alleles[0]);
                                                if (substr($codon, $j, 1) eq $rc_nt) {
                                                    substr($mod_codon, $j, 1) = &make_revcomp($alleles[0]);
                                                } else {
                                                    substr($mod_codon, $j, 1) = $rc_nt;
                                                }
                                                if ($codon_aminoacid{$mod_codon} eq $codon_aminoacid{$codon}) {
                                                    $ds++;
                                                } else {
                                                    $dn++;
                                                }
                                            } else {
                                                $ps1 = $ps1 + $sample_locus_syn{$sample1}{$locus};
                                                $ps2 = $ps2 + $sample_locus_syn{$sample2}{$locus};
                                                $pn1 = $pn1 + $sample_locus_nonsyn{$sample1}{$locus};
                                                $pn2 = $pn2 + $sample_locus_nonsyn{$sample2}{$locus};
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    $sample_sample_gene_ni{$sample1}{$sample2}{$gene} = "NA";
                    $sample_sample_gene_ni{$sample2}{$sample1}{$gene} = "NA";
                    if ($ps1*$dn*$ds > 0) {
                        $sample_sample_gene_ni{$sample1}{$sample2}{$gene} = ($pn1/$ps1)/($dn/$ds);
                    }
                    if ($ps2*$dn*$ds > 0) {
                        $sample_sample_gene_ni{$sample2}{$sample1}{$gene} = ($pn2/$ps2)/($dn/$ds);
                    }
                    #print"$gene $sample1 $sample2 pn1:$pn1 ps1:$ps1 pn2:$pn2 ps2:$ps2 dn:$dn ds:$ds ni1:$sample_sample_gene_ni{$sample1}{$sample2}{$gene} ni2:$sample_sample_gene_ni{$sample2}{$sample1}{$gene}\n";
                }
            }
        }
    }
}


sub make_revcomp {
    local($seq) = $_[0];
    local($modseq) = "";
    local($i);
    local($nt);
    local($i);
    for ($i = 0; $i < length($seq); $i++) {
        $nt = substr($seq, $i, 1);
        if ($nt eq "A") {
            $modseq = T.$modseq;
        } elsif ($nt eq "T") {
            $modseq = A.$modseq;
        } elsif ($nt eq "C") {
            $modseq = G.$modseq;
        } elsif ($nt eq "G") {
            $modseq = C.$modseq;
        } else {
            $modseq = $nt.$modseq;
        }
    }
    return ($modseq);
}

sub print_output_to_file {
    print"$outprefix.logfile.txt\n";
    open (OUT, ">$outprefix.logfile.txt");
    print OUT "$logtext";
    close(OUT);
    ###
    local($sample_norm_pi);
    print"$outprefix.intradiv.txt\n";
    open (OUT, ">$outprefix.intradiv.txt");
    print OUT "Sample\tIntra_pi\tNorm_intra_pi\tEst_genome_cov\tNum_loci\tTotal_num_alleles\tAverage_read_depth\n";
    foreach $sample (@samples_plus) {
        if (defined $sample_estcov{$sample}) {
            $sample_norm_pi = $sample_pi{$sample}*$est_num_loci/$num_loci{$sample};
            print OUT "$sample\t$sample_pi{$sample}\t$sample_norm_pi\t$sample_estcov{$sample}\t$num_loci{$sample}\t$sample_totalleles{$sample}\t$sample_avcount{$sample}\n";
            #print "#$sample $sample_estcov{$sample} $sample_pi{$sample} $sample_norm_pi $num_loci{$sample} $est_num_loci\n";
        } else {
            print OUT "$sample\t$sample_pi{$sample}\tNA\tNA\t$num_loci{$sample}\t$sample_totalleles{$sample}\t$sample_avcount{$sample}\n";
        }
    }
    close(OUT);
    if ($pi_only) {
        print"\n### Finished pogenom succesfully (--pi_only mode) ###\n\n";
        exit();
    }
    ###
    print"$outprefix.intradiv-per-locus.txt\n";
    open (OUT, ">$outprefix.intradiv-per-locus.txt");
    print OUT "Contig\tPosition";
    foreach $sample (@samples_plus) {
        print OUT "\t$sample Intra_locus_pi\t$sample Locus_Depth";
    }
    print OUT "\n";
    foreach $contig (keys %contig_pos) {
        @positions = (keys %{$contig_pos{$contig}});
        @positions = sort {$a <=> $b} @positions;
        foreach $pos (@positions) {
            print OUT "$contig\t$pos";
            $locus = $contig."|".$pos;
            foreach $sample (@samples_plus) {
                if (defined $sample_locus_pi{$sample}{$locus}) {
                    print OUT "\t$sample_locus_pi{$sample}{$locus}\t$sample_locus_totcount{$sample}{$locus}";
                } else {
                    print OUT "\tNA\tNA";
                }
            }
            print OUT "\n";
        }
    }
    close(OUT);
    ###
    print"$outprefix.allele-freqs.txt\n";
    @nts = ('A', 'T', 'C', 'G');
    open (OUT, ">$outprefix.allele-freqs.txt");
    print OUT "Contig\tPosition";
    foreach $sample (@samples_plus) {
        print OUT "\t$sample A\t$sample T\t$sample C\t$sample G";
    }
    print OUT "\n";
    foreach $contig (keys %contig_pos) {
        @positions = (keys %{$contig_pos{$contig}});
        @positions = sort {$a <=> $b} @positions;
        foreach $pos (@positions) {
            next if ($contig_pos{$contig}{$pos} == 0);
            #print "$contig\t$pos\t$contig_pos{$contig}{$pos}\n";
            print OUT "$contig\t$pos";
            $locus = $contig."|".$pos;
            foreach $sample (@samples_plus) {
                if (defined $sample_locus_pi{$sample}{$locus}) {
                    foreach $nt (@nts) {
                        $counts = 0;
                        if (defined $sample_locus_allel_counts{$sample}{$locus}{$nt}) {
                            $counts = $sample_locus_allel_counts{$sample}{$locus}{$nt};
                        }
                        print OUT "\t$counts";
                    }
                } else {
                    print OUT "\tNA\tNA\tNA\tNA";
                }
            }
            print OUT "\n";
        }
    }
    close(OUT);
    ###
    if ($gff_file) {
        print"$outprefix.aminoacid-freqs.txt\n";
        @aas = ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*');
        open (OUT, ">$outprefix.aminoacid-freqs.txt");
        print OUT "Contig\tGene\tAminoAcidPosition";
        foreach $sample (@samples_plus) {
            foreach $aminoacid (@aas) {
                print OUT "\t$sample $aminoacid";
            }
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                @codon_positions = (keys %{$sample_gene_codon_aminoacid_counts{'All_samples_combined'}{$gene}});
                @codon_positions = sort {$a <=> $b} @codon_positions;
                #print "@codon_positions\n";
                foreach $codon_pos (@codon_positions) {
                    $output_codon_pos = $codon_pos + 1;
                    print OUT "$contig\t$gene\t$output_codon_pos"; # 1 is first codon (start codon)
                    foreach $sample (@samples_plus) {
                        if (defined $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}) {
                            foreach $aminoacid (@aas) {
                                if (defined $sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid}) {
                                    print OUT "\t$sample_gene_codon_aminoacid_counts{$sample}{$gene}{$codon_pos}{$aminoacid}";
                                } else {
                                    print OUT "\t0";
                                }
                            }
                        } else {
                            foreach $aminoacid (@aas) {
                                print OUT "\tNA";
                            }
                        }
                    }
                    print OUT "\n";
                }
            }
        }
    }
    ###
    if ($gff_file) {
        print"$outprefix.intradiv-per-gene.txt\n";
        open (OUT, ">$outprefix.intradiv-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        foreach $sample (@samples_plus) {
            #print OUT "\t$sample Intra_gene_pi\t$sample Intra_gene_pi_p-value";
            print OUT "\t$sample Intra_gene_pi";
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                foreach $sample (@samples_plus) {
                    #print OUT "\t$sample_gene_pi{$sample}{$gene}\tNA";
                    print OUT "\t$sample_gene_pi{$sample}{$gene}";
                }
                print OUT "\n";
            }
        }
        close(OUT);
    }
    ###
    if ($gff_file and $genetic_code_file) {
        print"$outprefix.intradiv-aminoacid-per-gene.txt\n";
        open (OUT, ">$outprefix.intradiv-aminoacid-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        foreach $sample (@samples_plus) {
            #print OUT "\t$sample Intra_gene_aminoacid_pi\t$sample Intra_gene_aminoacid_pi_p-value";
            print OUT "\t$sample Intra_gene_aminoacid_pi";
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                foreach $sample (@samples_plus) {
                    #print OUT "\t$sample_gene_aminoacid_pi{$sample}{$gene}\tNA";
                    print OUT "\t$sample_gene_aminoacid_pi{$sample}{$gene}";
                }
                print OUT "\n";
            }
        }
        close(OUT);
        
        print"$outprefix.pNpS-per-gene.txt\n";
        open (OUT, ">$outprefix.pNpS-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        foreach $sample (@samples_plus) {
            #print OUT "\t$sample Intra_gene_aminoacid_pi\t$sample Intra_gene_aminoacid_pi_p-value";
            print OUT "\t$sample pNpS";
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            #print "$contig\n";
            @genes = @{ $contig_genes{$contig} };
            #print "@genes\n"; die;
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                foreach $sample (@samples_plus) {
                    if (defined $sample_gene_pNpS{$sample}{$gene}) {
                        print OUT "\t$sample_gene_pNpS{$sample}{$gene}";
                    } else {
                        print OUT "\tNA";
                    }
                }
                print OUT "\n";
            }
        }
        close(OUT);

    }
    ###
    if (@samples > 1) {
        print"$outprefix.fst.txt\n";
        open (OUT, ">$outprefix.fst.txt");
        foreach $sample (@samples) {
            print OUT "\t$sample";
        }
        print OUT "\n";
        for ($ix1 = 0; $ix1 < @samples; $ix1++) {
            print OUT "$samples[$ix1]";
            for ($ix2 = 0; $ix2 < @samples; $ix2++) {
                if ($ix1 == $ix2) {
                    print OUT "\tNA";
                } else {
                    print OUT "\t$sample_sample_fst{$samples[$ix1]}{$samples[$ix2]}";
                }
            }
            print OUT "\n";
        }
        close(OUT);
    }
    ###
    if ((@samples > 1) and $gff_file) {
        print"$outprefix.fst-per-gene.txt\n";
        open (OUT, ">$outprefix.fst-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        for ($ix1 = 0; $ix1 < @samples; $ix1++) {
            for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                print OUT "\t$samples[$ix1] $samples[$ix2]";
            }
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                    for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                        #print"$ix1 $ix2\n";
                        if (defined $sample_sample_gene_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}) {
                            print OUT "\t$sample_sample_gene_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}";
                        } else {
                            print OUT "\tNA";
                        }
                    }
                }
                print OUT "\n";
            }
        }
    }
    ###
    if ((@samples > 1) and $gff_file and $n_fst_permutations) {
        print"$outprefix.$n_fst_permutations.permuted-fst-per-gene.txt\n";
        open (OUT, ">$outprefix.$n_fst_permutations.permuted-fst-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        for ($i = 0; $i < $n_fst_permutations; $i++) {
            $perm = $i + 1;
            for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                    print OUT "\t$samples[$ix1] $samples[$ix2] #$perm";
                }
            }
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                for ($i = 0; $i < $n_fst_permutations; $i++) {
                    for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                        for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                            #print"$ix1 $ix2\n";
                            if (defined $sample_sample_gene_fst_perm{$samples[$ix1]}{$samples[$ix2]}{$gene}{$i}) {
                                print OUT "\t$sample_sample_gene_fst_perm{$samples[$ix1]}{$samples[$ix2]}{$gene}{$i}";
                            } else {
                                print OUT "\tNA";
                            }
                        }
                    }
                }
                print OUT "\n";
            }
        }
    }
    ###
    if ((@samples > 1) and $gff_file and $genetic_code_file) {
        print"$outprefix.fst-aminoacid-per-gene.txt\n";
        open (OUT, ">$outprefix.fst-aminoacid-per-gene.txt");
        print OUT "Contig\tGene\tLength\tStart\tEnd\tStrand\tNum_loci";
        for ($ix1 = 0; $ix1 < @samples; $ix1++) {
            for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                print OUT "\t$samples[$ix1] $samples[$ix2]";
            }
        }
        print OUT "\n";
        foreach $contig (@contigs) {
            @genes = @{ $contig_genes{$contig} };
            foreach $gene (@genes) {
                $num_loci = (keys %{$gene_locus{$gene}});
                print OUT "$contig\t$gene\t$gene_length{$gene}\t$gene_start{$gene}\t$gene_end{$gene}\t$gene_strand{$gene}\t$num_loci";
                for ($ix1 = 0; $ix1 < @samples; $ix1++) {
                    for ($ix2 = ($ix1 + 1); $ix2 < @samples; $ix2++) {
                        #print"$ix1 $ix2\n";
                        if (defined $sample_sample_gene_aminoacid_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}) {
                            print OUT "\t$sample_sample_gene_aminoacid_fst{$samples[$ix1]}{$samples[$ix2]}{$gene}";
                        } else {
                            print OUT "\tNA";
                        }
                    }
                }
                print OUT "\n";
            }
        }
    }
}


