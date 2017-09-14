# ComputeDeltaESRSeq

A simple Perl script to compute DeltaESRSeq to predict exon skipping induced by exonic variants

## Installation

To run the script you need Perl and the CPAn module REST::Client, and an active internet connection (the scripts uses togows.org REST service)

## Run

```bash
perl -T ComputeDeltaESRSeq.pl -l path/to/variant_list.txt -g genome_version
```

* -l txt file, list of variants

input format must be a file of format:

*	either tabulated:

Chr	Pos	Strand

*	or variants in HGVS genomic format plus the strand of interest, e.g.:

chr1:g.87466579C>T	+

*	Please note that for variants on strand -, nucleotides provided can be on strand - or +, in the latter case they will be converted on the minus strand, however WT and mutant nucleotides must be provided on the same strand.

i.e.:

chrX    32481660        C       A	-

And

chrX    32481660        G       T	-

are the same.

* -g genome version, hg19/hg38

## Output

Output is a tabulated file with variants (sorted), ESR WT, ESR mutant, Delta ESR

## WARNING:

Providing exonic positions is your responsibility!!!

Works only for substitutions

## Credits

This script uses togoWS web services, http://togows.org/

Toshiaki Katayama, Mitsuteru Nakao and Toshihisa Takagi: TogoWS: integrated SOAP and REST APIs for interoperable bioinformatics Web services. Nucleic Acids Research 2010, 38:W706-W711. doi:10.1093/nar/gkq386, PMID:20472643 (NAR Web Server Issue 2010)

The method is adapted from [di Giacomo et al., 2013](https://www.ncbi.nlm.nih.gov/pubmed/23983145) and the raw data (hexamers scores) from [Ke et al., 2011](https://www.ncbi.nlm.nih.gov/pubmed/21659425).