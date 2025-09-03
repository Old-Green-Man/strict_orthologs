This program uses reciprocal best blast hits to link genes across two
genomes. For alignments with multiple HSPs the bit score is estimated from
non-overlapping segments. An ortholog group must be self consistent in that all
members are reciprocal best hit to all other members. Such a strict criteria may
be useful for phylogenetic analyses where paralogous groupings are to be
avoided. The downside is that a single gene can split an ortholog group into two
groups if it has a single spurious reciprocal best hit. I've seen this happen
with misannotated genes where two closely related genomes were assembled
together. This could happen with legitimate sequences also, but should be
relatively rare, but highly dependent on the genomes under consideration and
their phylogenetic distribution. The goal here is to have reliable clusters for
doing phylogenetic analyses, even if some clusters may be split etc (higher
specificity; lower sensitivity).

The input blast outputs should be in tabular format as specified below. In addition, all HSPs between genome A and genome B should be together in the blast output. This is assumed for efficient processing. The output should include hits between genes in the same genome.

The gene "names" in columns 1 and 2 should be in the following format: genomeA_geneX. The underscore is the default separator, but can be changed with the -sep option.

You can test the program on the example.blastp and example.len files in the example/ directory. For example

```
mkdir tmp
cd tmp
strict_orthologs.pl example.len example.blastp > ortho.out 2> ortho.err
```

The main output will be in ortho.out. A file called rbh.dat will be written that can be used to rerun the program without the blast output using the -rbh_file option.

A line in the main output looks like this:
`>1 size=6 bits_median=109 pcnt_median=59.28 1_1417:73.7:128:75 2_867:73.7:128:75 3_288:73.7:128:75 4_47:54.3:103:81 4_538:56:103:73 5_1029:73.7:128:76`
The \>1 is the ortholog index.
size=6 means there are 6 genes in this group
The bit and percent medians are the all vs all medians within the group.
Each gene is annotated like: 1_1417:73.7:128:75, where
1_1417 means gene 1417 in genome 1
73.7 is the median percent identity between this gene and the others
128 is the same for the bit scores
75 is the sequence length

**The following help information is printed if the program is run without arguments.**

usage: ./strict_orthologs.pl gene-lengths-file [tabular-blast-output-file tabular-blast-output-file2 ...] [options]

This reads tabular blast output. It assumes that blast was run with the
following -outfmt.  -outfmt="6 qacc sacc bitscore nident qstart qend sstart
send"

When run the first time on blastoutput it will write a file rbh.dat.  This can be read in using -rbh_file=rbh.dat without needing the blast outputs and speeding things up.

gene-lengths-file should contain the lengths of all genes (format: len`<WHITE_SPACE>`id)

options
-------
-rbh_file preparsed reciprocal best hits including paralog links.  This bypasses
 running groupParalogs and should include all relevant reciprocal best hits
 needed for clustering.

-sep (default: _) separator between taxon name and gene in the query and sbjct ids.
 original blast search.  Each line should have length and id separated by
 space.

-useParalogs (default: true Include paralogs that are specific to a
 taxon into the same ortholog group.  These paralogs will be considered a single
 gene internally that share the various reciprocal hits.

-minParalogPcnt (default: 50) if n_ident / max(q,s) * 100 is less than this then not counted as paralogs

-singlets (default: true) print singlets as their own group

-minPcnt  (default: 0) minimum value of 100 * identities / max( qlen, slen ) to keep
  NOT CURRENTLY USED

-minPcnt2 (default: 0) minimum value of 100 * identities / qcov to keep
 This second (normal percent identity) is used linking CRISPR arrays where we
 don't want penalize so much for subsequence matches as long as a sufficient
 amount is aligned as determined by -minAlen
  NOT CURRENTLY USED

-minBits  (default: 40) minimum bit value to keep

-minAlen (default: 10) minimum per HSP alignment length.  This is to
 capture alignments of adjacent repeats and spacers.  The length should be at least
 greater than the maximum expected length of either repeat or spacer.

-diffPcnt (default: 20) if an HSP differs in percent identity from the
 first HSP for a given hit by more than this exclude it from calculating overall
 pcntId, bits, and coverage.

-help

-verbose (default: false) verbose logging

-debug   (default: false) Run extra debug checks
