## Description
In this experiment, fly samples with Mettl3 modification were sequenced together with control samples (3 replicates in each condition). The goal is to evaluate the change in A-to-I RNA editing after the Mettl3 modification overall and individually in a given set of genes.

The amount of editing in flies can be estimated by measuring the level of A-to-I RNA editing at a given set of known genomic positions. In each given genomic position, the editing level can be computed from sequencing data as a fraction of reads with an edited base (inosine I is read as guanine G by the sequencer, so the ratio of G count to the read depth).

- Note that individual nucleotide counts are computed from the reads mapped to the DNA genomic reference and so the direction of the transcript must be taken into account.
- Note that editing levels in individual genomic positions might vary highly, but still the consistent smaller shift of editing between the conditions might be statistically significant.

## Input files
- Tables with individual nucleotide counts (ACGT) for a given set of edited sites defined by genomic coordinates for each sample.
- Table of genes where editing should be estimated with genomic coordinates.
