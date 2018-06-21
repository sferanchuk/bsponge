# bsponge

A pipeline for open-reference clustering and annotation of amplicon sequencing data

## Description

### Prerequisites:

1. Input sequence libraries in fasta format
2. Reference database in fasta format, with taxonomy included in the sequence description or as a separate file
3. (optional and undocumented) Aligned reference sequences in fasta format (to be able to join alignments with different references)

### Steps of the pipeline:

1. For each input sample, align reads to reference database, with a single best hit, to obtain the alignments in .sam format.

2. Run open-reference clustering using the provided software (test_openref.cpp). The resulted clusters may include reads from all input samples.

3. Convert the formats and/or add the taxonomy associated with the reference database sequences

## Installation

The software was tested on Ubuntu Linux, with *usearch* software to run the alignments.




