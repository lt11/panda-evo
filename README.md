# panda-evo

## Overview

The `panda-evo` workflow is a CDS-oriented derivative of the [`panda`](https://github.com/lt11/panda) pipeline. The general pangenome-annotation workflow was adapted from `panda`, and the original repository remains the main reference for the broader method and background documentation.

The main difference is scope. While `panda` focuses on building and inspecting a feature-level pangenome, `panda-evo` focuses on CDS annotations to build a collection of valid coding sequences from a pangenome. These sequence collections are intended for downstream evolutionary analyses, for example `dN/dS` estimation.

The core workflow lives in `mrk/panda-evo.Rmd`. It stages GFF annotations, converts CDS annotations to BED, runs `impg` against a PAF alignment, reduces overlapping annotations into nonredundant CDS sub-blocks, and then builds one FASTA per CDS block by combining annotated CDS sequences with CDS reconstructed from genome sequence when needed.

The same logic is also split across the single-chunk scripts in `mrk/`, especially `mrk/c1.sh`, `mrk/c2.R`, `mrk/c3.R`, and `mrk/c4.R`.

For the original feature-level pipeline, repository structure, and additional documentation, refer to [`panda`](https://github.com/lt11/panda).

## Installation

This repository is not yet distributed as a packaged software. Instead, the workflow is executed in a local environment where all required tools and R packages are pre-installed.

All the following step are documented below:

- clone the repository
```
git clone https://github.com/lt11/panda-evo.git
cd panda
```

- install dependencies (e.g. `impg`) and make sure they are accessible in your `$PATH`

- install R packages, e.g.:
```
Rscript -e 'install.packages(c(
  "data.table",
  "rmarkdown"
))'
```

- configure external inputs

- check `impg` version compatibility

## Features

- Extraction of CDS from pangenome data
- Validation of coding sequences
- Generation of datasets tailored for evolutionary analyses
- Integration with the panda pipeline ecosystem

## Repository Layout

Top-level directories:

- `ids/ids-ps.txt` input lists of strain-haplotype identifiers used to define the analysis set and output column order.
- `aln/` input alignments in PAF format.
- `anno/` staged GFF and BED annotations produced during the workflow.
- `mrk/` the R Markdown workflow and single-chunk scripts.
- `png/` CDS pangenome tables produced by the pipeline.
- `seqs/` CDS FASTA outputs derived from the pangenome, currently under `seqs/cds/`.

Key source files:

- `mrk/panda-evo.Rmd`: end-to-end notebook.
- `mrk/c1.sh`: stages and normalises input GFF files.
- `mrk/c2.R`: converts GFF annotations to CDS BED intervals.
- `mrk/c3.R`: computes the CDS pangenome from `impg` output.
- `mrk/c4.R`: builds per-block CDS FASTA files and labels them as `non-evolvable`, `core`, or `dispensable`.
- `mrk/run-profvis.R`: helper to profile an R script with `profvis`.

## Inputs

The workflow expects these main inputs:

1. A single PAF alignment in `aln/`. E.g. you can use a PAF file generated from [`PGGB`](https://github.com/pangenome/pggb) or directly from [`wfmash`](https://github.com/waveygang/wfmash).
2. A list of strain-haplotype identifiers in (e.g. `ids/ids-ps.txt`). These must be the IDs of the sequences you used to build the input PAF alignment.
3. A directory of source GFF annotations outside this repository, configured in `mrk/c1.sh` and in `mrk/panda-evo.Rmd`.
4. CDS FASTA files and genome FASTA files outside this repository, used by `mrk/c4.R` to recover valid coding sequences.
5. A list of essential non-evolvable genes, also configured in `mrk/c4.R`. One for the S. cerevisiae is provided (`aux/essentiality.txt`)

Naming conventions matter:

- Contigs in the PAF use PanSN-style names with `#`, for example `SGDref#0#chrI`.
- Annotation files are named with `-`, for example `SGDref-0-features.gff`.
- `ids/ids-ps.txt` defines both the haplotypes (e.g. for example `SGDref-0` to analyse and the order of the output columns.

## Workflow

The pipeline runs in four main stages.

1. Stage GFF files.  
   `mrk/c1.sh` copies the requested annotation files into `anno/gff/`, optionally merging mitochondrial annotations.
2. Convert GFF to CDS BED.  
   `mrk/c2.R` filters the annotations to `CDS`, normalises coordinates, trims redundant suffixes from feature names, and writes BED files to `anno/bed/`.
3. Build the CDS pangenome.  
   `mrk/c3.R` runs `impg` for each BED file against the PAF alignment, groups overlapping hits into homology blocks, reduces them into nonredundant CDS sub-blocks, and writes the CDS pangenome tables to `png/`.
4. Build valid coding-sequence collections.  
   `mrk/c4.R` reads the CDS pangenome table and, for each CDS sub-block, writes a FASTA file in `seqs/cds/`. It first retrieves annotated CDS sequences from the CDS database; when a sequence is missing or too short, it tries to reconstruct a valid CDS directly from the genome sequence. Candidate sequences are accepted only if they satisfy coding constraints such as start codon, stop codon, and reading-frame consistency.

This last stage is the main purpose of `panda-evo`: to derive a curated collection of valid coding sequences from pangenome-defined CDS blocks for downstream comparative and evolutionary analyses.

## Sequence Classes

Each CDS block is labelled in `mrk/c4.R` as one of:

- `non-evolvable`: block contains an essential non-evolvable gene.
- `core`: block is present in more than `95%` of genomes.
- `dispensable`: block is below the core threshold.

These labels are embedded in the output FASTA filenames under `seqs/cds/`, for example `1-YDR381C-A-core.fa`.

## Dependencies

This repository is not packaged as an R package. The analysis assumes a local environment with:

- `bash`
- `Rscript`
- `rmarkdown` to render the notebook
- `impg`
- a PAF alignment, typically generated upstream from a pangenome/graph alignment workflow

R packages used in the scripts include:

- `data.table`
- `this.path`
- `scriptName`
- `GenomicRanges`
- `stringr`
- `Biostrings`
- `BiocParallel`
- `gtools`
- `tictoc`
- `profvis`

As in `panda`, the `impg` command may need to be adjusted to match the installed version:

- `impg 0.2.0`: `impg -I -p file.paf -b file.bed`
- `impg 0.2.3`: `impg query -I -p file.paf -b file.bed`

Check `mrk/c3.R` or the matching chunk in `mrk/panda-evo.Rmd` before running the workflow.

## Running an Analysis

There is no executable or standalone CLI entrypoint; the repository is driven from the notebook or by running a single chunk, as described below.

Before rendering:

- place exactly one input `.paf` file in `aln/`
- update `ids/ids-ps.txt` for the haplotypes to analyse
- set the external annotation, CDS, genome, and essential-gene paths in the notebook/scripts
- decide whether mitochondrial annotations should be included with `with_mito`
- confirm the `impg` command variant for the installed version
- review notebook chunks marked with `eval = FALSE`

Then render:

```bash
Rscript -e 'rmarkdown::render("mrk/panda-evo.Rmd")'
```

If you prefer to run the stages individually, the scripts map to the notebook sections:

```bash
bash mrk/c1.sh
Rscript mrk/c2.R
Rscript mrk/c3.R
Rscript mrk/c4.R
```

Profiling helper:

```bash
Rscript mrk/run-profvis.R c3.R
```

Run that command from `mrk/` if you want `run-profvis.R` to resolve the target script as written.

## Main Outputs

Primary CDS pangenome outputs are written to `png/`:

- `png/pan-features.txt`: CDS pangenome table. The first columns are `Class_id` and `Features_id`; remaining columns are haplotypes such as `SGDref#0` and `CMF#1`, containing CDS interval strings or `MA` when absent.
- `png/generators.txt`: the `class:feature#strand` entries that generated each output row.
- `pan-features.RData`: R object with the same CDS pangenome table.

Sequence outputs are written to `seqs/cds/`:

- one FASTA file per CDS sub-block
- filenames encode a numeric row index, a representative gene identifier, and the block class (`non-evolvable`, `core`, or `dispensable`)
- each FASTA collects the valid coding sequences associated with that pangenome CDS block

These FASTA collections are the key downstream product of `panda-evo`. They are meant to be used for analyses on homologous coding sequences across the pangenome, including examples such as `dN/dS` estimation.

## Notes

- `panda-evo` was derived from `panda`; for the original feature-level workflow and additional documentation, refer to the upstream repository: <https://github.com/lt11/panda>.
- The CDS validation step in `mrk/c4.R` enforces basic coding constraints and can reconstruct CDS candidates from genome sequence when the annotation-derived CDS FASTA is missing.
- `mrk/panda-evo.Rmd` documents the known `impg 0.2.0` warning related to asymmetric alignments.
