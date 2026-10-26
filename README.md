# CoReID

**CoReID** is a long-read sequencing workflow for identifying extrachromosomal circular DNA (eccDNA) from tandem-repeat evidence produced by rolling-circle amplification. It integrates circular consensus sequence (CCS) construction with raw-read alignment patterns to reconstruct candidate circular structures, then clusters concordant candidates across reads.

The project also provides a universal multi-sample integration script that standardizes and aggregates eccDNA calls from CoReID, CReSIL, NanoCircle, Flec, and FLED.

## Overview

For each long read, CoReID:

1. detects tandem-repeat structure and constructs a CCS with `pyccs`;
2. maps the raw read to a reference genome with mappy/minimap2;
3. screens for overlapping primary alignments consistent with repeated circular units;
4. reconstructs cyclic fragment order from the raw read whenever possible;
5. uses alignment of `CCS × 2` as a fallback when raw-read reconstruction is unresolved;
6. normalizes cyclic topology, clusters candidates with similar genomic structures, merges subset structures, and reports representative boundaries.

The current implementation does not use genome annotation or a pre-existing eccDNA database during calling.

## Installation

From a cloned repository:

```bash
python -m pip install .
```

For development, use an editable installation:

```bash
python -m pip install -e .
```

## Repository layout

The source files remain flat in the repository root. During installation, `pyproject.toml` maps this root directory to the Python package name `coreid`; no `ECCTRY/` or `coreid/` source subdirectory is required.

```text
CoReID/
├── __init__.py
├── main.py
├── align.py
├── cluster.py
├── env.py
├── find_ccs.py
├── find_eccdna.py
├── logger.py
├── utils.py
├── merge_eccdna.py
├── pyproject.toml
└── README.md
```

## Requirements

- Python 3.9 or later
- `pyccs`
- `mappy`
- `pysam`
- `pyfaidx`
- `numpy`
- `pandas`

All Python dependencies are installed automatically by `pip install .`.

```bash
pip install pyccs mappy pysam pyfaidx numpy
```

The reference genome should be FASTA indexed. If necessary, create its index before running:

```bash
samtools faidx reference.fa
```

## Input

CoReID accepts FASTA or FASTQ long-read files, either uncompressed or gzip-compressed:

- `.fa`, `.fasta`, `.fa.gz`, `.fasta.gz`
- `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`

The workflow requires a reference genome FASTA file. Reads without a successfully constructed CCS are not analyzed by the current calling workflow.

## Commands

Installation provides two commands:

```text
coreid
coreid-merge
```

### Per-sample eccDNA calling

```bash
coreid call_ECC \
  --in reads.fastq.gz \
  --out results/sample_01 \
  --prefix sample_01 \
  --threads 16 \
  --ref reference.fa \
  --mapq 30
```

Short-option equivalent:

```bash
coreid call_ECC \
  -i reads.fastq.gz \
  -o results/sample_01 \
  -p sample_01 \
  -t 16 \
  -r reference.fa \
  -q 30
```

### Parameters

| Option | Description | Default |
| --- | --- | --- |
| `-i`, `--in` | Input long-read FASTA/FASTQ file | Required |
| `-o`, `--out` | Output directory | Required |
| `-p`, `--prefix` | Sample/output prefix | `ECCTRY` |
| `-t`, `--threads` | Number of worker processes | All available CPUs |
| `-r`, `--ref` | Reference genome FASTA | Required |
| `-q`, `--mapq` | Minimum mapping quality used during raw-read fragment reconstruction | `30` |

## Output

The final result is written to:

```text
<output_directory>/<prefix>_eccdna_information.txt
```

| Column | Description |
| --- | --- |
| `id` | Sequential eccDNA candidate identifier. |
| `merge_region` | Representative genomic interval or comma-separated interval series, formatted as `chromosome:start-end_strand`. |
| `merge_len` | Sum of the representative interval lengths. |
| `num_region` | Number of genomic intervals in the reported structure. |
| `eccdna_status` | Candidate type; currently reported as `ccs_cyclic`. |
| `numreads` | Number of supporting reads assigned to the final clustered structure. |

The `tmp/` directory retains only two reusable CCS cache files:

```text
tmp/<prefix>.ccs.fa
tmp/<prefix>.raw.fa
```

When both files are present, a subsequent run with the same output directory and prefix reloads them instead of recalculating CCS sequences. Other candidate, filtering, tree, and serialized intermediate files are not written.

## Multi-sample integration

`merge_eccdna.py` merges eccDNA calls across samples and methods. It accepts one method per run:

- CoReID
- CReSIL
- NanoCircle
- Flec
- FLED

The script parses the method-specific result format, normalizes circular coordinates, clusters structurally similar calls across samples, and produces either an eccDNA-by-sample coordinate matrix or a gene-by-sample read-support matrix.

### Sample manifest

The sample manifest has **no header**. Each non-comment row contains three tab-separated fields:

```text
SampleID    file_role    path
```

Paths may be absolute or relative to the directory containing the manifest.

#### CoReID, CReSIL, and Flec

These methods use one result file per sample. Set `file_role` to `main`:

```text
ERR12717804    main    /storage/project/CoReID/ERR12717804_eccdna_information.txt
ERR12717805    main    /storage/project/CoReID/ERR12717805_eccdna_information.txt
ERR12717806    main    /storage/project/CoReID/ERR12717806_eccdna_information.txt
```

For CReSIL or Flec, replace the paths with that method's result files and use the corresponding `--method` value. If Flec files do not contain a header, add `--header F` to the command.

#### NanoCircle and FLED

Each sample must be listed twice: once with `simple` and once with `complex`:

```text
ERR12717804    simple     /storage/project/FLED/ERR12717804.simple.txt
ERR12717804    complex    /storage/project/FLED/ERR12717804.complex.txt
ERR12717805    simple     /storage/project/FLED/ERR12717805.simple.txt
ERR12717805    complex    /storage/project/FLED/ERR12717805.complex.txt
```

Do not combine outputs from different calling methods within one manifest. Run `coreid-merge` once per method.

### Commands

#### eccDNA-by-sample coordinate matrix

```bash
coreid-merge \
  --method CoReID \
  --input coreid_samples.tsv \
  --level eccDNA \
  --overlap 0.9 \
  --threads 8 \
  --prefix cohort \
  --out-dir results
```

This produces:

```text
results/cohort_CoReID_eccDNA_matrix.tsv
```

#### Gene-by-sample read-support matrix

```bash
coreid-merge \
  --method CoReID \
  --input coreid_samples.tsv \
  --level gene \
  --annotation /storage/reference/hg38_geneSymbol.bed \
  --overlap 0.9 \
  --threads 8 \
  --prefix cohort \
  --out-dir results
```

This produces:

```text
results/cohort_CoReID_gene_matrix.tsv
```

For NanoCircle or FLED, use `--method NanoCircle` or `--method FLED` with a manifest containing paired `simple`/`complex` files. Other command-line options are unchanged.

## Notes

- Coordinates and representative boundaries are inferred from clustered read-level alignments.
- The workflow is designed for long reads containing detectable tandem repeat units. Molecules without usable tandem-repeat/CCS evidence are outside the current detection scope.
- The current code is research software. Please validate calls with orthogonal evidence and appropriate controls before biological or clinical interpretation.

## Citation

If you use CoReID, please cite the accompanying manuscript. Citation details will be added upon publication.

