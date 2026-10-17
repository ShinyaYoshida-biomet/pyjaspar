# CLI Reference


## Motif retrieval

### motif-by-id

Fetch a single motif by JASPAR matrix ID.

```bash
pyjaspar motif-by-id MA0095.2

# Different output formats
pyjaspar motif-by-id MA0095.2 --motif-format pfm
pyjaspar motif-by-id MA0095.2 --motif-format transfac

# Metadata instead of matrix
pyjaspar motif-by-id MA0095.2 --metadata
pyjaspar motif-by-id MA0095.2 --metadata --format json

# Specific release
pyjaspar motif-by-id MA0095.2 --release 2020
```

### motifs-by-name

Fetch motifs by transcription factor name.

```bash
pyjaspar motifs-by-name CTCF
pyjaspar motifs-by-name CTCF --metadata --format json
```

### motifs

Fetch motifs with advanced filters.

```bash
# CORE vertebrate motifs
pyjaspar motifs --collection CORE --tax-group Vertebrates

# Save to file
pyjaspar motifs --collection CORE --tax-group Vertebrates -o motifs.jaspar

# PFM format
pyjaspar motifs --collection CORE --motif-format pfm

# Filter by minimum information content
pyjaspar motifs --min-ic 12 --min-length 8
```

### metadata

Export motif metadata.

```bash
pyjaspar metadata --tax-group Vertebrates
pyjaspar metadata --tax-group Vertebrates --format json
pyjaspar metadata --tax-group Vertebrates --format json -o metadata.json
```

## Analysis commands

Requires `pip install pyjaspar[analysis]`.

### scan

Scan a DNA sequence for motif occurrences.

```bash
# Literal DNA sequence
pyjaspar scan ACGTACGTACGTACGT --motif-id MA0139.1 --threshold 0.5

# FASTA file input
pyjaspar scan sequences.fasta --motif-id MA0139.1

# JSON output
pyjaspar scan ACGTACGT --motif-id MA0139.1 --format json

# Forward strand only
pyjaspar scan ACGTACGT --motif-id MA0139.1 --no-reverse
```

### similarity

Compare two motifs using similarity metrics.

```bash
# Pearson correlation (default)
pyjaspar similarity MA0001.1 MA0002.1

# Euclidean distance
pyjaspar similarity MA0001.1 MA0002.1 --metric euclidean

# KL divergence
pyjaspar similarity MA0001.1 MA0002.1 --metric kl

# Find best alignment (Pearson only)
pyjaspar similarity MA0001.1 MA0002.1 --best

# With specific offset
pyjaspar similarity MA0001.1 MA0002.1 --offset 3

# JSON output
pyjaspar similarity MA0001.1 MA0002.1 --best --format json
```

### align

Align two motifs and show the aligned view, rather than just a score. The default method is Matrix Align, the method of the JASPAR website; `--method pearson` uses the column-wise Pearson correlation instead.

```bash
# Matrix Align (default)
pyjaspar align MA0001.1 MA0002.1

# Pearson correlation, with at least 5 overlapping columns
pyjaspar align MA0001.1 MA0002.1 --method pearson --min-overlap 5

# Don't try the reverse complement
pyjaspar align MA0001.1 MA0002.1 --no-reverse-complement

# JSON output (score/offset/orientation only, not the rendered view)
pyjaspar align MA0001.1 MA0002.1 --format json
```

### search

Find the profiles most similar to a motif, like the "Search similar profiles" page of the JASPAR website.

```bash
# Top 10 CORE profiles (latest versions)
pyjaspar search MA0139.1

# Vertebrate profiles only, top 5
pyjaspar search MA0139.1 --tax-group Vertebrates --top 5

# Pearson correlation, all versions, JSON output
pyjaspar search MA0139.1 --method pearson --all-versions --format json
```

### enrichment

Test motif enrichment in foreground vs background sequences.

```bash
pyjaspar enrichment \
    --foreground fg_sequences.fasta \
    --background bg_sequences.fasta \
    --motif-ids MA0001.1,MA0002.1,MA0139.1 \
    --threshold 0.8

# JSON output
pyjaspar enrichment \
    -fg fg.fasta -bg bg.fasta \
    --motif-ids MA0001.1 \
    --format json
```

### infer

Predict which JASPAR profiles a protein binds from its amino acid sequence. The search runs on a JASPAR server (2020, 2022 or 2024 release), so it needs network access.

```bash
pyjaspar infer MSNEQTVLSSLIRDEDGNPAAPRK...

# FASTA file of proteins
pyjaspar infer proteins.fasta --format json
```

## Visualization commands

Requires `pip install pyjaspar[viz]`.

### logo

Draw the sequence logo of a motif to an image file; the format follows the extension.

```bash
pyjaspar logo MA0139.1 -o ctcf.png
pyjaspar logo MA0139.1 -o ctcf.svg
```

## Deep Learning (DL) commands

### dl profile

Fetch a DL profile by ID.

```bash
pyjaspar dl profile DL0001.1
pyjaspar dl profile DL0001.1 --format json
```

### dl model

Fetch a DL model by ID.

```bash
pyjaspar dl model BP000001.1
pyjaspar dl model BP000001.1 --format json
```

### dl motif

Fetch a DL motif pattern by MO* ID.

```bash
pyjaspar dl motif MO000001.1
pyjaspar dl motif MO000001.1 --format json
```

### dl profiles

Search DL profiles.

```bash
pyjaspar dl profiles --tf-name REST
pyjaspar dl profiles --tax-group vertebrates --format json
```

### dl models

Search DL models.

```bash
pyjaspar dl models --tf-name REST
pyjaspar dl models --cell-line K562
pyjaspar dl models --tf-name REST --cell-line K562 --format json
```

## General commands

### releases

List available JASPAR releases.

```bash
pyjaspar releases
pyjaspar releases --latest
```

### collections

List JASPAR collections available in a release. The set of collections
varies by release (e.g. `CORE` and `UNVALIDATED` only from 2022 onward).

```bash
pyjaspar collections
pyjaspar collections -r 2020
```

### cite

Show citation information.

```bash
pyjaspar cite
```

## Output formats

Most commands support `--format json` or `--format tsv` (default).

- **TSV**: Tab-separated values, suitable for piping to `cut`, `awk`, etc.
- **JSON**: Structured output, suitable for `jq` or programmatic use.

### Scripting examples

```bash
# Get motif IDs as a list
pyjaspar metadata --tax-group Vertebrates | cut -f1 | tail -n +2

# Count vertebrate motifs
pyjaspar metadata --format json --tax-group Vertebrates | jq '.count'

# Pipe to downstream tools
pyjaspar motif-by-id MA0095.2 --motif-format pfm > YY1.pfm
```
