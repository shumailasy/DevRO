# DevRO v5.0.0

DevRO (Deviant Read Orientation) is a structural variant (SV) discovery toolkit for paired-end whole-genome sequencing data.

This repository now includes a **publication-ready command-line workflow** that wraps and modernizes the legacy DevRO callers with:

- A unified CLI (`devro.py`) with reproducible run metadata.
- Parameterized window size, mapQ, insert size distribution, and read length.
- Scalable chunked processing for large genomes and large-structure variant scans.
- Parallel chunk execution for high-throughput runs.
- Backward-compatible core callers for duplication and inversion signatures.

---

## What changed in v5

### Legacy compatibility + modern controls
Legacy callers are still used for the core algorithmic logic:

- `VariantCaller_dup.pl`
- `VariantCaller_inv.pl`

But they can now be configured at runtime via environment variables set by the new CLI:

- `DEVRO_WINDOW_SIZE`
- `DEVRO_READ_LENGTH`
- `DEVRO_MIN_MAPQ`
- `DEVRO_MEAN_INSERT`
- `DEVRO_INSERT_SIGMA`
- `DEVRO_OUTPUT_DIR`

This enables reproducible and tunable runs without editing source code.

### Large-structure variant readiness
To support large genomes and large structural-variant projects:

- Region files can be split into chunks (`--chunk-lines`).
- Chunks can be run in parallel (`--threads`).
- Every run writes JSON metadata suitable for Methods sections and supplementary materials.

---

## Quick start

### 1) Inputs

- **Region file** (`--regions`): tab-separated with 3 columns:
  - `chromosome`, `start`, `chromosome_size`
- **Population config** (`--config`): tab-separated with 2 columns:
  - `group_id`, `bam_path`

Example config row:

```text
1	/path/sampleA.dedup.bam
```

### 2) Run duplication discovery

```bash
python3 devro.py call \
  --mode dup \
  --regions genome.windows.tsv \
  --config populations.tsv \
  --prefix cohortA \
  --output-dir results \
  --window-size 1000 \
  --min-mapq 10 \
  --read-length 150 \
  --mean-insert 400 \
  --insert-sigma 130 \
  --chunk-lines 5000 \
  --threads 8
```

### 3) Run inversion discovery

```bash
python3 devro.py call \
  --mode inv \
  --regions genome.windows.tsv \
  --config populations.tsv \
  --prefix cohortA \
  --output-dir results \
  --chunk-lines 5000 \
  --threads 8
```

---

## Output layout

```text
results/
  dup/
    <prefix>_DupSignatures.Allpop.txt
  inv/
    <prefix>_INVSignatures.Allpop.txt
  <prefix>.<mode>.run-metadata.json
```

`run-metadata.json` includes command settings and timestamps for provenance.

---

## Notes for publication-quality analyses

- Estimate insert-size parameters from each sequencing library (or a representative subset) before final production runs.
- Keep region chunking and parallelism fixed between cohorts for direct comparability.
- Archive metadata JSON files with output tables and downstream scripts.
- Use parser/annotation scripts in this repository for post-processing as needed.

---

## Legacy scripts

The original scripts remain available for historical compatibility:

- `VariantCaller_dup.pl`
- `VariantCaller_inv.pl`
- `VC_dels_ins_refdel.pl`
- `Parser_*.pl`
- `Annotate.get.AvDP.Norm.1K.fixedwindows.wrt.GenomeAvDepth.pl`

Recommended path for new analyses: use `devro.py` as the primary entry point.
