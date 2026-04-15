# fastoma-wrapper

A lightweight wrapper around the [FastOMA](https://github.com/DessimozLab/FastOMA) Nextflow pipeline for orthology inference from proteome data.

**No other repositories required.** This repo contains one script (`run_fastoma.sh`) that handles input pre-processing, LUCA database management, and FastOMA execution end-to-end.

---

## What this does

FastOMA infers orthology groups (HOGs) across species from protein sequences. This wrapper adds:

- Protein length filtering before FastOMA runs
- GFF ID prefixing to prevent cross-species collisions
- Automated LUCA database download and caching (`~/.fastoma/omamerdb.h5`)
- Optional NCBI data download via `ncbi-datasets-cli`
- Nextflow config generation tuned to available resources

FastOMA itself runs via `nextflow run DessimozLab/FastOMA` — no local clone of FastOMA is needed.

---

## Installation

### 1. Clone this repo

```bash
git clone https://github.com/chazhyseni/fastoma-wrapper.git
cd fastoma-wrapper
chmod +x run_fastoma.sh
```

### 2. Create the conda environment

```bash
mamba env create -f environment.yml
mamba activate fastoma
```

### 3. (Optional) Pre-download the LUCA database

The script downloads it automatically on first run (~7 GB), but you can pre-download:

```bash
mkdir -p ~/.fastoma
wget -O ~/.fastoma/omamerdb.h5 https://omabrowser.org/All/LUCA.h5
```

---

## Dependencies

| Tool | Version | Notes |
|------|---------|-------|
| Nextflow | 25.04.4+ | Included in `environment.yml` |
| Java | 11+ | Included in `environment.yml` |
| Python | 3.8+ | Included in `environment.yml` |
| biopython | any | Included in `environment.yml` |
| pandas | any | Included in `environment.yml` |
| wget | any | System package |
| ncbi-datasets-cli | any | Only needed with `--auto-download` |

---

## Quick start

### Using pre-existing protein and GFF files

```bash
./run_fastoma.sh \
    --config example_configs/primates.json \
    --proteins-dir /path/to/proteomes \
    --gff-dir /path/to/annotations \
    --output-dir results/primates_v1 \
    --threads 32 --memory 256G
```

Protein FASTAs must be named `<species_key>.fa` (e.g. `human.fa`).
GFF3 files must be named `<species_key>.gff`.
Species keys must match the keys in `SPECIES_LIST` in your config.

### Using CLI flags instead of a config file

```bash
./run_fastoma.sh \
    --species "human:Homo sapiens,chimp:Pan troglodytes,mouse:Mus musculus" \
    --tree "((human:6,chimp:6):75,mouse:81);" \
    --proteins-dir /path/to/proteomes \
    --gff-dir /path/to/annotations \
    --output-dir results/my_run
```

### Auto-downloading data from NCBI

```bash
./run_fastoma.sh \
    --config example_configs/primates.json \
    --auto-download \
    --output-dir results/primates_v1
```

Requires `ncbi-datasets-cli`: `conda install -c bioconda ncbi-datasets-cli`

---

## Config file format

```json
{
    "ANALYSIS_NAME": "primates_v1",
    "SPECIES_LIST": {
        "human":   "Homo sapiens",
        "chimp":   "Pan troglodytes",
        "gorilla": "Gorilla gorilla",
        "mouse":   "Mus musculus"
    },
    "SPECIES_TREE": "(((human:6,chimp:6):2,gorilla:8):72,mouse:80);"
}
```

- **`SPECIES_LIST`** keys must match filenames in `--proteins-dir` and `--gff-dir`
- **`SPECIES_TREE`** must be valid Newick format ending with `;`, with branch lengths

See [`example_configs/primates.json`](example_configs/primates.json) for a working example.

---

## All options

```
INPUT:
    --config FILE             JSON config file
    --species "k:Name,..."    Comma-separated "key:Scientific name" pairs
    --tree NEWICK             Newick tree with branch lengths

DATA SOURCES:
    --proteins-dir DIR        Per-species protein FASTAs (species.fa)
    --gff-dir DIR             Per-species GFF3 files (species.gff)
    --auto-download           Download from NCBI (requires ncbi-datasets-cli)

OUTPUT:
    --output-dir DIR          Required. Where results are written.
    --work-dir DIR            Temp files (default: OUTPUT_DIR/.work)
    --luca-db PATH            LUCA database path (default: ~/.fastoma/omamerdb.h5)
    --log-file PATH           Log file (default: OUTPUT_DIR/run.log)

RESOURCES:
    --threads N               CPU threads (default: nproc)
    --memory SIZE             Memory, e.g. 128G (default: 128G)
    --min-protein-length N    Discard proteins shorter than N aa (default: 50)

FASTOMA / NEXTFLOW:
    --min-sequence-length N   default: 60
    --filter-gap-row F        default: 0.4
    --filter-gap-col F        default: 0.6
    --nr-repr-per-hog N       default: 3
    --fastoma-revision REV    FastOMA git tag or branch (default: main)

EXECUTION:
    --analysis-name NAME      Label for this run
    --dry-run                 Validate inputs, exit without running
    --force-reprocess         Re-run even if output already exists
    --verbose                 Print species tree and proteome list before run
```

---

## Output

```
OUTPUT_DIR/
├── OrthologousGroups.tsv     Main result — one HOG per row, one species per column
├── run_summary.txt           Run statistics
├── run.log                   Full execution log
└── logs/
    ├── timeline.html         Nextflow task timeline
    ├── report.html           Nextflow run report
    └── trace.txt             Per-task resource usage
```

### OrthologousGroups.tsv format

```
OrthologyGroupID    human:ENST00000123    chimp:ENST00000456    mouse:ENST00000789
HOG:0000001         protein_A             protein_B             protein_C
HOG:0000002         protein_D             -                     protein_E
```

Dash (`-`) indicates the species has no member in that group.

---

## LUCA database

The LUCA database (`omamerdb.h5`, ~7 GB) is the orthology reference used by FastOMA. The script downloads it on first run and caches it at `~/.fastoma/omamerdb.h5`. All subsequent runs reuse the cached copy.

To use a different location: `--luca-db /your/preferred/path/omamerdb.h5`

---

## System requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| OS | Linux (Ubuntu 20.04+) or macOS | Linux |
| RAM | 16 GB | 128 GB |
| Storage | 10 GB + 7 GB (LUCA) | 100 GB |
| CPU | 4 cores | 32+ cores |

---

## Troubleshooting

| Problem | Fix |
|---------|-----|
| `Java heap space` error | `export NXF_OPTS="-Xms2G -Xmx8G"` |
| Nextflow resume fails | Delete `OUTPUT_DIR/.work/work/` and rerun |
| Low protein count after filtering | Use `--min-protein-length 30` |
| LUCA download fails or incomplete | Manually download: `wget -O ~/.fastoma/omamerdb.h5 https://omabrowser.org/All/LUCA.h5` |
| Species tree validation error | Tree must be valid Newick ending with `;` and include branch lengths |
| `biopython` not found | `pip install biopython pandas` |
| `ncbi-datasets-cli` not found | `conda install -c bioconda ncbi-datasets-cli` |

---

## Citation

If you use FastOMA in your research, please cite:

> Altenhoff AM et al. (2024). FastOMA: scalable orthology inference using the OMA Hierarchical Orthologous Groups framework. *bioRxiv*. https://github.com/DessimozLab/FastOMA

---

## License

MIT
