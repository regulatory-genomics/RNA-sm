# RNA-seq Analysis Pipeline

A comprehensive Snakemake-based RNA-seq data analysis pipeline for bulk RNA sequencing data.

## Overview

This pipeline processes paired-end RNA-seq data from raw FASTQ files through quality control, trimming, alignment, and quantification. It is designed to handle multiple technical replicates (runs) per biological sample and supports biological replicates for reproducibility analysis.

## Features

- **Quality Control**: FastQC and fastp for adapter trimming and quality filtering
- **Alignment**: STAR aligner with gene counting
- **Quantification**: Gene-level expression quantification from STAR ReadsPerGene output
- **Technical Replicate Handling**: Automatic merging of multiple runs per sample
- **Biological Replicate Support**: Support for replicate groups (via `replicate_name` column)
- **Comprehensive Reporting**: MultiQC report aggregating all QC metrics

## Pipeline Structure

```
RNA-sm/
├── config/
│   ├── config.yaml              # Main configuration file
│   └── test_info.csv            # Sample information file
├── workflow/
│   ├── Snakefile                # Main workflow file
│   ├── rules/
│   │   ├── common.smk           # Common functions and setup
│   │   ├── qc.smk               # Quality control rules
│   │   ├── align.smk            # Alignment rules
│   │   └── count.smk            # Quantification rules
│   ├── scripts/
│   │   └── count-matrix.py      # Count matrix generation script
│   └── env/
│       ├── fastp.yaml           # fastp environment
│       └── pandas.yaml          # pandas environment
└── schemas/
    ├── config.schema.yaml       # Config validation schema
    └── samples.schema.yaml      # Sample table validation schema
```

## Requirements

- Snakemake >= 8.0.0
- Conda/Mamba for environment management
- Python >= 3.8

## Sample Information File

The pipeline uses a CSV file (`config/test_info.csv`) with the following columns:

| Column         | Required | Description                                    |
|----------------|----------|------------------------------------------------|
| sample_name    | Yes      | Unique sample identifier                       |
| run            | Yes      | Run number (technical replicate), integer      |
| R1             | Yes      | Path to forward read FASTQ file                |
| R2             | Yes      | Path to reverse read FASTQ file                |
| replicate_name | No       | Biological replicate group name                |
| strandedness   | No       | Library strandedness: none/forward/reverse     |
| passqc         | No       | QC pass flag: 1=pass, 0=fail                   |

### Example Sample File

```csv
sample_name,run,R1,R2,replicate_name
test1,1,/path/to/test1_1_1.fq,/path/to/test1_1_2.fq,test1
test1,2,/path/to/test1_2_1.fq,/path/to/test1_2_2.fq,test1
test2,1,/path/to/test2_1_1.fq,/path/to/test2_1_2.fq,test2
```

### Key Features of Sample Handling

- **Technical Replicates**: Multiple runs per sample (same `sample_name`, different `run` numbers) are processed independently and then merged
- **Biological Replicates**: Samples sharing the same `replicate_name` are grouped for reproducibility analysis
- **QC Filtering**: If `passqc` column exists, only samples with at least one run passing QC (passqc=1) are included in downstream analysis

## Configuration

Edit `config/config.yaml` to customize:

```yaml
# Path to sample information CSV file
samples: config/test_info.csv

# Output directory path
result_path: results

# Reference genome settings
ref:
  species: homo_sapiens
  release: 100
  build: GRCh38

# Resource paths
resources:
  star_index: /path/to/star/index
  gtf: /path/to/annotation.gtf
  fasta: /path/to/genome.fa

# STAR alignment parameters
params:
  star:
    align:
      extra: "--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"
    index:
      extra: ""
```

## Running the Pipeline

### Dry Run (Check Execution Plan)

```bash
cd workflow
snakemake --use-conda -n
```

### Execute Pipeline

```bash
cd workflow
snakemake --use-conda --cores 8
```

### Execute with Cluster Support

```bash
cd workflow
snakemake --use-conda --cluster "sbatch -p partition -c {threads}" --jobs 10
```

## Output Structure

```
results/
├── qc/
│   ├── fastqc/                  # FastQC reports
│   ├── fastp/                   # Fastp trimming reports
│   └── multiqc_report.html      # Aggregated QC report
├── trimmed/                     # Trimmed FASTQ files (temp)
├── align/
│   ├── {sample_run}_sortedByCoord.out.bam    # Per-run alignments
│   ├── {sample_run}_ReadsPerGene.out.tab     # Per-run counts
│   ├── {sample_run}_Log.final.out            # STAR log files
│   ├── merged/                               # Merged BAMs (if multiple runs)
│   └── unmapped/                             # Unmapped reads
└── counts/
    └── all_counts.tsv           # Gene expression count matrix
```

## Pipeline Workflow

1. **Quality Control (FastQC)**: Initial quality assessment of raw FASTQ files
2. **Trimming (fastp)**: Adapter trimming and quality filtering
3. **Alignment (STAR)**: Map reads to reference genome with gene counting
4. **Merging**: Combine multiple runs per sample (if applicable)
5. **Quantification**: Generate gene count matrix from STAR output
6. **Reporting (MultiQC)**: Aggregate all QC metrics into comprehensive report

## Key Functions in common.smk

- `get_runs_for_sample(sample_name)`: Get all run identifiers for a sample
- `get_units_fastqs(wildcards)`: Get FASTQ file paths for a specific run
- `get_all_fastqs_for_sample(sample_name)`: Get all FASTQ files for a sample
- `get_replicate_samples(replicate_name)`: Get samples in a replicate group
- `get_samples_passing_qc()`: Filter samples by QC status
- `get_strandedness(sample_run)`: Get library strandedness information

## Technical Details

### STAR Alignment Parameters

Default STAR parameters:
- `--outSAMtype BAM SortedByCoordinate`: Output coordinate-sorted BAM
- `--quantMode GeneCounts`: Generate per-gene read counts
- `--sjdbGTFfile`: Use GTF for splice junction detection
- `--outReadsUnmapped Fastx`: Output unmapped reads

### Count Matrix Generation

The count matrix is generated from STAR's `ReadsPerGene.out.tab` files:
- Column selection based on strandedness setting:
  - `none`: Column 2 (unstranded)
  - `forward`: Column 3 (first-strand)
  - `reverse`: Column 4 (second-strand, typical for Illumina TruSeq)
- Technical replicates are automatically collapsed by summing counts

## Troubleshooting

### Issue: Duplicate sample_run identifiers

**Error**: `Duplicate (sample_name + run) pairs found in sample sheet!`

**Solution**: Ensure each combination of `sample_name` and `run` is unique in your CSV file.

### Issue: Missing FASTQ files

**Error**: Input files not found

**Solution**: 
- Verify file paths in CSV are correct and absolute
- Ensure files are accessible from the execution environment
- Check file permissions

### Issue: STAR index not found

**Solution**: 
- Build STAR index first (see `star_index` rule)
- Update `resources.star_index` path in config.yaml

## Citation

If you use this pipeline, please cite:

- **Snakemake**: Mölder et al. (2021), Sustainable data analysis with Snakemake
- **STAR**: Dobin et al. (2013), STAR: ultrafast universal RNA-seq aligner
- **fastp**: Chen et al. (2018), fastp: an ultra-fast all-in-one FASTQ preprocessor
- **FastQC**: Andrews (2010), FastQC: a quality control tool for high throughput sequence data
- **MultiQC**: Ewels et al. (2016), MultiQC: summarize analysis results for multiple tools

## License

This pipeline is provided as-is for research purposes.

## Contact

For questions or issues, please open an issue on the repository or contact the pipeline maintainer.