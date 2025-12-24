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

## License

This pipeline is provided as-is for research purposes.

## Contact

For questions or issues, please open an issue on the repository or contact the pipeline maintainer.