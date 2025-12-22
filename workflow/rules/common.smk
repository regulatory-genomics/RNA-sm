# ===================================================
#                Common Functions and Setup
# ===================================================

import os
import sys
import pandas as pd
import re
from snakemake.utils import validate

# Validate config schema
validate(config, schema="../../schemas/config.schema.yaml")

# ===================================================
#           Load Sample Annotation Table
# ===================================================

# Load sample information from CSV file
# Expected columns: sample_name, run, R1, R2, replicate_name
annot = pd.read_csv(config["samples"], dtype={"sample_name": str, "run": str})

# Filter out empty rows (rows where sample_name is NaN or empty)
annot = annot.dropna(subset=['sample_name'])
annot = annot[annot['sample_name'].astype(str).str.strip() != '']

# Ensure required columns exist
required_columns = {"sample_name", "run", "R1", "R2"}
missing_required = required_columns.difference(annot.columns)
if missing_required:
    raise ValueError(f"Sample table missing required columns: {', '.join(sorted(missing_required))}.")

# Strip whitespace from path columns
if 'R1' in annot.columns:
    annot['R1'] = annot['R1'].astype(str).str.strip()
if 'R2' in annot.columns:
    annot['R2'] = annot['R2'].astype(str).str.strip()

# Convert data types before validation
# Ensure run column is integer
annot['run'] = pd.to_numeric(annot['run'], errors='raise').astype(int)

# Ensure passqc column is integer if it exists (before validation)
if 'passqc' in annot.columns:
    annot['passqc'] = pd.to_numeric(annot['passqc'], errors='coerce').fillna(0).astype(int)

# Validate sample table schema (after type conversions)
validate(annot, schema="../../schemas/samples.schema.yaml")

# Create a unique identifier for each run: sample_name_run
annot['sample_name'] = annot['sample_name'].astype(str)
annot['sample_run'] = annot['sample_name'] + '_' + annot['run'].astype(str)

# Check for duplicate sample_run identifiers
duplicates = annot[annot.duplicated(subset=['sample_run'], keep=False)]
if not duplicates.empty:
    sys.stderr.write("\nError: Duplicate (sample_name + run) pairs found in sample sheet!\n")
    sys.stderr.write("The following rows are not unique:\n")
    sys.stderr.write(str(duplicates[['sample_name', 'run', 'sample_run']]) + "\n\n")
    sys.exit(1)

# Set sample_run as index for easy lookup
annot = annot.set_index('sample_run', drop=False)

# Get unique sample names (for downstream analysis and merging runs)
samples = annot.reset_index(drop=True).drop_duplicates(subset='sample_name', keep='first').set_index("sample_name").to_dict(orient="index")

# ===================================================
#              Helper Functions (defined early for use in setup)
# ===================================================

def get_runs_for_sample(sample_name):
    """Get all sample_run identifiers for a given sample name"""
    sample_runs = annot[annot['sample_name'] == sample_name].index.tolist()
    return sample_runs if len(sample_runs) > 0 else []

def get_replicate_samples(replicate_name):
    """Get all unique sample names for a given replicate_name"""
    if 'replicate_name' not in annot.columns:
        return []
    samples_in_rep = annot['sample_name'][annot['replicate_name'] == replicate_name].unique().tolist()
    return samples_in_rep if len(samples_in_rep) > 1 else []

# Get replicate sample groups (for reproducibility analysis if replicate_name is provided)
if 'replicate_name' in annot.columns:
    all_replicate_names = (
        annot['replicate_name']
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )
    # Only keep replicate names with >1 associated sample
    replicate_samples = [rep for rep in all_replicate_names if len(get_replicate_samples(rep)) > 1]
else:
    replicate_samples = []

# Set wildcard constraints
wildcard_constraints:
    sample="|".join([re.escape(s) for s in samples.keys()]),
    sample_run="|".join([re.escape(sr) for sr in annot.index.tolist()]),

# ===================================================
#              Additional Helper Functions
# ===================================================

def get_units_fastqs(wildcards):
    """Get fastq files for a sample_run"""
    u = annot.loc[wildcards.sample_run]
    fq1 = str(u["R1"]).strip()
    fq2 = str(u["R2"]).strip() if not pd.isna(u["R2"]) else None
    
    # Handle empty strings or 'nan' strings
    if fq1 == '' or fq1.lower() == 'nan':
        raise ValueError(f"R1 path is empty for sample_run {wildcards.sample_run}")
    
    # Ensure absolute paths
    if fq1 and not os.path.isabs(fq1):
        fq1 = os.path.abspath(fq1)
    if fq2 and fq2 != '' and fq2.lower() != 'nan':
        if not os.path.isabs(fq2):
            fq2 = os.path.abspath(fq2)
    else:
        fq2 = None
    
    # Validate that files exist (provide helpful error if not)
    if not os.path.exists(fq1):
        # Check for common alternative extensions
        alt_extensions = ['.fq.gz', '.fastq', '.fastq.gz', '.fq']
        base_path = fq1.rsplit('.', 1)[0] if '.' in fq1 else fq1
        alternatives = [base_path + ext for ext in alt_extensions if base_path + ext != fq1]
        existing_alternatives = [alt for alt in alternatives if os.path.exists(alt)]
        
        error_msg = f"Input file not found for sample_run {wildcards.sample_run}:\n"
        error_msg += f"  Expected: {fq1}\n"
        if existing_alternatives:
            error_msg += f"  Found alternatives: {', '.join(existing_alternatives)}\n"
            error_msg += f"  Please update the CSV file to use one of these paths."
        else:
            error_msg += f"  File does not exist. Please check the path in the CSV file."
        raise FileNotFoundError(error_msg)
    
    if fq2 and not os.path.exists(fq2):
        # Check for common alternative extensions
        alt_extensions = ['.fq.gz', '.fastq', '.fastq.gz', '.fq']
        base_path = fq2.rsplit('.', 1)[0] if '.' in fq2 else fq2
        alternatives = [base_path + ext for ext in alt_extensions if base_path + ext != fq2]
        existing_alternatives = [alt for alt in alternatives if os.path.exists(alt)]
        
        error_msg = f"Input file not found for sample_run {wildcards.sample_run}:\n"
        error_msg += f"  Expected: {fq2}\n"
        if existing_alternatives:
            error_msg += f"  Found alternatives: {', '.join(existing_alternatives)}\n"
            error_msg += f"  Please update the CSV file to use one of these paths."
        else:
            error_msg += f"  File does not exist. Please check the path in the CSV file."
        raise FileNotFoundError(error_msg)
    
    return [fq1, fq2]

def get_all_fastqs_for_sample(sample_name):
    """Get all R1 and R2 fastq files for a sample (all runs combined)"""
    sample_runs = get_runs_for_sample(sample_name)
    r1_files = []
    r2_files = []
    for sr in sample_runs:
        u = annot.loc[sr]
        r1 = u["R1"]
        if not pd.isna(r1):
            r1_files.append(str(r1))
        r2 = u["R2"]
        if not pd.isna(r2):
            r2_files.append(str(r2))
    return r1_files, r2_files

def normalize_strandedness(val):
    """Normalize strandedness value to expected format"""
    if pd.isna(val) or val == '':
        return "none"
    strand_val = str(val).strip().lower()
    # Map common variations to expected values for count-matrix.py
    if strand_val in ['', 'nan', 'none', 'n', 'unstranded', 'no']:
        return "none"
    elif strand_val in ['forward', 'f', 'yes', 'y', 'first']:
        return "yes"
    elif strand_val in ['reverse', 'r', 'second']:
        return "reverse"
    else:
        # Return as-is if it's already a valid value
        return strand_val

def get_strandedness(sample_run=None):
    """
    Get strandedness information for RNA-seq samples.
    Returns 'none', 'yes', or 'reverse' (matching count-matrix.py expectations).
    If strandedness column exists in annot, use it; otherwise default to 'none'.
    """
    
    if "strandedness" in annot.columns:
        if sample_run:
            strand_val = annot.loc[sample_run, "strandedness"]
            return normalize_strandedness(strand_val)
        else:
            # Return list of strandedness values for all sample_runs
            return [normalize_strandedness(annot.loc[sr, "strandedness"]) for sr in annot.index]
    else:
        config_strandedness = config.get("params", {}).get("star", {}).get("strandedness", "yes")
        # Default to 'none' if column doesn't exist
        return config_strandedness if sample_run else [config_strandedness] * len(annot.index)

def get_samples_passing_qc():
    """Get sample names that have at least one run passing QC (passqc=1).
    If passqc column doesn't exist, returns all samples (default: all pass QC)."""
    if 'passqc' not in annot.columns:
        # Default: all samples pass QC if column doesn't exist
        return list(samples.keys())
    
    # Get samples where at least one run has passqc == 1
    passing_samples = set()
    for sample_name in samples.keys():
        sample_runs = get_runs_for_sample(sample_name)
        for sample_run in sample_runs:
            passqc_value = annot.loc[sample_run, 'passqc']
            # Handle int, string, or boolean values
            if passqc_value == 1 or str(passqc_value).strip() == '1':
                passing_samples.add(sample_name)
                break  # At least one run passes, include the sample
    
    return list(passing_samples) if passing_samples else []

# ===================================================
#              Output Target Functions
# ===================================================

def get_final_output():
    """Define all final output files for the pipeline"""
    result_path = config.get("result_path", "results")
    
    outputs = []
    
    # Important processed - BAM files for all sample_runs
    for sample_run in annot.index:
        outputs.append(os.path.join(result_path, "Important_processed", "Bam", f"{sample_run}_sortedByCoord.out.bam"))
    
    # Merged BAM files per sample (if multiple runs exist)
    for sample_name in samples.keys():
        runs = get_runs_for_sample(sample_name)
        if len(runs) > 1:
            outputs.append(os.path.join(result_path, "Important_processed", "Bam", f"{sample_name}.bam"))
    
    # Downstream results - Count matrix
    outputs.append(os.path.join(result_path, "downstream_res", "counts", "all_counts.tsv"))
    
    # Report - MultiQC report
    outputs.append(os.path.join(result_path, "Report", "multiqc_report.html"))
    
    return outputs

def get_track_strandedness(sample_run):
    """
    Get strandedness for track generation.
    Returns 'Stranded' or 'Unstranded' for STAR --outWigStrand parameter.
    """
    strand_val = get_strandedness(sample_run)
    # STAR expects 'Stranded' or 'Unstranded' (capitalized)
    if strand_val in ['yes', 'reverse']:
        return "Stranded"
    else:
        return "Unstranded"

def is_stranded(sample_run):
    """Check if sample is stranded"""
    strand_val = get_strandedness(sample_run)
    return strand_val in ['yes', 'reverse']

def get_star_bg_filename(wildcards):
    """
    Maps the requested BigWig suffix to the specific raw STAR output file.
    Implements the logic from the Python script:
    - minusUniq -> str1 (if stranded)
    - plusUniq  -> str2 (if stranded)
    - uniq      -> str1 (if unstranded)
    """
    mapping = {
        # Stranded mappings (matches python script logic)
        "minusUniq": "Signal.Unique.str1.out.bg",
        "minusAll":  "Signal.UniqueMultiple.str1.out.bg",
        "plusUniq":  "Signal.Unique.str2.out.bg",
        "plusAll":   "Signal.UniqueMultiple.str2.out.bg",
        # Unstranded mappings
        "uniq":      "Signal.Unique.str1.out.bg",
        "all":       "Signal.UniqueMultiple.str1.out.bg"
    }
    return mapping[wildcards.suffix]


def get_final_bigwigs(wildcards):
    """
    Aggregation function: determine which BigWig files to build
    for all samples based on their strandedness.

    Returned paths are used by rule `all` in the Snakefile.
    """
    final_files = []

    for sample_name in samples.keys():
        # Get strandedness from first run of the sample
        sample_runs = get_runs_for_sample(sample_name)
        if sample_runs:
            strand_val = get_strandedness(sample_runs[0])
        else:
            strand_val = "none"

        if strand_val in ["yes", "reverse"]:
            # Stranded libraries → 4 tracks
            suffixes = ["minusUniq", "minusAll", "plusUniq", "plusAll"]
        else:
            # Unstranded libraries → 2 tracks
            suffixes = ["uniq", "all"]

        for suff in suffixes:
            final_files.append(
                os.path.join(
                    result_path,
                    "Important_processed",
                    "Track",
                    f"{sample_name}_{suff}.bw",
                )
            )

    return final_files


def get_track_strandedness(sample_name_or_run):
    """
    Get strandedness for track generation.
    Returns 'Stranded' or 'Unstranded' for STAR --outWigStrand parameter.
    Can accept either sample_name or sample_run.
    """
    # If this is a sample_run, use directly
    if sample_name_or_run in annot.index:
        strand_val = get_strandedness(sample_name_or_run)
    else:
        # Treat as sample_name: look up first run
        sample_runs = get_runs_for_sample(sample_name_or_run)
        if sample_runs:
            strand_val = get_strandedness(sample_runs[0])
        else:
            strand_val = "none"

    # STAR expects 'Stranded' or 'Unstranded' (capitalized)
    if strand_val in ["yes", "reverse"]:
        return "Stranded"
    else:
        return "Unstranded"


def get_bam_for_tracks(sample_name):
    """
    Get the BAM file to use for track generation.
    Prefers merged BAM if available, otherwise uses first sample_run BAM.
    """
    sample_runs = get_runs_for_sample(sample_name)
    if len(sample_runs) > 1:
        # Use merged BAM if multiple runs exist
        return os.path.join(result_path, "Important_processed", "Bam", f"{sample_name}.bam")
    else:
        # Use the single run BAM
        return os.path.join(
            result_path,
            "Important_processed",
            "Bam",
            f"{sample_runs[0]}_sortedByCoord.out.bam",
        )