import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),
    
def get_units_fastqs(wildcards):
    u = units.loc[(wildcards.sample, wildcards.unit)]
    return [u["fq1"], u["fq2"]]

def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]
# ---- FastQC helpers ----
from os.path import basename

def _strip_ext(path):
    """
    Strips read pair identifiers and extensions, then replaces the
    final underscore with a hyphen.
    e.g., C9n1_D56_B_1_1.fq.gz -> C9n1_D56_B-1
    """
    name = basename(path)
    
    # --- Part 1: Strip read identifiers and extensions ---
    pattern = r"_R?[12](\.fq|\.fastq)\.gz$"
    stripped_name = re.sub(pattern, "", name)

    if stripped_name != name:
        clean_basename = stripped_name
    else:
        # Fall back to just stripping extensions for non-paired files
        for ext in (".fq.gz", ".fastq.gz"):
            if name.endswith(ext):
                clean_basename = name[:-len(ext)]
                break
        else:  # This runs if the for loop doesn't break
            if "." in name:
                clean_basename = name.rsplit(".", 1)[0]
            else:
                clean_basename = name
    
    # --- Part 2: Replace the last underscore with a hyphen ---
    if '_' in clean_basename:
        parts = clean_basename.rsplit('_', 1)
        return '-'.join(parts)
    else:
        return clean_basename

_fq1 = units["fq1"].dropna().tolist()
_fq2 = units["fq2"].dropna().tolist()

fastq_by_basename = { _strip_ext(p): p for p in (_fq1 + _fq2) }

fastqc_targets = []
for sample, unit in units.index:
    fastqc_targets.append(f"results/fastqc/{sample}_{unit}_1_fastqc.html")
    fastqc_targets.append(f"results/fastqc/{sample}_{unit}_2_fastqc.html")

bam_targets = [
    f"results/align/{sample}-{unit}_sortedByCoord.out.bam"
    for sample, unit in units.index
]
def get_final_output():
    return fastqc_targets + bam_targets + ["results/counts/all.tsv"]
