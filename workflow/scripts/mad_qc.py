import sys
import argparse
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

def organize_exp(file1, file2, exp_col_idx, id_col_idx):
    """
    Reads two tables, checks uniqueness of IDs, and matches rows (inner join)
    based on the ID column.
    
    Args:
        exp_col_idx (int): 0-based index for expression value (R script used 7 -> Python 6)
        id_col_idx (int): 0-based index for gene/transcript ID (R script used 1 -> Python 0)
    """
    # Read files (assuming tab-separated based on R's read.delim)
    df1 = pd.read_csv(file1, sep="\t")
    df2 = pd.read_csv(file2, sep="\t")

    # 1. Validate Uniqueness (replicating R: stop("not unique"))
    # We access columns by integer position using iloc
    id1 = df1.iloc[:, id_col_idx]
    id2 = df2.iloc[:, id_col_idx]

    if not id1.is_unique:
        sys.exit(f"Error: IDs in {file1} are not unique.")
    if not id2.is_unique:
        sys.exit(f"Error: IDs in {file2} are not unique.")

    # 2. Extract specific columns and rename for merging
    # We create temporary dataframes with just ID and Exp
    sub1 = df1.iloc[:, [id_col_idx, exp_col_idx]].copy()
    sub2 = df2.iloc[:, [id_col_idx, exp_col_idx]].copy()
    
    # Standardize column names
    sub1.columns = ["id", "rep1"]
    sub2.columns = ["id", "rep2"]

    # 3. Align tables (replicating R's match logic)
    # R script implies strict matching. An inner join ensures we only keep
    # genes present in both files, aligned by ID.
    merged = pd.merge(sub1, sub2, on="id", how="inner")
    
    # R script checks: if(!identical(geneid1, geneid2)) -> implied by merge
    return merged

def main():
    parser = argparse.ArgumentParser(
        description="Calculate MAD QC metrics between two RSEM gene expression files"
    )
    parser.add_argument("file1", help="First RSEM genes.results file")
    parser.add_argument("file2", help="Second RSEM genes.results file")
    parser.add_argument(
        "--json",
        help="Output JSON file path (default: stdout)",
        default=None
    )
    parser.add_argument(
        "--plot",
        help="Output plot file path (default: MAplot.png)",
        default="MAplot.png"
    )
    
    args = parser.parse_args()
    
    file1 = args.file1
    file2 = args.file2
    json_out = args.json
    plot_out = args.plot

    # --- Configuration matches R script defaults ---
    # R used col 7 for Data (FPKM) -> Python index 6
    # R used col 1 for ID (Gene)   -> Python index 0
    EXP_COL = 6 
    ID_COL = 0
    A_CUTOFF = 0

    # 1. Organize Expression Data
    try:
        reps = organize_exp(file1, file2, EXP_COL, ID_COL)
    except Exception as e:
        sys.exit(str(e))

    # 2. Filter Zeros
    # R: nozero <- which(reps$rep1!=0 | reps$rep2!=0)
    # Keep rows where at least one replicate is non-zero
    reps_part = reps[(reps["rep1"] != 0) | (reps["rep2"] != 0)].copy()

    if reps_part.empty:
        sys.exit("Error: No non-zero expression data found.")

    # 3. Log Transformation
    # R: logrep1 <- log2(reps$rep1[nozero])
    # Note: If a value is 0, log2(0) = -inf. Numpy will warn, but we suppress it
    # because we filter these out later via the A_CUTOFF.
    with np.errstate(divide='ignore'):
        logrep1 = np.log2(reps_part["rep1"])
        logrep2 = np.log2(reps_part["rep2"])

    # 4. Calculate M and A
    # A = Average Log Expression
    # M = Log Ratio
    A = (logrep1 + logrep2) / 2
    M = logrep1 - logrep2

    # 5. Apply Cutoff for Statistics
    # R: M[A > Acutoff]
    # Note: Since Acutoff=0, any -inf values (from 0 counts) are automatically filtered out here.
    mask = A > A_CUTOFF
    
    M_filtered = M[mask]
    logrep1_filtered = logrep1[mask]
    logrep2_filtered = logrep2[mask]

    # 6. Compute Statistics
    stats = {}
    
    if len(M_filtered) > 0:
        # MAD of log ratios
        # Formula: median(|M|) * 1.4826
        # Note: This formula assumes the median of M is roughly 0.
        mad_val = np.median(np.abs(M_filtered)) * 1.4826
        stats["MAD of log ratios"] = round(mad_val, 3)

        # Pearson correlation
        pearson = np.corrcoef(logrep1_filtered, logrep2_filtered)[0, 1]
        stats["Pearson correlation"] = pearson # R script prints raw float

        # Spearman correlation
        spearman = logrep1_filtered.corr(logrep2_filtered, method='spearman')
        stats["Spearman correlation"] = spearman # R script prints raw float

        # SD of log ratios
        # Formula: sqrt(mean(M^2))
        # Note: This is Root Mean Square, not standard deviation (which centers the mean).
        # We use this exact formula to match the R script.
        sd_val = np.sqrt(np.mean(M_filtered**2))
        stats["SD of log ratios"] = round(sd_val, 3)
    else:
        stats = {"warning": "No genes passed the expression cutoff (A > 0)"}

    # 7. Output JSON
    json_str = json.dumps(stats, indent=4)
    if json_out:
        with open(json_out, 'w') as f:
            f.write(json_str)
    else:
        print(json_str)

    # 8. Generate MA Plot
    # R: bitmap("MAplot.png"); plot(A,M)
    # We must filter out -inf values for matplotlib to work
    plot_mask = np.isfinite(A) & np.isfinite(M)
    
    plt.figure(figsize=(6, 6))
    plt.scatter(A[plot_mask], M[plot_mask], s=1, alpha=0.4, color='black')
    plt.title("MA Plot")
    plt.xlabel("A (Average Log2 Expression)")
    plt.ylabel("M (Log2 Fold Change)")
    
    # Optional: Add red line at 0 (common in MA plots)
    plt.axhline(0, color='red', linewidth=1, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(plot_out, dpi=150)
    plt.close()  # Close the figure to free memory

if __name__ == "__main__":
    main()