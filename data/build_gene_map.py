"""
Build a unified gene ID mapping table from three sources:
  1. GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf  -> gtf_gene_id, gtf_gene_name, ncbi_gene_id
  2. ncbi_gene2ensembl_nfurzeri.csv              -> ncbi_gene_id <-> ensembl_gene_id
  3. query_to_atlas_gene_mapping.csv             -> ensembl_gene_id -> atlas_gene

Output: gene_id_map.csv  (one row per unique ensembl_gene_id)
"""

import re
import pandas as pd

GTF      = "GCF_043380555.1_NfurGRZ-RIMD1_genomic.gtf"
G2E_CSV  = "ncbi_gene2ensembl_nfurzeri.csv"
ATL_CSV  = "query_to_atlas_gene_mapping.csv"
OUT_CSV  = "gene_id_map.csv"

# ── 1. Parse GTF gene records ─────────────────────────────────────────────────
gtf_rows = []
with open(GTF) as f:
    for line in f:
        if line.startswith('#') or '\tgene\t' not in line:
            continue
        attrs  = line.strip().split('\t')[8]
        gid    = re.search(r'gene_id "([^"]+)"', attrs)
        gname  = re.search(r'(?<!\w)gene "([^"]+)"', attrs)
        dbxref = re.search(r'db_xref "GeneID:(\d+)"', attrs)
        if gid and gname and dbxref:
            gtf_rows.append({
                "gtf_gene_id":   gid.group(1),
                "gtf_gene_name": gname.group(1),
                "ncbi_gene_id":  int(dbxref.group(1)),
            })

gtf = pd.DataFrame(gtf_rows).drop_duplicates()

# ── 2. Load ncbi <-> ensembl mapping ─────────────────────────────────────────
g2e = pd.read_csv(G2E_CSV, usecols=["GeneID", "Ensembl_gene_identifier"])
g2e.columns = ["ncbi_gene_id", "ensembl_gene_id"]
g2e = g2e[g2e["ensembl_gene_id"] != "-"].drop_duplicates()

# ── 3. Load ensembl -> atlas_gene mapping ────────────────────────────────────
atl = pd.read_csv(ATL_CSV)
atl.columns = ["ensembl_gene_id", "atlas_gene"]
atl = atl.drop_duplicates()

# ── 4. Merge: gtf -> g2e -> atl on ncbi_gene_id / ensembl_gene_id ────────────
# Start from gtf (fullest NCBI annotation), merge in ensembl IDs, then atlas names
merged = (
    gtf
    .merge(g2e, on="ncbi_gene_id", how="left")
    .merge(atl, on="ensembl_gene_id", how="left")
)

# Also include any ensembl IDs in g2e that weren't captured via GTF
g2e_only = g2e[~g2e["ncbi_gene_id"].isin(gtf["ncbi_gene_id"])].copy()
g2e_only["gtf_gene_id"]   = pd.NA
g2e_only["gtf_gene_name"] = pd.NA
g2e_only = g2e_only.merge(atl, on="ensembl_gene_id", how="left")

# And atlas entries not reached via GTF or g2e
atl_only = atl[~atl["ensembl_gene_id"].isin(merged["ensembl_gene_id"].dropna())
               & ~atl["ensembl_gene_id"].isin(g2e_only["ensembl_gene_id"])].copy()
atl_only["gtf_gene_id"]   = pd.NA
atl_only["gtf_gene_name"] = pd.NA
atl_only["ncbi_gene_id"]  = pd.NA

full_map = pd.concat([
    merged[["gtf_gene_id","gtf_gene_name","ncbi_gene_id","ensembl_gene_id","atlas_gene"]],
    g2e_only[["gtf_gene_id","gtf_gene_name","ncbi_gene_id","ensembl_gene_id","atlas_gene"]],
    atl_only[["ensembl_gene_id","atlas_gene"]].assign(gtf_gene_id=pd.NA, gtf_gene_name=pd.NA, ncbi_gene_id=pd.NA),
], ignore_index=True).drop_duplicates()

# ── 5. Save ───────────────────────────────────────────────────────────────────
full_map.to_csv(OUT_CSV, index=False)

print(f"Saved {OUT_CSV}: {len(full_map)} rows")
print(f"\nCoverage (non-null per column):")
for col in full_map.columns:
    print(f"  {col:<25}: {full_map[col].notna().sum():>7}")

# Quick sanity: unique values per column
print(f"\nUnique values per column:")
for col in full_map.columns:
    print(f"  {col:<25}: {full_map[col].nunique():>7}")
