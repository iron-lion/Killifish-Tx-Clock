"""
Gene ID mapping: query Ensembl ENSNFUG IDs → Atlas NCBI gene names.

Background
----------
- Atlas uses NCBI RefSeq annotation: LOC107XXXXXX (unannotated) + lowercase gene symbols
- Query data uses Ensembl annotation: ENSNFUG00015XXXXXX IDs + mixed-case gene names
- Three-layer mapping (covers ~12,500 genes):
    1. Direct: lowercase(query gene_name) → Atlas gene name
    2. BioMart fallback: Ensembl ENSNFUG → external_gene_name → Atlas (for unmatched named genes)
    3. LOC: NCBI gene2ensembl ENSNFUG → GeneID 107XXXXXX → 'LOC107XXXXXX' in Atlas

Pre-built mapping file: data_matrices/query_to_atlas_gene_mapping.csv
To rebuild: python gene_mapping.py  (or GeneMapper.build_and_save(...))
"""

import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
MAPPING_FILE = REPO_ROOT / "data_matrices" / "query_to_atlas_gene_mapping.csv"
QUERY_DIR = REPO_ROOT / "query_data"
ATLAS_COUNTS_FILE = REPO_ROOT / "data_matrices" / "GSE308970_Counts_Atlas_allbatches_merged_v3.csv"
GENE2ENSEMBL_FILE = REPO_ROOT / "data_matrices" / "ncbi_gene2ensembl_nfurzeri.csv"


class GeneMapper:
    """Convert query ENSNFUG gene IDs to Atlas gene names.

    Parameters
    ----------
    mapping_file : Path
        Pre-built CSV with columns [ensembl_gene_id, atlas_gene].
        Generate with GeneMapper.build_and_save() if missing.
    """

    def __init__(self, mapping_file: Path = MAPPING_FILE):
        df = pd.read_csv(mapping_file)
        self._map = df.set_index("ensembl_gene_id")["atlas_gene"]

    # ------------------------------------------------------------------
    # Conversion
    # ------------------------------------------------------------------

    def convert(self, counts: pd.DataFrame) -> pd.DataFrame:
        """Map ENSNFUG index to Atlas gene names; keep only mapped genes.

        Parameters
        ----------
        counts : pd.DataFrame
            genes × samples, index = ensembl_gene_id (ENSNFUG...).

        Returns
        -------
        pd.DataFrame
            genes × samples, index = Atlas gene names.
            Genes not in the mapping are dropped.
            Duplicate Atlas names (rare) keep the first occurrence.
        """
        mapped_idx = self._map.reindex(counts.index).dropna()
        result = counts.loc[mapped_idx.index].copy()
        result.index = mapped_idx.values
        result = result[~result.index.duplicated(keep="first")]
        return result

    def summary(self) -> None:
        print(f"Gene mapping: {len(self._map):,} ENSNFUG IDs → Atlas gene names")

    # ------------------------------------------------------------------
    # Build mapping from scratch
    # ------------------------------------------------------------------

    @staticmethod
    def build_and_save(
        query_dir: Path = QUERY_DIR,
        atlas_counts_path: Path = ATLAS_COUNTS_FILE,
        gene2ensembl_path: Path = GENE2ENSEMBL_FILE,
        output_path: Path = MAPPING_FILE,
    ) -> pd.DataFrame:
        """Build and save the ENSNFUG → Atlas gene name mapping.

        Requires internet access (BioMart query for Method 2).

        Returns
        -------
        pd.DataFrame
            Mapping with columns [ensembl_gene_id, atlas_gene].
        """
        from pybiomart import Server

        # Load all unique ENSNFUG IDs across query files
        all_genes = pd.DataFrame()
        for f in sorted(Path(query_dir).glob("*.xlsx")):
            df = pd.read_excel(f)[["ensembl_gene_id", "gene_name"]]
            all_genes = pd.concat([all_genes, df])
        all_genes = all_genes.drop_duplicates("ensembl_gene_id").reset_index(drop=True)
        print(f"Unique ENSNFUG IDs in query files: {len(all_genes)}")

        # Load Atlas gene index
        atlas = pd.read_csv(atlas_counts_path, index_col=0, nrows=0)
        atlas_full = pd.read_csv(atlas_counts_path, index_col=0)
        atlas_lower = pd.Series(atlas_full.index.tolist(), index=atlas_full.index.str.lower())

        # Method 1: direct lowercase gene_name
        named = all_genes[all_genes["gene_name"].notna()].copy()
        named["gene_name_lower"] = named["gene_name"].str.lower()
        direct = named[named["gene_name_lower"].isin(atlas_lower.index)].copy()
        direct["atlas_gene"] = direct["gene_name_lower"].map(atlas_lower)
        print(f"Method 1 (direct lowercase): {len(direct)}")

        # Method 2: BioMart for unmatched named genes
        matched_ids = set(direct["ensembl_gene_id"])
        unmatched = named[~named["ensembl_gene_id"].isin(matched_ids)]
        bm_df = (
            Server("http://www.ensembl.org")["ENSEMBL_MART_ENSEMBL"]["nfurzeri_gene_ensembl"]
            .query(attributes=["ensembl_gene_id", "external_gene_name"])
        )
        bm_df.columns = ["ensembl_gene_id", "gene_name_bm"]
        bm_df["gene_name_bm_lower"] = bm_df["gene_name_bm"].str.lower()
        bm_df = bm_df.dropna(subset=["gene_name_bm"])
        unmatched_bm = unmatched.merge(bm_df[["ensembl_gene_id", "gene_name_bm_lower"]], on="ensembl_gene_id", how="left")
        bm_m = unmatched_bm[unmatched_bm["gene_name_bm_lower"].isin(atlas_lower.index)].dropna(subset=["gene_name_bm_lower"]).copy()
        bm_m["atlas_gene"] = bm_m["gene_name_bm_lower"].map(atlas_lower)
        print(f"Method 2 (BioMart fallback): {len(bm_m)}")

        # Method 3: NCBI gene2ensembl LOC107 → ENSNFUG
        g2e = pd.read_csv(gene2ensembl_path)
        loc107 = g2e[(g2e["GeneID"] > 107_000_000) & (g2e["GeneID"] < 108_000_000)].copy()
        loc107["atlas_gene"] = "LOC" + loc107["GeneID"].astype(str)
        loc107 = loc107[loc107["atlas_gene"].isin(atlas_full.index)]
        matched_ids2 = matched_ids | set(bm_m["ensembl_gene_id"])
        loc_m = all_genes.merge(loc107[["Ensembl_gene_identifier", "atlas_gene"]], left_on="ensembl_gene_id", right_on="Ensembl_gene_identifier", how="inner")
        loc_m = loc_m[~loc_m["ensembl_gene_id"].isin(matched_ids2)]
        print(f"Method 3 (LOC107 via gene2ensembl): {len(loc_m)}")

        # Combine and save
        final = pd.concat([
            direct[["ensembl_gene_id", "atlas_gene"]],
            bm_m[["ensembl_gene_id", "atlas_gene"]],
            loc_m[["ensembl_gene_id", "atlas_gene"]],
        ]).drop_duplicates("ensembl_gene_id")
        print(f"Total mapped: {len(final)} / {len(all_genes)} query genes")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        final.to_csv(output_path, index=False)
        print(f"Saved to {output_path}")
        return final


if __name__ == "__main__":
    print("Building gene ID mapping...")
    GeneMapper.build_and_save()
