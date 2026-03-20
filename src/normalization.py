"""
Normalizers for aging clock input preparation.

FrequencyNormalizer  — BayesAge 2.0: divide each gene by total reads per sample.
DESeq2Normalizer     — EN / PCR: size-factor normalization via pydeseq2.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference


class FrequencyNormalizer:
    """Frequency count normalization for BayesAge 2.0.

    Transforms raw counts into relative frequencies by dividing each gene's
    count by the total read count for that sample.
    """

    def normalize(self, counts: pd.DataFrame) -> pd.DataFrame:
        """Return frequency-normalized DataFrame (same shape as input).

        Parameters
        ----------
        counts : pd.DataFrame
            Raw counts, genes × samples (index=genes, columns=samples).

        Returns
        -------
        pd.DataFrame
            Relative frequencies in [0, 1], genes × samples.
        """
        totals = counts.sum(axis=0)  # total reads per sample
        freq = counts.div(totals, axis=1)
        return freq


class DESeq2Normalizer:
    """DESeq2 size-factor normalization via pydeseq2.

    Used as input for Elastic Net and PCR clocks.

    Parameters
    ----------
    design_factors : list of str
        Column names in metadata used as design factors (e.g., ['age_days']).
    """

    def __init__(self, design_factors: list[str] = None):
        self.design_factors = design_factors or ["age_days"]
        self._size_factors = None

    def normalize(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        save_path: Path | None = None,
    ) -> pd.DataFrame:
        """Run DESeq2 normalization and return normalized counts.

        Parameters
        ----------
        counts : pd.DataFrame
            Raw integer counts, genes × samples.
        metadata : pd.DataFrame
            Sample metadata indexed by sample ID (same columns as counts).
        save_path : Path, optional
            If provided, save normalized matrix as CSV.

        Returns
        -------
        pd.DataFrame
            DESeq2-normalized counts, genes × samples.
        """
        # pydeseq2 expects samples × genes
        counts_t = counts.T.copy()
        # Ensure integer type
        counts_t = counts_t.astype(int)

        # Align metadata to same sample order
        meta_aligned = metadata.loc[counts_t.index].copy()

        # Keep only design factor columns (preserve original types — pydeseq2 handles both numeric and categorical)
        meta_sub = meta_aligned[self.design_factors].copy()

        inference = DefaultInference(n_cpus=1)
        dds = DeseqDataSet(
            counts=counts_t,
            metadata=meta_sub,
            design_factors=self.design_factors,
            refit_cooks=False,
            inference=inference,
        )
        dds.fit_size_factors()
        self._size_factors = dds.obs["size_factors"]

        # Retrieve normalized counts (samples × genes), transpose back
        norm_counts_t = dds.layers["normed_counts"]
        norm_df = pd.DataFrame(
            norm_counts_t,
            index=counts_t.index,
            columns=counts_t.columns,
        ).T  # back to genes × samples

        if save_path is not None:
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            norm_df.to_csv(save_path)
            print(f"Saved normalized counts to {save_path}")

        return norm_df

    @property
    def size_factors(self):
        return self._size_factors
