# Frozen five-gene score

Genes: `PLOD2`, `MSN`, `TPM4`, `PPFIBP1`, and `SEC31A`.

For each cohort:

1. Collapse duplicated gene symbols by their mean expression where required.
2. Z-standardize each gene across samples.
3. Calculate the equal-weight mean of the five standardized genes for each
   sample.
4. Z-standardize that mean across samples in the cohort.

The score is not calculated when any frozen member is unavailable; remaining
genes are not reweighted. Cox-model hazard ratios correspond to a one-standard-
deviation increase in the final cohort-standardized score.

The gene set and equal weights were frozen before post-freeze testing in
GSE57495 and GSE183795. No outcome-based retuning was performed in those
cohorts.
