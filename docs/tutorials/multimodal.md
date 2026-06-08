# Multi-modal: intrinsic vs extrinsic variance

The flagship analysis (paper Fig 4): given matched **DNA** (subclones) and **RNA**
(cell types) on the same tissue, decompose each gene's spatial expression variance
into **genetic** (which clone a cell belongs to — intrinsic) vs **microenvironment**
(local tumour / immune density — extrinsic).

This is what makes spatial multi-modal data more than two side-by-side panels: a
single regression attributes expression to *what a cell is* vs *where it sits*.

```python
import spatial_omics as so

# 1. DNA side — clones (as in Getting started)
dna = so.datasets.simulate(coverage=0.4, n_clones=2)
so.pp.bin_qc(dna); so.pp.correct_bias(dna); so.pp.normalize(dna)
so.tl.dual_smooth(dna); so.tl.call_clones(dna)

# 2. RNA side — matched spatial expression with cell types
#    (synthetic here; on real data load your RNA AnnData with obs['cell_type'])
rna = so.datasets.simulate_rna(dna, n_genes=200)

# 3. Decompose
so.tl.variance_decomposition(dna, rna,
                             clone_key="clone", celltype_key="cell_type",
                             tumor_types=("Tumor",), immune_types=("Immune",))

rna.uns["spatial_omics_vardecomp"]
# {'mean_frac_genetic': ..., 'mean_frac_microenv': ..., 'n_genes_genetic_dominant': ...}

so.pl.variance_decomposition(rna)   # stacked genetic/microenv/residual per gene
```

Each gene gets `var['vd_genetic']`, `var['vd_microenv']`, `var['vd_residual']` and a
`var['vd_dominant']` call. For the synthetic data the recovered split matches the
ground-truth `var['true_driver']`.

!!! note "On real data"
    Provide an RNA `AnnData` aligned to the same DBiT grid (`obs['x_id']`,
    `obs['y_id']`, `obsm['spatial']`) with a cell-type annotation in
    `obs['cell_type']` (e.g. from RCTD deconvolution). Spots are matched to the DNA
    assay by grid coordinate.
