# Installation

## From PyPI (planned)

```bash
pip install spatial_omics
```

## From source (now)

```bash
git clone https://github.com/chxhstatis/spatial_omics
cd spatial_omics
pip install -e ".[dev,docs]"
```

`spatial_omics` needs Python ≥ 3.9 and `numpy, pandas, scipy, scikit-learn, anndata,
matplotlib` (installed automatically).

## Verify

```bash
python examples/quickstart.py     # runs the whole pipeline on synthetic data
pytest                            # runs the test suite
```

## Optional

- `pip install "spatial_omics[viz]"` — adds `squidpy` for richer spatial plots.
- Upstream (FASTQ → matrix) needs `bowtie2, samtools, bedtools, bbmap` — see
  [From FASTQ](tutorials/from-fastq.md).
