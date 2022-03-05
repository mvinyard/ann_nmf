# ![ann_nmf.logo](docs/imgs/ann_nmf.logo.svg)

[![PyPI pyversions](https://img.shields.io/pypi/pyversions/ann_nmf.svg)](https://pypi.python.org/pypi/ann_nmf/)
[![PyPI version](https://badge.fury.io/py/ann_nmf.svg)](https://badge.fury.io/py/ann_nmf)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

AnnData wrapper of the [**ARD-NMF module**](https://github.com/broadinstitute/SignatureAnalyzer-GPU/blob/master/ARD_NMF.py) from [**SignatureAnalyzer-GPU**](https://github.com/broadinstitute/SignatureAnalyzer-GPU).

## Installation

Install using [pip](https://pypi.org/project/ann-nmf/)
```BASH
pip install ann_nmf
```

Alternatively, install the development version:
```BASH
git clone https://github.com/mvinyard/ann_nmf.git; cd ann_nmf;
pip install -e .
```

## API overview

#### Import libraries and get some data
```python
import ann_nmf
import scanpy as sc

adata = sc.datasets.pbmc3k()
ann_nmf.ut.preprocess_raw_counts(adata)
```

#### Key class:
```python
nmf = ann_nmf.NMF(adata, outdir="nmf_results/pbmc3k") # saves .h5 file
nmf.run(n_runs=10, K0=20, max_iter=2000)
```

#### SignatureAnalyzer visualization:
```python
nmf.cluster()
nmf.signatures()
```
<a><img src="https://user-images.githubusercontent.com/47393421/156896801-a5816aaf-3ecb-466c-8b8b-e92ff20f02bc.png"  height=325 hspace=20></a><a><img src="https://user-images.githubusercontent.com/47393421/156896796-3037dabb-000b-4593-8eb1-161c2d22c62f.png" height=325 hspace=20></a>

```python
nmf.plot_umap_signatures(adata)
```
![image](https://user-images.githubusercontent.com/47393421/156896803-d835c4e0-9fa5-460a-a194-cfc63af2fa77.png)

```python
ann_nmf.ut.fetch_signature(nmf.nmf_genes, 1)
```
<img width="1081" alt="gene-sig-table-example" src="https://user-images.githubusercontent.com/47393421/156896865-a2c0a29e-3b97-4050-ba60-2b75911fa4a1.png">


## Conceptual background and foundational work

* ARD-NMF theory [Arxiv](https://arxiv.org/pdf/1111.6085.pdf)
* SignatureAnalyzer ([GitHub](https://github.com/getzlab/getzlab-SignatureAnalyzer))
* SignatureAnalyzer-GPU ([GitHub](https://github.com/broadinstitute/SignatureAnalyzer-GPU))

## Acknowledgements

* Most of the code to wrap SignatureAnalyzer in an AnnData-friendly API was borrowed directly (and shamelessly) from **Shankara Anand** (**[@shankara-a](https://github.com/shankara-a)**) with only slight refactoring for more flexibility with fewer dependencies on install.
