from setuptools import setup
import re, os, sys

setup(
    name="ann_nmf",
    version="0.0.1",
    python_requires=">3.6.0",
    author="Michael E. Vinyard - Harvard University - Massachussetts General Hospital - Broad Institute of MIT and Harvard",
    author_email="mvinyard@broadinstitute.org",
    url="https://github.com/mvinyard/ann_nmf",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    description="ann_nmf - AnnData wrapper of the ARD-NMF module from SignatureAnalyzer",
    packages=[
        "ann_nmf",
        "ann_nmf._supporting_functions",
        "ann_nmf._utilities",
    ],
    install_requires=[
        "anndata>=0.7.8",
        "matplotlib>=3.5.1",        
        "scanpy>=1.8.2",
        "signatureanalyzer>=0.0.7",
        "pandas>=1.3.5",
        "pydk>=0.0.4",
        "vinplots>=0.0.43",
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license="MIT",
)
