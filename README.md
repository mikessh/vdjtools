[![Build Status](https://travis-ci.org/mikessh/vdjtools.svg?branch=master)](https://travis-ci.org/mikessh/vdjtools)

## VDJtools

A comprehensive analysis framework for T-cell and B-cell repertoire sequencing data.
Compiled binaries are available from [here](https://github.com/mikessh/vdjtools/releases/latest).
The software is cross-platform and requires Java v1.8 to run and R to perform plotting.

Easy installation on **MacOS/Linux** via [Homebrew](http://brew.sh/) or [Linuxbrew](http://linuxbrew.sh/):
```bash
brew tap mikessh/repseq
brew install vdjtools
vdjtools CalcBasicStats ...
```

For **Windows** just use the latest VDJtools bundle marked with ``.win.zip`` suffix.

List of features and detailed documentation can be found at [ReadTheDocs](http://vdjtools-doc.readthedocs.org/en/latest/).

Example datasets and shell scripts are provided in a separate [repository](https://github.com/mikessh/vdjtools-examples).

### Please cite VDJtools as:

Shugay M et al. VDJtools: Unifying Post-analysis of T Cell Receptor Repertoires. [PLoS Comp Biol 2015; 11(11):e1004503-e1004503](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004503).

### VDJtools was used in the following publications:

- Feng Y et al. A mechanism for expansion of regulatory T-cell repertoire and its role in self-tolerance. [Nature 2015; doi:10.1038/nature16141](http://www.nature.com/nature/journal/vaop/ncurrent/full/nature16141.html)
- Britanova OV et al. Dynamics of Individual T Cell Repertoires: From Cord Blood to Centenarians. [J Immunol 2016; doi:10.4049/​jimmunol.1600005](http://www.jimmunol.org/content/196/12/5005.short)
- Joachims ML et al. Single-cell analysis of glandular T cell receptors in Sjögren’s syndrome. [JCI Insight 2016; doi:10.1172/jci.insight.85609](https://insight.jci.org/articles/view/85609)
