[![Build Status](https://travis-ci.org/mikessh/vdjtools.svg?branch=master)](https://travis-ci.org/mikessh/vdjtools)
[![Licence](https://img.shields.io/hexpm/l/plug.svg)](http://www.apache.org/licenses/LICENSE-2.0)
[![RepSeq](http://statsarray.com/wp-content/uploads/2014/03/omictools-logo.png)](http://omictools.com/rep-seq-c424-p1.html)

# VDJ tools

A comprehensive framework for post-analysis of immune repertoire sequencing data.
Compiled binaries are available from [here](https://github.com/mikessh/vdjtools/releases/latest).
The software is cross-platform and requires Java v1.7+ to run and R to generate high-quality graphics.

## Documentation

* See project wiki for detailed documentation: https://github.com/mikessh/vdjtools/wiki

* Examples: in `examples/` folder together with corresponding shell scripts.

* An example snippet to perform rarefaction, Variable segment usage analysis and repertoire clustering for RepSeq samples coming from donors of different ages

```bash
git clone https://github.com/mikessh/vdjtools-examples
VDJTOOLS="java -Xmx6G -jar vdjtools-1.0-SNAPSHOT.jar"
PARAMS="-S mitcr -m ./vdjtools-examples/aging_lite/samples/metadata.txt"
PLOT_PARAMS="-p -f age -n"
$VDJTOOLS RarefactionPlot $PARAMS $PLOT_PARAMS ./vdjtools/
$VDJTOOLS CalcSegmentUsage $PARAMS $PLOT_PARAMS ./vdjtools/
$VDJTOOLS BatchIntersectPair $PARAMS ./vdjtools/
$VDJTOOLS BatchIntersectPairPlot $PLOT_PARAMS ./vdjtools/
```

## Compiling from source

Clone the repository and compiling using [Maven](maven.apache.org)

```bash
git clone https://github.com/mikessh/vdjtools
cd vdjtools && mvn clean install
```

### Prerequisites

The following steps should be performed to compile VDJtools from sources:

* Make sure you are compiling under Java v1.7+

* [VDJdb](https://github.com/mikessh/vdjdb) dependency should be manually installed:

```bash
git clone https://github.com/mikessh/vdjdb.git
cd vdjdb
mvn clean install
```

* Install required R packages by running 
```bash
java -jar vdjtools-1.0-SNAPSHOT.jar RInstall
```