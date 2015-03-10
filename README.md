[![Build Status](https://travis-ci.org/mikessh/vdjtools.svg?branch=master)](https://travis-ci.org/mikessh/vdjtools)

# VDJ tools

A comprehensive framework for post-analysis of immune repertoire sequencing data.
Compiled binaries are available from [here](https://github.com/mikessh/vdjtools/releases/latest).
The software is cross-platform and requires Java v1.7+ to run.

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

## Documentation

* Wiki: https://github.com/mikessh/vdjtools/wiki

* Examples: in `examples/` folder together with corresponding shell scripts.

* Javadocs: http://mikessh.github.io/vdjtools-doc/ **(under development)**