# VDJ tools

A comprehensive framework for post-analysis of immune repertoire sequencing data.
Compiled binaries are available from [here](https://github.com/mikessh/vdjtools/releases/latest).

## Installing

Install by cloning the repository and compiling using [Maven](maven.apache.org)
```bash
git clone https://github.com/mikessh/vdjtools
cd vdjtools && mvn clean install
```

## Prerequisites

The following steps should be performed to compile VDJtools from sources:

* Make sure you are compiling under Java v1.7

* [VDJdb](https://github.com/mikessh/vdjdb) and [MiLib](https://github.com/milaboratory/milib) dependency should be manually installed:

```bash
git clone --branch 1.0 --depth 1 https://github.com/milaboratory/milib.git
cd milib && mvn clean install && cd ..
git clone https://github.com/mikessh/vdjdb.git
cd vdjdb
mvn clean install
```

* Install required R packages by running 
```
java -jar vdjtools-1.0-SNAPSHOT.jar RInstall
```

## Documentation

Please see the [wiki page](https://github.com/mikessh/vdjtools/wiki) for more information.

The examples are provided in `examples/` together with corresponding shell scripts.