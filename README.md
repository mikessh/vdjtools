# VDJ tools

A comprehensive framework for post-analysis of immune repertoire sequencing data

## Installing

Install by cloning the repository and compiling using Maven 
```bash
git clone https://github.com/mikessh/vdjtools
cd vdjtools && mvn clean install
```

## Prerequisites

* VDJdb dependency should be manually installed:

```bash
git clone https://github.com/mikessh/vdjdb
cd vdjdb && mvn clean install
```

* Install required R packages by running 
```
java -jar vdjtools-1.0-SNAPSHOT.jar RInstall
```

## Documentation

Please see the [wiki page](https://github.com/mikessh/vdjtools/wiki) for more information.

The examples are provided in `examples/` together with corresponding shell scripts.