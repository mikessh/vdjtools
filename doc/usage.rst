Usage
=====

Command line usage
------------------

General way to execute VDJtools
`routines <https://github.com/mikessh/vdjtools/wiki/Modules>`__ would be
the following,

::

    java -Xmx16G -jar vdjtools.jar RoutineName [arguments] -S software -m metadata.txt output/prefix

Here ``-S software`` argument specifies the input
`format <https://github.com/mikessh/vdjtools/wiki/Formats>`__ type.

The ``-m metadata.txt`` argument specifies
`metadata <https://github.com/mikessh/vdjtools/wiki/Metadata>`__ file
containing relative sample paths, sample names and any other information
that could be used later.

Alternatively, ``-m`` argument could be substituted with a
space-separated list of files, e.g.

::

    java -Xmx16G -jar vdjtools.jar RoutineName -S software sample1.txt[.gz] sample2.txt[.gz] ... output/prefix

Whether not explicitly used (such as in *BatchIntersectPairPlot*
routine) and applicable, plotting is turned on with ``-p`` argument.

The ``-h`` argument will bring up help message for specified routine.

Output prefix could be either an output directory name (if ended with
``/``) or an output file prefix, which would be appended with an
intuitive suffix and extension.

**Note on memory usage:** > Some routines could be memory demanding,
especially when running sample intersection routines with many large
(~1,000,000 clonotypes) datasets. In this case Java Virtual Machine
(JVM) could drop with ``java.lang.OutOfMemoryError: Java heap space``
exception. The memory limit for JVM is set manually using ``-Xmx?G``
parameter. Setting this parameter to 20-60Gb of memory should be enough
for most purposes, e.g. 100 samples with 500,000 clonotypes on average.

--------------

Usage examples
--------------

There are several data bundles and shell scripts that cover most of
VDJtools usage scenarios. All of the examples contain a folder with
clonotype abundance tables (``samples/``) and a
`metadata <https://github.com/mikessh/vdjtools/wiki/Metadata>`__ file,
as well as a shell script ``run.sh`` that contains a line-by-line
execution of various VDJtools routines. The user should modify the shell
script to ensure that ``$VDJTOOLS`` variable points to executable JAR
file.

--------------

Aging
~~~~~

[[/images/age-logo.jpg]]

The aging experiment involving 39 healthy donors of various ages and
both genders (see this
`paper <http://www.jimmunol.org/cgi/pmidlookup?view=long&pmid=24510963>`__
for details). This example allows to have a look at how a diverse set of
repertoire characteristics changes as we age.

The full dataset could be downloaded using this
`link <>`__\ [**TODO**\ ], while the lite version is provided in
VDJtools
`repository <https://github.com/mikessh/vdjtools/tree/master/examples>`__
for the sake of convenience. One can run the full set of analysis
routines on this dataset by first setting the path to VDJtools
executable jar and sample metadata:

.. code:: bash

    # assuming you compiled from source
    VDJTOOLS="java -Xmx6G -jar ../../target/vdjtools-1.0-SNAPSHOT.jar"
    # the data was processed using MiTCR
    PARAMS="-S mitcr -m samples/metadata.txt"

and then running routines sequentially

.. code:: bash

    # Basic
    $VDJTOOLS CalcBasicStats $PARAMS ./out/0
    $VDJTOOLS CalcSpectratype $PARAMS ./out/1
    # -p for plotting, -f specifies metadata column for coloring, 
    # -n tells that factor is continuous
    $VDJTOOLS CalcSegmentUsage $PARAMS -p -f age -n ./out/2
    # the following routines run on a single sample
    $VDJTOOLS PlotFancySpectratype -S mitcr ./samples/A4-i125.txt.gz ./out/3
    $VDJTOOLS PlotSpectratypeV -S mitcr ./samples/A4-i125.txt.gz ./out/4
    $VDJTOOLS PlotFancyVJUsage -S mitcr ./samples/A4-i125.txt.gz ./out/5

    # Diversity
    $VDJTOOLS PlotQuantileStats -S mitcr ./samples/A4-i125.txt.gz ./out/6
    # Compute the resampling-based diversity estimates by selecting half
    # of the reads (all samples contain 10000 reads)
    $VDJTOOLS CalcDiversityStats $PARAMS -x 5000 ./out/7
    # -l specifies metadata column used as label
    $VDJTOOLS RarefactionPlot $PARAMS -f age -n -l sample.id ./out/8

    # Intersect
    $VDJTOOLS IntersectPair -S mitcr -p ./samples/A4-i189.txt.gz ./samples/A4-i190.txt.gz ./out/9
    # computes various metrics characterizing divergence between repertoires
    $VDJTOOLS BatchIntersectPair $PARAMS ./out/10
    # plotting routine is separated from time-consuming batch intersection
    # sample clustering is performed on this stage.
    # Here we use relative sample overlap as metric and age as continuous factor
    $VDJTOOLS BatchIntersectPairPlot -f age -n -l sample.id ./out/10 ./out/10.age
    # here we use Variable segment Jensen-Shannon divergence and sex as discrete factor
    $VDJTOOLS BatchIntersectPairPlot -m vJSD -f sex -l sample.id ./out/10 ./out/10.sex

    # Annotation
    # you can use flexible filter for scanning dataset that accepts regexp syntax, 
    # '__' marks the column in annotation database aka VDJdb
    $VDJTOOLS ScanDatabase $PARAMS -f --filter "__origin__=~/EBV/" ./out/11

Below is an example of ``CalcSegmentUsage`` graphical output:

[[/images/age-vusage.png]]

--------------

Hematopoietic stem cell transfer (HSCT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[[/images/hsct-logo.jpg]]

HSCT is a great model for clonotype tracking and studying how the
diversity of immune repertoire restores following myeloablation. An
example dataset with four time points is available in VDJtools
`repository <https://github.com/mikessh/vdjtools/tree/master/examples>`__.

As the samples were built using various TCR analysis software, they were
all formatted in the same way to preserve minimal necessary information
(``simple`` format). Therefore parameters should be set as following:

.. code:: bash

    PARAMS="-S simple -m samples/metadata.txt"

The following analysis routines are

.. code:: bash

    # Basic
    $VDJTOOLS CalcBasicStats $PARAMS ./out/0
    $VDJTOOLS CalcSpectratype $PARAMS ./out/1
    $VDJTOOLS CalcSegmentUsage $PARAMS -p -f "Time post HSCT, months" -n ./out/2

    # Diversity
    # Note that selecting the factor having spaces in its name requires using double quotes
    $VDJTOOLS CalcDiversityStats $PARAMS ./out/3
    $VDJTOOLS RarefactionPlot $PARAMS -f "Time post HSCT, months" -n -l sample.id ./out/4

    # Intersect
    # this routine by default detects clonotypes that are present in 2 or more samples
    # and builds a time course for them, 
    # but here we trace clonotypes from first time point setting -x 0
    $VDJTOOLS IntersectPair -S simple -p ./samples/minus48months.txt.gz ./samples/4months.txt.gz ./out/5
    $VDJTOOLS IntersectSequential $PARAMS -f "Time post HSCT, months" -x 0 -p ./out/6

    # Annotation
    # can also use Groovy/Java syntax in filter
    $VDJTOOLS ScanDatabase $PARAMS -f --filter \
    "__origin__.contains('CMV')||__origin__.contains('EBV')" ./out/7

Rarefaction plot shows how repertoire diversity is lost and restored
during post-HSCT period. The output of ``ScanDatabase`` displays that
CMV- and EBV-specific clonotypes start to dominate in the repertoire:
they comprise ~4% of repertoire prior to HSCT, but increase more than
2-fold in post-HSCT period. Stackplot showing time course for the
abundance of top 100 clonotypes is displayed below:

[[/images/hsct-stackplot.png]]

--------------

Multiple sclerosis (MS)
~~~~~~~~~~~~~~~~~~~~~~~

[[/images/ms-logo.jpg]]

A usage example involving MS study will appear here when the VDJtools
paper is published :)
