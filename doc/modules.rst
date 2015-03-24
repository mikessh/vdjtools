Modules
=======

General
-------

VDJtools software package contains a comprehensive set of immune
repertoire post-analysis routines, which are subdivided into several
analysis modules. Each module's section provides command line usage
syntax and parameter descriptions for each of the routines, as well as
output example and description.

`Basic <https://github.com/mikessh/vdjtools/wiki/Modules-basic>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Summary statistics, spectratyping, etc

-  `CalcBasicStats <https://github.com/mikessh/vdjtools/wiki/Modules-basic#calcbasicstats>`__
   Computes summary statistics for samples: read counts, mean clonotype
   sizes, number of non-functional clonotypes, etc
-  `CalcSegmentUsage <https://github.com/mikessh/vdjtools/wiki/Modules-basic#calcsegmentusage>`__
   Computes Variable (V) and Joining (J) segment usage vectors
-  `CalcSpectratype <https://github.com/mikessh/vdjtools/wiki/Modules-basic#calcspectratype>`__
   Computes spectratype, the distribution of clonotype abundance by CDR3
   sequence length
-  `PlotFancySpectratype <https://github.com/mikessh/vdjtools/wiki/Modules-basic#plotfancyspectratype>`__
   Plots spectratype explicitly showing top N clonotypes
-  `PlotFancyVJUsage <https://github.com/mikessh/vdjtools/wiki/Modules-basic#plotfancyvjusage>`__
   Plots the frequency of different V-J pairings
-  `PlotSpectratypeV <https://github.com/mikessh/vdjtools/wiki/Modules-basic#plotspectratypev>`__
   Plots distribution of V segment abundance by resulting CDR3 sequence
   length

`Diversity <https://github.com/mikessh/vdjtools/wiki/Modules-diversity>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Repertoire richness and diversity

-  `PlotQuantileStats <https://github.com/mikessh/vdjtools/wiki/Modules-diversity#plotquantilestats>`__
   Visualizes repertoire clonality
-  `RarefactionPlot <https://github.com/mikessh/vdjtools/wiki/Modules-diversity#rarefactionplot>`__
   Plots rarefaction analysis results
-  `CalcDiversityStats <https://github.com/mikessh/vdjtools/wiki/Modules-diversity#calcdiversitystats>`__
   Computes repertoire diversity estimates

`Intersection <https://github.com/mikessh/vdjtools/wiki/Modules-intersection>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clonotype overlap between samples

-  `IntersectPair <https://github.com/mikessh/vdjtools/wiki/Modules-intersection#intersectpair>`__
   Computes intersection between a pair of samples
-  `BatchIntersectPair <https://github.com/mikessh/vdjtools/wiki/Modules-intersection#batchintersectpair>`__
   Computes pairwise intersections for a list of samples
-  `BatchIntersectPairPlot <https://github.com/mikessh/vdjtools/wiki/Modules-intersection#batchintersectpairplot>`__
   Plots results of batch intersection
-  `IntersectSequential <https://github.com/mikessh/vdjtools/wiki/Modules-intersection#intersectsequential>`__
   Intersects an ordered list of samples
-  `PoolSamples <https://github.com/mikessh/vdjtools/wiki/Modules-intersection#poolsamples>`__
   Pools several samples together

`Manipulation <https://github.com/mikessh/vdjtools/wiki/Modules-manipulation>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Filtering and resampling

-  `FilterNonFunctional <https://github.com/mikessh/vdjtools/wiki/Modules-manipulation#filternonfunctional>`__
   Filters non-functional clonotypes
-  `DownSample <https://github.com/mikessh/vdjtools/wiki/Modules-manipulation#downsample>`__
   Performs down-sampling
-  `ApplySampleAsFilter <https://github.com/mikessh/vdjtools/wiki/Modules-manipulation#applysampleasfilter>`__
   Given a list of samples filters clonotypes that are present in a
   specified sample
-  `Decontaminate <https://github.com/mikessh/vdjtools/wiki/Modules-manipulation#decontaminate>`__
   Filters possible cross-sample contamination

`Annotation <https://github.com/mikessh/vdjtools/wiki/Modules-annotation>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clonotype annotation

-  `ScanDatabase <https://github.com/mikessh/vdjtools/wiki/Module-annotations#scandatabase>`__
   Queries a database containing clonotype of known antigen specificity

--------------

Common parameters
-----------------

There are several parameters that are commonly used among analysis
routines:

+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                                                                                                                                                                                                                       |
+=============+========================+============+===================================================================================================================================================================================================================================================================================================================================================+
| ``-S``      | ``--software``         | string     | Software used to process RepSeq data. Currently supported: ``mitcr``,\ ``migec``,\ ``igblast`` and ``simple``. See `Formats <https://github.com/mikessh/vdjtools/wiki/Input#formats>`__ section                                                                                                                                                   |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path       | Path to metadata file. Should point to a tab-delimited file with the first two columns containing sample path and sample id respectively, and the remaining columns containing user-specified data. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section                                                            |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | `Intersection type <https://github.com/mikessh/vdjtools/wiki/Modules#intersection-type>`__, that specifies which clonotype features (CDR3 sequence, V/J segments, hypermutations) will be compared when checking if two clonotypes match. Allowed values: ``strict``,\ ``nt``,\ ``ntV``,\ ``ntVJ``,\ ``aa``,\ ``aaV``,\ ``aaVJ`` and ``aa!nt``.   |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Brings up the help message for selected routine                                                                                                                                                                                                                                                                                                   |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-p``      | ``--plot``             |            | (*plotting*) Enable plotting for routines that supports it. See `this <https://github.com/mikessh/vdjtools/wiki/Installation#plotting-routines>`__ section.                                                                                                                                                                                       |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--factor``           | string     | (*plotting*) Name of the sample metadata column that should be treated as factor. If the name contains spaces, the argument should be surrounded with double quotes, e.g. ``-f "Treatment type"``                                                                                                                                                 |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-n``      | ``--factor-numeric``   |            | (*plotting*) Treat the factor as numeric?                                                                                                                                                                                                                                                                                                         |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-l``      | ``--label``            | string     | (*plotting*) Name of the sample metadata column that should be treated as label. If the name contains spaces, the argument should be surrounded with double quotes, e.g. ``-l "Patient id"``                                                                                                                                                      |
+-------------+------------------------+------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Intersection type
-----------------

Some of VDJtools routines require to define clonotype matching strategy.
For example a common situation is when one is interested in calculating
overlap between samples or estimating the extent of convergent
recombination, which is the number of distinct nucleotide CDR3 sequences
per one CDR3 amino acid sequence. The list of strategies is given below.

+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   | Rule                                        | Note                                                                                                                                  |
+=============+=============================================+=======================================================================================================================================+
| strict      | *CDR3nt*\ (AND)*V*\ (AND)*J*\ (AND)*SHMs*   | Require full match for receptor nucleotide sequence                                                                                   |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| nt          | *CDR3nt*                                    |                                                                                                                                       |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| ntV         | *CDR3nt*\ (AND)*V*                          |                                                                                                                                       |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| ntVJ        | *CDR3nt*\ (AND)*V*\ (AND)*J*                |                                                                                                                                       |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| aa          | *CDR3aa*                                    |                                                                                                                                       |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| aaV         | *CDR3aa*\ (AND)*V*                          |                                                                                                                                       |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| aaVJ        | *CDR3aa*\ (AND)*V*\ (AND)*J*                |                                                                                                                                       |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+
| aa!nt       | *CDR3aa*\ (AND)((NOT)*CDR3nt*)              | Removes nearly all contamination bias from overlap results. Should not be used for samples from the same donor/tracking experiments   |
+-------------+---------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------+

Basic
-----

CalcBasicStats
~~~~~~~~~~~~~~

This routine computes a set of basic sample statistics, such as read
counts, number of clonotypes, etc.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar CalcBasicStats \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                             |
+=============+=======================+============+=========================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__            |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                    |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

The following table with ``.basicstats.txt`` suffix is generated,

+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column                 | Description                                                                                                                                                                                          |
+========================+======================================================================================================================================================================================================+
| sample\_id             | Sample unique identifier                                                                                                                                                                             |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ...                    | Metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section                                                                                                 |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| count                  | Number of reads in a given sample                                                                                                                                                                    |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| diversity              | Number of clonotypes in a given sample                                                                                                                                                               |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| mean\_frequency        | Mean clonotype frequency                                                                                                                                                                             |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| geomean\_frequency     | Geometric mean of clonotype frequency                                                                                                                                                                |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| nc\_diversity          | Number of non-coding clonotypes                                                                                                                                                                      |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| nc\_frequency          | Frequency of reads that belong to non-coding clonotypes                                                                                                                                              |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| mean\_cdr3nt\_length   | Mean length of CDR3 nucleotide sequence. Weighted by clonotype frequency                                                                                                                             |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| mean\_insert\_size     | Mean number of inserted random nucleotides in CDR3 sequence. Characterizes V-J insert for receptor chains without D segment, or a sum of V-D and D-J insert sizes. Weighted by clonotype frequency   |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| mean\_ndn\_size        | Mean number of nucleotides that lie between V and J segment sequences in CDR3. Weighted by clonotype frequency                                                                                       |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| convergence            | Mean number of unique CDR3 nucleotide sequences per a single CDR3 amino acid sequence they are translated to                                                                                         |
+------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Graphical output**

none

--------------

CalcSegmentUsage
~~~~~~~~~~~~~~~~

This routine computes Variable (V) and Joining (J) segment usage
vectors, i.e. the frequency of associated reads for each of V/J segments
present in sample(s). If plotting is on, will also perform clustering
for V/J usage vectors and samples *à la* gene expression analysis.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar CalcSegmentUsage \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                                               |
+=============+=======================+============+===========================================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                              |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                     |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-u``      | ``--unweighted``      |            | Instead of computing read frequency, will compute the number of unique clonotypes with specific /J segments                               |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-p``      | ``--plot``            |            | Turns on plotting. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                         |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--factor``          | string     | Specifies plotting factor. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                 |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-n``      | ``--numeric``         |            | Specifies if plotting factor is numeric. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-l``      | ``--label``           | string     | Specifies label used for plotting. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__         |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                                      |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

The following tables with
``.segments.[unwt or wt depending on -u parameter].[V or J].txt`` suffix
are generated,

+-----------------------------------+--------------------------------------------------------------------------------------------------------+
| Column                            | Description                                                                                            |
+===================================+========================================================================================================+
| sample\_id                        | Sample unique identifier                                                                               |
+-----------------------------------+--------------------------------------------------------------------------------------------------------+
| ...                               | Metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section   |
+-----------------------------------+--------------------------------------------------------------------------------------------------------+
| Segment name, e.g. TRBJ1-1        | Segment frequency in a given sample                                                                    |
+-----------------------------------+--------------------------------------------------------------------------------------------------------+
| Next segment name, e.g. TRBJ1-2   | ...                                                                                                    |
+-----------------------------------+--------------------------------------------------------------------------------------------------------+
| ...                               | ...                                                                                                    |
+-----------------------------------+--------------------------------------------------------------------------------------------------------+

**Graphical output**

Images, having the same name as tables, with the exception of ``.pdf``
extension, are created if plotting is on. They display segment usage
heatmap and hierarchical clustering for samples and segment.

[[/images/modules/basic-segmentusage.png]]

--------------

CalcSpectratype
~~~~~~~~~~~~~~~

Calculates
`spectratype <http://www.jimmunol.org/content/152/10/5109.full.pdf+html>`__,
that is, histogram of read counts by CDR3 nucleotide length. The
spectratype is useful to detect pathological and highly clonal
repertoires, as the spectratype of non-expanded T- and B-cells has a
symmetric gaussian-like distribution.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar CalcSpectratype \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                             |
+=============+=======================+============+=========================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__            |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-u``      | ``--unweighted``      |            | Instead of computing read frequency, will compute the number of unique clonotypes with specific CDR3 length             |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-a``      | ``--amino-acid``      |            | Will use CDR3 amino acid sequences for calculation instead of nucleotide ones                                           |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                    |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

The following table with
``.spectratype.[aa or nt  depending on -a parameter].[unwt or wt depending on -u parameter].txt``
suffix is generated,

+------------------------+--------------------------------------------------------------------------------------------------------+
| Column                 | Description                                                                                            |
+========================+========================================================================================================+
| sample\_id             | Sample unique identifier                                                                               |
+------------------------+--------------------------------------------------------------------------------------------------------+
| ...                    | Metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section   |
+------------------------+--------------------------------------------------------------------------------------------------------+
| CDR3 length, e.g. 22   | Frequency of reads with a given CDR3 length in a given sample                                          |
+------------------------+--------------------------------------------------------------------------------------------------------+
| Next CDR3 length, 23   | ...                                                                                                    |
+------------------------+--------------------------------------------------------------------------------------------------------+
| ...                    | ...                                                                                                    |
+------------------------+--------------------------------------------------------------------------------------------------------+

**Graphical output**

none

--------------

PlotFancySpectratype
~~~~~~~~~~~~~~~~~~~~

Plots a spectratype that also displays CDR3 lengths for top N clonotypes
in a given sample. This plot allows to detect the highly-expanded
clonotypes.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar PlotFancySpectratype [options] sample.txt output_prefix

**Parameters**

+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                    |
+=============+=======================+============+================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-t``      | ``--top``             | int        | Number of top clonotypes to visualize. Should not exceed 20, default is 10                                     |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                           |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+

**Tabular output**

Following table with ``.fancyspectra.txt`` prefix is generated,

+-------------------------------------+----------------------------------------------------------------------+
| Column                              | Description                                                          |
+=====================================+======================================================================+
| Len                                 | Length of CDR3 nucleotide sequence                                   |
+-------------------------------------+----------------------------------------------------------------------+
| Other                               | Frequency of clonotypes with a given CDR3 length, other than top N   |
+-------------------------------------+----------------------------------------------------------------------+
| Clonotype#N, e.g. CASRLLRAGSTEAFF   | Clonotype frequency, at the corresponding CDR3 length                |
+-------------------------------------+----------------------------------------------------------------------+
| Clonotype#N-1                       | ...                                                                  |
+-------------------------------------+----------------------------------------------------------------------+
| ...                                 | ...                                                                  |
+-------------------------------------+----------------------------------------------------------------------+

**Graphical output**

The following image file with ``.fancyspectra.pdf`` suffix,

[[/images/modules/basic-fancyspectra.png]]

--------------

PlotFancyVJUsage
~~~~~~~~~~~~~~~~

Plots a `circos <http://circos.ca/>`__-style V-J usage plot displaying
the frequency of various V-J junctions.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar PlotFancyVJUsage [options] sample.txt output_prefix

**Parameters**

+-------------+-----------------------+------------+-----------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                     |
+=============+=======================+============+=================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__    |
+-------------+-----------------------+------------+-----------------------------------------------------------------------------------------------------------------+
| ``-u``      | ``--unweighted``      |            | Instead of computing read frequency, will compute the number of unique clonotypes with specific V-J junctions   |
+-------------+-----------------------+------------+-----------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                            |
+-------------+-----------------------+------------+-----------------------------------------------------------------------------------------------------------------+

**Tabular output**

A matrix with rows corresponding to different J segments and columns
corresponding to different V segments. Each cells contains the frequency
of a give V-J junction. The file has
``.fancyvj.[unwt or wt depending on -u parameter].txt`` suffix.

**Graphical output**

An image having the same name as the output table, with the exception of
``.pdf`` extension, is generated. Arcs correspond to different V and J
segments, scaled to their frequency in sample. Ribbons represent V-J
pairings and their size is scaled to the pairing frequency.

[[/images/modules/basic-fancyvj.png]]

--------------

PlotSpectratypeV
~~~~~~~~~~~~~~~~

Plots a detailed spectratype containing additional info displays CDR3
length distribution for clonotypes from top N Variable segment families.
This plot is useful to detect type 1 and type 2 repertoire
`biases <http://www.nature.com/nri/journal/v6/n12/fig_tab/nri1977_T1.html>`__,
that could arise under pathological conditions.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar PlotSpectratypeV [options] sample.txt output_prefix

**Parameters**

+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                    |
+=============+=======================+============+================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-t``      | ``--top``             | int        | Number of top (by frequency) V segments to visualize. Should not exceed 12 default is 12                       |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-u``      | ``--unweighted``      |            | Instead of counting read frequency, will count the number of unique clonotypes                                 |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                           |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+

**Tabular output**

Following table with
``.spectraV.[unwt or wt depending on -u parameter].txt`` prefix is
generated,

+----------------------------+--------------------------------------------------------------------------------------------+
| Column                     | Description                                                                                |
+============================+============================================================================================+
| Len                        | Length of CDR3 nucleotide sequence                                                         |
+----------------------------+--------------------------------------------------------------------------------------------+
| Other                      | Frequency of clonotypes with a given CDR3 length, having V segments other than the top N   |
+----------------------------+--------------------------------------------------------------------------------------------+
| Segment#N, e.g. TRBV10-1   | Frequency of clonotypes with a given V segment at the corresponding CDR3 length            |
+----------------------------+--------------------------------------------------------------------------------------------+
| Segment#N-1                | ...                                                                                        |
+----------------------------+--------------------------------------------------------------------------------------------+
| ...                        | ...                                                                                        |
+----------------------------+--------------------------------------------------------------------------------------------+

**Graphical output**

The following image file with
``.spectraV.[unwt or wt depending on -u parameter].pdf`` suffix,

[[/images/modules/basic-spectrav.png]]

Diversity
---------

PlotQuantileStats
~~~~~~~~~~~~~~~~~

Plots a three-layer donut chart to visualize the repertoire clonality.
\* First layer ("set") includes the frequency of singleton ("1", met
once), doubleton ("2", met twice) and high-order ("3+", met three or
more times) clonotypes. Singleton and doubleton frequency is an
important factor in estimating the total repertoire diversity, e.g.
Chao1 diversity estimator (see `Colwell *et
al* <http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf>`__).
We have also recently
`shown <http://www.ncbi.nlm.nih.gov/pubmed/24510963>`__ that in whole
blood samples, singletons have very nice correlation with the number of
naive T-cells, which are the backbone of immune repertoire diversity. \*
The second layer ("quantile"), displays the abundance of top 20% ("Q1"),
next 20% ("Q2"), ... (up to "Q5") clonotypes for clonotypes from "3+"
set. In our experience this quantile plot is a simple and efficient way
to display repertoire clonality. \* The last layer ("top") displays the
individual abundances of top N clonotypes.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar PlotQuantileStats [options] sample.txt output_prefix

**Parameters**

+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                    |
+=============+=======================+============+================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-t``      | ``--top``             | int        | Number of top clonotypes to visualize. Should not exceed 10, default is 5                                      |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                           |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------------------+

**Tabular output**

Following table with ``.qstat.txt`` prefix is generated,

+----------+--------------------------------------------------------+
| Column   | Description                                            |
+==========+========================================================+
| Type     | Detalization level: ``set``, ``quantile`` or ``top``   |
+----------+--------------------------------------------------------+
| Name     | Variable name: "1", "Q1", "CASSLAPGATNEKLFF", etc      |
+----------+--------------------------------------------------------+
| Value    | Corresponding relative abundance                       |
+----------+--------------------------------------------------------+

**Graphical output**

Following plot with ``.qstat.pdf`` prefix is generated,

[[/images/modules/diversity-qstat.png]]

--------------

RarefactionPlot
~~~~~~~~~~~~~~~

Plots rarefaction curves for specified list of samples, that is, the
dependencies between sample diversity and sample size. Those curves are
interpolated from 0 to the current sample size and then extrapolated up
to the size of the largest of samples, allowing comparison of diversity
estimates. Interpolation and extrapolation are based on multinomial
models, see `Colwell *et
al* <http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf>`__
for details.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar RarefactionPlot \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                                                                                                        |
+=============+========================+============+====================================================================================================================================================================================================================================+
| ``-S``      | ``--software``         | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                       |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                              |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Set the intersection type used to collapse clonotypes before computing diversity. Defaults to ``strict`` (don't collapse at all). See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-s``      | ``--steps``            | integer    | Set the total number of points in the rarefaction curve, default is ``101``                                                                                                                                                        |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--factor``           | string     | Specifies plotting factor. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                          |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-n``      | ``--numeric``          |            | Specifies if plotting factor is numeric. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                            |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-l``      | ``--label``            | string     | Specifies label used for plotting. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                  |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--wide-plot``        |            | Set wide plotting area                                                                                                                                                                                                             |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                                                                                                                                                               |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

The following table with
``rarefaction.[intersection type shorthand].txt`` is generated:

+--------------+---------------------------------------------------------------------------------------------------------------+
| Column       | Definition                                                                                                    |
+==============+===============================================================================================================+
| sample\_id   | Sample unique identifier                                                                                      |
+--------------+---------------------------------------------------------------------------------------------------------------+
| ...          | Sample metadata columns, see `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section   |
+--------------+---------------------------------------------------------------------------------------------------------------+
| x            | Subsample size, reads                                                                                         |
+--------------+---------------------------------------------------------------------------------------------------------------+
| mean         | Mean diversity at given size                                                                                  |
+--------------+---------------------------------------------------------------------------------------------------------------+
| ciL          | Lower bound of 95% confidence interval                                                                        |
+--------------+---------------------------------------------------------------------------------------------------------------+
| ciU          | Upper bound of 95% confidence interval                                                                        |
+--------------+---------------------------------------------------------------------------------------------------------------+
| type         | Data point type: ``0=interpolation``, ``1=exact``, ``2=extrapolation``                                        |
+--------------+---------------------------------------------------------------------------------------------------------------+

**Graphical output**

A figure with the same suffix as output table and ``.pdf`` extension is
provided. Solid and dashed lines mark interpolated and extrapolated
regions of rarefaction curves respectively, points mark exact sample
size and diversity. Shaded areas mark 95% confidence intervals.

[[/images/modules/diversity-rarefaction.png]]

--------------

CalcDiversityStats
~~~~~~~~~~~~~~~~~~

Computes a set of diversity statistics, including \* Observed diversity
\*
`Chao <http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf>`__
and `Efron-Thisted <www.jstor.org/stable/2335721>`__ lower bound total
diversity (LBTD) estimates \*
`Shannon-Weaver <http://www.esajournals.org/doi/abs/10.2307/1934352>`__
and `Inverse
Simpson <http://www.esajournals.org/doi/abs/10.2307/1934352>`__
diversity indices \* `Extrapolated Chao diversity
estimate <http://viceroy.eeb.uconn.edu/estimates/EstimateSPages/EstSUsersGuide/References/ColwellEtAl2012.pdf>`__
(``chaoE``).

Diversity stats are computed in two modes: using original data and via
several re-sampling steps (usually down-sampling to the size of smallest
dataset).

-  The estimates computed on original data could be biased by uneven
   sampling depth (sample size), of those only ``chaoE`` is properly
   normalized to be compared between samples. While not good for
   between-sample comparison, the LBTD estimates provided for original
   data are most useful for studying the fundamental properties of
   repertoires under study, i.e. to answer the question how large the
   repertoire diversity of an entire organism could be.
-  Estimates computed using re-sampling are useful for between-sample
   comparison, e.g. we have successfully used the re-sampled
   (normalized) observed diversity to measure the repertoire aging
   trends (see `this <http://www.ncbi.nlm.nih.gov/pubmed/24510963>`__
   paper).

In our recent experience the ``chaoE`` estimate and LBTD estimates
compared on re-sampled data provide best results for between-sample
comparisons.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar CalcDiversityStats \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                                                                                                        |
+=============+========================+============+====================================================================================================================================================================================================================================+
| ``-S``      | ``--software``         | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                       |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                              |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Set the intersection type used to collapse clonotypes before computing diversity. Defaults to ``strict`` (don't collapse at all). See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-x``      | ``--downsample-to``    | integer    | Set the sample size to interpolate the diversity estimates via resampling. Default = size of smallest sample. Applies to diversity estimates stored in ``.resampled.txt`` table                                                    |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-X``      | ``--extrapolate-to``   | integer    | Set the sample size to extrapolate the diversity estimates. Default = size of largest sample. Currently, only applies to ``chaoE`` diversity estimate.                                                                             |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                                                                                                                                                               |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Two tables with ``diversity.[intersection type shorthand].txt`` and
``diversity.[intersection type shorthand].resampled.txt`` are generated,
containing diversity estimates computed on original and down-sampled
datasets respectively.

Note that ``chaoE`` estimate is only present in the table generated for
original samples. Both tables contain means and standard deviations of
diversity estimates. Also note that standard deviation and mean values
for down-sampled datasets are computed based on N=3 re-samples.

Here is an example column layout, similar between both output tables

+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+
| Column                                                                                                                                                    | Definition                                                                                                    |
+===========================================================================================================================================================+===============================================================================================================+
| sample\_id                                                                                                                                                | Sample unique identifier                                                                                      |
+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+
| ...                                                                                                                                                       | Sample metadata columns, see `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section   |
+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+
| reads                                                                                                                                                     | Number of reads in the sample                                                                                 |
+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+
| diversity                                                                                                                                                 | Diversity of the original sample (after collapsing to unique clonotypes according to ``-i`` parameter)        |
+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+
| extrapolate\_reads / resample\_reads                                                                                                                      | The reads used to extrapolate or re-sample in order to compute present diversity estiamtes                    |
+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+
| <*name*\ >\_mean\|Mean value of the diversity estimate <*name*\ > \|<*name*\ >\_std\|Standard deviation of the diversity estimate <*name*\ > \|...\|...   |                                                                                                               |
+-----------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------+

**Graphical output**

none

Intersection
------------

IntersectPair
~~~~~~~~~~~~~

Intersects clonotype lists from a pair of samples.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar IntersectPair [options] sample1.txt sample2.txt output_prefix

**Parameters**

+-------------+------------------------+------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                         |
+=============+========================+============+=====================================================================================================================================================+
| ``-S``      | ``--software``         | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                        |
+-------------+------------------------+------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Sample intersection rule. Defaults to ``strict``. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__    |
+-------------+------------------------+------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-t``      | ``--top``              | int        | Number of top clonotypes to visualize explicitly on stack are plot and provide in the collapsed joint table. Should not exceed 100, default is 20   |
+-------------+------------------------+------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-p``      | ``--plot``             |            | Turns on plotting. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                   |
+-------------+------------------------+------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                                                                                |
+-------------+------------------------+------------+-----------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Two joint clonotype abundance tables with
``paired.[intersection type shorthand].table.txt`` and
``paired.[intersection type shorthand].table.collapsed.txt`` suffices
are generated. The latter one is collapsed up to top N clonotypes. See
**tabular output** in
`IntersectSequential <https://github.com/mikessh/vdjtools/wiki/Modules#intersectsequential>`__
section for detailed description of table fields.

A summary table (``paired.[intersection type shorthand].summary.txt``
suffix) containing information on sample overlap size, etc, is also
provided. See tabular output in
`BatchIntersectPair <https://github.com/mikessh/vdjtools/wiki/Modules#batchintersectpair>`__
section for details.

**Graphical output**

A composite plot having
``paired.[intersection type shorthand].scatter.pdf`` suffix is
generated. It contains a scatterplot of clonotype abundances that show
overlapping clonotypes and a linear regression. Point size is scaled to
clonotype abundance. The plot also contains two marginal histograms each
composed of the overlapping (red) and total clonotype (grey) abundance
distributions in corresponding sample. Histograms are weighted by
clonotype abundance, i.e. they show read distribution by clonotype size.

[[/images/modules/intersect-pair-scatter.png]]

The second plot file with
``.paired.[intersection type shorthand].table.collapsed.pdf`` suffix
contains a clonotype abundance stack area plot. It shows details for top
N clonotypes, as well as collapsed ("NotShown") and non-overlapping
("NonOverlapping") clonotypes. Clonotype CDR3 amino acid sequence is
plotted against the sample where the clonotype reaches maximum
abundance.

[[/images/modules/intersect-pair-stack.png]]

--------------

BatchIntersectPair
~~~~~~~~~~~~~~~~~~

Performs clonotype list intersections between all possible pairs of
provided samples. At least 3 samples should be provided. Note that this
is one of most memory-demanding routines, as it will load all samples
into memory at once (unless used with ``--low-mem`` option). While this
tool provides only a tabular output, it should be used together with
`BatchIntersectPairPlot <https://github.com/mikessh/vdjtools/wiki/Modules#batchintersectpairplot>`__
which performs cluster analysis and visualization.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar BatchIntersectPair \
    [options] [sample1.txt sample2.txt sample3.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                    |
+=============+========================+============+================================================================================================================================================+
| ``-S``      | ``--software``         | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                   |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                          |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Sample intersection rule. Defaults to ``aa``. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--low-mem``          |            | Low memory mode, will keep only a pair of samples in memory during execution, but run much slower.                                             |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                                                                           |
+-------------+------------------------+------------+------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

A table suffixed
``intersect.batch.[intersection type shorthand].summary.txt`` with a
comprehensive information on sample pair intersections is generated.
This table is non-redundant: it contains ``N * (N - 1) / 2`` rows
corresponding to upper diagonal of matrix of possible pairs ``(i,j)``.
Table layout is given below.

\|General info \|-----

+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| Column          | Description                                                                                                                 |
+=================+=============================================================================================================================+
| 1\_sample\_id   | First sample unique identifier                                                                                              |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| 2\_sample\_id   | Second sample unique identifier                                                                                             |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| div1            | Total number of clonotypes in the first sample after identical clonotypes are collapsed based on intersection type ``-i``   |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| div2            | Same as above, second sample                                                                                                |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| div12           | Number of overlapping clonotypes                                                                                            |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| div21           | Same as above                                                                                                               |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| count1          | Total number of reads in the first sample                                                                                   |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| count2          | ...                                                                                                                         |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| count12         | For clonotypes **overlapping** between two samples: total number of reads they have in the **first** sample                 |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| count21         | ...                                                                                                                         |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| freq1           | Total clonotype relative abundance for the first sample (should be 1.0 if sample is unaltered)                              |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| freq2           | ...                                                                                                                         |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| freq12          | For clonotypes **overlapping** between two samples: their sum of relative abundances in the **first** sample                |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+
| freq21          | ...                                                                                                                         |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------+

\|Overlap metrics \|-----

+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column   | Description                                                                                                                                                       |
+==========+===================================================================================================================================================================+
| R        | Correlation between relative abundances of overlapping clonotypes in first and second samples                                                                     |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| D        | The diversity of overlap, equals to ``div12 / (div1 * div2)``                                                                                                     |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| F        | Relative abundance of overlapping clonotypes, equals to ``sqrt(freq1*freq2)``                                                                                     |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| F2       | Given fi1 and fi2 are the frequency of i'th overlapping clonotype in sample 1 and 2 respectively, equals to ``sum(sqrt(fi1 * fi2), i=1..div12)``                  |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| vJSD     | `Jensen-Shannon divergence <https://www.cise.ufl.edu/~anand/sp06/jensen-shannon.pdf>`__ of Variable segment usage distribution between two samples                |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| vjJSD    | Same as above for the concatenated Variable and Joining segment usage distributions, i.e. histograms with ``#V unique segments + #J unique segments`` bins each   |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| vj2JSD   | Same as above for the Variable-Joining junction usage matrix, i.e. histograms with ``#V unique segments * #J unique segments`` bins each                          |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sJSD     | Jensen-Shannon divergence of spectratypes                                                                                                                         |
+----------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+

\|Sample metadata \|-----

+----------+---------------------------------------------------------------------------------------------------------------------+
| Column   | Description                                                                                                         |
+==========+=====================================================================================================================+
| 1\_...   | First sample metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section   |
+----------+---------------------------------------------------------------------------------------------------------------------+
| 2\_...   | Second sample metadata columns                                                                                      |
+----------+---------------------------------------------------------------------------------------------------------------------+

**Graphical output**

none

--------------

BatchIntersectPairPlot
~~~~~~~~~~~~~~~~~~~~~~

This routine provides cluster analysis and plotting for
BatchIntersectPair output. Note that this routine requires 1) setting
input file prefix same as the output prefix of BatchIntersectPair and 2)
setting the same ``-i`` argument.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar BatchIntersectPairPlot \
    [options] batch_intersect_pair_output_prefix [output_prefix]

**Parameters**

+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                                                                                                                      |
+=============+========================+============+==================================================================================================================================================================================================================================================+
| ``-m``      | ``--measure``          | string     | Specifies which sample overlap metric to use. Defaults to ``F2``. Allowed values: ``R``,\ ``D``,\ ``F``,\ ``F2``,\ ``vJSD``,\ ``vjJSD``,\ ``vj2JSD`` and ``sJSD``. See tabular output of BatchIntersectPair for details.                         |
+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Intersection type, should be the same as used in BatchIntersectPair. Defaults to ``aa``. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                          |
+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--factor``           | string     | Specifies plotting factor. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                        |
+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-n``      | ``--numeric``          |            | Specifies if plotting factor is numeric. Also determines how the post-hoc tests for relation between factor value and clustering are performed. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-l``      | ``--label``            | string     | Specifies label used for plotting. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                |
+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                                                                                                                                                                             |
+-------------+------------------------+------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Three auxillary tables are generated: \* Table suffixed
``mds.coords.[value of -i argument].[value of -m argument].txt``
contains coordinates of samples computed using multi-dimensional scaling
(MDS) \* In case ``-n`` is not specified a table suffixed
``perms.[value of -i argument].[value of -m argument].txt`` contains
results of within- and between-cluster distances generated using sample
label permutation \* In case ``-k`` is specified, a table suffixed
``hc.cut$k.[value of -i argument].[value of -m argument].txt`` is
generated. It contains cluster labels for samples obtained using
dendrogram cutting

**Graphical output**

Hierarchical clustering output is stored in a file suffixed
``hc.[value of -i argument].[value of -m argument].pdf``. Clustering is
performed using ``hcl`` util in R with default parameters. Distances are
scaled as ``-log10(.)`` and ``(1-.)/2`` for relative overlap size and
correlation respectively; in case of Jensen-Shannon divergence no
scaling is performed. Node colors correspond to factor value.

[[/images/modules/intersect-batch-dendro.png]]

Multi-dimensional scaling is performed using ``isoMDS`` function from
``MASS`` R package with number of dimensions set as ``k=2``. The file is
suffixed
``mds.coords.[value of -i argument].[value of -m argument].pdf``.

[[/images/modules/intersect-batch-mds.png]]

A plot showing the significance of sample distances within- and
between-groups is generated in case the factor is non-numeric (``-n``).
It contains a histogram of distances obtained using permutations with
red vertical line indicating the observed distance and P-value. The file
is suffixed ``perms.[value of -i argument].[value of -m argument].pdf``.

[[/images/modules/intersect-batch-perms.png]]

--------------

IntersectSequential
~~~~~~~~~~~~~~~~~~~

This routine performs an all-vs-all intersection between an ordered list
of samples for clonotype tracking purposes. Users can specify clonotypes
from which sample to trace, e.g. the pre-therapy sample. Alternatively,
the output will contain all clonotypes present in at lease 2+ samples.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar IntersectSequential \
    [options] [sample1.txt sample2.txt sample3.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument          | Description                                                                                                                                                                                                                                                                                                                                        |
+=============+========================+===================+====================================================================================================================================================================================================================================================================================================================================================+
| ``-S``      | ``--software``         | string            | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                                                                                                                                       |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path              | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                                                                                                                              |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string            | Sample intersection rule. Defaults to ``strict``. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                                                                                                   |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--factor``           | string            | Specifies factor that should be treated as time variable. Factor values should be numeric. Defaults to 'time'. If such column is not present in metadata, time points are taken either from values provided with ``-s`` argument or sample order. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-x``      | ``--track-sample``     | integer           | A zero-based index of time point to track. If not provided, will consider all clonotypes that were detected in 2+ samples                                                                                                                                                                                                                          |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-s``      | ``--sequence``         | ``[t1,t2,...]``   | Time point sequence. Unused if -m is specified. If not specified, either time values from metadata, or sample indexes (as in command line) are used.                                                                                                                                                                                               |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-t``      | ``--top``              | int               | Number of top clonotypes to visualize explicitly on stack are plot and provide in the collapsed joint table. Should not exceed 100, default is 200                                                                                                                                                                                                 |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-p``      | ``--plot``             |                   | Turns on plotting. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                                                                                                                                                  |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |                   | Display help message                                                                                                                                                                                                                                                                                                                               |
+-------------+------------------------+-------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Summary table suffixed ``sequential.[value of -i argument].summary.txt``
is created with the following columns.

+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column          | Description                                                                                                                                                                                                                                                                                               |
+=================+===========================================================================================================================================================================================================================================================================================================+
| 1\_sample\_id   | First sample unique identifier                                                                                                                                                                                                                                                                            |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2\_sample\_id   | Second sample unique identifier                                                                                                                                                                                                                                                                           |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| value           | Value of the intersection metric                                                                                                                                                                                                                                                                          |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| metric          | Metric type: ``diversity``, ``frequency`` or ``count``. Metrics correspond to the number of unique clonotypes, total frequency and total read count for clonotypes overlapping between first and second sample. In case tracking is on (``-x``), only clonotypes present in tracked sample are counted.   |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1\_time         | Time value for the first sample                                                                                                                                                                                                                                                                           |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2\_time         | Time value for the second sample                                                                                                                                                                                                                                                                          |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 1\_...          | First sample metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section                                                                                                                                                                                         |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| 2\_...          | Second sample metadata columns                                                                                                                                                                                                                                                                            |
+-----------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Two joint clonotype abundance tables with
``sequential.[intersection type shorthand].table.txt`` and
``sequential.[intersection type shorthand].table.collapsed.txt``
suffices are generated. The latter one is collapsed up to top N
clonotypes. Those tables contain the following columns.

    NOTE: When several clonotype variants are present in samples that
    correspond to the same clonotype under ``-i`` conditions (e.g.
    several Variable segment variants when ``-i nt`` is set), only the
    most frequent form is taket to final output.

+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column          | Description                                                                                                                                                  |
+=================+==============================================================================================================================================================+
| count           | Clonotype count, normalized so that clonotypes with smallest frequency have count of ``1``                                                                   |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| freq            | Clonotype frequency, computed as geometric mean of clonotype frequencies in intersected samples. If clonotype is missing, its frequency is set to ``1e-9``   |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| cdr3nt          | CDR3 nucleotide sequence, see `Input <https://github.com/mikessh/vdjtools/wiki/Input>`__ section                                                             |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| cdr3aa          | CDR3 amino acid sequence                                                                                                                                     |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| v               | Variable segment                                                                                                                                             |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| d               | Diversity segment                                                                                                                                            |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| j               | Joining segment                                                                                                                                              |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| peak            | Index of a time point at which given clonotype reaches its maximum frequency                                                                                 |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| <sample name>   | Frequency of a given clonotype at corresponding sample                                                                                                       |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ...             |                                                                                                                                                              |
+-----------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Graphical output**

Summary table is visualized in a plot file suffixed
``sequential.[value of -i argument].summary.pdf``.

[[/images/modules/intersect-seq-summary.png]]

A plot file with
``.sequential.[intersection type shorthand].stackplot.pdf`` suffix
contains a clonotype abundance stack area plot. It shows details for top
N clonotypes, as well as collapsed ("NotShown") and non-overlapping
("NonOverlapping") clonotypes. Clonotype CDR3 amino acid sequence is
plotted against the sample where the clonotype reaches maximum
abundance. Clonotypes are colored by the peak position of their
abundance profile.

[[/images/modules/intersect-seq-stackplot.png]]

Clonotype abundance for top N clonotypes is also visualized using
heatmap (``.sequential.[intersection type shorthand].heatplot.pdf``). It
also includes a dendrogram showing the clustering of clonotype abundance
profiles. suffix contains a clonotype abundance stack area plot.
Clonotypes that are missing in a given sample are shown with grey.

[[/images/modules/intersect-seq-heatplot.png]]

--------------

PoolSamples
~~~~~~~~~~~

<*UNDER DEVELOPMENT*\ >

Sample manipulation
------------

FilterNonFunctional
~~~~~~~~~~~~~~~~~~~

Filters non-functional (non-coding) clonotypes, i.e. the ones that
contain a stop codon or frameshift in their receptor sequence. Those
clonotypes do not have any functional role, but they are useful for
dissecting and studying the V-(D)-J recombination machinery as they do
not pass thymic selection.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar FilterNonFunctional \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                             |
+=============+=======================+============+=========================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__            |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-e``      | ``--negative``        |            | Negative filtering, i.e. only non-functional clonotypes are retained                                                    |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                                                                                            |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                    |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``ncfilter:[retain or remove based on -e option]`` to ``..filter..``
metadata column.

Creates a filter summary file with a ``ncfilter.summary.txt`` suffix
containing info on the number of unique clonotypes that passed the
filtering process, their total frequency and count.

**Graphical output**

none

--------------

Downsample
~~~~~~~~~~

Down-samples a list of clonotype abundance tables by randomly selecting
a pre-defined number of reads. This routine could be useful for a)
normalizing samples for further highly-sensitive comparison b) speeding
up computation / decreasing file size.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar Downsample \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                             |
+=============+=======================+============+=========================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__            |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-x``      | ``--num-reads``       | integer    | Number of reads to take. **Required**                                                                                   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                                                                                            |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                    |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``ds:[-x value]`` to ``..filter..`` metadata column.

**Graphical output**

none

--------------

ApplySampleAsFilter
~~~~~~~~~~~~~~~~~~~

Retains/filters out all clonotypes found in a given sample *S* from
other samples. Useful when *S* contains some specific cells of interest
e.g. tumor-infiltrating T-cells or sorted tetramer+ T-cells.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar ApplySampleAsFilter \
    [options] [sample1.txt sample2.txt ... if -m is not specified] filter_sample output_prefix

**Parameters**

+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                                                                        |
+=============+========================+============+====================================================================================================================================================+
| ``-S``      | ``--software``         | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                       |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                              |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Sample intersection rule. Defaults to ``strict``. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__   |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-e``      | ``--negative``         |            | Negative filtering, i.e. only clonotypes absent in sample *S* are retained                                                                         |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``         |            | Compress output sample files                                                                                                                       |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                                                                               |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``asaf:[- if -e, + otherwise]:[-i value]`` to ``..filter..`` metadata
column.

**Graphical output**

none

--------------

Decontaminate
~~~~~~~~~~~~~

DNA contamination can occur at library prep stage, for example sample
barcode swithing resulting from PCR chimeras. Those could lead to a high
number of artificial shared clonotypes for samples sequenced in the same
batch. If no sophisticated library prep method (e.g. paired-end
barcoding) is applied, it is highly recommended to filter those before
performing any kind of cross-sample analysis.

This routine filters out all clonotypes that have a matching clonotype
in a different sample which is ``-r`` times more abundant.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar ApplySampleAsFilter \
    [options] [sample1.txt sample2.txt ... if -m is not specified] filter_sample output_prefix

**Parameters**

+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                              |
+=============+=======================+============+==========================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__             |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__    |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-r``      | ``--ratio``           | numeric    | Parent-to-child clonotype frequency ratio for contamination filtering. Defaults to ``20``                                |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
|             | ``--low-mem``         |            | Will process all sample pairs sequentially, avoiding loading all of them into memory. Slower but memory-efficient mode   |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                                                                                             |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                     |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``dec:[-r value]`` to ``..filter..`` metadata column.

**Graphical output**

none

Sample annotation
-----------------

ScanDatabase
~~~~~~~~~~~~

Annotates a set of samples using immune receptor database based on
V-(D)-J junction matching. By default uses
`VDJdb <https://github.com/mikessh/vdjdb>`__, which contains CDR3
sequences, Variable and Joining segments of known specificity obtained
using literature mining. This routine supports user-provided databases
and allows flexible filtering of results based on database fields. The
output of ScanDatabase includes both detailed (clonotype-wise)
annotation of samples and summary statistics. Only amino-acid CDR3
sequences are used in database querying.

**Command line usage**

::

    java -Xmx4G -jar vdjtools.jar ScanDatabase \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

**Parameters**

+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument         | Description                                                                                                                                                                                                                |
+=============+=======================+==================+============================================================================================================================================================================================================================+
| ``-S``      | ``--software``        | string           | Input format. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                               |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path             | Path to metadata file. See `Common parameters <https://github.com/mikessh/vdjtools/wiki/Modules#common-parameters>`__                                                                                                      |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-D``      | ``--database``        | path             | Path to an external database file. Will use built-in VDJdb if not specified.                                                                                                                                               |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-d``      | ``--details``         |                  | Will provide a detailed output for each sample with annotated clonotype matches                                                                                                                                            |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--fuzzy``           |                  | Will query database allowing at most 2 substitutions, 1 deletion and 1 insertion but no more than 2 mismatches simultaneously. If not set, only exact matches will be reported                                             |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--filter``          | ``expression``   | Logical filter on database columns. Supports Regex, .contains(), .startsWith(), etc. Database field names should be surrounded with ``__`` and will be substituted during execution. For example ``"__origin__=~/EBV/"``   |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--v-match``         |                  | V segment must to match                                                                                                                                                                                                    |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--j-match``         |                  | J segment must to match                                                                                                                                                                                                    |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |                  | Display help message                                                                                                                                                                                                       |
+-------------+-----------------------+------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

**Tabular output**

A summary table suffixed ``annot.[database name].summary.txt`` is
generated. First header line marked with ``##FILTER`` contains filtering
expression that was used. The table contains the following columns:

+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column                           | Description                                                                                                                                                                                                                                                                                      |
+==================================+==================================================================================================================================================================================================================================================================================================+
| sample\_id                       | Sample unique identifier                                                                                                                                                                                                                                                                         |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ...                              | Sample metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section                                                                                                                                                                                      |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| diversity                        | Number of clonotypes in sample                                                                                                                                                                                                                                                                   |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| match\_size                      | Number of matches between sample and database. In case ``--fuzzy`` mode is on, all matches will be counted. E.g. if clonotype ``a`` in the sample matches clonotypes ``A`` and ``B`` in the database and clonotype ``b`` in the sample matches clonotype B the value in this column will be 3.   |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sample\_diversity\_in\_matches   | Number of unique clonotypes in the sample that matched clonotypes from the database                                                                                                                                                                                                              |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| db\_diversity\_in\_matches       | Number of unique clonotypes in the database that matched clonotypes from the sample                                                                                                                                                                                                              |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| sample\_freq\_in\_matches        | Overall frequency of unique clonotypes in the sample that matched clonotypes from the database                                                                                                                                                                                                   |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| mean\_matched\_clone\_size       | Geometric mean of frequency of unique clonotypes in the sample that matched clonotypes from the database                                                                                                                                                                                         |
+----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Detailed database query results will be also reported for each sample if
``-d`` is specified. Those tables are suffixed
``annot.[database name].[sample id].txt`` and contain the following
columns.

+-------------------+-----------------------------------------------------------------------+
| Column            | Description                                                           |
+===================+=======================================================================+
| score             | CDR3 sequence alignment score                                         |
+-------------------+-----------------------------------------------------------------------+
| query\_cdr3aa     | Query CDR3 amino acid sequence                                        |
+-------------------+-----------------------------------------------------------------------+
| query\_v          | Query Variable segment                                                |
+-------------------+-----------------------------------------------------------------------+
| query\_j          | Query Joining segment                                                 |
+-------------------+-----------------------------------------------------------------------+
| subject\_cdr3aa   | Subject CDR3 amino acid sequence                                      |
+-------------------+-----------------------------------------------------------------------+
| subject\_v        | Subject Variable segment                                              |
+-------------------+-----------------------------------------------------------------------+
| subject\_j        | Subject Joining segment                                               |
+-------------------+-----------------------------------------------------------------------+
| v\_match          | ``true`` if Variable segments of query and subject clonotypes match   |
+-------------------+-----------------------------------------------------------------------+
| j\_match          | ``true`` if Joining segments of query and subject clonotypes match    |
+-------------------+-----------------------------------------------------------------------+
| mismatches        | Comma-separated list of query->subject mismatches                     |
+-------------------+-----------------------------------------------------------------------+
| ...               | Database fields corresponding to subject clonotype                    |
+-------------------+-----------------------------------------------------------------------+

**Graphical output**

none
