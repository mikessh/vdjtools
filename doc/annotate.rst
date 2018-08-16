.. _annotate:

Annotation
----------

.. _SegmentsToFamilies:

SegmentsToFamilies
^^^^^^^^^^^^^^^^^^

Will replace V and J segment IDs in samples with segment 'family' IDs. Here, 'families' are defined as clusters of V/J 
sequences built using hierarchial clustering of pairwise amino acid sequence alignment distances. Thus, two segments are 
assigned to the same family if they have homologous sequence. The actual table of segment <> family conversions can 
be accessed `here <https://github.com/mikessh/vdjtools/blob/master/src/main/resources/vj_families.txt>`__.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS SegmentsToFamilies \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+--------------------+----------------------------------------------------+
| Shorthand   |      Long name        | Argument           | Description                                        |
+=============+=======================+====================+====================================================+
| ``-m``      | ``--metadata``        | path               | Path to metadata file. See :ref:`common_params`    |
+-------------+-----------------------+--------------------+----------------------------------------------------+
| ``-s``      | ``--species``         | name               | [Required] Species name: ``human`` or ``mouse``.   |
+-------------+-----------------------+--------------------+----------------------------------------------------+
| ``-h``      | ``--help``            |                    | Display help message                               |
+-------------+-----------------------+--------------------+----------------------------------------------------+
| ``-c``      |                       |                    | Compressed output.                                 |
+-------------+-----------------------+--------------------+----------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Samples are returned as is, with the content of ``v`` and ``j`` columns replaced by families.

A metadata file will be created for resulting samples with ``segm2fam`` 
appended to the ``..filter..`` metadata column.


Graphical output
~~~~~~~~~~~~~~~~

none


--------------

.. _CalcDegreeStats:

CalcDegreeStats
^^^^^^^^^^^^^^^

Performs a TCR neighborhood enrichment test (TCRNET), testing each sample for clonotypes 
that have more neighbours (higher **degree** in a graph), i.e. clonotypes with similar CDR3 amino acid sequences, than would be expected 
by chance according to some control dataset. User can specify the actual **search scope** (i.e. 
number of allowed CDR3 mismatches), whether to only compare clonotypes with same V/J, and the 
control sample. If control sample is not provided, a pooling (see :ref:`PoolSamples`) of all provided samples is used. 
Note that this test, if supplied with real samples and a control pooled using ``-i strict`` option 
will account for the number of neighbours with the same CDR3 amino acid sequence, but distinct nucleotide 
sequences. If this is not desired, all input samples and control should be pre-pooled with ``-i aa`` or 
``-i aaVJ`` to collapse variants coding for the amino acid CDR3 sequence.

.. note:: 
    
    Running this routine will not return the actual clonotype graph for you, just annotate input samples. 
    To build the graph, one should refer to `VDJmatch <https://github.com/antigenomics/vdjmatch>`__ software 
    and its ``Cluster`` routine. Make sure the search scope option is the same as ``-o`` used for ``CalcDegreeStats`` 
    and that all scoring/filtering is turned off. Next, one should retain only the edges that connect pairs of 
    enriched clonotypes and enriched clonotypes with their neighbours.


Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS CalcDegreeStats \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument           | Description                                                                                                                                                |
+=============+=======================+====================+============================================================================================================================================================+
| ``-m``      | ``--metadata``        | path               | Path to metadata file. See :ref:`common_params`                                                                                                            |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-b``      | ``--background``      | path               | Path to the background (control) sample, used to compute expected statistics/P-values. If not provided, will pool input samples and uses them as control.  |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-o``      | ``--search-scope``    | s,i,d              | Search scope: number of substitutions (s), indels (id) and total number of mismatches (t) allowed. Default is ``1,0,1``                                    |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-g``      | ``--grouping``        | string             | Primary grouping type, limits set of clonotype comparisons: 'dummy' (no grouping, default), 'vj' (same V and J) or 'vjl' (same V, J and CDR3 length).      |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-g2``     | ``--grouping2``       | string             | Secondary grouping, used for computing statistics, accepts same values as ``-g``. By default will select 'vjl' if no indels allowed and 'vj' otherwise.    |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |                    | Display help message                                                                                                                                       |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. note:: 
    
    There are two possible schemes for running the algorithm. Firstly, one can select, 
    say a search scope of ``1,0,1`` allowing no indels, and ``-g vjl`` to only allow comparisons
    between clonotypes that match in V, J and CDR3 length. Then, one should 
    only consider ``p.value.g`` in the output and disregard all columns with ``g2/group2``.
    On the other hand, if one wants to allow comparison of clonotypes with different V/J, 
    and/or comparisons with indels, the option ``-g dummy`` should be used. If one thinks there 
    might be certain biases in V/J frequencies between control/background sample and input samples, 
    and one wants to control for them, he should select ``-g2 vj``, then observed degree values 
    will be provided as is (i.e. not limiting clonotype comparisons to a fixed V/J), 
    but the expected degree will be corrected to account for V/J usage difference 
    between input sample and control. One should only consider ``p.value.g2`` 
    in this case. See below for more explaination on output columns.

Tabular output
~~~~~~~~~~~~~~

Processed samples will have additional annotation columns appended to VDJtools clonotype 
table columns. These columns are the following:

+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Column          | Description                                                                                                                                                                                           |
+=================+=======================================================================================================================================================================================================+
| degree.s        | Degree (number of neighbours) of a given clonotype in sample. The degree is the number of unique clonotypes (incl. nucleotide variants) that match a given clonotype under specified search scope.    |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| group.count.s   | Number of unique clonotypes that match the group, defined by primary grouping (``-g``), of a given clonotype in sample, say have the same V and J.                                                    |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| group2.count.s  | Same as above, but the group is defined by secondary grouping ``-g2``.                                                                                                                                |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| degree.c        | Degree (number of neighbours) of a given clonotype in the control sample.                                                                                                                             |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| group.count.c   | Number of unique clonotypes in the control sample that match the group of given clonotype as defined by primary grouping (``-g``).                                                                    |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| group2.count.c  | Same as above, but the group is defined by secondary grouping ``-g2``.                                                                                                                                |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| p.value.g       | P-value for the neighbour (degree) enrichment of a given clonotype according to primary grouping. The P-value is computed as ``Pbinom(n=degree.s|p=degree.c/group.count.c, N=group.count.s)``.        |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| p.value.g2      | P-value for the neighbour (degree) enrichment of a given clonotype according to secondary grouping. The P-value is computed as ``Ppoisson(n=degree.s|lambda=group.count.s*degree.c/group.count.c)``.  |
+-----------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

A metadata file will be created for resulting samples with ``degstat`` 
appended to the ``..filter..`` metadata column.


Graphical output
~~~~~~~~~~~~~~~~

none


--------------


.. _CalcCdrAAProfile:

CalcCdrAAProfile
^^^^^^^^^^^^^^^^

Generates amino acid physical properties profile of CDR3. Amino acids are 
first grouped to corresponding CDR3 sub-regions and then binned by position 
within the sub-region. Amino acids in a given bin is scored according to 
its physical properties, sums of those scores and total number of amino acids
is reported for each sample/sub-region/bin/property combination.

For example under the **polarity** property amino acids are marked as polar (``1``) 
and non-polar (``0``) and the sum of these values is returned. When divided by 
the total number of amino acids one will get the fraction of polar amino acids 
in a given sample/sub-region. For **volume** the same operation will return the 
average volume of amino acids.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS CalcCdrAAProfile \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument           | Description                                                                                                                                                |
+=============+=======================+====================+============================================================================================================================================================+
| ``-m``      | ``--metadata``        | path               | Path to metadata file. See :ref:`common_params`                                                                                                            |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-w``      | ``--weighted``        |                    | If set, will weight amino acid property values by clonotype frequency.                                                                                     |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-n``      | ``--normalize``       |                    | If set, will normalize amino acid property values by dividing them by corresponding CDR3 sub-region size.                                                  |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-r``      | ``--region-list``     | region1,...        | List of CDR3 sub-regions to count statistics for, default is ``"CDR3-full,VJ-junc,V-germ,J-germ``                                                          |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-o``      | ``--property-list``   | property1,...      | List of amino acid physicochemical properties to use, see below for allowed value. Uses all amino acid properties from list below by default.              |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |                    | Display help message                                                                                                                                       |
+-------------+-----------------------+--------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------+

Supported CDR3 sub-regions:

+-----------------+--------------------------------------------------------------------------+
| Name            | Description                                                              |
+=================+==========================================================================+
| ``CDR3-full``   | Complete CDR3 region                                                     |
+-----------------+--------------------------------------------------------------------------+
| ``CDR3-center`` | Central 5 amino acids of CDR3                                            |
+-----------------+--------------------------------------------------------------------------+
| ``V-germ``      | Germline part of CDR3 region corresponding to Variable segment           |
+-----------------+--------------------------------------------------------------------------+
| ``D-germ``      | Germline part of CDR3 region corresponding to Diversity segment          |
+-----------------+--------------------------------------------------------------------------+
| ``J-germ``      | Germline part of CDR3 region corresponding to Joining segment            |
+-----------------+--------------------------------------------------------------------------+
| ``VD-junc``     | Variable-Diversity segment junction, applicable when D segment is mapped |
+-----------------+--------------------------------------------------------------------------+
| ``DJ-junc``     | Diversity-Joining segment junction, applicable when D segment is mapped  |
+-----------------+--------------------------------------------------------------------------+
| ``VJ-junc``     | Variable-Joining segment junction, including D segment if it is mapped   |
+-----------------+--------------------------------------------------------------------------+

Supported amino acid physical properties (see `full table <https://github.com/mikessh/vdjtools/blob/master/src/main/resources/profile/aa_property_table.txt>`__ for raw values):

+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| Name              | Description                                                                                                     | Reference                                                       |
+===================+=================================================================================================================+=================================================================+
| ``alpha``         | Preference to appear in alpha helices                                                                           | Stryer L et al. Biochemistry, 5th edition. ISBN 978-0716746843  |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``beta``          | Preference to appear in beta sheets                                                                             | Stryer L et al. Biochemistry, 5th edition. ISBN 978-0716746843  |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``turn``          | Preference to appear in turns                                                                                   | Stryer L et al. Biochemistry, 5th edition. ISBN 978-0716746843  |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``surface``       | Residues that have unchanged accessibility area when PPI partner is present                                     | `PMID:22559010 <http://www.ncbi.nlm.nih.gov/pubmed/22559010>`__ |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``rim``           | Residues that have changed accessibility area, but no atoms with zero accessibility in PPI interfaces           | `PMID:22559010 <http://www.ncbi.nlm.nih.gov/pubmed/22559010>`__ |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``core``          | Residues that have changed accessibility area and at least one atom with zero accessibility in PPI interfaces   | `PMID:22559010 <http://www.ncbi.nlm.nih.gov/pubmed/22559010>`__ |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``disorder``      | Intrinsic structural disorder-promoting, order-promoting and neutral amino acids                                | `PMID:11381529 <http://www.ncbi.nlm.nih.gov/pubmed/11381529>`__ |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``charge``        | Charged/non-charged amino acids                                                                                 | `Wikipedia <https://en.wikipedia.org/wiki/Amino_acid>`__        |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``pH``            | Amino acid pH level                                                                                             | `Wikipedia <https://en.wikipedia.org/wiki/Amino_acid>`__        |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``polarity``      | Polar/non-polar amino acids                                                                                     | `Wikipedia <https://en.wikipedia.org/wiki/Amino_acid>`__        |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``hydropathy``    | Amino acid hydropathy                                                                                           | `Wikipedia <https://en.wikipedia.org/wiki/Amino_acid>`__        |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``volume``        | Amino acid volume                                                                                               | `Wikipedia <https://en.wikipedia.org/wiki/Amino_acid>`__        |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``strength``      | Strongly-interacting amino acids / amino acids depleted by purifying selection in thymus                        | `PMID:18946038 <http://www.ncbi.nlm.nih.gov/pubmed/18946038>`__ |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``mjenergy``      | Mean value of MJ statistical potential for each amino acid, used to derive 'strength'                           | `PMID:8604144 <https://www.ncbi.nlm.nih.gov/pubmed/8604144>`__  |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
| ``kf1``..``kf10`` | Values of 10 Kidera factors summarizing physicochemical properties of amino acids                               | unpublished                                                     |
+-------------------+-----------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------+
    
Tabular output
~~~~~~~~~~~~~~

A summary table with averaged amino acid property values is generated, 
suffixed ``cdr3aa.profile.[wt or unwt based on -u].txt``. The table contains 
the following columns:

+---------------+---------------------------------------------------------------------------------------------------------------+
| Column        | Description                                                                                                   |
+===============+===============================================================================================================+
| sample\_id    | Sample unique identifier                                                                                      |
+---------------+---------------------------------------------------------------------------------------------------------------+
| ...           | Sample metadata columns. See `Metadata <https://github.com/mikessh/vdjtools/wiki/Input#metadata>`__ section   |
+---------------+---------------------------------------------------------------------------------------------------------------+
| region        | Current CDR3 sub-region, see above                                                                            |
+---------------+---------------------------------------------------------------------------------------------------------------+
| property      | Amino acid physical property name, see above                                                                  |
+---------------+---------------------------------------------------------------------------------------------------------------+
| mean          | Mean property value                                                                                           |
+---------------+---------------------------------------------------------------------------------------------------------------+

Graphical output
~~~~~~~~~~~~~~~~

none


--------------

.. _Annotate2:

Annotate
^^^^^^^^

This routine will compute a set of properties for each clonotype's CDR3 sequence and 
append them to resulting clonotype table. For example, number of added N-nucleotides 
and the sum of polar amino acids in CDR3. The main difference from :ref:`CalcCdrAAProfile` 
is that the former computes sample-level average while this routine performs calculation 
on clonotype level.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS Annotate \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+--------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument           | Description                                                                                                                                                                                                                                                               |
+=============+=======================+====================+===========================================================================================================================================================================================================================================================================+
| ``-m``      | ``--metadata``        | path               | Path to metadata file. See :ref:`common_params`                                                                                                                                                                                                                           |
+-------------+-----------------------+--------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-b``      | ``--base``            | param1,param2,...  | Comma-separated list of basic clonotype features to calculate and append to resulting clonotype tables. See below for allowed values. Default: ``cdr3Length,ndnSize,insertSize``                                                                                          |
+-------------+-----------------------+--------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-a``      | ``--aaprop``          | property1,...      | Comma-separated list of amino acid properties. Amino acid property value sum will be calculated for CDR3 sequence (blank annotations will be generated for non-coding clonotypes). See below for allowed values. Default: ``hydropathy,charge,polarity,strength,contact`` |
+-------------+-----------------------+--------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |                    | Display help message                                                                                                                                                                                                                                                      |
+-------------+-----------------------+--------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

List of basic annotation properties:

+----------------+--------------------------------------------------------------------------------------------------+
| Name           | Description                                                                                      |
+================+==================================================================================================+
| ``cdr3Length`` | Length of CDR3 region                                                                            |
+----------------+--------------------------------------------------------------------------------------------------+
| ``NDNSize``    | Number of nucleotides between last base of V germline and first base of J germline parts of CDR3 |
+----------------+--------------------------------------------------------------------------------------------------+
| ``insertSize`` | Number of added N-nucleotides                                                                    |
+----------------+--------------------------------------------------------------------------------------------------+
| ``VDIns``      | Number of added N-nucleotides in V-D junction or ``-1`` if D segment is undefined                |
+----------------+--------------------------------------------------------------------------------------------------+
| ``DJIns``      | Number of added N-nucleotides in D-J junction or ``-1`` if D segment is undefined                |
+----------------+--------------------------------------------------------------------------------------------------+

See :ref:`CalcCdrAAProfile` for the list of amino acid properties available for annotation. 
Sum of specified amino acid property values across all amino acids of CDR3 will be computed. 
It can be divided by ``cdr3Length / 3`` basic property value to get the average.
    
Tabular output
~~~~~~~~~~~~~~

Processed samples will have additional annotation columns appended to VDJtools clonotype 
table columns. Those columns will be prefixed with ``base.`` for basic CDR3 properties 
and ``aaprop.`` for CDR3 amino acid composition properties.

A metadata file will be created for resulting samples with ``annot:[-b value]:[-a value]`` 
appended to the ``..filter..`` metadata column.

Graphical output
~~~~~~~~~~~~~~~~

none

----------------

.. _ScanDatabase:

ScanDatabase (DEPRECATED since v1.0.5, use `VDJmatch <https://github.com/antigenomics/vdjmatch>`__)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotates a set of samples using immune receptor database based on
V-(D)-J junction matching. By default uses
`VDJdb <https://github.com/antigenomics/vdjdb-db>`__, which contains CDR3
sequences, Variable and Joining segments of known specificity obtained
using literature mining. This routine supports user-provided databases
and allows flexible filtering of results based on database fields. The
output of ScanDatabase includes both detailed (clonotype-wise)
annotation of samples and summary statistics. Only amino-acid CDR3
sequences are used in database querying.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS ScanDatabase \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument         | Description                                                                                                                                                                       |
+=============+=======================+==================+===================================================================================================================================================================================+
| ``-m``      | ``--metadata``        | path             | Path to metadata file. See :ref:`common_params`                                                                                                                                   |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-D``      | ``--database``        | path             | Path to an external database file. Will use built-in VDJdb if not specified.                                                                                                      |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-d``      | ``--details``         |                  | Will provide a detailed output for each sample with annotated clonotype matches                                                                                                   |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--fuzzy``           |                  | Will query database allowing at most 2 substitutions, 1 deletion and 1 insertion but no more than 2 mismatches simultaneously. If not set, only exact matches will be reported    |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--filter``          | ``expression``   | Logical pre-filter on database columns. See below                                                                                                                                 |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--v-match``         |                  | V segment must to match                                                                                                                                                           |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--j-match``         |                  | J segment must to match                                                                                                                                                           |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |                  | Display help message                                                                                                                                                              |
+-------------+-----------------------+------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. note:: 
    
    Database filter is a logical expression that contains
    reference to input table columns. Database column name references should 
    be surrounded with double underscores (``__``). Syntax supports Regex and 
    standard Java/Groovy functions such as ``.contains()``, ``.startsWith()``, 
    etc. Here are some examples:
    
    .. code-block:: groovy    
        
        __origin__=~/EBV/
        !(__origin__=~/CMV/)
        
    Note that the expression should be quoted: ``--filter "__origin__=~/HSV/"``

Tabular output
~~~~~~~~~~~~~~

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

Graphical output
~~~~~~~~~~~~~~~~

none