.. _preprocess:

Pre-processing
--------------

.. note::

    Most of routines specified in this section will output processed clonotype tables and re-normalize 
    individual clonotype frequencies by dividing their read count by the total read count in resulting (filtered/processed) sample.
    For some of the routines this behavior can be disabled with ``--save-freqs`` option. In this case original clonotype frequencies 
    will be carried over from input samples and they will likely not sum to ``1.0`` in the resulting clonotype table.

.. _Correct:

Correct
^^^^^^^

Performs frequency-based correction to eliminate erroneous clonotypes. Searches the sample for 
clonotype pairs that differ by one, two ... (up to specified depth) mismatches. In case 
the ratio of smallest to largest clonotype sizes is lower than the threshold specified 
as ``ratio ^ number_of_mismatches`` correction is performed. Largest clonotype in pair 
increases its size by the read count of the smaller one and the smaller 
one is discarded. Note that the original sample is not changed during correction, so 
all comparisons are performed with original count values and erroneous clonotypes are only 
removed after search procedure is finished. It is also possible to restrict correction to 
clonotypes with identical V/J segments using ``-a`` option.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS Correct \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+
| Shorthand | Long name           | Argument | Description                                                                                                         |
+===========+=====================+==========+=====================================================================================================================+
| ``-m``    | ``--metadata``      | path     | Path to metadata file. See :ref:`common_params`                                                                     |
+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+
| ``-d``    | ``--depth``         | 1+       | Maximum number of mismatches allowed between clonotypes being compared. Default is 2                                |
+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+
| ``-r``    | ``--ratio``         | [0, 1)   | Child-to-parent clonotype size ratio threshold under which child clonotype is considered erroneous. Default is 0.05 |
+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+
| ``-a``    | ``--match-segment`` |          | Check for erroneous clonotypes only among those that have identical V and J assignments                             |
+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+
| ``-c``    | ``--compress``      |          | Compress output sample files                                                                                        |
+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+
| ``-h``    | ``--help``          |          | Display help message                                                                                                |
+-----------+---------------------+----------+---------------------------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs corrected samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``corr:[-d option value]:[-r option value]:['vjmatch' or 'all' based on -a option]`` to 
``..filter..`` metadata column.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _Decontaminate:

Decontaminate
^^^^^^^^^^^^^

Cross-sample contamination can occur at library prep stage, for example sample
barcode swithing resulting from PCR chimeras. Those could lead to a high
number of artificial shared clonotypes for samples sequenced in the same
batch. If no sophisticated library prep method (e.g. paired-end
barcoding) is applied, it is highly recommended to filter those before
performing any kind of cross-sample analysis.

This routine filters out all clonotypes that have a matching clonotype
in a different sample which is ``-r`` times more abundant. Clonotype fractions 
within samples are considered, which is good for dealing with FACS-related contaminations.
In case of dealing with cross-sample contaminations in samples coming from the same 
sequencing lane use ``--read-based`` option that tells the routine to compare read counts.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS Decontaminate \
    [options] [sample1.txt sample2.txt ... if -m is not specified] filter_sample output_prefix

Parameters
~~~~~~~~~~

+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                              |
+=============+=======================+============+==========================================================================================================================+
| ``-S``      | ``--software``        | string     | Input format. See :ref:`common_params`                                                                                   |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
|             | ``--read-based``      | string     | If set will compare clonotype read counts. Clonotype fractions in corresponding samples are compared by default.         |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See :ref:`common_params`                                                                          |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-r``      | ``--ratio``           | numeric    | Parent-to-child clonotype frequency ratio for contamination filtering. Defaults to ``20``                                |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                                                                                             |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                                     |
+-------------+-----------------------+------------+--------------------------------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``dec:[-r value]`` to ``..filter..`` metadata column.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _DownSample:

DownSample
^^^^^^^^^^

Down-samples a list of clonotype abundance tables by randomly selecting
a pre-defined number of reads or clonotypes. This routine could be useful for

-  normalizing samples to remove certain biases for depth-dependent statistics
-  speeding up computation / decreasing file size and memory footprint.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS DownSample \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+------------+--------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                      |
+=============+=======================+============+==================================================+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See :ref:`common_params`  |
+-------------+-----------------------+------------+--------------------------------------------------+
| ``-x``      | ``--size``            | integer    | Number of reads/clonotypes to take. **Required** |
+-------------+-----------------------+------------+--------------------------------------------------+
| ``-u``      | ``--unweighted``      |            | Will not weight clonotypes by frequency          |
+-------------+-----------------------+------------+--------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                     |
+-------------+-----------------------+------------+--------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                             |
+-------------+-----------------------+------------+--------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs sub-samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``ds:[-x value]`` to ``..filter..`` metadata column.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _FilterNonFunctional:

FilterNonFunctional
^^^^^^^^^^^^^^^^^^^

Filters non-functional (non-coding) clonotypes, i.e. the ones that
contain a stop codon or frameshift in their receptor sequence. Those
clonotypes do not have any functional role, but they are useful for
dissecting and studying the V-(D)-J recombination machinery as they do
not pass thymic selection.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS FilterNonFunctional \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                        |
+=============+=======================+============+====================================================================================================+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See :ref:`common_params`                                                    |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-e``      | ``--negative``        |            | Negative filtering, i.e. only non-functional clonotypes are retained                               |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-e``      | ``--negative``        |            | Negative filtering, i.e. only non-functional clonotypes are retained                               |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
|             | ``--save-freqs``      |            | Don't re-calculate clonotype frequencies and use those from original sample (no re-normalization)  |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                               |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``ncfilter:[retain or remove based on -e option]`` to ``..filter..``
metadata column.

Creates a filter summary file with a ``ncfilter.summary.txt`` suffix
containing info on the number of unique clonotypes that passed the
filtering process, their total frequency and count.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _SelectTop:

SelectTop
^^^^^^^^^

Selects top N clonotypes from the sample. Useful for studying exapanded clonotypes 
and clonotypes with strong convergent recombination bias, as well as robust computing 
of unweighted statistics.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS SelectTop \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                        |
+=============+=======================+============+====================================================================================================+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See :ref:`common_params`                                                    |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-x``      | ``--top``             | integer    | Number of top clonotypes to take. **Required**                                                     |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
|             | ``--save-freqs``      |            | Don't re-calculate clonotype frequencies and use those from original sample (no re-normalization)  |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                                                                       |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                               |
+-------------+-----------------------+------------+----------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs sub-samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``top:[-x value]`` to ``..filter..`` metadata column.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _FilterByFrequency:

FilterByFrequency
^^^^^^^^^^^^^^^^^

Selects clonotypes that either have a frequency above the specified threshold and/or 
constitute more than a specified percent of reads (e.g. quantile threshold of 0.25 will 
top N clonotypes that in total contain 25% of reads in the sample). Those two filters 
can be used together or separately by setting either frequency threshold to 0 or 
quantile threshold to 1.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS FilterByFrequency \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name           | Argument    | Description                                                                                                                                               |
+=============+==========================+=============+===========================================================================================================================================================+
| ``-m``      | ``--metadata``           | path        | Path to metadata file. See :ref:`common_params`                                                                                                           |
+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-f``      | ``--freq-threshold``     | ``0.0-1.0`` | Clonotype frequency threshold. Default is ``0.01``                                                                                                        |
+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-q``      | ``--quantile-threshold`` | ``0.0-1.0`` | Quantile threshold. Will retain a set of top N clonotypes so that their total frequency is equal or less to the specified threshold. Default is ``0.25``  |
+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
|             | ``--save-freqs``         |             | Don't re-calculate clonotype frequencies and use those from original sample (no re-normalization)                                                         |
+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``           |             | Compress output sample files                                                                                                                              |
+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``               |             | Display help message                                                                                                                                      |
+-------------+--------------------------+-------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``freqfilter:[-f value]:[-q value]`` to ``..filter..`` metadata column.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _ApplySampleAsFilter:

ApplySampleAsFilter
^^^^^^^^^^^^^^^^^^^

Retains/filters out all clonotypes found in a given sample **S** from
other samples. Useful when **S** contains some specific cells of interest
e.g. tumor-infiltrating T-cells or sorted tetramer+ T-cells.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS ApplySampleAsFilter \
    [options] [sample1.txt sample2.txt ... if -m is not specified] filter_sample output_prefix

Parameters:

+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument   | Description                                                                                        |
+=============+========================+============+====================================================================================================+
| ``-m``      | ``--metadata``         | path       | Path to metadata file. See :ref:`common_params`                                                    |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-i``      | ``--intersect-type``   | string     | Sample intersection rule. Defaults to ``strict``. See :ref:`common_params`                         |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-e``      | ``--negative``         |            | Negative filtering, i.e. only clonotypes absent in sample *S* are retained                         |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+
|             | ``--save-freqs``       |            | Don't re-calculate clonotype frequencies and use those from original sample (no re-normalization)  |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``         |            | Compress output sample files                                                                       |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``             |            | Display help message                                                                               |
+-------------+------------------------+------------+----------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``asaf:[- if -e, + otherwise]:[-i value]`` to ``..filter..`` metadata
column.

Graphical output
~~~~~~~~~~~~~~~~

none

--------------

.. _FilterBySegment:

FilterBySegment
^^^^^^^^^^^^^^^

Filters clonotypes that have V/D/J segments that match a specified segment set.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS FilterBySegment \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix

Parameters:

+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name        | Argument   | Description                                                                                                 |
+=============+=======================+============+=============================================================================================================+
| ``-m``      | ``--metadata``        | path       | Path to metadata file. See :ref:`common_params`                                                             |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| ``-n``      | ``--negative``        |            | Retain only clonotypes that lack specified V/D/J segments.                                                  |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| ``-v``      | ``--v-segments``      | v1,v2,...  | A comma-separated list of Variable segment names. Non-matching incomplete names will be partially matched.  |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| ``-d``      | ``--d-segments``      | d1,d2,...  | A comma-separated list of Diversity segment names. Non-matching incomplete names will be partially matched. |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| ``-j``      | ``--j-segments``      | j1,j2,...  | A comma-separated list of Joining segment names. Non-matching incomplete names will be partially matched.   |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
|             | ``--save-freqs``      |            | Don't re-calculate clonotype frequencies and use those from original sample (no re-normalization)           |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``        |            | Compress output sample files                                                                                |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+
| ``-h``      | ``--help``            |            | Display help message                                                                                        |
+-------------+-----------------------+------------+-------------------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs filtered samples to the path specified by output prefix and
creates a corresponding metadata file. Will also append
``segfilter:[retain or remove based on -e option]:[-v value]:[-d value]:[-j value]`` 
to ``..filter..`` metadata column.

Creates a filter summary file with a ``segfilter.summary.txt`` suffix
containing info on the number of unique clonotypes that passed the
filtering process, their total frequency and count.

Graphical output
~~~~~~~~~~~~~~~~

none
