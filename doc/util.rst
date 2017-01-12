.. _util:

Utilities
---------

.. _SplitMetadata:

SplitMetadata
^^^^^^^^^^^^^

Splits metadata file into separate metadata files according to the set of values in specified column(s). 
Can be handly for implementing pipelines using VDJtools.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS SplitMetadata [options] metadata.txt output_dir
    
Parameters:

+-------------+------------------------+---------------------+-----------------------------------------------------------------+
| Shorthand   |      Long name         | Argument            | Description                                                     |
+=============+========================+=====================+=================================================================+
| ``-c``      | ``--columns``          | string1,string2,... | A comma separated list of column name(s) to split metadata by.  |
+-------------+------------------------+---------------------+-----------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Output resulting metadata files to specified folder. Unique combinations of metadata entries in specified columns will be appended to names of corresponding metadata files,
relative sample paths will be handled appropriately.

-------------

.. _FilterMetadata:

FilterMetadata
^^^^^^^^^^^^^^

Filters metadata by evaluating expression over values in specified metadata columns, e.g.:

.. code-block:: java

    "__chain__=~/TR[AB]/"
    "__chain__=='TRA'||__chain__=='TRB'"
    "__chain__.contains('TRA')"
    "!__condition__.startsWith('control')"

Both Java and Groovy syntax are supported, column names should be marked by double underscores before and after the name.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS FilterMetadata [options] metadata.txt output_dir output_suffix
    
Parameters:

+-------------+------------------------+--------------+-------------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument     | Description                                                                                                       |
+=============+========================+==============+===================================================================================================================+
| ``-f``      | ``--filter``           | expression   | Filter expression, should be surrounded with quotation marks, metadata column names should be marked with ``__``. |
+-------------+------------------------+--------------+-------------------------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Filtered metadata table with corresponding suffix will be created in the specified folder, relative sample paths will be handled appropriately.

-------------

.. _Convert:

Convert
^^^^^^^

Converts datasets from an arbitrary supported format to :ref:`vdjtools_format`.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS Convert \
    [options] [sample1.txt sample2.txt ... if -m is not specified] output_prefix
    
Parameters:

+-------------+------------------------+-----------+-------------------------------------------------------------------------------------------------------------+
| Shorthand   |      Long name         | Argument  | Description                                                                                                 |
+=============+========================+===========+=============================================================================================================+
| ``-S``      | ``--software``         | path      | Format to convert from, see the :ref:`supported_input` section                                              |
+-------------+------------------------+-----------+-------------------------------------------------------------------------------------------------------------+
| ``-m``      | ``--metadata``         | path      | Path to metadata file. See :ref:`common_params`                                                             |
+-------------+------------------------+-----------+-------------------------------------------------------------------------------------------------------------+
| ``-c``      | ``--compress``         |           | Compressed output for clonotype table. See :ref:`common_params`                                             |
+-------------+------------------------+-----------+-------------------------------------------------------------------------------------------------------------+

Tabular output
~~~~~~~~~~~~~~

Outputs converted samples to the path specified by output prefix and creates a 
corresponding metadata file. Will also append ``conv:[-S value]`` to ``..filter..`` 
metadata column.

-------------

.. _Rinstall:

RInstall
^^^^^^^^

Prints the list of required R packages and installs dependencies into a local library 
(`RPackages` folder) which is placed in the parent folder of VDJtools jar. 
If this routine does not return with "PASSED" message, manual installation of 
packages that failed to deploy is required.

Command line usage
~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    $VDJTOOLS RInstall
