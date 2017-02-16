Usage
-----

Command line usage
^^^^^^^^^^^^^^^^^^

General way to execute VDJtools routines would be the following,

.. code-block:: bash

    java -Xmx16G -jar vdjtools.jar RoutineName [arguments] -m metadata.txt output/prefix
    
Output prefix could be either an output directory name (if ended with
``/``) or an output file prefix. Most VDJtools routines will append 
the prefix with an intuitive suffix and extension.

The ``-m metadata.txt`` argument specifies a metadata file with relative sample paths, 
sample names and any other information to provide this information later in analysis.
For more details, see the :ref:`metadata` section.

Alternatively, ``-m`` argument could be substituted with a
space-separated list of files, e.g.

.. code-block:: bash

    java -Xmx16G -jar vdjtools.jar RoutineName sample1.txt[.gz] sample2.txt[.gz] ... output/prefix

Whether not explicitly used (such as in "...Plot" routines) and applicable, 
plotting is turned on with ``-p`` argument.

The ``-h`` argument will bring up help message for specified routine.

.. warning:: 

    Consider allocating sufficient memory for Java Virtual Machine
    when running the pipeline. To do so, execute the java with the 
    ``-Xmx`` argument, e.g.: 
    
    .. code-block:: bash
    
        java -Xmx16G -jar vdjtools.jar RoutineName [arguments] 
    
    If insufficient amount memory is allocated, the Java Virtual Machine
    could drop with a *Java Heap Space Out of Memory* error.

.. tip::

    Some routines could be memory demanding, especially when running sample 
    intersection/joining/pooling with a high number of large (~1,000,000 clonotypes)
    datasets. Setting the ``-Xmx`` argument to 20-60Gb of memory should be enough
    for most purposes, e.g. 100 samples with 500,000 clonotypes on average.

    Another way to work this around is to down-sample datasets to ~100,000 reads
    each using the :ref:`downsample` routine.

Importing clonotype tables
^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to proceed with VDJtools analysis datasets should be converted to
VDJtools format (see :ref:`vdjtools_format`). To do this run either of the following commands:

.. code-block:: bash

    java -Xmx16G -jar vdjtools.jar Convert -S software -m metadata.txt ... output_dir/
    
or

.. code-block:: bash

    java -Xmx16G -jar vdjtools.jar Convert -S software sample1.txt[.gz] sample2.txt[.gz] ... output_dir/
    
An additional ``-c`` flag could be set to compress output files.
