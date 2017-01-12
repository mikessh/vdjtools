Input
-----

Clonotype tables
^^^^^^^^^^^^^^^^

The processing stage of RepSeq analysis starts with mapping of Variable, 
Diversity and Joining segments. Mapped reads are then assembled into clonotypes
and stored as a clonotype abundance tables.

.. _clonotype_spec:

Clonotype
~~~~~~~~~

VDJtools **clonotype** specification includes the following fields:

-  Variable (*V*) segment name.

-  Diversity (*D*) segment name for some of the receptor chains (TRB,
   TRD and IGH). Set to `.` if not aplicable or D segment was not
   identified.

-  Joining (*J*) segment name.

-  Complementarity determining region 3 nucleotide sequence (*CDR3nt*).
   CDR3 starts with Variable region reference point (conserved Cys residue) 
   and ends with Joining segment reference point (conserved Phe\Try).

-  Translated CDR3 sequence (*CDR3aa*).

-  Somatic hypermutations (*SHMs*) in the variable segment (antibody only, **planned**).

.. important::
   For ambiguous segment assignments encoded by a comma separated list 
   of segment names only the first one is selected.

.. hint::
   In case of non-coding CDR3 sequences, the convention is to
   translate in both directions: upstream from V segment 
   reference point and downstream from J segment reference point.
   The resulting sequence (e.g. ``CASSLA_TNEKFF``) 
   is linked by a ``_`` symbol that marks the incomplete codon.

Clonotype **abundance** data is represented by *count* and *frequency* fields:

-  *Count*: number of reads or cDNA/DNA molecules in case
   `UMIs <https://github.com/mikessh/migec#migec-molecular-identifier-group-based-error-correction-pipeline>`__
   are used.

-  *Frequency*: the share of clonotype in the sample. While seemingly
   redundant, this property is left for compatibility with cases when
   the sample represents a subset of another one, e.g. clonotypes from
   PBMCs filtered by intersection with lymph node clonotypes.

The following fields are optional, but are used for computing various
statistics and visualization:

-  *Vend*, *Dstart*, *Dend* and *Jstart* - marking V, D and J segment
   boundaries within CDR3 nucleotide sequence (inclusive)

.. tip::
   VDJtools accepts `gzip <http://www.gzip.org/>`__-compressed
   files, such files should have an ``.gz`` suffix. Input data 
   should be provided in a form of tab-delimited table.

.. _vdjtools_format:

VDJtools format
^^^^^^^^^^^^^^^

This is a core tabular format for VDJtools. All datasets 
should be converted to this format using the :ref:`convert` routine 
prior to analysis. Columns 8-10 are optional.

+-----------+-------------+---------------------------+------------------+------------+-----------+-----------+------------+-----------+-----------+-----------+
| column1   | column2     | column3                   | column4          | column5    | column6   | column7   | column8    | column9   | column10  | column11  |
+===========+=============+===========================+==================+============+===========+===========+============+===========+===========+===========+
| count     | frequency   | CDR3nt                    | CDR3aa           | V          | D         | J         | Vend       | Dstart    | Dend      | Jstart    |
+-----------+-------------+---------------------------+------------------+------------+-----------+-----------+------------+-----------+-----------+-----------+
| 1176      | 9.90E-02    | TGTGCCAGC...AAGCTTTCTTT   | CAST...EAFF      | TRBV12-4   | TRBD1     | TRBJ1-1   | 11         | 14        | 16        | 23        |
+-----------+-------------+---------------------------+------------------+------------+-----------+-----------+------------+-----------+-----------+-----------+

All additional columns after column 10 will be considered as clonotype annotations 
and carried over unmodified during most stages of VDJtools analysis. This is especially 
useful when processing results of :ref:`Annotate` and :ref:`ScanDatabase` routines.

.. _supported_input:

Formats supported for conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MiTCR
~~~~~

Output from MiTCR software (`executable jar <http://files.milaboratory.com/mitcr/1.0.3.3-beta/mitcr.jar>`__, 
`documentation <http://files.milaboratory.com/mitcr/Manual.pdf>`__) in 
``full`` mode can be used without any pre-processing. Corresponding 
table should start with **two header lines** (default MiTCR output 
stores processing options and version in the first line), followed by a clonotype 
list.

Run :ref:`convert` routine with ``-S mitcr`` argument to prepare datasets 
in this format for VDJtools analysis.

MiGEC
~~~~~

`MiGEC <https://github.com/mikessh/migec>`__ is a software for V/D/J mapping and CDR3 
extraction that relies on BLAST algorithm for running alignments. MIGEC software 
additionally implements processing of unique molecular identifier (UMI)-tagged libraries 
for error correction and dataset normalization. Default output of MIGEC software 
can be directly used with VDJtools.

Run :ref:`convert` routine with ``-S migec`` argument to prepare datasets 
in this format for VDJtools analysis.

IgBlast (MIGMAP)
~~~~~~~~~~~~~~~~

As IgBlast doesn't compute a canonical clonotype abundance table, 
VDJtools supports output of `MIGMAP <https://github.com/mikessh/igblastwrp>`__, 
a versatile IgBlast wrapper. Note that currently no somatic hypermutation (SHM) 
information is imported by VDJtools, neither there are any dedicated VDJtools 
routines to analyze SHM profiles, but you check out `post-analysis provided by MIGMAP <https://github.com/mikessh/migmap/tree/master/post>`__.

Run :ref:`convert` routine with ``-S migmap`` argument to prepare datasets 
in this format for VDJtools analysis.

ImmunoSEQ
~~~~~~~~~

One of the most commonly used RepSeq data format, more than 90% of recently published studies  
were performed using `immunoSEQ <http://www.adaptivebiotech.com/content/immunoseq-0>`__ 
assay. We have implemented a parser for clonotype tables as provided by 
`Adaptive Biotechnologies <http://www.adaptivebiotech.com/>`__.

-  The resulting datasets for most studies that use ImmunoSEQ technology can be accessed and exported using the 
   `ImmunoSEQ Analyzer <https://clients.adaptivebiotech.com/>`__.

-  Example datasets in this format could be found in the 
   `Supplementary Data <http://ard.bmj.com/content/suppl/2014/12/11/annrheumdis-2014-206226.DC1/annrheumdis-2014-206226supp_tcr-primary-data.zip>`__ 
   section of `Spreafico R et al. Ann Rheum Dis. 2014 <http://ard.bmj.com/content/early/2014/12/11/annrheumdis-2014-206226.full>`__.

-  Column header information was taken from **page 24** of the immunoSEQ Analyzer 
   `manual <https://clients.adaptivebiotech.com/assets/downloads/immunoSEQ_AnalyzerManual.pdf>`__

-  VDJtools will use V/J segment information only at the family level, as many of the clonotypes miss 
   segment (`-X`) and allele (`-X*0Y`) information. 
   The clonotype table is then collapsed to handle unique V/J/CDR3 entries.

-  Raw clonotype tables in this format do not contain CDR3 nucleotide sequence. 
   Instead, an entire sequencing read (first column) is provided. Therefore, we have 
   implemented additional algorithms for CDR3 extraction and "virtual" translation 
   to tell out-of-frame clonotypes from partially read ones.

.. attention::
   Some of the clonotype entries will dropped during conversion as they contain an incomplete 
   CDR3 sequence (lacking J segment), which is due to short reads used in immunoSEQ assay, 
   see this `blog post <http://www.immunoseq.com/comparing-adaptive-data-and-imgt-data-on-cdr3-region-amino-acid-sequences/>`__ 
   for details.
   
Run :ref:`convert` routine with ``-S immunoseq`` argument to prepare datasets 
in this format for VDJtools analysis. This option should be used in case you have selected 
``Export samples`` option in the ImmunoSEQ analyzer. In case you have used the ``Export samples v2`` 
option you should pass the ``-S immunoseqv2`` argument to VDJtools Convert routine.
   
IMGT/HighV-QUEST 
~~~~~~~~~~~~~~~~

Another commonly used RepSeq processing tool is the 
`IMGT/HighV-QUEST <http://www.imgt.org/IMGTindex/IMGTHighV-QUEST.html>`__ web server.

Please refer to the official `documentation <http://www.imgt.org/HighV-QUEST/help.action?section=doc>`__ 
to see the description of output files and their formats.

.. tip:: 
    The output for each submission consists of several files and only 
    
    .. code:: bash
    
        3_Nt-sequences_${chain}_${sx}_${date}.txt
    
    should be used as an input for VDJtools :ref:`convert` routine. 
    
Run :ref:`convert` routine with ``-S imgthighvquest`` argument to prepare datasets 
in this format for VDJtools analysis.

VDJdb
~~~~~

VDJtools has native support for the analysis of clonotype tables annotated 
with `VDJdb <https://github.com/antigenomics/vdjdb-standalone>`__ software. 
Note that as those tables can list the same clonotype several times with 
different annotation, they should not be used directly in most VDJtools 
routines (e.g. diversity statistics), check out 
`VDJdb README <https://github.com/antigenomics/vdjdb-standalone#some-notes>`__ 
for corresponding guidelines and workarounds.

Vidjil
~~~~~~

VDJtools supports parsing output Json files produced by the 
`Vidjil <http://www.vidjil.org/>`__ software. VDJtools will only use 
top clonotypes which have V/D/J detalization in the output.

RTCR
~~~~

VDJtools supports parsing the ``results.tsv`` table with clonotype list 
generated by the `RTCR <https://github.com/uubram/RTCR>`__ software.

Run :ref:`convert` routine with ``-S rtcr`` argument to prepare datasets 
in this format for VDJtools analysis.

MiXCR
~~~~~

Output from `MiXCR <https://github.com/milaboratory/mixcr>`__ software ``export`` routine 
in ``full`` (default) mode can be used without any pre-processing. 

Run :ref:`convert` routine with ``-S mixcr`` argument to prepare datasets 
in this format for VDJtools analysis.

IMSEQ
~~~~~

Output from `IMSEQ <http://www.imtools.org/>`__ software can be used 
if results are collapsed to nucleotide-level clonotypes using
``-on`` argument with IMSEQ. 

Run :ref:`convert` routine with ``-S imseq`` argument to prepare datasets 
in this format for VDJtools analysis.

.. _metadata:

Metadata
^^^^^^^^

Most VDJtools routines could be run with a sample batch. In this case
paths to input files could be provided via command line (space separated), 
but a more elegant solution is to specify a metadata file via ``-m`` option.
The primary purpose of a metadata file is to organize and annotate datasets.

.. note::
   -  VDJtools will append metadata fields to its output tables to
      facilitate the exploration of analysis results.
      
   -  Metadata entries are used as a factor in some analysis routines and
      most plotting routines.

   -  When performing tasks that involve modifying clonotype abundance
      tables themselves, such as down-sampling, VDJtools will also provide
      a copy of metadata file pointing to newly generated samples.

   -  Newly generated metadata file would contain an additional
      ``..filter..`` column, which has a comma-separated list of filters
      that were applied. For example the :ref:`downsample` routine run with
      ``-n 50000`` will append ``ds:50000`` to the ``..filter..`` column.
      Note that this column name is reserved and should not be modified.

Below are the basic guidelines for creating a metadata file.

-  Metadata file should be a tab-delimited table, e.g.

    +-----------------+--------------+-------------+-------+
    | #file.name      | sample.id    | col.name    | ...   |
    +=================+==============+=============+=======+
    | sample\_1.txt   | sample\_1    | A           | ...   |
    +-----------------+--------------+-------------+-------+
    | sample\_2.txt   | sample\_2    | A           | ...   |
    +-----------------+--------------+-------------+-------+
    | sample\_3.txt   | sample\_3    | B           | ...   |
    +-----------------+--------------+-------------+-------+
    | sample\_4.txt   | sample\_4    | C           | ...   |
    +-----------------+--------------+-------------+-------+
    | ...             | ...          | ...         | ...   |
    +-----------------+--------------+-------------+-------+

-  Header is mandatory, first two columns should be named **file\_name**
   and **sample\_id**. Names of the remaining columns will be later used
   to specify metadata variable name

-  First two columns should contain the file name and sample id
   respectively.
   
   -  The file name should be either an absolute path
      (e.g. ``/Users/username/somedir/file.txt``) or a path relative to the
      parent directory of metadata file (e.g. ``../file.txt``)
   
   -  Sample IDs should be unique

-  Columns after **sample.id** are treated as metadata entries. There
   are also several cases when info from metadata is used during
   execution:
   
   -  VDJtools plotting routines could be directed to use metadata fields
      for naming samples and creating intuitive legends. If column name
      contains spaces it should be quoted, e.g. ``-f "patient id"``

   -  Metadata fields are categorized as factor (contain only strings),
      numeric (contain only numbers) and semi-numeric (numbers and
      strings). Numeric and semi-numeric fields could be used for
      gradient coloring by plotting routines.
