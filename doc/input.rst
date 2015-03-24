Input
=====

Formats
-------

Raw reads with mapped V-(D)-J junctions should assembled into clonotypes
by the RepSeq processing software prior to post-analysis using VDJtools.
VDJtools supports input files in form of tab-delimited clonotype
abundance tables. To fully characterize a **clonotype**, the following
fields are required:

-  Variable (*V*) segment name

-  Diversity (*D*) segment name for some of the receptor chains (TRB,
   TRD and IGH)

-  Joining (*J*) segment name

    NOTE: in case a comma separated list of segment names is provided,
    only the first one is selected

-  Complementarity determining region 3 nucleotide sequence (*CDR3nt*)

-  CDR3 amino acid sequence (*CDR3aa*)

-  In case of non-coding nucleotide sequences, the convention is to
   translate both starting from V segment and backwards from J segment
   up to the CDR3 mid point, producing sequence like ``CASSLAgg?TNEKFF``
   where ``gg?`` marks the incomplete codon (MiGEC). This could also be
   replaced by ``~`` symbol (MiTCR).

Clonotype **abundance** data is provided in the following fields:

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
   boundaries within CDR3 nucleotide sequence

VDJtools also supports parsing somatic hypermutations (*SHMs*) in
antibody sequences, as hypermutation analysis modules are planned to be
added in future versions.

VDJtools could also read `gzip <http://www.gzip.org/>`__-compressed
files, such files should have an ``.gz`` suffix. VDJtools currently
supports four input types, specified by a mandatory ``-S`` argument,
with the following table layouts.

MiTCR
~~~~~

Output from `MiTCR <mitcr.milaboratory.com>`__ software in ``full`` mode
can be used without any pre-processing. The table shoul start with **two
header lines** (default MiTCR output stores processing options and
version in the first line), which are skipped by VDJtools. Clonotypes
are then listed one at a line. Fields marked as <unused> are skipped by
VDJtools but are necessary for correct table layout.

Table example:

+-----------+-----------+------------------------------+
| column1   | column2   | column3                      |
+===========+===========+==============================+
| count     | freq      | CDR3nt                       |
+-----------+-----------+------------------------------+
| 84        | 0.0084    | TGCGCCAGCAG...GAGCTGTTTTTT   |
+-----------+-----------+------------------------------+

+------------+------------+--------------------+------------+
| column4    | column5    | column6            | column7    |
+============+============+====================+============+
| <unused>   | <unused>   | CDR3aa             | <unused>   |
+------------+------------+--------------------+------------+
| .          | .          | CASSQEGTGYSGELFF   | .          |
+------------+------------+--------------------+------------+

+-----------+------------+------------+------------+------------+
| column8   | column9    | column10   | column11   | column12   |
+===========+============+============+============+============+
| V         | <unused>   | J          | <unused>   | D          |
+-----------+------------+------------+------------+------------+
| TRBV4-1   | .          | TRBJ2-2    | .          | TRBD1      |
+-----------+------------+------------+------------+------------+

+------------+------------+------------+------------+
| column13   | column14   | column15   | column16   |
+============+============+============+============+
| Vend       | Dstart     | Dend       | Jstart     |
+------------+------------+------------+------------+
| 16         | 18         | 26         | 31         |
+------------+------------+------------+------------+

To use this input type execute VDJtools routines with ``-S mitcr``
argument.

MiGEC
~~~~~

Default output of `MiGEC <https://github.com/mikessh/migec>`__ software
can be directly used with VDJtools. The table should start with a
**single header line**. The table should have the following columns:

+-----------+--------------+--------------------------------------+
| column1   | column2      | column3                              |
+===========+==============+======================================+
| Count     | Percentage   | CDR3 nucleotide sequence             |
+-----------+--------------+--------------------------------------+
| 198       | 0.0295       | TGTGCGAAAGAGAG...ACGGTATGGACGTCTGG   |
+-----------+--------------+--------------------------------------+

+----------------------------+------------------------------+--------------+--------------+
| column4                    | column5                      | column6      | column7      |
+============================+==============================+==============+==============+
| CDR3 amino acid sequence   | V segments                   | J segments   | D segments   |
+----------------------------+------------------------------+--------------+--------------+
| CAKERYSSSPDYYYYGMDVW       | IGHV3-NL1\ *01\|IGHJ6*\ 01   | .            |              |
+----------------------------+------------------------------+--------------+--------------+

+------------------------------+-------------------------------+------------------------------+-------------------------------+
| column8                      | column9                       | column10                     | column11                      |
+==============================+===============================+==============================+===============================+
| Last V nucleotide position   | First D nucleotide position   | Last D nucleotide position   | First J nucleotide position   |
+------------------------------+-------------------------------+------------------------------+-------------------------------+
| 10                           | .                             | .                            | 31                            |
+------------------------------+-------------------------------+------------------------------+-------------------------------+

To use this input type execute VDJtools routines with ``-S migec``
argument.

IgBlast
~~~~~~~

As IgBlast doesn't provide means for clonotype abundance table
generation, I wrote a small wrapper to handle this task (see
`IgBlastWrapper <https://github.com/mikessh/igblastwrp>`__ repository).
IgBlastWrapper should be run with ``-l 2`` argument to generate the full
table.

+------------+------------+--------------+----------------+
| column1    | column2    | column3      | column4        |
+============+============+==============+================+
| <unused>   | <unused>   | mig\_count   | mig\_freq      |
+------------+------------+--------------+----------------+
| .          | .          | 7            | 0.0103857567   |
+------------+------------+--------------+----------------+

+-----------------------+----------------------+------------------------------+
| column5               | column6              | column7                      |
+=======================+======================+==============================+
| cdr1nt                | cdr2nt               | cdr3nt                       |
+-----------------------+----------------------+------------------------------+
| GGGCACG...ATCTTATCC   | TTTGACC...GCGAAATG   | TGTGCAAGAGG...TTCGACCCCTGG   |
+-----------------------+----------------------+------------------------------+

+------------+------------+------------------------+
| column8    | column9    | column10               |
+============+============+========================+
| cdr1aa     | cdr2aa     | cdr3aa                 |
+------------+------------+------------------------+
| GHGDTILS   | FDPEDGEM   | CARGGSLNVVPPAANWFDPW   |
+------------+------------+------------------------+

+------------+------------+------------+
| column11   | column12   | column13   |
+============+============+============+
| inFrame    | noStop     | complete   |
+------------+------------+------------+
| true       | true       | true       |
+------------+------------+------------+

+--------------------------------+-------------+------------+
| column14                       | column15    | column16   |
+================================+=============+============+
| vSegment                       | dSegment    | jSegment   |
+--------------------------------+-------------+------------+
| IGHV1-24\ *01\|IGHD3-16*\ 01   | IGHJ5\*01   |            |
+--------------------------------+-------------+------------+

+------------+------------+------------+-------------------------------------+
| column17   | column18   | column19   | column20                            |
+============+============+============+=====================================+
| <unused>   | <unused>   | <unused>   | mutations                           |
+------------+------------+------------+-------------------------------------+
| .          | .          | .          | 2:1E0:2:1E0,13:T>A,4:V>E,FW2\|...   |
+------------+------------+------------+-------------------------------------+

See the IgBlastWrapper repository for details on hypermutation encoding.
To use this input type execute VDJtools routines with ``-S igblast``
argument.

Simple
~~~~~~

Simple format is a lightweight clonotype abundance table format that
allows to use data that was generated with software not currently
supported by VDJtools / in-house processed data / incomplete data.

+-----------+-------------+---------------------------+------------------+
| column1   | column2     | column3                   | column4          |
+===========+=============+===========================+==================+
| count     | frequency   | CDR3nt                    | CDR3aa           |
+-----------+-------------+---------------------------+------------------+
| 1176      | 9.90E-02    | TGTGCCAGC...AAGCTTTCTTT   | CASTVDSLDTEAFF   |
+-----------+-------------+---------------------------+------------------+

+------------+-----------+-----------+
| column5    | column6   | column7   |
+============+===========+===========+
| V          | D         | J         |
+------------+-----------+-----------+
| TRBV12-4   | .         | TRBJ1-1   |
+------------+-----------+-----------+

To use this input type execute VDJtools routines with ``-S simple``
argument.

Metadata
--------

Most VDJtools routines could be run with a sample batch. In this case
the paths for the sample files to be analyzed could be provided via
command line, but a more elegant solution is to specify a metadata file
via ``-m`` option. Basic guidelines for creating metadata file are the
following.

-  Metadata file should be a tab-delimited table, e.g.

+-----------------+--------------+-------------+-------+
| #file\_name     | sample\_id   | col\_name   | ...   |
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
-  The file name should be an absolute path
   (``/Users/username/somedir/file.txt``) or a path relative to the
   directory of parent metadata file (``./file.txt``)
-  Sample IDs should be unique

-  Columns after **sample\_id** are treated as metadata entries. There
   are also several cases when info from metadata is used during
   execution:
-  VDJtools plotting routines could be directed to use metadata fields
   for naming samples and creating intuitive legends

   -  Metadata fields are categorized as factor (contain only strings),
      numeric (contain only numbers) and semi-numeric (numbers and
      strings). Numeric and semi-numeric fields could be used for
      gradient coloring by plotting routines.

-  The **time** column is used to specify time points in sequential
   intersection

Note
~~~~

Metadata is introduced to allow easy organizing files to be analyzed.
VDJtools will also append metadata fields to its output tables to
facilitate the exploration of analysis results.

-  When performing tasks that involve modifying clonotype abundance
   tables themselves, such as down-sampling, VDJtools will also provide
   a copy of metadata file pointing to newly generated samples.

-  Newly generated metadata file would contain an additional
   ``..filter..`` column, which has a comma-separated list of filters
   that were applied. For example the downsample routine run with
   ``-n 50000`` will append ``ds:50000`` to the ``..filter..`` column.
   Note that this column name is reserved and should not be modified.
