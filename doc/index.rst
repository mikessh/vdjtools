.. vdjtools documentation master file, created by
   sphinx-quickstart on Tue Mar 24 17:33:45 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

VDJtools: a framework for post-analysis of repertoire sequencing data
=====================================================================

VDJtools is an open-source Java/Groovy-based framework designed to
facilitate analysis of immune repertoire sequencing
(`RepSeq <http://www.ncbi.nlm.nih.gov/pubmed/22043864>`__) data.
VDJtools computes a wide set of statistics and is able to perform
various forms of cross-sample analysis. Both comprehensive tabular
output and publication-ready plots are provided.

For the period of VDJtools development, there were no other software tools able to 
perform a comprehensive RepSeq post-analysis. Therefore most of the analysis of
this kind was done using in-house scripts, which definitely leads to
"re-inventing the bicycle" problem and loss of analysis reproducibility.

The main aims of the **VDJtools Project** are:

-  To ensure consistency between post-analysis methods and results
-  To save the time of bioinformaticians analyzing RepSeq data
-  To create an API framework facilitating development of new RepSeq
   analysis applications
-  To provide a simple enough command line tool so it could be used by
   immunologists and biologists with little computational background
   
VDJtools source code and binaries are located `here <https://github.com/mikessh/vdjtools>`__.

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   
   intro
   install
   usage
   examples
   input
   modules
   basic
   diversity
   overlap
   preprocess
   operate
   annotate
   util
   vdjviz
