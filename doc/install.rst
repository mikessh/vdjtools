.. _install:

Installing VDJtools
-------------------

Installing binaries
^^^^^^^^^^^^^^^^^^^

First make sure that you have installed Java Runtime Environment (JRE) v1.8 by running
``java -version``.  Any recent Linux distribution will provide it via its
package manager.  If not, or if your system is running MacOSX or Windows,
download the JRE from `Oracle <http://java.com/en/download/>`__.

Then download and unpack the VDJtools binaries from the `latest
release <https://github.com/mikessh/vdjtools/releases/latest>`__.

The program is then run by executing the following line:

.. code:: bash

    java -jar path-to-vdjtools-X.X.X.jar

where ``X.X.X`` stands for the VDJtools version (omitted further
for simplicity). This will bring up the list of available routines. To
see the details (parameters, etc) for a specific routine execute

.. code:: bash

    java -jar vdjtools.jar RoutineName -h

Windows
~~~~~~~

Dedicated VDJtools bundle can be downloaded from the
`release <https://github.com/mikessh/vdjtools/releases/latest>`__ section
and is marked with ``.win.zip`` suffix.

Linux
~~~~~

A VDJtools bundle can be downloaded from the
`release <https://github.com/mikessh/vdjtools/releases/latest>`__ section
which includes the required vdjtools.jar file.

All plotting is handled by R and will require several R packages some of
which will be available via your distribution package manager.
See :ref:`install-plotting` below.

MacOS
~~~~~

Installation can be performed using `Homebrew <http://brew.sh/>`__ package manager:

.. code:: bash

    brew tap homebrew/science
    brew tap mikessh/repseq
    brew install vdjtools

Note that this sets ``vdjtools`` as a shortcut for ``java -jar vdjtools-X.X.X.jar``. JVM arguments
such as ``-Xmx`` can be still passed to the script, e.g. ``vdjtools -Xmx20G CalcBasicStats ...``.

.. _install-plotting:

Setting up plotting routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All plotting in VDJtools framework is performed via running R scripts.
Therefore one needs to install `R <http://www.r-project.org/>`__
programming language and several of its packages.  Make sure that

.. code:: bash

    Rscript --version

runs successfully. Note that all R scripts were tested under R version 3.1.0.

The pre-compiled ```*.win.zip`` includes all the required R packages
and the homebrew installation will install them automatically.
In all other cases the required packages need to be manually installed.

These are the required packages:

============  ===================
CRAN package  Debian package
============  ===================
ape           r-cran-ape
circlize
FField
ggplot2       r-cran-ggplot2
gplots        r-cran-gplots
grid
gridExtra
MASS          r-cran-mass
plotrix       r-cran-plotrix
RColorBrewer  r-cran-rcolorbrewer
reshape       r-cran-reshape
reshape2      r-cran-reshape2
scales        r-cran-scales
VennDiagram
============  ===================

If your Linux distribution includes pre-packaged versions of a package,
those should be prefered.  The following will install the existing for
Debian and Debian based distributions such as Ubuntu and Mint:

.. code:: bash

    apt-get install r-cran-ape r-cran-ggplot2 r-cran-gplots r-cran-mass \
      r-cran-plotrix r-cran-rcolorbrewer r-cran-reshape r-cran-reshape2 \
      r-cran-scales

while the other packages will have to be installed via R itself:

.. code:: r

    install.packages(c("circlize", "grid", "gridExtra", "VennDiagram"))

Alternatively, VDJtools has a ref:`Rinstall` routine:

.. code:: bash

    java -jar vdjtools.jar Rinstall

This would also print the list of required R modules, so in case
``Rinstall`` fails, they could be installed manually by running the following
command in R:

.. code:: r

    install.packages(c("reshape2", "FField", "reshape", "gplots",
                       "gridExtra", "circlize", "ggplot2", "grid",
                       "VennDiagram", "ape", "MASS", "plotrix",
                       "RColorBrewer", "scales"))

Note that most issues with package installation can be resolved by switching to correct CRAN mirror.

Dedicated windows binaries already have all R packages bundled, and the options summarized above
should be considered only when troubleshooting R script execution issues.

Compiling from source
^^^^^^^^^^^^^^^^^^^^^

VDJtools could be compiled from source code using `Apache
Maven <http://maven.apache.org/>`__. Compilation should be performed
under JRE v1.8 by running the following commands:

.. code:: bash

    git clone https://github.com/mikessh/vdjtools.git
    cd vdjtools/
    mvn clean install

Binaries could then be found under the ``vdjtools/target/`` folder.
