Install
=======

Getting binaries
----------------

Core package
~~~~~~~~~~~~

First make sure that Java Runtime Environment (JRE) v1.7+, available
`here <http://www.oracle.com/technetwork/java/javase/downloads/java-se-jre-7-download-432155.html>`__,
is installed on your system. Then download and unpack the binaries from
the `latest
release <https://github.com/mikessh/vdjtools/releases/latest>`__.

The program is then run by executing the following line

::

    java -jar vdjtools-X.X.X.jar

where ``X.X.X`` stands for the latest VDJtools version (omitted further
for simplicity). This will bring up the list of available routines. To
see the details (parameters, etc) for a specific routine execute

::

    java -jar vdjtools.jar RoutineName -h

Additional argument (``-Xmx``) that sets the memory limit for Java
Virtual Machine should be set for most cases. For example,

::

    java -Xmx16G -jar vdjtools.jar RoutineName [arguments]

Plotting routines
~~~~~~~~~~~~~~~~~

All plotting in VDJtools framework is performed via running R scripts.
Therefore one needs to install `R <http://www.r-project.org/>`__
programming language and make sure that the

::

    Rscript --version

runs successfully. All R scripts were tested under R version 3.1.0.

Then install all dependencies to a local library automatically using

::

    java -jar vdjtools.jar Rinstall

This would also print the list of required R modules, so in case
``Rinstall`` fails, they could be installed manually. To skip this step
run the following command in R:

.. code:: r

    install.packages(c("ggplot2","reshape2","grid",
                       "ape","MASS","plotrix",
                       "RColorBrewer","FField","reshape",
                       "gplots","gridExtra","circlize"))

Compiling from source
---------------------

VDJtools could be compiled from source code using `Apache
Maven <http://maven.apache.org/>`__. Compilation should be performed
under Java v1.7.

The only dependencies that need to be manually installed are
`VDJdb <https://github.com/mikessh/vdjdb>`__ and
`MiLib <https://github.com/milaboratory/milib>`__, used for clonotype
annotation purposes:

.. code:: bash

    git clone --branch 1.0 --depth 1 https://github.com/milaboratory/milib.git
    cd milib && mvn clean install && cd ..
    git clone https://github.com/mikessh/vdjdb.git
    cd vdjdb
    mvn clean install

Then proceed with compiling VDJtools itself:

.. code:: bash

    git clone https://github.com/mikessh/vdjtools.git
    cd vdjtools
    mvn clean install

Binaries could then be found at ``vdjtools/target/`` directory
