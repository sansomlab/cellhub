Installation
============


Dependencies
------------

Core dependencies include:

- Cellranger (from 10X Genomics) >= 6.0
- Python3
- Various Python packages (see python/requirements.txt)
- R >= 4.0
- Various R libraries (see R/install.packages.R)
- Latex
- The provided cellhub R library


Installation
------------

1. Install the cgat-core pipeline system following the instructions here `https://github.com/cgat-developers/cgat-core/ <https://github.com/cgat-developers/cgat-core/>`_.

2. Clone and install the cellhub-devel repository e.g.

.. code-block:: Bash
     
     git clone https://github.com/COMBATOxford/cellhub-devel.git
     cd cellhub-devel
     python setup.py develop

.. note:: Running "python setup.py develop" is necessary to allow pipelines to be launched via the "cellhub" command.

3. In the same virtual or conda environment as cgat-core install the required python packages::

     pip install -r cellhub-devel/python/requirements.txt

4. To install the required R packages (the "BiocManager" and "devtools" libraries must be pre-installed)::

     Rscript cellhub-devel/R/install.packages.R
     
5. Install the cellhub R library::

     R CMD INSTALL R/cellhub

.. note:: On some systems, automatic detection of the HDF5 library by the R package "hdf5r" (a dependency of loomR) is problematic. This can be worked around by explicitely passing the path to your HDF5 library h5cc or h5pcc binary, e.g.

.. code-block:: Bash

		install.packages("hdf5r", configure.args="--with-hdf5=/apps/eb/dev/skylake/software/HDF5/1.10.6-gompi-2020a/bin/h5pcc")
