CodingGuidelines
================

In the guidelines below "xxx" denotes the name of a pipeline such as e.g. "cell_qc".

1. Paths should never be hardcoded. They should be read from the configuation file.
2. Yaml configuration files should be named pipeline_xxx.yml
3. The output of individual pipelines should be written to a subfolder name "xxx.dir" to keep the root directory clean (it should only contain these directories and the yml configuration files!).
4. Python code should pass pep8 checks, this will be enforced at some point in the future.
5. Pipelines should be documented using the sphinx "autodocs" module
6. For the autodocs documentation system to work pipelines should not read or write files outside of ruffus tasks (see below).


Writing pipelines
-----------------

The pipelines live in the "pipelines" folder/module.

Auxillary task functions live in the "pipelines/task" folder/module.

If you need to read or write files outside of ruffus tasks, for example in generator functions it is essential to test that the script has not been imported e.g.::

  if __name__ == "__main__":
     # read file
     yield(inputs, outputs)
  else:
     # don't read file
     yield(None, None)

The reason for this is that the sphinx autodocs module needs to import the piplines to build the documentation and this will fail if the pipeline attempts to read or write files from/to non-existent paths when imported.


Yaml configuration files
------------------------

The cgat-core system only supports configuration files name "pipeline.yml".

We work around this by overriding the cgat-core functionality using a helper function in pipelines/task/control.py as follows::

  import Pipeline as P
  import tasks.control as C

  # Override function to collect config files
  P.control.write_config_files = C.write_config_files

Default yml files must be located at the path pipelines/pipeline_xxx/pipeline_xxx.yml


Writing documentation
---------------------

The sources files for the documentation are found in::

  docsrc

The documentation source files for the pipelines can be found in::

  docsrc/pipelines

To build the documentation cd to the docsrc folder and run::

  make github

This will build the documentation and copy the latex output to the "docs" folder. You then need to cd to the "docs" folder and run::

  make

To compile the pdf.

When the repo is made public we will switch to using html documentation on readthedocs. Unfortunately there is no straightforward solution for private html hosting.
