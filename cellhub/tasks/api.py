'''
API
===

Overview
--------

The primary and secondary analysis pipelines define and register their outputs via a common api. The api comprises of an "api" folder into which pipeline outputs are symlinked (by the "register_dataset" "api" class method).

This module contains the code for registering and accessing pipeline outputs
from a common location.

There are classes that provide methods for:

(1) registering pipeline outputs to the common service endpoint
(2) discovering the information available from the service endpoint (not yet written)
(3) accessing information from the service endpoint (not yet written)

The service endpoint is the folder "api". We use a rest-like syntax for providing access to the pipline outputs.

Usage
-----

Registering outputs on the service endpoint
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- All matrices registered on the API that hold per-cell statistics must have "library_id" and "barcode" columns. The library identifiers must correspond with those given in the pipeline_cellranger_multi.yml file. The barcodes field should contain the untouched Cellranger barcodes.

Please see :doc:`pipeline_cellranger_multi.py</pipelines/pipeline_cellranger_multi>` or :doc:`pipeline_cell_qc.py</pipelines/pipeline_cell_qc>` source code for examples of how to register results on the API.

As an example the code used for registering the qcmetrics outputs is reproduced with some comments here: ::

  import cellhub.tasks.api as api

  file_set={}

  ...

  # the set of files to be registered is defined as a dictionary
  # the keys are arbitrary and will not appear in the api

  file_set[library_id] = {"path": tsv_path,
                          "description":"qcmetric table for library " +\
                                        library_id,
                          "format":"tsv"}

  # an api object is created, passing the pipeline name
  x = api.api("cell.qc")

  # the dataset to be deposited is added
  x.define_dataset(analysis_name="qcmetrics",
                   data_subset="filtered",
                   file_set=file_set,
                   analysis_description="per library tables of cell GEX qc statistics",
                   file_format="tsv")

  # the dataset is linked in to the API
  x.register_dataset()


Discovering available datasets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At present the API can be browsed on the command line. Programmatic access is expected in a future update.


Accessing datasets
^^^^^^^^^^^^^^^^^^

At present datasets are accessed directly via the "api" endpoint.


Class and method documentation
------------------------------

'''

import yaml
import os
import shutil
import re
import copy
from pprint import pprint

class api():
    '''
    A class for defining and registering datasets on the cellhub api.

    When initialising an instance of the class, the pipeline name
    is passed e.g.::

      x = cellhub.tasks.api.register("cell_qc")

    .. note:: pipeline names are sanitised to replace spaces, underscores and hypens with periods.
    '''

    def __init__(self, pipeline = None, endpoint="api"):

        if pipeline is None or pipeline == "":
            raise ValueError("a pipeline name must be specified")

        self.pipeline = re.sub("[ \-_]",".",pipeline)
        self.endpoint = endpoint
        self.dataset_defined = False

    def define_dataset(self,
                    analysis_name = None,
                    analysis_description = None,
                    data_subset = None,
                    data_id = None,
                    data_format = None,
                    file_set = None,
    ):
        '''
        Define the dataset.

        The "data_subset", "data_id" and "data_format" parameters are optional.

        The file_set is a dictionary that contains the files to be registered: ::

          { "name": { "path": "path/to/file",
                      "format": "file-format",
                      "link_name": "api link name", # optional
                      "description": "free-text" }

        the top level "name" keys are arbitrary and not exposed in the API

        e.g. for cell ranger output the file_set dictionary looks like this: ::

          {"barcodes": {"path":"path/to/barcodes.tsv",
                        "format": "tsv",
                        "description": "cell barcode file"},
          {"features": {"path":"path/to/features.tsv",
                        "format": "tsv",
                        "description": "features file"},
          {"matrix": {"path":"path/to/matrix.mtx",
                        "format": "market-matrix",
                        "description": "Market matrix file"}
          }

        '''

        if analysis_name is None:
            raise ValueError("The analysis name must be specified")

        if analysis_description is None:
            raise ValueError("The analysis description must be specified")

        if file_set is None:
            raise ValueError("The file_set  must be specified")

        self.data_subset = data_subset # e.g. for cell ranger full|filtered
        self.data_id = data_id  # e.g. for cell ranger this is the library_id
        self.data_format = data_format # e.g. "mtx", "bam", "h5" as necessary

        # check the files exist
        for file_name in file_set.keys():
            file_path = file_set[file_name]["path"]
            if not os.path.exists(file_path):
                raise ValueError("file_set file : " + file_name + " does not "
                                 "exist at path: " + file_path)

        self.file_set = file_set # a dictionary of the files to be registered

        self.analysis_name = analysis_name

        self.dataset_defined = True


    def register_dataset(self):
        '''
        Register the dataset on the service endpoint. The method:

        1. creates the appropriate folders in the "api" endpoint folder
        2. symlinks the source files to the target location
        3. constructs and deposits the manifest.yml file

        The location at which datasets will be registered is defined as: ::

          api/pipeline.name/analysis_name/[data_subset/][data_id/][data_format/]

        (data_subset, data_id and data_format are [optional])

        '''

        if not self.dataset_defined:

            raise ValueError("A dataset must be defined (register_dataset) before"
                             " the register_datast method is called")

        # 1. create the output folders
        endpoint_location = os.path.join(self.endpoint,
                                         self.pipeline,
                                         self.analysis_name)

        if self.data_subset is not None:

            endpoint_location = os.path.join(endpoint_location,
                                             self.data_subset)

        if self.data_id is not None:

            endpoint_location = os.path.join(endpoint_location,
                                             self.data_id)

        if self.data_format is not None:

            endpoint_location = os.path.join(endpoint_location,
                                             self.data_format)

        # ensure the endpoint location is clean
        if os.path.exists(endpoint_location):
            shutil.rmtree(endpoint_location)

        os.makedirs(endpoint_location)

        # 2. construct and write the manifest
        out_file_set = copy.deepcopy(self.file_set)

        for x in out_file_set.keys():

            out_file_set[x]["path"] = os.path.basename(
                out_file_set[x]["path"])

        manifest = { "pipeline": self.pipeline,
                     "analysis": self.analysis_name,
                     "analysis_description": self.analysis_name,
                     "data_subset": self.data_subset,
                     "identifier": self.data_id,
                     "files" : out_file_set }

        manifest_file = os.path.join(endpoint_location, "manifest.yml")
        with open(manifest_file, "w") as yml:
            yaml.dump(manifest, yml, allow_unicode=True)


        # 3. link in the files
        for output_file in self.file_set.keys():

            source_file = os.path.abspath(self.file_set[output_file]["path"])

            if "link_name" in self.file_set[output_file].keys():
                link_location = os.path.join(endpoint_location,
                                             self.file_set[output_file]["link_name"])
            else:
                link_location = os.path.join(endpoint_location,
                                             os.path.basename(
                                             self.file_set[output_file]["path"]))

            os.symlink(os.path.relpath(source_file,
                                       os.path.dirname(link_location)),
                       link_location)


    def show(self):
        '''
        Print the api object for debugging.
        '''

        pprint(vars(self))


    def reset_endpoint(self):
        '''
        Clean the dataset endpoint
        '''

        endpoint_location = os.path.join(self.endpoint,
                                         self.pipeline)

        if os.path.exists(endpoint_location):

            shutil.rmtree(endpoint_location)
