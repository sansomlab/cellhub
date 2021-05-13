'''
Overview
========

This module contains the code for registering and accessing pipeline outputs
from a common location.

There are classes that provide methods for:

(1) registering pipeline outputs to the common service endpoint
(2) discovering the information avaliable from the service endpoint
(3) accessing information from the service endpoint

The service endpoint is the folder "api". We use a rest-like syntax for providing access to piple outputs. The general structure is: ::

api/pipline.name/analysis.name/data.subset.or.slice/data.files
api/pipline.name/analysis.name/data.subset.or.slice/manifest.yml

Where the manifest.yml file contains a description of dataset and file formats.

For example the layout of the api for pipeline_cellranger_multi.py looks like this ::

api/cellranger.multi/gex/filtered/library_id/
api/cellranger.multi/gex/full/library_id/
api/cellranger.multi/adt/filtered/library_id/
api/cellranger.multi/adt/full/library_id/
api/cellranger.multi/hto/filtered/library_id/
api/cellranger.multi/hto/full/library_id/
api/cellranger.multi/vdj/filtered/library_id/
api/cellranger.multi/vdj/full/library_id/


Usage
=====

Registering pipline ouputs on the service endpoint
--------------------------------------------------

Pipelines must register their output datsets on the API using a subclass of the "register" class. Practically this involves:

(1) creation of the appropriate folders in the "api" endpoint folder
(2) symlinking of files to the appropriate folders
(3) construction and deposition of the manifest.yml  file
(4) sanity checking.


Discovering avaliable datasets
------------------------------

TBD.


Accessing datasets
------------------

Datasets can be accessed directly via the endpoint. However it is recommend to access them via a subclass of the "read" class" which should provide basic sanity checking.
'''

import yaml
import os
import shutil
import re
import copy
from pprint import pprint

class register():
    '''Prototype class for registering outputs on the service endpoint'''

    def __init__(self, pipeline = None, endpoint="api"):

        if pipeline is None or pipeline == "":
            raise ValueError("a pipeline name must be specified")

        self.pipeline = re.sub("[ -_]",".",pipeline)
        self.endpoint = endpoint

    def reset(self):
        '''
        Clean the pipeline endpoint
        '''

        endpoint_location = os.path.join(self.endpoint,
                                         self.pipeline)

        if os.path.exists(endpoint_location):

            shutil.rmtree(endpoint_location)

    def dataset(self,
                analysis_name,
                file_set,
#                description,
                file_format,
                analysis_description = None,
                data_subset = None,
                data_id = None
                ):




        # self.description = description # free text description
        self.file_format = file_format # e.g. tsv, market-matrix etc



        self.data_subset = data_subset # e.g. for cell ranger full|filtered
        self.data_id = data_id  # e.g. for cell ranger this is the library_id


        # The file_set is a dictionary that contains the files to be registered
        #
        # { "name": { "path": "path/to/file",
        #             "format": "file-format",
        #             "description": "free-text" }
        #
        # the top level "name" keys are arbitrary and not exposed in the API
        #
        # e.g. for cell ranger output
        #
        # {"barcodes": {"path":"path/to/barcodes.tsv",
        #               "format": "tsv",
        #               "description": "cell barcode file"},
        # {"features": {"path":"path/to/features.tsv",
        #               "format": "tsv",
        #               "description": "features file"},
        # {"matrix": {"path":"path/to/matrix.mtx",
        #               "format": "market-matrix",
        #               "description": "Market matrix file"}
        # }

        # check the files exist

        print("*************************")
        print(file_set)
        for file_name in file_set.keys():
            file_path = file_set[file_name]["path"]
            if not os.path.exists(file_path):
                raise ValueError("file_set file : " + file_name + " does not "
                                 "exist at path: " + file_path)

        self.file_set = file_set # a dictionary of the files to be registered

        self.analysis_name = analysis_name


    def report(self):

        pprint(vars(self))

    def deposit(self):

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


            print("<<<<<<<<<<>>>>>>>>>>")
            print(self.file_set[output_file])

            source_file = os.path.abspath(self.file_set[output_file]["path"])
            link_location = os.path.join(endpoint_location,
                                         os.path.basename(
                                             self.file_set[output_file]["path"]))


            os.symlink(os.path.relpath(source_file,
                                       os.path.dirname(link_location)),
                       link_location)
