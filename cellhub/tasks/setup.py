"""
setup.py
========

A parent class to hold routine task variables. Can be extended to meet the needs of the different 
pipelines.

Variables can be accessed by name from class instances or from the var dictionary.

"""

import os
import math
import types
import pandas as pd

class setup():
    '''
    A class for routine setup of pipeline task by:
    
    - defining commonly used variables
    - creating essential folders
    
    '''
    
    def set_resources(self, PARAMS, memory="4G", cpu=1):
        '''
        calculate the resource requirements and return a
        dictionary that can be used to update the local variables
        '''

        if not memory.endswith("G"):
            raise ValueError("Memory must be specified as XXG")

        gb_requested = int(memory[:-1])
        
        mem_gb = int(math.ceil(gb_requested / float(cpu) ))


        if "resources_mempercore" in PARAMS.keys():
            mpc = PARAMS["resources_mempercore"]

            if not mpc:
            
                ncpu = cpu
            
            else:
                # memory is requested via the core count
                # the number of GB per core is given in the resources_mempercore
                # (note that the schedule will simply ignore the job_memory)
                if not isinstance(mpc, int):
                    raise ValueError("resources_mempercore must be False or integer")

                cpu_needed_for_mem_request = math.ceil(gb_requested / mpc)
                ncpu = max(cpu, cpu_needed_for_mem_request)
        
        else:
            ncpu = cpu

        self.job_memory = str(mem_gb) + "G"
        self.job_threads = ncpu
        self.r_memory = gb_requested * 1000
        self.resources = {"job_memory": self.job_memory,
                          "job_threads": self.job_threads}


    def __init__(self, infile, outfile, PARAMS,
                 memory="4G", cpu=1,
                 make_outdir=True):
    
        self.set_resources(PARAMS, memory=memory, cpu=cpu)

        self.var = {}

        # self.outfile = os.path.relpath(outfile)
        # self.var["outfile"] = self.outfile

        self.outdir = os.path.dirname(outfile) 
        self.var["outdir"] = self.outdir
         
        self.outname = os.path.basename(outfile)
        self.var["outname"] = self.outname


        if not infile is None:
        
            # self.infile = os.path.relpath(infile)
            # self.var["infile"] = self.infile
            
            self.indir = os.path.dirname(infile)
            self.var["indir"] = self.indir    
            
            self.inname = os.path.basename(infile)
            self.var["inname"] = self.inname

        if make_outdir:

         # take care of making the output directory.
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)

        if outfile.endswith(".sentinel"):
        
            self.log_file = outfile.replace(".sentinel", ".log")
            self.var["log_file"] = self.log_file

      

