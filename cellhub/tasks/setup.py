"""
setup.py
========

A parent class to help setup pipeline tasks. It can be extended to meet the needs of the different 
pipelines. The class is used to obtain a task object that:

* defines job resource requirements
* provides access to variables (by name or via a .var dictionary)
* creates an outfolder based on the outfile name

"""

import os
import math
import pandas as pd

class setup():
    '''
    A class for routine setup of pipeline tasks.
    
    Args:
        infile: The task infile path or None
        outfile: The task outfile path (typically ends with ".sentinel")
        memory: The total memory needed for execution of the task. If no
            unit is given, gigabytes are assumed. Recognised units are
            "M" for megabyte and "G" for gigabytes. 4 gigabytes can be
            requested by passing "4", "4GB" or "4096M". Default = "4GB".
        cpu: The number of cpu cores required (used to populate job_threads)
        make_outdir: True|False. Default = True.
        expose_var: True|False. 
            Should the self.var dictionary be created from self.__dict__.
            Default = True.

    Attributes:
        job_threads: The number of threads that will be requested
        job_memory: The amount of memory that will be requested per thread
        resources: A dictionary with keys "job_threads" and "job_memory" for populating
            the P.run() kwargs, e.g. ``P.run(statement, **t.resources)``
        outname: The os.path.basename of outfile
        outdir: The os.path.dirname of outfile
        indir: If an infile path is given, the os.path.dirname of the  infile.
        inname: If an infile path is given, the os.path.basename of the infile.
        log_file: If the outfile path ends with ".sentinel"

    '''
    
    def parse_mem(self, memory):
        '''
        Return an integer that represents the amount of memory
        needed by the task in gigabytes.
        '''
    
        if memory in [None, "None", "none", 
                      False, "False", "false", 
                      ""]:
        
            # set the default
            G = 4
        
        elif isinstance(memory, int):
        
            G = int(memory)
        
        elif memory.endswith("G"):
        
            G = int(memory[:-1])
            
        elif memory.endswith("M"):
        
            G = int(memory[:-1]) / 1000
            
        else:
            raise ValueError(
                'Memory request not recognised. Please specify the memory '
                'required in gigabytes (G) or megabytes (M), e.g. "4G" or '
                '"4000M". If a unit is not specified, G will be assumed')
            
        return(G)
    
    
    def set_resources(self, PARAMS, memory="4G", cpu=1):
        '''
        calculate the resource requirements and return a
        dictionary that can be used to update the local variables
        '''

        gb_requested = self.parse_mem(memory)

        # CGAT-core expects memory to be specified per core
        mem_gb = int(math.ceil(gb_requested / float(cpu) ))

        if "resources_mempercore" in PARAMS.keys():
        
            mpc = self.parse_mem(PARAMS["resources_mempercore"])

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
                 make_outdir=True,
                 expose_var=True):
    
        self.set_resources(PARAMS, memory=memory, cpu=cpu)

        # self.outfile = os.path.relpath(outfile)
        # self.var["outfile"] = self.outfile

        self.outdir = os.path.dirname(outfile)          
        self.outname = os.path.basename(outfile)

        if not infile is None:
        
            # self.infile = os.path.relpath(infile)
            # self.var["infile"] = self.infile
            
            self.indir = os.path.dirname(infile) 
            
            self.inname = os.path.basename(infile)

        if make_outdir:

         # take care of making the output directory.
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)

        if outfile.endswith(".sentinel"):
        
            self.log_file = outfile.replace(".sentinel", ".log")
        
        if expose_var:
            self.var = self.__dict__

      

