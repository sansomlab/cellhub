import os
import math
import types
import pandas as pd

# def get(memory="4G", cpu=1):
#     '''calculate the resource requirements and return a
#        dictonary that can be used to update the local variables'''

#     if not memory.endswith("G"):
#         raise ValueError("Memory must be specified as XXG")

#     gb_requested = int(memory[:-1])

#     mem_gb = int(math.ceil(gb_requested / float(cpu) ))

#     spec = {"job_memory": str(mem_gb) + "G",
#             "job_threads": cpu,
#             "r_memory": gb_requested * 1000}

#     return spec


def get_resources(memory="4G", cpu=1, PARAMS={"resources_mempercore":False}):
    '''calculate the resource requirements and return a
       dictionary that can be used to update the local variables'''

    if not memory.endswith("G"):
        raise ValueError("Memory must be specified as XXG")

    gb_requested = int(memory[:-1])
    
    mpc = PARAMS["resources_mempercore"]

    mem_gb = int(math.ceil(gb_requested / float(cpu) ))

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


    spec = {"job_memory": str(mem_gb) + "G",
            "job_threads": ncpu,
            "r_memory": gb_requested * 1000}

    return(spec["job_threads"],
           spec["job_memory"],
           spec["r_memory"])