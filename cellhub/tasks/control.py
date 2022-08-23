from cgatcore import experiment as E
from pathlib import Path
import shutil
import os
import sys

def write_config_files(pipeline_path, general_path):
    '''
    Retrieve default yaml configuration file from:
    cellhub/yaml/pipeline_name.yml
    '''

    paths = [os.path.join(Path(pipeline_path).parents[0], "yaml")]
    config_files = [ os.path.basename(pipeline_path) + ".yml" ]

    for dest in config_files:
        if os.path.exists(dest):
            E.warn("file `%s` already exists - skipped" % dest)
            continue

        for path in paths:
            src = os.path.join(path, dest)
            if os.path.exists(src):
                shutil.copyfile(src, dest)
                E.info("created new configuration file `%s` " % dest)
                break
        else:
            raise ValueError(
                "default config file `%s` not found in %s" %
                (config_files, paths))


def get_parameter_file(pipeline_path):
    '''
    Return the local yml file path if the pipeline is being executed,
    otherwise return the location of the yml file in the repo.
    '''

    pipeline_name = os.path.basename(pipeline_path)


    if len(sys.argv) > 1:
    
        if sys.argv[1] == "make":
        
            # We require a local configuration file to exist

            yml_file = pipeline_name.replace(".py", ".yml")
            E.info("Using local yml file: " + yml_file)

            if not os.path.exists(yml_file):

                cmd = pipeline_name.replace("pipeline_","").split(".")[0]

                raise ValueError('local configuration file missing. Please e.g. run '
                                '"cellhub ' + cmd +  ' config" to check'
                                'out a local copy of the default file')
                
                
        elif sys.argv[1] in ["config", "-M"] :
            
            # necessary to:
            # - allow the file to be parsed by the interpreter (config)
            # - allow sphinx autodocs call to importlib.import_module (-M)     
            
            E.info("Using the default configuration file")
            
            yml_file = os.path.join(os.path.dirname(pipeline_path),
                                  "yaml",
                               pipeline_name.replace(".py", ".yml"))
            
            if not os.path.exists(yml_file):
                raise ValueError("default configuration file missing")
            
        else:
        
            raise ValueError('pipeline command not recognised: ' + sys.argv[1])
            

    else:
    
        # Catch all.
        
        yml_file = os.path.join(os.path.dirname(pipeline_path),
                                "yaml",
                                pipeline_name.replace(".py", ".yml"))
        E.warn("Using default yml file: " + yml_file)

        if not os.path.exists(yml_file):
            raise ValueError("default configuration file missing")
            

    return(yml_file)





