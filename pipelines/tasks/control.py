from cgatcore import experiment as E
import shutil
import os
import sys

def write_config_files(pipeline_path, general_path):
    '''create default configuration files in `path`.
    '''

    paths = [pipeline_path, general_path]
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


def get_parameter_file(pipeline_path, name):
    '''
    Return the local yml file path if the pipeline is being executed,
    otherwise return the location of the yml file in the repo.
    '''

    pipeline_name = os.path.basename(pipeline_path)

    if name == "__main__" and sys.argv[1] != "config":
        return(pipeline_name.replace(".py", ".yml"))

    else:
        return(os.path.splitext(pipeline_path)[0] + "/" + pipeline_name.replace(".py", ".yml"))
