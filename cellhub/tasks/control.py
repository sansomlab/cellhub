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

    if sys.argv[1] == "make":
        yml_file = pipeline_name.replace(".py", ".yml")

        if not os.path.exists(yml_file):

            cmd = pipeline_name.replace("pipeline_","").split(".")[0]

            raise ValueError('local configuration file missing. Please e.g. run '
                             '"cellhub ' + cmd +  ' config" to check'
                             'out a local copy of the default file')

    else:
        yml_file = os.path.splitext(pipeline_path)[0] + "/" + pipeline_name.replace(".py", ".yml")

        if not os.path.exists(yml_file):
            raise ValueError("default configuration file missing")

    return(yml_file)
