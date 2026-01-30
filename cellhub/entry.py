"""
cellhub.py - single cell data warehousing and analysis pipelines
================================================================

:Tags: Single-cell

To use a specific workflow, type::

    cellhub <workflow> [workflow options] [workflow arguments]

For this message and a list of available keywords type::

    cellhub --help

To get help for a specify workflow, type::

    cellhub <workflow> --help
"""

import os
import subprocess
import sys
import re
import glob
# import imp

# import pipelines


def print_list_in_columns(list2print, ncolumns):
    """output list *l* in *ncolumns*."""
    list_len = len(list2print)
    assert list_len > 0, "attempt to print an empty list."

    # compute the alignment width
    max_width = max([len(x) for x in list2print]) + 3

    # compute the total number of rows
    nrows = list_len // ncolumns
    if list_len % ncolumns != 0:
        nrows += 1

    # build columns
    columns = [list2print[x * nrows : x * nrows + nrows] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(list2print) % ncolumns)):
        columns[-(x + 1)].append("")

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = "%-" + str(max_width) + "s"
    pattern = " ".join([p for x in range(ncolumns)])

    # put it all together
    return "\n".join([pattern % row for row in rows])


def main_cgat(proj_dir, cmdargs=None):
    command = cmdargs[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)

    if cmdargs[2] == "profile":
        import cellhub.tasks.profile as p

        p.profile(pipeline + ".log", show_fields=False)
        return

    # remove 'cellhub' from sys.argv
    del cmdargs[0]

    # specify a named logfile
    sys.argv.append("--pipeline-logfile=" + pipeline + ".log")

    # (file, pathname, description) = imp.find_module(pipeline, [proj_dir])

    # module = imp.load_module(pipeline, file, pathname, description)
    # module.main(cmdargs)


def main_smk(proj_dir, cmdargs):
    command, subcommand, addit_ops = cmdargs[1], cmdargs[2], cmdargs[3:]
    snakefile = os.path.join(proj_dir, "Snakefile")
    if not "--cores" in addit_ops:
        addit_ops = ["--cores=1"] + addit_ops
    if not "--jobs" in addit_ops:
        addit_ops = ["--jobs=1"] + addit_ops
    cmd = [
        "snakemake",
        "-s",
        snakefile,
        "--config",
        f"target={command}",
        f"mode={subcommand}",
    ] + addit_ops
    print(" ".join(cmd))
    subprocess.run(cmd)


def main(argv=None):
    cmdargs = sys.argv
    print(cmdargs)

    # paths to look for pipelines:
    path = os.path.abspath(os.path.dirname(__file__))

    if len(cmdargs) == 1 or cmdargs[1] == "--help" or cmdargs[1] == "-h":
        pipelines = []
        pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print((globals()["__doc__"]))
        print("The list of available pipelines are:\n")
        cmdlist = print_list_in_columns(
            sorted(
                [os.path.basename(x)[len("pipeline_") : -len(".py")] for x in pipelines]
            ),
            3,
        )
        print(f"{cmdlist}\n")
        return

    if cmdargs[1] in ["annotation", "cluster"]:
        assert sys.version_info >= (
            3,
            11,
        ), f"Python >= 3.11 required, found {sys.version}."
        main_smk(proj_dir=path, cmdargs=cmdargs)
    else:
        assert sys.version_info < (
            3,
            11,
        ), f"Python < 3.11 required, found {sys.version}."
        main_cgat(proj_dir=path, cmdargs=cmdargs)


if __name__ == "__main__":
    sys.exit(main())
