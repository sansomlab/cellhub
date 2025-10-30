#!/usr/bin/env bash

usage() { 
    echo "Usage: $0 <CMD> <SUBCMD>" >&2
    echo "\tCMD: annotation, cluster, etc." >&2
    echo "\tSUBCMD: genconfig or execute." >&2
    exit 1
}

if [[ $# -lt 2 ]]; then
    usage
fi

cmd="$1"
subcmd="$2"

shift 2

CELLHUB_HOME=/users/sansom/vlw740/devel/cellhub/ # to set this when installing cellhub

# To run with snakemake: `cellhub.bash annotation execute --cores 1 -p --jobs 5 --executor slurm --default-resources mem_mb=16000`
if [[ "$cmd" = "annotation" || "$cmd" = "cluster" ]]; then
# if [[ "$cmd" = "cluster" ]]; then
    snakemake -s $CELLHUB_HOME/cellhub/Snakefile --config target="$cmd" mode="$subcmd" "$@"
else
    if [[ "$subcmd" == "genconfig" ]]; then
        cellhub $cmd config
    elif [[ "$subcmd" == "execute" ]]; then
        cellhub $cmd make full "$@"
    else
        echo "Unknown sub-command: ${subcmd}."
        usage
    fi
fi
