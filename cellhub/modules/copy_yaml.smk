import os

TARGET = config["target"]


rule target:
    input:
        f"{workflow.basedir}/yaml/pipeline_{TARGET}.yml",
    output:
        f"pipeline_{TARGET}.yml",
    shell:
        "cp {input} {output}"
