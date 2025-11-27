import os

TARGET = config["target"]


rule full:
    input:
        f"{workflow.basedir}/yaml/config_{TARGET}.yml",
    output:
        f"config_{TARGET}.yml",
    shell:
        "cp {input} {output}"
