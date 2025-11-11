import os

PARAMS = config["annotation"]
SPECIES = PARAMS["species"]
ENSEMBL_RELEASE = PARAMS["ensembl_release"]
ENSEMBL_HOST = (
    ""
    if PARAMS["ensembl_host"] == "default"
    else f' --ensemblhost="{PARAMS["ensembl_host"]}"'
)

RSCRIPT_DIR = f"{workflow.basedir}/{os.pardir}/R/scripts/"
ANNOTATION_DIR = "annotation.dir"
API_DIR = "api/annotation"


rule full:
    input:
        f"{API_DIR}/ensembl/ensembl.gene_name.map.tsv.gz",
        f"{API_DIR}/ensembl/ensembl.to.entrez.tsv.gz",
        f"{API_DIR}/kegg/kegg.pathways.rds",


rule fetch_ensembl:
    output:
        gene_map=f"{ANNOTATION_DIR}/ensembl.gene_name.map.tsv.gz",
        to_entrez=f"{ANNOTATION_DIR}/ensembl.to.entrez.tsv.gz",
    log:
        f"{ANNOTATION_DIR}/fetch_ensembl.log",
    params:
        script=f"{RSCRIPT_DIR}/annotation_fetch_ensembl.R",
        outdir=ANNOTATION_DIR,
        species=SPECIES,
        ensembl_release=ENSEMBL_RELEASE,
        ensembl_host=ENSEMBL_HOST,
    threads: 1
    shell:
        """
        Rscript "{params.script}" \
            --ensemblversion="{params.ensembl_release}" \
            {params.ensembl_host} \
            --species="{params.species}" \
            --outdir="{params.outdir}" \
            &> "{log}"
        """


rule fetch_kegg:
    output:
        rds=f"{ANNOTATION_DIR}/kegg.pathways.rds",
    log:
        f"{ANNOTATION_DIR}/kegg.pathways.log",
    params:
        script=f"{RSCRIPT_DIR}/annotation_fetch_kegg.R",
        outdir=ANNOTATION_DIR,
        species=SPECIES,
    threads: 1
    shell:
        """
        Rscript "{params.script}" \
            --species="{params.species}" \
            --outfile="{params.outdir}/kegg.pathways.rds" \
            &> "{log}"
        """


rule register:
    input:
        ensembl1=f"{ANNOTATION_DIR}/ensembl.gene_name.map.tsv.gz",
        ensembl2=f"{ANNOTATION_DIR}/ensembl.to.entrez.tsv.gz",
        kegg=f"{ANNOTATION_DIR}/kegg.pathways.rds",
    output:
        ensembl1=f"{API_DIR}/ensembl/ensembl.gene_name.map.tsv.gz",
        ensembl2=f"{API_DIR}/ensembl/ensembl.to.entrez.tsv.gz",
        kegg=f"{API_DIR}/kegg/kegg.pathways.rds",
    params:
        api_dir=API_DIR,
    threads: 1
    shell:
        """
        ln -sf "$(realpath {input.ensembl1})" "{output.ensembl1}"
        ln -sf "$(realpath {input.ensembl2})" "{output.ensembl2}"
        ln -sf "$(realpath {input.kegg})" "{output.kegg}"
        """
