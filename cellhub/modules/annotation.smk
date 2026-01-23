import os

from cellhub.parser.annotation_setup import AnnotationSetup


print(config)
annot = AnnotationSetup(config)


RSCRIPT_DIR = f"{workflow.basedir}/{os.pardir}/R/scripts/"


rule full:
    input:
        os.path.join(annot.api_dir, "ensembl", "ensembl.gene_name.map.tsv.gz"),
        os.path.join(annot.api_dir, "ensembl", "ensembl.to.entrez.tsv.gz"),
        os.path.join(annot.api_dir, "kegg", "kegg.pathways.rds"),


rule fetch_ensembl:
    output:
        ensembl1=os.path.join(annot.out_dir, "ensembl.gene_name.map.tsv.gz"),
        ensembl2=os.path.join(annot.out_dir, "ensembl.to.entrez.tsv.gz"),
    log:
        os.path.join(annot.out_dir, "fetch_ensembl.log"),
    params:
        script=f"{RSCRIPT_DIR}/annotation_fetch_ensembl.R",
        outdir=annot.out_dir,
        species=annot.species,
        ensembl_release=annot.ensembl_release,
        ensembl_host=annot.ensembl_host,
    threads: annot.threads
    resources:
        mem_mb=annot.mem_mb,
        time=annot.time,
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
        kegg=os.path.join(annot.out_dir, "kegg.pathways.rds"),
    log:
        os.path.join(annot.out_dir, "kegg.pathways.log"),
    params:
        script=f"{RSCRIPT_DIR}/annotation_fetch_kegg.R",
        outdir=annot.out_dir,
        species=annot.species,
    threads: annot.threads
    resources:
        mem_mb=annot.mem_mb,
        time=annot.time,
    shell:
        """
        Rscript "{params.script}" \
            --species="{params.species}" \
            --outfile="{params.outdir}/kegg.pathways.rds" \
            &> "{log}"
        """


rule register:
    input:
        ensembl1=rules.fetch_ensembl.output.ensembl1,
        ensembl2=rules.fetch_ensembl.output.ensembl2,
        kegg=rules.fetch_kegg.output.kegg,
    output:
        ensembl1=os.path.join(annot.api_dir, "ensembl", "ensembl.gene_name.map.tsv.gz"),
        ensembl2=os.path.join(annot.api_dir, "ensembl", "ensembl.to.entrez.tsv.gz"),
        kegg=os.path.join(annot.api_dir, "kegg", "kegg.pathways.rds"),
    threads: annot.threads
    resources:
        mem_mb=annot.mem_mb,
        time=annot.time,
    shell:
        """
        ln -sf "$(realpath {input.ensembl1})" "{output.ensembl1}"
        ln -sf "$(realpath {input.ensembl2})" "{output.ensembl2}"
        ln -sf "$(realpath {input.kegg})" "{output.kegg}"
        """
