stopifnot(
    require(devtools),
    require(BiocManager))


install.package <- function(x, source="CRAN", github_repo="")
{
    pkg <- x

    if (!require(pkg, character.only = TRUE))
    {
        message("installing missing package: ", pkg)

        if(source=="CRAN")
            {
                install.packages(pkg, dep=TRUE)

            } else if (source=="bioconductor")
            {
                BiocManager::install(pkg)

            } else if (source=="github")
            {
                install_github(github_repo)
            }

        if(!require(pkg,character.only = TRUE)) stop("Package not found")

    }
}

cran_packages <- c("Cairo",
                   "circlize",
                   "clustree",
                   "colormap",
                   "cowplot",
                   "data.table",
                   "docopt",
                   "dplyr",
                   "futile.logger",
                   "ggExtra",
                   "ggplot2",
                   "ggrepel",
                   "ggwordcloud",
                   "gplots",
                   "grDevices",
                   "grid",
                   "gridExtra",
                   "knitr",
                   "magrittr",
                   "Matrix",
                   "methods",
                   "optparse",
                   "openxlsx",
                   "pheatmap",
                   "RColorBrewer",
                   "reshape2",
                   "R.utils",
                   "S4Vectors",
                   "scales",
                   "stringr",
                   "tibble",
                   "tidyr",
                   "tidyverse",
                   "yaml")

bioconductor_packages <- c("ComplexHeatmap",
                           "DropletUtils",
                           "scater",
			   "SingleR",
			   "celldex",
               "ensembldb")

github_packages <- c("loomR"="mojaveazure/loomR",
                     "gsfisher"="sansomlab/gsfisher")


message("installing cran packages")
for(x in cran_packages)
{
    install.package(x, source="CRAN")
}

message("installing bioconductor packages")
for(x in bioconductor_packages)
{
    install.package(x, source="bioconductor")
}

message("installing github packages")
for(x in names(github_packages))
{
    install.package(x, source="github", github_repo=github_packages[x])
}
