FROM bioconductor/bioconductor_docker:RELEASE_3_20

RUN R -e 'BiocManager::install("GenomicRanges")'
RUN R -e 'install.packages("devtools")'
RUN R -e 'devtools::install_github("dzhang32/ggtranscript")'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("Polychrome")'
RUN R -e 'install.packages("data.table")'
RUN R -e 'install.packages("R.utils")'
RUN R -e 'install.packages("cowplot")'
RUN R -e 'install.packages("svglite")'