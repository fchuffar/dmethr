conda create --clone epimedtools_env -n dmethr_env 
conda activate dmethr_env

conda install -c bioconda r-openxlsx r-writexls
conda install -c conda-forge r-mclust

# conda install -c conda-forge umap-learn r-umap r-caret r-randomforest
# conda install -c bioconda bioconductor-champ snakemake bedtools bioconductor-annotatr





# under R
# devtools::install_github("fchuffar/dnamaging")


# conda install  -c anaconda libopenblas
#
# conda install -c r \
#   r-base=3.6.1     \
#   r-rmarkdown      \
#   r-devtools       \
#   r r-xml
#
#   install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.tar.gz")
#
#
# # conda install --no-update-deps -c r r-xml
#
# conda install -c bioconda \
#   bioconductor-Biobase    \
#   bioconductor-affy       \
#   bioconductor-GEOquery
#
# # conda install --no-update-deps -c bioconda            \
# #   bioconductor-illuminahumanmethylation27kanno.ilmn12.hg19    \
# #   bioconductor-illuminahumanmethylation450kanno.ilmn12.hg19   \
# #   bioconductor-illuminahumanmethylationepicanno.ilm10b4.hg19
#
#
# #conda install -c r r-locfit
# conda install -c r r-glmnet
#
#
#
# conda install -c bioconda bioconductor-rhdf5lib
# conda install -c bioconda bioconductor-rhdf5
# conda install -c bioconda bioconductor-hdf5array
# conda install -c bioconda snakemake
#
#
# BiocManager::install("HDF5Array")
# BiocManager::install("minfi")
# BiocManager::install("IlluminaHumanMethylation27kanno.ilmn12.hg19")
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
# BiocManager::install("EpiDISH")
#
#
#
#
# # https://github.com/fchuffar/epimedtools
# Sys.setenv(TAR = "/bin/tar") ; devtools::install_github("fchuffar/epimedtools")
# Sys.setenv(TAR = "/bin/tar") ; devtools::install_github("fchuffar/dnamaging")

# BiocManager::install("impute")


# install.packages("dbplyr")
# install.packages("tidyr")
# install.packages("memoise")
# install.packages("glmnet")
# install.packages("glmnetUtils")
# install.packages("beeswarm")
# install.packages("openxlsx")

#
#
#
# # bioconductor-methylclockdata                                \
# # bioconductor-methylclock                                    \
#
#
#
#
