# `dmethr`: a 3-step pipeline to explore DNA methylation controlled gene expression from omics data

DNA methylation, the addition of a methyl group to the DNA molecule, is one of several epigenetic mechanisms that cells use to control gene expression. The most studied DNA methylation mark is that of a Cytosine preceding a Guanine (CpG dinucleotides). Since DNA methylation levels have been found lower in many cancers compared to most mammalian somatic tissues [Nordor 2017] [Zhang 2021], reseachers have tried to investigate the implication of DNA methylation as a key regulator of gene expression contributing to cancer development [Jones 2012] [Shi 2020].
Although methylation of CpG-rich regions (called CpG islands) has been associated with stable gene repression, the rules directing the relationship between DNA methylation and gene expression have not yet been fully established. In order to address this question, researchers need to explore -omics data of various origins, including genome-wide genomic features, methylomes (DNA methylation data) and transcriptomic data. Previous studies have established 2/3 of human genes are associated with low methylated CpG-rich promoter ( figure 1) region, regardless of their expression levels, a small subset of genes with high tissue specificity of expression and whose repression in most normal non-germline (somatic) tissues is associated with highly methylated CpG-rich promoter region in non-cancer.
The bio-informatic pipelines developed here aimed at identifying these genes and test their expression variation upon demethylation in experimental models and in cancer samples.
The first objective of this study is to classify genes according to the three following characteristics: i) CpG density of their promoter region, ii) level of methylation of the promoter region in normal tissues, iii) level of expression in a demethylating context. The aim is to identify a set of genes whose promoter region is CpG-rich and widely methylated in most tissues and whose expression is up-regulated in a demethylating context. The relationship between the expression of these genes and the methylation of their promoter region is then explored in the context of cancer. For this purpose
we are developing `dmethr`.



# Installation

```
devtools::install_github("fchuffar/dmethr")
```

# Conda env

```
source ~/conda_config.sh
# conda create -n dmethr_env
conda activate dmethr_env
# mamba install -c anaconda -c bioconda -c conda-forge -c r r-base libopenblas bioconductor-geoquery bioconductor-affy bioconductor-biobase r-seqinr r-rcpparmadillo r-devtools r-fastmap r-matrix r-kernsmooth r-catools r-gtools r-nortest r-survival r-beanplot r-gplots r-sass jquery r-nlme r-bslib r-sourcetools r-fontawesome r-xtable r-httpuv r-dbi r-igraph
# devtools::install_github("fchuffar/epimedtools")
# gse='GSE41037'; rmarkdown::render('01_build_study_generic.Rmd', output_file=paste0('01_build_study_',gse,'.html'));
# gse='GSE119617'; rmarkdown::render('01_build_study_generic.Rmd', output_file=paste0('01_build_study_',gse,'.html'));
# gse='GSE40279'; rmarkdown::render('01_build_study_generic.Rmd', output_file=paste0('01_build_study_',gse,'.html'));
# gse='GSE42861'; rmarkdown::render('01_build_study_generic.Rmd', output_file=paste0('01_build_study_',gse,'.html'));
# dim(s$data) # problem with GSE42861 : only 374449 probes. not the case with R3.6.1 
snakemake -k --cores 1 -s 00_custom_studies_wf.py -pn

```

# Run pipeline

```
rmarkdown::render("00_dmethr_pipeline.Rmd")
rmarkdown::render("01_tss_cpg_status.Rmd")
rmarkdown::render("02_mapreduce_mean.Rmd")
rmarkdown::render("03_da_GSE45332.Rmd")
rmarkdown::render("03_da_GSE5816.Rmd")
rmarkdown::render("03_da_GSE14315.Rmd")
rmarkdown::render("04_results.Rmd")
rmarkdown::render("00_dmethr_poster.Rmd")
rmarkdown::render("00_dmethr_poster.Rmd", output_format="pdf_document")
```

# Run pipeline on the cluster

## Cluster prerequisites

https://github.com/fchuffar/practicle_sessions/blob/master/ciment_infrastructure/config_cluster.sh

## Package requirements

```
# BSgenome.Hsapiens.UCSC.hg38
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

cd ~/projects
ln -s ~/../chuffarf/projects/tcga_studies
ln -s ~/../chuffarf/projects/genes
ln -s ~/../chuffarf/projects/study_gpl570
ln -s ~/../chuffarf/projects/all_human_tissues
ln -s ~/../chuffarf/projects/hadaca
```

## Launch pipeline

```
ssh dahu
chandler
oarsub --project epimed  -l /nodes=1,core=32,walltime=00:30:00 -I -t devel
oarsub --project epimed  -l /nodes=1,core=4,walltime=6:00:00 -I 
oarsub --project epimed -p "network_address='luke41'" -l /core=2,walltime=168:00:00 -I

mkdir -p ~/projects/
cd ~/projects/
git clone https://github.com/fchuffar/dmethr.git
cd ~/projects/dmethr/vignettes
git pull

R
# source("00_run_pipeline.R")
# rmarkdown::render("00_meta_analysis.Rmd")
```

## Knit meta_analysis on laptop

```
cd ~/projects/dmethr/vignettes/
rsync -auvP dahu:~/projects/dmethr/vignettes/*.xlsx ~/projects/dmethr/vignettes/.
# rsync -auvP dahu:~/../chuffarf/projects/dmethr/vignettes/*.xlsx ~/projects/dmethr/vignettes/.
R
# rmarkdown::render("00_meta_analysis.Rmd")


```


# Development

```
devtools::document(); devtools::install(); devtools::check(build_args="--no-build-vignettes")
```

