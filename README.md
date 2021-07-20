# dmethr
Hunting DNA methylation repression of transcriptional activity and its reactivation in a demethylating context

# Installation

```
devtools::install_github("fchuffar/dmethr")
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

## knit meta_analysis on laptop

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

