# dmethr
Hunting DNA methylation repression of transcriptional activity and its reactivation in a demethylating context



# Run pipeline on the cluster

https://github.com/fchuffar/practicle_sessions/blob/master/ciment_infrastructure/config_cluster.sh

```
ssh dahu
chandler
oarsub --project epimed  -l /nodes=1,core=32,walltime=00:30:00 -I -t devel
oarsub --project epimed  -l /nodes=1,core=4,walltime=6:00:00 -I 
oarsub --project epimed -p "network_address='luke41'" -l /core=2,walltime=168:00:00 -I

mkdir -p ~/projects/
cd ~/projects/
git clone https://github.com/fchuffar/dmethr.git
cd ~/projects/dmethr/results
git pull

R
# source("00_run_pipeline.R")
# rmarkdown::render("meta_analysis.Rmd")




```

# Run pipeline on the cluster

```
# knit meta_analysis on laptop
cd ~/projeccd ~/projects/dmethr/results/
rsync -auvP dahu:~/projects/dmethr/results/*.xlsx ~/projects/dmethr/results/.
# rsync -auvP dahu:~/../chuffarf/projects/dmethr/results/*.xlsx ~/projects/dmethr/results/.
R
# rmarkdown::render("meta_analysis.Rmd")


```


# Installation

## Requirements

```
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

cd ~/projects
ln -s ~/../chuffarf/projects/tcga_studies
ln -s ~/../chuffarf/projects/genes

```

