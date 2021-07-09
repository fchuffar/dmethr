# dmethr
Hunting DNA methylation repression of transcriptional activity and its reactivation in a demethylating context



# Launch on the cluster

https://github.com/fchuffar/practicle_sessions/blob/master/ciment_infrastructure/config_cluster.sh

```
ssh dahu
oarsub --project epimed  -l /nodes=1,core=32,walltime=00:30:00 -I -t devel
oarsub --project epimed  -l /nodes=1,core=4,walltime=6:00:00 -I 

mkdir -p ~/projects/
cd ~/projects/
git clone https://github.com/fchuffar/dmethr.git
cd ~/projects/dmethr/results

R
source("00_run_pipeline.R")
```


# Installation

##Requirements

```
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

