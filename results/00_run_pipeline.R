# stp = 100
# ud_str = 2500
# feattrmt = "raw"
# reducer = "mean"
# prefix3 = paste0(reducer, "_", feattrmt, "_", ud_str, "_", step)
# rmarkdown::render("01_tss_cpg_status.Rmd", output_file=paste0("01_tss_cpg_status_", prefix3, ".html"))
# rmarkdown::render("02_mapreduce_mean.Rmd", output_file=paste0("02_mapreduce_mean_", prefix3, ".html"))
# rmarkdown::render("03_da_GSE45332.Rmd"   , output_file=paste0("03_da_GSE45332_"   , prefix3, ".html"))
# rmarkdown::render("03_da_GSE5816.Rmd"    , output_file=paste0("03_da_GSE5816_"    , prefix3, ".html"))
# rmarkdown::render("04_results.Rmd"       , output_file=paste0("04_results_"       , prefix3, ".html"))


# nb_rnd_feat = 500
# feature_pretreatment = "cen"
# prefix3 = paste0(feature_pretreatment, "_", nb_rnd_feat)
# rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))



# 
# nb_rnd_feat = 2000
# feature_pretreatment = "raw"
# prefix3 = "centered"
# prefix3 = paste0(feature_pretreatment, "_", nb_rnd_feat)
# rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))
# 
nb_rnd_feat = 0
ud_strs = c(2500, 1000, 500)
feature_pretreatments = c("cen", "raw")
reducer_func2_name = "max"

stats = NULL
for (feature_pretreatment in feature_pretreatments) {
  print(feature_pretreatment)
  for (ud_str in ud_strs) {
    prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
    rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))
    tmp_list = list(
      feature_pretreatment=feature_pretreatment, 
      ud_str=ud_str, 
      nb_rnd_feat=nb_rnd_feat, 
      nb_cpg_rich=length(idx_rich),
      nb_methpp=length(gs_methplusplus),
      nb_genes=length(corefeats), 
      #pval_surv=1,
      #pval_momik=1,
      NULL
    )
    if (is.null(stats)) {
      stats = tmp_list
    } else{
      stats = rbind(stats, tmp_list)
    }
  }
}

stats

#beta_reg = c()

layout(matrix(1:3, 1), respect=TRUE)
for(ud_str in ud_strs){
  feature_pretreatment = "raw"
  prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
  featsout_raw = openxlsx::read.xlsx (paste0("feats_", prefix3, ".xlsx"))
  feature_pretreatment = "cen"
  prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
  featsout_cen = openxlsx::read.xlsx (paste0("feats_", prefix3, ".xlsx"))
  plot(jitter(featsout_raw$nb_probes), jitter(featsout_cen$nb_probes), main=ud_str)
  abline(a=0, b=1)
  
  X = featsout_cen$nb_probes
  Y = featsout_raw$nb_probes
  reg = lm(X~Y)
  plot(X, Y, main=ud_str, xlab="raw", ylab="center")
  abline(reg)
  beta_reg = c(beta_reg, reg$coefficients[[2]])
}

beta_reg

# rmarkdown::render("01_tss_cpg_status.Rmd")
# rmarkdown::render("02_mapreduce_mean.Rmd")
# rmarkdown::render("03_da_GSE45332.Rmd"   )
# rmarkdown::render("03_da_GSE5816.Rmd"    )
# rmarkdown::render("04_results.Rmd"       )



