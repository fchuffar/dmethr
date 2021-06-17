# step = 100
# ud_str = 2500
# feattrmt = "raw"
# reducer = "mean"
# prefix3 = paste0(reducer, "_", feattrmt, "_", ud_str, "_", step)
# rmarkdown::render("01_tss_cpg_status.Rmd", output_file=paste0("01_tss_cpg_status_", prefix3, ".html"))
# rmarkdown::render("02_mapreduce_mean.Rmd", output_file=paste0("02_mapreduce_mean_", prefix3, ".html"))
# rmarkdown::render("03_da_GSE45332.Rmd"   , output_file=paste0("03_da_GSE45332_"   , prefix3, ".html"))
# rmarkdown::render("03_da_GSE5816.Rmd"    , output_file=paste0("03_da_GSE5816_"    , prefix3, ".html"))
# rmarkdown::render("04_results.Rmd"       , output_file=paste0("04_results_"       , prefix3, ".html"))


nb_rnd_feat = 2000
feature_pretreatment = "cen"
prefix3 = paste0(feature_pretreatment, "_", nb_rnd_feat)
rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))

nb_rnd_feat = 2000
feature_pretreatment = "raw"
prefix3 = "centered"
prefix3 = paste0(feature_pretreatment, "_", nb_rnd_feat)
rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))


ud_strs = c(2500, 1000, 500)
feature_pretreatments = c("cen", "raw")



for (ud_str in ud_strs) {
  for (ud_str in ud_strs) {
    nb_rnd_feat = 2000
    feature_pretreatment = "raw"
    prefix3 = "centered"
    prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat)
    rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))
  }
}

# rmarkdown::render("01_tss_cpg_status.Rmd")
# rmarkdown::render("02_mapreduce_mean.Rmd")
# rmarkdown::render("03_da_GSE45332.Rmd"   )
# rmarkdown::render("03_da_GSE5816.Rmd"    )
# rmarkdown::render("04_results.Rmd"       )



