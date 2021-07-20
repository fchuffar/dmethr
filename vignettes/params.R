# run models
tcga_projects = c("TCGA-LUSC")
#, "TCGA-LUAD")

nb_rnd_feat = 0
ud_strs = c(250, 500, 1000, 2500)
# feature_pretreatments = c("raw", "cen")
# feature_pretreatments = c("cen")
feature_pretreatments = c("raw")

reducer_func2_names = c("mean", "max", "min")
reducer_func2_names = c("mean", "max")

models = c(
  GSE45332 = "DNMT DKO vs. WT (ref.)", 
  GSE14315 = "5-aza vs. Ctrl (ref.)",
  GSE5816  = "5-aza 1000 nM vs. DMSO (ref.)"#,
  #GSE25427 = "5-aza vs. Mock (ref.)",
  #GSE22250 = "5-aza vs. WT (ref.)"
)
gses = names(models)


