# run models
tcga_projects = c("TCGA-LUSC", "TCGA-LUAD")

nb_rnd_feat = 0
ud_strs = c(2500, 1000, 500, 250)
# feature_pretreatments = c("cen", "raw")
feature_pretreatments = c("raw")
reducer_func2_names = c("mean", "max")

# Compute 00_dmethr_pipeline_XXX.html and feats_XXX.xlsx files
for (feature_pretreatment in feature_pretreatments) {
  print(feature_pretreatment)
  for (ud_str in ud_strs) {
    print(ud_str)
    for (reducer_func2_name in reducer_func2_names) {
      prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
      rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=paste0("00_dmethr_pipeline_", prefix3, ".html"))
    }
  }
}


# Meta analysis start here
if (!exists("mread.xlsx")) {mread.xlsx = memoise::memoise(openxlsx::read.xlsx)}
    
stats = NULL
for (feature_pretreatment in feature_pretreatments) {
  print(feature_pretreatment)
  for (ud_str in ud_strs) {
    print(ud_str)
    for (reducer_func2_name in reducer_func2_names) {
      print(reducer_func2_name)
      
      prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
      print(prefix3)
      featsout = mread.xlsx (paste0("feats_", prefix3, ".xlsx"))
      
      
      tmp_list = list(
        feature_pretreatment = feature_pretreatment, 
        ud_str = ud_str, 
        reducer_func2_name = reducer_func2_name,
        nb_rnd_feat = nb_rnd_feat, 
        nb_cpg_rich = sum(featsout$cpg_status %in% "rich"),
        nb_methpp =   sum(featsout$methplusplus),
        NULL
      )

      for (gse in gses) {
        prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
        print(prefix3)
        feats = mread.xlsx (paste0("feats_", prefix3, ".xlsx"))
  
        gene_set = rownames(feats)[feats[["methplusplus"]]]
        expression_vector = feats[,paste0("l2fc_", gse)]
        names(expression_vector) = rownames(feats)
        expression_vector = expression_vector[!is.na(expression_vector)]
        utest = wilcox.test(expression_vector~ifelse(names(expression_vector)%in%gene_set, "methplusplus", "others"), las=2)
        tmp_list[[paste0("utest_pval_", gse)]] = utestutest$p.value
      }

      
      if (is.null(stats)) {
        stats = tmp_list
      } else{
        stats = rbind(stats, tmp_list)
      }
    }
  }
}
stats



plot(stats$ud_str, -log10(stats$upval), col=as.numeric(as.factor(stats$gse)), pch=ifelse(stats$reducer_func2_name=="mean", 1, 2))




beta_reg = c()

layout(matrix(1:3, 1), respect=TRUE)
for(ud_str in ud_strs){
  
    reducer_func2_name = "mean"
    feature_pretreatment = "raw"
    prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
    featsout_raw = openxlsx::read.xlsx (paste0("feats_", prefix3, ".xlsx"))
    feature_pretreatment = "cen"
    prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
    featsout_cen = openxlsx::read.xlsx (paste0("feats_", prefix3, ".xlsx"))
    #plot(jitter(featsout_raw$nb_probes), jitter(featsout_cen$nb_probes), main=ud_str)
    #abline(a=0, b=1)
    
    X = featsout_cen$nb_probes
    Y = featsout_raw$nb_probes
    reg = lm(X~Y)
    plot(X, Y, main=ud_str, xlab="raw", ylab="center")
    abline(reg, col = 2)
    abline(a=0, b=1, lty=2, col="grey")
    beta_reg = c(beta_reg, reg$coefficients[[2]])
  
}
beta_reg



