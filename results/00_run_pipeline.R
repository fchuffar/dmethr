source("common.R")

# run models
tcga_projects = c("TCGA-LUSC")
#, "TCGA-LUAD")

nb_rnd_feat = 0
ud_strs = c(250, 500, 1000, 2500)
# feature_pretreatments = c("raw", "cen")
# feature_pretreatments = c("cen")
feature_pretreatments = c("raw")
reducer_func2_names = c("max", "mean", "min")

gses = c(
  "GSE45332", 
  "GSE5816", 
  "GSE14315"#, 
  # "GSE25427", 
  # "GSE22250"
  )
models = c(
  GSE45332 = "DNMT DKO vs. WT (ref.)", 
  GSE5816  = "5-aza 1000 nM vs. DMSO (ref.)",
  GSE14315 = "5-aza vs. Ctrl (ref.)"#,
  #GSE25427 = "5-aza vs. Mock (ref.)",
  #GSE22250 = "5-aza vs. WT (ref.)"
)

source("common.R")
# Compute 00_dmethr_pipeline_XXX.html and feats_XXX.xlsx files
for (feature_pretreatment in feature_pretreatments) {
  for (ud_str in ud_strs) {
    for (reducer_func2_name in reducer_func2_names) {
      prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
      feats_filename = paste0("feats_", prefix3, ".xlsx")
      feats_filename = paste0("feats_", prefix3, ".xlsx")
      html_filename = paste0("00_dmethr_pipeline_", prefix3, ".html")
      if (!file.exists(feats_filename)) {
        print(paste0(feats_filename))
        rmarkdown::render("00_dmethr_pipeline.Rmd", output_file=html_filename)
      }
    }
  }
}

# stop("EFN")


# Meta analysis start here
    
stats = NULL
for (feature_pretreatment in feature_pretreatments) {
  print(feature_pretreatment)
  for (ud_str in ud_strs) {
    print(ud_str)
    for (reducer_func2_name in reducer_func2_names) {
      print(reducer_func2_name)
      
      prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
      feats_filename = paste0("feats_", prefix3, ".xlsx")
      print(feats_filename)
      feats = mread.xlsx (feats_filename)
        
      tmp_list = list(
        feature_pretreatment = feature_pretreatment, 
        ud_str = ud_str, 
        reducer_func2_name = reducer_func2_name,
        nb_rnd_feat = nb_rnd_feat, 
        nb_cpg_rich = sum(feats$cpg_status %in% "rich"),
        nb_methpp =   sum(feats$methplusplus)
      )

      for (gse in gses) {
        print(gse)
        gene_set = rownames(feats)[feats[["methplusplus"]]]
        expression_vector = feats[,paste0("l2fc_", gse)]
        names(expression_vector) = rownames(feats)
        expression_vector = expression_vector[!is.na(expression_vector)]
        utest = wilcox.test(expression_vector~ifelse(names(expression_vector)%in%gene_set, "methplusplus", "others"), las=2)
        pv = utest$p.value
        print(pv)
        tmp_list[[paste0("utest_pval_", gse)]] = pv
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
stats[,1]
stats = data.frame(lapply(data.frame(stats, stringsAsFactors=FALSE), unlist), stringsAsFactors=FALSE)
stats[,1]




layout(matrix(1:2,1), respect = TRUE)
plot(0,0,col=0, xlim=range(stats$ud_str), ylim=c(0,50), xlab="ud_str", ylab="-log10(pval_utest)")
for (i in 1:length(gses)) {
  gse = gses[i]
  points(stats$ud_str, -log10(stats[[paste0("utest_pval_",gse)]]), 
       col=i, pch=ifelse(stats[["reducer_func2_name"]]=="mean", 1, 2))
  for (reducer_func2_name in reducer_func2_names) {
    lines(stats[stats$reducer_func2_name==reducer_func2_name,]$ud_str, -log10(stats[stats$reducer_func2_name==reducer_func2_name, paste0("utest_pval_",gse)]), 
           col=i, pch=ifelse(stats[["reducer_func2_name"]]=="mean", 1, 2))
  }
} 
plot(stats$ud_str, stats$nb_methpp, pch=ifelse(stats[["reducer_func2_name"]]=="mean", 1, 2), xlab = "ud_str", ylab = "nbgenes methplusplus")
for (reducer_func2_name in reducer_func2_names) {
    lines(stats[stats$reducer_func2_name==reducer_func2_name,]$ud_str, stats[stats$reducer_func2_name==reducer_func2_name,]$nb_methpp)
}






# layout(matrix(1:4, 4), respect=TRUE, byrow=TRUE)
# beta_reg = c()
# for(ud_str in ud_strs){
#     reducer_func2_name = "mean"
#     feature_pretreatment = "raw"
#     prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
#     feats_raw = openxlsx::read.xlsx (paste0("feats_", prefix3, ".xlsx"))
#     feature_pretreatment = "cen"
#     prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
#     feats_cen = openxlsx::read.xlsx (paste0("feats_", prefix3, ".xlsx"))
#     #plot(jitter(feats_raw$nb_probes), jitter(feats_cen$nb_probes), main=ud_str)
#     #abline(a=0, b=1)
# 
#     X = feats_cen$nb_probes
#     Y = feats_raw$nb_probes
#     reg = lm(X~Y)
#     plot(X, Y, main=ud_str, xlab="raw", ylab="center")
#     abline(reg, col = 2)
#     abline(a=0, b=1, lty=2, col="grey")
#     beta_reg = c(beta_reg, reg$coefficients[[2]])
# }
# beta_reg


