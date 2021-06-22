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



# run models
tcga_project = "TCGA-LUSC"

nb_rnd_feat = 0
ud_strs = c(2500, 1000, 500)
feature_pretreatments = c("cen", "raw")
reducer_func2_names = c("mean", "max")

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

stop("ici")

# meta analisis strat here

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
        feature_pretreatment=feature_pretreatment, 
        ud_str=ud_str, 
        reducer_func2_name= featsout$cp,
        nb_rnd_feat=nb_rnd_feat, 
        nb_cpg_rich=sum(featsout$cpg_status %in% "rich"),
        nb_methpp=sum(featsout$methplusplus),
        NULL
      )
      if (is.null(stats)) {
        stats = tmp_list
      } else{
        stats = rbind(stats, tmp_list)
      }
    }
  }
}
stats













gses = c("GSE5816", "GSE45332")

for (gse in gses) {
  #plot(0,0, col=0, xlim=c(0,100); ylim=c(0,10))
  for (feature_pretreatment in feature_pretreatments) {
    print(feature_pretreatment)
    for (ud_str in ud_strs) {
      print(ud_str)
      for(reducer_func2_name in reducer_func2_names) {
        prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
        featsout = mread.xlsx (paste0("feats_", prefix3, ".xlsx"))
        rownames(featsout) = featsout[,4]
        gene_symbols = rownames(featsout[featsout$methplusplus,])[order(featsout[featsout$methplusplus,][[paste0("l2fc_", gse)]], decreasing=TRUE)][1:100]
        pvals = sapply(gene_symbols, function(gene_symbol) {
          print(gene_symbol)
          mdata = try(mget_multiomic_data(gene_symbols=gene_symbol, tcga_project=tcga_project)) 
          if (class(mdata) != "try-error"){
            if ((!is.null(mdata$d)) & (gene_symbol %in% colnames(mdata$d))) {
              mdata$d 
              head(mdata$d)
              Y=mdata$d[,gene_symbol]
              X=apply(mdata$d[,mdata$probes],1,mean)
              m= lm(X~Y)
              sm=summary(m)
              fstat=sm$fstatistic
              pval = pf(fstat[[1]], fstat[[2]], fstat[[3]], lower.tail = FALSE)
              return(pval)
            } else {
              warning(paste0("probleme with ", gene_symbol))
              return(NULL)
            }
          } else {
            warning(paste0("probleme with ", gene_symbol))
            return(NULL)
          }
        })  
        #plot(density(unlist(-log10(pvals)), na.rm=TRUE))
        #lines(density(-log10(featsout$pval_omic), na.rm=TRUE), col=which(ud_str==ud_strs))
      }
    }
  }
}















pvals =c()

for (gene_symbol in gene_symbols[1:length(gene_symbols)]) {
  print(gene_symbol)
  mdata = try(mget_multiomic_data(gene_symbols=gene_symbol, tcga_project=tcga_project)) 
  if (class(mdata)!= "try-error"){
    if (!is.null(mdata$d)) {
      mdata$d 
      head(mdata$d)
      Y=mdata$d[,gene_symbol]
      X=apply(mdata$d[,mdata$probes],1,mean)
      m= lm(X~Y)
      sm=summary(m)
      fstat=sm$fstatistic
      pval = pf(fstat[[1]], fstat[[2]], fstat[[3]], lower.tail = FALSE)
      pvals = c(pvals, pval)
    } else {
      warning(paste0("probleme with ", gene_symbol))
    }
  }
  #try(momic_pattern(gene_symbol, tcga_project))
}
layout(1, respect=TRUE)
plot(density(-log10(pvals)))

















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

# rmarkdown::render("01_tss_cpg_status.Rmd")
# rmarkdown::render("02_mapreduce_mean.Rmd")
# rmarkdown::render("03_da_GSE45332.Rmd"   )
# rmarkdown::render("03_da_GSE5816.Rmd"    )
# rmarkdown::render("04_results.Rmd"       )




#plot(0,0, col=0, xlim=???; xlim=???)
i=1
for (feature_pretreatment in feature_pretreatments) {
  print(feature_pretreatment)
  for (ud_str in ud_strs) {
    print(ud_str)
    for(reducer_func2_name in reducer_func2_names) {
      prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
      featsout = mread.xlsx (paste0("feats_", prefix3, ".xlsx"))
      #lines(density(-log10(featsout$pval_omic), na.rm=TRUE), col=i)
      lines(density(-log10(featsout$pval_omic), na.rm=TRUE), col=which(ud_str==ud_strs))
      i = i+1
    }
  }
}

#plot(0,0, col=0, xlim=???; xlim=???)
i=1
for (feature_pretreatment in feature_pretreatments) {
  print(feature_pretreatment)
  for (ud_str in ud_strs) {
    print(ud_str)
    for(reducer_func2_name in reducer_func2_names) {
      prefix3 = paste0(feature_pretreatment, "_", ud_str, "_", nb_rnd_feat, "_", reducer_func2_name)
      featsout = mread.xlsx (paste0("feats_", prefix3, ".xlsx"))
      #lines(density(-log10(featsout$pval_omic), na.rm=TRUE), col=i)
      lines(density(-log10(featsout$pval_omic), na.rm=TRUE), col=which(feature_pretreatment==feature_pretreatments))
      i = i+1
    }
  }
}




corefeatsout_raw_2500_0_mean = openxlsx::read.xlsx ("corefeats_raw_2500_0_mean.xlsx")
corefeatsout_raw_1000_0_mean = openxlsx::read.xlsx ("corefeats_raw_1000_0_mean.xlsx")
corefeatsout_raw_500_0_mean = openxlsx::read.xlsx ("corefeats_raw_500_0_mean.xlsx")
corefeatsout_raw_2500_0_max = openxlsx::read.xlsx ("corefeats_raw_2500_0_max.xlsx")
corefeatsout_raw_1000_0_max = openxlsx::read.xlsx ("corefeats_raw_1000_0_max.xlsx")
corefeatsout_raw_500_0_max = openxlsx::read.xlsx ("corefeats_raw_500_0_max.xlsx")
corefeatsout_cen_2500_0_mean = openxlsx::read.xlsx ("corefeats_cen_2500_0_mean.xlsx")
corefeatsout_cen_1000_0_mean = openxlsx::read.xlsx ("corefeats_cen_1000_0_mean.xlsx")
corefeatsout_cen_500_0_mean = openxlsx::read.xlsx ("corefeats_cen_500_0_mean.xlsx")
corefeatsout_cen_2500_0_max = openxlsx::read.xlsx ("corefeats_cen_2500_0_max.xlsx")
corefeatsout_cen_1000_0_max = openxlsx::read.xlsx ("corefeats_cen_1000_0_max.xlsx")
corefeatsout_cen_500_0_max = openxlsx::read.xlsx ("corefeats_cen_500_0_max.xlsx")
# 
# 
 plot(density(-log10(corefeatsout_raw_2500_0_mean$pvals), col=1))
lines(density(-log10(corefeatsout_raw_1000_0_mean$pvals)), col=corefeatsout_raw_1000_0_mean$pvals+1)
lines(density(-log10(corefeatsout_raw_500_0_mean$pvals)), col=corefeatsout_raw_1000_0_mean$pvals+2)
lines(density(-log10(corefeatsout_raw_2500_0_max$pvals)), col=corefeatsout_raw_2500_0_max$pvals+3)
lines(density(-log10(corefeatsout_raw_1000_0_max$pvals)), col=corefeatsout_raw_1000_0_max$pvals+4)
lines(density(-log10(corefeatsout_raw_500_0_max$pvals)), col=corefeatsout_raw_500_0_max$pvals+5)
# lines(density(-log10(corefeatsout_cen_2500_0_mean$pvals)), col=corefeatsout_cen_2500_0_mean$pvals+6)
# lines(density(-log10(corefeatsout_cen_1000_0_mean$pvals)), col=corefeatsout_cen_1000_0_mean$pvals+7)
# lines(density(-log10(corefeatsout_cen_500_0_mean$pvals)), col=corefeatsout_cen_500_0_mean$pvals+8)
# lines(density(-log10(corefeatsout_cen_2500_0_max$pvals)), col=corefeatsout_cen_2500_0_max$pvals+9)
# lines(density(-log10(corefeatsout_cen_1000_0_max$pvals)), col=corefeatsout_cen_1000_0_max$pvals+10)
# lines(density(-log10(corefeatsout_cen_500_0_max$pvals)), col=corefeatsout_cen_500_0_max$pvals+11)


# pvals
# 
# for (ud_str in ud_strs){
#   for(feature_pretreatment in feature_pretreatments){
#     for(reducer_func2_name in reducer_func2_names){
#       
# }
# }
# }