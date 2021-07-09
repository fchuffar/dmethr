if (!exists("mread.xlsx")) {mread.xlsx = memoise::memoise(openxlsx::read.xlsx)}


fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)
  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


expression_in_gtex = function(gene_symbols) {
  
  study_gtex_filename = "~/projects/gtex/results/study_gtex.rds"
  s = mreadRDS(study_gtex_filename)  
  if (length(gene_symbols)!=0) {
    tmp_idx_features = intersect(gene_symbols, rownames(s$data))

    # ```{r label="sort data by sample"}
    tissues = rev(names(sort(table(s$exp_grp$tissue_group_level1))))
    s$exp_grp = s$exp_grp[s$exp_grp$tissue_group_level1 %in% tissues,]
    s$exp_grp$tissue_group_level1 = factor(s$exp_grp$tissue_group_level1, levels=tissues)
    idx_samples = rownames(s$exp_grp)[order(s$exp_grp$tissue_group_level1)]
    s$data = s$data[,idx_samples]
    s$exp_grp = s$exp_grp[idx_samples,]
    # ```
    # ```{r label="normalization"}
    data = normalization(s$data[tmp_idx_features,], "zscore_rows")
    data = data[!apply(is.na(data), 1, any),]
    # ```
    # ```{r label="clustering each tissue group"}
    outputs = lapply(tissues, function (tissue) {
      # tissue = "prostate"
      print(tissue)
      idx_sample = rownames(s$exp_grp)[s$exp_grp$tissue_group_level1 == tissue]
      d = data[,idx_sample]
      sum(is.na(d))
      d = d[apply(d, 1, any),]
      foo = plot_expr_hm(
        data=d                           ,
        rsc=NULL                         , 
        csc=NULL                         , 
        nb_grp_row=4                     ,
        nb_grp_col=4                     , 
        main=tissue                      , 
        hcmeth_cols="eucl_dist"          , 
        hcmeth_rows="eucl_dist"          , 
        normalization=FALSE              , 
        ordering_func=median             , 
        colors=c("cyan", "black", "red") , 
        PCA=FALSE                        ,
        PLOT_MAIN_HM=FALSE
      )
      dendro = as.dendrogram(foo$hc_col)  
      sample = foo$hc_col$labels[foo$hc_col$order][ceiling(length(foo$hc_col$order)/2)]
      tissue = tissue
      list(dendro=dendro, sample=sample, tissue=tissue)
    })
    dendro_cols = do.call(merge, lapply(outputs, "[[", "dendro"))
    labels = sapply(outputs, "[[", "tissue")
    names(labels) = sapply(outputs, "[[", "sample")
    # ```
    # ```{r}
    # data = s$data[tmp_idx_features,]

    colors = c("cyan", "black", "red")
    cols = colorRampPalette(colors)(20)
    main=""

    tmp_tab = table(s$exp_grp$tissue_group_level1)
    ColSideColors = unlist(apply(cbind(data.frame(tmp_tab), col=rep(c("white", "grey"), length(tmp_tab))[1:length(tmp_tab)]), 1, function (l) {
      rep(l[[3]], as.numeric(l[[2]]))
    }))

    # colnames by tissues
    tmp_data = data
    tmp_cn = colnames(tmp_data)
    names(tmp_cn) = tmp_cn
    tmp_cn[] = NA
    tmp_cn
    tmp_cn[names(labels)] = labels
    colnames(tmp_data) = tmp_cn
    # colnames(tmp_data)


    tmp_data[tmp_data < -4 ] = -4
    tmp_data[tmp_data > 4 ] = 4
    foo = plot_expr_hm(
      data=tmp_data                    ,
      rsc=NULL                         , 
      csc=ColSideColors                , 
      Colv=dendro_cols                  , 
      main=""                          , 
      hcmeth_cols=FALSE                , 
      hcmeth_rows="cor"                , 
      normalization=FALSE              , 
      ordering_func=median             , 
      colors=c("cyan", "black", "red") , 
      PCA=FALSE                        ,
      PLOT_MAIN_HM=TRUE                ,
      cexRow=0.6                       , 
      cexCol=0.6
    )

    # ```
    #
    # ```{r}

    Rowv = as.dendrogram(foo$hc_row)

    tmp_data = t(apply(data, 1, function(l) {
      # l = data [1,]
      # bp = boxplot(l~s$exp_grp[colnames(data),]$tissue_group_level1, las=2)
      m = lm(l~s$exp_grp[colnames(data),]$tissue_group_level1)
      tmp_coef =  m$coefficients

      oo = tmp_coef[1]
      tmp_coef = tmp_coef  + tmp_coef[1]
      tmp_coef[1] =  oo 
      names(tmp_coef) = c(levels(s$exp_grp$tissue_group_level1)[1], do.call(rbind, strsplit(names(tmp_coef)[-1], "tissue_group_level1"))[,2])
      # points(tmp_coef, col=2)
      tmp_coef
    }))

    RowSideColors = rep("white", nrow(data))

    tmp_data[tmp_data < -4 ] = -4
    tmp_data[tmp_data > 4 ] = 4
    foo = plot_expr_hm(
      data=tmp_data                    ,
      rsc=RowSideColors                , 
      csc=NULL                         , 
      Rowv=Rowv                        ,
      # Colv=dendro_cols                  ,
      main=""                          , 
      hcmeth_cols=FALSE                , 
      hcmeth_rows="cor"                , 
      normalization=FALSE              , 
      ordering_func=median             , 
      colors=c("cyan", "black", "red") , 
      PCA=FALSE                        ,
      PLOT_MAIN_HM=TRUE                ,
      cexRow=0.6                       , 
      cexCol=0.6
    )
  }
}

truncate_survival = function(s, censoring_time) {
  ## ending survival study at XX months
  censoring_time = 60
  s$exp_grp$dead = as.logical(s$exp_grp$dead)
  idx = which(s$exp_grp$os_months>censoring_time)
  if (length(idx) > 0) {
    s$exp_grp[idx,]$os_months = censoring_time         # replace by truncation value
    s$exp_grp[idx,]$dead = FALSE               # replace deaths by censors  
  }
  s$exp_grp$os = survival::Surv(s$exp_grp$os_months, s$exp_grp$dead)  
  return(s)
}

et_gsea_plot = function(expression_vector, gene_set, prefix, nperm=1000) {
  
  
  utest = wilcox.test(expression_vector~ifelse(names(expression_vector)%in%gene_set, "methplusplus", "others"), las=2)
  boxplot(rank(expression_vector)~ifelse(names(expression_vector)%in%gene_set, "methplusplus", "others"), main=paste0("Mann-Whitney U test (pval=", signif(utest$p.value, 3), ")"), ylab=paste0("rank(log2FoldChange)"), xlab="", col=adjustcolor(c("red", "grey"), alpha.f=.5))  
  if (i==1) {fig_label("B", cex=3)}
  return(utest)
  
  
  # den_vec_name = "DEN14TuvsDEN14NT"
  # n_top = 500
  # enrichement = "ENRICHED"
  # i = idx_kc[1]

  top = gene_set
  gs_filename = paste0(prefix, ".grp")
  write.table(top, gs_filename, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  gsea_input = data.frame(names(expression_vector), expression_vector)
  gsea_input = gsea_input[order(gsea_input[,2]),]
  gsea_input_filename = paste0("gsea_input.rnk")
  print(paste("gsea_input were exported in", gsea_input_filename, "file."))
  write.table(gsea_input, gsea_input_filename, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

  gsea_input = read.table(gsea_input_filename, row.names=1)
  stoda_input = gsea_input[,1]
  names(stoda_input) = rownames(gsea_input)
  # stoda::gseaplot(v=stoda_input, gs=top)

  dir.create("gsea_out", showWarnings=FALSE)
  outdir = paste0("gsea_out/", prefix, "_", i, "")
  unlink(outdir, recursive=TRUE)

  cmd = "/Applications/GSEA_4.1.0/gsea-cli.sh"
  args = paste0(" GSEAPreranked -gmx ", gs_filename, " -collapse No_Collapse -mode Max_probe -norm meandiv -nperm ", nperm, " -rnk ", gsea_input_filename, " -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 100000 -set_min 1 -zip_report false -out ", outdir, "")
  print(paste(cmd, args))
  system2(cmd, args)

  endplot_filename = paste0(outdir, "/", list.files(outdir)[1], "/enplot_", gs_filename, "_1.png")
  endplot_filename

  # grid::grid.raster(png::readPNG("gsea_out/DEN14TuvsDEN14NT_500_DEPLETED_TCGA-3K-AAZ8-01A/my_analysis.GseaPreranked.1623067518202/enplot_gs_DEN14TuvsDEN14NT_500_DEPLETED.grp_1.png"), x=1, y=1, width=1)

  addImg <- function(
    obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
    x = NULL, # mid x coordinate for image
    y = NULL, # mid y coordinate for image
    width = NULL, # width of image (in x coordinate units)
    interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
  ){
    if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
    USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
    PIN <- par()$pin # The current plot dimensions, (width, height), in inches
    DIM <- dim(obj) # number of x-y pixels for the image
    ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
    WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
    HEIi <- WIDi * ARp # height in inches
    HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
    rasterImage(image = obj, 
      xleft = x-(width/2), xright = x+(width/2),
      ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
      interpolate = interpolate)
  }
  par(mar=c(0, 0, 0, 0))
  plot.new()  
  par(mar=c(0, 0, 0, 0))
  # endplot_filename = "enplot_gs_DEN14TuvsDEN14NT_500_DEPLETED.grp_1.png"
  addImg(png::readPNG(endplot_filename), x=0.5,y=0.5,width=1)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  return(NULL)
}


normalization = function(data, normalization) {
  # normalization
  # colnames(data) = s$exp_grp[colnames(data),]$tissue_group_level1
  if (normalization=="zscore_rows" | normalization==TRUE) {
    data = data - apply(data, 1, mean)
    data = data / apply(data, 1, sd)    
  } else if (normalization=="zscore_cols") {
    data = t(data)
    data = data - apply(data, 1, mean)
    data = data / apply(data, 1, sd)    
    data = t(data)
  } else if (normalization=="rank_cols") {
    data = apply(data, 2, rank)
  } else if (normalization=="qqnorm_cols") {
    data = apply(data, 2, function(c) {
      qqnorm(c, plot.it=FALSE)$x
    })
  } else {
    # "do noting"
  }
  data
}


plot_expr_hm = function(data, 
  rsc=NULL                         , 
  csc=NULL                         , 
  Colv=NULL                        ,
  Rowv=NULL                        ,
  # hc_row=NULL,
  dendrogram=NULL                  ,
  nb_grp_row=4                     ,
  nb_grp_col=4                     , 
  main=""                          , 
  hcmeth_cols="eucl_dist"          , 
  hcmeth_rows="cor"                , 
  normalization=FALSE              , 
  ordering_func=median             , 
  colors=c("cyan", "black", "red") , 
  PCA=FALSE                        , 
  PLOT_MAIN_HM=TRUE                ,
  RowSideColors,
  ...  
) {
  # Remove rows with no variation (needed to clustering rows according to cor)
  # data = data[apply(data, 1, function(l) {length(unique(l))})>1, ]
  
  # normalization
  data = normalization(data, normalization)
  # RowSideColors = rsc
  # hc_row=NULL

  # clustering samples...
  if (is.null(Colv)) {
    if (hcmeth_cols != FALSE) {
      tmp_d = t(data)
      if (hcmeth_cols == "cor") {
        # ... based on correlation
        tmp_d = tmp_d[!apply(is.na(tmp_d), 1, any), ]
        d = dist(1 - cor(t(tmp_d), method="pe"))
        hc_col = hclust(d, method="complete")
        Colv = as.dendrogram(hc_col)
      } else if (hcmeth_cols == "ordered") {
        # ... ordered by median
        data = data[,order(apply(tmp_d, 1, ordering_func, na.rm=TRUE), decreasing=TRUE)]
        hc_col = Colv = NULL      
      } else {
        # ... based on eucl. dist.
        d = dist(tmp_d)
        hc_col = hclust(d, method="complete")
        Colv = as.dendrogram(hc_col)
      }
    } else {
      hc_col = Colv = NULL          
    }
    ColSideColors = rep("white", ncol(data))
    names(ColSideColors) = colnames(data)    
  } else {
    ColSideColors = csc    
    hc_col = Colv
  }

  # clustering features...
  if (is.null(Rowv)) {
    if (hcmeth_rows != FALSE) {
      tmp_d = data
      if (hcmeth_rows == "eucl_dist") {
        # ... based on eucl. dist.
        d = dist(tmp_d)
        hc_row = hclust(d, method="complete")
        Rowv = as.dendrogram(hc_row)
      } else if (hcmeth_rows == "ordered") {
        # ... ordered by median
        data = data[order(apply(tmp_d, 1, ordering_func, na.rm=TRUE), decreasing=TRUE),]
        hc_row = Rowv = NULL      
      } else {
        # ... bases on correlation
        tmp_d = tmp_d[!apply(is.na(tmp_d), 1, any), ]
        d = dist(1 - cor(t(tmp_d), method="pe"))
        hc_row = hclust(d, method="single")
        Rowv = as.dendrogram(hc_row)      
      }
    } else {
      hc_row = Rowv = NULL    
    }
    RowSideColors = rep("white", nrow(data))
    names(RowSideColors) = rownames(data)
  } else {
    RowSideColors = rsc
    hc_row = Rowv
  }


  print(RowSideColors)

  if (!is.null(rsc)) {
    RowSideColors = rsc
  }



  if (!is.null(csc)) {
    ColSideColors = csc
  }


  
  if (PCA) {
    tmp_d = t(data)
    
    nb_clusters = c()
    scores = c()
    best_score = NULL
    for (nb_cluster in 2:10) {
      for (i in 1:20) {
        k = kmeans(tmp_d, centers=nb_cluster)

        com1 = k$cluster + 10000
        com2 = as.numeric(as.factor(names(k$cluster)))
        names(com1) = names(com2) = paste0("id", 1:length(com1))
        # score = igraph::compare(com1, com2, method="nmi")
        score = -igraph::compare(com1, com2, method="vi")
        score
        
        nb_clusters = c(nb_clusters, nb_cluster)
        scores = c(scores, score)

        if (is.null(best_score)) {
          best_k = k
          best_score = score
          best_nb_cluster = nb_cluster
        } else if (score > best_score) {
          best_k = k
          best_score = score
          best_nb_cluster = nb_cluster
        }
      }
    }
    
    k = best_k
    nb_cluster = best_nb_cluster
    score = best_score

    ColSideColors = palette(RColorBrewer::brewer.pal(n=8, "Dark2"))[k$cluster[colnames(data)]]    
    names(ColSideColors) = colnames(data)    

    if (!is.null(csc)) {
      ColSideColors = csc
    }

    # PCA on tissues
    pca = prcomp(tmp_d, scale=FALSE)
    PLOT_SAMPLE_LABELS = length(unique(rownames(pca$x))) < nrow(pca$x)
    if (PLOT_SAMPLE_LABELS) {
      sample_labels = t(sapply(unique(rownames(pca$x)), function(t) {        
        idx= which(rownames(pca$x)==t)
        if (length(idx)==1) {
          pca$x[idx,]          
        } else {
          apply(pca$x[idx,], 2, mean)                    
        }
      }))      
    }
    v = pca$sdev * pca$sdev
    p = v / sum(v) * 100
    layout(matrix(1:6,2, byrow=FALSE), respect=TRUE)
    barplot(p)
    i=3
    j=2
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    if (PLOT_SAMPLE_LABELS) text(sample_labels[,i], sample_labels[,j], rownames(sample_labels))
    i=1
    j=3
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    if (PLOT_SAMPLE_LABELS) text(sample_labels[,i], sample_labels[,j], rownames(sample_labels))
    i=1
    j=2
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    if (PLOT_SAMPLE_LABELS) text(sample_labels[,i], sample_labels[,j], rownames(sample_labels))
    i=4
    j=5
    plot(pca$x[,i], pca$x[,j], xlab=paste0("PC", i, "(", signif(p[i], 3), "%)"), ylab=paste0("PC", j, "(", signif(p[j], 3), "%)"), col=ColSideColors[rownames(pca$x)])
    if (PLOT_SAMPLE_LABELS) text(sample_labels[,i], sample_labels[,j], rownames(sample_labels))

    plot(jitter(nb_clusters), scores)
    points(nb_cluster, score, col=2)

    # stop("EFN")    
  # } else {
  #   ColSideColors = rep("white", ncol(data))
  #   names(RowSideColors) = colnames(data)
  #   hc_col = Colv = NULL


  ColSideColors = palette(RColorBrewer::brewer.pal(n=length(unique(colnames(data))), "Dark2"))[as.factor(colnames(data))]    
  names(ColSideColors) = colnames(data)    

  }

  # stop("EFN")







  # if (is.null(rsc)) {
  #   grps = list()
  #   ct = cutree(hc_row, nb_grp_row)
  #   for (i in 1:nb_grp_row) {
  #     grps[[palette()[i]]] = names(ct)[ct==i]
  #   }
  #   # print(grps)
  #   RowSideColors = palette()[ct[rownames(data)]]
  #   names(RowSideColors) = rownames(data)
  # } else {
  #   RowSideColors = rep("white", nrow(data))
  #   names(RowSideColors) = rownames(data)
  #   idx = intersect(rownames(data), names(rsc))
  #   RowSideColors[idx] = rsc[idx]
  # }


  if (is.null(dendrogram)) {
    if (!is.null(Colv) & !is.null(Rowv)) {
      dendrogram="both"
    } else if (!is.null(Rowv)) {
      dendrogram="row"
    } else if (!is.null(Colv)) {
      dendrogram="col"
    } else {
      dendrogram="none"
    }    
  }

  # colors = c("green", "black", "red")
  # colors = c("blue", "yellow", "red")
  # colors = rev(RColorBrewer::brewer.pal(n=11, "RdYlBu"))
  cols = colorRampPalette(colors)(20)
  if (PLOT_MAIN_HM) {
    print(RowSideColors)
    print(Rowv)
    print(dim(data))
    foo = gplots::heatmap.2(data, Rowv=Rowv, Colv=Colv, dendrogram=dendrogram, trace="none", col=cols, main=paste0(main, " (", nrow(data), " genes x ", ncol(data), " samples)"), mar=c(10,5), useRaster=TRUE, RowSideColors=RowSideColors, ColSideColors=ColSideColors, ...)      
  }
  return(list(rsc=RowSideColors, csc=ColSideColors, hc_row=hc_row, hc_col=hc_col))  
}



plot_volcano = function (feats, gse, gs_name, ...) {
  plot(feats[,paste0("l2fc_", gse)], -log10(feats[,paste0("pval_", gse)]), col="grey", xlab="log2FoldChange", ylab="adjusted pval", ...)
  idx = rownames(feats)[feats[[gs_name]]]
  points(feats[idx,][,paste0("l2fc_", gse)], -log10(feats[idx,][,paste0("pval_", gse)]), col=2, pch=1)
  # text(feats[idx,][,paste0("l2fc_", gse)], -log10(feats[idx,][,paste0("pval_", gse)]), idx)
  legend("topleft", pch=1, col=2, paste0(gs_name, " (", length(idx), ")"))  
}


if (!exists("mreadRDS")) mreadRDS = memoise::memoise(readRDS)

get_genes = function(tcga_project, gene_filename="~/projects/genes/bed_grch38_epimeddb.rds") {
  genes = mreadRDS(gene_filename)
  if (!missing(tcga_project)) {
    s_cnv   = mreadstudyRDS(paste0("tcga_studies/study_", tcga_project, "_cnv.rds"))
    s_trscr = mreadstudyRDS(paste0("tcga_studies/study_", tcga_project, "_trscr.rds"))    
    genes = genes[rownames(genes) %in% intersect(rownames(s_trscr$data), rownames(s_cnv$data)),]      
  }
  ## index meth probes by chr
  pf_chr_colname = colnames(genes)[1]
  chrs = as.character(unique(genes[[pf_chr_colname]]))
  chrs_indexed_genes = lapply(chrs, function(chr) {
    # print(chr)
    idx = rownames(genes)[!is.na(genes[[pf_chr_colname]]) & genes[[pf_chr_colname]]==chr]  
    ret = genes[idx,]
    return(ret)
  })
  names(chrs_indexed_genes) = chrs

  return(list(genes=genes, chrs_indexed_genes=chrs_indexed_genes))
}
if (!exists("mget_genes")) { mget_genes = memoise::memoise(get_genes)}

get_neighborhood_genes = function(bed, interaction_range=2500, nb_genes_max=1, START_TO_START=TRUE) {
  genes_singleton = mget_genes()
  genes = genes_singleton$genes

  # bed = genes[1,]
  # bed[[2]] = genes[1,]
  # bed[[3]] = ed_grch38_epimeddbgenes[1,]

  chr = as.character(bed[[1]])
  start = as.numeric(bed[[2]])
  end =   as.numeric(bed[[3]])
  strand = as.character(bed[[6]])
  true_start = ifelse(strand%in%"+", start, end)
  # gene_symbols = rownames(genes[as.character(genes[,1])==chr & genes[,3]>=start-1000 & genes[,2] <= end+1000,1:6])
  tmp_sub_genes = genes_singleton$chrs_indexed_genes[[chr]]
  tmp_sub_genes$tss = tmp_sub_genes[,2]
  tmp_sub_genes[tmp_sub_genes[,6]=="-","tss"] =   tmp_sub_genes[tmp_sub_genes[,6]=="-",3]
  tmp_sub_genes = tmp_sub_genes[tmp_sub_genes[,3]>=start-interaction_range & tmp_sub_genes[,2]<= end+interaction_range,]
  tmp_sub_genes$dist = (tmp_sub_genes$tss - true_start)
  if (START_TO_START) {
    tmp_sub_genes = tmp_sub_genes[abs(tmp_sub_genes$dist) <= interaction_range,]    
  }
  if (nb_genes_max>0) {
    tmp_sub_genes = tmp_sub_genes[order(abs(tmp_sub_genes$dist)),][1:nb_genes_max,]    
  }
  tmp_sub_genes
  tmp_sub_genes =   tmp_sub_genes[!is.na(tmp_sub_genes[,1]),]
  return(tmp_sub_genes)
  ## stress test
  # bed = genes[1,]
  # chr = as.character(bed[[1]])
  # tmp_sub_genes = genes_singleton$chrs_indexed_genes[[chr]]
  # set.seed(1)
  # foo = lapply(round(runif(1000, 1, max(tmp_sub_genes[,3]))), function(i) {
  #   # print(i)
  #   bed[[2]] = i
  #   bed[[3]] = bed[[2]]+1
  #   get_neighborhood_genes(bed, nb_genes_max=3)
  # })
  # sapply(foo, nrow)

  # foo = apply(feat_cpg_binmax[1:50,], 1, get_neighborhood_genes)
  # table(sapply(foo, nrow))
  # foo = do.call(rbind, foo)
  # head(foo)
  # dim(foo)
}


plot_feat_den = function (x, y, xlim, ylim, ...) {
  if (missing(xlim)) {
    xlim = range(x)    
  }
  if (missing(ylim)) {
    ylim = range(y)    
  }
  layout(matrix(c(2,1,1,2,1,1,4,3,3), 3), respect=TRUE)
  # plot(x, y, pch=".", xlim=xlim, ylim=ylim)
  smoothScatter (x, y, xlim=xlim, ylim=ylim, ...)
  den_x = density(x,)
  plot(den_x$x, den_x$y, type="l", xlim=xlim, xlab="")
  den_y = density(y)
  # den_y = density(y, breaks=c(-1, seq(min(ylim), max(ylim)), max(y)+1))
  plot(den_y$y, den_y$x, type="l", ylim=ylim, ylab="")  
}

get_matbin_from_seq = function(seq, up_str, dwn_str, step) {
  foo = seq(1,up_str+dwn_str, step)
  bins = cbind(foo, foo+step-1)
  mat = t(epimedtools::monitored_apply(t(t(seq)), 1, function(s) {
    apply(bins, 1, function(b) {
      # b = bins[1,]
      str = substr(s, b[[1]], b[[2]])
      # l = nchar(str)
      lnn = nchar(gsub("N", "", str))
      # denom = l - lnn
      if (lnn==0) {
        return(NA)
      } else {
        return(stringr::str_count(str, "CG") / lnn)
      }
    })
  }))
  offset = (bins[,2] - ud_str - step/2)/1000
  colnames(mat) = paste0("binmax ", ifelse(offset<0, "-", "+"), abs(offset), "kb")
  return(mat)  
}
if (!exists("mget_matbin_from_seq")) {mget_matbin_from_seq = memoise::memoise(get_matbin_from_seq)}


get_seq_from_bed = function(bed, up_str, dwn_str, chrom_sizes, genome=BSgenome.Hsapiens.UCSC.hg38::Hsapiens) {
  seq = epimedtools::monitored_apply(mod=10, bed, 1, function(gene) {
    chr = gene[[1]]
    strand = gene[[6]]
    if (strand == "+") {
      tss = as.numeric(gene[[2]])
      bef = up_str
      aft = dwn_str
    } else {
      tss = as.numeric(gene[[3]])
      bef = dwn_str
      aft = up_str
    }  
    beg = max(1,tss - bef)
    end = min(chrom_sizes[chr,2], tss + aft)
    str = as.character(BSgenome::getSeq(genome, chr, beg, end))
    return(str)
  })
  return(seq)
}


build_feats_epic_grch38 = function() {
  # params
  extend_region_dist = 1000
  # meth pf
  pf_chr_colname = "seqnames"
  pf_pos_colname = "start"
  pf_orig = mreadRDS("~/projects/datashare/platforms/EPIC.hg38.manifest.full.fch.rds")
  pf_orig = pf_orig[pf_orig[,pf_pos_colname]>0,]    
  pf_orig = pf_orig[order(pf_orig[[pf_chr_colname]],pf_orig[[pf_pos_colname]]), ]
  ## index meth probes by chr
  chrs = unique(pf_orig[[pf_chr_colname]])
  chrs_indexed_methpf = lapply(chrs, function(chr) {
    print(chr)
    idx = rownames(pf_orig)[!is.na(pf_orig[[pf_chr_colname]]) & pf_orig[[pf_chr_colname]]==chr]  
    ret = pf_orig[idx,]
    return(ret)
  })
  names(chrs_indexed_methpf) = chrs

  fat_feat = lapply(unique(pf_orig[,1]), function(chr) {
    d = pf_orig[pf_orig[,1]==chr,]
    i = intervals::Intervals(c(d[,2], d[,3]), type="Z")
    # enlarge your fat feat
    l = extend_region_dist
    c = intervals::close_intervals( intervals::contract( intervals::reduce(intervals::expand(i, l)), l) )
    dim(c)
    df = data.frame(chr, c[,1], c[,2])
    return(df)
  })
  fat_feat = do.call(rbind, fat_feat)
  dim(fat_feat)
  fat_feat[,4] = paste0(fat_feat[,1], ":", fat_feat[,2], "-", fat_feat[,3])
  fat_feat[,5] = fat_feat[,3] - fat_feat[,2]
  fat_feat[,6] = "+"
  fat_feat = fat_feat[fat_feat[,5]>1,]
  rownames(fat_feat) = fat_feat[,4]
  colnames(fat_feat)  = c("chr", "start", "end", "id", "score", "strand")
  dim(fat_feat)
  head(fat_feat)

  ## index probes by feat name
  print("# indexing probes by feat name")
  feat_indexed_probes = epimedtools::monitored_apply(fat_feat, 1, function(feat) {
    # feat = fat_feat[3,]
    # print(feat)
    chr = feat[[1]]
    len = as.numeric(feat[[5]])
    meth_platform = chrs_indexed_methpf[[chr]]
    ret = dmprocr::get_probe_names(feat, meth_platform, pf_chr_colname, pf_pos_colname, 0, len) 
    # meth_platform[ret,1:3]
    # feat
    return(ret)
  })
  
  nb_probes = sapply(feat_indexed_probes, length)
  fat_feat$score = nb_probes[rownames(fat_feat)]
  saveRDS(fat_feat, "feats_epic_grch38.rds")
  # return("feats_epic_grch38.rds")
  # layout(matrix(1:2, 1), respect=TRUE)
  # barplot(table(sapply(feat_indexed_probes, length)), las=2, main="#probes")
  # plot(fat_feat$end - fat_feat$start, fat_feat$score, main="#probes")
  
  # options(scipen=999)
  # probes = unique(unlist(feat_indexed_probes))
  # pf = pf_orig[probes, c(pf_chr_colname, pf_pos_colname)]
  # pf[,3] = pf[,2] + 1
  # pf[,4] = probes
  # pf[,5] = 1
  # pf[,6] = "+"
  # write.table(pf  , file="pf.bed"  , sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)
  return(feat_indexed_probes)
}









# from ~/project/01_momik/results/common.R

get_feat_indexed_probes = function(feats_bed6, probes_bed2, up_str, dwn_str) {
  pf_bed = probes_bed2
  pf_bed[,3] = pf_bed[,2] + 1
  pf_bed[,4] = rownames(pf_bed)
  pf_bed[,5] = 0
  pf_bed[,6] = "+"
  pf_bed = pf_bed[!pf_bed[,1]%in%"*",]
  head(pf_bed)
  # options(scipen=999)
  write.table(pf_bed, file="illumina_450k_hg38.bed", sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)    
  ## index meth probes by chr
  chrs = unique(feats_bed6[,1])
  chrs_indexed_methpf = lapply(chrs, function(chr) {
    print(chr)
    idx = rownames(pf_bed)[!is.na(pf_bed[,1]) & pf_bed[,1]==chr]  
    ret = pf_bed[idx,]
    return(ret)
  })
  names(chrs_indexed_methpf) = chrs
  ## index probes by gene name
  print("# indexing probes by gene name")
  feat_indexed_probes = epimedtools::monitored_apply(feats_bed6, 1, function(gene) {
    # gene = randomall_genes[1,]genes=readRDS("~/fchuffar/projects/genes/bed_grch38_epimeddb.rds")
    # print(gene)
    chr = gene[[1]]
    meth_platform = chrs_indexed_methpf[[chr]]
    ret = dmprocr::get_probe_names(gene, meth_platform, 1, 2, up_str, dwn_str) 
    return(ret)
  })
  barplot(table(sapply(feat_indexed_probes, length)))
  return(feat_indexed_probes)
}
if (!exists("mget_feat_indexed_probes")) mget_feat_indexed_probes = memoise::memoise(get_feat_indexed_probes)



get_multiomic_data = function(gene_symbols, tcga_project, feat_indexed_probes, region_id, interaction_range=2500) {
  # # warning: feat_indexed_probes is a global variable
  s_cnv   = mreadstudyRDS(paste0("tcga_studies/study_", tcga_project, "_cnv.rds"))
  s_meth  = mreadstudyRDS(paste0("tcga_studies/study_", tcga_project, "_meth.rds"))
  s_trscr = mreadstudyRDS(paste0("tcga_studies/study_", tcga_project, "_trscr.rds"))
  genes_singleton = mget_genes(tcga_project)
  genes = genes_singleton$genes
  if (missing(feat_indexed_probes)) {
    # params
    pf_chr_colname = "Chromosome"
    pf_pos_colname = "Start"
    if (missing(region_id)) {
      feat = genes[gene_symbols,]
      chr = feat[[1]]
      # meth pf
      if (!exists("pf_orig")) {
        pf_orig = s_meth$platform
        pf_orig = pf_orig[pf_orig[,pf_pos_colname]>0,]    
        pf_orig = pf_orig[order(pf_orig[[pf_chr_colname]],pf_orig[[pf_pos_colname]]), ]
      }  
      meth_platform = pf_orig[pf_orig[[pf_chr_colname]]%in%feat[[1]], ]
      head(meth_platform)
      probes = dmprocr::get_probe_names(feat, meth_platform, pf_chr_colname, pf_pos_colname, interaction_range, interaction_range) 

      # print(feat_indexed_probes)
      tss = ifelse(feat[[6]]=="+", feat[[2]], feat[[3]])
      region_id = paste0(feat[[1]], ":", tss-interaction_range, "-", tss+interaction_range)      
      feat_indexed_probes = list()
      feat_indexed_probes[[gene_symbols]] = feat_indexed_probes[[region_id]] = probes
      feat_indexed_probes
    } else {
      tmp_pf = strsplit(region_id, ":|-", perl=TRUE)[[1]]
      feat = genes[1,1:6]
      feat[[1]] = tmp_pf[[1]]
      feat[[2]] = as.numeric(tmp_pf[[2]])
      feat[[3]] = as.numeric(tmp_pf[[3]])
      feat[[4]] = "foo"
      feat[[5]] = 0
      chr = feat[[1]]
      # meth pf
      if (!exists("pf_orig")) {
        pf_orig = s_meth$platform
        pf_orig = pf_orig[pf_orig[,pf_pos_colname]>0,]    
        pf_orig = pf_orig[order(pf_orig[[pf_chr_colname]],pf_orig[[pf_pos_colname]]), ]
      }  
      meth_platform = pf_orig[pf_orig[[pf_chr_colname]]%in%feat[[1]], ]
      head(meth_platform)
      probes = dmprocr::get_probe_names(feat, meth_platform, pf_chr_colname, pf_pos_colname, 0, feat[[3]]-feat[[2]])       
      feat_indexed_probes = list()
      feat_indexed_probes[[region_id]] = probes
      feat_indexed_probes
    }
    # print(feat_indexed_probes)
    # print(region_id)
  }
  if (missing(gene_symbols)) {
    tmp_reg = strsplit(region_id, ":|-", perl=TRUE)[[1]]
    chr = tmp_reg[1]
    start = as.numeric(tmp_reg[2])
    end =   as.numeric(tmp_reg[3])
    # gene_symbols = rownames(genes[as.character(genes[,1])==chr & genes[,3]>=start-1000 & genes[,2] <= end+1000,1:6])
    tmp_sub_genes = genes_singleton$chrs_indexed_genes[[chr]]
    gene_symbols = rownames(tmp_sub_genes[tmp_sub_genes[,3]>=start- interaction_range & tmp_sub_genes[,2] <= end+interaction_range,1:6])
    if (length(gene_symbols) == 0) {
      return(NULL)
    }
  }
  if (missing(region_id)) {
    gene_symbol = gene_symbols[1]
    tmp_pf = get_pf_from_feat_indexed_probes(feat_indexed_probes)
    tmp_pf = tmp_pf[as.character(tmp_pf$chr)==as.character(genes[gene_symbol, 1]), ]
    tss = ifelse(genes[gene_symbol,6]=="+", genes[gene_symbol,2], genes[gene_symbol,3])
    reg = tmp_pf[tmp_pf[,2]<tss&tmp_pf[,3]>tss, ]
    region_id = paste0(reg[1], ":", reg[2], "-", reg[3])
  }
  # debugged and optimized version of dmprocr::trscr_meth_analysis https://github.com/bcm-uga/dmprocr
  preproc_omics_data = function(region_id, gene_symbols, s_cnv, s_meth, s_trscr, feat_indexed_probes) {
    # meth_data
    # meth_probe_idx = intersect(feat_indexed_probes[[region_id]], rownames(s_meth$data))
    meth_probe_idx = feat_indexed_probes[[region_id]]
    if (length(meth_probe_idx) <= 1) {
        return(NULL)
    }
    tmp_reg = strsplit(region_id, ":|-", perl=TRUE)[[1]]
    chr = tmp_reg[1]
    meth_data = s_meth$stuffs$chrs_indexed_data[[chr]][meth_probe_idx, ]
    meth_data = meth_data[, apply(is.na(meth_data), 2, sum)/nrow(meth_data) < 0.5]
    meth_data = meth_data[apply(is.na(meth_data), 1, sum)/ncol(meth_data) < 0.5, ]
    # dim(meth_data)
    meth_probe_idx = rownames(meth_data)
    # idx_sample according to s_cnv if needed
    if (length(gene_symbols)==1) {
      tmp_gene_symbols = c(gene_symbols, gene_symbols)
    } else {
      tmp_gene_symbols = gene_symbols
    }
    gene_symbol = gene_symbols[1]
    FAST = FALSE
    if (FAST) {
      idx_sample = intersect(colnames(s_trscr$data), colnames(meth_data))
    } else {
      if (!is.null(s_cnv)) {
        idx_sample = intersect (
          intersect(
            colnames(s_trscr$data)[order(s_trscr$data[gene_symbol,])],
            colnames(meth_data)
          ),
          colnames(s_cnv$data)[apply(abs(s_cnv$data[tmp_gene_symbols, ]) < 0.2, 2, all)]
        )
      } else {
        idx_sample = intersect(
          colnames(s_trscr$data)[order(s_trscr$data[gene_symbol,])],
          colnames(meth_data)
        )
      }
      if (length(idx_sample) <= 1) {
        return(NULL)
      }      
    }
    d = data.frame(t(meth_data[, idx_sample]))
    tmp_trscr_data = data.frame(t(data.frame(s_trscr$data[tmp_gene_symbols, idx_sample])))
    genes_to_keep = sapply(gene_symbols, function(g) {
      if (length(unique(tmp_trscr_data[[g]])) == 1) {
        return(FALSE)
      } else {
        return(TRUE)
      }
    })
    gene_symbols = gene_symbols[genes_to_keep]
    for (g in gene_symbols) {
      d[[g]] = tmp_trscr_data[[g]]
    }
    ret = list(
      d=d,
      probes=meth_probe_idx,
      gene_symbols=gene_symbols,
      region_id=region_id,
      tcga_project=tcga_project
    )
    return(ret)
  }
  ret = preproc_omics_data(region_id, gene_symbols, s_cnv, s_meth, s_trscr, feat_indexed_probes)
  return(ret)
}
if (!exists("mget_multiomic_data")) {mget_multiomic_data = memoise::memoise(get_multiomic_data)}


momic_pattern = function(gene_symbols, tcga_project, interaction_range=2500, LAYOUT=TRUE, ...) {    
  data = get_multiomic_data(gene_symbols=gene_symbols, tcga_project=tcga_project, interaction_range=interaction_range, ...)
  if (LAYOUT) {layout(matrix(c(2,1,1,1,1), 1))}
  # layout(matrix(c(1, 1, 2, 2, 2, 2), 2), respect=TRUE)
  # transcriptome
  # par(mar=c(10, 4.1, 4.1, 2.1))
  # methylome
  colors = c("cyan", "black", "red")
  cols = colorRampPalette(colors)(20)
  breaks = seq(0, 1, length.out = length(cols) + 1)
  main = paste0(tcga_project, " ", data$gene_symbol, " expression and TSS+/-", interaction_range/1000, "kb methylation")
  # par(mar=c(10, 4.1, 4.1, 2.1))

  par(mar=c(5.7, 0, 4.1, 2.1))
  image(t(data$d[,data$probes]), col=cols, breaks=breaks, xaxt="n", yaxt="n", main=main)
  axis(1, (1:nrow(t(data$d[,data$probes])) - 1)/(nrow(t(data$d[,data$probes])) - 1), rownames(t(data$d[,data$probes])), las = 2)

  par(mar=c(5.7, 4.1, 4.1, 0))
  plot(data$d[,gene_symbols[1]], (seq(0,1,length=nrow(data$d)+1)[-1]) - .5/nrow(data$d), 
    main="", 
    xlab="expr.", 
    ylab=paste0(nrow(data$d), " samples"), 
    yaxt="n",
    ylim=c(0,1), 
    type="l",
    lwd=3,
    yaxs = "i"
  ) 
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  return(data)
}

readstudyRDS = function(rds_file){
  s = readRDS(rds_file)
  rownames(s$data) = gsub("/", "_", gsub("-", "_", rownames(s$data)))
  ## index platform by chr
  pf_chr_colname = colnames(s$platform)[1]
  if (pf_chr_colname == "Chromosome") { # Fix it!!
    chrs = as.character(unique(s$platform[[pf_chr_colname]]))
    s$stuffs$chrs_indexed_platform = lapply(chrs, function(chr) {
      # print(chr)
      idx = rownames(s$platform)[!is.na(s$platform[[pf_chr_colname]]) & s$platform[[pf_chr_colname]]==chr]  
      ret = s$platform[idx,]
      return(ret)
    })
    names(s$stuffs$chrs_indexed_platform) = chrs

    s$stuffs$chrs_indexed_data = lapply(chrs, function(chr) {
      # print(chr)
      idx = rownames(s$stuffs$chrs_indexed_platform[[chr]])
      ret = s$data[idx,]
      return(ret)
    })
    names(s$stuffs$chrs_indexed_data) = chrs      
  }
  return(s)
}

if (!exists("mreadstudyRDS")) {mreadstudyRDS = memoise::memoise(readstudyRDS)}




compute_anadiff = function(s, idx_ctl, idx_ttt) {
  d = data.frame(treatment_5aza=s$exp_grp[c(idx_ctl,idx_ttt), "treatment_5aza"], cell_line=s$exp_grp[c(idx_ctl,idx_ttt), "cell_line" ])
  result = epimedtools::monitored_apply(s$data, 1, function(l){
    #l = s$data[1,]
    d$expr = l[c(idx_ctl, idx_ttt)]
    m = lm(expr~treatment_5aza+cell_line, d)
    beta = m$coefficients[[2]]
    pval = anova(m)[1,5]
    ret = c(beta=beta, pval=pval)
    return(ret)
  })
  return(result)
}
if (!exists("mcompute_anadiff")) {mcompute_anadiff = memoise::memoise(compute_anadiff)}

