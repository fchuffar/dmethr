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


momic_pattern = function(gene_symbols, tcga_project, interaction_range=2500, ...) {    
  data = get_multiomic_data(gene_symbols=gene_symbols, tcga_project=tcga_project, interaction_range=interaction_range, ...)
  par(mar=c(10, 4.1, 4.1, 2.1))
  layout(matrix(c(1,2,2,2,2), 1))
  # layout(matrix(c(1, 1, 2, 2, 2, 2), 2), respect=TRUE)
  # transcriptome
  # par(mar=c(10, 4.1, 4.1, 2.1))
  plot(data$d[,gene_symbols[1]], 1:length(data$d[,gene_symbols[1]]), 
    main=paste(data$gene_symbol, tcga_project), 
    xlab="log2(normalized expression)", 
    ylab=paste0(nrow(data$d), " samples"), 
    yaxt="n"
  ) 
  # methylome
  colors = c("cyan", "black", "red")
  cols = colorRampPalette(colors)(20)
  breaks = seq(0, 1, length.out = length(cols) + 1)
  main = paste0("methylome TSS+/-", interaction_range)
  # par(mar=c(10, 4.1, 4.1, 2.1))
  image(t(data$d[,data$probes]), col=cols, breaks=breaks, xaxt="n", 
      yaxt="n", main=main)
  axis(1, (1:nrow(t(data$d[,data$probes])) - 1)/(nrow(t(data$d[,data$probes])) - 1), rownames(t(data$d[,data$probes])), 
      las = 2)
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

