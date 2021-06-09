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
