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
