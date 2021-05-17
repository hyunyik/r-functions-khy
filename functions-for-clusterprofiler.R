run_ego <- function(gene, ont, pvalueCutoff = 0.05, nCluster = null) {
  require(GOSemSim)
  require(enrichplot)
  require(clusterProfiler)
  require(ggplot2)
  
  d <- godata('org.Mm.eg.db', ont = ont)
  ego <- enrichGO(gene = names(deg.sig.gene.list), OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = ont, pvalueCutoff = 0.05)
  ego@result[ego@result$pvalue %>% order,]
  ego <- pairwise_termsim(ego, semData = d, showCategory = dim(ego@result)[1])
  p <- emapplot_cluster(ego, showCategory = dim(ego@result)[1], repel = T, cluster_num_label = T, cex_label_group = 2)
  p <- p + xlim(min(p$data$x)-1,max(p$data$x)+1) +ylim(min(p$data$y)-1,max(p$data$y)+1)
  p
  
  go_cls_list <- list()
  for (i in sort(unique(p$data$color))) {
    pp <- p$data[p$data$color==i,][,c(3, 4, 5)]
    colnames(pp) <- c(paste0("Cluster ",i," GO terms"), "Size", "p-Value")
    pp <- pp[order(pp$Size, decreasing = T), ]
    rownames(pp) <- NULL
    print(pp)
    go_cls_list[paste0("Cluster ",i," GO terms")] <- pp
    write.csv(pp, file = paste0(ont, "_GO_cluster_",i ,".csv"))
  }
  
  r <- list()
  r[["enrichGO"]] <- ego
  r[["p.GOcluster"]] <- p
  r[["list.GOcluster"]] <- go_cls_list
  return(r)
}

run_emap_w_ego_termsim <- function(ego, d, pvalueCutoff = 0.05, nCluster = NULL, showCategory = NULL) {
  require(GOSemSim)
  require(enrichplot)
  require(clusterProfiler)
  require(ggplot2)
  
  if(is.null(showCategory)) { showCategory <- dim(ego@result)[1] }
  
  ego <- pairwise_termsim(ego, semData = d, showCategory = showCategory)
  p <- emapplot_cluster(ego, showCategory = dim(ego@result)[1], repel = T, cluster_num_label = T, cex_label_group = 2, nCluster = nCluster)
  p <- p + xlim(min(p$data$x)-1,max(p$data$x)+1) +ylim(min(p$data$y)-1,max(p$data$y)+1)
  p
  
  go_cls_list <- list()
  for (i in sort(unique(p$data$color))) {
    pp <- p$data[p$data$color==i,][,c(3, 4, 5)]
    colnames(pp) <- c(paste0("Cluster ",i," GO terms"), "Size", "p-Value")
    pp <- pp[order(pp$Size, decreasing = T), ]
    rownames(pp) <- NULL
    print(pp)
    go_cls_list[paste0("Cluster ",i," GO terms")] <- pp
    #write.csv(pp, file = paste0(d@ont, "_GO_cluster_",i ,".csv"))
  }
  
  r <- list()
  r[["enrichGO"]] <- ego
  r[["p.GOcluster"]] <- p
  r[["list.GOcluster"]] <- go_cls_list
  return(r)
}

draw_gene_network <- function(d, deg.sig.gene.list, cluster.num, plot.margin = 2, ) {
  go_cls_list <- r.BP$list.GOcluster
  
  d  %>% slice(grep(paste(go_cls_list[cluster.num] %>% unlist, collapse = "|"), d@result$Description)) %>% cnetplot(circular = F, colorEdge = T, cex_label_gene = 1.2, node_label = "gene", layout = "dh", foldChange = deg.sig.gene.list, cex_category = 1.5, showCategory = 5) -> cp
  cp <- cp + xlim(min(cp$data$x)-plot.margin, max(cp$data$x)+plot.margin) +ylim(min(cp$data$y)-plot.margin, max(cp$data$y)+plot.margin) + 
    ggtitle(paste0("Gene Network of GO Cluster ", cluster.num)) + 
    theme(plot.title = element_text(size = 30, hjust = 0.5)) +
    theme(legend.direction = "vertical", legend.box = "horizontal", legend.position = "bottom")
  cp <- cp + theme(legend.title = element_text(size = 20, face = "bold"), legend.margin = margin(r = 10, l = 10), legend.text = element_text(size = 18))
  cp
}