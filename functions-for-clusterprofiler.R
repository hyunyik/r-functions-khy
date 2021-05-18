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
  require(clusterProfiler.dplyr)
  
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

draw_gene_network <- function(run_empa_r, deg.sig.gene.list, cluster.num, plot.margin = 2) {
  go_cls_list <- run_empa_r$list.GOcluster
  d <- run_empa_r$enrichGO
  d %>% clusterProfiler.dplyr::slice(grep(paste(go_cls_list[cluster.num] %>% unlist, collapse = "|"), d@result$Description)) %>% cnetplot(circular = F, colorEdge = T, cex_label_gene = 1.2, node_label = "gene", layout = "dh", foldChange = deg.sig.gene.list, cex_category = 1.5, showCategory = 5) -> cp
  cp <- cp + xlim(min(cp$data$x)-plot.margin, max(cp$data$x)+plot.margin) +ylim(min(cp$data$y)-plot.margin, max(cp$data$y)+plot.margin) + 
    ggtitle(paste0("Gene Network of GO Cluster ", cluster.num)) + 
    theme(plot.title = element_text(size = 30, hjust = 0.5)) +
    theme(legend.direction = "vertical", legend.box = "horizontal", legend.position = "bottom")
  cp <- cp + theme(legend.title = element_text(size = 20, face = "bold"), legend.margin = margin(r = 10, l = 10), legend.text = element_text(size = 18))
  cp
}

## https://github.com/YuLab-SMU/DOSE/blob/master/R/gsea.R
## GSEA algorithm (Subramanian et al. PNAS 2005)
## INPUTs to GSEA
## 1. Expression data set D with N genes and k samples.
## 2. Ranking procedure to produce Gene List L.
## Includes a correlation (or other ranking metric)
## and a phenotype or profile of interest C.
## 3. An exponent p to control the weight of the step.
## 4. Independently derived Gene Set S of N_H genes (e.g., a pathway).
## Enrichment Score ES.
## 2. Evaluate the fraction of genes in S ("hits") weighted
## by their correlation and the fraction of genes not in S ("miss")
## present up to a given position i in L.
gseaScores <- function(geneList, geneSet, leading.edge, exponent=1, fortify=FALSE) {
  ###################################################################
  ##    geneList                                                   ##
  ##                                                               ##
  ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
  ##    according to the correlation, r(g_j)=r_j,                  ##
  ##    of their expression profiles with C.                       ##
  ##                                                               ##
  ###################################################################
  
  ###################################################################
  ##    exponent                                                   ##
  ##                                                               ##
  ## An exponent p to control the weight of the step.              ##
  ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
  ##   the standard Kolmogorov-Smirnov statistic.                  ##
  ##   When p = 1, we are weighting the genes in S                 ##
  ##   by their correlation with C normalized                      ##
  ##   by the sum of the correlations over all of the genes in S.  ##
  ##                                                               ##
  ###################################################################
  
  ## genes defined in geneSet should appear in geneList.
  ## this is a must, see https://github.com/GuangchuangYu/DOSE/issues/23
  geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  leading.genes <- unlist(strsplit(leading.edge, "/"))
  print(length(leading.genes))
  #leading.genes <- ifelse(length(leading.genes)<40, leading.genes, leading.genes[1:40])
  #print(leading.genes)
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits),
                   gene.name=names(geneList),
                   leading.edge=names(geneList) %in% leading.genes[1:40]
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}

## https://github.com/YuLab-SMU/enrichplot/blob/master/R/gseaplot.R
##' extract gsea result of selected geneSet
##'
##'
##' @title gsInfo
##' @param object gseaResult object
##' @param geneSetID gene set ID
##' @return data.frame
##' @author Guangchuang Yu
## @export
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  leading.edge <- object@result[geneSetID,]$core_enrichment
  df <- gseaScores(geneList, geneSet, leading.edge, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
