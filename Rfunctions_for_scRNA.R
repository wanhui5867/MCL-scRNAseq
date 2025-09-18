library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
library(Matrix)
library(DoubletFinder)
library(scRepertoire)
library(harmony)
library(clustree)
library(scPred)
library(RColorBrewer)
library(ggpubr)
library(eoffice)
library(grid)
library(fgsea)
library(EnhancedVolcano)
library(clusterProfiler)
library(pheatmap)
gene.lymphchip = read_lines('~/OneDrive - Karolinska Institutet//Mac/Project/MCL50_2022/1.driver/_Ref/lymphChip.gene.list')

myheatmap <- function(mat,...) {
  bk = c(seq(-2,2,0.1))
  pheatmap(mat,
           color=c(colorRampPalette(colors = c("blue","black"))(length(bk)/2),colorRampPalette(colors = c("black","yellow"))(length(bk)/2)),
           show_rownames=T,
           show_colnames=T,
           scale="row",
           # cluster_rows = T, 
           # cluster_cols = T,
           breaks = bk,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "complete",
           ...
  )
  
} 
#myheatmap(avgexp)



# add lable
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  heatmap <- pheatmap$gtable
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}




#  # example for heatmap -----

# MATRIX
# mat_ave_cluster <- AverageExpression(q.all.nb, return.seurat = T)@assays$RNA@data
# mat_cluster     <- mat_ave_cluster[top5$gene,]

# LABEL
# markers_genes <- FindAllMarkers(q.all.nb, logfc.threshold = log(1.5), test.use = "wilcox", 
#                                 min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE,  random.seed = 1234, 
#                                 assay = "RNA") %>% 
#   group_by(cluster) %>% filter(p_val_adj <= 0.05) 
# top5 <- markers_genes %>% group_by(cluster) %>% top_n(-5, avg_log2FC)

# HEATMAP
# myheatmap(mat_cluster)
# heat <- myheatmap(mat_ave_cluster[unique(markers_genes$gene), ])  # affected by cluster size
# 
# # according to cluster order to order genes
# heat.order = markers_genes %>% 
#   arrange(match(cluster, c( 'Bplasma-1', 'Bplasma-2',  'Bnaive-1', 'Bnaive-2', 'plasmablast', 'Bmemory','B1/Breg')), -avg_log2FC) 
# heat <- myheatmap(mat_ave_cluster[unique(heat.order$gene), heat.order$cluster ], cluster_rows = F, cluster_cols = F)  # affected by cluster size

# pdf(file = 'q.b.cells/heatmap.markers.cluster.pdf', h = 7, w = 5)
# add.flag(heat, kept.labels = top5$gene,  repel.degree = 0.2)
# dev.off()



# total cells for each cell types -----------
plotbar.total.cells <-function(meta){
  ggplot(meta %>% select(HTO_classification, cluster ) %>% 
           mutate(cluster = factor(cluster, levels = order_cluster )) , 
         aes( x = cluster,  fill = HTO_classification)) +
    geom_bar( width = 0.8) +
    #facet_wrap( ~ cluster) +
    scale_fill_manual(values = color_sample) +
    ylab('No. of cells') +
    theme_pubr() +   rremove('xlab')
} 



myVolcano <- function(res, ident.1, ident.2,  FCcutoff = log2(1.5), pCutoff = 1e-02, ...) {
  # FCcutoff = log(1.2)
  # pCutoff = 1e-02
  # ident.1 = "Bnaive-1"
  # ident.2 = "Bnaive-2"
  # show.genes = markers.b.activated
  
  res = res %>% 
    mutate(color = if_else(avg_log2FC > FCcutoff & p_val_adj < pCutoff, 'red2', if_else(avg_log2FC < -FCcutoff & p_val_adj < pCutoff, 'royalblue', 'grey60')),
           change = if_else(avg_log2FC > FCcutoff & p_val_adj < pCutoff, ident.1 , if_else(avg_log2FC < -FCcutoff & p_val_adj < pCutoff, ident.2, 'NS') ))
  point.color.vec = res$color
  names(point.color.vec) = res$change
  
  EnhancedVolcano(res,  #MODIFY
                  lab = rownames(res), #MODIFY
                  #selectLab = show.genes,
                  x = 'avg_log2FC',
                  y = 'p_val_adj',
                  border = 'full',
                  pCutoff = pCutoff,
                  FCcutoff = FCcutoff,  # MODIFY
                  cutoffLineWidth = 0.7,
                  cutoffLineType = "dashed",
                  pointSize = 2.0,
                  labSize = 4.0,
                  # point color default
                  # col = c("grey60", "grey60", "grey60", "red2"),
                  # legendLabels = c("NS", "avg_logFC", "p_val_adj", "p_val_adj & avg_logFC"),
                  # point color by group
                  colCustom = point.color.vec,
                  colAlpha = 1,
                  widthConnectors = 0.5,
                  drawConnectors = TRUE,
                  colConnectors = "black",
                  title = paste(ident.1, 'vs', ident.2) , # MODIFY
                  titleLabSize = 16, subtitle = NULL,
                  legendLabSize = 16, legendPosition = 'right',
                  caption = paste0("Log2 fold change cutoff: ", round(FCcutoff,2), "; p.ajust cutoff: ", pCutoff),
                  ylab = bquote(~-Log[10]~"adj-p"), ...) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black', size = 1),
          plot.background = element_rect(fill = 'white'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
}

# markers_genes.naive.all <- FindMarkers(q.b.naive, ident.1 = "Bnaive-1", ident.2 = "Bnaive-2", min.pct = 0.1, logfc.threshold = 0)  # input res data
# myVolcano(markers_genes.naive.all, ident.1 = "Bnaive-1", ident.2 = "Bnaive-2", FCcutoff = log2(1.2), pCutoff = 1e-02 ) # label all changed genes
# myVolcano(markers_genes.naive.all, ident.1 = "Bnaive-1", ident.2 = "Bnaive-2", selectLab = markers.b.activated , xlim = c(-1.5, 1.5), ylim = c(0,10)) # label selected genes & ajust x and y limitation
# ggsave(p, filename =  'q.b.cells/Volcano.Bnaive.selectgenes.pdf', w = 9, h = 6)

# pct cells  for each cell types
plotbar.pct.cells <-function(meta){
  ggplot(meta %>% select(HTO_classification, cluster ) %>% 
           mutate(cluster = factor(cluster, levels = order_cluster )) , 
         aes( x = cluster,  fill = HTO_classification)) +
    geom_bar(position = 'fill', width = 0.8) +
    #facet_wrap( ~ cluster) +
    scale_fill_manual(values = color_sample) +
    ylab('Pct. of cells') +
    theme_pubr() +   rremove('xlab')
} 


# total cells for each sample
plotbar.total.samples <-function(meta){
  ggplot(meta %>% select(HTO_classification, cluster ) %>% 
           mutate(cluster = factor(cluster, levels = order_cluster )) , 
         aes( fill = cluster,  x = HTO_classification)) +
    geom_bar( width = 0.8) +
    #facet_wrap( ~ cluster) +
    scale_fill_brewer(palette="Set2") +
    ylab('No. of cells') +
    theme_pubr() +  rremove('xlab')
} 


# pct cells for each sample
plotbar.pct.samples <-function(meta){
  ggplot(meta %>% select(HTO_classification, cluster ) %>% 
           mutate(cluster = factor(cluster, levels = order_cluster )) , 
         aes( fill = cluster,  x = HTO_classification)) +
    geom_bar(position = 'fill', width = 0.8) +
    #facet_wrap( ~ cluster) +
    scale_fill_brewer(palette="Set2") +
    ylab('No. of cells') +
    theme_pubr() +  rremove('xlab')
} 




# proport of cells for each sample
plotbar.sample.cell <-function(meta){
  
  df.sample.cell = meta %>% 
    select(HTO_classification, cluster ) %>% 
    group_by(cluster, HTO_classification) %>% 
    summarise(N_cells = n()) %>% 
    ungroup() %>% group_by(HTO_classification) %>%   
    mutate(`% of cells` = N_cells/sum(N_cells)*100) %>% 
    mutate(cluster = factor(cluster, levels = order_cluster ))
  
  
  ggbarplot(df.sample.cell, x = 'HTO_classification', y = '% of cells', 
            position = position_identity(),
            facet.by =  'cluster', nrow = 1, 
            #scale = 'free_y',
            fill = 'HTO_classification', palette = color_sample) +
    rremove('legend') + rremove('xlab')
}


# proport of cells for each sample
plotbar.cell.sample <-function(meta){
  
  df.sample.cell = meta %>% 
    select(HTO_classification, cluster ) %>% 
    group_by(cluster, HTO_classification) %>% 
    summarise(N_cells = n()) %>% 
    ungroup() %>% group_by(HTO_classification) %>%   
    mutate(`% of cells` = N_cells/sum(N_cells)*100) %>% 
    mutate(cluster = factor(cluster, levels = order_cluster ))
  
  
  ggbarplot(df.sample.cell, x = 'cluster', y = '% of cells', 
            position = position_identity(),
            facet.by =  'HTO_classification', nrow = 1, 
            #scale = 'free_y',
            fill = 'cluster') +
    rremove('legend') + rremove('xlab')
}


# proport of cells for each sample
plotbar.sample.gene <-function(meta){
  
  ggboxplot(meta, x = 'HTO_classification', y = 'nFeature_RNA', 
            outlier.shape = NA, 
            ylab = 'No. of genes per cell',
            order = order_sample,
            #position = position_identity(),
            facet.by =  'cluster', nrow = 1, 
            #scale = 'free_y',
            fill = 'HTO_classification', palette = color_sample) +
    rremove('legend') + rremove('xlab')
}


# calculate the proportion of cells expressed specific gene(s)
PrctCellExpringGene <- function(object, genes, group.by = "all", assay = "RNA", datatype = "counts", threshold = 0){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object, assay = assay, datatype=datatype, threshold=threshold))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    results = lapply(list, PrctCellExpringGene, genes=genes, assay = assay, datatype=datatype, threshold=threshold)
    results %>% reduce(full_join, by="Markers") %>% select(any_of("Markers")) -> genelist
    results %>% reduce(full_join, by="Markers") %>% select(!any_of("Markers")) %>% "colnames<-"(names(results)) -> percentages
    combined <- cbind(genelist,percentages)
    return(combined)
  }
}

calc_helper <- function(object,genes,assay,datatype,threshold){
  counts = slot(object[[assay]],datatype)
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    round((sum(counts[genes,]>threshold)/ncells)*100,1)
  }else{return(NA)}
}

# GSVA geneset collection
geneSetFunc <- function(listIn){
  sets <- Map(GeneSet, 
              listIn, 
              setName=names(listIn),
              MoreArgs=list(geneIdType=SymbolIdentifier()))
  GeneSetCollection(sets)
}
