#------------------- package -------
library(Seurat)
library(patchwork)
library(tidyverse)
library(cowplot)
library(Matrix)
library(DoubletFinder)
library(scRepertoire)
library(harmony)
library(clustree)
#library(scPred)
library(RColorBrewer)
library(ggpubr)
library(eoffice)
#library(SeuratData)6+9*
library(Azimuth)
library(djvdj)
library(pheatmap)
library(paletteer)
library(ggnewscale)

source(file = '_src/Rfunctions_for_scRNA.R')

#----------------- project setting ------------
phe = read_tsv('sample_info.xls')


#------------------- color setting -------------
color_tf = structure(brewer.pal(9, 'Set1')[c(1,9)], names = c('TRUE', 'FALSE'))

color_clone1 = structure(brewer.pal(9, 'Set1')[c(1,3,9)], names = c('Dominant', 'Minor','NA'))

isotype_list = c('IGHD', 'IGHM', "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4")
col_isotype = structure( brewer.pal(8, 'Paired'),  names = isotype_list)

# order.cancer_type = c('MCL', 'DLBCL', 'Health')
# col.cancer_type = structure(  c(get_palette("Dark2", 2), "#868686FF" )  , names = order.cancer_type)

order.cancer_type = c('MCL', 'Health')
col.cancer_type = structure(  c(get_palette("Dark2", 1), "#868686FF" )  , names = order.cancer_type)

order.tissue = c(  "LN" , "BM", "Tonsil", "Pleural fluid",  "Spleen", "Intestine" )
col.tissue = structure( brewer.pal(6, 'Set2'), names = order.tissue)

color_cellcycle = structure( brewer.pal(9, 'Set1')[c(1,5,9)], names = c('G2M', 'S','G1'))


order.pr = c('Primary', 'Relapse', 'Normal'  )
col.pr = structure( get_palette("jco", 3), names = order.pr)

order.samtime = c('Normal', 'Primary', 'Pretreatment','Relapse'  )
col.samtime = structure( c("#868686FF","#0073C2FF", "#78B3E1", "#EFC000FF"), names = order.samtime)


order.cell_major = c("Malignant",  "B"  , "T" , "Myeloid",   "NK" ,  "Eryth" , "LMPP", "HSC", "Others" , "Doublets" )
col.cell_major = structure(c(brewer.pal(9, 'Set1'), 'grey80'), names = order.cell_major)

order.patient = c('MC3', 'MC202', 'MC203', 'MC205', 'MC206', 'MC201','MC4',  'MC1', 'MC2',  'MC204',  'MC207',   
                  "NC01", "NC02", "NC03","NC04", "NC05", "NC06",
                  "NC07", "NC08", "NC09")
col.patient = structure(c(brewer.pal(4, 'PuRd')[4],brewer.pal(12, 'Paired')[seq(2,13,2)] ,  brewer.pal(8, 'Dark2')[c(2,3,4,7)], rep("#666666",6),rep("#999999",3)),
                        names = order.patient)
col.donor = col.patient


order.sample = c(
                 "MC3P1", "MC3P2","MC3Rb", "MC3Ri", 
                 "MC202P",  "MC202R", 
                 "MC203P1" , "MC203P2" ,  
                 "MC205P",  "MC205R", 
                 "MC206P",  "MC206R", 
                 "MC201P1", "MC201P2",
                 "MC4R1","MC4R2",
                 "MC1R", "MC2P", "MC204P", "MC207R",
                 "NC01", "BM1","BM2","BM3" , "BM4", "BM5",
                 "RLN1", "RLN2","RLN3")
col.sample = structure(c( brewer.pal(4, 'PuRd'),
               brewer.pal(12, 'Paired') , # paired
               brewer.pal(8, 'Dark2')[c(2,3,4,7)], # single patietn
               rep("#666666",6),rep("#999999",3)),  # NC 
               names = order.sample)

order.treatment = c('ASCT', 'R-Chemo', 'Untreated')
col.treatment = structure(c('#a32a2a', '#2929a9', 'grey60'),  names = order.treatment)


col.dataset = structure(c('#6D8325FF',  '#BD5630FF',  '#8785B2FF'),
                        names = c('This study', 'PMID: 35549406', 'PMID: 3568781'))
col.dataset4 = structure(c('#6D8325FF',  '#BD5630FF',  '#8785B2FF', '#FF8C42FF'),
                        names = c('This study', 'PMID: 35549406', 'PMID: 3568781','Bone marrow with B cell enrichment'))

# comparisons
pairs_DR = c(  'MC202', 'MC205', 'MC206')
pairs_PR = c( 'MC3' )
pairs_PP = c('MC3', 'MC201', 'MC203')
pairs_RR = c('MC4')

#--------------------- Markers -------------------
markers.mcl         <- c("CCND1", "SOX11",  "TP53")
markers.B           <- c("CD19","MS4A1","CD79A")
markers.T            <- c("CD3D","CD3E","CD3G")
#markers.NK          <- c("XCL2", "NKG7", "GNLY")
markers.NK          <- c('KLRC1', 'KLRD1', 'KLRF1') 
markers.NKsub          <- c( 'KLRF1','NCAM1', 'FCGR3A') #NKp80 (KLRF1) CD56 (NCAM1) and CD16 (FCGR3A)

markers.Myeloid    <-c("CD68", "CD33","CST3")
markers.Mono        <- c('LYZ', 'S100A9', 'FCN1')
markers.macro       <- c("LYZ","C1QA","C1QB")
markers.pDC         <- c("NRP1","CLEC4C")
markers.Eryth       <- c('HBQ1', 'HBM', 'GYPA')
markers.mag         <- c('GATA2', 'FCER1A', 'GATA1')
markers.Eryth.early <- c('APOC1', 'TESPA1', 'GATA2')
markers.Eryth.late <- c('GYPA', 'GYPB', 'HBA1')
markers.NMP <- c('LYZ', 'MPO', 'SERPINB1')
markers.ELP <- c('FLT3', 'LTB', 'RUNX2') #early lymphoid progenitors??
markers.T.sub            <- c("CD3E","CD4","CD8A")
makrers.prolif    <- c('TOP2A', 'TU88', 'MKI67', 'UBE2C')
makrers.epith  <- c('EPCAM', 'KRT18', 'KRT19') #epithelial cells

# ------------------ functions ----------------------
plot_package <- function(data.filt, mydir){
  
  dir.create(mydir)
  
  # feature plot
  plot_grid(ncol = 3, 
            FeaturePlot(data.filt, features = "nFeature_RNA"),
            FeaturePlot(data.filt, features = "percent_hb"),
            FeaturePlot(data.filt, features = "percent_mito"),
            FeaturePlot(data.filt, features = "percent_ribo"),
            FeaturePlot(data.filt, features = "riboLT5"),
            DimPlot(data.filt, group.by = 'Phase'))
  ggsave(last_plot(), filename = paste0(mydir, 'feature_plot.raw.pdf'), h = 8, w = 12)
  
  # Groups
  plot_grid(ncol = 2,
            DimPlot(data.filt, group.by = "Cancer_type", reduction = 'umap', cols = col.cancer_type) +  ggtitle('Cancer type'),  
            #DimPlot(data.filt, group.by = "Batch", reduction = 'umap', cols = get_palette('Set2', n_distinct(phe$Batch))) +  ggtitle('umap'),
            DimPlot(data.filt, group.by = "Tissue", reduction = 'umap', cols = col.tissue) +  ggtitle('Tissue'),
            DimPlot(data.filt, group.by = "Donor", reduction = 'umap', cols = col.donor)+  ggtitle('Donor'), 
            DimPlot(data.filt, group.by = "Sampling_time", reduction = 'umap', cols = col.samtime) +  ggtitle('Sample')
            
  )
  
  ggsave(last_plot(), filename = paste0(mydir, 'dimplot.Groups.pdf'), h = 16, w = 24)
  
  #split samples
  
  DimPlot(data.filt, group.by = "orig.ident", split.by ="orig.ident" , reduction = 'umap', ncol = 6) +  ggtitle('umap') & NoLegend()
  ggsave(last_plot(), filename = paste0(mydir, '/dimplot.SampleSplit.pdf'), h = 20, w = 20)
  
  
  # markers
  FeaturePlot(data.filt, reduction = 'umap', dims = 1:2, 
              features = c(markers.B,markers.T.sub ,markers.NK,markers.Myeloid , markers.Eryth,markers.mcl),
              ncol = 6,   order = T) & NoLegend()
  
  ggsave(last_plot(), filename =  paste0(mydir,'/featureplot.markers.pdf'), h = 10, w= 24)
  
  
}

my.featureplot <- function(data, gene, mydir, w = 6, h = 6) {
  ggsave( FeaturePlot(data, reduction = "umap", features = gene, raster = T) &
            NoAxes() &  NoLegend() & 
            labs(title = gene) &
            # scale_color_manual(values = color_tf) +
            theme(panel.background = element_rect(colour = 'black', size = 2)) ,
          filename =  paste0(mydir, '/mark.', gene, '.png'), w = w, h = h)
  
}


my.Vlnplot <- function(data, gene, x, color){
  VlnPlot(data, features = gene,   pt.size = 0, group.by =  x,  stack = F, same.y.lims = F) & NoLegend()  & 
    stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)  &
    scale_fill_manual(values = color)
} 


plot_umap <- function(q.p, mydir, w = 6, h = 6, reduction = 'umap', sel.clust = "SCT_snn_res.0.5"){
  # Dimplot
  # w = 6
  # h = 6
  # reduction = 'umap'
  # sel.clust = "SCT_snn_res.0.5"
  
  # by samples
  ggsave( DimPlot(q.p, reduction = reduction, group.by = "Sample", repel = T, label = F, label.size = 7) + NoAxes()  +
            labs(title = paste0('by Sample') ) +
            scale_color_manual(values = col.sample) +
            theme(panel.background = element_rect(colour = 'black', size = 2))  ,
          filename =  paste0(mydir, reduction, '.Sample.T.png'), w = w+2, h = h)
  
  ggsave( DimPlot(q.p, reduction = reduction, group.by = "Sample", repel = T, label = T, label.size = 7) + NoAxes()  +
            labs(title = paste0('by Sample') ) +  NoLegend() +
            scale_color_manual(values = col.sample) +
            theme(panel.background = element_rect(colour = 'black', size = 2))  ,
          filename =  paste0(mydir,reduction, '.Sample.T.L.png'), w = w, h = h)
  
  
  # by cell cycle
  ggsave( DimPlot(q.p, reduction = reduction, group.by = 'Phase', repel = T, label = F, label.size = 7) + NoAxes()  +
            labs(title = paste0('by Cell cycle') ) +
            scale_color_manual(values = color_cellcycle) +
            theme(panel.background = element_rect(colour = 'black', size = 2))  ,
          filename =  paste0(mydir,reduction, '.CellCycle.T.png'), w = w+1, h = h)
  
  
  # plot igK/igL
  ggsave(FeaturePlot(q.p, reduction = reduction,  features  = "igl_kl_ratio", repel = T,
                     raster = F, label = F, label.size = 7) + NoAxes()  +
           scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")), limits = c(0,1)) +
           labs(title = paste0('by IGK/IKL ratio') ) +
           theme(panel.background = element_rect(colour = 'black', size = 2))  ,
         filename =  paste0(mydir,reduction, '.IgK_L.T.png'), w = w+1, h = h)
  
  # plot CNVscore
  ggsave( FeaturePlot(q.p, reduction = reduction, features = "CNV_score", repel = T,
                      label = F, label.size = 7, raster = T) + NoAxes()  +
            labs(title = paste0('by CNV score') ) +
            #scale_colour_gradient2(midpoint = CNV.cutoff , limits = c(0, CNV.max) , low =   "#438ec2", mid = "white", high =  "#c2438e") +
            theme(panel.background = element_rect(colour = 'black', size = 2))  ,
          filename =  paste0(mydir,reduction, '.CNVscore.T.png'), w = w+1, h = h)
  
  
  
  # plot SHM
  ggsave( FeaturePlot(q.p, reduction = reduction, features  = "SHM_VH", repel = T,
                      max.cutoff = quantile(q.p$SHM_VH, na.rm = T, probs = 0.99),
                      # cols =  c("#00AFBB", "#E7B800", "#FC4E07"),
                      cols =  c( "#00AFBB", "#FC4E07"),
                      label = F, label.size = 7) + NoAxes()  +
            #scale_color_manual(values = color_clone1) +
            labs(title = paste0('by SHM') ) +
            theme(panel.background = element_rect(colour = 'black', size = 2))  ,
          filename =  paste0(mydir,reduction, '.VDJ.T.png'), w = w+1, h = h)
  
  
  # plot Cluster
  ggsave( DimPlot(q.p, reduction = reduction, group.by = sel.clust, repel = T, label = T, label.size = 7) + NoAxes()  +
            labs(title = paste0('by Cluster') ) +
            # scale_color_manual(values = color_cellcycle) +
            theme(panel.background = element_rect(colour = 'black', size = 2))  ,
          filename =  paste0(mydir,reduction, '.Cluster.T.png'), w = w, h = h)
  
  
  
  
}


## enrichment getsets
stromal_gs <- read.gmt("_ref/Signatures from Lenz 2008 NEJM.gmt.txt")
kegg_gs <- read.gmt("_ref/c2.cp.kegg.v2023.1.Hs.symbols.gmt")
hm_gs <- read.gmt("_ref/h.all.v2023.1.Hs.symbols.gmt")
react_gs <- read.gmt("_ref/c2.cp.reactome.v2023.2.Hs.symbols.gmt")

## checkpoints
gene.ckp.T <- c('CTLA4', 'TIGIT', 'PDCD1', 'ICOS', 'TNFRSF18', 'TNFRSF4', 'LAG3', 'TNFRSF9', 'CD40LG', 'CD27', 'HAVCR2', 'CD274')
gene.ckp.B <- c('CD70',  'CD80', 'CD86', 'PD1', 'CD40', 'ICAM1','CD274')


# color for TME
col.mo <- structure(scales::hue_pal()(7), names = c('CD14 Mono', 'GMP', 'CD16 Mono', 'pDC', 'cDC', 'GMP-prolif','Macro'))
color.T1 = structure(c("#8DD3C7", "#FDB462","#BEBADA" ,"#FB8072"), names = c("CD4", "CD8", 'gdT', 'NK'))
color.T2 = structure( c(paletteer_d('pals::tol', n = 12)), 
                      names = c("CD4.Tm" , "CD8.Tn" , "CD56.NK", "CD4.Tfh" ,"CD4.Tn" , "CD8.Teff" ,"CD8.Tex"  ,"CD4.Treg",  "CD16.NK",    "CD8.Tm",  "CD8.Tprolif" ,  "gdT"  ))

# color.T3 = structure( c(paletteer_d("pals::stepped", n = n_distinct(q.t$celltype.T.l3))), names = levels(q.t$celltype.T.l3))


color.major_celltype.l2 <- structure(c("#E41A1C", "#377EB8", "#4DAF4A","#8DD3C7", "#FDB462","#BEBADA" , "#FF7F00",  "#984EA3",  "grey80" ), 
                                     names = c("Malignant",   "B", "T" ,  "CD4",  "CD8", "gdT" ,  "NK" ,    "Myeloid" ,      "Others"   ))
