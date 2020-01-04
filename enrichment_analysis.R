############################################################################################################################################

# Analysis for Carolina's paper

############################################################################################################################################
setwd('~/Documents/GeneExprBlood_Psych_scripts/')

library(data.table)
library(readxl)
source('./fisher_overlap.R')
library(ggplot2)

genesetanalysis = read_xlsx('Cell_type_fordalila.xlsx', sheet=2)
celltype = read_xlsx('Cell_type_fordalila.xlsx', sheet=3)
background = read_xlsx('Cell_type_fordalila.xlsx', sheet=4)

lake = transpose(read.table('gene sets/list_celltypes_lake.txt', fill=TRUE,sep='\t',col.names = paste0("V",seq_len(841))))
galatro = transpose(read.table('gene sets/list_celltypes_galatro2017-1296g_microglia.txt', fill=TRUE,sep='\t'))

setnames(lake,unlist(lake[1,]))
setnames(galatro,unlist(galatro[1,]))

lake = lake[-1,]
darmanislake = darmanislake[-1,]
velmesh = velmesh[-1,]
galatro = galatro[-1,]

## Part A. Cell Type enrichments
enrichmentCellType = data.frame()
genelists = c('Blood','Bottlenecks','GandalFDR<0.05 cross disorders ∩ Blood',
              'Gandaloverlap down-regulated genes across disorders ∩ Blood',
              'Gandaloverlap up-regulated genes across disorders ∩ Blood',
              'GWAS∩Blood','Blood∩Exome')
for(i in genelists) {
  print(i)
  genes = celltype[,i][[1]]
  genes = genes[!is.na(genes)]
  # Lake
  for(c in colnames(lake)) {
    f = as.numeric(ORA(genes, lake[,c], background$'genesymbol_cortex&blood>=0.5', background$'genesymbol_cortex&blood>=0.5','greater'))
    enrichmentCellType = rbind(enrichmentCellType, 
                        data.frame(GeneList=i,Dataset="Lake", CellType=c, OR=f[[1]], p=f[[2]]))
  }
  # Galatro
  for(c in colnames(galatro)) {
    f = as.numeric(ORA(genes, galatro[,c], background$'genesymbol_cortex&blood>=0.5', background$'genesymbol_cortex&blood>=0.5','greater'))
    enrichmentCellType = rbind(enrichmentCellType, 
                               data.frame(GeneList=i,Dataset="galatro", CellType=c, OR=f[[1]], p=f[[2]]))
  }
}

# FDR correction
# Method 1: FDR correct in bulk (NOTE: don't use this method, see below instead)
#enrichmentCellType$p[enrichmentCellType$OR < 1] = 1  # This line was included in Mike's original code used for CD paper
# enrichmentCellType$fdr = p.adjust(enrichmentCellType$p, method = "fdr")
# enrichmentCellType$log10fdr = -log10(enrichmentCellType$fdr) 
# enrichmentCellType$log10fdr[is.infinite(enrichmentCellType$log10fdr)] = max(enrichmentCellType$log10fdr[!is.infinite(enrichmentCellType$log10fdr)])
# enrichmentCellType$text = "*"

# Method 2: FDR correct by gene list (NOTE: use this method)
for (g in genelists){
  enrichmentCellType$fdr[enrichmentCellType$GeneList==g] = p.adjust(enrichmentCellType$p[enrichmentCellType$GeneList==g], method = "fdr")
  enrichmentCellType$log10fdr[enrichmentCellType$GeneList==g] = -log10(enrichmentCellType$fdr[enrichmentCellType$GeneList==g]) 
  enrichmentCellType$log10fdr[is.infinite(enrichmentCellType$log10fdr[enrichmentCellType$GeneList==g])] = max(enrichmentCellType$log10fdr[!is.infinite(enrichmentCellType$log10fdr[enrichmentCellType$GeneList==g])])
#  enrichmentCellType$text[enrichmentCellType$GeneList==g] = as.character(signif(enrichmentCellType$OR[enrichmentCellType$GeneList==g],2))
}

# Only display OR for signficant results
enrichmentCellType$text=""
enrichmentCellType$text[enrichmentCellType$fdr < 0.05] = as.character(signif(enrichmentCellType$OR[enrichmentCellType$fdr<0.05],2))

# See significant rows @ FDR <.05
enrichmentCellType[enrichmentCellType$fdr<0.05,]

write.table(enrichmentCellType,'CellType_enrichmentResults_LakeGalatroOnly_FDRcorrectByGeneList.csv',sep=',',quote=F,row.names=T)

## Part B. Tissue enrichments
enrichmentGeneSets = data.frame()
genesetlists = c('Gandal FDR<0.05 cross disorders',
  'Gandal overlap down-regulated genes across disorders','Gandal overlap up-regulated genes across disorders',
  'GWAS','Exome LGD','Exome LGD & missense')
for (i in c('Blood','TF')){
  for(j in genesetlists) {
    genes = genesetanalysis[,i][[1]]
    genes = genes[!is.na(genes)]
    setgenes = genesetanalysis[,j][[1]]
    setgenes = setgenes[!is.na(setgenes)]
    f = as.numeric(ORA(genes, setgenes, background$'genesymbol_cortex&blood>=0.5', background$'genesymbol_cortex&blood>=0.5','greater'))
    enrichmentGeneSets = rbind(enrichmentGeneSets, 
                               data.frame(GeneList=i,Dataset=j, OR=f[[1]], p=f[[2]]))
  }
}

# FDR correction
# Method 1: FDR correct in bulk
#enrichmentGeneSets$p[enrichmentGeneSets$OR < 1] = 1  # This line was included in Mike's original code used for CD paper
# enrichmentGeneSets$fdr = p.adjust(enrichmentGeneSets$p, method = "fdr")
# enrichmentGeneSets$log10fdr = -log10(enrichmentGeneSets$fdr) 
# enrichmentGeneSets$log10fdr[is.infinite(enrichmentGeneSets$log10fdr)] = max(enrichmentGeneSets$log10fdr[!is.infinite(enrichmentGeneSets$log10fdr)])
# enrichmentGeneSets$text = as.character(signif(enrichmentGeneSets$OR,2))

# Method 2: FDR correct by gene list (NOTE: use this method)
for (g in c('Blood','TF')){
  enrichmentGeneSets$fdr[enrichmentGeneSets$GeneList==g] = p.adjust(enrichmentGeneSets$p[enrichmentGeneSets$GeneList==g], method = "fdr")
  enrichmentGeneSets$log10fdr[enrichmentGeneSets$GeneList==g] = -log10(enrichmentGeneSets$fdr[enrichmentGeneSets$GeneList==g]) 
  enrichmentGeneSets$log10fdr[is.infinite(enrichmentGeneSets$log10fdr[enrichmentGeneSets$GeneList==g])] = max(enrichmentGeneSets$log10fdr[!is.infinite(enrichmentGeneSets$log10fdr[enrichmentGeneSets$GeneList==g])])
  enrichmentGeneSets$text[enrichmentGeneSets$GeneList==g] = as.character(signif(enrichmentGeneSets$OR[enrichmentGeneSets$GeneList==g],2))
}

# See significant rows @ FDR <.05
enrichmentGeneSets[enrichmentGeneSets$fdr<0.05,]

# Only display OR for signficant results
enrichmentGeneSets$text=""
enrichmentGeneSets$text[enrichmentGeneSets$fdr < 0.05] = as.character(signif(enrichmentGeneSets$OR[enrichmentGeneSets$fdr<0.05],2))

write.table(enrichmentGeneSets,'GeneSet_enrichmentResults_FDRcorrectAll.csv',sep=',',quote=F,row.names=T)

## Part C. Make plots
# Cell type enrichments
enrichmentCellType_lake = enrichmentCellType[enrichmentCellType$Dataset=='Lake',]

g = ggplot(enrichmentCellType_lake, aes(x=CellType,y=GeneList, label=text)) +
  geom_tile(aes(fill=log10fdr),color="grey60") + 
  coord_fixed(ratio=2.5) +
  scale_fill_gradient(low = "white", high = "red","-log10 FDR", limits = c(0,3)) + 
  geom_text(size=3, color="black") + 
  ylab("") + 
  xlab("") + 
  theme(
    #plot.margin = margin(1,.8,1,.8, "cm"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1), 
    panel.background = element_blank(),
    plot.background = element_blank()
)

aspect_ratio = 1.8
height = 7
ggsave('cellTypeEnrichment_Lake.pdf',plot=g,height=height, width=height*aspect_ratio)

#Cell-Type enrichment
#g5=ggplot(enrichmentCellType[enrichmentCellType$Module==moduleColor,],aes(x=CellType,y=log10fdr, fill=Dataset)) + geom_bar(stat="identity") + coord_flip() + geom_hline(yintercept = -log10(0.05),lty=2) + xlab("") + ggtitle("Cell-Type Enrichment") +theme(axis.text.y = element_text(size=6))

# Blood/TF enrichments
g2 = ggplot(enrichmentGeneSets, aes(x=Dataset,y=GeneList, label=text)) +
  geom_tile(aes(fill=log10fdr),color="grey60") + 
  coord_fixed(ratio=2.5) +
  scale_fill_gradient(low = "white", high = "mediumpurple4","-log10 FDR", limits = c(0,5)) + 
  geom_text(size=3, color="black") + 
  ylab("") + 
  xlab("") + 
  theme(
    #plot.margin = margin(1,.8,1,.8, "cm"),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1), 
    panel.background = element_blank(),
    plot.background = element_blank()
  )

aspect_ratio = 1.8
height = 7
ggsave('BloodTFEnrichment.pdf',plot=g2,height=height, width=height*aspect_ratio)
