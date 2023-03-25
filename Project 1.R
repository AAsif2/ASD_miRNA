library('purrr')
library('ggplot2')
library('ggpubr')
library('VennDiagram')
library('RColorBrewer')
library('gsubfn')
library('dplyr')

Wucortex<-read.table('Wu cortex.csv', header=TRUE, sep=',')

Wucerebellum<-read.table('Wu cerebellum.csv', header=TRUE, sep=',')

Ander<-read.table('Ander.csv', header=TRUE, sep=',')

Mor<-read.table('Mor1.csv', header=TRUE, sep=',')



#download list of validated miRNAs from miRBase
alias<- read.table('aliases 2.txt', header=FALSE, sep='\t')

#rename alias columns
colnames(alias)[1] = 'miRNA'
colnames(alias)[2] = 'aliases'

#keep only human miRs
alias<- alias[grep("hsa", alias$aliases), ]

#capitalise mir>miR names in alias
alias$aliases = gsub('mir','miR',alias$aliases)



#add ID column to Wucortex
Wucortex$ID <- NA

#add ID names to Wucortex 
for(i in Wucortex$miRNA){
  subset = alias[grep(i, alias$aliases),]
  ID = subset[1,1]
  Wucortex$ID[Wucortex$miRNA==i] = ID
}

#replace ID for miR-1 to correct ID
Wucortex$ID[Wucortex$miRNA == 'hsa-miR-1']<-'MIMAT0000416'


# add ID column to Wucerebellum
Wucerebellum$ID <- NA

#add iD names to Wucerebellum
for(i in Wucerebellum$miR.ID){
  subset = alias[grep(i, alias$aliases),]
  ID = subset[1,1]
  Wucerebellum$ID[Wucerebellum$miR.ID==i] = ID
}

#replace ID for miR-1 to correct ID
Wucerebellum$ID[Wucerebellum$miR.ID == 'hsa-miR-1']<-'MIMAT0000416'



# add ID column to Ander
Ander$ID <- NA

#add iD names to Ander
for(i in Ander$miRNA){
  subset = alias[grep(i, alias$aliases),]
  ID = subset[1,1]
  Ander$ID[Ander$miRNA== i] = ID
}

#replace NA for miR-4709-3p to known ID
Ander$ID[which(is.na(Ander$ID))]<-'MI0019812'
#replace ID for miR-1 to correct ID
Ander$ID[Ander$miRNA == 'miR-1']<-'MIMAT0000416'


#add ID column to Mor
Mor$ID <- NA

#add ID names to Mor
for(i in Mor$miRNA){
  subset = alias[grep(i, alias$aliases),]
  ID = subset[1,1]
  Mor$ID[Mor$miRNA==i] = ID
}



#subset for upregulated Wucortex
Wucortexup <- subset(Wucortex, log2.fold.change. > 0 & FDR.adjusted.Pval.diagnosis <= 0.05,
                     select=c(ID))

#subset for downregulated Wucortex
Wucortexdown <- subset(Wucortex, log2.fold.change. < 0 & FDR.adjusted.Pval.diagnosis <= 0.05,
                     select=c(ID))

#subset for all dysregulated Wucortex
Wucortexdys <- subset(Wucortex, log2.fold.change. !=0 & FDR.adjusted.Pval.diagnosis <= 0.05,
                       select=c(ID))


#subset for upreg Wucerebellum
Wucerebellumup <- subset(Wucerebellum, Beta.diagnosis..log2.fold.change..ASD.vs..CTL. > 0 & FDR.adjusted.Pval.diagnosis <= 0.05,
                       select=c(ID))

#subset for downregulated Wucerebellum
Wucerebellumdown <- subset(Wucerebellum, Beta.diagnosis..log2.fold.change..ASD.vs..CTL. < 0 & FDR.adjusted.Pval.diagnosis <= 0.05,
                       select=c(ID))

#subset for all dysregulated Wucerebellum
Wucerebellumdys <- subset(Wucerebellum, Beta.diagnosis..log2.fold.change..ASD.vs..CTL. != 0 & FDR.adjusted.Pval.diagnosis <= 0.05,
                          select=c(ID))


#subset for upreg Ander (don't need to subset for p.value as list is already confirmed stat signif DE miRs)
Anderup <- subset(Ander, Fold.change > 0,
                         select=c(ID))

#subset for downregulated Ander (don't need to subset for p.value as list is already confirmed stat signif DE miRs)
Anderdown <- subset(Ander, Fold.change < 0,
                           select=c(ID))

#subset for dysregulated Ander (don't need to subset for p.value as list is already confirmed stat signif DE miRs)
Anderdys <- subset(Ander, Fold.change !=0,
                    select=c(ID))


#subset for upreg Mor (don't need to subset for p.value as list is already confirmed stat signif DE miRs)
Morup <- subset(Mor, Fold.change > 1,
                  select=c(ID))

#subset for downregulated Mor (don't need to subset for p.value as list is already confirmed stat signif DE miRs)
Mordown <- subset(Mor, Fold.change < 1,
                    select=c(ID))
Mordown <- Mordown$ID

#subset for dysregulated Mor (don't need to subset for p.value as list is already confirmed stat signif DE miRs)
Mordys <- subset(Mor, Fold.change !=0,
                   select=c(ID))



#Remove NAs due to experimentally unvalidated miRNAs
Wucortexup= Wucortexup[complete.cases(Wucortexup), ]
Wucerebellumup= Wucerebellumup[complete.cases(Wucerebellumup), ]
Anderup= Anderup[complete.cases(Anderup), ]
Morup= Morup[complete.cases(Morup), ]

Wucortexdown= Wucortexdown[complete.cases(Wucortexdown), ]
Wucerebellumdown= Wucerebellumdown[complete.cases(Wucerebellumdown), ]
Anderdown= Anderdown[complete.cases(Anderdown), ]

Wucortexdys= Wucortexdys[complete.cases(Wucortexdys),]
Wucerebellumdys= Wucerebellumdys[complete.cases(Wucerebellumdys),]
Anderdys= Anderdys[complete.cases(Anderdys),]
Mordys= Mordys[complete.cases(Mordys), ]


Wucerebellum_list = as.character(Wucerebellum$ID)

Wucerebellum_list = Wucerebellum_list[complete.cases(Wucerebellum_list)]



# Venn diagram upregulated miRNAs
venn.diagram(
  x = list(Wucortexup, Wucerebellumup, Anderup, Morup),
  main = "Number of significantly upregulated miRNAs",
  main.pos = c(0.5, 0.72),
  main.fontface = "bold",
  category.names = c("Wu - cerebrum", "Wu - cerebellum", "Ander - cerebrum", "Mor - cerebrum"),
  lty = rep ('blank'), 
  fill = c("#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854"),
  margin = 0.6,
  filename = '#Upreg_brain_miRNA.png',
  output=TRUE
  

)


#Upregulated miRNAs 
Upregulated<-list(Wucortexup, Wucerebellumup, Anderup, Morup)


#Identify overlapping upregualted miRNAs
Wuupoverlap<- calculate.overlap(x=list('Wucortexup'=Wucortexup,'Wucerebellumup'=Wucerebellumup))
WucerebellumMorupoverlap <- calculate.overlap(x=list('Wucerebellumup'=Wucerebellumup, 'Morup'=Morup))
WucortexMorupoverlap <- calculate.overlap(x=list('Wucortexup'=Wucortexup, 'Morup'=Morup))

#Fisher's matrix for Wucortex and Wucerebllum
Wuupreg_Fishers_table = matrix(c(4,43,28,439), nrow = 2, ncol = 2, dimnames = list(c('Cortex yes','Cortex no'), c('Cerebellum yes','Cerebellum no')))
fisher.test(Wuupreg_Fishers_table)



# Venn diagram downregulated miRNAs
venn.diagram(
  x = list(Wucortexdown, Wucerebellumdown, Anderdown, Mordown),
  main = "Number of significantly downregulated miRNAs",
  main.pos = c(0.5, 0.72),
  main.fontface = "bold",
  category.names = c("Wu - cerebrum", "Wu - cerebellum", "Ander - cerebrum", "Mor - cerebrum"), 
  lty = rep ('blank'), 
  fill = c("#66C2A5", "#8DA0CB", "#E78AC3", "#A6D854"),
  margin = 0.6,
  filename = '#Downreg_brain_miRNA.png',
  output=TRUE
  
  
)


#Downregulated miRNAs
Downregulated<-list(Wucortexdown, Wucerebellumdown, Anderdown, Mordown)

#Identify overlapping downregulated miRNAs
Wudownoverlap<- calculate.overlap(x=list('Wucortexdown'=Wucortexdown,'Wucerebellumdown'=Wucerebellumdown))


#Fisher's matrix for downreg Wucortex and Wucerebllum
Wudownreg_Fishers_table = matrix(c(1,22,4,487), nrow = 2, ncol = 2, dimnames = list(c('Cortex yes','Cortex no'), c('Cerebellum yes','Cerebellum no')))
fisher.test(Wudownreg_Fishers_table)



# Venn diagram dysregulated miRNAs
venn.diagram(
  x = list(Wucortexdys, Wucerebellumdys, Anderdys, Mordys),
  main = "Number of significantly dysregulated miRNAs",
  main.pos = c(0.5, 0.72),
  main.fontface = "bold",
  category.names = c("Wu - cerebrum", "Wu - cerebellum", "Ander - cerebrum", "Mor - cerebrum"),
  lty = rep ('blank'), 
  fill = c("#66C2A5",  "#8DA0CB", "#E78AC3", "#A6D854"),
  margin = 0.6,
  filename = '#Dysreg_brain_miRNA.png',
  output=TRUE
  
  
)


#Dysregulated miRNAs
Dysregulated<-list(Wucortexdys, Wucerebellumdys, Anderdys, Mordys)

#Identify overlapping dysregulated miRNAs
Wudysopverlap <- calculate.overlap(x=list('Wucortexdys'=Wucortexdys,'Wucerebellumdys'=Wucerebellumdys))
WucortexMordysoverlap <- calculate.overlap(x=list('Wucortexdys'=Wucortexdys, 'Mordys'=Mordys))
WucerebellumMordysoverlap <- calculate.overlap(x=list('Wucerebellumdys'=Wucerebellumdys, 'Mordys'=Mordys))

  
#Fisher's matrix for Wucortex and Wucerebellum
Wudysreg_Fishers_table = matrix(c(5,65,32,412), nrow = 2, ncol = 2, dimnames = list(c('Cortex yes','Cortex no'), c('Cerebellum yes','Cerebellum no')))
fisher.test(Wudysreg_Fishers_table)


#Fisher's matrix for Wucortex and Wucerebellum without downreg miRNA
Wudysreg2_Fishers_table = matrix(c(4,65,32,413), nrow = 2, ncol = 2, dimnames = list(c('Cortex yes','Cortex no'), c('Cerebellum yes','Cerebellum no')))
fisher.test(Wudysreg2_Fishers_table)



#create scatterplot for all miRNAs in Wucortex + Wucerebellum (Beta.diagnosis..log2.fold.change..ASD.vs..CTL.)
Wu<- merge(Wucortex, Wucerebellum, by = 'ID', all=TRUE) 
Wuplot_df <- unique(Wu[,c(1,3,4,7,8)])
Wuplot_df= Wuplot_df[complete.cases(Wuplot_df),]  

Wuplot_df$Statistically.significant<- "Neither study"
Wuplot_df$Statistically.significant [Wuplot_df$FDR.adjusted.Pval.diagnosis.x<0.05]<- "Wu cerebrum"
Wuplot_df$Statistically.significant [Wuplot_df$FDR.adjusted.Pval.diagnosis.y<0.05]<- "Wu cerebellum"
Wuplot_df$Statistically.significant [Wuplot_df$FDR.adjusted.Pval.diagnosis.x<0.05 & Wuplot_df$FDR.adjusted.Pval.diagnosis.y<0.05]<- "Wu cerebrum & Wu cerebellum"


Wuplot<- ggplot(Wuplot_df, aes(x=log2.fold.change., y=Beta.diagnosis..log2.fold.change..ASD.vs..CTL., color = Statistically.significant)) +
  labs(x="log2 fold change in Wu cerebrum", y= "log2 fold change in Wu cerebellum",
       title="Fold change in miRNA expression in ASD subjects compared to controls", color = "p < 0.5")+
  geom_point() +scale_color_manual(breaks = c("Wu cerebrum", "Wu cerebellum", "Wu cerebrum & Wu cerebellum", "Neither study"), values=c("seagreen1","steelblue1","red","gray54"))

Wuplot + scale_x_continuous(limits = c(-1,1))+
  scale_y_continuous(limits = c(-1,1))+
theme_bw() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  stat_cor(method = "pearson", label.x = 3, label.y = 30)


#dataframe holding dysregulated miRNAs and their IDs
dysregulatedmiRNA_df <- data.frame (ID = c(Wudysopverlap$a3, WucerebellumMordysoverlap$a3, WucortexMordysoverlap$a3),
                                    miRNA = c(NA)
)

for(i in dysregulatedmiRNA_df$ID){
  subset = alias[grep(i, alias$miRNA),]
  miRNA = subset[1,2]
  dysregulatedmiRNA_df$miRNA[dysregulatedmiRNA_df$ID==i] = miRNA
}

dysregulatedmiRNA_df<- dysregulatedmiRNA_df[!grepl("MI0016088", dysregulatedmiRNA_df$ID),]

#NOT WORKING splitting dyregulatedmiRNA_df column into multiple columns using pattern matching so can select mirs to put into dysreg mir vector
#NOT WORKING  miRNApattern<- "(('.*'), )?(('.p'),)?"
#NOT WORKING  r <- read.pattern(text = as.character(dysregulatedmiRNA_df$miRNA), pattern = miRNApattern, as.is = TRUE)
#NOT WORKING  dysregulatedmiRNA_df2 <- cbind(miRNA = dysregulatedmiRNA_df$miRNA, r[c(FALSE, TRUE)], stringsAsFactors = FALSE)
#NOT WORKING  nc <- ncol(dysregulatedmiRNA_df2) 
#NOT WORKING  names(dysregulatedmiRNA_df2)[-1] <- paste0("miRNA_", 1:(nc-1))

#easier way to creat dysreg mir vector which works
dysregulatedmiRNAs<- c("hsa-miR-130b-5p","hsa-miR-148a-3p", "hsa-miR-3938", "hsa-miR-425-3p", "hsa-miR-19b-3p", "hsa-miR-155-5p", "hsa-miR-21-3p")


#Finding SFARI genes in miRNA target lists obtained from mienturnet

#read files
SFARI<-read.table('SFARI.csv', header=TRUE, sep=',')

miRNAtargetsmiRTB <- read.table('Mient_miRTarBase1Target.csv', header=TRUE, sep=',')

miRNAtargetsTS <- read.table('Mient_TS1Target.csv', header=TRUE, sep=',')

#subset target data to find stat. signif results
#subset miRTarbase genes for miRNA targets with FDR <0.05
miRNAtargetsmiRTB_0.05 <- subset(miRNAtargetsmiRTB, FDR <= 0.05,)

#subset miRTarbase genes for miRNA targets with p-value <0.1
miRNAtargetsmiRTB_0.1 <- subset(miRNAtargetsmiRTB, FDR < 0.1,)

#subset miRTarbase genes for miRNA targets with p-value <0.2
miRNAtargetsmiRTB_0.2 <- subset(miRNAtargetsmiRTB, FDR < 0.2,)

#subset TargetScan genes for miRNA targets with p-value <0.05
miRNAtargetsTS_0.05 <- subset(miRNAtargetsTS, FDR <= 0.05,)

#subset TargetScan genes for miRNA targets with p-value <0.2
miRNAtargetsTS_0.2 <- subset(miRNAtargetsTS, FDR <= 0.2,)


#Identify which miRTarBase genes targeted by the miRNAs are in SFARI
miRNAtargetsmiRTB_0.05SFARI <- subset(miRNAtargetsmiRTB_0.05, miRNAtargetsmiRTB_0.05$Gene.Symbol %in% SFARI$gene.symbol)

miRNAtargetsmiRTB_0.1SFARI <- subset(miRNAtargetsmiRTB_0.1, miRNAtargetsmiRTB_0.1$Gene.Symbol %in% SFARI$gene.symbol)

miRNAtargetsmiRTB_0.2SFARI <- subset(miRNAtargetsmiRTB_0.2, miRNAtargetsmiRTB_0.2$Gene.Symbol %in% SFARI$gene.symbol)

#Identify which TargetScan genes targeted by the miRNAs are in SFARI
miRNAtargetsTS_0.2SFARI <- subset(miRNAtargetsTS_0.2, miRNAtargetsTS_0.2$Gene.Symbol %in% SFARI$gene.symbol)



#Identify which genes the miRNAs target in Nowakowski
miRNANowakowski_df <- Nowakowski[Nowakowski$miR %in% miRNAtargets$microRNA.1 | Nowakowski$miR %in% miRNAtargets$microRNA.2 | Nowakowski$miR %in% miRNAtargets$microRNA.3 | Nowakowski$miR %in% miRNAtargets$microRNA.4 | Nowakowski$miR %in% miRNAtargets$microRNA.5 | Nowakowski$miR %in% miRNAtargets$microRNA.6, ]



for (i in dysregulatedmiRNA_df$ID){
  ID = i #just a label to remind what i is
  miR_targets = subset(alluniquemiRNAtargets, miRTarBase.ID == ID)
  for (j in scdata$Cluster){
    cluster = j  
    cell_type_markers = subset(scdata, cluster == Celltype)
    A = subset(scdata$Gene %in% alluniquemiRNAtargets$Target.Gene)
    A_size = length(unique(A$Gene))
    B = length(alluniquemiRNAtargets$Target.Gene) - A_size
    C = length(scdata$Gene) - A_size
    D = length(allhumangenes$Gene.name)- A_size - B - C
    celltype_fishersmatrix = matrix(c(A,B,C,D), nrow=2, ncol=2, dimnames = list(c('Cell marker', 'Not a cell marker'), c('Target gene', 'Not a target gene')))
    celltype_fishersresult = fisher.test(celltype_fishersmatrix)
    my_output <- data.frame(miRNA = ID, CellType = cluster,
                            pvalue=fishers_result$p.value[1], OddsRatio=fishers_result$estimate[1])
    celltype_fishersresultsall = rbind(celltype_fishersreults, my_output)
  }
}



#find cells ALL miRNA targets are enriched in, miRTarBase


allhumangenes<- read.table('human protein coding genes.txt', header=TRUE, sep='\t')
scdata<- read.table('Lake B.csv', header=TRUE, sep=',')[,c(1,4,7)]
colnames(scdata) <- c("Gene","logFC","Cluster")
AllmiRTarBase<-read.table('mirtarbase.csv', header=TRUE, sep=',')[,c(1,2,4,7,8)]

All_miRNAs<- unique(AllmiRTarBase$miRNA)
AllmiRTarBase_Targets <- unique(AllmiRTarBase$Target.Gene)

#for each cell type in sc data with different versions, reduce them down to their main type e.g. purk cell 1 and 2 both become purk cell 
excitatory_neurons <- unique(scdata$Cluster)[1:13]
granule_cell <- unique(scdata$Cluster) [14]
inhibitory_neuron <- unique(scdata$Cluster)[15:25]
purk_cell <- unique(scdata$Cluster)[26:27]

#create new column in sc data which includes the now fewer different types of cells
scdata$CellType[scdata$Cluster %in% excitatory_neurons] = "Excitatory Neuron"
scdata$CellType[scdata$Cluster %in% inhibitory_neuron] = "Inhibitory Neuron"
scdata$CellType[scdata$Cluster %in% purk_cell] = "Purkinje Cell"
scdata$CellType[scdata$Cluster == "Gran"] = "Granule cell"
scdata$CellType[scdata$Cluster == "End"] = "Endothelial cell"
scdata$CellType[scdata$Cluster == "Per"] = "Pericyte"
scdata$CellType[scdata$Cluster == "Ast"] = "Astrocyte"
scdata$CellType[scdata$Cluster == "Ast_Cer"] = "Cerebellar Astrocyte"
scdata$CellType[scdata$Cluster == "Oli"] = "Oligodendrocyte"
scdata$CellType[scdata$Cluster == "OPC"] = "Oligodendrocyte Precursor Cell"
scdata$CellType[scdata$Cluster == "OPC_Cer"] = "Cerebelllar Oligodendrocyte Precursor Cell"
scdata$CellType[scdata$Cluster == "Mic"] = "Microglia"


#prepare dataframe holding fishers results of miRNAs and cell types
AllmiRTarBase_fishersresults_all <- data.frame(miRNA=as.character(),
                                               CellType=as.character(),
                                               pvalue=as.numeric(), OddsRatio = as.numeric(),
                                               Total_Targets=as.numeric(), Total_Markers=as.numeric(),
                                               Number_of_Markers_Targeted=as.numeric(), Markers_Targeted=as.character())


cell_types = unique(scdata$CellType)
scdata_filtered <- subset(scdata, logFC>0.5)

for (i in All_miRNAs) {
  print(i)
  AllmiRTarBase_miR_targets = unique(subset(AllmiRTarBase, miRNA == i)[,3])
  AllmiRTarBase_miR_targets = as.character(AllmiRTarBase_miR_targets)
  for (j in cell_types) {
    cell_type_markers = as.character(unique(subset(scdata_filtered, scdata_filtered$CellType == j)[,1]))
    
    AllmiRTarBase_targets_and_markers = subset(cell_type_markers, cell_type_markers %in% AllmiRTarBase_miR_targets)
    AllmiRTarBase_targets_and_markers = unique(AllmiRTarBase_targets_and_markers)
    list_AllmiRTarBase_targets_and_markers = paste(AllmiRTarBase_targets_and_markers, collapse = ", ")
    
    A = length(AllmiRTarBase_targets_and_markers)
    B = length(AllmiRTarBase_miR_targets) - A 
    C = length(cell_type_markers) - A
    D = length(allhumangenes$Gene.name)- A - B - C
    
    fishersmatrix = matrix(c(A,B,C,D), nrow=2, ncol=2, dimnames = list(c('Cell marker', 'Not a cell marker'), c('Target gene', 'Not a target gene')))
    
    AllmiRTarBase_fishersresult = fisher.test(fishersmatrix)
    
    my_output_AllmiRTarBase <- data.frame(miRNA = i, CellType = j,
                                          pvalue=AllmiRTarBase_fishersresult$p.value[1], OddsRatio=AllmiRTarBase_fishersresult$estimate[1],
                                          Total_Targets=length(AllmiRTarBase_miR_targets), Total_Markers=length(cell_type_markers),
                                          Number_of_Markers_Targeted=A,
                                          Markers_Targeted = list_AllmiRTarBase_targets_and_markers)
    
    AllmiRTarBase_fishersresults_all = rbind(AllmiRTarBase_fishersresults_all, my_output_AllmiRTarBase)
  }
}



AllmiRTarBase_fishersresults_all$FDR <- p.adjust(AllmiRTarBase_fishersresults_all$pvalue, method = "fdr")
AllmiRTarBase_fishersresults_all$FDR.bonf <- p.adjust(AllmiRTarBase_fishersresults_all$pvalue, method = "bonferroni")
AllmiRTarBase_fishersresults_all_padj <- subset(AllmiRTarBase_fishersresults_all, AllmiRTarBase_fishersresults_all$FDR <0.5)


write.table(AllmiRTarBase_fishersresults_all, file="All miRNAs MirTarBase Fishers Enrichment Results logFC CellMarkers_0.5.txt",
            sep="\t",col.names = T,quote=F,row.names = F)



#scatterplot for all cell types

ggplot(AllmiRTarBase_fishersresults_all, aes(y=OddsRatio, x=FDR, colour=CellType))+geom_point()+  
  labs(x="FDR", y= "Odds Ratio",
       title="Cell type enrichment for all miRNA targets identified in miRTarBase", color = "Cell type")+
  theme_bw() +
  geom_vline(xintercept=0.05, linetype = 2)



#scatterplot for ASD/non-ASD miRNAs

AllmiRTarBase_fishersresults_all$ASDmiRNA <- "Non-ASD"
AllmiRTarBase_fishersresults_all$ASDmiRNA[AllmiRTarBase_fishersresults_all$miRNA %in% dysregulatedmiRNAs] <- "ASD"

ggplot(AllmiRTarBase_fishersresults_all, %>% arrange(ASDmiRNA) aes(y=OddsRatio, x=FDR, colour=ASDmiRNA))+geom_point()+ 
  labs(x="FDR", y= "Odds Ratio",
       title="Cell type enrichment for all miRNA targets identified in miRTarBase", color = "miRNA")+
  scale_color_manual(breaks = c("ASD", "Non-ASD"), values=c("red", "gray54"))

ggplot(AllmiRTarBase_fishersresults_all, aes(y=OddsRatio, x=FDR, colour=ASDmiRNA))+geom_point()+ 
  labs(x="FDR", y= "Odds Ratio",
       title="Cell type enrichment for all miRNA targets identified in miRTarBase", color = "miRNA")+
  scale_color_manual(breaks = c("ASD", "Non-ASD"), values=c("red", "gray54")) +
  theme_bw() +
  geom_vline(xintercept=0.05, linetype = 2)



##scatter showing cerebellum/cerebrum
AllmiRTarBase_fishersresults_all$brainregion <- "Cerebrum & Cerebellum"

AllmiRTarBase_fishersresults_all$brainregion[AllmiRTarBase_fishersresults_all$CellType == 'Excitatory Neuron' | AllmiRTarBase_fishersresults_all$CellType == 'Inhibitory Neuron'] <- "Cerebrum"

AllmiRTarBase_fishersresults_all$brainregion[AllmiRTarBase_fishersresults_all$CellType == 'Purkinje Cell' | AllmiRTarBase_fishersresults_all$CellType == 'Granule cell' | AllmiRTarBase_fishersresults_all$CellType == 'Cerebellar Astrocyte' |   AllmiRTarBase_fishersresults_all$CellType == 'Cerebelllar Oligodendrocyte Precursor Cell'] <- "Cerebellum"


ggplot(AllmiRTarBase_fishersresults_all, aes(y=OddsRatio, x=FDR, colour=brainregion))+geom_point()+ 
  labs(x="FDR", y= "Odds Ratio",
       title="Cell type enrichment for all miRNA targets identified in miRTarBase, by brain region", color = "Brain region")+
  theme_bw() +
  geom_vline(xintercept=0.05, linetype = 2)



#scatterplot for neurons/glia

AllmiRTarBase_fishersresults_all$NeuronGlia <- ifelse(AllmiRTarBase_fishersresults_all$CellType == 'Excitatory Neuron' | AllmiRTarBase_fishersresults_all$CellType == 'Inhibitory Neuron' | AllmiRTarBase_fishersresults_all$CellType == 'Purkinje Cell' | AllmiRTarBase_fishersresults_all$CellType == 'Granule cell', 'Neurons',
                                                      ifelse( AllmiRTarBase_fishersresults_all$CellType == 'Cerebellar Astrocyte' |   AllmiRTarBase_fishersresults_all$CellType == 'Cerebelllar Oligodendrocyte Precursor Cell' | AllmiRTarBase_fishersresults_all$CellType == 'Endothelial cell' | AllmiRTarBase_fishersresults_all$CellType == 'Pericyte' | AllmiRTarBase_fishersresults_all$CellType == 'Astrocyte' | AllmiRTarBase_fishersresults_all$CellType == 'Oligodendrocyte' | AllmiRTarBase_fishersresults_all$CellType == 'Oligodendrocyte Precursor Cell' | AllmiRTarBase_fishersresults_all$CellType == 'Microglia', 'Glia',
                                                              ifelse(AllmiRTarBase_fishersresults_all$CellType == NA, 'Other')))


ggplot(AllmiRTarBase_fishersresults_all, aes(y=OddsRatio, x=FDR, colour=NeuronGlia))+geom_point()+
  labs(x="FDR", y= "Odds Ratio",
       title="Cell type enrichment for all miRNA targets identified in miRTarBase", color = "Cell type")+ 
  theme_bw() +
  geom_vline(xintercept=0.05, linetype = 2)

