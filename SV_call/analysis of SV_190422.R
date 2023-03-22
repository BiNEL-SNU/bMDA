library(RColorBrewer)
library(stringr)
library(ggbeeswarm)
library(reshape2)
library(ggsignif)
library(ComplexHeatmap)
library(gplots)
library(ggplot2)
library(Rcpp)
library(ggrepel)
library(ggrastr)
library(MASS)
library(tidyverse)
library(broom)
library(tibble)
library(tidyr)
library(viridis)
library(dplyr)

display.brewer.all(colorblindFriendly = TRUE)


setwd("D:/BiNEL/01_Research/02_Single cell genomics/01_daily/21.10.14 NGS73")
groupID = "190422"; slideNo="31"
#groupID = "190513"; slideNo="34"


infoAll = read.table(paste("../20.12.29 NGS70/cluster_", slideNo, "_FinalCluster.csv", sep=""), header = T, stringsAsFactors = F, sep ="," ) %>%
        left_join(read.table(paste("../20.12.29 NGS70/CNA row order slide", slideNo, ".csv", sep=""), sep = ",", header=T)) %>%
        mutate(name = sapply(str_split(name, "\\.", n=2), function(x) { return(x[1]);})) %>%
        select(-clusterSNV) %>%
        select(-clusterPhyloAlpha, -clusterPhyloCNAorder, -clusterPhyloSNVorder) %>%
        left_join(read.table("../20.12.29 NGS70/slide31 phylogeny.csv", header = T, stringsAsFactors = F, sep =",") %>%
                          select(name, clusterPhyloAlpha,clusterPhyloCNA, clusterPhyloSNV,clusterPhyloCNAorder,clusterPhyloSNVorder,clusterPhyloSNV_WESorder) %>% rename(name2=name)) %>%
        mutate(barcodeNo = ifelse(name == "190422_Tumor", "Bulk", sapply(str_split(name, "_", n=4), function(x) { return(x[3]);})))


#data=read.table('confident.site.01.variant',header = T)
INPUT1 = "190422_SV_MANTAandGRIDSS_merged.filtered" #high confident variant calls across all samples
INPUT2 = "190422_SV_all_merged_lowConf_bulkFilter.filtered" #variants detected on bulk
anno1 = read.table(paste(groupID, " SV/", INPUT1, '.annotated.tsv', sep=""),header = T, sep = "\t") %>%
        select(ID, Annotation_mode, CytoBand, Gene_name, Gene_count, AnnotSV_ranking_score, AnnotSV_ranking_criteria, ACMG_class)
data1= read.table(paste(groupID, " SV/", INPUT1, '.variant', sep=""),header = T, sep = ",") %>%
        left_join(anno1) #%>% left_join(anno2)

anno2 = read.table(paste(groupID, " SV/", INPUT2, '.annotated.tsv', sep=""),header = T, sep = "\t") %>%
        select(ID, Annotation_mode, CytoBand, Gene_name, Gene_count, AnnotSV_ranking_score, AnnotSV_ranking_criteria, ACMG_class)
data2= read.table(paste(groupID, " SV/", INPUT2, '.variant', sep=""),header = T, sep = ",") %>%
        left_join(anno2) #%>% left_join(anno2)

#annovar
#bedpe = read.table(paste(groupID, " SV/", INPUT2, '.bedpe', sep=""),header = F, sep = "\t")
#colnames(bedpe) = c("CHR1", "START1", "END1", "CHR2", "START2", "END2", "ID")
#anno2 = read.table(paste(groupID, " SV/", INPUT2, '.hg19_multianno.txt', sep=""),header = T, sep = "\t") %>%
#        mutate(ID = bedpe$ID) %>% select (ID, Func.refGene, Gene.refGene, GeneDetail.refGene)

data_merged = data1 %>% filter(X190422_Tumor != 1) %>%
  bind_rows(data2) %>%
  mutate(idx = 1:n()) %>% group_by(idx) %>%
  mutate(Gene_name_trunc = ifelse(length(unlist(strsplit(Gene_name, ";"))) > 3, paste(c(unlist(strsplit(Gene_name, ";"))[1:3], "..."), collapse = ";"), Gene_name)) %>% 
  ungroup() %>% select(-idx) %>%
  mutate(Gene_count_trunc = ifelse(Gene_count >= 10, 10, Gene_count),
         significance = ifelse(ACMG_class %in% c(5) | AnnotSV_ranking_score >= 0.99, 1,
                               ifelse(ACMG_class %in% c(4) | AnnotSV_ranking_score >= 0.90, 0.5, 0)))


dataFiltered = data_merged %>%
  filter(FILTER == "PASS", SUPP >= 2,
         CHROM %in% c(as.character(1:24), "X", "Y"), 
         CHR2 %in% c(as.character(1:24), "X", "Y")) %>%
  filter(CHROM %in% c(5)) %>%
  #filter(CHROM %in% c(5,12), CHR2 %in% c(5,12)) %>%
  filter(SVTYPE == "TRA" | (SVTYPE %in% c("INV", "DUP", "DEL") & abs(SVLEN) > 1000))

if(F) {# find variants on each subclone
target_sample = (inputDF %>% filter(cluster_custom %in% c(1)) %>% mutate(name2 = paste("X", name,sep="")))$name2 %>% str_replace_all("201229-bMDA-BRCA", "201229.bMDA.BRCA")
off_target_sample = (inputDF %>% filter(cluster_custom %in% c(2,3, 4, 5)) %>% mutate(name2 = paste("X", name,sep="")))$name2 %>% str_replace_all("201229-bMDA-BRCA", "201229.bMDA.BRCA")

dataFiltered = dataFiltered %>%
  filter((rowSums(dataFiltered[,target_sample]) >=2) & (rowSums(dataFiltered[,off_target_sample]) <= 1))
}

dataFiltered %>% group_by(SVTYPE) %>% summarise(count = n())
dataFiltered$SUPP %>% hist(breaks = 20)
dataFiltered %>% select(CHROM,POS, SUPP)
dataFiltered$SVLEN %>% abs %>% hist(xlim = c(0,10000), breaks = 100000)
dataFiltered %>% arrange(-abs(SVLEN)) %>% select(CHROM,POS, CHR2,END, SVTYPE, SVLEN)


# export gene gene name for TRA
#data_annotated %>% select(CHROM, POS,SVLEN, SVTYPE, Gene_name, Gene_count, AnnotSV_ranking_score, Func.refGene, Gene.refGene) %>% filter(SVTYPE == "TRA") %>% View
dataFiltered %>% 
        filter(SVTYPE == "TRA") %>%
        filter(Gene_name != "") %>%
        mutate(Gene_name = ifelse(Gene_name == "LOC644135;PDE7B", "PDE7B", Gene_name)) %>%
        group_by(Gene_name) %>% summarise(CHROM = CHROM[1], POS = round(mean(POS)), Gene_name = Gene_name[1]) %>%
        arrange(as.numeric(CHROM), POS) %>%
        mutate(CHROM = paste("hs", CHROM, sep=""), POS2 = POS) %>%
        select(CHROM, POS, POS2, Gene_name) %>%
        write.table("190422_TRA_GeneName.txt", sep = "\t", row.names = F, col.names = F, quote = F)

data_for_mat=dataFiltered %>% select(starts_with("X"))#[,14:ncol(dataFiltered)] #dev
sampleNames = colnames(data_for_mat) %>% str_replace_all("X", "") %>% str_replace_all("201229.bMDA.BRCA","201229-bMDA-BRCA")
nSamples = length(sampleNames)
inputDF = as.data.frame(sampleNames, stringsAsFactors=FALSE) %>% rename(name = sampleNames) %>%
        left_join(infoAll)


mat = t(as.matrix(data_for_mat))
rownames(mat) = sampleNames
colnames(mat) = paste0(dataFiltered$CHROM, ":", dataFiltered$POS)



##### phylogeny based drawing
inputDF = inputDF %>% 
        mutate(cluster_custom = ifelse(clusterPhyloSNV_WESorder %in% 1:4, 1, 
                                       ifelse(clusterPhyloSNV_WESorder %in% 5:8, 2,
                                              ifelse(clusterPhyloSNV_WESorder %in% 9:13, 3,
                                                     ifelse(clusterPhyloSNV_WESorder %in% 14:18, 4,
                                                            ifelse(clusterPhyloSNV_WESorder %in% 19:24, 5, 0))))))


ha_sample <- rowAnnotation(phy = inputDF$clusterPhyloSNV_WESorder, #clusterPhyloCNAorder
                           col = list(phy = setNames(c("black", colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2")[1:3])(24)), 0:24)))
#c("#4888C8","#5054A0","#C94D98", "#F8DA30","#79B752","#CE3037")

ha_sample <- rowAnnotation(phy = inputDF$cluster_custom,
                           col = list(phy = setNames(c("black", colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2")[1:3])(5)), 0:5)))

ha_method = rowAnnotation(#c = variantClust,
                          m = inputDF$groupMethod, 
                          markingNo = anno_text(inputDF$markingNo %>% as.character, gp = gpar(fontsize = 9)),
                          #AUC = anno_text(sprintf("%.3f", inputDF$AUC), gp = gpar(fontsize = 9)),
                          col = list(#c = setNames(colCluster[1:(variantClust %>% unique %>% length)], variantClust %>% unique %>% sort),
                                     m = setNames(c("#CE3037", RColorBrewer::brewer.pal(1, "Dark2"))[1:(inputDF$groupMethod %>% unique %>% length)], 
                                                  inputDF$groupMethod %>% unique %>% sort %>% rev))) # just = "center"
ha_chrompos <- HeatmapAnnotation(
        SV_size = anno_barplot(dataFiltered$SVLEN %>% abs),
        #pos = anno_text(colnames(mat),location = 1, rot = 90, just = "right", gp = gpar(fontsize = 9)),
        geneCount = dataFiltered$Gene_count_trunc,
        col = list(geneCount = colorRamp2(c(0, 10), c("#FFFFFF", RColorBrewer::brewer.pal(8, "Set1")[3])))
        )


## cancer consensus genes
consensusFileName = "../22.03.05 NGS77/variant heatmap/annotation/Census_allMon Mar 28 15_16_57 2022.csv"
consensusDF = read.table(consensusFileName, header = T, stringsAsFactors = F, sep ="," ) %>% rename(Gene.refGene = Gene.Symbol) %>%
        select(Gene.refGene, Tier) %>%filter(Tier %in% c(1,2))

#all anno
df_geneName = dataFiltered %>% select(Gene_name_trunc) %>% mutate(index = 1:n()) %>% filter(Gene_name_trunc != "")

ha_gene <- HeatmapAnnotation(
        significance = dataFiltered$significance,
        gene = anno_mark(at = df_geneName$index, labels = df_geneName$Gene_name_trunc, side="bottom"))

#concensus anno
df_geneName = dataFiltered %>% select(SVTYPE, Gene_name) %>% mutate(index = 1:n()) %>% filter(Gene_name != "") %>%
        group_by(index) %>%
        filter(any(sapply(consensusDF$Gene.refGene, function(x) {return(x %in% unlist(strsplit(Gene_name, ";")))}))) %>%
        mutate(Gene_name_selected = paste(consensusDF$Gene.refGene[sapply(consensusDF$Gene.refGene, function(x) {return(x %in% unlist(strsplit(Gene_name, ";")))})], collapse= ";")) %>%
        mutate(Gene_name_selected = ifelse(length(unlist(strsplit(Gene_name_selected, ";"))) > 3, paste(c(unlist(strsplit(Gene_name_selected, ";"))[1:3], "..."), collapse = ";"), Gene_name_selected)) %>%
        group_by(SVTYPE, Gene_name) %>% filter(row_number()==1)
ha_gene <- HeatmapAnnotation(
        significance = dataFiltered$significance,
        idx = anno_text((dataFiltered %>% select(SVTYPE, Gene_name) %>% mutate(index = 1:n()))$index),
        gene = anno_mark(at = df_geneName$index, labels = df_geneName$Gene_name_selected, side="bottom"),
        col = list(significance = colorRamp2(c(0, 1), c("#FFFFFF", RColorBrewer::brewer.pal(8, "Set1")[5]))))


varCols = c("#FFFFFF", RColorBrewer::brewer.pal(9, "Blues")[6]) #c("#C7C7C7", "#1F77B4")
colors = structure(varCols, names = c("0", "1")) # black, red, green, 
library(circlize); col_fun = colorRamp2(c(0, 1), varCols)

#pdf("heatmap.pdf",width=30/2.54, height=nrow(mat)*0.55/2.54) #whole chr
pdf("heatmap.pdf",width=17/2.54, height=nrow(mat)*0.47/2.54) #chr12
pdf("heatmap.pdf",width=18/2.54, height=nrow(mat)*0.39/2.54) #chr5

Heatmap(mat,
        name = "Variant Detected", #title of legend
        col = col_fun,
        #col = colors, 
        na_col = "#C7C7C7",
        row_names_gp = gpar(fontsize = 7),
        show_row_names =F,
        cluster_columns = T,
        row_split = inputDF$groupMethod %>% factor(c("Tumor", "bMDA")),
        column_split = dataFiltered$SVTYPE,
        column_gap = unit(2, "mm"),
        use_raster = F,
        top_annotation = ha_chrompos,
        left_annotation = ha_sample,
        right_annotation = ha_method,
        bottom_annotation = ha_gene,
        #col = col_fun,
        cluster_rows = F,#dend, #T,
        #row_order = inputDF$rowRank %>% order,
        #row_order = with(inputDF, order(groupTemplate,-AUC)),
        row_order = with(inputDF, order(clusterPhyloSNV_WESorder)), #clusterPhyloSNVorder #clusterPhyloCNAorder
        #column_order = order(as.numeric(gsub("bin", "", colnames(cnv.matrix)))),
        show_column_dend = FALSE, show_column_names =F
)
dev.off()


## plot chr12 graphical view
library(chromoMap)
cytoBand = read.table("190422 SV/cytoBand.hg19.txt", header = F, stringsAsFactors = F, sep ="\t" )
colnames(cytoBand) = c("CHROM", "START", "END", "ARM", "BAND")
cytoBand = cytoBand %>% 
  mutate(chr_num = str_replace(CHROM, "chr", ""), 
         chr_num = str_replace(chr_num, "X", "23"),
         chr_num = str_replace(chr_num, "Y", "24"),
         chr_num = as.integer(chr_num)) %>%
  arrange(chr_num, START)

chromFile = cytoBand %>% group_by(CHROM) %>% summarise(START = 1, END = max(END))
centromere = cytoBand %>% filter(BAND %in% c("acen")) %>% group_by(CHROM, BAND) %>% 
  summarise(START = min(START), END = max(END)) %>% mutate(MID = as.integer((START+END)/2))
chromFile = chromFile %>% left_join(centromere %>% select(CHROM, MID) %>% rename(CENTROMERE = MID))

cyto_anno = cytoBand %>% mutate(ID = paste("id",as.character(1:n()), sep = "")) %>%
  filter(!(BAND %in% c("gneg", 'acen', 'gvar', 'stalk'))) %>%
  select(ID, CHROM, START, END, BAND, ARM)

sample_w_kata = (inputDF %>% filter(cluster_custom == 5) %>% mutate(name2 = paste("X", name,sep="")))$name2 %>% str_replace_all("201229-bMDA-BRCA", "201229.bMDA.BRCA")
sample_wo_kata = (inputDF %>% filter(!(cluster_custom %in% c(0,5))) %>% mutate(name2 = paste("X", name,sep="")))$name2 %>% str_replace_all("201229-bMDA-BRCA", "201229.bMDA.BRCA")


data_chromPlot = data_merged %>%
  filter(FILTER == "PASS", SUPP >= 2,
         CHROM %in% c(as.character(1:24), "X", "Y"), 
         CHR2 %in% c(as.character(1:24), "X", "Y")) %>%
  arrange(-abs(SVLEN)) %>% 
  filter(SVTYPE %in% c("INV", "TRA", "DUP", "DEL")) %>% #DUP, DEL
  filter(CHROM %in% c(5,12), CHR2 %in% c(5,12)) %>%
  filter(SVTYPE == "TRA" | (SVTYPE %in% c("INV", "DUP", "DEL") & abs(SVLEN) > 100000))

  
data_chromPlot_target = data_chromPlot %>% 
  select(CHROM,POS, CHR2,END, SVTYPE, SVLEN)  %>% 
  filter(rowSums(data_chromPlot[,sample_w_kata]) >= 1) %>%
  #filter(rowSums(data_chromPlot[,sample_wo_kata]) >= 1) %>%
  mutate(name1 = paste(CHROM, ":", POS, sep=""), name2 = paste(CHR2, ":", END, sep=""), ploidy1 = 1, ploidy2 = 1)

link_data = data_chromPlot_target %>% select(name1, ploidy1, name2, ploidy2, SVTYPE) %>%
  rename(V1 = name1, V2 = ploidy1, V3 = name2, V4 = ploidy2, V5 = SVTYPE)


for_cyto_anno = 
  bind_rows(data_chromPlot_target %>% select(CHROM, POS, name1) %>% rename(ID = name1),
            data_chromPlot_target %>% select(CHR2, END, name2) %>% rename(CHROM = CHR2, POS = END, ID = name2)) %>%
  unique %>% rename(START = POS) %>%
  mutate(ARM = "NotDefined", BAND = "SV", END = START + 1, CHROM = paste("chr", CHROM, sep = ""))

chromFile %>% filter(CHROM %in% c("chr12", "chr5")) %>% #select(-CENTROMERE) %>%
  write.table("chromosome_file.txt", sep = "\t", quote = F, row.names = F, col.names = F)
cyto_anno %>% filter(CHROM %in% c("chr12", "chr5")) %>% 
  write.table("annotation_file.txt", sep = "\t", quote = F, row.names = F, col.names = F)

cyto_anno %>% filter(CHROM %in% c("chr12", "chr5")) %>% 
  bind_rows(for_cyto_anno) %>% arrange(CHROM, START, END) %>%
  write.table("annotation_file_withLink.txt", sep = "\t", quote = F, row.names = F, col.names = F)

varCols = c("#FFFFFF", RColorBrewer::brewer.pal(8, "Dark2")[3])
varCols = c("#FFFFFF", RColorBrewer::brewer.pal(8, "Greys")[6])

colors = structure(varCols, names = c("0", "1")) # black, red, green, 
library(circlize); col_fun = colorRamp2(c(0, 1), varCols)

chromoMap(ch.files = c("chromosome_file.txt"),
          data.files = c("annotation_file.txt"),
          segment_annotation = T,
          n_win.factor = 1,
          data_based_color_map = T,
          win.summary.display=T,
          data_type = "categorical",
          #data_colors = list(c("orange","yellow","orange","yellow")),
          data_colors = list(c(col_fun(0.25), col_fun(0.75), col_fun(1), col_fun(0.5))),
          ch_gap = 50,
          top_margin = 100,
          legend= T, 
          lg_x = 100,
          export.options = T,
)



link_cols = RColorBrewer::brewer.pal(12, "Paired")
chromoMap(ch.files = c("chromosome_file.txt"),
          data.files = c("annotation_file_withLink.txt"),
          segment_annotation = T,
          n_win.factor = 1,
          data_based_color_map = T,
          data_type = "categorical",
          show.links = T,
          loci_links = link_data,
          links.colors = c(link_cols[4], link_cols[2], link_cols[6], link_cols[4]),
          data_colors = list(c(col_fun(0), col_fun(0.25), col_fun(0.75), col_fun(1), col_fun(0.5))),
          ch_gap = 50,
          top_margin = 100,
          legend= T, 
          lg_x = 100,
          links.lg_x = 200,
          export.options = T,
)

          