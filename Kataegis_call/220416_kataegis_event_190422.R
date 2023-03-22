rm(list= ls())
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(KataegisPortal)
library(stringr)
library(ggplot2)
library(ggrastr)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyr)
library(dplyr)

root_dir = "D:/BiNEL/01_Research/02_Single cell genomics/01_daily/21.10.14 NGS73/ketaegis"
data_path = "bMDA_WGS_somatic_variants_190422_6_20_30"
setwd(root_dir)

source("readVCF.r")

dir_list = list.files(data_path)
sample_name = c()
sample_name[1] = 'Bulk_varscan'
sample_name[2] = 'Bulk_GATK'
for (i in 3:length(dir_list)){
  sample_name[i] = str_split(dir_list[i], pattern = "_")[[1]][3]
}

ff <- function(input, TargetFolder){
  CountsFile = paste(TargetFolder, input, sep="/")
  data = as.matrix(readVCF(file = CountsFile))
}

print("Read File")
vcfs = lapply(dir_list, ff, TargetFolder = data_path)
print("Done")

temp <- vcfs[[1]]
temp <- temp[1:max(which(temp[,'#CHROM'] == 'X')),]
temp[,1] <- str_replace_all(temp[,1], 'chr', '')
temp <- as.data.frame(temp)[1:4]
temp[,2] <- as.integer(temp[,2])
colnames(temp) <- c('chr', 'pos', 'ref', 'alt')
temp <- temp %>% mutate(chr=paste('chr', chr, sep=""))
temp <- temp %>% filter(alt %in% c("A", "G", "T", "C"))
inputs <- list()
inputs[[1]] = temp

for (i in 2:length(vcfs)){
  temp <- vcfs[[i]]
  temp <- temp[1:max(which(temp[,'#CHROM'] == 'X')),]
  temp[,1] <- str_replace_all(temp[,1], 'chr', '')
  temp <- as.data.frame(temp)[1:4]
  temp[,2] <- as.integer(temp[,2])
  colnames(temp) <- c('chr', 'pos', 'ref', 'alt')
  temp <- temp %>% mutate(chr=paste('chr', chr, sep=""))
  temp <- temp %>% filter(alt %in% c("A", "G", "T", "C"))
  inputs[[i]] <- temp
}

katas = lapply(inputs, mutSNP.input, chr = "chr", pos = "pos", ref = "ref", alt = "alt", build = "hg19")
katpoints = list()
for (i in 1:length(sample_name)){
  temp <- katPoint(katas[[i]], sample = sample_name[i], txdb = TxDb.Hsapiens.UCSC.hg19.knownGene)
  katpoints[[i]] <- temp
}

results <- bind_rows(katpoints)
results <- results %>% mutate(chrom = str_replace(chrom, 'chr', ''))
results <- results %>% mutate(chrom = ifelse(chrom=='X', 23, chrom))
results <- results %>% mutate(chrom = as.integer(chrom), start = as.integer(start), end = as.integer(end))
results <- arrange(results, chrom, start, end)
bed_results <- results %>% dplyr::select(chrom, start, end)
write.table(bed_results, '220416_kataegis_result.bed', sep = '\t', row.names = FALSE, col.names = FALSE)

merge_result <- read.table("220416_kategis_merge_result.bed")
kata_event_counts <- list()
for(i in 1:nrow(merge_result)){
  merge_start <- merge_result[i,2]
  merge_end <- merge_result[i,3]
  merge_chrom <- merge_result[i,1]
  temp <- results %>% filter(chrom == merge_chrom, start >= merge_start, end <= merge_end) %>% 
    dplyr::select(-start, -end) %>% mutate(kataegis_event = 1,start=merge_start, end=merge_end)
  kata_event_counts[[i]] <- temp
}
kata_event_count <- bind_rows(kata_event_counts) %>%
  arrange(chrom, start, end) %>% dplyr::select(sample, chrom, start, end, everything() )

bed_anno = read.table("../kataegis_220416.annotated.bed", header=F, stringsAsFactors = F, sep="\t")
colnames(bed_anno) <- c("chrom", "start", "end", "gene")
kata_event_count = kata_event_count %>% left_join(bed_anno)
kata_event_count %>% dplyr::select(geneName, gene) %>% unique
kata_event_count %>% dplyr::select(chrom, start, end, geneName) %>% unique %>% 
  write.table('220416_kategis_merge_result.annotated_new.bed', sep = '\t', row.names = FALSE, col.names = FALSE, quote = F)

target_kataegis = kata_event_count %>% group_by(chrom, start, end, geneName, gene) %>%
  summarise(n_sample = n(), conf0_cnt = sum(confidence == 0), conf1_cnt = sum(confidence == 1), conf2_cnt = sum(confidence == 2), conf3_cnt = sum(confidence == 3)) %>%
  filter(n_sample > 3, (n_sample - conf0_cnt) >= 3) %>% View

kata_event_count %>% filter(chrom == 12, start == 1460520)
kata_event_count %>% filter(chrom == 19, start == 11575721)


kata_event_count %>% filter(chrom == 1, start == 17231327)
kata_event_count %>% filter(chrom == 1, start == 30196204)
kata_event_count %>% filter(chrom == 2, start == 102618954)
kata_event_count %>% filter(chrom == 2, start == 113103985)
kata_event_count %>% filter(chrom == 3, start == 68303283)
kata_event_count %>% filter(chrom == 5, start == 15248631)
kata_event_count %>% filter(chrom == 5, start == 43464626)
kata_event_count %>% filter(chrom == 8, start == 42644818)
kata_event_count %>% filter(chrom == 12, start == 2828979)
kata_event_count %>% filter(chrom == 21, start == 44755928)
kata_event_count %>% filter(chrom == 23, start == 104252671)





target_kataegis = kata_event_count %>% group_by(chrom, start, end) %>%
  summarise(n_sample = n(), conf0_cnt = sum(confidence == 0), conf1_cnt = sum(confidence == 1), conf2_cnt = sum(confidence == 2), conf3_cnt = sum(confidence == 3)) %>%
  #filter(n_sample > 1, (n_sample - conf0_cnt) >= 1) %>% 
  filter(n_sample > 3, (n_sample - conf0_cnt) >= 3) 
kata_event_selected = target_kataegis %>%
  dplyr::select(chrom,start,end) %>% mutate(sel = TRUE)

final_data <- kata_event_count %>% 
  left_join(kata_event_selected) %>% filter(sel == TRUE) %>% dplyr::select(-sel) %>% ##filtering
  dplyr::select(sample, chrom, start, end, kataegis_event) %>%
  spread(key = sample, value = kataegis_event)
final_data[is.na(final_data)] = 0
write.csv(final_data, "../220416_kataegis_event_persample.csv", row.names=FALSE)

pdf(paste("all chromosome.pdf", sep = ""), width=18/2.54, height=10/2.54)
for (i in c(1:length(sample_name))) {
  mutDis.plot(plot.data = katas[[i]], sample=sample_name[i])
}
dev.off()


dfs = list()
for (n in 1:nrow(target_kataegis)) {
  thisChrom = target_kataegis$chrom[n]
  thisData = kata_event_count %>% filter(chrom == thisChrom, start == target_kataegis$start[n])
  thisSamples = thisData$sample
  idx_kata = match(thisSamples, sample_name)
  idx_notkata = c(1:length(sample_name))[!(c(1:length(sample_name)) %in% idx_kata)]
  thisAnno = paste(thisChrom, "_", target_kataegis$start[n], "_", thisData$geneName[1], sep = "")
  df_kata = data.frame(sample_name[idx_kata], rep(RColorBrewer::brewer.pal(9, "Blues")[6], length(idx_kata)), rep(thisAnno, length(idx_kata))) 
  colnames(df_kata) = c("sampleIdx", "color", "geneName")
  df_notkata = data.frame(sample_name[idx_notkata], rep("#C7C7C7", length(idx_notkata)), rep(thisAnno, length(idx_notkata)))
  colnames(df_notkata) = c("sampleIdx", "color", "geneName")
  dfs[[n]] = bind_rows(df_kata, df_notkata)
  if (T) {
    pdf(paste(thisAnno, "_kata", ".pdf", sep = ""), width=18/2.54, height=10/2.54)
    for (i in idx_kata){
      mutDis.plot(plot.data = katas[[i]], sample=sample_name[i], chr = thisChrom)
      #(plot.data = katas[[i]], sample=sample_name[i])
    }
    dev.off()
    pdf(paste(thisAnno, "_notkata", ".pdf", sep = ""), width=18/2.54, height=10/2.54)
    for (i in idx_notkata){
      mutDis.plot(plot.data = katas[[i]], sample=sample_name[i], chr = thisChrom)
      #(plot.data = katas[[i]], sample=sample_name[i])
    }
    dev.off()
  }
  
}



## spatial marking
groupID = "190422"; slideNo="31"
infoAll = read.table(paste("../../20.12.29 NGS70/cluster_", slideNo, "_FinalCluster.csv", sep=""), header = T, stringsAsFactors = F, sep ="," ) %>%
  left_join(read.table(paste("../../20.12.29 NGS70/CNA row order slide", slideNo, ".csv", sep=""), sep = ",", header=T)) %>%
  mutate(name = sapply(str_split(name, "\\.", n=2), function(x) { return(x[1]);})) %>%
  dplyr::select(-clusterSNV) %>%
  dplyr::select(-clusterPhyloAlpha, -clusterPhyloCNAorder, -clusterPhyloSNVorder) %>%
  left_join(read.table("../../20.12.29 NGS70/slide31 phylogeny.csv", header = T, stringsAsFactors = F, sep =",") %>%
              dplyr::select(name, clusterPhyloAlpha,clusterPhyloCNA, clusterPhyloSNV,clusterPhyloCNAorder,clusterPhyloSNVorder,clusterPhyloSNV_WESorder) %>% dplyr::rename(name2=name)) %>%
  mutate(barcodeNo = ifelse(name == "190422_Tumor", "Bulk", sapply(str_split(name, "_", n=4), function(x) { return(x[3]);}))) %>%
  mutate(sampleIdx = sapply(str_split(name, "_"), function(x) { return(x[3]);})) %>% 
  dplyr::select(name, sampleIdx, isolationRemark) %>%
  dplyr::rename(NOTE = isolationRemark) %>%
  left_join(read.table("../../20.12.29 NGS70/region marking on H&E/31 - 190422_outputFinal.txt", stringsAsFactors = F, header=T))


df_for_spatial = infoAll %>% left_join(bind_rows(dfs)) %>%
  arrange(geneName) %>%
  filter(!is.na(color))

df_for_spatial %>% write.table("31 - 190422_kataegis.csv", quote = F, row.names = F, sep=",")





### CNA plots -----
bin.count = "50k"; undo.sd = 0.8
name.file.postfix = paste("uniqMap", bin.count, "k150.data.txt", sep=".")
baseDir = "CNA plot"

df_for_cna = infoAll %>% right_join(bind_rows(dfs)) %>%
  arrange(geneName) %>% filter(sampleIdx != "Bulk_varscan") %>% mutate(name = ifelse(is.na(name), "190422_Tumor.3G", name))

#for(thisGene in unique(df_for_cna$geneName)) {
  
thisGene = "12_1460520_LINC00942"

## MergeLevel Preprocessing - perform everytime if input list changes (e.g. removal of bad quality samples)
inputDF = df_for_cna %>% filter(geneName == thisGene) %>% #kata+
  filter(color == "#4292C6")

inputDF = df_for_cna %>% filter(geneName == thisGene) %>% #kata-
  filter(color != "#4292C6" | sampleIdx == "Bulk_GATK" ) %>%
  filter(sampleIdx %in% c("46", "43", "41", "32", "13", "22", "Bulk_GATK"))

library(aCGH)
inputPath = baseDir
outputPath=paste(inputPath, "/MergeLevels/", sep="")

if(F){
  cnv <- list();
  pdf("MergeLevels.pdf",width=16, height=6)
  for(n in 1:nrow(inputDF)) {
    name = inputDF$name[n]
    fileName = paste(inputDF$name[n], ".", name.file.postfix, sep="")
    cnv[[n]] <- read.table(paste(inputPath, "/", fileName, sep=""),header = T, stringsAsFactors = F)
    data = cnv[[n]]
    
    data.log2=data$log2ratio
    data.log2.seg.mean=data$log2seg.mean
    data.merged=mergeLevels(data.log2,data.log2.seg.mean,ansari.sign=0.1) ##ansari.sign is not default value
    seg.mean.LOWESS.MergeLevels=2^(data.merged$vecMerged)
    cn.seg.mean.MergeLevels=seg.mean.LOWESS.MergeLevels*(mean(data$cn.ratio)/mean(data$lowratio))
    data.new=cbind(data[,1:7],seg.mean.LOWESS.MergeLevels ,data[,9], cn.seg.mean.MergeLevels, 
                   round(cn.seg.mean.MergeLevels),data[,12], log2(seg.mean.LOWESS.MergeLevels), log2(round(seg.mean.LOWESS.MergeLevels)), data[,15])
    # replaces : seg.mean.LOWESS, cn.seg.mean, copy.number, log2seg.mean, log2copy.number
    # invariants : cn.ratio, log2ratio
    colnames(data.new)=colnames(data)
    plot(data.new$log2ratio,pch=20,col='gray')
    points(data.log2.seg.mean,col="red",pch=20)#old seg mean
    points(data.new$log2seg.mean,col="blue",pch=20) #new 
    title(name, outer=F)
    
    write.table(data.new,file=paste(outputPath,fileName,sep=''),quote = F, sep='\t',row.names = F,col.names = T)
  }
  dev.off()
}

## draw CNA plots
inputPath = baseDir #use default CNV result
#inputPath = paste(baseDir, "/MergeLevels/", sep="") #use MergeLevel result

##setting parameters
plot_seg = T
plot_cn = F
col.data="#CCCCCC"
col.cn="black"
col.seg="blue"#"blue" #"black"

y.at <- c(1.000, 2.000, 4.000,8)
y.labels <- c("1", "2", "4", "8") #c("0.0625","0.125","0.25","0.5","1", "2", "4", "8", "16","32","64")

cnv <- list();
for(n in 1:nrow(inputDF)) {
  cnv[[n]] <- read.table(paste("../", inputPath, "/SD", as.character(undo.sd), "/", inputDF$name[n], ".", name.file.postfix, sep=""),header = TRUE, stringsAsFactors = FALSE)
}
names(cnv) = inputDF$name
dataForPlot = bind_rows(cnv, .id = 'name') %>% left_join(inputDF)


## prepare plot
chr <- cnv[[1]]$chrom
chr.shift <- c(chr[-1], chr[length(chr)])
vlines <- c(1, cnv[[1]]$abspos[which(chr != chr.shift) + 1], cnv[[1]]$abspos[nrow(cnv)])
chr.text <- c(1:22, "X", "Y")
vlines.shift <- c(vlines[-1], 4*10^9)
chr.at <- vlines + (vlines.shift - vlines) / 2
chr.at <- chr.at[c(1:24)]
hlines <- c(1.0, 2.0, 4.0)

## plot start
pal = brewer.pal(n = 8, name = "Dark2")
col.use = c("black", pal)
palette(col.use)

inputDF$name
yrange = c(0.25,20)
## plot using ggplots

targetChrom = str_split(thisGene, "_")[[1]][1]
targetPos = as.integer(str_split(thisGene, "_")[[1]][2])

dataForPlot = dataForPlot %>% filter(chrom == targetChrom)


p = ggplot(data = dataForPlot, aes(group = name))+
  geom_point_rast(size = 0.3, aes(x = chrompos, y = cn.ratio), colour = col.data) + 
  geom_line(size = 0.5, aes(x = chrompos, y = cn.ratio), colour = col.data) + 
  geom_point_rast(size = 1, aes(x = chrompos, y = cn.seg.mean, colour = color)) + #cn.seg.mean
  geom_line(size = 1, aes(x = chrompos, y = cn.seg.mean, colour = color)) + #cn.seg.mean
  xlab("Genome Position (Gbp)") + ylab("Copy Number") + 
  facet_grid(rows = vars(sampleIdx), scales = "free") + 
  scale_colour_manual(values=c("black", "black")) +
  #scale_colour_manual(values=c(brewer.pal(n = 8, name = "Dark2")[3]), guide=FALSE) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        #strip.text.y = element_blank(), #facet_grid lables
        axis.line = element_line(colour = "black")) + 
  #scale_x_continuous(breaks = chr.at, labels = chr.text) + geom_vline(xintercept = vlines) + #geom_hline(yintercept = hlines) + #label chr position instead of absolute position
  geom_vline(xintercept = targetPos) +
  #geom_vline(xintercept = 1099675) +
  #geom_vline(xintercept = 1605099)+
  #scale_x_continuous(breaks = x.at, labels = x.labels) +
  scale_y_continuous(limits = yrange, breaks = y.at, trans = "log10"); print(p)
#expand_limits(x=0, y=0) + 
#scale_x_continuous(expand = c(0,0), breaks = x.at, labels = x.labels) +
#scale_y_continuous(expand = c(0,0), limits = yrange, breaks = y.at, trans = "log10"); print(p)
ggsave(paste("../", thisGene, "_CNA.pdf", sep=""), width=20, height=length(inputDF$name) * 1.3, units = c("cm"), dpi = 300)



##### heatmap drawing
##setting parameters
threshold.amp = log(1.01, base=2)
threshold.del = log(0.99, base=2)
cn.ref = 2.0
bin.count = "7k"
undo.sd = 0.8
sampling.Count = 479087
#pal = brewer.pal(n = 8, name = "Set1")
#col.use = c("black", pal[2], pal[3], pal[4])
color.vector = c("#C94D98","#5054A0","#4888C8","#79B752","#F8DA30","#CE3037","#D77AB2","#7C7FB8","#76A6D6",
                 "#9BC97D","#FAE364", "#DA6469")
col.use = c("black", color.vector[3], color.vector[6], color.vector[4], color.vector[9], color.vector[12])
palette(col.use)

set.seed(3)
inputDF = df_for_cna %>% filter(geneName == thisGene) %>%
  mutate(kata = ifelse(color == "#4292C6", "kata", "notkata"))


targetChrom = str_split(thisGene, "_")[[1]][1]
targetPos = as.integer(str_split(thisGene, "_")[[1]][2])

name.file.postfix = paste("uniqMap", bin.count, "k150.data.txt", sep=".")
cnv.raw <- list()
for(n in 1:nrow(inputDF)) {
  cnv.raw[[n]] <- read.table(paste("../", inputPath, "/SD", as.character(undo.sd), "/",inputDF$name[n], ".", name.file.postfix, sep=""), stringsAsFactors=F, header=T) %>%
    mutate(cn.gain.loss.log = log(cn.seg.mean / cn.ref,base=2)) %>%
    filter(chrom == targetChrom)
  #cnv.raw[[n]]$cn.heatmap = cnv.raw[[n]]$cn.gain.loss.log
  cnv.raw[[n]]$cn.heatmap = sapply(cnv.raw[[n]]$cn.gain.loss.log, function(x) max(-3.595, x)) #remove data with too low value
  cnv.raw[[n]]$cn.heatmap = sapply(cnv.raw[[n]]$cn.heatmap, function(x) if((x>0)&&(x<threshold.amp)) {0} else {x}) #only report significantly amplifed value
  cnv.raw[[n]]$cn.heatmap = sapply(cnv.raw[[n]]$cn.heatmap, function(x) if((x<0)&&(x>threshold.del)) {0} else {x}) #only report significantly amplifed value
  if(n == 1) {
    cnv.matrix=t(data.matrix(cnv.raw[[n]]$cn.heatmap))
  } else {
    cnv.matrix = rbind(cnv.matrix, t(data.matrix(cnv.raw[[n]]$cn.heatmap)))
  }
}

rownames(cnv.matrix) = inputDF$name 
colnames(cnv.matrix) = paste0("bin", seq_len(ncol(cnv.matrix)))

##prepare for plot
chr <- cnv.raw[[1]]$chrom
chr.shift <- c(chr[-1], chr[length(chr)])
vlines <- c(1, which(chr != chr.shift) + 1, nrow(cnv.raw[[1]]))
chr.text <- c(1:22, "X", "Y")
vlines.shift <- c(vlines[-1], 4*10^9)
chr.at <- vlines + (vlines.shift - vlines) / 2
chr.at <- chr.at[c(1:24)]
hlines <- c(1.0, 2.0, 4.0)


## draw complex heatmap


ha_sample <- rowAnnotation(k = inputDF$kata,
                           col = list(k = setNames(c(RColorBrewer::brewer.pal(9, "Blues")[6], "#C7C7C7"), 
                                                          inputDF$kata %>% unique)))

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf("heatmap.pdf",width=20/2.54, height=11/2.54) 
Heatmap(cnv.matrix,
        name = "Fold Amplification", #title of legend
        row_names_gp = gpar(fontsize = 7),show_row_names =F,
        column_split = cnv.raw[[1]]$chrom %>% str_replace_all(c("X" = "23", "Y" = "24")) %>% as.numeric(),
        column_gap = unit(1, "mm"),
        cluster_columns = F,
        split = inputDF$kata,
        use_raster = T,
        left_annotation = ha_sample,
        col = col_fun,
        #cluster_rows = T,
        #row_order = with(inputDF, order(groupTemplate,-AUC)),
        column_order = order(as.numeric(gsub("bin", "", colnames(cnv.matrix)))),
        show_column_dend = FALSE, show_column_names =F
)
dev.off()





}

