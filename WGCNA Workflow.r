library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()

setwd('E:/WGCNA')

### Read data
data<- read.csv('Counts.csv', row.names = 1)
Sample_info<- read.csv('Sample_info.csv', row.names = 1)


### Quality control
#### Step1: Detect outlier genes with WGCNA package
gg<- goodSamplesGenes(t(data))
summary(gg)
gg$allOK

#Since it indicates tat all genes are not good genes let's quantify the outlier genes
#We'll do it for samples as well
table(gg$goodGenes)
table(gg$goodSamples)

#Around 4475 outlier genes but no outlier samples. Let's remove those genes
data<- data[gg$goodGenes=='TRUE',]

#### Step2: Identify outlier samples with hierachical clustering
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

#### Step3: Identify outlier samples with PCA
pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

#N.B. We need to remove batch effect if there's any at this stage before we proceed further

#In both these methods we found Sample 3, 15, 16, 21 and 22 are distantly related 
#So we'll exclude them

discard_samples<- c('Sample5')
data.final<- data[,!colnames(data) %in% discard_samples]
Sample_info_final<- Sample_info[!row.names(Sample_info) %in% discard_samples,]

### Perform VST normalization with DESeq2 since it's a count matrix
#N.B. FPKM or RPKM matrices need to be log transformed
#Before performing normalization make sure if the rownames in matrix and colnames in
#Sample info are identical

all(rownames(Sample_info_final)%in%colnames(data.final)) #If present
all(rownames(Sample_info_final)==colnames(data.final)) #If present in same order

#### Create deseq2 object

dds <- DESeqDataSetFromMatrix(countData = data.final,
                              colData = Sample_info_final,
                              design = ~ 1) #Specifying no design
nrow(dds)
#### Remove samples that have <15 counts in more than 75% of the samples
#Recommended by WGCNA
# keep<- keep <- dds[rowSums(counts(dds) >= 15)>=35,]  
#Since this filtering generates extremely low read counts, we'll follow modified filter.
#Moreover, the downstream analysis with higher # of genes  will depend largely on 
#the available ram of PC. So, we'll consider gene numbers based on our PC capacity.

keep <- dds[rowSums(counts(dds) >= 15)>=39,] #75% of 52 samples is 52*0.75=39
nrow(keep) #5075 genes

#### Perform VST transformation
dds.vst<- vst(keep) 
dds.vst<- dds.vst[1:8000,] #We reduced the # of genes to 8000 because of low RAM
nrow(dds.vst)

#### get normalized counts
norm.counts <- assay(dds.vst) %>% 
  t()

### Constructing network
#### Set softhreshold power
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

#### Set the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)
#N.B. unsigned -> nodes with positive & negative correlation are treated equally 
#N.B. signed -> nodes with negative correlation are considered *unconnected*, treated as zero

sft.data <- sft$fitIndices

#### visualization to pick power

Plot1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


Plot2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(Plot1, Plot2, nrow = 2)

#### convert matrix to numeric and assign power
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor


#### Estimate w.r.t blocksize memory depends on available ram of PC
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 8000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)
# grey module = all genes that doesn't fall into other modules were assigned to the grey module


### Finding the association between module and traits
#### Prepare and manipulate data for association analysis
#### Converting categorical variable into binary group variables
type <- Sample_info_final %>% 
        mutate(type = ifelse(grepl('Lung cancer', Group), 1, 0)) %>% 
        select(4)



Sample_info_final$Status <- factor(Sample_info_final$Status, 
                                   levels = c("None", "Newly diagnosed", 
                                   "Post Chemotherapy", "On Doxorubicin"))

Status<- binarizeCategoricalColumns(Sample_info_final$Status,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)


Status.final<- cbind(type, Status)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

status.module.corr <- cor(module_eigengenes, Status.final, use = 'p')
status.module.corr.pvals <- corPvalueStudent(status.module.corr , nSamples)


### visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, Status.final, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
   column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[14:17],
             y = names(heatmap.data)[1:13],
             col = c("blue1", "green", "pink", "red"))



module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'blue') %>% 
  rownames()



### Identifying significant genes against different features of data

# Calculating the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:10,1:10]


### Calculating the gene significance and associated p-values against a particular trait

gene.signf.corr <- cor(norm.counts, traits$`data.Post Chemotherapy.vs.all`, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

#Sub-setting based on significant association: X1:corr, X2:pval
Sig_chemo<-data.frame(cbind(gene.signf.corr, gene.signf.corr.pvals))
Sig_chemo_final<- filter(Sig_chemo, Sig_chemo$X1> 0.3 | Sig_chemo$X1< -0.3
                         & Sig_chemo$X2<0.05)  

write.csv(Sig_chemo_final, 'sig_genes_chemo.csv')


