#in version 4.4.1 R
rm(list = ls())

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(plyr)
library(pheatmap)
library(ggh4x)
library(ggrepel)

library(readxl)
library(readr)
library("QFeatures")
#library(edgeR)
library("limma")
library("here")
library("NormalyzerDE")
library("factoextra")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
#library("patchwork")
library(MsCoreUtils)
library(tidyr)
library(ggVennDiagram)

theme_flo <-function() {
  theme_classic() +
    theme(
      axis.line=element_line(size=0.4),
      plot.title=element_text(hjust=0.5, face = 'bold', size = 10),
      axis.title=element_text(face = 'bold',  size = 10),
      axis.text =element_text(color = 'black', face= 'plain', size = 8),
      legend.margin = margin(0.2),
      legend.text = element_text(size = 8),
      legend.key.height = (unit(0.5, 'cm')),
      legend.key.width = (unit(0.5, 'cm')),
      legend.title=element_text(face = 'bold', size = 8),
      legend.direction = ("vertical"),
      legend.box = ("vertical"),
      legend.position = c(0.9, 0.9),
      axis.ticks.length= unit(0.1, "cm"),
      axis.ticks = element_line(size = 0.6),
      panel.spacing = unit(0.7, "lines"),
      strip.background = element_rect(size = 0.8, colour= NA, fill=NA),
      strip.placement = ("outside"),
      strip.text.x= element_text(hjust=0.05, size = 10, face="bold", margin = margin(l = 0.1))
    )
}

#import the tsv file with the intensity results here
df<-readr::read_tsv("C:\\Users\\Florence\\OneDrive - University of Cambridge\\Proteomics\\IP_results\\IPsamples.pg_matrix.tsv")
df1 <- df
#Remove columns with duplicated samples
df1_Rep2 <- dplyr::select(df1, -c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60))
df1_Rep2 <- df1_Rep2[!duplicated(df1_Rep2$Genes), ] #remove duplicated proteins

#update column names
colnames(df1_Rep2)[5:164] <- substring(names(df1_Rep2[5:164]), 63) 
colnames(df1_Rep2) <- gsub(".raw", "", names(df1_Rep2)) 
colnames(df1_Rep2) <- gsub("STR", "PUT", names(df1_Rep2)) #original samples were labelled as striatum, updated to putamen
df2 <- df1_Rep2 

colnames(df2)[5:164] <-  sub("_1A", "_1", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_2025", "_REP", colnames(df2)[5:164])

colnames(df2)[5:164] <-  sub("_3", "_PD5", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_4", "_PD13", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_5", "_PD14", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_6", "_CON15", colnames(df2)[5:164]) # 6 is T: order change during MS run
colnames(df2)[5:164] <-  sub("_7", "_PD1", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_8", "_PD7", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_9", "_PD9", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_10", "_PD10", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_11", "_PD12", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_12", "_CON4", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_13", "_CON5", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_14", "_CON6", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_15", "_CON7", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_16", "_CON8", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_17", "_CON10", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_18", "_CON12", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_19", "_CON13", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_20", "_CON14", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_1", "_PD2", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_2", "_PD3", colnames(df2)[5:164])

#Rep 2 only
colnames(df2)[5:164] <-  sub("PD3B", "PD3", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD5C", "PD5", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD13D", "PD13", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD14E", "PD14", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON15T", "CON15", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD1F", "PD1", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD7G", "PD7", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD9H", "PD9", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD10I", "PD10", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("PD12J", "PD12", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON4K", "CON4", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON5L", "CON5", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON6M", "CON6", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON7N", "CON7", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON8O", "CON8", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON10P", "CON10", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON12Q", "CON12", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON13R", "CON13", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("CON14S", "CON14", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321191235", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321193040", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321194904", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321200711", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321202520", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321204331", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321210143", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321211952", "", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("_REP0321213801", "", colnames(df2)[5:164])

colnames(df2)[5:164] <-  sub("1_", "PUT_F1_Syn211_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("2_", "PUT_F1_IgG_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("3_", "PUT_F2_Syn211_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("4_", "PUT_F2_IgG_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("5_", "FC_F1_Syn211_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("6_", "FC_F1_IgG_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("7_", "FC_F2_Syn211_", colnames(df2)[5:164])
colnames(df2)[5:164] <-  sub("8_", "FC_F2_IgG_", colnames(df2)[5:164])

df4 <- df2

#write.csv(df4, "D:/Newgraph/CoIP_rawdata.csv") #export raw data for supplementary

#add the dataframe to a QFeatures object, only the quantitative columns- each of the sample columns
cc_qf_IP <- readQFeatures(assayData = df4,
                       quantCols = 5:164, #for original
                       #quantCols = 7:164, #for rerun
                       name = "proteins_raw")

#look at the rowdata to the QFeatures object: this is the non-quantitative data ie the gene descriptions. These are also called features
cc_qf_IP[["proteins_raw"]] %>%
  rowData() %>%
  names()

#look at the column data: the names of the quantitative columns
cc_qf_IP[["proteins_raw"]] %>%
  colData() %>% 
  rownames()

#Import metadata
my_meta_data <- read_excel("C:/Users/Florence/OneDrive - University of Cambridge/Proteomics/IP_results/IPmeta_data.xlsx", 
                           sheet = "metadata_rep2")
metadata_df <- data.frame(my_meta_data)

## Annotate colData with condition information
cc_qf_IP$Case <- metadata_df$Case
cc_qf_IP$Fraction <- metadata_df$Fraction
cc_qf_IP$Sample_Type <- metadata_df$Sample_Type
cc_qf_IP$Region <- metadata_df$Region
cc_qf_IP$Antibody <- metadata_df$Antibody
colData(cc_qf_IP)
colData(cc_qf_IP[["proteins_raw"]]) <- colData(cc_qf_IP)

## Filtering ##
## Extract a copy of the raw protein data
raw_data_copy <- cc_qf_IP[["proteins_raw"]] 

## Re-add the assay to our QFeatures object with a new name
cc_qf_IP <- addAssay(x = cc_qf_IP, 
                  y = raw_data_copy, 
                  name = "proteins_filtered")

#Make a dataframe with the raw numbers and a column for the Genes
df_proteins <- 
  cc_qf_IP[["proteins_filtered"]] %>% 
  assay() %>%
  as.data.frame() %>%
  mutate(Genes = rowData(cc_qf_IP[["proteins_filtered"]])$Genes) 

dim(df_proteins)

#see how many missing values there are before filtering out
#see how many missing values there are before filtering out
mv_filtered <- nNA(cc_qf_IP, i = "proteins_filtered")
mv_filtered_df <- mv_filtered$nNAcols %>% as_tibble() 
mv_filtered_df <- mv_filtered_df %>% separate(name, c("Region", "Fraction", "Antibody", "Case"), sep = "_")
mv_filtered_df$Fraction <- factor(mv_filtered_df$Fraction, levels = c("F1", "F2"))
mv_filtered_df$Sample_Type <- ifelse(grepl("PD",mv_filtered_df$Case),"PD","CON")

#bar plot version
missing <- mv_filtered_df %>%
  ggplot(aes(y = Case, x = pNA,  fill = Sample_Type, alpha = Antibody)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  facet_wrap2( ~ Fraction + Region, ncol = 1, nrow = 4, axes = "all", strip.position = "top") + 
  labs(y = "Sample", x = "Proportion missing values") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by =0.2)) +
  scale_fill_manual(values = c("black", "#d40000ff")) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  scale_alpha_manual(values = c(0.3, 0.8)) +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 5.5))

ggsave(missing, filename = file.path('D:/Newgraph', 'missing_prefilter_CoIP_F1.pdf'), units = 'mm', height = 250, width =180, dpi = 300)

#grouped box plot version
missing <- mv_filtered_df %>%
  ggplot(aes(y = interaction(Antibody, Sample_Type), x = pNA, alpha = Antibody, color = Sample_Type)) +
  geom_boxplot(position = position_dodge(0.9), width = 0.7) +
  geom_point(size= 1.2, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9)) +
  facet_wrap2( ~ Fraction + Region, ncol = 1, nrow = 4, axes = "all", strip.position = "top") + 
  labs(y = "Proportion missing values", x = "") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by =0.2)) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  #scale_color_manual(values = c( '#A9A9A9',  "#276ab3")) +
  #scale_alpha_manual(values = c(0.3, 0.8)) +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8))

ggsave(missing, filename = file.path('D:/Newgraph', 'missing_prefilter_CoIP.pdf'), units = 'mm', height = 150, width =89, dpi = 300)


df_proteins_copy <-df_proteins

#Step 1: separate the total dataframe into chunks to remove highly missing proteins
#Filter out proteins that are missing in > 6/10 samples per group
dfa <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("PUT_F1_Syn211_PD")))) >= 4)   
dfb <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("PUT_F2_Syn211_PD")))) >= 4)   
dfc <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("PUT_F1_Syn211_CON")))) >= 4)    
dfd <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("PUT_F2_Syn211_CON")))) >= 4)     
dfe <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("FC_F1_Syn211_PD")))) >= 4)   
dff <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("FC_F2_Syn211_PD")))) >= 4)     
dfg <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("FC_F1_Syn211_CON")))) >= 4)   
dfh <- df_proteins_copy %>% filter(rowSums(!is.na(across(starts_with("FC_F2_Syn211_CON")))) >= 4)   

df_proteins_BP <- dfa %>% 
  merge(dfb, all = T) %>% 
  merge(dfc, all =T) %>%
  merge(dfd, all = T) %>%
  merge(dfe, all = T) %>%
  merge(dff, all = T) %>%
  merge(dfg, all = T) %>%
  merge(dfh, all = T)

df_proteins_BP <- df_proteins_BP[!duplicated(df_proteins_BP$Genes), ] #checks for duplicated proteins
df_MAPT <- subset(df_proteins_BP, Genes == "MAPT;cRAP-MAPT") #keep MAPT in the dataset
df_proteins_BP <- df_proteins_BP %>% filter(!grepl("KRT", Genes)) 
df_proteins_BP <- df_proteins_BP %>% filter(!grepl("cRAP", Genes)) 
df_proteins_BP <- df_proteins_BP %>% filter(!grepl("IG", Genes)) #likely to be contaminant proteins from the Co-IP antibody
df_proteins_BP <- df_proteins_BP %>% filter(!grepl("B4GALT1", Genes)) #this was higher in Syn211 vs IgG, from the azido modification kit used to biotinylate Syn211
df_proteins_BP <- df_proteins_BP %>% filter(!is.na(Genes))

df_proteins_filt <- rbind(df_proteins_BP, df_MAPT)

dim(df_proteins_filt) #filtered df
#df_proteins_filt$Genes <- gsub("MAPT;cRAP-MAPT", "MAPT", df_proteins_filt$Genes)

#Join this to rowData
df_proteins_filt['marker'] = 1 #Make a column to show which genes are present in unfiltered rowData

#Put the rowdata in a df- it includes the Gene names
protein_rowdata <- 
  cc_qf_IP[["proteins_filtered"]] %>% 
  rowData() %>% 
  as.data.frame()

#combine the rowdata with the filtered protein data frame, the rows will match up based on the Gene names and all of the proteins
#that have been filtered out on the above criteria will have a NA in their Marker column
mergedrowData <- left_join(x = protein_rowdata,
                           y = df_proteins_filt,
                           by = c("Genes" = "Genes"))

names(mergedrowData)
## Now add the Marker column to the proteins_filtered assay and filter out all of the proteins that don't have a 1 Marker
rowData(cc_qf_IP[["proteins_filtered"]])$Missing <- mergedrowData$marker

cc_qf_IP <- cc_qf_IP %>%
  filterFeatures(~ Missing == 1, ##keep features with this feature 
                 i = "proteins_filtered", na.rm = TRUE)

head(rowData(cc_qf_IP[["proteins_filtered"]]))

mv_filtered <- nNA(cc_qf_IP, i = "proteins_filtered")
mv_filtered_df <- mv_filtered$nNAcols %>% as_tibble() 
mv_filtered_df <- mv_filtered_df %>% separate(name, c("Region", "Fraction", "Antibody", "Case"), sep = "_")
mv_filtered_df$Fraction <- factor(mv_filtered_df$Fraction, levels = c("F1", "F2"))
mv_filtered_df$Sample_Type <- ifelse(grepl("PD",mv_filtered_df$Case),"PD","CON")

#grouped box plot version
missing <- mv_filtered_df %>%
  ggplot(aes(y = interaction(Antibody, Sample_Type), x = pNA, alpha = Antibody, color = Sample_Type)) +
  geom_boxplot(position = position_dodge(0.9), width = 0.7) +
  geom_point(size= 1.2, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), shape = 16) +
  facet_wrap2( ~ Fraction + Region, ncol = 1, nrow = 4, axes = "all", strip.position = "top") + 
  labs(y = "Proportion missing values", x = "") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by =0.2)) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  #scale_color_manual(values = c( '#A9A9A9',  "#276ab3")) +
  #scale_alpha_manual(values = c(0.3, 0.8)) +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8))

ggsave(missing, filename = file.path('D:/Newgraph', 'missing_postfilter_CoIP.pdf'), units = 'mm', height = 150, width =89, dpi = 300)


cc_qf_IP <- impute(object = cc_qf_IP,
                method = "MinDet", 
                i = "proteins_filtered",
                MARGIN = 1, # Intentionally keeping as 1 to capture stochastic, absent/present per protein = 1
                name = "proteins_imputed")


cc_qf_IP[["proteins_imputed"]] %>% assay() %>% as_tibble()


cc_qf_IP[["proteins_imputed"]] %>%
  assay() %>%
  longFormat() %>%
  ggplot(aes(x = log2(value))) +
  geom_histogram() + 
  theme_bw() +
  xlab("(Log2) Abundance")

cc_qf_IP <- logTransform(object = cc_qf_IP, 
                      base =2, 
                      i = "proteins_imputed", 
                      name = "log_proteins")


cc_qf_IP <- normalize(cc_qf_IP, 
                   i = "log_proteins", 
                   name = "log_proteins_centered",
                   method = "center.median") 

df<- cc_qf_IP[["log_proteins_centered"]] %>% assay() %>% as_tibble()
df<- cbind(rowData(cc_qf_IP[["log_proteins_centered"]]), df)
write.csv(df, "D:/Newgraph/CoIP_filteredpros.csv") #export data

#Plot LFQ per case
prenorm_df <- cc_qf_IP[["log_proteins"]] %>%
  assay() %>%
  longFormat()

#prenorm_df$Group <- ifelse(grepl("CON", prenorm_df$colname), "CON", "PD")
#prenorm_df$Antibody <- ifelse(grepl("IgG", prenorm_df$colname), "IgG", "Syn211")
#prenorm_df$Region <- ifelse(grepl("FC", prenorm_df$colname), "FC", "PUT")
prenorm_df$Sample_Type <- ifelse(grepl("CON", prenorm_df$colname), "CON", "PD")
prenorm_df <- prenorm_df %>% separate(colname, c("Region", "Fraction", "Antibody", "Case"), sep = "_")


pre_norm <- prenorm_df %>%
  ggplot(aes(y = interaction(Antibody, Case, Region, Fraction), x = value, color = Sample_Type)) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  #scale_alpha_manual(values = c(0.4, 1)) +
  #facet_wrap2( ~ Fraction, ncol = 1, nrow = 2, axes = "all", strip.position = "top") + 
  geom_boxplot(outlier.size = 1.2, position = position_dodge(0.9), width = 0.7) +
  labs(y = "Sample", x = "log2 protein abundance")+
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 5))
pre_norm

ggsave(pre_norm, filename = file.path('D:/Newgraph', 'prenorm_LFQ.pdf'), units = 'mm', height = 220, width =89, dpi = 300)

post_norm_df <- cc_qf_IP[["log_proteins_centered"]] %>%
  assay() %>%
  longFormat()

#prenorm_df$Region <- ifelse(grepl("FC", prenorm_df$colname), "FC", "PUT")
post_norm_df$Sample_Type <- ifelse(grepl("CON", post_norm_df$colname), "CON", "PD")
post_norm_df <- post_norm_df %>% separate(colname, c("Region", "Fraction", "Antibody", "Case"), sep = "_")

post_norm <- post_norm_df %>%
  ggplot(aes(y = interaction(Antibody, Case, Region, Fraction), x = value, color = Sample_Type)) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  #scale_alpha_manual(values = c(0.4, 1)) +
  #facet_wrap2( ~ Fraction, ncol = 1, nrow = 2, axes = "all", strip.position = "top") + 
  geom_boxplot(outlier.size = 1.2, position = position_dodge(0.9), width = 0.7) +
  labs(y = "Sample", x = "log2 protein abundance")+
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 5))
post_norm

ggsave(post_norm, filename = file.path('D:/Newgraph', 'post_norm_LFQ.pdf'), units = 'mm', height = 220, width =89, dpi = 300)


########### Block 1: identify proteins that are significantly increased in 211 vs IgG, start with both fractions for one-sided volcano plot #######
all_proteins <- cc_qf_IP[["log_proteins_centered"]]
row_names <-  all_proteins %>% rowData() %>% as_tibble() %>% pull(Genes) %>% as.list()

all_proteins <- "rownames<-"(all_proteins, row_names)
all_proteins$Sample_Type <- factor(all_proteins$Sample_Type , levels = c("CON", "PD"))
all_proteins$Region <- factor(all_proteins$Region , levels = c("FC", "PUT"))
all_proteins$Antibody <- factor(all_proteins$Antibody , levels = c("IgG", "211"))

Case <- all_proteins$Case
Sample_Type <- all_proteins$Sample_Type 
Fraction <- all_proteins$Fraction
Antibody <- all_proteins$Antibody
Region <- all_proteins$Region

#make a joint block with Region, Fraction, Case
df <- data.frame(Case= Case, Fraction = Fraction, Region = Region)
block_vars_ch <- paste(df$Case, df$Fraction, df$Region, sep = "_") 

#see if there is a duplicate correlation for case
Antibody <- factor(all_proteins$Antibody)
design <- model.matrix(~ 0 + Antibody) 
colnames(design) <- c("IgG", "S211")
dupcor <- duplicateCorrelation(assay(all_proteins),design,block= block_vars_ch)
dupcor$consensus.correlation #small positive correlation

fit <- lmFit(assay(all_proteins), design, block = block_vars_ch, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(S211 - IgG, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)
summary(decideTests(final_model, method="global"))

## Format results
limma_results <- topTable(fit = final_model,
                          coef = c(1), 
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

limma_results <- limma_results %>%
  mutate(direction = ifelse(logFC > 0, "Up", "Down"),
         significance = ifelse(adj.P.Val < 0.05, "sig", "not.sig"),
         fillcol = ifelse(direction == "Up" & significance == "sig", "Syn211",
                          ifelse(direction == "Up" & significance == "not.sig", "No Change",
                                 ifelse(direction == "Down" & significance == "sig", "IgG",
                                        "No Change"))))

enrichedinSyn211 <- subset(limma_results, logFC > 0.5) %>% dplyr::select(Protein, logFC, adj.P.Val) %>% as_tibble()
enrichedinIgG <- subset(limma_results, logFC < -0.5) %>% dplyr::select(Protein, logFC, adj.P.Val) %>% as_tibble()
limma_results$fillcol <- factor(limma_results$fillcol, c("No Change", "IgG", "Syn211"))

### calculate the mean expression for each fraction and merge with Limma results
all_proteins_df <- cc_qf_IP[["log_proteins_centered"]] %>% assay() %>% as_tibble()

# Select columns containing "weight"
cols <- grepl("Syn211", names(all_proteins_df))
df_F1<- all_proteins_df[, cols]
df_F1$F1_mean <- rowMeans(df_F1[,1:80])

cols <- grepl("IgG", names(all_proteins_df))
df_F2<- all_proteins_df[, cols]
df_F2$F2_mean <- rowMeans(df_F2[,1:80])

genes <-  all_proteins %>% rowData() %>% as_tibble() 
genes <- dplyr::select(genes, c("Genes"))
mean_intensity <- cbind(genes, df_F1$F1_mean,df_F2$F2_mean )
mean_intensity$Genes <- as.character(mean_intensity$Genes)
limma_results$Protein <- as.character(limma_results$Protein)
merged <- left_join(mean_intensity, limma_results, by = c("Genes" = "Protein"))

write.csv(merged, "D:/Newgraph/CoIP_limma_results.csv") #export results for supplementary

#for 211 vs IgG: 
S211_volcano <- limma_results %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = fillcol)) +
  ylab("")+ #-log10 adjusted p-value
  xlab('') + #log2 fold change
  scale_x_continuous(expand = c(0,0), limits = c(0, 3), breaks = seq(0, 3, by =0.5))+
  scale_y_continuous(expand = c(0,0), limits=c(0, 30), breaks = seq(0, 30, by =5)) +
  scale_alpha_manual(values = c(0.4,1)) +
  geom_point(data = subset(limma_results, adj.P.Val >= 0.05 & logFC > -0.5 | logFC < 0.5 ), size = 1.5, alpha = 0.4, color = "grey", shape = 16) +
  geom_hline(yintercept = 1.301, color = "black", linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "black", , linetype = "dashed", alpha = 0.4) +
  geom_point(data = subset(enrichedinIgG, adj.P.Val <= 0.05), size = 1.5, color = "#276ab3", fill = "#276ab3",alpha = 0.4, shape = 21) + #shape = 16
  geom_point(data = subset(enrichedinSyn211, adj.P.Val <= 0.05), size = 1.5, color = "#276ab3", fill = "#276ab3", alpha = 0.4, shape = 21) +
  geom_text_repel(data = subset(enrichedinSyn211, adj.P.Val <= 0.05), aes(label = Protein), size = 2.5, color = "#276ab3") +
  guides(color= guide_legend(title = '', override.aes = aes(label = "")), alpha = "none") +
  theme_flo()+
  theme(legend.box = ("vertical"), legend.position = c(0.88,0.92)) 

ggsave(S211_volcano, filename = file.path('D:/Newgraph', 'S211_volcano_nocentre.pdf'), units = 'mm', height = 70, width = 90, dpi = 300)


############### BLOCK 2: compare binding partners in F1 and F2 ###################
#now look at the fractions for the venn diagram
all_proteins <- cc_qf_IP[["log_proteins_centered"]]
row_names <-  all_proteins %>% rowData() %>% as_tibble() %>% pull(Genes) %>% as.list()

all_proteins <- "rownames<-"(all_proteins, row_names)
all_proteins$Sample_Type <- factor(all_proteins$Sample_Type , levels = c("CON", "PD"))
all_proteins$Region <- factor(all_proteins$Region , levels = c("FC", "PUT"))
all_proteins$Antibody <- factor(all_proteins$Antibody , levels = c("IgG", "211"))

Case <- all_proteins$Case
Sample_Type <- all_proteins$Sample_Type 
Fraction <- all_proteins$Fraction
Antibody <- all_proteins$Antibody
Region <- all_proteins$Region

#First do at the Fraction level only
Antibody <- factor(all_proteins$Antibody)
design <- model.matrix(~ 0 + Antibody:Fraction) 
colnames(design) <- c("IgG_F1", "S211_F1", "IgG_F2", "S211_F2")
dupcor <- duplicateCorrelation(assay(all_proteins), design, block = block_vars_ch)
dupcor$consensus.correlation #small positive correlation

fit <- lmFit(assay(all_proteins), design, block = block_vars_ch, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(F1 = (S211_F1 - IgG_F1), 
                           F2 = (S211_F2 - IgG_F2), levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)
vennums <- summary(decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5)) # use separate = same as using each coefficient on its own
vennDiagram(decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5), include=c("up"), circle.col = c("#276ab3", "#276ab3"))

m <- decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5)
#write.csv(m, "D:/Newgraph/CoIP_filteredpros.csv") #export data for supplementary

#Next do with Fraction:Region level
Antibody <- factor(all_proteins$Antibody)
design <- model.matrix(~ 0 + Antibody:Fraction:Region) 
colnames(design) <- c("IgG_F1_FC", "S211_F1_FC", "IgG_F2_FC", "S211_F2_FC", "IgG_F1_PUT", "S211_F1_PUT", "IgG_F2_PUT", "S211_F2_PUT")
dupcor <- duplicateCorrelation(assay(all_proteins), design, block = block_vars_ch)
dupcor$consensus.correlation #small positive correlation


fit <- lmFit(assay(all_proteins), design, block = block_vars_ch, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(
                            F1_PUT = (S211_F1_PUT - IgG_F1_PUT), 
                            F1_FC = (S211_F1_FC - IgG_F1_FC),
                            F2_PUT = (S211_F2_PUT - IgG_F2_PUT), 
                            F2_FC = (S211_F2_FC - IgG_F2_FC) , levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)
vennums <- summary(decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5)) # use separate = same as using each coefficient on its own
vennDiagram(decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5), include=c("up"), circle.col = c("#276ab3", "#276ab3"))

m <- decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5)
#write.csv(m, "D:/Newgraph/CoIP_filteredpros.csv") #export data for supplementary

#get the 47 proteins with overlap in both fractions and regions
pros_47 <- (decideTests(final_model, method="global", adjust.method="BH",p.value=0.05, lfc = 0.5)) %>% as.data.frame()
pros_47 <- subset(pros_47, F1_PUT == 1 & F1_FC == 1 & F2_PUT == 1 & F2_FC == 1) 
pros_47$ID <- as.character(rownames(pros_47))
BP <- dplyr::select(pros_47, c(ID))

limma_results <- topTable(fit = final_model,
                          coef = c(1, 2,3,4), 
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

pros_47_FC<- subset(limma_results, (limma_results$Protein %in% pros_47$ID))
#write.csv(pros_47_FC, "D:/Newgraph/top47BPs.csv") #export data for supplementary

#GO enrichment analysis of the 47 proteins, start by using all proteins as the background
all_proteins <- cc_qf_IP[["log_proteins_centered"]]
protein_info <- all_proteins %>% 
  rowData() %>%
  as_tibble %>%
  select(ID = Genes,
         Protein_numb = Protein.Group,
         Protein_description = First.Protein.Description)

protein_info <- protein_info[!duplicated(protein_info$ID), ] #checks for duplicated proteins

BP_des <- BP %>%
  left_join(protein_info, by = "ID")

ego_BP <- enrichGO(gene = BP_des$Protein_numb,               # list of down proteins
                   universe = protein_info$Protein_numb,      # all proteins 
                   OrgDb = org.Hs.eg.db,                  # database to query
                   keyType = "UNIPROT",                   # protein ID encoding 
                   qvalueCutoff = 0.05,
                   ont = "BP",                            # can be CC, MF, BP, or ALL
                   readable = TRUE)

cnet_top47 <- 
  cnetplot(ego_BP, 
           color_category = "#276ab3",
           color_gene = "black",
           showCategory = 5, #there are 112 enriched terms
           #node_label = "category",
           alpha_gene = 0.4,
           cex_label_category = 0.8,
           #cex_label_gene = 0.6 ,
           shadowtext = "none",
           layout = "dh") +
  theme_flo() +
  labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  guides(size = "none") 
cnet_top47

ggsave(cnet_top47, filename = file.path('D:/Newgraph', 'cnet_47pros_nocentre.pdf'), units = 'mm', height = 80, width =180, dpi = 300)

ego_BP_df <- ego_BP %>% as_tibble()
ego_BP_df <- ego_BP_df[order(ego_BP_df$Count, decreasing = TRUE),]  
ego_BP_df$Description <- factor(ego_BP_df$Description , levels=unique(ego_BP_df$Description))
ego_BP_df$GeneRatio_num <- as.numeric(sapply(strsplit(ego_BP_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(ego_BP_df$GeneRatio, "/"), "[", 2))

#write.csv(ego_BP_df, "D:/Newgraph/CoIP_47BP_GO.csv") #export results for supplementary

## Make a plot for the top 47 binding partners
top47_refs <- read_excel("C:/Users/Florence/OneDrive - University of Cambridge/proteomicspaper/ResultsFigures/top47_refs.xlsx")
top47_interactions <- dplyr::select(top47_refs, c("Gene_ID", "Novel","Human_brain", "PD_genetic","PD", "model", "Genetic", "Functional"))

top47_interactions[is.na(top47_interactions)] <- 0

top47_interactions_long <- pivot_longer(top47_interactions, cols = c("Human_brain", "Novel","PD", "model", "Genetic", "Functional"), names_to = "group")
summary <- ddply(top47_interactions_long, c("group"), summarise,
      sum = sum(value))

summary$group <- factor(summary$group, levels = c("Human_brain", "model", "Functional", "PD","Genetic", "Novel"))

bar_diffs <- ggplot(summary, aes(x=group, y=sum)) +
  labs(x="Interaction with ", 
       y= 'No. proteins') +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.6, fill = "#276ab3", color = "#276ab3") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 30)) +
  geom_text(stat = "identity", label = summary$sum, vjust = -0.5, size = 2.5) +
  theme_flo()

ggsave(bar_diffs, filename = file.path('D:/Newgraph', 'heatmap_diffs.pdf'), units = 'mm', height = 35, width =88, dpi = 300)

################# Block 3: take list of enriched proteins in 211, subset DF for these, look at presence/absence 
#now get a list of all binding partners for all regions/fractions to continue analysing
BPs <- subset(limma_results, F1_PUT >= 0.5 | F1_FC >= 0.5 | F2_PUT >= 0.5  | F2_FC >= 0.5 )
BPs <- subset(BPs, adj.P.Val <= 0.05) 
BPs$BP <- 1
BPs <- BPs[!duplicated(BPs$Protein), ] #checks for duplicated proteins

centre_data_copy <- cc_qf_IP[["log_proteins_centered"]] 
cc_qf_IP <- addAssay(x = cc_qf_IP, 
                     y = centre_data_copy, 
                     name = "log_proteins_BP")

protein_rowdata2 <- 
  cc_qf_IP[["log_proteins_BP"]] %>% 
  rowData() %>% 
  as.data.frame()

mergedrowData <- left_join(x = protein_rowdata2,
                           y = BPs,
                           by = c("Genes" = "Protein"))

rowData(cc_qf_IP[["log_proteins_BP"]])$BP_marker <- mergedrowData$BP

cc_qf_IP <- cc_qf_IP %>%
  filterFeatures(~ BP_marker == 1, ##keep features with this feature 
                 i = "log_proteins_BP", na.rm = TRUE)

all_proteins <- cc_qf_IP[["log_proteins_BP"]]
row_names <-  all_proteins %>% rowData() %>% as_tibble() %>% pull(Genes) %>% as.list()
all_proteins <- "rownames<-"(all_proteins, row_names)

Antibody <- factor(all_proteins$Antibody)
design <- model.matrix(~ 0 + Antibody:Fraction:Region:Sample_Type) 
colnames(design) <- c("S211_F1_FC_CON", "IgG_F1_FC_CON", "S211_F2_FC_CON", "IgG_F2_FC_CON", "S211_F1_PUT_CON", "IgG_F1_PUT_CON",
                      "S211_F2_PUT_CON", "IgG_F2_PUT_CON", "S211_F1_FC_PD", "IgG_F1_FC_PD", "S211_F2_FC_PD", "IgG_F2_FC_PD",
                      "S211_F1_PUT_PD", "IgG_F1_PUT_PD", "S211_F2_PUT_PD", "IgG_F2_PUT_PD")
dupcor <- duplicateCorrelation(assay(all_proteins), design, block = block_vars_ch)
dupcor$consensus.correlation #accounting for the same sample in diff fractions/antibodies etc


fit <- lmFit(assay(all_proteins), design, block = block_vars_ch, correlation = dupcor$consensus.correlation)
contrasts <- makeContrasts(
  F1_PUT_CON = (S211_F1_PUT_CON - IgG_F1_PUT_CON), 
  F1_FC_CON = (S211_F1_FC_CON - IgG_F1_FC_CON),
  F2_PUT_CON = (S211_F2_PUT_CON - IgG_F2_PUT_CON), 
  F2_FC_CON = (S211_F2_FC_CON - IgG_F2_FC_CON) ,
  F1_PUT_PD = (S211_F1_PUT_PD - IgG_F1_PUT_PD), 
  F1_FC_PD = (S211_F1_FC_PD - IgG_F1_FC_PD),
  F2_PUT_PD = (S211_F2_PUT_PD - IgG_F2_PUT_PD), 
  F2_FC_PD = (S211_F2_FC_PD - IgG_F2_FC_PD) ,levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)

limma_results <- topTable(fit = final_model,
                          coef = c(1, 2,3,4,5,6,7,8), 
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

summary(decideTests(final_model, method="global"))
limma_diffs <- limma_results

##### replace intensity values with keys for present/absent in PD and controls. First look at Log FC
# Controls: if log FC < 0.5, it's absent, give value = -2
# If log FC > 0.5, it's present, give value = 1
for (i in colnames(limma_diffs[,2:5])) { #CON columns
  limma_diffs[[i]] = case_when(
    limma_diffs[[i]]   < 0.5  ~ -2,
    limma_diffs[[i]]   > 0.5  ~ 1,
    .default = limma_diffs[[i]] 
  )
}

# PD: if log FC < 0.5, it's absent, give value = -3
# If log FC > 0.5, it's present, give value = 2
for (i in colnames(limma_diffs[,6:9])) { #PD columns
  limma_diffs[[i]] = case_when(
    limma_diffs[[i]]   < 0.5  ~ -3,
    limma_diffs[[i]]   > 0.5  ~ 2,
    .default = limma_diffs[[i]] 
  )
}

#now filter out the non-significant differences
limma_diffs_sig <- subset(limma_diffs, adj.P.Val < 0.05) #1212 proteins were true BPs

#make codes for whether proteins were uniquely present in one condition or another
# 3 = present in both
# -5 = absent in both
# -2 = present in CON, absent in PD
# 0 = present in PD, absent in CON
limma_diffs_sig$F1_FC <- limma_diffs_sig$F1_FC_PD + limma_diffs_sig$F1_FC_CON
limma_diffs_sig$F2_FC <- limma_diffs_sig$F2_FC_PD + limma_diffs_sig$F2_FC_CON
limma_diffs_sig$F1_PUT <- limma_diffs_sig$F1_PUT_PD + limma_diffs_sig$F1_PUT_CON
limma_diffs_sig$F2_PUT <- limma_diffs_sig$F2_PUT_PD + limma_diffs_sig$F2_PUT_CON

#summarise results with a bar plot
PDCON_pres <- colSums(dplyr::select(limma_diffs_sig, c("F1_FC",  "F2_FC", "F1_PUT", "F2_PUT")) == 3) %>% as.data.frame()
names(PDCON_pres) <- "PDCON_pres"

PDCON_abs <- colSums(dplyr::select(limma_diffs_sig, c("F1_FC",  "F2_FC", "F1_PUT", "F2_PUT")) == -5) %>% as.data.frame()
names(PDCON_abs) <- "PDCON_abs"

unique <- subset(limma_diffs_sig, F1_FC != 3 & F2_FC != 3 & F1_PUT != 3 & F2_PUT != 3)
PD <- subset(unique, F1_FC == 0 | F2_FC == 0 | F1_PUT == 0 | F2_PUT == 0)
PD <- subset(PD, F1_FC != -2 & F2_FC != -2  & F1_PUT != -2  & F2_PUT != -2  )
PDonly <- colSums(dplyr::select(PD, c("F1_FC",  "F2_FC", "F1_PUT", "F2_PUT")) == 0) %>% as.data.frame()
names(PDonly) <- "PDonly"

CON <- subset(unique, F1_FC == -2 | F2_FC == -2 | F1_PUT == -2 | F2_PUT == -2)
CON <- subset(CON, F1_FC != 0 & F2_FC != 0  & F1_PUT != 0  & F2_PUT != 0)
CONonly <- colSums(dplyr::select(CON, c("F1_FC",  "F2_FC", "F1_PUT", "F2_PUT")) == -2) %>% as.data.frame()
names(CONonly) <- "CONonly"

summary<- t(cbind(PDCON_pres, PDonly, CONonly)) %>% as.data.frame()
summary$status <- row.names(summary)
longsum <- pivot_longer(summary, cols = c("F1_FC", "F2_FC", "F1_PUT", "F2_PUT"), names_to = "group", values_to = "count")
longsum$status <- factor(longsum$status, levels = c("PDCON_pres", "CONonly", "PDonly"))
longsum$group <- factor(longsum$group, levels = c("F1_FC", "F1_PUT", "F2_FC", "F2_PUT"))

#write.csv(longsum, "D:/Newgraph/BPsummary_sum.csv") #export results for supplementary

barplot <- ggplot(longsum, aes(x=group, y=count,  group = status, fill = status, color = status)) + 
  theme_flo() + 
  ylab("Number of proteins") +
  xlab('') +
  scale_y_continuous(expand = c(0,0), limits=c(0, 500), breaks = seq(0, 500, by =100)) +
  scale_fill_manual(values = c("#276ab3", "black", "#d40000ff")) +
  scale_color_manual(values = c("#276ab3", "black", "#d40000ff")) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7, alpha = 0.6) + 
  theme(legend.position = c(0.8,0.9)) +
  guides(colour= "none") 

ggsave(barplot, filename = file.path('D:/Newgraph', 'barplot_diffs.pdf'), units = 'mm', height = 40, width =74, dpi = 300)

## heatmap for only uniquely altered proteins
PD_heatmap <- dplyr::select(PD, c("Protein","F1_FC",  "F2_FC", "F1_PUT", "F2_PUT")) #73 proteins
CON_heatmap <- dplyr::select(CON, c("Protein","F1_FC",  "F2_FC", "F1_PUT", "F2_PUT")) #151 proteins

summary<- rbind(PD_heatmap, CON_heatmap) 
m <- data.matrix(dplyr::select(summary, c("F1_FC", "F2_FC", "F1_PUT", "F2_PUT")))
rownames(m)<-summary$Protein

heatmap_diffsunique <- pheatmap(m, cluster_rows = T, cluster_cols = T, show_rownames = F, border_color = NA,
                          fontsize_col = 10, treeheight_col = 7, treeheight_row = 20,
                          angle_col = 0, color = c("#ecececff", "black", "#d40000ff")) #, #ffe6aaff, color = myColor,  breaks=myBreaks

ggsave(heatmap_diffsunique, filename = file.path('D:/Newgraph', 'heatmap_diffs_unique_nocentre.pdf'), units = 'mm', height = 80, width =90, dpi = 300)

### GO terms associated with differentially present proteins. 
#Using 'unique', which excludes proteins that are present in PD and CON in any fraction/region combo
PD_ids <- dplyr::select(subset(PD), c("Protein")) 
CON_ids <- dplyr::select(subset(CON), c("Protein")) 

all_proteins <- cc_qf_IP[["log_proteins_BP"]]
protein_info <- all_proteins %>% 
  rowData() %>%
  as_tibble %>%
  dplyr::select(ID = Genes,
         Protein_numb = Protein.Group,
         Protein_description = First.Protein.Description)

protein_info <- protein_info[!duplicated(protein_info$ID), ] #checks for duplicated proteins

diff_PD_des <- PD_ids %>%
  left_join(protein_info, c("Protein" = "ID"))

diff_CON_des <- CON_ids %>%
  left_join(protein_info, c("Protein" = "ID"))

# no sig enriched terms in the 73 proteins across both fractions in PD
ego_PD<- enrichGO(gene = diff_PD_des$Protein_numb,               # list of down proteins
                         universe = protein_info$Protein_numb,      # all proteins 
                         OrgDb = org.Hs.eg.db,                  # database to query
                         keyType = "UNIPROT",                   # protein ID encoding 
                         qvalueCutoff = 0.05,
                         ont = "BP",                            # can be CC, MF, BP, or ALL
                         readable = TRUE)

ego_CON <- enrichGO(gene = diff_CON_des$Protein_numb,               # list of down proteins
                  universe = protein_info$Protein_numb,      # all proteins 
                  OrgDb = org.Hs.eg.db,                  # database to query
                  keyType = "UNIPROT",                   # protein ID encoding 
                  qvalueCutoff = 0.05,
                  ont = "BP",                            # can be CC, MF, BP, or ALL
                  readable = TRUE)

simplify(ego_CON)

cnet_CON <- 
  cnetplot(simplify(ego_CON), 
           color_category = "#808080ff",
           color_gene = "black",
           showCategory = 8,
           #node_label = "category",
           #color_gene = "red",
           alpha_gene = 0.4,
           cex_label_category = 0.6,
           cex_label_gene = 0.4 ,
           shadowtext = "none",
           layout = "gem") +
  theme_flo() +
  labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  guides(size = "none") 
cnet_CON

ggsave(cnet_CON, filename = file.path('D:/Newgraph', 'cnet_allFsRs_rerun.pdf'), units = 'mm', height = 60, width =120, dpi = 300)

ego_CON_df <- ego_CON %>% as_tibble()
ego_CON_df <- ego_CON_df[order(ego_CON_df$Count, decreasing = TRUE),]  
ego_CON_df$Description <- factor(ego_CON_df$Description , levels=unique(ego_CON_df$Description))
ego_CON_df$GeneRatio_num <- as.numeric(sapply(strsplit(ego_CON_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(ego_CON_df$GeneRatio, "/"), "[", 2))
#write.csv(ego_CON_df, "D:/Newgraph/CONF1_GO.csv") #export results

############## Controls and PD: GO terms only in Fraction 1
PDonlydf <- subset(unique, F1_PUT == 0 | F1_FC == 0)
PDonlydf <- subset(PDonlydf, F2_FC != 0 & F2_PUT != 0 )
PD <- subset(PDonlydf, F1_FC != -2 & F2_FC != -2  & F1_PUT != -2  & F2_PUT != -2  )
PD_ids <- select(subset(PD), c("Protein")) 

CONonlydf <- subset(unique, F1_FC == -2 | F1_PUT == -2)
CONonlydf <- subset(CONonlydf, F2_FC != -2 &  F2_PUT != -2)
CON <- subset(CONonlydf, F1_FC != 0 & F2_FC != 0  & F1_PUT != 0  & F2_PUT != 0 )
CON_ids <- select(subset(CON), c("Protein")) 

all_proteins <- cc_qf_IP[["log_proteins_BP"]]
protein_info <- all_proteins %>% 
  rowData() %>%
  as_tibble %>%
  select(ID = Genes,
         Protein_numb = Protein.Group,
         Protein_description = First.Protein.Description)

protein_info <- protein_info[!duplicated(protein_info$ID), ] #checks for duplicated proteins

diff_PD_des <- PD_ids %>%
  left_join(protein_info, c("Protein" = "ID"))

diff_CON_des <- CON_ids %>%
  left_join(protein_info, c("Protein" = "ID"))

ego_PD<- enrichGO(gene = diff_PD_des$Protein_numb,               # list of down proteins
                  universe = protein_info$Protein_numb,      # all proteins 
                  OrgDb = org.Hs.eg.db,                  # database to query
                  keyType = "UNIPROT",                   # protein ID encoding 
                  qvalueCutoff = 0.05,
                  ont = "BP",                            # can be CC, MF, BP, or ALL
                  readable = TRUE)

ego_CON <- enrichGO(gene = diff_CON_des$Protein_numb,               # list of down proteins
                    universe = protein_info$Protein_numb,      # all proteins 
                    OrgDb = org.Hs.eg.db,                  # database to query
                    keyType = "UNIPROT",                   # protein ID encoding 
                    qvalueCutoff = 0.05,
                    ont = "BP",                            # can be CC, MF, BP, or ALL
                    readable = TRUE)


egoPD_df <- simplify(ego_PD) %>% as.data.frame()
egoCON_df <- simplify(ego_CON) %>% as.data.frame()

egoPD_df$Sample_Type <- "PD"
egoCON_df$Sample_Type <- "CON"

egoPD_df <- egoPD_df[order(egoPD_df$Count, decreasing = TRUE),]  
egoPD_df$Description <- factor(egoPD_df$Description , levels=unique(egoPD_df$Description))
egoPD_df$GeneRatio_num <- as.numeric(sapply(strsplit(egoPD_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(egoPD_df$GeneRatio, "/"), "[", 2))

egoCON_df <- egoCON_df[order(egoCON_df$Count, decreasing = TRUE),]  
egoCON_df$Description <- factor(egoCON_df$Description , levels=unique(egoCON_df$Description))
egoCON_df$GeneRatio_num <- as.numeric(sapply(strsplit(egoCON_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(egoCON_df$GeneRatio, "/"), "[", 2))

ego_BP_CONPD<- rbind(egoPD_df, egoCON_df)
#write.csv(ego_BP_CONPD, "D:/Newgraph/CONF1_PD_GO.csv") #export results

ego_BP_doptplot <- ggplot(ego_BP_CONPD, aes(x = Count, y = fct_rev(Description), alpha = p.adjust, color = Sample_Type))+
  ylab("GO Biological Process") +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8)) +
  geom_point(aes(size = GeneRatio_num),  shape = 16) +
  scale_alpha(range = c(1, 0.4)) +
  scale_y_discrete(labels = function(y) str_wrap(y, width=35)) +
  scale_color_manual(values = c("black", "red")) +
  #scale_color_gradient(guide = guide_colorbar(reverse = TRUE, title = "", order = 1), low = "#53acf0ff", high = "#010048") +
  guides(size = guide_legend(title = 'Gene Ratio', order = 2), fill = "none") +
  theme(legend.position = "right") +
  theme(plot.margin = margin(0.3,0.1,0.1,0.1, "cm")) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.minor = element_line(color = "grey92", size = 0.4, linetype = 2), panel.grid.major= element_line(color = "grey92", size = 0.4, linetype = 2)) 

ggsave(ego_BP_doptplot, filename = file.path('D:/Newgraph', 'CC_dotplotbinding.pdf'), units = 'mm', height = 105, width =105, dpi = 300)




