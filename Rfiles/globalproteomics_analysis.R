
rm(list = ls())
library(ggplot2)
library(tidyverse)
library(readxl)
library(ggpubr)
library(plyr)
library(readr)
library(pheatmap)
library("QFeatures")
library("limma")
library("here")
library("QFeatures")
library("NormalyzerDE")
library("limma")
library("factoextra")
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")
library("patchwork")
library(MsCoreUtils)
library(ggrepel)
library(ggh4x)
#library(missForest)


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

#Import raw data
df <- read.csv("C:/Users/Florence/OneDrive - University of Cambridge/proteomicspaper/ResultsFigures/Globalproteomics_rawdata.csv", sep = ",")
df <- df[!duplicated(df$Genes), ] #remove duplicated proteins

#add the dataframe to a QFeatures object, only the quantitative columns- each of the sample columns
cc_qf <- readQFeatures(assayData = df,
                       quantCols = 5:44, 
                       name = "proteins_raw")

#look at the rowdata to the QFeatures object: this is the non-quantitative data ie the gene descriptions. These are also called features
cc_qf[["proteins_raw"]] %>%
  rowData() %>%
  names()

#look at the column data: the names of the quantitative columns
cc_qf[["proteins_raw"]] %>%
  colData() %>% 
  rownames()

#Import metadata
metadata_df <- read_excel("C:/Users/Florence/OneDrive - University of Cambridge/proteomicspaper/ResultsFigures/global_meta_data.xlsx")

# Annotate colData with condition information
cc_qf$Case <- metadata_df$Case
cc_qf$Fraction <- metadata_df$Fraction
cc_qf$Sample_Type <- metadata_df$Sample_Type
colData(cc_qf[["proteins_raw"]]) <- colData(cc_qf)

######## Filtering #######
## Extract a copy of the raw protein data
raw_data_copy <- cc_qf[["proteins_raw"]] 

## Re-add the assay to our QFeatures object with a new name
cc_qf <- addAssay(x = cc_qf, 
                  y = raw_data_copy, 
                  name = "proteins_filtered")

#Make a dataframe with the raw numbers and a column for the Genes
df_proteins <- 
  cc_qf[["proteins_filtered"]] %>% 
  assay() %>%
  as.data.frame() %>%
  mutate(Genes = rowData(cc_qf[["proteins_filtered"]])$Genes) 

dim(df_proteins) 

#see how many missing values there are before filtering out
mv_filtered <- nNA(cc_qf, i = "proteins_filtered")
mv_filtered_df <- mv_filtered$nNAcols %>% as_tibble() 
#mv_filtered_df <- mv_filtered_df  %>%mutate(Fraction = colData(cc_qf)$Fraction) 
mv_filtered_df <- mv_filtered_df %>% separate(name, c("Case", "Fraction"), sep = "_")
mv_filtered_df$Fraction <- factor(mv_filtered_df$Fraction, levels = c("F1", "F2"))
mv_filtered_df$Sample_Type <- ifelse(grepl("PD",mv_filtered_df$Case),"PD","CON")


missing <- 
  ggplot(mv_filtered_df, aes(y = Case, x = pNA,  fill = Sample_Type)) +
  geom_bar(stat = "identity",position = position_dodge(0.9), width = 0.7, alpha = 0.6) +
  facet_wrap2( ~ Fraction, ncol = 1, nrow = 2, axes = "all", strip.position = "top") + 
  labs(y = "Sample", x = "Proportion missing values") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by =0.2)) +
  scale_fill_manual(values = c("black", "#d40000ff")) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 6))

ggsave(missing, filename = file.path('D:/Newgraph', 'missing_prefilter_Strap.pdf'), units = 'mm', height = 100, width =90, dpi = 300)


#Filter out proteins that are missing in > 5/10 samples per group
dfa <- df_proteins %>% filter(rowSums(!is.na(across(contains("CON") & contains("F1")))) >= 5) 
dfb <- df_proteins %>% filter(rowSums(!is.na(across(contains("PD") & contains("F1")))) >= 5)
dfc <- df_proteins %>% filter(rowSums(!is.na(across(contains("CON") & contains("F2")))) >= 5) 
dfd <- df_proteins %>% filter(rowSums(!is.na(across(contains("PD") & contains("F2")))) >= 5)
df_proteins_noNA <- rbind(dfa,dfb,dfc, dfd) 

df_proteins_noNA_nodups <- df_proteins_noNA[!duplicated(df_proteins_noNA$Genes), ] #removes duplicated proteins

#additional filtering for contaminants
df_proteins_noNA_nodups_noKR <- df_proteins_noNA_nodups %>% filter(!grepl("KRT", Genes)) 
df_proteins_noNA_nodups_noNA <- df_proteins_noNA_nodups_noKR %>% filter(rowSums(!is.na(across(contains("Genes")))) >= 1) 
df_proteins_MAPT <- subset(df_proteins_noNA_nodups_noNA, Genes == "MAPT;cRAP-MAPT") #keep MAPT
df_proteins_noNA_nodups_nocRAP <- df_proteins_noNA_nodups_noNA %>% filter(!grepl("cRAP", Genes)) 
df_proteins_noNA_nodups_nocRAP_MAPT <- rbind(df_proteins_noNA_nodups_nocRAP, df_proteins_MAPT)

dim(df_proteins_noNA_nodups_nocRAP_MAPT) #filtered df
#df_proteins_noNA_nodups_nocRAP_MAPT$Genes <- gsub("MAPT;cRAP-MAPT", "MAPT", df_proteins_noNA_nodups_nocRAP_MAPT$Genes)

#Join this to rowData
df_proteins_noNA_nodups_nocRAP_MAPT['marker'] = 1 #Make a column to show which genes are present in unfiltered rowData

#Put the rowdata in a df- it includes the Gene names
protein_rowdata <- 
  cc_qf[["proteins_filtered"]] %>% 
  rowData() %>% 
  as.data.frame()

#combine the rowdata with the filtered protein data frame, the rows will match up based on the Gene names and all of the proteins
mergedrowData <- left_join(x = protein_rowdata,
                           y = df_proteins_noNA_nodups_nocRAP_MAPT,
                           by = c("Genes" = "Genes"))

## Now add the Marker column to the proteins_filtered assay and filter out all of the proteins that don't have a 1 Marker
rowData(cc_qf[["proteins_filtered"]])$Missing <- mergedrowData$marker

cc_qf <- cc_qf %>%
  filterFeatures(~ Missing == 1, ##keep features with this feature 
                 i = "proteins_filtered", na.rm = TRUE)

# check missing values after filtering
mv_filtered <- nNA(cc_qf, i = "proteins_filtered")
mv_filtered_df <- mv_filtered$nNAcols %>% as_tibble() 
#mv_filtered_df <- mv_filtered_df  %>%mutate(Fraction = colData(cc_qf)$Fraction) 
mv_filtered_df <- mv_filtered_df %>% separate(name, c("Case", "Fraction"), sep = "_")
mv_filtered_df$Fraction <- factor(mv_filtered_df$Fraction, levels = c("F1", "F2"))
mv_filtered_df$Sample_Type <- ifelse(grepl("PD",mv_filtered_df$Case),"PD","CON")


missing <- 
  ggplot(mv_filtered_df, aes(y = Case, x = pNA,  fill = Sample_Type)) +
  geom_bar(stat = "identity",position = position_dodge(0.9), width = 0.7, alpha = 0.6) +
  facet_wrap2( ~ Fraction, ncol = 1, nrow = 2, axes = "all", strip.position = "top") + 
  labs(y = "Sample", x = "Proportion missing values") + 
  scale_x_continuous(expand = c(0,0), limits = c(0, 1), breaks = seq(0, 1, by =0.2)) +
  scale_fill_manual(values = c("black", "#d40000ff")) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 6))

ggsave(missing, filename = file.path('D:/Newgraph', 'missing_postfilter_Strap.pdf'), units = 'mm', height = 100, width =90, dpi = 300)


####### imputation, log transformation, and center-median normalisation ##########
cc_qf <- impute(object = cc_qf,
                method = "MinDet",
                i = "proteins_filtered",
                MARGIN = 2,
                name = "proteins_imputed")

cc_qf[["proteins_imputed"]] %>% assay() %>% as_tibble()

cc_qf[["proteins_imputed"]] %>%
  assay() %>%
  longFormat() %>%
  ggplot(aes(x = log2(value))) +
  geom_histogram() + 
  theme_bw() +
  xlab("(Log2) Abundance")

cc_qf <- logTransform(object = cc_qf, 
                      base = 2, 
                      i = "proteins_imputed", 
                      name = "log_proteins")

cc_qf <- normalize(cc_qf, 
                   i = "log_proteins", 
                   name = "log_proteins_centered",
                   method = "center.median")

df<- cc_qf[["log_proteins_centered"]] %>% assay() %>% as_tibble()
df<- cbind(rowData(cc_qf[["log_proteins_centered"]]), df)
#write.csv(df, "D:/Newgraph/filteredglobalpros.csv") #export data

#plot LFQ per sample, pre and post normalising

prenorm_df <- cc_qf[["log_proteins"]] %>%
  assay() %>%
  longFormat()

#prenorm_df$Fraction <- ifelse(grepl("F1", prenorm_df$colname), "F1", "F2")
prenorm_df <- prenorm_df %>% separate(colname, c("Case", "Fraction"), sep = "_")
prenorm_df$Sample_Type <- ifelse(grepl("PD",prenorm_df$Case),"PD","CON")


pre_norm <- prenorm_df %>%
  ggplot(aes(y = Case, x = value, color = Sample_Type)) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  geom_boxplot() +
  facet_wrap2( ~ Fraction, ncol = 1, nrow = 2, axes = "all", strip.position = "top") + 
  labs(y = "Sample", x = "log2 protein abundance")+
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 6))
pre_norm

ggsave(pre_norm, filename = file.path('D:/Newgraph', 'prenorm_LFQ.pdf'), units = 'mm', height = 100, width =89, dpi = 300)

post_norm_df <- cc_qf[["log_proteins_centered"]] %>%
  assay() %>%
  longFormat()

post_norm_df <- post_norm_df %>% separate(colname, c("Case", "Fraction"), sep = "_")
post_norm_df$Sample_Type <- ifelse(grepl("PD",post_norm_df$Case),"PD","CON")


post_norm <- post_norm_df %>%
  ggplot(aes(y = Case, x = value, color = Sample_Type)) +
  scale_color_manual(values = c("black", "#d40000ff")) +
  geom_boxplot() +
  facet_wrap2( ~ Fraction, ncol = 1, nrow = 2, axes = "all", strip.position = "top") + 
  labs(y = "Sample", x = "log2 protein abundance")+
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 6))
post_norm

ggsave(post_norm, filename = file.path('D:/Newgraph', 'post_norm_LFQ.pdf'), units = 'mm', height = 100, width =89, dpi = 300)

protein_pca <- cc_qf[["log_proteins_centered"]] %>% 
  assay() %>%
  t() %>%
  prcomp(scale = TRUE, center = FALSE)

summary(protein_pca)

S_Trap_PCA <- protein_pca$x %>%
  as_tibble() %>%
  mutate(Fraction= cc_qf[["log_proteins_centered"]]$Fraction) %>% 
  mutate(Case = cc_qf[["log_proteins_centered"]]$Case) %>%
  mutate(Sample_Type = cc_qf[["log_proteins_centered"]]$Sample_Type) %>%
  ggplot(aes(x = PC1, y = PC2, colour = factor(Sample_Type), shape = Fraction)) +
  geom_point(size = 1.5, alpha = 0.6) + 
  #scale_y_continuous(limits = c(-70, 82))+
  scale_color_manual(values = c("black", "#d40000ff")) +
  scale_shape_manual(values = c(16,1), labels= c('Fraction 1', 'Fraction 2')) +
  ylab("PC2 (11.6%)") +
  xlab('PC1 (67.9%)') +
  guides(colour = guide_legend(title = ''), shape = guide_legend(title = 'Fraction')) +
  theme_flo() +
  theme(legend.position = c(0.9, 0.85), legend.box = "vertical" ) +
  theme(panel.grid.minor = element_line(color = "grey92", size = 0.4, linetype = 2), panel.grid.major= element_line(color = "grey92", size = 0.4, linetype = 2)) 

ggsave(S_Trap_PCA, filename = file.path('D:/Newgraph', 'Strap_PCA.pdf'), units = 'mm', height = 60, width = 90, dpi = 300)

######## Limma analysis: Membrane fraction (F2) vs Cytosolic fraction (F1) #######
all_proteins <- cc_qf[["log_proteins_centered"]]
row_names <-  all_proteins %>% rowData() %>% as_tibble() %>% pull(Genes) %>% as.list()
all_proteins <- "rownames<-"(all_proteins, row_names)
all_proteins$Sample_Type <- factor(all_proteins$Sample_Type , levels = c("CON", "PD"))

Case <- all_proteins$Case
Sample_Type <- all_proteins$Sample_Type 
Fraction <- all_proteins$Fraction

#calculate duplicate correlation for case
Fraction <- factor(all_proteins$Fraction)
design <- model.matrix(~0+Fraction)
colnames(design) <- levels(Fraction)
dupcor <- duplicateCorrelation(assay(all_proteins),design,block=all_proteins$Case)
dupcor$consensus.correlation 

#now do pairwise comparisons between the fractions, controlling for correlations between the donors 
fit <- lmFit(assay(all_proteins), design, block=all_proteins$Case, correlation=dupcor$consensus.correlation)
contrasts <- makeContrasts(F2-F1, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)
summary(decideTests(final_model, method="global"))

## Format results
limma_results <- topTable(fit = final_model,
                          coef = 1,
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

## Verify
head(limma_results)

## Add direction and significance information
limma_results <- limma_results %>%
  mutate(direction = ifelse(logFC > 0.5, "Membrane", ifelse(logFC < -0.5, "Cytosol", "No change")),
         significance = ifelse(adj.P.Val < 0.05 & direction != "No Change", "sig", "not.sig"),
         fillcol = ifelse(direction == "Membrane" & significance == "sig", "sig",
                          ifelse(direction == "Membrane" & significance == "not.sig", "No Change",
                                 ifelse(direction == "Cytosol" & significance == "sig", "Cytosol",
                                        "No Change"))))

limma_results$direction = ifelse(limma_results$logFC < -0.5, "Cytosol", limma_results$direction)

enrichedinF2 <- subset(limma_results, logFC > 0.5) %>% select(Protein, logFC, adj.P.Val) %>% as_tibble()
enrichedinF1 <- subset(limma_results, logFC < -0.5) %>% select(Protein, logFC, adj.P.Val) %>% as_tibble()
limma_results$fillcol <- factor(limma_results$fillcol, c("No Change", "Down", "Up"))

### calculate the mean expression for each fraction and merge with Limma results
all_proteins_df <- cc_qf[["log_proteins_centered"]] %>% assay() %>% as_tibble()

# Select columns containing "weight"
cols <- grepl("F1", names(all_proteins_df))
df_F1<- all_proteins_df[, cols]
df_F1$F1_mean <- rowMeans(df_F1[,1:20])

cols <- grepl("F2", names(all_proteins_df))
df_F2<- all_proteins_df[, cols]
df_F2$F2_mean <- rowMeans(df_F2[,1:20])

genes <-  all_proteins %>% rowData() %>% as_tibble() 
genes <- dplyr::select(genes, c("Genes"))
mean_intensity <- cbind(genes, df_F1$F1_mean,df_F2$F2_mean )
merged <- left_join(mean_intensity, limma_results, by = c("Genes" = "Protein"))

#write.csv(merged, "D:/Newgraph/GlobalProt_limma_results.csv") #export results

#for Cytosol (F1) versus Membrane (F2)
F1_F2_volcano <- limma_results %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = fillcol)) +
  ylab("")+ #-log10 adjusted p-value
  xlab('') + #log2 fold change
  scale_x_continuous(limits = c(-6.5, 12), breaks = seq(-6, 12, by =3))+
  scale_y_continuous(expand = c(0,0), limits=c(0, 35), breaks = seq(0, 35, by =5)) +
  scale_alpha_manual(values = c(0.4,1)) +
  geom_point(data = subset(limma_results, adj.P.Val >= 0.05 & logFC > -0.5 | logFC < 0.5 ), size = 1.5, alpha = 0.4, color = "grey", shape = 16) +
  geom_hline(yintercept = 1.301, color = "black", linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "black", , linetype = "dashed", alpha = 0.4) +
  geom_point(data = subset(enrichedinF1, adj.P.Val <= 0.05), size = 1.5, color = "#276ab3", alpha = 0.4, shape = 16) +
  geom_point(data = subset(enrichedinF2, adj.P.Val <= 0.05), size = 1.5, color = "#276ab3", alpha = 0.4, shape = 1) +
  guides(color= guide_legend(title = '', override.aes = aes(label = "")), alpha = "none") +
  theme_flo()+
  theme(legend.box = ("vertical"), legend.position = c(0.88,0.92)) 

ggsave(F1_F2_volcano, filename = file.path('D:/Newgraph', 'STRAP_F1F2volcano.pdf'), units = 'mm', height = 60, width = 90, dpi = 300)

######################## GO analysis for each fraction ###############################
protein_info <- all_proteins %>% 
  rowData() %>%
  as_tibble %>%
  select(ID = Genes,
         Protein_numb = Protein.Group,
         Protein_description = First.Protein.Description)

#left join this to the significant proteins for their GO enrichment
limma_results_des <- limma_results %>%
  left_join(protein_info, by = c("Protein" = "ID"))

limma_results_des_changing <- limma_results_des %>% 
  as_tibble() %>%
  filter(significance =="sig")

sig_up <- limma_results_des_changing %>%
  filter(logFC  > 0.5 )

sig_down <- limma_results_des_changing %>%
  filter(logFC < -0.5)

ego_up_CC <- enrichGO(gene = sig_up$Protein_numb,               # list of up proteins: membrane fraction
                      universe = limma_results_des$Protein_numb,      # all proteins 
                      OrgDb = org.Hs.eg.db,                  # database to query
                      keyType = "UNIPROT",                   # protein ID encoding 
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      ont = "CC",                            #Cellular component
                      readable = TRUE)

ego_down_CC <- enrichGO(gene = sig_down$Protein_numb,               # list of down proteins: cytosol fraction
                        universe = limma_results_des$Protein_numb,      # all proteins 
                        OrgDb = org.Hs.eg.db,                  # database to query
                        keyType = "UNIPROT",                   # protein ID encoding 
                        qvalueCutoff = 0.05,
                        ont = "CC",                            #Cellular component
                        readable = TRUE)

#put results in DF to plot
ego_up_CC_df <- ego_up_CC %>% as_tibble()
ego_up_CC_df <- ego_up_CC_df[order(ego_up_CC_df$Count, decreasing = TRUE),]  
ego_up_CC_df$Description <- factor(ego_up_CC_df$Description , levels=unique(ego_up_CC_df$Description))
ego_up_CC_df$GeneRatio_num <- as.numeric(sapply(strsplit(ego_up_CC_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(ego_up_CC_df$GeneRatio, "/"), "[", 2))
ego_up_CC_df$Fraction <- "Membrane"

ego_down_CC_df <- ego_down_CC %>% as_tibble()
ego_down_CC_df <- ego_down_CC_df[order(ego_down_CC_df$Count, decreasing = TRUE),]  
ego_down_CC_df$Description <- factor(ego_down_CC_df$Description , levels=unique(ego_down_CC_df$Description))
ego_down_CC_df$GeneRatio_num <- as.numeric(sapply(strsplit(ego_down_CC_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(ego_down_CC_df$GeneRatio, "/"), "[", 2))
ego_down_CC_df$Fraction <- "Cytosolic"

merged <- rbind(ego_up_CC_df, ego_down_CC_df)
write.csv(merged, "D:/Newgraph/GlobalProt_Fraction_GO.csv") #export results

#top 5 in membrane fraction
Dotplot_CC_up <- ggplot(ego_up_CC_df[1:5,], aes(x = Count, y = fct_rev(Description), color = p.adjust))+
  ylab("GO Cellular component") +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8)) +
  geom_point(aes(size = GeneRatio_num), shape = 1, alpha = 0.8, stroke = 1.1) +
  scale_y_discrete(labels = function(y) str_wrap(y, width=25)) +
  scale_color_gradient(guide = guide_colorbar(reverse = TRUE, title = "", order = 1), low = "#53acf0ff", high = "#010048") +
  guides(size = guide_legend(title = 'Gene Ratio', order = 2)) +
  theme(legend.position = "right") +
  theme(plot.margin = margin(0.3,0.1,0.1,0.1, "cm")) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.minor = element_line(color = "grey92", size = 0.4, linetype = 2), panel.grid.major= element_line(color = "grey92", size = 0.4, linetype = 2)) 

ggsave(Dotplot_CC_up, filename = file.path('D:/Newgraph', 'STRAP_F2_CC5.pdf'), units = 'mm', height = 55, width =89, dpi = 300)

#Top 20 for supplementary
Dotplot_CC_up <- ggplot(ego_up_CC_df[1:20,], aes(x = Count, y = fct_rev(Description), color = p.adjust))+
  ylab("GO Cellular component") +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8)) +
  geom_point(aes(size = GeneRatio_num), shape = 1, alpha = 0.8, stroke = 1.1) +
  scale_y_discrete(labels = function(y) str_wrap(y, width=25)) +
  scale_color_gradient(guide = guide_colorbar(reverse = TRUE, title = "", order = 1), low = "#53acf0ff", high = "#010048") +
  guides(size = guide_legend(title = 'Gene Ratio', order = 2)) +
  theme(legend.position = "right") +
  theme(plot.margin = margin(0.3,0.1,0.1,0.1, "cm")) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.minor = element_line(color = "grey92", size = 0.4, linetype = 2), panel.grid.major= element_line(color = "grey92", size = 0.4, linetype = 2)) 

ggsave(Dotplot_CC_up, filename = file.path('D:/Newgraph', 'STRAP_F2_CC20.pdf'), units = 'mm', height = 150, width =89, dpi = 300)

#top 5 in cytosol for main text
Dotplot_CC_down <- ggplot(ego_down_CC_df[1:5,], aes(x = Count, y = fct_rev(Description), color = p.adjust))+
  ylab("GO Cellular component") +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8)) +
  geom_point(aes(size = GeneRatio_num), alpha = 0.8, shape = 16) +
  scale_y_discrete(labels = function(y) str_wrap(y, width=25)) +
  scale_color_gradient(guide = guide_colorbar(reverse = TRUE, title = "", order = 1), low = "#53acf0ff", high = "#010048") +
  guides(size = guide_legend(title = 'Gene Ratio', order = 2)) +
  theme(legend.position = "right") +
  theme(plot.margin = margin(0.3,0.1,0.1,0.1, "cm")) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.minor = element_line(color = "grey92", size = 0.4, linetype = 2), panel.grid.major= element_line(color = "grey92", size = 0.4, linetype = 2)) 

ggsave(Dotplot_CC_down, filename = file.path('D:/Newgraph', 'STRAP_F1_CC.pdf'), units = 'mm', height = 55, width =89, dpi = 300)

#top 20 for supplementary
Dotplot_CC_down <- ggplot(ego_down_CC_df[1:20,], aes(x = Count, y = fct_rev(Description), color = p.adjust))+
  ylab("GO Cellular component") +
  theme_flo() +
  theme(axis.text =element_text(color = 'black', face= 'plain', size = 8)) +
  geom_point(aes(size = GeneRatio_num), alpha = 0.8, shape = 16) +
  scale_y_discrete(labels = function(y) str_wrap(y, width=25)) +
  scale_color_gradient(guide = guide_colorbar(reverse = TRUE, title = "", order = 1), low = "#53acf0ff", high = "#010048") +
  guides(size = guide_legend(title = 'Gene Ratio', order = 2)) +
  theme(legend.position = "right") +
  theme(plot.margin = margin(0.3,0.1,0.1,0.1, "cm")) +
  coord_cartesian(clip = "off") +
  theme(panel.grid.minor = element_line(color = "grey92", size = 0.4, linetype = 2), panel.grid.major= element_line(color = "grey92", size = 0.4, linetype = 2)) 

ggsave(Dotplot_CC_down, filename = file.path('D:/Newgraph', 'STRAP_F1_CC20.pdf'), units = 'mm', height = 150, width =89, dpi = 300)

########### Limma with contrasts to look at PD vs control differences in each fraction ######
#this analysis has the same effect as following the steps in Section 9.5.2 and 9.5.3 of the Limma instruction manual on factorial designs
#which is achieved by using model.matrix with interactions only and not main effects

all_proteins <- cc_qf[["log_proteins_centered"]]
row_names <-  all_proteins %>% rowData() %>% as_tibble() %>% pull(Genes) %>% as.list()
all_proteins <- "rownames<-"(all_proteins, row_names)
all_proteins$Sample_Type <- factor(all_proteins$Sample_Type , levels = c("CON", "PD"))

Case <- all_proteins$Case
Sample_Type <- all_proteins$Sample_Type 
Fraction <- all_proteins$Fraction

Fraction <- factor(all_proteins$Fraction)
design <- model.matrix(~0 + Sample_Type:Fraction)
colnames(design) <- c("CON_F1", "PD_F1", "CON_F2", "PD_F2")
fit <- lmFit(assay(all_proteins), design)
contrasts <- makeContrasts(F1 = PD_F1-CON_F1, F2 = PD_F2-CON_F2, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)
summary(decideTests(final_model, method="global"))

## Limma results for each fraction
limma_results <- topTable(fit = final_model,
                          coef = c(1),
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 


## Add direction and significance information
limma_results <- limma_results %>%
  mutate(direction = ifelse(logFC > 0, "PD", "CON"),
         significance = ifelse(adj.P.Val < 0.05, "sig", "not.sig"),
         fillcol = ifelse(direction == "PD" & significance == "sig", "PD",
                          ifelse(direction == "PD" & significance == "not.sig", "No Change",
                                 ifelse(direction == "CON" & significance == "sig", "CON",
                                        "No Change"))))

limma_results$Fraction <- "Cytosolic"
write.csv(limma_results, "D:/Newgraph/GlobalProt_limma_results_PDCONF1.csv") #export results

enrichedinPD <- subset(limma_results, logFC > 0.5) %>% select(Protein, logFC, adj.P.Val) %>% as_tibble()
enrichedinCON <- subset(limma_results, logFC < -0.5) %>% select(Protein, logFC, adj.P.Val) %>% as_tibble()

#list of F1 enriched proteins here
PD <- subset(enrichedinPD, adj.P.Val <= 0.05)
con <- subset(enrichedinCON, adj.P.Val <= 0.05)
PDCON <- rbind(PD, con)
#write.csv(PDCON, "D:/Newgraph/PDCONenriched_F1.csv")

limma_results$fillcol <- factor(limma_results$fillcol, c("No Change", "CON", "PD"))

CON_PD_volcanoF1 <- limma_results %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = fillcol)) +
  ylab("")+
  xlab('') + #log2 fold change
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by =3)) + 
  scale_y_continuous(expand = c(0,0), limits=c(NA, 5)) +
  scale_alpha_manual(values = c(0.4,1)) +
  geom_point(size = 1.5, alpha = 0.6, color = "grey", shape = 16) + #change shape for F
  geom_hline(yintercept = 1.301, color = "black", linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "black", , linetype = "dashed", alpha = 0.4) +
  geom_point(data = subset(enrichedinCON, adj.P.Val <= 0.05), size = 1.5, color = "black", alpha = 0.8, shape = 16) +
  geom_text_repel(data = subset(enrichedinCON, adj.P.Val <= 0.05), aes(label = Protein), size = 2.5, color = "black") +
  geom_point(data = subset(enrichedinPD, adj.P.Val <= 0.05), size = 1.5, color = "#d40000ff", alpha = 0.4, shape = 16) +
  geom_text_repel(data = subset(enrichedinPD, adj.P.Val <= 0.05), aes(label = Protein), size = 2.5, color = "#d40000ff") +
  guides(color= guide_legend(title = '', override.aes = aes(label = "")), alpha = "none") +
  theme_flo()+
  theme(legend.box = ("vertical"), legend.position = c(0.88,0.92)) 

ggsave(CON_PD_volcanoF1, filename = file.path('D:/Newgraph', 'STRAP_F1_PDCONvolcano.pdf'), units = 'mm', height = 60, width = 90, dpi = 300)

#GO enrichment analysis for F1 enriched proteins
all_proteins_F1 <- all_proteins
protein_info_F1 <- all_proteins_F1 %>% 
  rowData() %>%
  as_tibble %>%
  select(Protein = Genes,
         Protein_numb = Protein.Group,
         Protein_description = First.Protein.Description)

#left join this to the significant proteins for their GO enrichment
limma_results_des_F1 <- limma_results %>%
  left_join(protein_info_F1, by = "Protein")

limma_results_des_changingF1 <- limma_results_des_F1 %>% 
  as_tibble() %>%
  filter(significance =="sig")

sig_upF1 <- limma_results_des_changingF1 %>%
  filter(logFC > 0.5)

sig_downF1 <- limma_results_des_changingF1 %>%
  filter(logFC < 0.5)

sig_downF1$Protein_numb <- gsub("P0DP23;P0DP24;P0DP25", "P0DP25", sig_downF1$Protein_numb)

#ego up = enriched in PD. No significant terms
ego_up_F1 <- enrichGO(gene = sig_upF1$Protein_numb,               # list of down proteins
                        universe = limma_results_des_F1$Protein_numb,      # all proteins 
                        OrgDb = org.Hs.eg.db,                  # database to query
                        keyType = "UNIPROT",                   # protein ID encoding 
                        qvalueCutoff = 0.05,
                        ont = "MF",                            # can be CC, MF, BP, or ALL
                        readable = TRUE)

#ego down = enriched in controls. Significant terms plotted with CNET below
ego_down_F1 <- enrichGO(gene = sig_downF1$Protein_numb,               # list of down proteins
                        universe = limma_results_des_F1$Protein_numb,      # all proteins 
                        OrgDb = org.Hs.eg.db,                  # database to query
                        keyType = "UNIPROT",                   # protein ID encoding 
                        qvalueCutoff = 0.05,
                        ont = "MF",                            # can be CC, MF, BP, or ALL
                        readable = TRUE)

cnet_F1 <- 
  cnetplot(ego_down_F1, 
           color_category = "#808080ff",
           color_gene = "black",
           alpha_gene = 0.4,
           cex_label_category = 0.5,
           cex_label_gene = 0.5 ,
           shadowtext = "none",
           layout = "gem") +
  theme_flo() +
  labs(x = "", y = "") +
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  guides(size = "none") 
cnet_F1

ggsave(cnet_F1, filename = file.path('D:/Newgraph', 'cnet_F1.pdf'), units = 'mm', height = 65, width =90, dpi = 300)

ego_down_df <- ego_down_F1 %>% as_tibble()
ego_down_df <- ego_down_df[order(ego_down_df$Count, decreasing = TRUE),]  
ego_down_df$Description <- factor(ego_down_df$Description , levels=unique(ego_down_df$Description))
ego_down_df$GeneRatio_num <- as.numeric(sapply(strsplit(ego_down_df$GeneRatio, "/"), "[", 1)) / as.numeric(sapply(strsplit(ego_down_df$GeneRatio, "/"), "[", 2))
ego_down_df$Fraction <- "Cytosol"
write.csv(ego_down_df, "D:/Newgraph/ego_down_df.csv") #export results

## F2
limma_results <- topTable(fit = final_model,
                          coef = 2,
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
rownames_to_column("Protein") 

## Verify
head(limma_results)

## Add direction and significance information
limma_results <- limma_results %>%
  mutate(direction = ifelse(logFC > 0, "PD", "CON"),
         significance = ifelse(adj.P.Val < 0.05, "sig", "not.sig"),
         fillcol = ifelse(direction == "PD" & significance == "sig", "PD",
                          ifelse(direction == "PD" & significance == "not.sig", "No Change",
                                 ifelse(direction == "CON" & significance == "sig", "CON",
                                        "No Change"))))

limma_results$Fraction <- "Membrane"
write.csv(limma_results, "D:/Newgraph/GlobalProt_limma_results_PDCONF2.csv") #export results

enrichedinPD <- subset(limma_results, logFC > 0.5) %>% select(Protein, logFC, adj.P.Val) %>% as_tibble()
enrichedinCON <- subset(limma_results, logFC < -0.5) %>% select(Protein, logFC, adj.P.Val) %>% as_tibble()
limma_results$fillcol <- factor(limma_results$fillcol, c("No Change", "CON", "PD"))

PD <- subset(enrichedinPD, adj.P.Val <= 0.05)
con <- subset(enrichedinCON, adj.P.Val <= 0.05)
PDCON_F2 <- rbind(PD, con)
#write.csv(PDCON_F2, "D:/Newgraph/PDCONenriched_F2.csv")

CON_PD_volcanoF2 <- limma_results %>%
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = fillcol)) +
  ylab("")+
  xlab('') + #log2 fold change
  scale_x_continuous(limits = c(-6, 6), breaks = seq(-6, 6, by =3)) + 
  scale_y_continuous(expand = c(0,0), limits=c(NA, 5)) +
  scale_alpha_manual(values = c(0.4,1)) +
  geom_point(size = 1.5, alpha = 0.6, color = "grey", shape = 1) + #change shape for F
  geom_hline(yintercept = 1.301, color = "black", linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "black", , linetype = "dashed", alpha = 0.4) +
  geom_point(data = subset(enrichedinCON, adj.P.Val <= 0.05), size = 1.5, color = "black", alpha = 0.8, shape = 1) +
  geom_text_repel(data = subset(enrichedinCON, adj.P.Val <= 0.05), aes(label = Protein), size = 2.5, color = "black") +
  geom_point(data = subset(enrichedinPD, adj.P.Val <= 0.05), size = 1.5, color = "#d40000ff", alpha = 0.4, shape = 1) +
  geom_text_repel(data = subset(enrichedinPD, adj.P.Val <= 0.05), max.overlaps = 20, aes(label = Protein), size = 2.5, color = "#d40000ff") +
  guides(color= guide_legend(title = '', override.aes = aes(label = "")), alpha = "none") +
  theme_flo()+
  theme(legend.box = ("vertical"), legend.position = c(0.88,0.92)) 

ggsave(CON_PD_volcanoF2, filename = file.path('D:/Newgraph', 'STRAP_F2_PDCONvolcano.pdf'), units = 'mm', height = 60, width = 90, dpi = 300)

#GO enrichment analysis
all_proteins_F2 <- all_proteins
protein_info_F2 <- all_proteins_F2 %>% 
  rowData() %>%
  as_tibble %>%
  select(Protein = Genes,
         Protein_numb = Protein.Group,
         Protein_description = First.Protein.Description)
PDCON$Fraction <- 1
PDCON_F2$Fraction <- 2
PDCON_tot <-rbind(PDCON, PDCON_F2)
pronames <- PDCON_tot %>%
  left_join(protein_info_F2, by = "Protein")
write.csv(pronames, "D:/Newgraph/pronames.csv")

#left join this to the significant proteins for their GO enrichment
limma_results_des_F2 <- limma_results %>%
  left_join(protein_info_F2, by = "Protein")

limma_results_des_changingF2 <- limma_results_des_F2 %>% 
  as_tibble() %>%
  filter(significance =="sig")

sig_upF2 <- limma_results_des_changingF2 %>%
  filter(logFC > 0.5)

sig_downF2 <- limma_results_des_changingF2 %>%
  filter(logFC < 0.5)

#ego down: enriched in Controls: no significantly enriched terms
ego_down_F2 <- enrichGO(gene = sig_downF2$Protein_numb,               # list of down proteins
                        universe = limma_results_des_F2$Protein_numb,      # all proteins 
                        OrgDb = org.Hs.eg.db,                  # database to query
                        keyType = "UNIPROT",                   # protein ID encoding 
                        qvalueCutoff = 0.05,
                        ont = "MF",                            # can be CC, MF, BP, or ALL
                        readable = TRUE)

#ego up: enriched in PD: no significantly enriched terms
ego_up_F2 <- enrichGO(gene = sig_upF2$Protein_numb,               # list of down proteins
                        universe = limma_results_des_F2$Protein_numb,      # all proteins 
                        OrgDb = org.Hs.eg.db,                  # database to query
                        keyType = "UNIPROT",                   # protein ID encoding 
                        qvalueCutoff = 0.05,
                        ont = "MF",                            # can be CC, MF, BP, or ALL
                        readable = TRUE)

####### Now compare F2- F1 difference in PD vs CON
Fraction <- factor(all_proteins$Fraction)
design <- model.matrix(~0+Fraction)
colnames(design) <- levels(Fraction)
dupcor <- duplicateCorrelation(assay(all_proteins), design,block = all_proteins$Case)
dupcor$consensus.correlation 

Fraction <- factor(all_proteins$Fraction)
design <- model.matrix(~0 + Sample_Type:Fraction)
colnames(design) <- c("CON_F1", "PD_F1", "CON_F2", "PD_F2")
fit <- lmFit(assay(all_proteins), design, block= all_proteins$Case, correlation=dupcor$consensus.correlation)
contrasts <- makeContrasts(Diff = ((PD_F2-PD_F1)) - ((CON_F2-CON_F1)), PD = (PD_F2-PD_F1), CON = (CON_F2-CON_F1), levels=design)
fit2 <- contrasts.fit(fit, contrasts)
final_model <- eBayes(fit2, trend=TRUE, robust = T)
summary(decideTests(final_model, method="global"))

##Difference between fractions for PD and controls: calculate logFC together for correlation plot
limma_results <- topTable(fit = final_model,
                          coef = c(2,3), # 1 for Diff, 2 for PD, 3 for CON
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

write.csv(limma_results, "D:/Newgraph/GlobalProt_logFC.csv") #export results

corr_F1F2 <- ggplot(limma_results, aes(x = PD, y = CON)) + 
  labs(x="PD Fraction 2 - Fraction 1", 
       y= 'CON Fraction 2 - Fraction 1') +
  stat_cor(method = "pearson", size = 3, label.x.npc = 0.05, show.legend = FALSE, p.accuracy = 0.0001, aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  geom_point(size = 1.2, alpha = 0.3, shape = 17, colour = "black") +
  geom_hline(yintercept = c(-0.5, 0.5), color = "black", linetype = "dashed", alpha = 0.4) +
  geom_vline(xintercept = c(-0.5, 0.5), color = "black", , linetype = "dashed", alpha = 0.4) +
  scale_y_continuous(expand = c(0, 0), limits=c(-8, 12), breaks = seq(-8, 12, by =4)) + 
  scale_x_continuous(expand = c(0, 0), limits=c(-8,12), breaks = seq(-8, 12, by =4)) + 
  geom_smooth(method='lm', formula= y~x, show.legend = FALSE, se = FALSE) +
  guides(fill = 'none', colour = guide_legend(title = '')) +
  theme_flo() +
  theme(legend.position = c(0.96, 0.85))

ggsave(corr_F1F2, filename = file.path('D:/Newgraph', 'STRAP_F1F2_corr.pdf'), units = 'mm', height = 52, width = 110, dpi = 300)

###### Heatmap: calculate the log FC for CON and PD individually to get appropriate p-value
limma_results <- topTable(fit = final_model,
                          coef = c(2,3), # 1 for Diff, 2 for PD, 3 for CON
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

##### specifically investigate only the differentially expressed proteins from the above limma analysis in each fraction
##Difference between fractions for PD and controls
corrdiffs <- subset(limma_results, (limma_results$Protein %in% PDCON$Protein) | (limma_results$Protein %in% PDCON_F2$Protein))
#write.csv(corrdiffs, "D:/Newgraph/diffpros_F1F2.csv") #export results

PDCON$Effect <- ifelse(PDCON$logFC < 0, "F1_CON", "F1_PD")
PDCON_F2$Effect <- ifelse(PDCON_F2$logFC < 0, "F2_CON", "F2_PD")
PDCON_tot <-rbind(PDCON_F2, PDCON)

diffs <- left_join(corrdiffs, select(PDCON_tot, c("Protein", "Effect")), by = "Protein")
diffs$sig <- ifelse(diffs$adj.P.Val < 0.05, "sig", "not")
diffs$sig <- ifelse(diffs$adj.P.Val < 0.05, "*", "not")
diffs$sig <- ifelse(diffs$adj.P.Val < 0.01, "**", diffs$sig)
diffs$sig <- ifelse(diffs$adj.P.Val < 0.001, "***", diffs$sig)
diffs$sig <- ifelse(diffs$adj.P.Val < 0.001, "****", diffs$sig)

diffs$Effect <-factor(diffs$Effect, levels = c("F1_PD", "F1_CON", "F2_PD","F2_CON"))
diffs = diffs[order(diffs$Effect), ]

#write.csv(diffs, "D:/Newgraph/GlobalProt_CONloFC.csv") #export results
annotation_row = data.frame(
  Effect = factor(diffs$Effect, c("F1_PD", "F1_CON", "F2_PD","F2_CON")),
  Sig = factor(diffs$sig, levels = c("sig", "not"))
)
rownames(annotation_row) = diffs$Proteins

ann_colors = list(
  Effect = c(F1_PD = "red", F1_CON = "black", F2_PD = "maroon",F2_CON = "grey"),
  Sig = c(sig = "#276ab3", not = "white")
)

matrix <- data.matrix(select(diffs, c("PD", "CON")))
rownames(matrix)<- diffs$Protein 

paletteLength <- 100
myColor <- colorRampPalette(c("navy", "#ecececff", "#d40000ff"))(paletteLength)
myBreaks <- c(seq(min(matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(matrix)/paletteLength, max(matrix), length.out=floor(paletteLength/2)))


heatmap_diffs <- pheatmap(matrix, cluster_rows = F, cluster_cols = F, show_rownames = T, 
                          fontsize_col = 10, fontsize_row = 8, treeheight_row = 15, treeheight_col = 0, 
                          angle_col = 0, annotation_row = annotation_row, 
                          annotation_colors = ann_colors, breaks=myBreaks) #, color = myColor

ggsave(heatmap_diffs, filename = file.path('D:/Newgraph', 'heatmap_diffs.pdf'), units = 'mm', height = 115, width =90, dpi = 300)

#need to redo limma with the difference as the coefficient to get the appropriate p-value
#for the difference between PD and CON
limma_results <- topTable(fit = final_model,
                          coef = c(1), # 1 for Diff, 2 for PD, 3 for CON
                          adjust.method = "BH",    # Method for multiple hypothesis testing
                          number = Inf) %>%        # Print results for all proteins
  rownames_to_column("Protein") 

corrdiffs <- subset(limma_results, (limma_results$Protein %in% PDCON$Protein) | (limma_results$Protein %in% PDCON_F2$Protein))
PDCON$Effect <- ifelse(PDCON$logFC < 0, "F1_CON", "F1_PD")
PDCON_F2$Effect <- ifelse(PDCON_F2$logFC < 0, "F2_CON", "F2_PD")
PDCON_tot <-rbind(PDCON_F2, PDCON)
diffs <- left_join(corrdiffs, select(PDCON_tot, c("Protein", "Effect")), by = "Protein")
diffs$sig <- ifelse(diffs$adj.P.Val < 0.05, "sig", "not")
diffs$sig <- ifelse(diffs$adj.P.Val < 0.05, "*", "not")
diffs$sig <- ifelse(diffs$adj.P.Val < 0.01, "**", diffs$sig)
diffs$sig <- ifelse(diffs$adj.P.Val < 0.001, "***", diffs$sig)
diffs$sig <- ifelse(diffs$adj.P.Val < 0.001, "****", diffs$sig)
diffs$Effect <-factor(diffs$Effect, levels = c("F1_PD", "F1_CON", "F2_PD","F2_CON"))
diffs = diffs[order(diffs$Effect), ]

write.csv(diffs, "D:/Newgraph/GlobalProt_logFCfractionsdiff.csv") #export results



