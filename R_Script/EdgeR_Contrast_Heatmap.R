
library("tximport")
library("limma")
library("locfit")
library("edgeR")
library("tidyverse")
library("dplyr")

#import files from salmon to R 1st step
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon

#creating files from the salmon_quant outputs file paths
files <- file.path(getwd(), "sturgeon_quant", list.files("./sturgeon_quant/"), "quant.sf")
names(files) <- paste0(list.files("./sturgeon_quant/"))

#actually reading in the files from that path
txi.salmon <- tximport(files, type = "salmon", txOut = TRUE)
head(txi.salmon$counts)

files

#gene models to estimates differential expression (all computation) impress ken!!!

cts <- txi.salmon$counts
head(cts)
normMat <- txi.salmon$length
head(normMat)
# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)

# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

#creating design matrix
metadata <- metadata %>%
  mutate(
    state = factor(state, levels = c("Fed", "Fast")),
    time = factor(time, levels = c("24hrs", "48hrs", "168hrs", "336hrs")),
    treatment = factor(
      treatment,
      levels = c(
        "Fed_24hrs", "Fast_24hrs",
        "Fed_48hrs", "Fast_48hrs",
        "Fed_168hrs", "Fast_168hrs",
        "Fed_336hrs", "Fast_336hrs")))

write_tsv(x = metadata, file = "metadata.txt")

# filtering using the design information
design <- model.matrix(~0 + treatment, data = metadata)
keep <- filterByExpr(y, design)
head(keep)
y <- y[keep, ]

#For creating a matrix of CPMs within edgeR, the following code chunk can be used: since we already removed zero across transcripts and samples 

cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)

#portion of creating for negative binomial distribution estimateDisp
y<-estimateDisp(y,design)
y

#Take a look at a PCA
MDS_full <- plotMDS(y, top = nrow(y$counts), labels = metadata$file_ID, dim.plot =c(2,3))

#where you get the variance from P1 and P2 for plot
MDS_full$var.explained

MDS_data_tibble <- tibble(x_coord = MDS_full$x, y_coord = MDS_full$y, Time = metadata$time, State = metadata$state)

MDS_data_tibble$Time <- factor(MDS_data_tibble$Time, levels = c("24hrs", "48hrs", "168hrs", "336hrs"))

PCA_plot12 <- ggplot(data = MDS_data_tibble, aes(x = x_coord, y = y_coord, color = Time, shape = State)) +
  geom_point(size = 4.5, alpha = 0.9) +
  
  xlab(bquote("Dimension 1 (6.9%)")) +
  ylab(bquote("Dimension 2 (4.9%)")) +
  #fishualize::scale_colour_fish_d(option = "Naso_lituratus", direction = 1) +
  scale_colour_manual(values = c("#d73027", "#f46d43", "#74add1","#4575b4" )) +
  scale_shape_discrete(labels = c("Fasted", "Fed")) +
  theme_bw(base_size = 20) +
  theme(legend.background = element_blank(), legend.key = element_blank(), legend.key.size = unit(0.4, 'cm'), legend.spacing.y = unit(0.1, 'cm'), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
PCA_plot12
ggsave(filename = "Dimensions_PCA12.pdf", plot = PCA_plot12, dpi = 3000, scale = 1.34)

PCA_plot <- ggplot(data = MDS_data_tibble, aes(x = x_coord, y = y_coord, color = Time, shape = State)) +
  geom_point(size = 4.5, alpha = 0.9) +
  
  xlab(bquote("Dimension 2 (4.9%)")) +
  ylab(bquote("Dimension 3 (3.3%)")) +
  #fishualize::scale_colour_fish_d(option = "Naso_lituratus", direction = 1) +
  scale_colour_manual(values = c("#d73027", "#f46d43", "#74add1","#4575b4" )) +
  scale_shape_discrete(labels = c("Fasted", "Fed")) +
  theme_bw(base_size = 20) +
  theme(legend.background = element_blank(), legend.key = element_blank(), legend.key.size = unit(0.4, 'cm'), legend.spacing.y = unit(0.1, 'cm'), axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())
PCA_plot


ggsave(filename = "Dimensions_PCA.pdf", plot = PCA_plot, dpi = 3000, scale = 1.34)

#installing fishualize package
install.packages("fishualize")

#quality control estimate dispersion CPM logs per millions by you dont want dispersion by CPM (want to be relatively flat)
plotBCV(y)

#fit the GLM using the edgeR data and the design
y_fit <- glmQLFit(y, design, robust = TRUE)
result <- glmQLFTest(y_fit, coef=2)
plotQLDisp(y_fit)

#comparisons for differntial gene expression
treatment_contrast <- makeContrasts(FedvsFast_24 = treatmentFast_24hrs - treatmentFed_24hrs,
                                    FedvsFast_48 = treatmentFast_48hrs - treatmentFed_48hrs,
                                    FedvsFast_168 = treatmentFast_168hrs - treatmentFed_168hrs,
                                    FedvsFast_336 = treatmentFast_336hrs - treatmentFed_336hrs,
                                    levels = design)
treatment_contrast
#pull results from contrasts
quasi_test_FedvsFast_24 <- glmQLFTest(y_fit, contrast = treatment_contrast[,"FedvsFast_24"])

quasi_test_FedvsFast_48 <- glmQLFTest(y_fit, contrast = treatment_contrast[,"FedvsFast_48"])

quasi_test_FedvsFast_168 <- glmQLFTest(y_fit, contrast = treatment_contrast[,"FedvsFast_168"])

quasi_test_FedvsFast_336 <- glmQLFTest(y_fit, contrast = treatment_contrast[,"FedvsFast_336"])

#FedvsFast_24
nrow(topTags(quasi_test_FedvsFast_24, p.value = 0.05, n = nrow(y_fit)))
quasi_test_FedvsFast_24

plotMD(quasi_test_FedvsFast_24)

FedvsFast_24_res <- as.data.frame(topTags(quasi_test_FedvsFast_24, p.value = 0.05, n = nrow(y_fit)))

summary(FedvsFast_24_res$logFC)

hist(FedvsFast_24_res$logFC)
#looking at the sum of up and down regulated genes for fedvsfast_24hrs
upregulated <- sum(FedvsFast_24_res$logFC > logFC_threshold)
upregulated # 73 upregulated


downregulated <- sum(FedvsFast_24_res$logFC < logFC_threshold)
downregulated # 194 downregulated

#FedvsFast_48
nrow(topTags(quasi_test_FedvsFast_48, p.value = 0.05, n = nrow(y_fit)))
quasi_test_FedvsFast_48

plotMD(quasi_test_FedvsFast_48)

FedvsFast_48_res <- as.data.frame(topTags(quasi_test_FedvsFast_48, p.value = 0.05, n = nrow(y_fit)))

summary(FedvsFast_48_res$logFC)

hist(FedvsFast_48_res$logFC)

#looking at the sum of up and down regulated genes for fedvsfast_48hrs
upregulated <- sum(FedvsFast_48_res$logFC > logFC_threshold)
upregulated # 351 upregulated

downregulated <- sum(FedvsFast_48_res$logFC < logFC_threshold)
downregulated # 463 downregulated

#FedvsFast_168
nrow(topTags(quasi_test_FedvsFast_168, p.value = 0.05, n = nrow(y_fit)))
quasi_test_FedvsFast_168

plotMD(quasi_test_FedvsFast_168)

FedvsFast_168_res <- as.data.frame(topTags(quasi_test_FedvsFast_168, p.value = 0.05, n = nrow(y_fit)))

summary(FedvsFast_168_res$logFC)

hist(FedvsFast_168_res$logFC)
#looking at the sum of up and down regulated genes for fedvsfast_168hrs
upregulated <- sum(FedvsFast_168_res$logFC > logFC_threshold)
upregulated # 1065 upregulated

downregulated <- sum(FedvsFast_168_res$logFC < logFC_threshold)
downregulated # 6427 downregulated

#FedvsFast_336
nrow(topTags(quasi_test_FedvsFast_336, p.value = 0.05, n = nrow(y_fit)))
quasi_test_FedvsFast_336

plotMD(quasi_test_FedvsFast_336)

FedvsFast_336_res <- as.data.frame(topTags(quasi_test_FedvsFast_336, p.value = 0.05, n = nrow(y_fit)))

summary(FedvsFast_336_res$logFC)

hist(FedvsFast_336_res$logFC)
#looking at the sum of up and down regulated genes for fedvsfast_336hrs
upregulated <- sum(FedvsFast_336_res$logFC > logFC_threshold)
upregulated # 2431 upregulated

downregulated <- sum(FedvsFast_336_res$logFC < logFC_threshold)
downregulated #9739 downregulated

# looking at upregulated genes at 336hrs
upregulated_genes_FedvsFast_336 <- FedvsFast_336_res[FedvsFast_336_res$logFC > logFC_threshold, ]

#FedvsFast_24
write_tsv(x = FedvsFast_24_res %>% rownames_to_column(var = "transcript_ID"), file = "FedvsFast_24.txt")

#FedvsFast_48
write_tsv(x = FedvsFast_48_res %>% rownames_to_column(var = "transcript_ID"), file = "FedvsFast_48.txt")

#FedvsFast_168
write_tsv(x = FedvsFast_168_res %>% rownames_to_column(var = "transcript_ID"), file = "FedvsFast_168.txt")

#FedvsFast_336
write_tsv(x = FedvsFast_336_res %>% rownames_to_column(var = "transcript_ID"), file = "FedvsFast_336.txt")

#removing cpms with filter for differential gene expression 
cpms_export <- data.frame(TrinityID = rownames(y),cpms,check.names = FALSE)

write.csv(cpms_export,"CPM_expression.csv", row.names = FALSE)
################################################################################################################################################################
################################################################################################################################################################
library("tximport")
library("limma")
library("locfit")
library("edgeR")
library("tidyverse")
library(dplyr)
library(stringr)
library(ggplot2)
library(ggtext)

annotations <- read_tsv("LS_OE_Clean_eggNOG_annotation.txt")

annotations <- annotations %>% dplyr::select(TrinityID, gene_ID, Preferred_name, evalue, score, Description)

cpms <- read_csv("CPM_expression.csv")
cpms <- as.data.frame(cpms)

# Define IP3/olfactory pathway keywords


olfactory_IP3 <- c(
  "v2r", "v1r", "vomeronasal", "olfc",
  "gnaq", "gna11", "gnao", "gna14",
   "plcb3", "plcb4",
   "trpc2", "itpr3",
  "olfactory receptor C family", "olfactory receptor class A")

pattern_genes_IP3 <- paste0("\\b(",paste(olfactory_IP3, collapse = "|"),")\\b")


# Identify IP3-related annotated transcripts

olfactory_master_list_IP3 <- annotations %>%
  filter(grepl(pattern_genes_IP3, Preferred_name, ignore.case = TRUE) |grepl(pattern_genes_IP3, Description, ignore.case = TRUE)) %>%
  rowwise() %>%
  mutate(QueryTerm = paste(olfactory_IP3[sapply(olfactory_IP3, function(x) {grepl(x, Preferred_name, ignore.case = TRUE) |
  grepl(x, Description, ignore.case = TRUE)}) ],
  collapse = "; ")) %>%
  ungroup() %>%
  distinct(TrinityID, .keep_all = TRUE) %>%
  mutate(evalue = as.numeric(evalue),
    score = as.numeric(score),
    GeneLabel = ifelse(
      is.na(Preferred_name) |
        Preferred_name == "" |
        Preferred_name == "-",
      Description,
      Preferred_name))

# keeping ALL isoforms Best e-value, highest score, TrinityID

olfactory_master_list_IP3 <- olfactory_master_list_IP3 %>%
  group_by(gene_ID) %>%
  arrange(evalue, desc(score), TrinityID, .by_group = TRUE) %>%
  mutate(IsoformNumber = row_number(),n_isoforms = n()) %>%
  ungroup()

olfactory_master_list_IP3 <- olfactory_master_list_IP3 %>%
  group_by(GeneLabel) %>%
  mutate(
    n_geneIDs = n_distinct(gene_ID),
    GeneGroup = dense_rank(gene_ID)
  ) %>%
  ungroup()

olfactory_master_list_IP3 <- olfactory_master_list_IP3 %>%
  mutate(GeneLetter = ifelse(n_geneIDs > 1,LETTERS[GeneGroup],""),
GeneDisplay = case_when(n_geneIDs > 1 ~ paste0(GeneLabel, " (",GeneLetter,IsoformNumber,")"),
n_isoforms > 1 ~ paste0(GeneLabel, " (",IsoformNumber,")"), TRUE ~ GeneLabel))

#  Select IP3 genes from CPM data

IP3_CPM <- cpms %>%
  filter(TrinityID %in% olfactory_master_list_IP3$TrinityID) %>%
  left_join(
    olfactory_master_list_IP3 %>%
      dplyr::select(TrinityID,gene_ID,GeneLabel,GeneDisplay, Preferred_name,Description,evalue,score),by = "TrinityID")

# Create expression matrix

IP3_expression <- IP3_CPM %>%
  dplyr::select(
    -GeneLabel,
    -GeneDisplay,
    -Preferred_name,
    -Description,
    -evalue,
    -score,
    -gene_ID,
    TrinityID
  ) %>%
  column_to_rownames("TrinityID")

# Read and organize metadata

metadata <- read_tsv("metadata.txt")

metadata <- metadata[match(colnames(IP3_expression), metadata$file_ID),]

metadata$state <- factor(metadata$state,levels = c("Fed", "Fast"))

metadata$time <- factor(metadata$time,levels = c("24hrs","48hrs","168hrs","336hrs"))

metadata_ordered <- metadata %>%
  arrange(state, time, sample)

sample_order <- metadata_ordered$file_ID

# Reorder expression matrix

IP3_expression_ordered <- IP3_expression[,sample_order,drop = FALSE]

# Create State × Time groups

group <- paste(
  metadata_ordered$state,
  metadata_ordered$time,
  sep = "_")

# Average logCPM across biological replicates

IP3_logCPM <- log2(IP3_expression_ordered + 1)

IP3_logCPM_mean <- sapply(unique(group),function(g) {rowMeans(IP3_logCPM[, group == g, drop = FALSE], na.rm = TRUE)})

IP3_logCPM_mean <- as.data.frame(IP3_logCPM_mean)

# gene-wise Z-score across the 8 group means


IP3_mean_z <- t(scale(t(IP3_logCPM_mean)))

#  Add GeneLabels

gene_labels <- IP3_CPM$GeneDisplay[match(rownames(IP3_mean_z),IP3_CPM$TrinityID)]

gene_labels <- make.unique(gene_labels)

gene_order <- order(gene_labels,na.last = TRUE)

IP3_mean_z <- IP3_mean_z[gene_order,,drop = FALSE]

rownames(IP3_mean_z) <- gene_labels[gene_order]

# Convert to long format for ggplot

IP3_heatmap_df <- as.data.frame(IP3_mean_z) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene,names_to = "Group",values_to = "Zscore") %>%
  mutate(State = ifelse(grepl("^Fed", Group),"Fed","Fast"),
    Time = case_when(
      grepl("24", Group)  ~ "24",
      grepl("48", Group)  ~ "48",
      grepl("168", Group) ~ "168",
      grepl("336", Group) ~ "336"),
    Group = factor(Group,
      levels = c(
        "Fed_24hrs",
        "Fed_48hrs",
        "Fed_168hrs",
        "Fed_336hrs",
        "Fast_24hrs",
        "Fast_48hrs",
        "Fast_168hrs",
        "Fast_336hrs")),
    Gene = factor(Gene,levels = rev(sort(unique(Gene)))))

# Heatmap

bold_genes_IP3 <- c(
  "Vomeronasal 2 receptor (C1)",
  "Vomeronasal type-1 receptor (A1)",
  "Vomeronasal type-1 receptor (B1)")

IP3_heatmap <- ggplot(IP3_heatmap_df, aes(x = Group,y = Gene,fill = Zscore)) +
  geom_tile( width = 0.95, height = 0.95) +
   scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = c(-2.5, -1.25, 0, 1.25, 2.5),
    name = "Z-score") +
  scale_x_discrete(labels = c("24", "48", "168", "336", "24", "48", "168", "336")) +
  scale_y_discrete(labels = function(x) {ifelse(x %in% bold_genes_IP3,paste0("<b>", x, "</b>"),x)})+
  annotate("text",x = 2.5, y = n_distinct(IP3_heatmap_df$Gene) + 2,label = "Fed",fontface = "bold",size = 12) +
  annotate( "text",x = 6.5,y = n_distinct(IP3_heatmap_df$Gene) + 2,label = "Fasted",fontface = "bold",size = 12) +
  annotate( "text",x = 7,y = 5.6, label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 7,y = 3.6,label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 8,y = 2.6,label = "*",fontface = "bold",size =10) +
  
  geom_vline(xintercept = 4.5,linewidth = 1) +
   labs( x = "Time (hrs)",y = "Transcript Annotation") +
  theme_classic(base_size = 24) +
   theme(axis.title.x = element_text(face = "bold", size = 26),
    axis.title.y = element_text(face = "bold",size = 26),
    axis.text.x = element_text(angle = 0,hjust = 0.5,face = "bold",size = 20),
    axis.text.y = ggtext::element_markdown(size = 15),
    axis.ticks.x = element_blank(),
    legend.title = element_text(face = "bold", size = 24),
    legend.text = element_text(face = "bold",size = 15),
    legend.position = "right",
    panel.grid = element_blank(),
    plot.margin = margin(40, 15, 15, 15)) +
  coord_cartesian(ylim = c(0,n_distinct(IP3_heatmap_df$Gene) + 2),clip = "off")

IP3_heatmap

ggsave(filename = "IP3_heatmap.pdf", plot = IP3_heatmap, dpi = 3000, scale = 3.5)
####################################################################
####################################################################
####################################################################
####################################################################

  # Define strict gene symbols using word boundaries (\\b) so "gnal" doesn't match "signal"
 Olfactory_cAMP <- c("adcy3", "gnal", "cnga2", "cnga4", "cngb1", "ano2", "omp", "taar1", "OR6N1","OR52D1","Odorant Receptor", "trace amine", "PDE1C", "CALM1","CAMK2A","CAMK2B", "CAMK2D", "CAMK2G")

pattern_genes_cAMP <- paste0("\\b(", paste( Olfactory_cAMP, collapse = "|"), ")\\b")

# Dynamic selection fix for 'Preferred_name' variations
target_col_cAMP <- colnames(annotations)[grep("preferred_name", colnames(annotations), ignore.case = TRUE)]

# Identify cAMP-related annotated transcripts

olfactory_master_list_cAMP <- annotations %>%
  filter(grepl(pattern_genes_cAMP, Preferred_name, ignore.case = TRUE) 
         |grepl(pattern_genes_cAMP, Description, ignore.case = TRUE)) %>%
  rowwise() %>%
  mutate(QueryTerm = paste(Olfactory_cAMP[sapply(Olfactory_cAMP, function(x) {grepl(x, Preferred_name, ignore.case = TRUE) |
      grepl(x, Description, ignore.case = TRUE)}) ],
      collapse = "; ")) %>%
  ungroup() %>%
  distinct(TrinityID, .keep_all = TRUE) %>%
  mutate(evalue = as.numeric(evalue),
         score = as.numeric(score),
         GeneLabel = ifelse(
           is.na(Preferred_name) |
             Preferred_name == "" |
             Preferred_name == "-",
           Description,
           Preferred_name))

olfactory_master_list_cAMP <- olfactory_master_list_cAMP %>%
  group_by(gene_ID) %>%
  arrange(evalue, desc(score), TrinityID, .by_group = TRUE) %>%
  mutate(IsoformNumber = row_number(),n_isoforms = n()) %>%
  ungroup()

olfactory_master_list_cAMP <- olfactory_master_list_cAMP %>%
  group_by(GeneLabel) %>%
  mutate(
    n_geneIDs = n_distinct(gene_ID),
    GeneGroup = dense_rank(gene_ID)
  ) %>%
  ungroup()

olfactory_master_list_cAMP <- olfactory_master_list_cAMP %>%
  mutate(GeneLetter = ifelse(n_geneIDs > 1,LETTERS[GeneGroup],""),
         
         GeneDisplay = case_when(n_geneIDs > 1 ~ paste0(GeneLabel, " (",GeneLetter,IsoformNumber,")"),
                                 n_isoforms > 1 ~ paste0(GeneLabel, " (",IsoformNumber,")"), TRUE ~ GeneLabel))

# Select cAMP genes from CPM data

cAMP_CPM <- cpms %>%
  filter(TrinityID %in% olfactory_master_list_cAMP$TrinityID) %>%
  left_join(
    olfactory_master_list_cAMP %>%
      dplyr::select(TrinityID,gene_ID,GeneLabel,GeneDisplay, Preferred_name,Description,evalue,score),by = "TrinityID")

# Create expression matrix

cAMP_expression <- cAMP_CPM %>%
  dplyr::select(
    -GeneLabel,
    -GeneDisplay,
    -Preferred_name,
    -Description,
    -evalue,
    -score,
    -gene_ID,
    TrinityID
  ) %>%
  column_to_rownames("TrinityID")


# Reorder expression matrix

cAMP_expression_ordered <- cAMP_expression[,sample_order,drop = FALSE]


# Average logCPM across biological replicates

cAMP_logCPM <- log2(cAMP_expression_ordered + 1)

cAMP_logCPM_mean <- sapply(unique(group),function(g) {rowMeans(cAMP_logCPM[, group == g, drop = FALSE], na.rm = TRUE)})

cAMP_logCPM_mean <- as.data.frame(cAMP_logCPM_mean)

# Gene-wise Z-score across the 8 group means

cAMP_mean_z <- t(scale(t(cAMP_logCPM_mean)))

# 13. Add GeneLabels

gene_labels <- cAMP_CPM$GeneDisplay[match(rownames(cAMP_mean_z),cAMP_CPM$TrinityID)]

gene_labels <- make.unique(gene_labels)

gene_order <- order(gene_labels,na.last = TRUE)

cAMP_mean_z <- cAMP_mean_z[gene_order,,drop = FALSE]

rownames(cAMP_mean_z) <- gene_labels[gene_order]

# Convert to long format for ggplot

cAMP_heatmap_df <- as.data.frame(cAMP_mean_z) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene,names_to = "Group",values_to = "Zscore") %>%
  mutate(State = ifelse(grepl("^Fed", Group),"Fed","Fast"),
         Time = case_when(
           grepl("24", Group)  ~ "24",
           grepl("48", Group)  ~ "48",
           grepl("168", Group) ~ "168",
           grepl("336", Group) ~ "336"),
         Group = factor(Group,
                        levels = c(
                          "Fed_24hrs",
                          "Fed_48hrs",
                          "Fed_168hrs",
                          "Fed_336hrs",
                          "Fast_24hrs",
                          "Fast_48hrs",
                          "Fast_168hrs",
                          "Fast_336hrs")),
         Gene = factor(Gene,levels = rev(sort(unique(Gene)))))



# Heatmap

bold_genes_cAMP <- c(
  "ADCY3",
  "CAMK2G",
  "CNGA2 (1)",
  "CNGA2 (2)",
  "Trace amine-associated receptor (B1)",
  "Odorant receptor, family H, subfamily 134, member 1 (F2)"
  )
cAMP_heatmap <- ggplot(cAMP_heatmap_df, aes(x = Group,y = Gene,fill = Zscore)) +
  geom_tile( width = 0.95, height = 0.95) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = c(-2.5, -1.25, 0, 1.25, 2.5),
    name = "Z-score") +
  scale_x_discrete(labels = c("24", "48", "168", "336", "24", "48", "168", "336")) +
  scale_y_discrete(labels = function(x) {ifelse(x %in% bold_genes_cAMP,paste0("<b>", x, "</b>"),x)})+
  annotate("text",x = 2.5, y = n_distinct(cAMP_heatmap_df$Gene) + 2,label = "Fed",fontface = "bold",size = 12) +
  annotate( "text",x = 6.5,y = n_distinct(cAMP_heatmap_df$Gene) + 2,label = "Fasted",fontface = "bold",size = 12) +
  annotate( "text",x = 7,y = 40.6, label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 7,y = 35.6,label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 7,y = 34.6,label = "*",fontface = "bold",size =10) +
  annotate( "text",x = 7,y = 33.6,label = "*",fontface = "bold",size =10) +
  annotate( "text",x = 8,y = 40.6, label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 8,y = 34.6,label = "*",fontface = "bold",size =10) +
  annotate( "text",x = 8,y = 33.6,label = "*",fontface = "bold",size =10) +
  annotate( "text",x = 8,y = 9.65,label = "*",fontface = "bold",size =10) +
  annotate( "text",x = 8,y = 0.65,label = "*",fontface = "bold",size =10) +
  
  
  geom_vline(xintercept = 4.5,linewidth = 1) +
  labs( x = "Time (hrs)",y = "Transcript Annotation") +
  theme_classic(base_size = 24) +
  theme(axis.title.x = element_text(face = "bold", size = 26),
        axis.title.y = element_text(face = "bold",size = 26),
        axis.text.x = element_text(angle = 0,hjust = 0.5,face = "bold",size = 20),
        axis.text.y = ggtext::element_markdown(size = 15),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold", size = 24),
        legend.text = element_text(face = "bold",size = 15),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(40, 15, 15, 15)) +
  coord_cartesian(ylim = c(0,n_distinct(cAMP_heatmap_df$Gene) + 2),clip = "off")

cAMP_heatmap

ggsave(filename = "cAMP_heatmap.pdf", plot = cAMP_heatmap, dpi = 3000, scale = 3.5)
########################################################################################################################################
########################################################################################################################################

# Define strict gene symbols using word boundaries (\\b) so "gnal" doesn't match "signal"
Biophysical_properties <- c("Sodium channel, voltage-gated", "voltage-dependent sodium", "Potassium voltage-gated", "voltage-dependant potassium", "voltage-gated calcium channel")

pattern_genes_Biophysical <- paste0("\\b(", paste(Biophysical_properties, collapse = "|"), ")\\b")

# Dynamic selection fix for 'Preferred_name' variations
target_col_Biophysical <- colnames(annotations)[grep("preferred_name", colnames(annotations), ignore.case = TRUE)]

#  Identify iophysical proterties annotated transcripts

olfactory_master_list_Biophysical <- annotations %>%
  filter(grepl(pattern_genes_Biophysical, Preferred_name, ignore.case = TRUE) 
         |grepl(pattern_genes_Biophysical, Description, ignore.case = TRUE)) %>%
  rowwise() %>%
  mutate(QueryTerm = paste(Biophysical_properties[sapply(Biophysical_properties, function(x) {grepl(x, Preferred_name, ignore.case = TRUE) |
      grepl(x, Description, ignore.case = TRUE)}) ],
      collapse = "; ")) %>%
  ungroup() %>%
  distinct(TrinityID, .keep_all = TRUE) %>%
  mutate(evalue = as.numeric(evalue),
         score = as.numeric(score),
         GeneLabel = ifelse(
           is.na(Preferred_name) |
             Preferred_name == "" |
             Preferred_name == "-",
           Description,
           Preferred_name))

olfactory_master_list_Biophysical <- olfactory_master_list_Biophysical[-87, ]

# keeping ALL isoform, Best e-value, highest score, TrinityID 

olfactory_master_list_Biophysical <- olfactory_master_list_Biophysical %>%
  group_by(gene_ID) %>%
  arrange(evalue, desc(score), TrinityID, .by_group = TRUE) %>%
  mutate(IsoformNumber = row_number(),n_isoforms = n()) %>%
  ungroup()

olfactory_master_list_Biophysical <- olfactory_master_list_Biophysical %>%
  group_by(GeneLabel) %>%
  mutate(
    n_geneIDs = n_distinct(gene_ID),
    GeneGroup = dense_rank(gene_ID)) %>%
  ungroup()

olfactory_master_list_Biophysical <- olfactory_master_list_Biophysical %>%
  mutate(GeneLetter = ifelse(n_geneIDs > 1,LETTERS[GeneGroup],""),
         GeneDisplay = case_when(n_geneIDs > 1 ~ paste0(GeneLabel, " (",GeneLetter,IsoformNumber,")"),
                                 n_isoforms > 1 ~ paste0(GeneLabel, " (",IsoformNumber,")"), TRUE ~ GeneLabel))

# Select Biophysical genes from CPM data

Biophysical_CPM <- cpms %>%
  filter(TrinityID %in% olfactory_master_list_Biophysical$TrinityID) %>%
  left_join(
    olfactory_master_list_Biophysical %>%
      dplyr::select(TrinityID,gene_ID,GeneLabel,GeneDisplay, Preferred_name,Description,evalue,score),by = "TrinityID")

# Create expression matrix

Biophysical_expression <- Biophysical_CPM %>%
  dplyr::select(
    -GeneLabel,
    -GeneDisplay,
    -Preferred_name,
    -Description,
    -evalue,
    -score,
    -gene_ID,
    TrinityID
  ) %>%
  column_to_rownames("TrinityID")


# Reorder expression matrix

Biophysical_expression_ordered <- Biophysical_expression[,sample_order,drop = FALSE]


# Average logCPM across biological replicates

Biophysical_logCPM <- log2(Biophysical_expression_ordered + 1)

Biophysical_logCPM_mean <- sapply(unique(group),function(g) {rowMeans(Biophysical_logCPM[, group == g, drop = FALSE], na.rm = TRUE)})

Biophysical_logCPM_mean <- as.data.frame(Biophysical_logCPM_mean)

# ene-wise Z-score across the 8 group means

Biophysical_mean_z <- t(scale(t(Biophysical_logCPM_mean)))

#Add GeneLabels

gene_labels <- Biophysical_CPM$GeneDisplay[match(rownames(Biophysical_mean_z),Biophysical_CPM$TrinityID)]

gene_labels <- make.unique(gene_labels)

gene_order <- order(gene_labels,na.last = TRUE)

Biophysical_mean_z <- Biophysical_mean_z[gene_order,,drop = FALSE]

rownames(Biophysical_mean_z) <- gene_labels[gene_order]

#Convert to long format for ggplot

Biophysical_heatmap_df <- as.data.frame(Biophysical_mean_z) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(cols = -Gene,names_to = "Group",values_to = "Zscore") %>%
  mutate(State = ifelse(grepl("^Fed", Group),"Fed","Fast"),
         Time = case_when(
           grepl("24", Group)  ~ "24",
           grepl("48", Group)  ~ "48",
           grepl("168", Group) ~ "168",
           grepl("336", Group) ~ "336"),
         Group = factor(Group,
                        levels = c(
                          "Fed_24hrs",
                          "Fed_48hrs",
                          "Fed_168hrs",
                          "Fed_336hrs",
                          "Fast_24hrs",
                          "Fast_48hrs",
                          "Fast_168hrs",
                          "Fast_336hrs")),
         Gene = factor(Gene,levels = rev(sort(unique(Gene)))))

Biophysical_heatmap_df

# heatmap
bold_genes_Biophysical <- c(
  "KCNA4",
  "KCNH2 (C1)",
  "KCNQ1 (5)",
  "KCNC2 (2)")
Biophysical_heatmap <- ggplot(Biophysical_heatmap_df, aes(x = Group,y = Gene,fill = Zscore)) +
  geom_tile( width = 0.95, height = 0.95) +
  scale_fill_gradient2(
    low = "navy",
    mid = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    breaks = c(-2.5, -1.25, 0, 1.25, 2.5),
    name = "Z-score") +
  scale_x_discrete(labels = c("24", "48", "168", "336", "24", "48", "168", "336")) +
  scale_y_discrete(labels = function(x) {ifelse(x %in% bold_genes_Biophysical,paste0("<b>", x, "</b>"),x)})+
  annotate("text",x = 2.5, y = n_distinct(Biophysical_heatmap_df$Gene) + 2,label = "Fed",fontface = "bold",size = 12) +
  annotate( "text",x = 6.5,y = n_distinct(Biophysical_heatmap_df$Gene) + 2,label = "Fasted",fontface = "bold",size = 12) +
  annotate( "text",x = 7,y = 21.6, label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 7,y = 44.6, label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 8,y = 61.6, label = "*",fontface = "bold",size = 10) +
  annotate( "text",x = 7,y = 73.6, label = "*",fontface = "bold",size = 10) +
  geom_vline(xintercept = 4.5,linewidth = 1) +
  labs( x = "Time (hrs)",y = "Transcript Annotation") +
  theme_classic(base_size = 24) +
  theme(axis.title.x = element_text(face = "bold", size = 26),
        axis.title.y = element_text(face = "bold",size = 26),
        axis.text.x = element_text(angle = 0,hjust = 0.5,face = "bold",size = 20),
        axis.text.y = ggtext::element_markdown(size = 15),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold", size = 24),
        legend.text = element_text(face = "bold",size = 15),
        legend.position = "right",
        panel.grid = element_blank(),
        plot.margin = margin(40, 15, 15, 15)) +
  coord_cartesian(ylim = c(0,n_distinct(Biophysical_heatmap_df$Gene) + 2),clip = "off")

Biophysical_heatmap

ggsave(filename = "Biophysical_heatmap2.pdf", plot = Biophysical_heatmap, dpi = 3000, scale = 3.5)

#Filter out the columns i want for downstream analysis
olfactory_master_list_Biophysical <- olfactory_master_list_Biophysical %>% dplyr::select(TrinityID, Preferred_name.x, Description.x, evalue, score)
######################################################################################################################################################################
######################################################################################################################################################################
library(ggpubr)
# Page/Figure 1: A + B
AB_plot <- ggarrange( IP3_heatmap,cAMP_heatmap, ncol = 1,nrow = 2,labels = c("A", "B"),font.label = list(size = 30, color = "black", face = "bold"),common.legend = TRUE,legend = "right",align = "v")

AB_plot

ggsave(filename = "AB_heatmaps.pdf", plot = AB_plot,dpi = 3000,scale = 5)


# Page/Figure 2: C
C_plot <- ggarrange( Biophysical_heatmap,labels = "C",font.label = list(size = 30, color = "black", face = "bold"),common.legend = TRUE, legend = "right")

C_plot

ggsave(filename = "C_heatmap.pdf",plot = C_plot,dpi = 3000,scale = 5)
