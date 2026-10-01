library("tximport")
library("limma")
library("locfit")
library("edgeR")
library("tidyverse")
library(dplyr)
library(ggplot2)
library(reshape2)
library("clusterProfiler")
library("GO.db")
library(KEGGREST)

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
normMat <- txi.salmon$length

dim(cts)

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

# Smoother, shorter, and won't break on future package updates
metadata <- tibble(file_ID = list.files("./sturgeon_quant/")) %>% 
  separate(col = file_ID, into = c("state", "time", "sample"), sep = "_", remove = FALSE) %>% 
  mutate(treatment = as.factor(paste0(state, "_", time)))

# filtering using the design information
design <- model.matrix(~0 + treatment, data = metadata)
keep <- filterByExpr(y, design)
y <- y[keep, ]

#For creating a matrix of CPMs within edgeR, the following code chunk can be used: since we already removed zero across transcripts and samples
cpms <- edgeR::cpm(y, offset = y$offset, log = FALSE)

#portion of creating for negative binomial distribution estimateDisp
y<-estimateDisp(y,design)
y

y_fit <- glmQLFit(y, design, robust = TRUE)

##################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################### time course analysis using a joint design matrix#####################################################################
#############################################################################################################################################################################################

# Match metadata rows to the order of count-matrix columns
metadata_time <- metadata_time[
  match(colnames(y), metadata_time$file_ID),]

#Subset count matrix 'y' to match the cleaned metadata
#y_clean <- y[, !is.na(metadata_time$time)]

# Subset count matrix using the SAME rows/samples
y_clean <- y[, metadata_time$file_ID %in% metadata_clean$file_ID]

#Generate polynomial terms based on the continuous time variable
time_poly <- poly(metadata_clean$time, df=3)

#Define the state factor
state_factor <- factor(metadata_clean$state, levels = c("Fed", "Fast"))

#Create a design matrix with interaction terms
# This fits a base curve for 'Fed' and calculates how 'Fast' deviates from it
design_joint <- model.matrix(~ state_factor * time_poly)
colnames(design_joint)
design_joint

#Estimate Dispersion and Fit the Model
# Estimate global and gene-wise dispersion
y_clean <- estimateDisp(y_clean, design_joint)

# Convert the row names into a proper data frame
background_df <- data.frame(gene = rownames(y_clean))

# Save 
write.csv(background_df, "global_background_universe.csv", row.names = FALSE)

# Fit the Quasi-Likelihood GLM
fit_joint <- glmQLFit(y_clean, design_joint, robust=TRUE)

#Extracting Unique and Common Genes Mathematically
# Test all interaction coefficients simultaneously
fit_unique_response <- glmQLFTest(fit_joint, coef=6:8)

# These genes change over time in a way that is unique to their group
summary(decideTests(fit_unique_response))
time_res_unique <- as.data.frame(topTags(fit_unique_response, p.value = 0.05, n = nrow(fit_joint)))

#Finding Genes with a common Response Over Time
#Test if time matters at all (overall time effect across both groups)
fit_any_time <- glmQLFTest(fit_joint, coef=c(3:5, 6:8))
de_time_genes <- rownames(topTags(fit_any_time, p.value = 0.05, n = nrow(fit_joint)))

#Extract genes that are significant for time but not significant for interaction
interaction_genes <- rownames(time_res_unique) # From step A above

common_time_genes <- setdiff(de_time_genes, interaction_genes)

#Extracting Log-CPM and Fit Values
time_res_joint.obs <- cpm(y_clean, log = TRUE, prior.count = fit_joint$prior.count)

# Fitted values based on the model
time_res_joint.fit <- cpm(fit_joint, log = TRUE)

#################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################### Prepare the True 4-Timepoint Profile Matrix trying clustering instead of the timeseries, since it is missing key points##
# creating the different clusters of transcript profiles
#Identify sample indices for every specific condition
idx_fed_24   <- which(metadata_clean$time == 24   & metadata_clean$state == "Fed")
idx_fed_48   <- which(metadata_clean$time == 48   & metadata_clean$state == "Fed")
idx_fed_168  <- which(metadata_clean$time == 168  & metadata_clean$state == "Fed")
idx_fed_336  <- which(metadata_clean$time == 336  & metadata_clean$state == "Fed")

idx_fast_24  <- which(metadata_clean$time == 24   & metadata_clean$state == "Fast")
idx_fast_48  <- which(metadata_clean$time == 48   & metadata_clean$state == "Fast")
idx_fast_168 <- which(metadata_clean$time == 168  & metadata_clean$state == "Fast")
idx_fast_336 <- which(metadata_clean$time == 336  & metadata_clean$state == "Fast")

# 2. Extract unique/interaction gene names edgeR step
interaction_genes <- rownames(time_res_unique)

# 3. Build a simplified matrix of mean fitted expression values across the timeline
mean_profiles <- data.frame(
  Fed_24h    = rowMeans(time_res_joint.fit[interaction_genes, idx_fed_24,   drop=FALSE]),
  Fed_48h    = rowMeans(time_res_joint.fit[interaction_genes, idx_fed_48,   drop=FALSE]),
  Fed_168h   = rowMeans(time_res_joint.fit[interaction_genes, idx_fed_168,  drop=FALSE]),
  Fed_336h   = rowMeans(time_res_joint.fit[interaction_genes, idx_fed_336,  drop=FALSE]),
  Fasted_24h = rowMeans(time_res_joint.fit[interaction_genes, idx_fast_24,  drop=FALSE]),
  Fasted_48h = rowMeans(time_res_joint.fit[interaction_genes, idx_fast_48,  drop=FALSE]),
  Fasted_168h= rowMeans(time_res_joint.fit[interaction_genes, idx_fast_168, drop=FALSE]),
  Fasted_336h= rowMeans(time_res_joint.fit[interaction_genes, idx_fast_336, drop=FALSE])
)

#Scale rows to Z-scores (normalizes baseline differences to focus strictly on trajectory shapes)
profiles_scaled <- t(scale(t(as.matrix(mean_profiles))))

#evaluating the amount of clusters needed to determine a response
install.packages("factoextra")
install.packages("rlang")
library(rlang)
library(factoextra)
#evaluating different values of K
fviz_nbclust(profiles_scaled,kmeans,method = "wss")

#comparing the silhouette
sil <- fviz_nbclust(profiles_scaled,kmeans,method = "silhouette")
sil
sil$data

#Run the K-Means Clustering Algorithm

# Set a random seed so cluster numbers stay identical every time the script is run
set.seed(123)

# Run K-means clustering to split genes into 4 shape groups
cluster_results <- kmeans(profiles_scaled, centers = 4, iter.max = 100, nstart = 25)

# Add the cluster assignments back to data as a tracking dataframe
gene_clusters <- data.frame(
  GeneID = names(cluster_results$cluster),
  Cluster = factor(cluster_results$cluster))

#Generate a Line Plot Grid (Visual Differentiation)
# 1. Set up a 2x2 grid layout for the 4 clusters
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), xpd = NA)

# Define our 4 true sequential timepoints on the x-axis
x_points <- 1:4

# Loop through each of the 4 clusters to draw their patterns correctly
for (c in 1:4) {
  # Extract the scaled profile data for genes in this specific cluster
  cluster_genes <- gene_clusters$GeneID[gene_clusters$Cluster == c]
  cluster_data  <- profiles_scaled[cluster_genes, , drop = FALSE]
  
  # Separate the data into the two distinct fish groups
  fed_matrix   <- cluster_data[, 1:4, drop = FALSE]   # Fed 24, 48, 168, 336
  fast_matrix  <- cluster_data[, 5:8, drop = FALSE]   # Fasted 24, 48, 168, 336
  
  # Calculate the average trajectory line for each group separately
  fed_mean  <- colMeans(fed_matrix)
  fast_mean <- colMeans(fast_matrix)
  
  # Create an empty plot canvas
  plot(x_points, fed_mean, type = "n", xaxt = "n",
       ylim = c(-2.5, 2.5), xlab = "Time (Hours)", ylab = "Z-score",
       main = paste("Cluster", c, "(n =", length(cluster_genes), "Transcripts)"))
  
  # Set the x-axis labels to show the true time intervals
  axis(1, at = 1:4, labels = c("24", "48", "168", "336"))
  
  # Draw the background lines for every individual gene
  for (i in 1:nrow(cluster_data)) {
    lines(x_points, fed_matrix[i, ],  col = rgb(0.8, 0.1, 0.1, alpha = 0.08)) # Red for Fed
    lines(x_points, fast_matrix[i, ], col = rgb(0.1, 0.4, 0.8, alpha = 0.08)) # Blue for Fasted}
  
  # Overlay trend lines
  lines(x_points, fed_mean,  col = "firebrick",  lwd = 4)   # Red line for Fed cohort
  lines(x_points, fast_mean, col = "dodgerblue", lwd = 4)   # Blue line for Fasted cohort
  
  # Place a legend on the first cluster plot
  if(c == 1) {
    legend("topleft",
           inset = c(-1, -1.5),
           legend = c("Fed", "Fast"),
           col = c("firebrick", "dodgerblue"),
           lwd = 3,
           bty = "n",
           cex = 0.8)}}
# Reset plotting layout back to standard 1x1 default
par(mfrow = c(1, 1))

# Extract target clusters
cluster1_IDs <- gene_clusters$GeneID[gene_clusters$Cluster == 1] 
cluster2_IDs <- gene_clusters$GeneID[gene_clusters$Cluster == 2] 
cluster3_IDs <- gene_clusters$GeneID[gene_clusters$Cluster == 3] 
cluster4_IDs <- gene_clusters$GeneID[gene_clusters$Cluster == 4] 
# Extract target clusters

#Clean the column names of clustering object
colnames(gene_clusters) <- c("GeneID", "Cluster")

#Merge cluster numbers with full edgeR time-series statistics
master_cluster_df <- merge(
  gene_clusters, 
  time_res_unique, 
  by.x = "GeneID", 
  by.y = "row.names", 
  all.x = TRUE
)
master_cluster_df <- as.data.frame(master_cluster_df)
write_tsv(x = master_cluster_df %>% rownames_to_column(var = "Transcript"), file = "master_cluster.txt")

#Loop through numbers 1 to 4 to filter and save each cluster separately
for (c in 1:4) {
  # Subset the master dataframe for just the current cluster
  cluster_subset <- master_cluster_df[master_cluster_df$Cluster == c, ]
  
  # Generate a dynamic file name for each cluster
  file_name <- paste0("Cluster_", c, "_Data.csv")
  
  # Export the file 
  write.csv(cluster_subset, file = file_name, row.names = FALSE)

#################################################################################################################################################################
#################################################################################################################################################################
#entap annotated transcriptome
annotation <- read_tsv("LS_OE_Clean_eggNOG_annotation.txt") 

#go terms from the annotated transcriptome
GoTerms <- read_tsv("eggnog_annotated_gene_ontology_terms.tsv")

#all genes in the transcriptome 
background_df <-read_csv("global_background_universe.csv")

GoTerms <- GoTerms %>%
  filter(!str_detect(go_term, regex("^obsolete", ignore_case = TRUE)))

# Extract unique Trinity IDs into a new data frame called background
Back_Ground <- data.frame(gene = unique(background_df$gene))

#Extract the column as a text vector for clusterProfiler
global_background_universe <- Back_Ground %>% pull(gene)

#Format EnTAP data frame specifically for clusterProfiler
# (Column 1 MUST be go_id, Column 2 MUST be the gene ID)
custom_term2gene <- GoTerms %>% 
  dplyr::select(go_id, query_sequence) %>%
  # Keep only GO IDs that exist in the active GO database
  filter(go_id %in% keys(GO.db, keytype = "GOID"))

Cluster1 <- read_csv("Cluster_1_Data.csv")
Cluster2 <- read_csv("Cluster_2_Data.csv")
Cluster3 <- read_csv("Cluster_3_Data.csv")
Cluster4 <- read_csv("Cluster_4_Data.csv")
Master_Cluster <- read_tsv("master_cluster.txt")

#Perform the Left Join matching GeneID to TrinityID
Master_Annotated <- Master_Cluster %>%
  left_join(annotation, by = c("GeneID" = "TrinityID"))

#Filter out all forms of missing annotations to count true matches
Master_Annotations_Clean <- Master_Annotated %>%
  filter(
    (!is.na(Preferred_name) & 
       trimws(Preferred_name) != "" &
       !tolower(trimws(Preferred_name)) %in% c("none", "-")) |
      (!is.na(GOs) & trimws(GOs) != "") |
      (!is.na(KEGG_ko) & trimws(KEGG_ko) != ""))

#Get the number of successfully annotated genes
Master_annotated <- nrow(Master_Annotations_Clean)

#Calculate the percentage of the cluster that got annotated
pct_annotated <- (Master_annotated / nrow(Master_Annotated)) * 100

#Print the summary 
cat("--- Master Cluster Annotation Summary ---\n")
cat("Total Genes master cluster:     ", nrow(Master_Annotated), "\n") #3427
cat("Successfully Annotated Genes: ", Master_annotated, "\n") #1002
cat("Annotation Rate:              ", round(pct_annotated, 2), "%\n") #29.24%

#Perform the Left Join matching GeneID to TrinityID
Cluster1_annotated <- Cluster1 %>%
  left_join(annotations_new, by = c("GeneID" = "TrinityID"))

#Filter out all forms of missing annotations to count true matches
Cluster1_Annotations_Clean <- Cluster1_annotated %>%
  filter(
    (!is.na(Preferred_name) & 
       trimws(Preferred_name) != "" &
       !tolower(trimws(Preferred_name)) %in% c("none", "-")) |
      (!is.na(GOs) & trimws(GOs) != "") |
      (!is.na(KEGG_ko) & trimws(KEGG_ko) != ""))

# Get the number of successfully annotated genes
Cluster1_total_annotated <- nrow(Cluster1_Annotations_Clean)

#Calculate the percentage of the cluster that got annotated
Cluster1_pct_annotated <- (Cluster1_total_annotated / nrow(Cluster1_annotated)) * 100

#Print the summary 
cat("--- Cluster 1 Annotation Summary ---\n")
cat("Total Genes in Cluster 1:     ", nrow(Cluster1_annotated), "\n") #540
cat("Successfully Annotated Genes: ", Cluster1_total_annotated, "\n") #91
cat("Annotation Rate:              ", round(Cluster1_pct_annotated, 2), "%\n") #16.85%

#Pull unique transcript IDs 
Cluster1_ID <- Cluster1 %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run enrichment tool (All)
Cluster1_enrich_results <- enricher(
  gene         = Cluster1_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = custom_term2gene,           
  pvalueCutoff = 0.05,                       
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 1000                        
)
head(Cluster1_enrich_results)

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster1_enrichr <- as.data.frame(Cluster1_enrich_results)%>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

# Write the ID and pvalue columns straight to the file
write.table(Cluster1_enrichr[, c("ID", "pvalue")], 
            file = "Cluster1_enrichr.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Perform the Left Join matching GeneID to TrinityID
Cluster2_annotated <- Cluster2 %>%
  left_join(annotations_new, by = c("GeneID" = "TrinityID"))

# Filter out all forms of missing annotations to count true matches
Cluster2_Annotations_Clean <- Cluster2_annotated %>%
  filter(
    (!is.na(Preferred_name) & 
       trimws(Preferred_name) != "" &
       !tolower(trimws(Preferred_name)) %in% c("none", "-")) |
      (!is.na(GOs) & trimws(GOs) != "") |
      (!is.na(KEGG_ko) & trimws(KEGG_ko) != "")
  )
# Get the number of successfully annotated genes
Cluster2_total_annotated <- nrow(Cluster2_Annotations_Clean)

# Calculate the percentage of the cluster that got annotated
pct_annotated <- (Cluster2_total_annotated / nrow(Cluster2_annotated)) * 100

# Print the summary report 
cat("--- Cluster 2 Annotation Summary ---\n")
cat("Total Genes in Cluster 2:     ", nrow(Cluster2_annotated), "\n") #1861
cat("Successfully Annotated Genes: ", Cluster2_total_annotated, "\n") #716
cat("Annotation Rate:              ", round(pct_annotated, 2), "%\n") #38.47%

#Pull unique transcript IDs 
Cluster2_ID <- Cluster2_annotated %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool (All)
Cluster2_enrich_results <- enricher(
  gene         = Cluster2_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = custom_term2gene,           
  pvalueCutoff = 0.05,                       
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 1000                         
)
head(Cluster2_enrich_results)

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster2_enrichr <- as.data.frame(Cluster2_enrich_results)%>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

# Write just the ID and pvalue columns straight to the file
write.table(Cluster2_enrichr[, c("ID", "pvalue")], 
            file = "Cluster2_enrichr.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Perform the Left Join matching GeneID to TrinityID
Cluster3_annotated <- Cluster3 %>%
  left_join(annotations_new, by = c("GeneID" = "TrinityID"))

#Filter out all forms of missing annotations to count true matches
Cluster3_Annotations_Clean <- Cluster3_annotated %>%
  filter(
    (!is.na(Preferred_name) & 
       trimws(Preferred_name) != "" &
       !tolower(trimws(Preferred_name)) %in% c("none", "-")) |
      (!is.na(GOs) & trimws(GOs) != "") |
      (!is.na(KEGG_ko) & trimws(KEGG_ko) != "")
  )
#Get the number of successfully annotated genes
Cluster3_total_annotated <- nrow(Cluster3_Annotations_Clean)

#Calculate the percentage of the cluster that got annotated
Cluster3_pct_annotated <- (Cluster3_total_annotated / nrow(Cluster3_annotated)) * 100

# Print the summary 
cat("--- Cluster 3 Annotation Summary ---\n")
cat("Total Genes in Cluster 2:     ", nrow(Cluster3_annotated), "\n") #572
cat("Successfully Annotated Genes: ", Cluster3_total_annotated, "\n") #162
cat("Annotation Rate:              ", round(Cluster3_pct_annotated, 2), "%\n") #29.02%

#Pull unique transcript IDs 
Cluster3_ID <- Cluster3_annotated %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool (All)
Cluster3_enrich_results <- enricher(
  gene         = Cluster3_ID,       
  universe     = global_background_universe, 
  TERM2GENE    = custom_term2gene,           
  pvalueCutoff = 0.05,                       
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 1000                         
)
head(Cluster3_enrich_results)

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster3_enrichr <- as.data.frame(Cluster3_enrich_results)%>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

Cluster3_enrichr <- Cluster3_enrichr %>%
  mutate(ID = if_else(ID == "GO:0062023","GO:0031012",ID))

# Write just the ID and pvalue columns straight to the file
write.table(Cluster3_enrichr[, c("ID", "pvalue")], 
            file = "Cluster3_enrichr.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#Perform the Left Join matching GeneID to TrinityID
 Cluster4_annotated <- Cluster4 %>%
  left_join(annotations_new, by = c("GeneID" = "TrinityID"))

# Filter out all forms of missing annotations to count true matches
Cluster4_Annotations_Clean <- Cluster4_annotated %>%
  filter(
    (!is.na(Preferred_name) & 
       trimws(Preferred_name) != "" &
       !tolower(trimws(Preferred_name)) %in% c("none", "-")) |
      (!is.na(GOs) & trimws(GOs) != "") |
      (!is.na(KEGG_ko) & trimws(KEGG_ko) != ""))
# Get the number of successfully annotated genes
Cluster4_total_annotated <- nrow(Cluster4_Annotations_Clean)

# Calculate the percentage of the cluster that got annotated
Cluster4_pct_annotated <- (Cluster4_total_annotated / nrow(Cluster4_annotated)) * 100

# Print the summary 
cat("--- Cluster 4 Annotation Summary ---\n")
cat("Total Genes in Cluster 4:     ", nrow(Cluster4_annotated), "\n") #454
cat("Successfully Annotated Genes: ", Cluster4_total_annotated, "\n") #29
cat("Annotation Rate:              ", round(Cluster4_pct_annotated, 2), "%\n") #6.39%

#Pull UNIQUE transcript IDs 
Cluster4_ID <- Cluster4 %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool (All)
Cluster4_enrich_results <- enricher(
  gene         = Cluster4_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = custom_term2gene,          
  pvalueCutoff = 0.05,                       
  qvalueCutoff = 0.05,
  minGSSize = 5,
  maxGSSize = 1000 
)
head(Cluster4_enrich_results)

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster4_enrichr <- as.data.frame(Cluster4_enrich_results)%>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

# Write just the ID and pvalue columns straight to the file
write.table(Cluster4_enrichr[, c("ID", "pvalue")], 
            file = "Cluster4_enrichr.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#############################################################################################################################################################################################################################################################################################################################################################################################################################################################                          takeing the revigo results inputting them back into r to make graphs       ###
#################################################################################################################################################################

#all enrichr files were for clusters were imported into revigo. Then revigo files were re imported into Rstudio in the different Gene otngeny categories (GO Terms)

#############################################################################################################################################################################################################################################################################################################################################################################################################################################################                          takeing the revigo results inputting them back into r to make graphs       ###
############################################################################################################################################################################################################################################################################################################################################################
#selecting the top 15 combined score go terms in each Revigo Cluster 1

#calling the revigo analysis
Cluster1_BP <- read_tsv("Cluster/Cluster#1_BP.tsv")

#joining the annotated file with the revigo file
Cluster1_BP_GO <- Cluster1_BP %>% left_join(Cluster1_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster1_BP_GO <- Cluster1_BP_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore,Count)

#calling the revigo analysis
Cluster1_MF <- read_tsv("Cluster/Cluster#1_MF.tsv")

#joinging the annotated file with the revigo file
Cluster1_MF_GO <- Cluster1_MF %>% left_join(Cluster1_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster1_MF_GO <- Cluster1_MF_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)

#calling the revigo analysis
Cluster1_CC <- read_tsv("Cluster/Cluster#1_CC.tsv")

#joinging the annotated file with the revigo file
Cluster1_CC_GO <- Cluster1_CC %>% left_join(Cluster1_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster1_CC_GO <- Cluster1_CC_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)


#selecting the top 15 combined score go terms in each Revigo file
top20_BP_Cluster1 <- Cluster1_BP_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

top20_MF_Cluster1 <- Cluster1_MF_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

top20_CC_Cluster1 <- Cluster1_CC_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

######################################################################################################################################################
#####################################################################################################################################################
#selecting the top 15 combined score go terms in each Revigo Cluster 2

#calling the revigo analysis
Cluster2_BP <- read_tsv("Cluster/Cluster#2_BP.tsv")

#joinging the annotated file with the revigo file
Cluster2_BP_GO <- Cluster2_BP %>% left_join(Cluster2_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster2_BP_GO <- Cluster2_BP_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)


#calling the revigo analysis
Cluster2_MF <- read_tsv("Cluster/Cluster#2_MF.tsv")

#joinging the annotated file with the revigo file
Cluster2_MF_GO <- Cluster2_MF %>% left_join(Cluster2_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster2_MF_GO <- Cluster2_MF_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)

#calling the revigo analysis
Cluster2_CC <- read_tsv("Cluster/Cluster#2_CC.tsv")

#joinging the annotated file with the revigo file
Cluster2_CC_GO <- Cluster2_CC %>% left_join(Cluster2_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster2_CC_GO <- Cluster2_CC_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)


#selecting the top 15 combined score go terms in each Revigo file
top20_BP_Cluster2 <- Cluster2_BP_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

top20_MF_Cluster2 <- Cluster2_MF_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

top20_CC_Cluster2 <- Cluster2_CC_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

################################################################################################################################################################################################################################################################################################################################################
#################################################################################################################################################################################
#selecting the top 15 combined score go terms in each Revigo Cluster 3
#calling the revigo analysis
Cluster3_BP <- read_tsv("Cluster/Cluster#3_BP.tsv")

#joinging the annotated file with the revigo file
Cluster3_BP_GO <- Cluster3_BP %>% left_join(Cluster3_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster3_BP_GO <- Cluster3_BP_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)

#calling the revigo analysis
Cluster3_MF <- read_tsv("Cluster/CLuster#3_MF.tsv")

#joinging the annotated file with the revigo file
Cluster3_MF_GO <- Cluster3_MF %>% left_join(Cluster3_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster3_MF_GO <- Cluster3_MF_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)

#calling the revigo analysis
Cluster3_CC <- read_tsv("Cluster/Cluster#3_CC.tsv")

#joinging the annotated file with the revigo file
Cluster3_CC_GO <- Cluster3_CC %>% left_join(Cluster3_enrichr, by = c("TermID" = "ID"))

#reducing the size of the file
Cluster3_CC_GO <- Cluster3_CC_GO %>% dplyr::select(TermID, Name, Value, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust, CombinedScore, Count)

Cluster3_CC_GO %>%
  filter(is.na(CombinedScore)) %>%
  dplyr::select(TermID, Name)

top20_BP_Cluster3 <- Cluster3_BP_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

top20_MF_Cluster3 <- Cluster3_MF_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)

top20_CC_Cluster3 <- Cluster3_CC_GO %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 15)
################################################################################################################################################################################################################################################################################################################################################
########################################################################################################################################################################################################################################################################################
#selecting the top 15 biological function genes from each cluster and putting them into graph
#calling the revigo analysis

#combined into a single data frame
Biological <- bind_rows(
  mutate(top20_BP_Cluster1, Treatment = "Cluster #1"),
  mutate(top20_BP_Cluster2, Treatment = "Cluster #2"),
  mutate(top20_BP_Cluster3, Treatment = "Cluster #3"))

# Reshape data for plotting
Biological_long <- melt(Biological, id.vars = c("TermID", "Name", "Treatment"), measure.vars = "CombinedScore")

# Set the levels of Treatment to ensure the legend reads Fast, Fed, Common
Biological$Treatment <- factor(Biological$Treatment, levels = c("Cluster #1", "Cluster #2", "Cluster #3"))

# Create a combined factor to order the names by treatment 
Biological$Name <- factor(Biological$Name, levels = Biological$Name[order(factor(Biological$Treatment, levels = c("Cluster #1", "Cluster #2", "Cluster #3")), Biological$Enrichment)])

#reverse the order of the go term plots so that it follows the legend
Biological$Name <- factor(
  Biological$Name,
  levels = rev(Biological$Name))

#lollipop graph 
 Biological_Graph_lollipop <- ggplot(Biological, aes(y = Name, x = CombinedScore, color = Treatment)) +
  geom_segment(aes(x = 0, xend = CombinedScore, y = Name, yend = Name), linewidth = 1) +
  geom_point(aes(size = Count)) +
  scale_color_manual(values = c(
    "Cluster #1" = "black",
    "Cluster #2" = "#d73027",
    "Cluster #3" = "#4575b4")) +
  scale_size_continuous(range = c(1, 10)) + 
  labs(
    y = "Biological Process",
    x = "Combined Score",
    color = "Dataset",
    size = "Transcripts") +
   xlim(0, 5000) +
  theme_classic(base_size = 20) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title = element_text(size = 26, face = "bold"),
    legend.title = element_text(face = "bold", size = 26),
    legend.text  = element_text(face = "bold", size = 22))

Biological_Graph_lollipop

# Save the plot
ggsave(filename = "Biological_Graph.pdf", plot = Biological_Graph, dpi = 3000, scale = 3.5)

################################################################################################################################################################################################################################################################################################################################################
#######################################################################################################################################################################################
#selecting the top 15 biological function genes fro meach cluster and putting them into graph
#calling the revigo analysis

#combined into a single data frame
Cellular <- bind_rows(
  mutate(top20_CC_Cluster1, Treatment = "Cluster #1"),
  mutate(top20_CC_Cluster2, Treatment = "Cluster #2"),
  mutate(top20_CC_Cluster3, Treatment = "Cluster #3"))

#Create an invisible unique ID column for plotting to keep the groups separate
Cellular <- Cellular %>%
  mutate(Plot_ID = row_number())

#Factor the Plot_ID so ggplot respects the correct sorted sequence
Cellular$Plot_ID <- factor(Cellular$Plot_ID, levels = unique(Cellular$Plot_ID))

#highest values at the top of the plot
Cellular$Plot_ID <- factor(Cellular$Plot_ID, levels = rev(levels(Cellular$Plot_ID)))

# Set the levels of Treatment to ensure the legend reads Fast, Fed, Common
Cellular$Treatment <- factor(Cellular$Treatment, levels = c("Cluster #1", "Cluster #2", "Cluster #3"))

# Create a combined factor to order the names by treatment 
Cellular$Name <- factor(Cellular$Name, levels = unique(Cellular$Name[order(factor(Cellular$Treatment, levels = c("Cluster #1", "Cluster #2", "Cluster #3")), Cellular$Enrichment)]))

#reverse the order of the go term plots so that it follows the legend
Cellular$Name <- factor(
  Cellular$Name,
  levels = rev(unique(Cellular$Name)))

#lollipop graph fro cellular response
Cellular_Graph_lollipop <- ggplot(Cellular, aes(y = Plot_ID, x = CombinedScore, color = Treatment)) +
  geom_segment(aes(x = 0, xend = CombinedScore, y = Plot_ID, yend = Plot_ID), linewidth = 1) +
  geom_point(aes(size = Count)) +
  scale_y_discrete(labels = setNames(as.character(Cellular$Name), Cellular$Plot_ID)) +
  scale_color_manual(values = c(
    "Cluster #2" = "#d73027",
    "Cluster #3" = "#4575b4",
    "Cluster #1" = "black" )) +
  scale_size_continuous(range = c(1, 10)) + 
  labs(
    y = "Cellular Components",
    x = "Combined Score",
    color = "Dataset",
    size = "Transcripts") +
  theme_classic(base_size = 20) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title = element_text(size = 26, face = "bold"),
    legend.title = element_text(face = "bold", size = 26),
    legend.text  = element_text(face = "bold", size = 22))

Cellular_Graph_lollipop
ggsave(filename = "Cellular_Graph_lollipop.pdf", plot = Cellular_Graph_lollipop, dpi = 3000, scale = 3.5)
################################################################################################################################################################################################################################################################################################################################################
##################################################################################################################################################################################################################
#selecting the top 15 biological function genes fro meach cluster and putting them into graph
#calling the revigo analysis

#combined into a single data frame
Molecular <- bind_rows(
  mutate(top20_MF_Cluster1, Treatment = "Cluster #1"),
  mutate(top20_MF_Cluster2, Treatment = "Cluster #2"),
  mutate(top20_MF_Cluster3, Treatment = "Cluster #3"))

#Create an invisible unique ID column for plotting to keep the groups separate
Molecular <- Molecular %>%
  mutate(Plot_ID = row_number())

#Factor the Plot_ID so ggplot respects the correct sorted sequence
Molecular$Plot_ID <- factor(Molecular$Plot_ID, levels = unique(Molecular$Plot_ID))

# highest values at the top of the plot
Molecular$Plot_ID <- factor(Molecular$Plot_ID, levels = rev(levels(Molecular$Plot_ID)))

# Set the levels of Treatment to ensure the legend reads Fast, Fed, Common
Molecular$Treatment <- factor(Molecular$Treatment, levels = c("Cluster #1", "Cluster #2", "Cluster #3"))

# Create a combined factor to order the names by treatment 
Molecular$Name <- factor(Molecular$Name, levels = unique(Molecular$Name[order(factor(Molecular$Treatment, levels = c("Cluster #1", "Cluster #2", "Cluster #3")), Molecular$Enrichment)]))

#reverse the order of the go term plots so that it follows the legend
Molecular$Name <- factor(
  Molecular$Name,
  levels = rev(unique(Molecular$Name)))

#lollipop graph fro Molecular response
Molecular_Graph_lollipop <- ggplot(Molecular, aes(y = Plot_ID, x = CombinedScore, color = Treatment)) +
  geom_segment(aes(x = 0, xend = CombinedScore, y = Plot_ID, yend = Plot_ID), 
               linewidth = 1) +
  geom_point(aes(size = Count)) +
  scale_y_discrete(labels = setNames(as.character(Molecular$Name), Molecular$Plot_ID)) +
  scale_color_manual(values = c(
    "Cluster #2" = "#d73027",
    "Cluster #3" = "#4575b4",
    "Cluster #1" = "black")) +
  scale_size_continuous(range = c(1, 10)) + 
  labs(
    y = "Molecular Components",
    x = "Combined Score",
    color = "Dataset",
    size = "Transcripts") +
  theme_classic(base_size = 20) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title = element_text(size = 26, face = "bold"),
    legend.title = element_text(face = "bold", size = 26),
    legend.text  = element_text(face = "bold", size = 22))

Molecular_Graph_lollipop
ggsave(filename = "Molecular_Graph_lollipop.pdf", plot = Molecular_Graph_lollipop, dpi = 3000, scale = 3.5)

Molecular_Graph_lollipop
Cellular_Graph_lollipop
Biological_Graph_lollipop

ggsave(filename = "Molecular_Graph_lollipop.pdf", plot = Molecular_Graph_lollipop, dpi = 3000, scale = 3.5)

# making a combined plot of all go term graphs
library(patchwork)
Combined_Plot <- (
  (Cellular_Graph_lollipop + 
     labs(x = NULL) + 
     guides(color = "none", size = "none") + 
     xlim(0, 7000)) / (Biological_Graph_lollipop + 
       labs(x = NULL) + 
       xlim(0, 7000)) / (Molecular_Graph_lollipop + 
       guides(color = "none", size = "none") + 
       xlim(0, 7000))) + 
  plot_layout(heights = c(nrow(Cellular), nrow(Biological), nrow(Molecular))
  ) & theme(
    axis.text.y = element_text(size = 23, color = "black"),
    axis.text.x = element_text(size = 23, color = "black"),
    axis.title  = element_text(size = 32, face = "bold", color = "black"),
    legend.title = element_text(size = 28, face = "bold", color = "black"))

Combined_Plot

ggsave(filename = "combined_lollipop.pdf", plot = Combined_Plot, dpi = 3000, scale = 5)

################################################################################################################################
################################################################################################################################

annotation <- read_tsv("../LS_OE_Clean_eggNOG_annotation.txt") 

###all KEGG pathways from ENTAP
KEGG <- read_tsv("../LS_OE_KEGG_eggNOG_annotation.txt")

kegg_clean_universe <- KEGG %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "" & KEGG_Pathway != "NA"& KEGG_Pathway != "-")

# Column 1 Pathway ID, Column 2 Transcript ID
kegg_term2gene <- kegg_clean_universe %>% 
  dplyr::select(KEGG_Pathway, TrinityID) %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "")

#Extract only the 'ko' and 'K' orthology terms 
ko_term2gene <- kegg_term2gene %>%
  filter(str_detect(KEGG_Pathway, "(?i)^ko")) 

ko_term2gene <- kegg_term2gene %>%
  filter(str_detect(KEGG_Pathway, "(?i)^ko")) %>%
  distinct(KEGG_Pathway, TrinityID)

# Extract only the 'map' terms 
map_term2gene <- kegg_term2gene %>%
  filter(str_detect(KEGG_Pathway, "(?i)^map"))

#kegg pathways list
kegg_pathway_list <- keggList("pathway", "ko")

kegg_descriptions <- data.frame(
  KEGG_ID      = names(kegg_pathway_list),
  Description = as.character(kegg_pathway_list),
  stringsAsFactors = FALSE)

# run enrichment and combine
Cluster1 <- read_csv("../Cluster_1_Data.csv")

Cluster1_annotated <- Cluster1 %>%
  left_join(annotation, by = c("GeneID" = "TrinityID"))

# Pull unique transcript IDs Cluster1 data frame
Cluster1_ID <- Cluster1_annotated %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool
Cluster1_KEGG_results<- enricher(
  gene         = Cluster1_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = ko_term2gene,           
  pvalueCutoff = 0.05,                       
  pAdjustMethod = "BH",
  minGSSize = 3)

# Convert to data frame
Cluster1_KEGG_results <- as.data.frame(Cluster1_KEGG_results)
head(Cluster1_KEGG_results)

# Add descriptions using the verified clean key
Cluster1_KEGG_results <- Cluster1_KEGG_results %>%
  left_join(kegg_descriptions, by = c("ID" = "KEGG_ID")) %>%
  dplyr::select(ID, Description.y, everything())

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster1_KEGG_results <- Cluster1_KEGG_results %>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

#reducing the size of the file
Cluster1_KEGG_results <- Cluster1_KEGG_results %>% dplyr::select(ID, Description.y, Enrichment, CombinedScore, GeneRatio,   BgRatio, RichFactor, FoldEnrichment,zScore, pvalue, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust)

#call data
Cluster2 <- read_csv("../Cluster_2_Data.csv")

Cluster2_annotated <- Cluster2 %>%
  left_join(annotation, by = c("GeneID" = "TrinityID"))

# Pull UNIQUE transcript IDs Cluster2 data frame
Cluster2_ID <- Cluster2_annotated %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool
Cluster2_KEGG_results<- enricher(
  gene         = Cluster2_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = ko_term2gene,           
  pvalueCutoff = 0.05,                       
  pAdjustMethod = "BH",
  minGSSize = 3)

# Convert to data frame
Cluster2_KEGG_results <- as.data.frame(Cluster2_KEGG_results)
head(Cluster2_KEGG_results)

# Add descriptions using the verified clean key
Cluster2_KEGG_results <- Cluster2_KEGG_results %>%
  left_join(kegg_descriptions, by = c("ID" = "KEGG_ID")) %>%
  dplyr::select(ID, Description.y, everything())

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster2_KEGG_results <- Cluster2_KEGG_results %>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

#reducing the size of the file
Cluster2_KEGG_results <- Cluster2_KEGG_results %>% dplyr::select(ID, Description.y, Enrichment, CombinedScore, Count, GeneRatio, BgRatio, RichFactor, FoldEnrichment,zScore, pvalue, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust)

Cluster2_KEGG_results <- Cluster2_KEGG_results %>%
  mutate(Description.y = if_else(ID == "ko01130","Biosynthesis of antibiotics",Description.y))

# run enrichment and combine
Cluster3 <- read_csv("../Cluster_3_Data.csv")

Cluster3_annotated <- Cluster3 %>%
  left_join(annotation, by = c("GeneID" = "TrinityID"))

# Pull unique transcript IDs from Cluster3 data frame
Cluster3_ID <- Cluster3_annotated %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool
Cluster3_KEGG_results<- enricher(
  gene         = Cluster3_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = ko_term2gene,           
  pvalueCutoff = 0.05,                       
  pAdjustMethod = "BH",
  minGSSize = 3)

# Convert to data frame
Cluster3_KEGG_results <- as.data.frame(Cluster3_KEGG_results)
head(Cluster3_KEGG_results)

# Add descriptions using the verified clean key
Cluster3_KEGG_results <- Cluster3_KEGG_results %>%
  left_join(kegg_descriptions, by = c("ID" = "KEGG_ID")) %>%
  dplyr::select(ID, Description.y, everything())

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster3_KEGG_results <- Cluster3_KEGG_results %>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

Cluster3_KEGG_results <- Cluster3_KEGG_results %>%
  mutate(Count = ifelse(row_number() == 1, 5, Count))

#reducing the size of the file
Cluster3_KEGG_results <- Cluster3_KEGG_results %>% dplyr::select(ID, Description.y, Enrichment, CombinedScore, Count,GeneRatio,   BgRatio, RichFactor, FoldEnrichment,zScore, pvalue, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust)

#  run enrichment and combine
Cluster4 <- read_csv("../Cluster_4_Data.csv")

Cluster4_annotated <- Cluster4 %>%
  left_join(annotation, by = c("GeneID" = "TrinityID"))

# Pull unique transcript IDs from Cluster4 data frame
Cluster4_ID <- Cluster4_annotated %>% 
  distinct(GeneID) %>% 
  pull(GeneID)

# Run the enrichment tool
Cluster4_KEGG_results<- enricher(
  gene         = Cluster4_ID,         
  universe     = global_background_universe, 
  TERM2GENE    = ko_term2gene,           
  pvalueCutoff = 0.05,                       
  pAdjustMethod = "BH",
  minGSSize = 3)

# Convert to data frame
Cluster4_KEGG_results <- as.data.frame(Cluster4_KEGG_results)
head(Cluster4_KEGG_results)

# Add descriptions using the verified clean key
Cluster4_KEGG_results <- Cluster4_KEGG_results %>%
  left_join(kegg_descriptions, by = c("ID" = "KEGG_ID")) %>%
  dplyr::select(ID, Description.y, everything())

#Turn the formal results object into a standard table / Calculate the score by multiplying statistical weight by effect size
Cluster4_KEGG_results <- Cluster4_KEGG_results %>%
  mutate(CombinedScore = -log(pvalue) * zScore) %>%
  arrange(desc(CombinedScore))

#reducing the size of the file
Cluster4_KEGG_results <- Cluster4_KEGG_results %>% dplyr::select(ID, Description.y, Enrichment, CombinedScore, Count, GeneRatio,   BgRatio, RichFactor, FoldEnrichment,zScore, pvalue, GeneRatio, BgRatio, FoldEnrichment, pvalue, p.adjust)

################################################################################################################################################################################################################################################################################################################################################
#selecting the top 23 combined score Cluster 2

#selecting the top 23 combined score go terms in each Revigo file
top20_KEGG_Cluster2 <- Cluster2_KEGG_results %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 23) %>%
  mutate(
    Description.y = factor(Description.y, levels = rev(Description.y)))

#selecting the top 20 combined score go terms in each Revigo file
top20_KEGG_Cluster3 <- Cluster3_KEGG_results %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 20) %>%
  mutate(Description.y = factor(Description.y, levels = rev(Description.y)))

#selecting the top 20 combined score go terms in each Revigo file
top20_KEGG_Cluster4 <- Cluster4_KEGG_results %>%
  arrange(desc(CombinedScore)) %>%
  slice_head(n = 20) %>%
  mutate(Description.y = factor(Description.y, levels = rev(Description.y)))

#combined into a single data frame
KEGG_Final <- bind_rows(
  mutate(top20_KEGG_Cluster2, Treatment = "Cluster #2"),
  mutate(top20_KEGG_Cluster3, Treatment = "Cluster #3"),
  mutate(top20_KEGG_Cluster4, Treatment = "Cluster #4"))

#Create an invisible unique ID column for plotting to keep the groups separate
KEGG_Final <- KEGG_Final %>%
  mutate(Plot_ID = row_number())

#Factor the Plot_ID so ggplot respects the correct sorted sequence
KEGG_Final$Plot_ID <- factor(KEGG_Final$Plot_ID, levels = unique(KEGG_Final$Plot_ID))

#Flip the levels for the highest values at the top of the plot
KEGG_Final$Plot_ID <- factor(KEGG_Final$Plot_ID, levels = rev(levels(KEGG_Final$Plot_ID)))

KEGG_Final$Treatment <- factor(KEGG_Final$Treatment, levels = c("Cluster #2", "Cluster #3", "Cluster #4"))

#lollipop graph fro KEGG_Final response
KEGG_lollipop <- ggplot(KEGG_Final, aes(y = Plot_ID, x = CombinedScore, color = Treatment)) +
  geom_segment(aes(x = 0, xend = CombinedScore, y = Plot_ID, yend = Plot_ID), 
               linewidth = 1) +
  geom_point(aes(size = Count)) +
  scale_y_discrete(labels = setNames(as.character(KEGG_Final$Description.y), KEGG_Final$Plot_ID)) +
  scale_color_manual(values = c(
    "Cluster #2" = "#d73027",
    "Cluster #3" = "#4575b4",
    "Cluster #4" = "purple")) +
  scale_size_continuous(range = c(1, 10)) + 
  labs(
    y = "KEGG Pathway",
    x = "Combined Score",
    color = "Dataset",
    size = "Transcripts") +
  theme_classic(base_size = 20) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title = element_text(size = 26, face = "bold"),
    legend.title = element_text(face = "bold", size = 26),
    legend.text  = element_text(face = "bold", size = 22))
KEGG_lollipop

# Save the plot
ggsave(filename = "KEGG_Lollipop.pdf", plot = KEGG_lollipop, dpi = 3000, scale = 2)



