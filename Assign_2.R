#Assignemnt 1 - Software Tools
#By Veedhi Solanki, last updated on October 26, 2023

#Theme: Unsupervised Machine Learning or Clustering of Sequences

#We are trying to figure out the gene clustering patterns within the genus Drosophila.Here I am predicting two contrasting hypothesis which will be proven true or rejected through the analysis here via sequnce alignment and clustering. There are contrasting information available in literature that support each of these hypothesis more or less. So, let's figure which hypothesis is supported by the data obtained from NCBI and our analysis on it. 

#Hypothesis (H1): The clustering patterns of gene sequences within the Drosophila genus will exhibit significant similarity, indicating a conserved evolutionary pressure and functional constraints across genes, thus suggesting their suitability as reliable markers for classification applications.

#Hypothesis (H2): The clustering patterns of gene sequences within the Drosophila genus will demonstrate substantial variability, suggesting diverse evolutionary pressures and functional roles among different genes. This variation may imply the existence of specific gene clusters more suitable as markers for classifying distinct subgroups within the genus.

#These contrasting hypothesis, provide an opportunity for intriguing exploration of the clustering patterns of different genes within the Drosophila genus. This approach allows to investigate the extent of similarity or dissimilarity between the genes, providing insights into their evolutionary relationships and potential implications for classification applications.

###Loading packages:
library(rentrez)
library(tidyverse)
library(Biostrings)
library(muscle)
library(ape)
library(DECIPHER)
library(phyloseq)
library(seqinr)
library(stringr)
library(RSQLite)
library(dplyr)
library(gplots)
library(Rtsne)
library(ggplot2)
library(cluster)
library(RColorBrewer)


#I am not labeling specifically Code section 1 as I am doing quality control and data exploration even throughout the script and not just at the beginning. 

###Obtaining data for genes COI and Adh:

entrez_db_searchable("nuccore")
#I am restricting gene lenths around theor actual length
search_term_COI <- "(Drosophila[Organism] AND COI[Gene] AND 600:700[SLEN])"
search_term_Adh <- "(Drosophila[Organism] AND Adh[Gene]AND 600:700[SLEN])"

search_COI <- entrez_search(db = "nuccore", term = search_term_COI, retmax = 100)
search_COI #3858 hits
search_Adh <- entrez_search("nuccore", term = search_term_Adh, retmax = 100)
search_Adh #43 hits

#Exploring both datasets
search_COI$ids
search_Adh$ids
class(search_COI$ids)
class(search_Adh$ids)
length(search_COI$ids)
length(search_Adh$ids)
maxHits_COI <- search_COI$count
maxHits_Adh <- search_Adh$count

#Set retmax argument to maxHits
search_COI_2 <- entrez_search(db = "nuccore", term = search_term_COI, retmax = maxHits_COI, use_history = T)
search_Adh_2 <- entrez_search("nuccore", term = search_term_Adh, retmax = maxHits_Adh, use_history = T)

#Exploting new search dataset
search_COI_2$ids
search_COI_2$retmax
search_COI_2$count
search_Adh_2$ids
search_Adh_2$retmax
search_Adh_2$count
length(search_COI_2$ids)
length(search_Adh_2$ids)

#Getting data with maxhits
# Fetch the sequences in smaller batches because the data is huge
batch_size <- 50

COI_batches <- list()
for (i in seq(0, maxHits_COI, by = batch_size)) {
  COI_batch <- entrez_fetch(db = "nuccore", web_history = search_COI_2$web_history, rettype = "fasta", retstart = i, retmax = batch_size)
  if (length(COI_batch) > 0) {
    COI_batches <- c(COI_batches, list(COI_batch))
  }
}

Adh_batches <- list()
for (i in seq(0, maxHits_Adh, by = batch_size)) {
  Adh_batch <- entrez_fetch(db = "nuccore", web_history = search_Adh_2$web_history, rettype = "fasta", retstart = i, retmax = batch_size)
  if (length(Adh_batch) > 0) {
    Adh_batches <- c(Adh_batches, list(Adh_batch))
  }
}

#Merge the batches
COI_sequences <- unlist(COI_batches)
Adh_sequences <- unlist(Adh_batches)

#Exploring the sequnces of both genes
class(COI_sequences)
class(Adh_sequences)
head(COI_sequences)
head(Adh_sequences)

#saving sequnces in fasta file
write(COI_sequences, "COI_sequences.fasta", sep = "\n") 
write(Adh_sequences, "Adh_sequences.fasta", sep = "\n") 

stringSet_COI <- readDNAStringSet("COI_sequences.fasta")
class(stringSet_COI)
head(names(stringSet_COI))
stringSet_Adh <- readDNAStringSet("Adh_sequences.fasta")
class(stringSet_Adh)
head(names(stringSet_Adh))

#Converting the sequnces from NCBI to data frame for each gene
dfCOI <- data.frame(COI_Title = names(stringSet_COI), COI_Sequence = paste(stringSet_COI))
View(dfCOI)
dfAdh <- data.frame(Adh_Title = names(stringSet_Adh), Adh_Sequence = paste(stringSet_Adh))
View(dfAdh)
#Adding the coloumns species name to the datatframes
dfCOI$Species_Name <- word(dfCOI$COI_Title, 2L, 3L)
dfCOI <- dfCOI[, c("COI_Title", "Species_Name", "COI_Sequence")]
View(dfCOI)
dfAdh$Species_Name <- word(dfAdh$Adh_Title, 2L, 3L)
dfAdh <- dfAdh[, c("Adh_Title", "Species_Name", "Adh_Sequence")]
View(dfAdh)

#Exploring some components of the dataframe
length(dfCOI$COI_Sequence) #3858
length(dfAdh$Adh_Sequence) #43

#Filtering sequnces in both dataframes
missing.data <- 0.01
length.var <- 50
chosen.model <- "JC"
clustering.threshold <- 0.03
clustering.method <- "average"

dfCOI_1 <- dfCOI %>%
  filter(!is.na(COI_Sequence)) %>%
  mutate(COI_Sequence = str_remove_all(COI_Sequence, "^N+|N+$|-")) %>%
  filter(str_count(COI_Sequence, "N") <= (missing.data * str_count(COI_Sequence))) %>%
  filter(str_count(COI_Sequence) >= median(str_count(COI_Sequence)) - length.var & str_count(COI_Sequence) <= median(str_count(COI_Sequence)) + length.var)
dfCOI_1 <- separate(dfCOI_1, COI_Title, into = c("Accession", "Voucher", "Gene_Description", "Sequence_Description"), sep = " ", extra = "merge")

dfAdh_1 <- dfAdh %>%
  filter(!is.na(Adh_Sequence)) %>%
  mutate(Adh_Sequence = str_remove_all(Adh_Sequence, "^N+|N+$|-")) %>%
  filter(str_count(Adh_Sequence, "N") <= (missing.data * str_count(Adh_Sequence))) %>%
  filter(str_count(Adh_Sequence) >= median(str_count(Adh_Sequence)) - length.var & str_count(Adh_Sequence) <= median(str_count(Adh_Sequence)) + length.var)
dfAdh_1 <- separate(dfAdh_1, Adh_Title, into = c("Accession", "Voucher", "Gene_Description", "Sequence_Description"), sep = " ", extra = "merge")

#Getting to know some data sizes for for genes in dataframe
length(dfCOI_1$COI_Sequence) #3798
length(dfAdh_1$Adh_Sequence) #34
length(unique(dfCOI_1$Species_Name)) #139
length(unique(dfAdh_1$Species_Name)) #11

summary(nchar(dfCOI_1$COI_Sequence))
summary(nchar(dfAdh_1$Adh_Sequence))

#Histogram to display the distribution of sequence lengths.

my_colors <- brewer.pal(8, "Pastel2")

# Histogram for COI Sequence
hist(nchar(dfCOI_1$COI_Sequence), col = my_colors[4], 
     main = "Histogram of COI Sequence Lengths",
     xlab = "Sequence Length",
     ylab = "Frequency")

# Histogram for Adh Sequence
hist(nchar(dfAdh_1$Adh_Sequence), col = my_colors[4], 
     main = "Histogram of Adh Sequence Lengths",
     xlab = "Sequence Length",
     ylab = "Frequency")

#Exploring the sequences it self (unaligned)
class(dfCOI_1)
dfCOI_1$COI_Sequence <- DNAStringSet(dfCOI_1$COI_Sequence)
class(dfCOI_1$COI_Sequence)
names(dfCOI_1$COI_Sequence) <- dfCOI_1$Accession
names(dfCOI_1$COI_Sequence)
dfCOI_1$COI_Sequence
BrowseSeqs(dfCOI_1$COI_Sequence)

class(dfAdh_1)
dfAdh_1$Adh_Sequence <- DNAStringSet(dfAdh_1$Adh_Sequence)
class(dfAdh_1$Adh_Sequence)
names(dfAdh_1$Adh_Sequence) <- dfAdh_1$Accession
names(dfAdh_1$Adh_Sequence)
dfAdh_1$Adh_Sequence
BrowseSeqs(dfAdh_1$Adh_Sequence)

#Sequence alignment for Adh and COI:
dfCOI_1.alignment <- DNAStringSet(muscle::muscle(dfCOI_1$COI_Sequence, maxiters = 3, gapopen = -10000), use.names = TRUE)
dfCOI_1.alignment
BrowseSeqs(dfCOI_1.alignment)
mean(str_count(as.character(dfCOI_1.alignment), "-"))
length(dfCOI_1.alignment[[1]])

dfAdh_1.alignment <- DNAStringSet(muscle::muscle(dfAdh_1$Adh_Sequence, maxiters = 3, gapopen = -10000 ), use.names = TRUE)
dfCOI_1.alignment
BrowseSeqs(dfAdh_1.alignment)
mean(str_count(as.character(dfAdh_1.alignment), "-"))
length(dfAdh_1.alignment[[1]])

BrowseSeqs(dfCOI_1.alignment)
BrowseSeqs(dfAdh_1.alignment)

writeXStringSet(dfCOI_1.alignment, file = "COI_alignment.fasta")
writeXStringSet(dfCOI_1.alignment, file = "Adh_alignment.fasta")

###Calculating pairwise alignment matrix: First alignment without gaps and then calculating matrix
#Pairwise distance for COI: Calculation 
sequence_data_COI <- readDNAStringSet("COI_alignment.fasta")
sequence_data_COI_XStringSet <- as(sequence_data_COI, "XStringSet")
# Remove gap characters from the sequences
sequence_data_COI_no_gaps <- gsub("-", "", sequence_data_COI_XStringSet)
# Convert back to DNAStringSet
sequence_data_COI_no_gaps <- DNAStringSet(sequence_data_COI_no_gaps)
# Perform sequence alignment on the data without gaps
alignment_COI <- AlignSeqs(sequence_data_COI_no_gaps)
BrowseSeqs(alignment_COI)

# Convert the aligned sequences to a matrix
sequence_matrix_COI <- as.DNAbin(alignment_COI)
# Calculate pairwise distances
dist_mat_COI <- dist.dna(sequence_matrix_COI)
# Convert the distance matrix to a data frame for better visualization
dist_df_COI <- as.data.frame(as.matrix(dist_mat_COI))
# Convert the distance matrix to a data frame
dist_df_COI <- as.data.frame(as.matrix(dist_mat_COI))
rownames(dist_df_COI) <- colnames(dist_df_COI)

#Pairwise distance for Adh: Calculation 
sequence_data_Adh <- readDNAStringSet("Adh_alignment.fasta")
sequence_data_Adh_XStringSet <- as(sequence_data_Adh, "XStringSet")
sequence_data_Adh_no_gaps <- gsub("-", "", sequence_data_Adh_XStringSet)
sequence_data_Adh_no_gaps <- DNAStringSet(sequence_data_Adh_no_gaps)
alignment_Adh <- AlignSeqs(sequence_data_Adh_no_gaps)
sequence_matrix_Adh <- as.DNAbin(alignment_Adh)
dist_mat_Adh <- dist.dna(sequence_matrix_Adh)
dist_df_Adh <- as.data.frame(as.matrix(dist_mat_Adh))
rownames(dist_df_Adh) <- colnames(dist_df_Adh)

# Compute Sum of Squared Errors (SSE) for different values of k : COI
set.seed(123)  #since this method randomly choose one centre to start with, this code was essential for reproductibility of same graph
sse_COI <- c()
for (i in 2:10) {
  kmeans_out_COI <- kmeans(distanceMatrix_COI, centers = i)
  sse_COI[i] <- kmeans_out_COI$tot.withinss
}
# Plot the elbow plot : COI
plot(1:10, sse_COI, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters (k)", ylab = "Sum of Squared Errors", 
     main = "Elbow Method for Optimal k for COI")

# Compute Sum of Squared Errors (SSE) for different values of k : Adh
set.seed(123)  
sse_Adh <- c(10)
for (i in 2:10) {
  kmeans_out_Adh <- kmeans(distanceMatrix_Adh, centers = i)
  sse_Adh[i] <- kmeans_out_Adh$tot.withinss
}
# Plot the elbow plot : Adh
plot(1:10, sse_Adh, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters (k)", ylab = "Sum of Squared Errors", 
     main = "Elbow Method for Optimal k for Adh")

#################################################################
#Dendrograms were not used a lot in analysis except the result for Adh was compared to the one obtained from elbow diagram. This was since the data of COI is huge.
missing.data <- 0.01
length.var <- 50
chosen.model <- "JC"
clustering.threshold <- 0.05
clustering.method <- "single"

dnaBin.COI <- as.DNAbin(alignment_COI)
class(dnaBin.COI)
distanceMatrix_COI <- dist.dna(dnaBin.COI, model = chosen.model, as.matrix = TRUE, pairwise.deletion = TRUE)
head(distanceMatrix_COI)
clusters.COI <- DECIPHER::TreeLine(myDistMatrix = distanceMatrix_COI,
                                   method = clustering.method,
                                   cutoff = clustering.threshold,
                                   showPlot = TRUE,
                                   type = "both",
                                   verbose = TRUE)

class(clusters.COI)
clusters.COI
length(clusters.COI)
length(unique(unlist(clusters.COI[[1]][1])))

dnaBin.Adh <- as.DNAbin(dfAdh_1.alignment)
class(dnaBin.Adh)
distanceMatrix_Adh <- dist.dna(dnaBin.Adh, model = chosen.model, as.matrix = TRUE, pairwise.deletion = TRUE)
head(distanceMatrix_Adh)
clusters.Adh <- DECIPHER::TreeLine(myDistMatrix = distanceMatrix_Adh,
                                   method = "complete",
                                   cutoff = clustering.threshold,
                                   showPlot = TRUE,
                                   type = "both",
                                   verbose = TRUE)
class(clusters.Adh)
clusters.Adh
length(clusters.Adh)
length(unique(unlist(clusters.Adh[[1]][1])))
#################################################################
# Visualizing matrix: COI
unique_dist_mat_COI <- unique(as.matrix(dist_mat_COI))
# Run t-SNE on the unique matrix
tsne_result_COI <- Rtsne(unique_dist_mat_COI)
clusters_COI <- kmeans(tsne_result_COI$Y, centers = 5)$cluster
# Plot the t-SNE result 
plot_COI <- ggplot(data = as.data.frame(tsne_result_COI$Y), aes(x = V1, y = V2, color = factor(clusters_COI))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE Plot of the Drosophila gene COI") +
  scale_color_brewer(palette = "Pastel2") + 
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),  
        panel.grid.major = element_line(color = "gray90"),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = "black"),  
        axis.text = element_text(color = "black"),  
        legend.title = element_text(size = 12, color = "black"),  
        legend.text = element_text(size = 10, color = "black"), 
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),  
        text = element_text(color = "black")) 
plot_COI


#Visualizing matrix: Adh
unique_dist_mat_Adh <- unique(as.matrix(dist_mat_Adh))
tsne_result_Adh <- Rtsne(unique_dist_mat_Adh)
clusters_Adh <- kmeans(tsne_result_Adh$Y, centers = 3)$cluster
plot_Adh <- ggplot(data = as.data.frame(tsne_result_Adh$Y), aes(x = V1, y = V2, color = factor(clusters_Adh))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE Plot of the Drosophila gene Adh") +
  scale_color_brewer(palette = "Pastel2") +  
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white"),  
        panel.grid.major = element_line(color = "gray90"),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(color = "black"),  
        legend.title = element_text(size = 12, color = "black"),  
        legend.text = element_text(size = 10, color = "black"),  
        plot.title = element_text(hjust = 0.5, size = 16, color = "black"),  
        text = element_text(color = "black"))  
plot_Adh


# Compute Silhouette index : COI
COIdistanceMatrix <- as.matrix(distanceMatrix_COI)
COIhclust_object <- hclust(as.dist(COIdistanceMatrix), method = "complete")
k <- 5
COIclusters <- cutree(COIhclust_object, k = k)
sil_COI <- silhouette(COIclusters, dist = COIdistanceMatrix)
summary_sil_COI <- summary(sil_COI)
par(mar = c(5, 4, 4, 8))
plot(sil_COI, col = rainbow(k), border = NA, main = "Silhouette Plot: COI", cex.names = 0.6)


# Compute Silhouette index : Adh
AdhdistanceMatrix <- as.matrix(distanceMatrix_Adh)
Adhhclust_object <- hclust(as.dist(AdhdistanceMatrix), method = "complete")
k <- 3
Adhclusters <- cutree(Adhhclust_object, k = k)
sil_Adh <- silhouette(Adhclusters, dist = AdhdistanceMatrix)
summary_sil_Adh <- summary(sil_Adh)
par(mar = c(5, 4, 4, 8))
plot(sil_Adh, col = rainbow(k), border = NA, main = "Silhouette Plot: Adh", cex.names = 0.6)

#END of Script