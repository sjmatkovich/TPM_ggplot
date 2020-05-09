# Scot Matkovich
# May 2020

## Input data
count_table <- read.table(file="GSE116250_featureCounts_matrix.txt", row.names=1) # integer values assigned to each ENSG gene for each sample
types <- read.table(file="GSE116250_unpaired_design-disease-gender-age-ethn.txt", row.names=1, header=TRUE, stringsAsFactors = FALSE) # sample grouping by disease and factors used in voom modeling
lengths <- read.csv(file="GRCh38-94_ENSG_lengths.csv", row.names=1) # exon lengths used for conversion to TPM
voom_log2_cpm <- read.table(file="GSE116250_disease-age-gender_TMM_voom_log2cpm_table.txt", row.names=1) # output from limma-voom linear modeling adjustment of input data, v$E
out.file.prefix <- c("GSE116250_disease-age-gender")

## Libraries (tidyverse)
library(tidyverse)

## Produce TPM from featureCounts

      # function from https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/ and https://www.biostars.org/p/307603/
      counts_to_tpm <- function(counts, featureLength) {
        # Ensure valid arguments.
        stopifnot(nrow(featureLength) == nrow(counts))
        # Compute effective lengths of features in each library.
        effLen <- featureLength
        # Process one column at a time.
        tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
          rate = log(counts[,i]) - log(effLen)
          denom = log(sum(exp(rate), na.rm = TRUE))
          exp(rate - denom + log(1e6))
        }))
        # Copy the row and column names from the original matrix.
        colnames(tpm) <- colnames(counts)
        rownames(tpm) <- rownames(counts)
        return(tpm)
      }
      
      tpm_matrix <- counts_to_tpm(count_table, lengths)

      # transpose matrix of TPM
      trans.tpm <- as.data.frame(t(tpm_matrix))
      sample.names <- rownames(trans.tpm)
      trans.tpm <- cbind(sample.names, trans.tpm)
      rm(sample.names)

## import metadata incl group designations and merge tables
      samples <- rownames(types)
      samples.disease <- cbind(samples, disease=types$disease)
      trans.tpm.types <- merge(samples.disease, trans.tpm, by.x="samples", by.y="sample.names")
      
## write out unadjusted TPM matrix
      write.csv(trans.tpm.types, file=paste0(out.file.prefix,"_unadjusted_TPM.csv"))

## obtain mean TPM for each ENSG in each disease group
      # trans.tpm.groups.split <- split(select(trans.tpm.types, -c("samples","disease")), trans.tpm.types$disease)
      # tpm.group.means <- lapply(trans.tpm.groups.split, colMeans)
      
## process voom_log2cpm table so that TPM can be adjusted to correspond to linear modelling statistical results
      trans.voom.cpm <- as.data.frame(2^(t(voom_log2_cpm)))
      sample.names.voom <- rownames(trans.voom.cpm)
      trans.voom.cpm <- cbind(sample.names.voom, trans.voom.cpm)
      rm(sample.names.voom)
      trans.voom.types <- merge(samples.disease, trans.voom.cpm, by.x="samples", by.y="sample.names.voom")
      trans.voom.groups.split <- split(select(trans.voom.types, -c("samples","disease")), trans.voom.types$disease)
      voom.group.means <- lapply(trans.voom.groups.split, colMeans)
      trans.voom.adj <- mapply('/', select(trans.voom.types,-c("samples","disease")), voom.group.means[["NF"]]) # select the correct group mean for normalizing data, generally nonHF, NF, control etc.
      trans.voom.adj <- cbind(samples=trans.voom.types$sample, disease=trans.voom.types$disease, as.data.frame(trans.voom.adj))
      
      # ensure trans.tpm.types and trans.voom.adj have ENSG entries (columns) in the same order. Since they were prepared by merging with the same samples.disease vector, the rows are already in the same order
      trans.tpm.types.select <- select(trans.tpm.types, -c("samples","disease"))
      trans.voom.adj.select <- select(trans.voom.adj, -c("samples","disease"))
      
## write out adj.tpm matrix
      adj.tpm <- trans.tpm.types.select[, order(names(trans.tpm.types.select))] * trans.voom.adj.select[, order(names(trans.voom.adj.select))]
      rm(trans.tpm.types.select)
      rm(trans.voom.adj.select)
      adj.tpm <- cbind(samples=trans.tpm.types$samples, disease=trans.tpm.types$disease, adj.tpm)
      write.csv(adj.tpm, file=paste0(out.file.prefix,"_voom_TPM.csv"))
      
## template for ggplot boxplot loop
genes <- c("ENSG00000248527","ENSG00000228794","ENSG00000188976","ENSG00000188290")
for (g in seq_along(genes)) {
  # print(ggplot(trans.tpm.types, aes(disease, y=.data[[genes[g]]])) + geom_boxplot())
  bx <- ggplot(trans.tpm.types, aes(disease, y=.data[[genes[g]]])) + geom_boxplot()
  # make sure that y axis starts at x=0; change theme for ggplot; add title; retrieve fc and FDR from a lookup table and add to plot; etc.
  ggsave(bx, file=paste0(genes[g],".png"))
}

