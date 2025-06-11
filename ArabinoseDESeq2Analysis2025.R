#Analysis of Arabinose RNAseq Data 2025; LMR, based on IM and JK workflow

# install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# change the library path as needed, to the first value returned by typing
#   .libPaths()
#   in the console
my_lib_path = .libPaths()[1]

# install packages
BiocManager::install(c('Rsubread'), lib=my_lib_path)
BiocManager::install(c('DESeq2'), lib=my_lib_path)
BiocManager::install(c('tximport'), lib=my_lib_path)
BiocManager::install(c("rhdf5"), lib=my_lib_path)
BiocManager::install(c("BSgenome"), lib=my_lib_path)
BiocManager::install(c("GenomicFeatures"), lib=my_lib_path)
BiocManager::install(c('apeglm'), lib=my_lib_path)
BiocManager::install(c('EnhancedVolcano'), lib=my_lib_path)
BiocManager::install(c("Matrix"), lib=my_lib_path)
BiocManager::install(c('ensembldb'), lib=my_lib_path)
BiocManager::install(c("AnnotationDbi"), lib=my_lib_path)
BiocManager::install(c("org.EcK12.eg.db"), lib=my_lib_path)
install.packages(c("readxl"), lib=my_lib_path)
install.packages(c("openxlsx"), lib=my_lib_path)
install.packages(c("stringr"), lib=my_lib_path)


# load packages
library(GenomicFeatures)
library(Rsubread)
library(DESeq2)
library(tximport)
library(rhdf5)
library(EnhancedVolcano)
library(apeglm)
library(Matrix)
library(ensembldb)
library(AnnotationDbi)
library(org.EcK12.eg.db)
library(readxl)
library(openxlsx)
library(stringr)



# create list of quantification file names
# your working directory should be in the directory with the LMR-P1-28-0 etc folders
# can set this in the console with setwd("<dirname>")
# there should be TWO of these potential working directories, one for planktonic and one for biofilm, run them separately
# your samples.txt,
# these .txt files need to be filled manually by you :)
# have a Filename column (ex LMR-P1-28-0/quant/quant.sf), a Sample column
#    (P1-28), and a Sugar column (treated vs untreated) 
samples <- read.table(file="samples.txt", header = TRUE)
filenames <- paste(samples$Filename)
files <- file.path(filenames)
names(files) <- paste0(1:6) # NOTE: change this between 6 & 24 depending on how many samples you're working with

#Check to make sure your files match the filenames in samples.txt
files
all(file.exists(files))

# create tx2gene 
# can get the gff3.gz file (renamed) from:
# https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-56/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.49.gff3.gz
txdb <- makeTxDbFromGFF("EColi_k12.gff3.gz")
k <- keys( txdb, keytype = "TXNAME" )
tx2gene <- select(txdb, k, "GENEID", "TXNAME") # this is "EcoGene ID, accession number"; change as you want to!
head(tx2gene)

#import quantification data to DESEQ2

txi_all <- tximport(files = files, type = "salmon", tx2gene = tx2gene)
names(txi_all)
head(txi_all$counts)

#check comparability
all(rownames(samples) %in% colnames(txi_all$counts))

#Get necessary columns from sample info files
# Sample is the column containing the descriptors, like P1-28, P2-37, etc
# Sugar is the column detailing whether the sample has the sugar or not.
#   Should be either "treated" or "untreated", or you'll run into problems a lil further down.

all_info <- samples[c("Sample", "BorP", "Temp","Sugar")]

# make DESEQ2 objects
dds_all <- DESeqDataSetFromTximport(txi_all, colData=all_info, design=~Sugar)

#pre filtering out low count obs
filt <- rowSums(counts(dds_all)) >= 10
dds_all <- dds_all[filt,]

# set factor level - this is detailing where the baseline is
# this is where it's important that you have the Sugar column as "treated" and "untreated"
# if you want to change that above, just make sure you change it here too
dds_all$Sugar <- relevel(dds_all$Sugar, ref = "untreated")

# DESeq
dds_all <- DESeq(dds_all)

# can change the alpha value as desired
res_all <- results(dds_all, alpha = 0.05)

# build the lookup_table of EcoGene IDs : Official Gene Symbols
gff <- read.delim("EColi_k12.gff3", header=FALSE, comment.char="#") # read in gff3 file
gff_length <- nrow(gff)
lookup_table <- data.frame(EcoGeneID = character(0), GeneSymbol = character(0)) # data frame to be filled
gff_line_no <- 1
while (gff_line_no < gff_length) {
  gff_type <- gff$V3[gff_line_no]
  gff_info <- gff$V9[gff_line_no]
  if (gff_type == "gene") {
    ecogene_id <- str_extract(gff_info, "(?<=gene_id=)[^;]+")
    gene_symbol <- str_extract(gff_info, "(?<=Name=)[^;]+")
    lookup_table <- rbind(lookup_table, data.frame(EcoGeneID = ecogene_id, GeneSymbol = gene_symbol))
  }
  gff_line_no <- gff_line_no + 1
}

# convert all EcoGeneIDs (b0001, b0002, etc) to Official Gene Symbols using the lookup table
b_numbers_in_res <- rownames(res_all)
res_length <- length(b_numbers_in_res)
iter <- 1
while (iter <= res_length) {
  if (b_numbers_in_res[iter] %in% lookup_table$EcoGeneID) {
    b_number <- subset(lookup_table, EcoGeneID == b_numbers_in_res[iter], select = GeneSymbol)[1,]
    rownames(res_all)[iter] <- b_number
  }
  print(iter)
  iter <- iter + 1
}

# view!
View(data.frame(res_all))

# Load necessary library
library(openxlsx)

# Assuming your DESeq2 results are in a data frame named `res_all`
# Convert DESeq2 results to a data frame if it isn't already
res_all_df <- as.data.frame(res_all)

# Add gene names as a column (assuming rownames are gene names)
res_all_df$Gene <- rownames(res_all_df)

# Reorder columns to move Gene to the first position
res_all_df <- res_all_df[, c("Gene", "baseMean", "log2FoldChange", "lfcSE", 
                             "stat", "pvalue", "padj")]

# Filter for significant rows: padj < 0.01 and abs(log2FoldChange) > 1
significant_genes <- res_all_df[res_all_df$padj < 0.01 & 
                                  abs(res_all_df$log2FoldChange) > 1, ]

# Further separate significant genes into upregulated and downregulated
upregulated_genes <- significant_genes[significant_genes$log2FoldChange > 2, ]
downregulated_genes <- significant_genes[significant_genes$log2FoldChange < -2, ]

# Create a new workbook
wb <- createWorkbook()

# Add the first sheet with all data
addWorksheet(wb, "All Data")
writeData(wb, "All Data", res_all_df)

# Add the second sheet with only significant genes
addWorksheet(wb, "Significant Genes")
writeData(wb, "Significant Genes", significant_genes)

# Add the third sheet with only upregulated genes
addWorksheet(wb, "Upregulated Genes")
writeData(wb, "Upregulated Genes", upregulated_genes)

# Add the fourth sheet with only downregulated genes
addWorksheet(wb, "Downregulated Genes")
writeData(wb, "Downregulated Genes", downregulated_genes)

# Save the workbook to an Excel file
saveWorkbook(wb, "DESeq2_Results_Extended.xlsx", overwrite = TRUE)

# Optional: Print a message to confirm success
cat("Excel file 'DESeq2_Results_Extended.xlsx' created successfully.\n")


# print out these unfiltered data frames to file, ordered by FC
unfiltered_by_FC <- res_all[order(res_all$log2FoldChange),]
write.xlsx(data.frame(unfiltered_by_FC), file = "unfiltered_by_FC.xlsx", rowNames=TRUE)
