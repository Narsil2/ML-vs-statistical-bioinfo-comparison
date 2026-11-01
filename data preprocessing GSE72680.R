Working_directory = "/home/student2024/Desktop/Project_Depression/GSE72680"# Replace with an appropriate directory
setwd(Working_directory)

# Importing packages (not all of them are required)
library(stringr)
library(dplyr)
library(openxlsx)
library(grid)
library(gdata)
library(magrittr)
library(data.table)
library(XML)
library(RCurl)
library(stringi)
library(httr)
library(FlowSorted.Blood.450k)
library(minfi)
library(wateRmelon)
library(lumi)
library(lmtest)
library(sva)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(GEOquery)
library(meffil)

# A function to quickly read data
smart_fread = function(x){
  x = as.data.frame(fread(x, nThread = 14))
  rownames(x) = x$V1
  x$V1 = NULL
  return(x)
}

# Importing betas from signals
RAW_geo_Signals = readGEORawFile(filename = "/home/student2024/Desktop/Project_Depression/GSE72680/GSE72680_signals.txt", sep = "\t")
# Default Uname = "Unmethylated signal" and Mname = "Methylated signal" work

# Saving beta values after signals
Beta_from_signals_raw = getBeta(RAW_geo_Signals, offset = 100) # Offset of 100 is Illumina's default
#write.table(Beta_from_signals_raw, "Beta_from_signals_raw_GSE72680.txt", sep = "\t")

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(Beta_from_signals_raw),]
all(rownames(Beta_from_signals_raw) == rownames(annot)) #True

# Need to obtain detection P-vals
RAW_geo_Signals_data = smart_fread("/home/student2024/Desktop/Project_Depression/GSE72680/GSE72680_signals.txt")
pvals = RAW_geo_Signals_data[, seq(from=1, to=ncol(RAW_geo_Signals_data), by = 3)]
pvals = pvals[rownames(Beta_from_signals_raw),]

# Renaming columns for pvals
colnames(pvals) = stri_replace_all_fixed(colnames(pvals), pattern = "_Detection PVal", replacement = "")

# Check
all(rownames(pvals) == rownames(Beta_from_signals_raw)) # TRUE
all(colnames(pvals) == colnames(Beta_from_signals_raw)) # TRUE
all(rownames(Beta_from_signals_raw) == rownames(annot)) # TRUE

#filtering samples to match the normalised data
cols_to_select_samples = colSums(pvals <= 0.00005) >= 0.75 * nrow(pvals) # Select samples with more than 75% of probes with detection P-value less than 0.00005

# Remove sites on X and Y chromosomes
indexes_to_select_X = which(annot[,'chr'] != 'chrX')
indexes_to_select_Y = which(annot[,'chr'] != 'chrY')
indexes_to_select_XY = intersect(indexes_to_select_X,indexes_to_select_Y)
lines_to_select_XY = rownames(annot[indexes_to_select_XY,])

# Remove probes that do not target CpGs (starting with "ch")
lines_to_select_ch = rownames(annot)[which(substr(rownames(annot), 1, 2) == 'cg')]

lines_to_select = intersect(lines_to_select_XY, lines_to_select_ch)

betas = Beta_from_signals_raw[lines_to_select, cols_to_select_samples]

betas = na.omit(betas)
m.val = beta2m(betas)


# Writing output
write.table(betas, file="blood_beta_combined_GSE72680.csv", row.names = TRUE)
write.table(m.val, file="M_values_combined_GSE72680.csv", row.names = TRUE) 
