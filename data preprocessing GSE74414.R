#setting up the working directory
Working_directory = "/home/student2024/Desktop/Project_Depression/GSE74414" # Replace with an appropriate directory
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

#reading files
GSE74414_methyl = smart_fread("/home/student2024/Desktop/Project_Depression/GSE74414/GSE74414_all_raw_methylated_signals.txt")
GSE74414_non_methyl = smart_fread("/home/student2024/Desktop/Project_Depression/GSE74414/GSE74414_all_raw_unmethylated_signals.txt")
GSE74414_det_p = smart_fread("/home/student2024/Desktop/Project_Depression/GSE74414/GSE74414_all_detection.p.values.txt")

# Resetting row names
rownames(GSE74414_methyl) = GSE74414_methyl$ID_REF
rownames(GSE74414_non_methyl) = GSE74414_non_methyl$ID_REF
rownames(GSE74414_det_p) = GSE74414_det_p$ID_REF
GSE74414_methyl$ID_REF = NULL
GSE74414_non_methyl$ID_REF = NULL
GSE74414_det_p$ID_REF = NULL
all(rownames(GSE74414_methyl) == rownames(GSE74414_non_methyl)) # TRUE
all(rownames(GSE74414_methyl) == rownames(GSE74414_det_p))  # TRUE
all(colnames(GSE74414_methyl) == colnames(GSE74414_non_methyl)) # TRUE
all(colnames(GSE74414_methyl) == colnames(GSE74414_det_p))  # TRUE

# Remove first rows
GSE74414_methyl = GSE74414_methyl[-1,]
GSE74414_non_methyl = GSE74414_non_methyl[-1,]
GSE74414_det_p = GSE74414_det_p[-1,]
ref_rownames = rownames(GSE74414_methyl)
ref_colnames = colnames(GSE74414_methyl)

# Set as numeric
GSE74414_methyl = as.data.frame(apply(GSE74414_methyl, 2, as.numeric))
GSE74414_non_methyl = as.data.frame(apply(GSE74414_non_methyl, 2, as.numeric))
GSE74414_det_p = as.data.frame(apply(GSE74414_det_p, 2, as.numeric))
rownames(GSE74414_methyl) = ref_rownames
rownames(GSE74414_non_methyl) = ref_rownames
rownames(GSE74414_det_p) = ref_rownames

colnames(GSE74414_methyl) = paste0(colnames(GSE74414_methyl), " Methylated signal")
colnames(GSE74414_non_methyl) = paste0(colnames(GSE74414_non_methyl), " Unmethylated signal")
merged_signals = cbind(GSE74414_non_methyl, GSE74414_methyl)
fwrite(merged_signals, "merged_signals.csv", row.names = TRUE, sep = ",")

# Making RAW GenomicMethylSet 
GSE74414_raw = readGEORawFile(filename = "merged_signals.csv",
                              Uname = "Unmethylated signal",
                              Mname = "Methylated signal",
                              sep = ",")

Beta_from_signals_raw = getBeta(GSE74414_raw, offset = 100)
fwrite(Beta_from_signals_raw, "Beta_from_signals_raw_GSE74414.csv", row.names = TRUE, sep = ",")

# Check
all(rownames(Beta_from_signals_raw) == rownames(GSE74414_methyl)) #False #not matching

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(Beta_from_signals_raw),]
all(rownames(Beta_from_signals_raw) == rownames(annot)) #True

# Need to obtain detection P-vals
pvals = GSE74414_det_p
pvals = pvals[rownames(Beta_from_signals_raw),]

# Check
all(colnames(Beta_from_signals_raw) == colnames(pvals)) # TRUE
all(rownames(Beta_from_signals_raw) == rownames(pvals)) # TRUE

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

write.csv(betas, "blood_beta_GSE74414.csv", row.names = TRUE, sep = ";")
write.csv(m.val, "M_values_GSE74414.csv", row.names = TRUE, sep = ";") 




