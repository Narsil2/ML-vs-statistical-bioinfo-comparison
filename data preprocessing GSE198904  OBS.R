Working_directory = "/home/student2024/Desktop/Project_Depression/GSE198904"# Replace with an appropriate directory
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

# Reading files
GSE198904_OBS_signals = smart_fread("/home/student2024/Desktop/Project_Depression/GSE198904/GSE198904_OBSERVEMDD0001_Matrix_signal_intensities.txt")
GSE198904_OBS_p = GSE198904_OBS_signals[,seq(from = 1, to = ncol(GSE198904_OBS_signals)-2, by=3)]
fwrite(GSE198904_OBS_signals, "GSE198904_OBS_signals.csv", row.names = TRUE, sep = ";")
rm(GSE198904_OBS_signals)
gc()

# Raw MethylSet
GSE198904_OBS_raw = readGEORawFile(filename = "GSE198904_OBS_signals.csv",
                                    Uname = "Unmethylated_Signal",
                                    Mname = "Methylated_Signal",
                                    array = "IlluminaHumanMethylationEPIC",
                                    annotation = "ilm10b2.hg19",
                                    sep = ";")

Beta_from_signals_raw = getBeta(GSE198904_OBS_raw, offset = 100)
fwrite(as.data.frame(Beta_from_signals_raw), "Beta_from_signals_raw_GSE198904_OBS.csv", row.names = TRUE, sep = ",")
rm(GSE198904_OBS_raw)
gc()

# Getting annotation and matching the rows
annot = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annot = as.data.frame(annot)
annot$Type = ifelse(annot$Type == "I",1,2)
annot = annot[rownames(Beta_from_signals_raw),]
all(rownames(Beta_from_signals_raw) == rownames(annot)) #True
gc()

# Need to obtain detection P-vals
pvals = GSE198904_OBS_p
pvals = pvals[rownames(Beta_from_signals_raw),]

# Check
all(colnames(Beta_from_signals_raw) == colnames(pvals)) # TRUE
all(rownames(Beta_from_signals_raw) == rownames(pvals)) # TRUE
rm(GSE198904_OBS_p)
gc()

#filtering samples to match the normalised data
cols_to_select_samples=colSums(pvals<=0.00005)>=0.75*nrow(pvals) #Select samples with more than 75% of probes with detection P-value less than 0.00005

# Remove sites on X and Y chromosomes
indexes_to_select_X=which(annot[,'chr']!='chrX')
indexes_to_select_Y=which(annot[,'chr']!='chrY')
indexes_to_select_XY=intersect(indexes_to_select_X,indexes_to_select_Y)
lines_to_select_XY=rownames(annot[indexes_to_select_XY,])

# Remove probes that do not target CpGs (starting with "ch")
lines_to_select_ch=rownames(annot)[which(substr(rownames(annot),1,2)=='cg')]

lines_to_select = intersect(lines_to_select_XY, lines_to_select_ch)

betas = Beta_from_signals_raw[lines_to_select, cols_to_select_samples]

betas = na.omit(betas)
m.val = beta2m(betas)

write.csv(betas, "betas_OBS.csv", row.names = TRUE, sep = ";")
gc()
write.csv(m.val, "Mval_OBS.csv", row.names = TRUE, sep = ";")
gc()
