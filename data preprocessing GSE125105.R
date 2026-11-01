#importing required libraries
library(stringi)
library(stringr)
library(minfi)
library(wateRmelon)
#GSE125105
#reading files
Working_directory = "/home/student2024/Desktop/Project_Depression/GSE125105" #Replace with an appropriate directory
setwd(Working_directory)

# Preparing target file
All_files = list.files("...")
All_files = All_files[stri_detect_fixed(All_files, pattern ="Red.idat")]
Full_file_name = stri_replace_all_fixed(All_files, pattern = "_Red.idat", replacement = "")
Participant = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
Sentrix.Barcode = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[2]
  return(x)
})
Sample.Section = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[3]
  return(x)
})
barcode = sapply(Full_file_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[-1]
  x = paste0(x, collapse = "_")
  return(x)
})
Targets_GSE125105 = data.frame(Full_file_name = Full_file_name, Participant = Participant, Sentrix.Barcode = Sentrix.Barcode, Sample.Section = Sample.Section, barcode = barcode)
rm(list = c("All_files", "Full_file_name","Participant", "Sentrix.Barcode", "Sample.Section", "barcode"))

# Reading IDAT files
IDAT_path = "/home/student2024/Desktop/Project_Depression/GSE125105"
Mset.raw = read.metharray.exp(base = IDAT_path, force=T, recursive = TRUE)


# Adjusting corresponding target file
Targets_GSE125105 = Targets_GSE125105[match(Mset.raw@colData@rownames, Targets_GSE125105$Full_file_name),]
Mset.raw@colData@rownames == Targets_GSE125105$Full_file_name

# Getting annotation
annot = as.data.frame(getAnnotation(Mset.raw))
write.csv(annot, "Annot_Illumina_450_K.csv")

# Background correction with dye-bias normalization ("noob" method)
Mset.A = preprocessNoob(Mset.raw)

# Background correction with Funnorm (may be required)
#Mset.B = preprocessFunnorm(Mset.raw)
#save(Mset.B, file = "Mset_Funnorm_normalized_bumphunt_combined_GSE125105.rda")

# Saving beta values after preprocessing with "noob"
Beta_Noob_Normalized = getBeta(Mset.A, offset = 100) # Offset of 100 is Illumina's default
#write.table(Beta_Noob_Normalized, "Beta_Noob_Normalized_combined_GSE125105.txt", sep = "\t")

# Obtaining beta value
#betQN = betaqn(Mset.A) # from the wateRmelon package
#betQN[betQN == 1] = 1 - min(betQN)

# Remove sites on X and Y chromosomes
#annot = as.data.frame(getAnnotation(Mset.raw))
indexes_to_select_X = which(annot[,'chr'] != 'chrX')
indexes_to_select_Y = which(annot[,'chr'] != 'chrY')
indexes_to_select_XY = intersect(indexes_to_select_X,indexes_to_select_Y)
lines_to_select_XY = rownames(annot[indexes_to_select_XY,])


# Remove probes that do not target CpGs (starting with "ch")
lines_to_select_ch = rownames(annot)[which(substr(rownames(annot), 1, 2) == 'cg')]

# Intersect all these probes that were selected, to only obtain the probes that pass all these filtering steps
lines_to_select = intersect(lines_to_select_XY,lines_to_select_ch)

betas = as.data.frame(Beta_Noob_Normalized)
pvals = detectionP(Mset.raw)
#samples have to match the normalized dataset
cols_to_select_samples = colSums(pvals <= 0.00005) >= 0.75 * nrow(pvals) #Select samples with more than 75% of probes with detection P-value less than 0.00005

betas = betas[lines_to_select,cols_to_select_samples]
betas = na.omit(betas)
M.val = beta2m(betas)

write.table(betas,file="adjusted_blood_beta_combined_GSE125105.txt", sep = "\t")
write.table(M.val,"M_values_adj_combined_GSE125105.txt", sep = "\t")
write.csv(Targets_GSE125105, "Targets_GSE125105.csv", row.names = FALSE)
