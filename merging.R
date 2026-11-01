#importing libraries
library(stringr)
library(dplyr)
library(data.table)
library(stringi)

# A function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 10, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}

#importing data
# GSE125105_MPIP
GSE125105_MPIP_mval = smart_fread('/home/student2024/Desktop/Project_Depression/GSE125105/M_values_adj_combined_GSE125105.txt')
GSE125105_MPIP_beta = smart_fread('/home/student2024/Desktop/Project_Depression/GSE125105/adjusted_blood_beta_combined_GSE125105.txt')

# GSE72680_GRADY
GSE72680_GRADY_mval = smart_fread('/home/student2024/Desktop/Project_Depression/GSE72680/M_values_combined_GSE72680.csv')
GSE72680_GRADY_beta = smart_fread('/home/student2024/Desktop/Project_Depression/GSE72680/blood_beta_combined_GSE72680.csv')

# GSE113725_RDE
GSE113725_RDE_mval = smart_fread('/home/student2024/Desktop/Project_Depression/GSE113725/M_values_adj_GSE113725_full.txt')
GSE113725_RDE_beta = smart_fread('/home/student2024/Desktop/Project_Depression/GSE113725/adjusted_blood_beta_GSE113725_full.txt')

# GSE198904_DHRC
GSE198904_DHRC_mval = smart_fread('/home/student2024/Desktop/Project_Depression/GSE198904/Mval_DHRC.csv')
GSE198904_DHRC_beta = smart_fread('/home/student2024/Desktop/Project_Depression/GSE198904/betas_DHRC.csv')

# GSE198904_OBS
GSE198904_OBS_mval = smart_fread('/home/student2024/Desktop/Project_Depression/GSE198904/Mval_OBS.csv')
GSE198904_OBS_beta = smart_fread('/home/student2024/Desktop/Project_Depression/GSE198904/betas_OBS.csv')

# GSE74414_MPIP2
GSE74414_MPIP2_mval = smart_fread('/home/student2024/Desktop/Project_Depression/GSE74414/M_values_GSE74414.csv')
GSE74414_MPIP2_beta = smart_fread('/home/student2024/Desktop/Project_Depression/GSE74414/blood_beta_GSE74414.csv')

combined_pheno = as.data.frame(read.csv('/home/student2024/Desktop/Project_Depression/combined_dataset_phenotypes.csv'))

#inspecting overlapping methylation
CpGs_GSE125105_MPIP_m = rownames(GSE125105_MPIP_mval)
CpGs_GSE72680_GRADY_m = rownames(GSE72680_GRADY_mval)
CpGs_GSE113725_RDE_m = rownames(GSE113725_RDE_mval)
CpGs_GSE198904_DHRC_m = rownames(GSE198904_DHRC_mval)
CpGs_GSE198904_OBS_m = rownames(GSE198904_OBS_mval)
CpGs_GSE74414_MPIP2_m = rownames(GSE74414_MPIP2_mval)

CpG_list = list(
  "GSE125105 MPIP" = CpGs_GSE125105_MPIP_m,
  "GSE72680 GRADY" = CpGs_GSE72680_GRADY_m,
  "GSE113725 RDE" = CpGs_GSE113725_RDE_m,
  "GSE198904 DHRC" = CpGs_GSE198904_DHRC_m,
  "GSE198904 OBS" = CpGs_GSE198904_OBS_m,
  "GSE74414 MPIP2" = CpGs_GSE74414_MPIP2_m
)

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

Intersecting_CpGs = Reduce(intersect, CpG_list)
sapply(CpG_list, function(x) all(Intersecting_CpGs %in% x)) # All -> TRUE

Combined_Mval = list(
  GSE125105_MPIP_mval[Intersecting_CpGs,],
  GSE72680_GRADY_mval[Intersecting_CpGs,],
  GSE113725_RDE_mval[Intersecting_CpGs,],
  GSE198904_DHRC_mval[Intersecting_CpGs,],
  GSE198904_OBS_mval[Intersecting_CpGs,],
  GSE74414_MPIP2_mval[Intersecting_CpGs,]
)

# Check
sapply(Combined_Mval, function(x) all(rownames(x) == Intersecting_CpGs)) # All -> TRUE

Combined_Mval = do.call(cbind, Combined_Mval)
Combined_Mval = Combined_Mval[,combined_pheno$ID]

#repeating for beta
#inspecting overlapping methylation
CpGs_GSE125105_MPIP_b = rownames(GSE125105_MPIP_beta)
CpGs_GSE72680_GRADY_b = rownames(GSE72680_GRADY_beta)
CpGs_GSE113725_RDE_b = rownames(GSE113725_RDE_beta)
CpGs_GSE198904_DHRC_b = rownames(GSE198904_DHRC_beta)
CpGs_GSE198904_OBS_b = rownames(GSE198904_OBS_beta)
CpGs_GSE74414_MPIP2_b = rownames(GSE74414_MPIP2_beta)

CpG_list = list(
  "GSE125105 MPIP" = CpGs_GSE125105_MPIP_b,
  "GSE72680 GRADY" = CpGs_GSE72680_GRADY_b,
  "GSE113725 RDE" = CpGs_GSE113725_RDE_b,
  "GSE198904 DHRC" = CpGs_GSE198904_DHRC_b,
  "GSE198904 OBS" = CpGs_GSE198904_OBS_b,
  "GSE74414 MPIP2" = CpGs_GSE74414_MPIP2_b
)

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

Intersecting_CpGs = Reduce(intersect, CpG_list)
sapply(CpG_list, function(x) all(Intersecting_CpGs %in% x)) # All -> TRUE

Combined_beta = list(
  GSE125105_MPIP_beta[Intersecting_CpGs,],
  GSE72680_GRADY_beta[Intersecting_CpGs,],
  GSE113725_RDE_beta[Intersecting_CpGs,],
  GSE198904_DHRC_beta[Intersecting_CpGs,],
  GSE198904_OBS_beta[Intersecting_CpGs,],
  GSE74414_MPIP2_beta[Intersecting_CpGs,]
)

# Check
sapply(Combined_beta, function(x) all(rownames(x) == Intersecting_CpGs)) # All -> TRUE

Combined_beta = do.call(cbind, Combined_beta)
Combined_beta = Combined_beta[,combined_pheno$ID]

write.csv(Combined_beta, file="/home/student2024/Desktop/Project_Depression/Combined_beta_non_adj.csv", row.names = TRUE) 
