library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sigCpGs <- readLines('/home/student2024/important features common.csv')
all_cpg <- readLines('/home/student2024/cpg sites.csv')
gst_go <- gometh(sig.cpg=sigCpGs, collection="GO", all.cpg = all_cpg,
              plot.bias=TRUE)
topGSA_go <- topGSA(gst)
