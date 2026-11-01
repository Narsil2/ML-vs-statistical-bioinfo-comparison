library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sigCpGs <- readLines('/home/student2024/important features common.csv')
png(file="final_plot.png",width = 600, height = 350)
all_cpg <- readLines('/home/student2024/cpg sites.csv')
gst <- gometh(sig.cpg=sigCpGs, collection="KEGG", all.cpg = all_cpg,
              plot.bias=TRUE)
dev.off()

topGSA <- topGSA(gst)

library(ggplot2)
library(RColorBrewer)
nb.cols <- 20
mycolors <- colorRampPalette(brewer.pal(8, "YlGnBu"))(nb.cols)
topgsa_plot <- ggplot(data =topGSA, aes(y=Description , x = DE, fill =(P.DE))) +
  geom_col() +
  xlab('Gene Names') +
  ylab('Differentially Methylated Genes') +
  ggtitle('Significantly Enriched Pathways') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
topgsa_plot
ggsave(topgsa_plot, filename = 'topgsa.png', device = 'png')

