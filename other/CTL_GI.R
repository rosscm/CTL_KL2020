library(ggplot2)
library(ggrepel)

dataDir <- "C:/Users/rumi/Desktop"
geneDir <- sprintf("%s/PhD/diffGI/diffGI", dataDir)


dir.create(sprintf("%s/RA_results_2", geneDir))
#dt <- format(Sys.time(), "%d%m%20%y%R")
dt <- format(Sys.Date(), "20%y%m%d")

outDir <- sprintf("%s/RA_results_2/%s_out_%s_plots", geneDir, dt, basename(geneDir))
if (!file.exists(outDir)) dir.create(outDir)

datF <- list.files(pattern="result_table_diffGI.*csv", path=geneDir, full.names=T)
dat <- lapply(datF, function(x) read.csv(x))

dat_names <- lapply(datF, function(x) unlist(strsplit(basename(x), "_"))[4])

for (i in 1:length(dat)) {

df <- dat[[i]]

sig_infile <- subset(df, df$fdr < 0.31)
wt_sub <- subset(df, abs(df$WT) > 5 & df$fdr > 0.98)
KO_sub <- subset(df, abs(df[, 3]) > 5 & df$fdr > 0.98)
sig_supp <- subset(sig_infile, sig_infile$diff > 0)
sig_sensi <- subset(sig_infile, sig_infile$diff < 0)


size_sp = -log(df$fdr)*0.8
size_sig = -log(sig_infile$fdr)*0.8
size_WT_sig = -log(wt_sub$fdr)*0.8
size_KO_sig = -log(KO_sub$fdr)*0.8
size_sig_supp <- -log(sig_supp$fdr)*0.8
size_sig_sensi <- -log(sig_sensi$fdr)*0.8

for_label <- rbind(sig_infile, wt_sub, KO_sub)
for_label <- for_label[-c(which(duplicated(for_label))), ]

ggplot(df, aes(WT, df[,3]))+
  geom_hline(yintercept = c(5,-5), linetype = "dashed", color = "gray")+
  geom_vline(xintercept = c(5,-5), linetype = "dashed", color = "gray")+
  geom_point(aes(WT, df[,3], col = "Insignificant Genes", fill = "Insignificant Genes", size = size_sp), pch = 21)+
  geom_point(data = sig_supp, aes(WT, sig_supp[,3], col = "Top Differential Suppressors", fill = "Top Differential Suppressors",
                                    size = size_sig_supp), pch = 21, alpha = 0.9)+
  geom_point(data = sig_sensi, aes(WT, sig_sensi[,3], col = "Top Differential Sensitizers", fill = "Top Differential Sensitizers",
                                  size = size_sig_sensi), pch = 21, alpha = 0.9)+
  geom_point(data = wt_sub, aes(WT, wt_sub[,3], col = "Top Non-Differential Genes", fill = "Top Non-Differential Genes", 
                                size = size_WT_sig), pch = 21, alpha = 0.9)+
  geom_point(data = KO_sub, aes(WT, KO_sub[,3], col = "Top Non-Differential Genes", fill = "Top Non-Differential Genes",
                                size = size_KO_sig), pch = 21, alpha = 0.9)+
  scale_color_manual(values = c("Insignificant Genes" = 'lightgray',"Top Differential Suppressors" = 'black',
                                "Top Differential Sensitizers" = 'black', "Top Non-Differential Genes" = "red"))+
  scale_fill_manual(name = "Type", values = c("Insignificant Genes" = 'lightgray',"Top Differential Suppressors" = 'yellow',
                                              "Top Differential Sensitizers" = 'steelblue', 
                                              "Top Non-Differential Genes" = "red"))+
  guides(Type = guide_legend(order=1),
         size = FALSE,
         color = FALSE)+
  geom_text_repel(data = for_label, aes(WT, for_label[,3]),
                  label = for_label$X, segment.size = 0.2)+
  xlab("RencaHA [DrugZ z score]")+ ylab(paste("RencaHA.", dat_names[i], " KO [DrugZ z score]", sep = ""))+
  scale_y_continuous(breaks=seq(-20,20,5))+
  scale_x_continuous(breaks=seq(-20,20,5))+
  theme_bw() +
  theme(text=element_text(family="sans", size=14),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        axis.line=element_line(colour="black"),
        strip.text=element_text(face="bold"),
        legend.position="bottom")


plot_name <- sprintf("%s/%s_FDR30_noinsig_GI_scatterplot.pdf", outDir, dat_names[i])

ggsave(filename = plot_name, plot = last_plot(), width = 10, height = 10, dpi = 300)

}
