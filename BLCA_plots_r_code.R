#==========
Section-1: Epigenetic regulator gene aberrations in bladder cancer 
#==========
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[1.1]. plot somatic mutation rate of epigenes in all cancer types
library(ggplot2)
library(reshape)
library(reshape2)
library(patchwork)
myinf1 ="/home/u251079/BLCA_code/BLCA_ms_plots_data/Mutation_Num_All_categories_each_cancer_TCGA.csv"
Mut_rate = read.csv(myinf1, row.names = 1, header = T)
class(Mut_rate)
Mut_rate1<-Mut_rate
colnames(Mut_rate1)[7]<-"Mutation rate"
class(Mut_rate1)

p1<-ggplot(Mut_rate1, aes(x = reorder(rownames(Mut_rate1), -`Mutation rate`), y = `Mutation rate`)) +
  geom_bar(stat = "identity", fill = "salmon", color = "white") +
  theme_classic() + xlab("") +ylab("Mutation rate")+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"))+
  coord_cartesian(expand = FALSE) # this can remove extra space on the y-axis starting point (Start from ZERO)

print(p1)
# Add row names as a new column
Mut_rate1$row_names <- rownames(Mut_rate1)

# Melt the data frame for plotting
data_melt <- melt(Mut_rate1, id.vars = c("row_names", "Mutation rate"))

p2 <- ggplot(data_melt, aes(x = reorder(row_names, -`Mutation rate`), y = reorder(variable, +value))) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2(low = "#2874A6", mid = "white", high = "red", midpoint = 0.35) +
  theme_classic() + xlab("")+ylab("")+
  labs(fill = "Mutation \n   Rate")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_text(size=8, vjust = 0.1))

print(p2)
library(patchwork)
p3<-(p1 / p2) + plot_layout(ncol = 1)

print(p3)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_1/Mut_TCGA.pdf", width = 9, height =5)
print(p3)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[1.2]. Somatic Mutation and CNV pattern pattern in BLCA
rm(list = ls())
# Reorder of x axis problem in this two plots
sm_cnv="/home/u251079/BLCA_code/BLCA_ms_plots_data/BLCA_Epi_sm_cnv.csv" 
sm_cnv= read.csv(sm_cnv, row.names = 1, header = T)
sm_cnv<-sm_cnv[order(sm_cnv$Num_SM, decreasing = T),]

sm_cnv1<-melt(sm_cnv)
names(sm_cnv1)

# aes(x = reorder(Gene, -(value)), y = value, fill=Category)
p4<-ggplot(sm_cnv1, aes(x =reorder(Gene, -value), y = value, fill=Category)) +
  geom_bar(data = subset(sm_cnv1, variable =="Num_SM"),aes(y=value,fill=Category), stat="identity", position = "identity")+
  geom_text(data = subset(sm_cnv1, variable =="Num_SM"),aes(label=value), hjust=0.5, vjust=-0.3, size=2.3, color="black")+
  theme_classic() + 
  xlab("")+
  ylab("Number of Mutations")+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.80, 0.60))+
  coord_cartesian(ylim = c(0, 110), expand = F) # this can remove extra space on the y-axis starting point (Start from ZERO)

print(p4)


p4<-ggplot(sm_cnv1, aes(x =Gene, y = value, fill=Category)) +
  geom_bar(data = subset(sm_cnv1, variable =="Num_SM"),aes(y=value,fill=Category), stat="identity", position = "identity")+
  geom_text(data = subset(sm_cnv1, variable =="Num_SM"),aes(label=value), hjust=0.5, vjust=-0.3, size=2.3, color="black")+
  theme_classic() + 
  xlab("")+
  ylab("Number of Mutations")+ 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.85, 0.70))+
  coord_cartesian(ylim = c(0, 110), expand = F) # this can remove extra space on the y-axis starting point (Start from ZERO)

print(p4)

p5<-ggplot(sm_cnv1, aes(x =Gene, y = value, fill=Category))+
  geom_bar(data = subset(sm_cnv1, variable =="amp.num"),aes(y=value), stat="identity", position = "identity")+
  geom_text(data = subset(sm_cnv1, variable =="amp.num"),aes(label=value), hjust=0.5, vjust=-0.4, size=2.3)+
  geom_bar(data = subset(sm_cnv1, variable =="del.num"),aes(y=-(value)), stat="identity", position = "identity")+
  geom_text(data = subset(sm_cnv1, variable =="del.num"),aes(y=-(value),label=value), hjust=0.5, vjust=1, size=2.3)+
  theme_classic() + 
  xlab("")+
  ylab("Num of CNVs")+
  scale_y_continuous(breaks=seq(-100,100,10),labels=abs(seq(-100,100,10))) + 
  geom_hline(yintercept = 0)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
        axis.text.y = element_text(color="black"),
        legend.position = "none")

print(p5)
# Plot p5 with the same reordered genes as p4
p5 <- ggplot(sm_cnv1, aes(x = reorder(Gene, -abs(value)), y = value, fill = Category)) +
  geom_bar(data = subset(sm_cnv1, variable == "amp.num"), aes(y = value), stat = "identity", position = "identity") +
  geom_text(data = subset(sm_cnv1, variable == "amp.num"), aes(label = value), hjust = 0.5, vjust = -0.4, size = 2.3) +
  geom_bar(data = subset(sm_cnv1, variable == "del.num"), aes(y = -value), stat = "identity", position = "identity") +
  geom_text(data = subset(sm_cnv1, variable == "del.num"), aes(y = -value, label = value), hjust = 0.5, vjust = 1, size = 2.3) +
  theme_classic() + 
  xlab("") +
  ylab("Number of CNVs") +
  scale_y_continuous(breaks = seq(-100, 100, 10), labels = abs(seq(-100, 100, 10))) + 
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none")

print(p5)

# for p1 need to add label to the top of the 
library(patchwork)
p6<-(p4/p5) + plot_layout(ncol = 1)
print(p6)
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_1/SM_CNV_TCGA_BLCA.pdf", width = 10, height =6)
print(p6)
dev.off()


#@@@@@@@@@@@@@@@@@@ New updated above figure
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
# Read and preprocess data
sm_cnv = "/home/u251079/BLCA_code/BLCA_ms_plots_data/BLCA_Epi_sm_cnv.csv" 
sm_cnv = read.csv(sm_cnv, row.names = 1, header = TRUE)
sm_cnv <- sm_cnv[order(sm_cnv$Num_SM, decreasing = TRUE),]
sm_cnv1 <- melt(sm_cnv)

# Extracting Num_SM values for each gene for reordering in p5
num_sm_values <- subset(sm_cnv1, variable == "Num_SM")[, c("Gene", "value")]
names(num_sm_values) <- c("Gene", "Num_SM")
sm_cnv1 <- merge(sm_cnv1, num_sm_values, by = "Gene")

sm_cnv1$value<-ifelse(sm_cnv1$value %in% NA, 0, sm_cnv1$value)
# Plot p4
p4 <- ggplot(sm_cnv1, aes(x = reorder(Gene, -Num_SM), y = value, fill = Category)) +
  geom_bar(data = subset(sm_cnv1, variable == "Num_SM"), aes(y = value, fill = Category), stat = "identity", position = "identity") +
  geom_text(data = subset(sm_cnv1, variable == "Num_SM"), aes(label = value), hjust = 0.5, vjust = -0.3, size = 2.3, color = "black") +
  theme_classic() +  xlab("") +ylab("Number of Somatic Mutations")+
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black"), 
        legend.position = c(0.85, 0.70)) +
  coord_cartesian(ylim = c(0, 110), expand = TRUE)

print(p4)

# Plot p5
p5 <- ggplot(sm_cnv1, aes(x = reorder(Gene, -Num_SM), y = value, fill = Category)) +
  geom_bar(data = subset(sm_cnv1, variable == "amp.num"), aes(y = value), stat = "identity", position = "identity") +
  geom_text(data = subset(sm_cnv1, variable == "amp.num"), aes(label = value), hjust = 0.5, vjust = -0.4, size = 2.3) +
  geom_bar(data = subset(sm_cnv1, variable == "del.num"), aes(y = -value), stat = "identity", position = "identity") +
  geom_text(data = subset(sm_cnv1, variable == "del.num"), aes(y = -value, label = value), hjust = 0.5, vjust = 1, size = 2.3) +
  theme_classic() + 
  xlab("") + 
  ylab("Number of CNVs") + 
  scale_y_continuous(breaks = seq(-100, 100, 10), labels = abs(seq(-100, 100, 10))) + 
  geom_hline(yintercept = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), axis.text.y = element_text(color = "black"), legend.position = "none")

print(p5)
# Combining p4 and p5 with patchwork
combined_plot <- p4 / p5 + plot_layout(ncol = 1)

# Print the combined plot
print(combined_plot)
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_1/SM_CNV_TCGA_BLCA.pdf", width = 10, height =6)
print(combined_plot)
dev.off()



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[1.3]. Volcano plot to show the prognostic association of them with ER gene mutation status of TCGA_BLCA patients
rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
mycox_epi<-"/home/u251079/BLCA_code/BLCA_ms_plots/BLCA_KM_plots/Epi_genes_survival_results.csv"

my_cox1<-read.csv(mycox_epi, row.names = 1, header = T)
names(my_cox1)
my_cox1$neg_log10_p<-c(-log10(my_cox1$PV))
row.names(my_cox1)<-my_cox1$name
my_cox1<-my_cox1[1:47,] # removed missing values
my_cox1$sig<-"Not_sig"
names(my_cox1)
my_cox1$sig[my_cox1$neg_log10_p >= 1.30103 | my_cox1$HR >= 4]<- "sig"

# to show top gene names on the plot we need to create gene name in seperate column
my_cox1$delab<-ifelse(my_cox1$sig == "sig", my_cox1$name, NA)
# add theme
p7<-ggplot(my_cox1, aes(HR, neg_log10_p,  col=my_cox1$sig, label=my_cox1$delab)) +
  geom_vline(xintercept =1, col="#bb0c00", linetype="dashed")+
  geom_hline(yintercept = 1.30103, col="gray", linetype="dashed")+
  geom_point(size=1)+
  scale_color_manual(values = c( "#8080ff","#ff80ff" ),
                     labels= c("Not Significant", "Significant"))+
  # coord_cartesian(ylim = c(min(my_cox1$neg_log10_p)-0.05,max(my_cox1$neg_log10_p)+0.2), ylim =c(0,5))+
  labs(color="Prognostic", 
       x=expression("Hazard Ratio"), y=expression("-log"[10]*"(p-value)"))+
  geom_text_repel(max.overlaps = Inf)+
  theme_classic()+
  theme(legend.position = c(0.85, 0.125),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 3),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key.height= unit(0.10, 'cm'),
        legend.key.width= unit(0.20, 'cm'),
        legend.background = element_rect(size=0.1, linetype="solid", colour ="black"))
print(p7)
# legend.background = element_rect(size=0.3, linetype="solid", colour ="black"),

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_1/1_3_Volcano_plot.pdf", width = 3.5, height =3)
print(p7)
dev.off()





#===========
Section-2: Gene signatures for driver ER gene aberrations 
#===========
[2.1.1-Supplementary figure] volcano plot between EpiRG signature scores to overall survival
rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(readr)
# This data contain adj.p values of surv and cox results using BH method and HR using coxPH
myin="BLCA_code/BLCA_ms_plots_data/Survival_Cox.csv"
Cox_uni_results<-read.csv(myin, row.names = 1, header = T)
my_cox<-Cox_uni_results[1:13,]
my_cox$sig<-"Not_sig"
names(my_cox)
my_cox$sig[my_cox$coxph.p <= 0.05 & my_cox$HR > 0.9]<- "sig"
my_cox$lab<- gsub("uni.noj__", "", rownames(my_cox))

# to show top gene names on the plot we need to create gene name in seperate column
my_cox$delab<-ifelse(my_cox$sig == "sig",my_cox$lab , NA)
my_cox$neg_Log10p<-c(-log10(my_cox$coxph.p))

names(my_cox)
# add theme
plot_Sig_cox<-ggplot(my_cox, aes(HR, neg_Log10p,  col=my_cox$sig, label=my_cox$delab)) +
  geom_vline(xintercept =1, col="#bb0c00", linetype="dashed")+
  geom_hline(yintercept = 1.30103, col="gray", linetype="dashed")+
  geom_point(size=1)+
  scale_color_manual(values = c( "#8080ff","#ff80ff" ),
                     labels= c("Non Significant", "Significant"))+
  # coord_cartesian(ylim = c(min(my_cox1$neg_log10_p)-0.05,max(my_cox1$neg_log10_p)+0.2), ylim =c(0,5))+
  labs(color="Prognostic", 
       x=expression("Hazard Ratio"), y=expression("-log"[10]*"(p-value)"))+
  geom_text_repel(max.overlaps = Inf)+
  theme_classic()+
  theme(legend.position = c(0.85, 0.40),
        legend.title = element_text(size = 5),
        legend.text = element_text(size = 3),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.key.height= unit(0.10, 'cm'),
        legend.key.width= unit(0.20, 'cm'),
        legend.background = element_rect(size=0.1, linetype="solid", colour ="black"))
print(plot_Sig_cox)
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/Supp_fig_TCGA_BLCA_Cox.pdf", width = 3.5, height =3)
print(plot_Sig_cox)
dev.off() 


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[2.1]. Schematic diagram
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[2.2]. Boxplot (Example) → The signatures can distinguish genes with mutations from without (Choi)
rm(list = ls())
library(ggpubr)
library(rstatix)
inpt="/home/u251079/BLCA_code/BLCA_ms_plots_data/chao_box.csv"
c_box<-read.csv(inpt, header = T)
table(c_box$Group~c_box$Gene)
p <- ggplot(c_box, aes(Gene, Value, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.colour = "gray") +  # Add position_dodge to create a gap
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.4, color = "gray") +  # Add color to jitter points
  coord_cartesian(ylim = c(min(c_box$Value), max(c_box$Value) + 10)) +
  labs(title = "Choi_signatures", x = NULL, y = "Signature Score") +
  stat_compare_means(method = "wilcox.test", label = "p", label.y = max(c_box$Value) + 5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black"),
        axis.text.y = element_text( color = "black"),
        legend.position = c(0.8, 0.08),
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(0.35, "cm"),  # Adjust the size of legend keys (symbols)
        legend.text = element_text(size = 4),  # Adjust the size of legend text
        legend.margin = margin(3, 3),
        legend.box.margin = margin(0.1, 0.1, 0.1, 0.1),
        legend.box.background = element_rect(colour = "black", size = 0.5),
        plot.title = element_text(hjust = 0.5, size = 12, color="black")) +
  scale_x_discrete(labels = c("FGFR3", "RB1", "TP53"))


print(p)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/Choi_sig_box_plot_ggplot.pdf", width = 3, height =4)
print(p)
dev.off() 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

[2.3]. Boxplot (Example) → The signatures can distinguish genes with mutations from without using CCLE data
rm(list = ls())
library(ggpubr)
inp="/home/u251079/BLCA_code/BLCA_ms_plots_data/2_C_CCLE_box.csv" 
ER_sig4=read.csv(inp,row.names = 1, header = T)

p1_value <- wilcox.test(Sig_EP300_mut ~ CREBBP_mut, data = dat1, alternative = "g")$p.value


p <- ggplot(ER_sig4, aes(Var2, value.1, fill = value)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.colour = "black") +  # Add position_dodge to create a gap
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.4, color = "blue") +  # Add color to jitter points
  coord_cartesian(ylim = c(min(ER_sig4$value.1), max(ER_sig4$value.1) + 1)) +
  labs(title = "ER_sig in CCLE bladder cancer cell lines", x = NULL, y = "ER_signature scores") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black"),
        axis.text.y = element_text( color = "black"),
        legend.position = c(0.80, 0.70), legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_x_discrete(labels = c("FGFR3", "KDM6A", "RB1", "TP53"))

print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/CCLE_ER_genes_Boxplots.pdf", width = 6, height =4)
print(p)
dev.off() 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

[2.4].Boxplot → CREBP and P300 signatures can cross- predict each others (Exp: Sig.score of CREBP significantly predict P300 Mut vs Wild?) 
rm(list = ls())
in1="/home/u251079/BLCA_code/BLCA_ms_plots_data/Box_sec_2_TCGA.csv"
# Load the necessary library
library(gridExtra)
dat= read.csv(in1, row.names = 1, header = T)
colnames(dat)
dat1=as.data.frame(ifelse(dat[,7:12] == 1, "Mut", "Wild"))
dat1=cbind(dat[,1:6], dat1)
colnames(dat1)
# Sig.score of CREBBP significantly predict P300 Mut vs Wild
p_value <- wilcox.test(Sig_CREBBP_mut ~ EP300_mut, data = dat1, alternative = "g")$p.value
p_value <-round(p_value,3)
p <- ggplot(dat1, aes(x = EP300_mut, y = Sig_CREBBP_mut)) +
  geom_boxplot(width = 0.7, lwd = 1, aes(color = EP300_mut), notch = TRUE, notchwidth = 0.5) +
  geom_jitter(width = 0.2, size = 0.5) +
  labs(x = "EP300", y = "Signature CREBBP") +
  ylim(min(dat1$Sig_CREBBP_mut), max(dat1$Sig_CREBBP_mut) + 10) +
  geom_text(aes(label = paste(signif(p_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_CREBBP_mut) + 4, size = 5, vjust = -1) +
  theme_classic() +
  scale_x_discrete(expand = c(0.25, 0.25)) +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10, color = "black"),
        text = element_text(family = "sans", face = "plain")  # Set font style to plain
  )

print(p)

#@@@@@@@@@@@@@ Jan 23 2024
p <- ggplot(dat1, aes(x = EP300_mut, y = Sig_CREBBP_mut, fill = EP300_mut)) +
  geom_boxplot(width = 0.7, lwd = 1, color = "black", notch = TRUE, notchwidth = 0.5) +  # Set black border here
  labs(x = "EP300", y = "Signature CREBBP") +
  ylim(min(dat1$Sig_CREBBP_mut), max(dat1$Sig_CREBBP_mut) + 10) +
  geom_text(aes(label = paste(signif(p_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_CREBBP_mut) + 4, size = 5, vjust = -1) +
  theme_classic() +
  scale_x_discrete(expand = c(0.25, 0.25)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10, color = "black"),
        text = element_text(family = "sans", face = "plain"))  # Set font style to plain

print(p)
#@@@@@@@@@@@@@


#Sig.score of EP300 significantly predict CREBBP Mut vs Wild
#__________________________________________________________
p1_value <- wilcox.test(Sig_EP300_mut ~ CREBBP_mut, data = dat1, alternative = "g")$p.value
p1_value <-round(p1_value,3)
# 
p1<- ggplot(dat1, aes(x=CREBBP_mut, y=Sig_EP300_mut))+
  geom_boxplot(width=0.7, lwd=1,aes(color=CREBBP_mut),notch = TRUE, notchwidth = 0.5)+
  geom_jitter(width = 0.2, size=0.5)+ labs(x="CREBBP", y="Signature EP300")+
  ylim(min(dat1$Sig_EP300_mut), max(dat1$Sig_EP300_mut) + 10)+
  geom_text(aes(label = paste(signif(p1_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_EP300_mut) + 4, size = 5, vjust = -1)+
  theme_classic()+
  scale_x_discrete(expand = c(0.25, 0.25))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position="none",
        plot.title = element_text(hjust = 0.5, size = 10))
print(p1)

#@@@@@@@@@@@@@ Jan 23 2024
p1 <- ggplot(dat1, aes(x=CREBBP_mut, y=Sig_EP300_mut, fill=CREBBP_mut))+
  geom_boxplot(width = 0.7, lwd = 1, color = "black", notch = TRUE, notchwidth = 0.5) +  # Set black border here
  labs(x="CREBBP", y="Signature EP300")+
  ylim(min(dat1$Sig_EP300_mut), max(dat1$Sig_EP300_mut) + 10)+
  geom_text(aes(label = paste(signif(p1_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_EP300_mut) + 4, size = 5, vjust = -1)+
  theme_classic()+
  scale_x_discrete(expand = c(0.25, 0.25))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position="none",
        plot.title = element_text(hjust = 0.5, size = 10))
print(p1)
#@@@@@@@@@@@@@




#Sig.score of Signature score of CHD6_mut can predicts with CREBBP amp vs wild
#__________________________________________________________
p3_value <- wilcox.test(Sig_CHD6_mut ~ CHD6_amp, data = dat1, alternative = "g")$p.value
p3_value <-round(p3_value,3)
# 
p3<-ggplot(dat1, aes(x=CHD6_amp, y=Sig_CHD6_mut))+
  geom_boxplot(width=0.7, lwd=1,aes(color=CHD6_amp),notch = TRUE, notchwidth = 0.5)+
  geom_jitter(width = 0.2, size=0.5)+ labs(x="CHD6", y="Signature CHD6_mut")+
  ylim(min(dat1$Sig_CHD6_mut), max(dat1$Sig_CHD6_mut) + 10)+
  geom_text(aes(label = paste(signif(p3_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_CHD6_mut) + 4, size = 5, vjust = -1)+
  theme_classic()+scale_x_discrete(expand = c(0.25, 0.25))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position="none",plot.title = element_text(hjust = 0.5, size = 10))+
  scale_x_discrete(labels=c("amp", "wild"))

print(p3)

#@@@@@@@@@@@@@ Jan 23 2024
p3 <- ggplot(dat1, aes(x=CHD6_amp, y=Sig_CHD6_mut, fill=CHD6_amp))+
  geom_boxplot(width = 0.7, lwd = 1, color = "black", notch = TRUE, notchwidth = 0.5) +  # Set black border here
  labs(x="CHD6", y="Signature CHD6_mut")+
  ylim(min(dat1$Sig_CHD6_mut), max(dat1$Sig_CHD6_mut) + 10)+
  geom_text(aes(label = paste(signif(p3_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_CHD6_mut) + 4, size = 5, vjust = -1)+
  theme_classic()+scale_x_discrete(expand = c(0.25, 0.25))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position="none",plot.title = element_text(hjust = 0.5, size = 10))+
  scale_x_discrete(labels=c("amp", "wild"))

print(p3)



##___________Plot and save as one big image----------
library(patchwork)
p4<-p+ p1+ p3+ plot_annotation(title= "TCGA_BLCA", tag_levels = "I")+plot_layout(nrow = 1)
print(p4)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/2_4_sig_TCGA_BLCA_box_plot_new.pdf", width = 6, height =5)
print(p4)
dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[2.5] Dot Plots
# Dot plot TP53
inf_data="/home/u251079/BLCA_code/BLCA_ms_plots_data/2_7_dot_tp53.csv"
data1= read.csv(inf_data, header = T, row.names = 1)
data1<-na.omit(data1)

p1 <- ggplot(data1, aes(reorder(protein_change, -TP53_sig_score), TP53_sig_score, color=Mut_type)) +
  geom_point(na.rm = T)+
  labs(title = "The signature scores of each individual TP53 mutation on protein change", x = NULL, y = "TP53 mut sig.score")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))

print(p1)
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/2_7_TCGA_dot_plot_tp53.pdf", width = 20, height =4)
print(p1)
dev.off()

# Dot plot KDM6A
inf_data1="/home/u251079/BLCA_code/BLCA_ms_plots_data/2_8_dot_kdm6a.csv"
data2= read.csv(inf_data1, row.names = 1, header = T)
p2 <- ggplot(data2, aes(reorder(protein_change, -KDM6A_sig_score), KDM6A_sig_score, color=Mut_type)) +
  geom_point(na.rm = T)+
  labs(title = "The signature scores of each individual KDM6A mutation on protein change", x = NULL, y = "KDM6A mut sig.score")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), plot.title = element_text(size = 10, hjust = 0.5))

print(p2)
library(ggpubr)
# p3<-ggarrange(p1, p2, ncol = 1, nrow = 2)
# print(p3)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/2_8_TCGA_dot_plot_kdm6a.pdf", width = 8, height =4)
print(p2)
dev.off()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Box plots of categories of Somatic Mutations in TP53 and KMD6A 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
#Box plots
inf1= "/home/u251079/BLCA_code/BLCA_ms_plots_data/2_7_dot_tp53.csv"
myinf1 = "/home/u251079/r_program/TCGA_BLCA__epiGene_iRAS_09_17.txt"
data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Convert all the signature scores in data to Z-scores
sig_mean = apply(data, 2, mean, na.rm=T)
sig_sd = apply(data, 2, sd, na.rm=T)

# Initialize a matrix for Z-datas
z_scores = matrix(0, nrow = nrow(data), ncol = ncol(data))
row.names(z_scores)<- rownames(data)
colnames(z_scores)<-colnames(data)
# Calculate Z-datas for each column
for (k in 1:ncol(data)) {
  cat("\r", k)
  z_scores[, k] = (data[, k] - sig_mean[k]) / sig_sd[k]
}
class(z_scores)
z_scores<-as.data.frame(z_scores)
data=z_scores



#@@@@@@@@@@@@@@@@@@
# Box plot TP53
data1= read.csv(inf1, header = T, row.names = 1)

which(data1$ID %in% rownames(data))
table(data1$Mut_type)
# Instead of Wild type extract the data only have Mut type then plot 
 data$TP53_Mut_type<- ifelse(rownames(data) %in%data1$ID , data1$Mut_type, "Wild")

data_tp53=data[,c("uni.noj__TP53_mut","TP53_Mut_type")]
library(reshape2)
library(ggpubr)
names(data_tp53)[1]<-"TP53_mut"
table(data_tp53$TP53_Mut_type)

se=mean(data_tp53$TP53_mut)
se1=sd(data_tp53$TP53_mut)
data_tp53$Z_score<-c((data_tp53$TP53_mut-se)/se1)

# se <- c("5'Flank","Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del")
# se1 <- which(data_tp53$TP53_Mut_type %in% se)
# data_tp53 <- data_tp53[!data_tp53$TP53_Mut_type %in% se, ]
se <- c("5'Flank", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del")
data_tp53$TP53_Mut_type <- ifelse(data_tp53$TP53_Mut_type %in% se, "other", data_tp53$TP53_Mut_type)
# remove Wild
se<-which(data_tp53$TP53_Mut_type%in% "Wild")
data_tp53<- data_tp53[-se,]


# Create a summary table with the counts of samples for each TP53_Mut_type
sample_counts <- table(data_tp53$TP53_Mut_type)

# Create the ggboxplot
p <- ggboxplot(data = data_tp53,
               x = "TP53_Mut_type",
               y = "Z_score",
               color = "TP53_Mut_type",
               palette = "jco",
               legend = "none") +
  ylim(min(data_tp53$Z_score), max(data_tp53$Z_score)) +
  #geom_hline(yintercept = mean(data_tp53$TP53_mut), linetype = 2) +
   # stat_compare_means(method = "anova") +
   # stat_compare_means(label = "p.signif", method = "wilcox.test",
   #                    ref.group = ".all.") +
  labs(title = "TCGA_BLCA_TP53", x = NULL, y = "Signature_TP53_mut (Z-Score)") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, color = "black"), 
        axis.title.y = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=8, color = "black"),
        axis.text.y = element_text(size = 10, color = "black")) +
  scale_x_discrete(labels = paste(names(sample_counts), "\n(N =", sample_counts, ")"))

print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/TP53_Box_SM_categories.pdf", width = 3, height =3, units = "in", res = 900)
print(p)
dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@
# Box plots KDM6A
myinf2= "/home/u251079/BLCA_code/BLCA_ms_plots_data/2_8_dot_kdm6a.csv"
data2= read.csv(myinf2,header = T, row.names = 1)
which(data2$ID %in% rownames(data))

data$KDM6A_Mut_type<- ifelse(rownames(data) %in%data2$ID , data2$Mut_type, "Wild")

data_kdm6A=data[,c("uni.noj__KDM6A_mut","KDM6A_Mut_type")]
library(reshape2)
library(ggpubr)
names(data_kdm6A)[1]<-"KDM6A_mut"
table(data_kdm6A$KDM6A_Mut_type)


# Create a summary table with the counts of samples for each Mut_type
sample_counts <- table(data_kdm6A$KDM6A_Mut_type)

# Remove Wild
se=which(data_kdm6A$KDM6A_Mut_type %in% "Wild")
data_kdm6A<- data_kdm6A[-se,]

# Create the ggboxplot
p1 <- ggboxplot(data = data_kdm6A,
                x = "KDM6A_Mut_type",
                y = "KDM6A_mut",
                color = "KDM6A_Mut_type",
                palette = "jco",
                legend = "none") +
  ylim(min(data_kdm6A$KDM6A_mut), max(data_kdm6A$KDM6A_mut)) +
  #geom_hline(yintercept = mean(data_kdm6A$KDM6A_mut), linetype = 2) +
  # stat_compare_means(method = "anova", label.y = 120) +
  # stat_compare_means(label = "p.signif", method = "wilcox.test",
  #                    ref.group = ".all.") +
  labs(title = "TCGA_BLCA_KDM6A", x = NULL, y = "Signature_KDM6A_mut (Z-score)") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, color = "black"), 
        axis.title.y = element_text(size = 10,color = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=8, color = "black"),
        axis.text.y = element_text(size = 10, color = "black")) +
  scale_x_discrete(labels = paste(names(sample_counts), "\n(N =", sample_counts, ")"))

print(p1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/KDM6A_Box_SM_categories.pdf", width = 3.5, height =3)
print(p1)
dev.off()

library(patchwork)
p2<- p+p1+plot_layout(ncol = 2)
print(p2)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/Box_plot_SM_categories.pdf", width = 6, height =4)
print(p2)
dev.off() 



#===============
Section-3:Association of driver gene signatures with prognosis 
#===============
[3.1].Volcano plot (Global) → X-HR vs. -log10(P) for all signatures

rm(list = ls())
library(ggrepel)
library(tidyverse)
library(dplyr)
inf="/home/u251079/BLCA_code/BLCA_ms_plots_data/HR_with_Z_scores.csv"
my_cox_global=read.csv(inf, row.names = 1, header = T)
rownames(my_cox_global)
dim(my_cox_global)
my_cox_global<-my_cox_global[1:21,]

se= c("uni.noj__TP53_del","uni.noj__CHD6_ampmut","uni.noj__CHD7_ampmut", "uni.noj__CREBBP_delmut","uni.noj__EP300CREBBP_mut","uni.noj__EP300CREBBP_delmut", "uni.noj__TP53_delmut")
# Remove the rows with names matching the values in 'se'
my_cox_global <- my_cox_global[!(rownames(my_cox_global) %in% se), ]
rownames(my_cox_global)<-gsub("uni.noj__","", rownames(my_cox_global))
my_cox_global$GSE13507_coxph.p_Log10<-c(-log10(my_cox_global$GSE13507_coxph.p))
my_cox_global$Sjodahl_coxph.p_Log10<-c(-log10(my_cox_global$GSE32894_coxph.p))

my_cox_global$GSE13507_level<-ifelse(my_cox_global$GSE13507_coxph.p_Log10 >=2, "Sig", NA)
my_cox_global$GSE13507_labs1<-ifelse(my_cox_global$GSE13507_coxph.p_Log10 >=2, rownames(my_cox_global),NA)

my_cox_global<-my_cox_global %>% 
  mutate(Exp_GSE13507= case_when(GSE13507_coxph.p_Log10 >=2 & GSE13507_HR <1 ~ "Down",
                        GSE13507_coxph.p_Log10 >=2 & GSE13507_HR >1 ~ "Up"))

pl1<-ggplot(my_cox_global, aes(GSE13507_HR, GSE13507_coxph.p_Log10)) +
  geom_vline(xintercept = 1, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed") +
  geom_point(aes(color=Exp_GSE13507) ,size=2)+
  coord_cartesian(xlim = c(0, 2.5), ylim = c(min(my_cox_global$GSE13507_coxph.p_Log10), 
                                             max(my_cox_global$GSE13507_coxph.p_Log10)+0.5))+
  geom_text_repel(aes(label=GSE13507_labs1),size=2)+
  theme_classic()+
  labs(title = "GSE13507 driver ER signatures",  y=expression("-log"[10]*"(p-value)"), x="Hazard Ratio")+
  scale_color_manual(values = c("dodgerblue3", "firebrick3"))+
  theme( plot.title = element_text(size=10, hjust = 0.5),
         legend.position = "none")
  
print(pl1)


#-------
my_cox_global$GSE32894_level<-ifelse(my_cox_global$Sjodahl_coxph.p_Log10 >=2, "Sig", NA)
my_cox_global$GSE32894_labs2<-ifelse(my_cox_global$Sjodahl_coxph.p_Log10 >=2, rownames(my_cox_global),NA)
# In inside the mutate dont use "<-" this symbol before case_when, it gives us NULL colname
my_cox_global<- my_cox_global %>% 
  mutate(Exp_GSE32894 = case_when(Sjodahl_coxph.p_Log10 >=2 & GSE32894_HR >1 ~ "Up",
                                   Sjodahl_coxph.p_Log10 >=2 & GSE32894_HR <1 ~ "Down"))

pl2<-ggplot(my_cox_global, aes(GSE32894_HR, Sjodahl_coxph.p_Log10)) +
  geom_vline(xintercept = 1, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 2, color = "darkred", linetype = "dashed") +
  geom_point(aes(color=Exp_GSE32894),size=2)+
  coord_cartesian(xlim = c(0, 3))+
  theme_classic()+
  labs(title = "GSE32894 driver ER signatures", y=expression("-log"[10]*"(p-value)"), x="Hazard Ratio")+
  scale_color_manual(values = c("dodgerblue3", "firebrick3"))+
  theme(legend.position = "none",
        plot.title = element_text(size=10, hjust = 0.5))+
  geom_text_repel(aes(label=GSE32894_labs2), size=2)
print(pl2)

library(patchwork)
p2<-p+p1+plot_layout(ncol = 2)
print(p2)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_1_Global_Vol_plot.pdf", width = 8, height =6)
print(p2)
dev.off() 

#@@@@@@@@@@@ COX HZ

mycox_HZ<-my_cox_global[,c(5,10)]

print(mycox_HZ)

                GSE13507_HR GSE32894_HR
ARID1A_mut        1.1294787   0.8649591
CHD6_mut          1.6940169   1.8821384
CHD7_mut          1.8365340   1.6887167
CREBBP_mut        1.6529939   2.1597317
EP300_mut         2.0594361   1.7753684
KDM6A_mut         0.6003847   0.4797199
TP53_mut          2.1802777   2.6557218
CHD6_amp          0.7304972   0.5742282
CHD7_amp          0.7657386   0.5984827
PRDM9_amp         1.8774413   1.9728366
CHD3_del          0.6331398   0.3912820
CREBBP_del        1.1946442   0.9875792
HDAC4_del         0.8648735   0.6198821
PHF23_del         0.6460149   0.3877488
EP300CREBBP_mut   1.8707262   1.8943059


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list = ls())
library(survival)
library(ggplot2)
library(survminer)
library(gridExtra)
# Sjodahl_GSE32894_data_for_KM with 224 patient informatyion 
input1= "/home/u251079/BLCA_code/BLCA_ms_plots_data/Sjodahl_GSE32894_data_for_KM.csv"
# GSE13507 data with 165 survival patients
input2= "/home/u251079/BLCA_code/BLCA_ms_plots_data/GSE13507_data_for_KM.csv"

dat1= read.csv(input1, row.names = 1, header = T)
dat2= read.csv(input2, row.names = 1, header = T)

[3.1.a]. Supplementary= KM plot (E) → EP300CREBBP_mut (in 2 datasets)
                 GSE13507_HR GSE32894_HR
EP300CREBBP_mut   1.8707262   1.8943059

# GSE32894 data
fit1<- survfit(Surv(t.dfs, e.dfs)~EP300CREBBP_mut, dat1)
surv_pvalue(fit1)

p1 <- ggsurvplot(
  fit1,
  data = dat1,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low", "High"),
  legend.title = "Signature",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold", hjust=0.5, vjust=0.5),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE32894 EP300CREBBP_mut score",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p1$plot <- p1$plot + annotate(
  "text",
  x = 30, y = 0.15,
  label = "HR = 1.89 \n p < 0.0001",
  size = 3.5,
  hjust = 0
)

print(p1)

# GSE13507 data
names(dat2)
fit2<- survfit(Surv(t.surv, e.surv)~EP300CREBBP_mut, dat2)
surv_pvalue(fit2)


p2 <- ggsurvplot(
  fit2,
  data = dat2,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low", "High"),
  legend.title = "Signature",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold"),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE13507 EP300CREBBP_mut score",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())
print(p2)

# Add HR and p-value to the plot
p2$plot <- p2$plot + annotate(
  "text",
  x = 45, y = 0.15,
  label = "HR = 1.87 \n p = 0.012",
  size = 3.5,
  hjust = 0
)
print(p2)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/SuppKM_EP300CBP_plot.pdf", width = 12, height =4)
print(grid.arrange(p1$plot, p$plot, ncol = 2, nrow=1) )
dev.off() 





3.2.KM plot (E) → TP53 (in 2 datasets)

                GSE13507_HR GSE32894_HR
TP53_mut          2.1802777   2.6557218

#$$$$$$$$$$$$$$$$ GSE32894 data
fit1<- survfit(Surv(t.dfs, e.dfs)~TP53_mut, dat1)
surv_pvalue(fit1)

p1 <- ggsurvplot(
  fit1,
  data = dat1,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (112)", "High (112)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold", hjust=0.5, vjust=0.5),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE32894 TP53_mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p1$plot <- p1$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 2.65 \n p < 0.0001",
  size = 3.5,
  hjust = 0
)

print(p1)


#$$$$$$$$$$$$$$$$$$$$ GSE13507 data HR :    2.1802777
fit2<- survfit(Surv(t.surv, e.surv)~TP53_mut, dat2)
surv_pvalue(fit2)

p2 <- ggsurvplot(
  fit2,
  data = dat2,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (82)", "High (83)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold"),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE13507 TP53_Mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p2$plot <- p2$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 2.18 \n p < 0.0001",
  size = 3.5,
  hjust = 0
)
print(p2)

#+++++++++++++++++++++++++++++++++++++++++++++++++
[3.3].KM plot (E) → CREBBP (in 2 datasets)


                GSE13507_HR GSE32894_HR
CREBBP_mut        1.6529939   2.1597317

# $$$$$$$$$$$$$$$$$.  GSE32894_data
fit3<- survfit(Surv(t.dfs, e.dfs)~CREBBP_mut, dat1)

surv_pvalue(fit3)

p3 <- ggsurvplot(
  fit3,
  data = dat1,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (112)", "High (112)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold", hjust=0.5, vjust=0.5),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE32894 CREBBP_mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p3$plot <- p3$plot + annotate(
  "text",
  x = 45, y = 0.15,
  label = "HR = 2.16 \n p < 0.0001",
  size = 3.5,
  hjust = 0
)

print(p3)


#$$$$$$$$$$$$$$$$$ GSE13507 data 
GSE13507_HR : 1.6529939

fit4<- survfit(Surv(t.surv, e.surv)~CREBBP_mut, dat2)
surv_pvalue(fit4)

p4 <- ggsurvplot(
  fit4,
  data = dat2,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (82)", "High (83)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold"),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE13507 CREBBP_Mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p4$plot <- p4$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 1.65 \n p = 0.012",
  size = 3.5,
  hjust = 0
)


# Print the modified plot
print(p4)

#+++++++++++++++++++++++++++++++++++++++++++++++++
[3.4].KM plot (E) → KDM6A (in 2 datasets)

                  GSE13507_HR GSE32894_HR
KDM6A_mut         0.6003847   0.4797199

#$$$$$$$$$$$$$$ GSE32894_data
fit5<- survfit(Surv(t.dfs, e.dfs)~KDM6A_mut, dat1)
surv_pvalue(fit5)
p5 <- ggsurvplot(
  fit5,
  data = dat1,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (112)", "High (112)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold", hjust=0.5, vjust=0.5),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE32894 KDM6A_mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p5$plot <- p5$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 0.48 \n p < 0.0001",
  size = 3.5,
  hjust = 0
)

print(p5)

#$$$$$$$$$$$$$$$$$$$ GSE13507 data 
GSE13507_HR:  0.6003847  

fit6<- survfit(Surv(t.surv, e.surv)~KDM6A_mut, dat2)
surv_pvalue(fit6)
# Create the Kaplan-Meier plot for
p6 <- ggsurvplot(
  fit6,
  data = dat2,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (82)", "High (83)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold"),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE13507 KDM6A_Mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme = theme_classic())

# Add HR and p-value to the plot
p6$plot <- p6$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 0.6 \n p = 0.015",
  size = 3.5,
  hjust = 0
)


# Print the modified plot
print(p6)


3.5.KM plot (E) → EP300 (in 2 datasets)

                 GSE13507_HR GSE32894_HR
EP300_mut         2.0594361   1.7753684

#$$$$$$$$$$$$$ GSE32894_data
fit7<- survfit(Surv(t.dfs, e.dfs)~EP300_mut, dat1)
surv_pvalue(fit7)

p7 <- ggsurvplot(
  fit7,
  data = dat1,
  pval = F, 
  conf.int = FALSE,
  size = 1,
  legend.labs = c("Low (112)", "High (112)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold", hjust=0.5, vjust=0.5),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE32894 EP300_mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme =  theme_classic())

# Add HR and p-value to the plot
p7$plot <- p7$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 1.78 \n p < 0.0001",
  size = 3.5,
  hjust = 0
)

print(p7)

#$$$$$$$$$$$$$$$$$$$ GSE13507 data 
GSE13507_HR: 2.0594361

fit8<- survfit(Surv(t.surv, e.surv)~EP300_mut, dat2)
surv_pvalue(fit8)
# Create the Kaplan-Meier plot for
p8 <- ggsurvplot(
  fit8,
  data = dat2,
  pval = F,
  conf.int = F,
  size = 1,
  legend.labs = c("Low (82)", "High (83)"),
  legend.title = "signature score",
  legend = c(0.2, 0.2),
  risk.table = FALSE,
  font.x = c(12, "black"),
  font.y = c(12, "black"),
  font.main = c(10, "black", "bold"),
  xlab = "Survival Time (months)",
  ylab = "Probability of Survival",
  title = "GSE13507 EP300_Mut",
  font.tickslab = c(12, "black"),
  palette = c("blue", "#C70039"),
  ggtheme =  theme_classic())

# Add HR and p-value to the plot
p8$plot <- p8$plot + annotate(
  "text",
  x = 55, y = 0.15,
  label = "HR = 2.06 \n p = 0.0066",
  size = 3.5,
  hjust = 0
)


# Print the modified plot
print(p8)


library(patchwork)
KM_plots<-p1$plot+p2$plot+p3$plot+p4$plot+p5$plot+p6$plot+p7$plot+p8$plot+plot_layout(ncol = 2, nrow = 4)
print(KM_plots)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/8_KM_plot.pdf", width =10 , height =15)
print(KM_plots)
dev.off() 

#@@@@@@@@@@ 01-09-2024

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3A_Test.pdf", width =3.5 , height =3.5)
print(pl2)
dev.off() 



pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3B_Test.pdf", width =3.5 , height =3.5)
print(p2)
dev.off() 

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3C_Test.pdf", width =3.5 , height =3.5)
print(p4)
dev.off() 

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3D_Test.pdf", width =3.5 , height =3.5)
print(p6)
dev.off() 



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[3.5].Forest plot → CREBBP after adjusting for clinical factors (in 1 dataset)
rm(list = ls())
library(survival)
library(survminer)
# GSE13507_data 
inp="/home/u251079/BLCA_code/BLCA_ms_plots_data/GSE13507_for_Box_Forest.csv"
info=read.csv(inp, row.names = 1, header=T)
colnames(info)<-gsub("uni.noj__", "", colnames(info))

names(info)
fit.coxph <- coxph(Surv(t.surv, e.surv) ~CREBBP_mut+ SEX + AGE + invasiveness+ Grade+ T_stage, data=info)
p<-ggforest(fit.coxph, info, main = "Hazard ratio_GSE13507")
print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_5_Forest_plot.pdf", width = 7, height =5)
print(p)
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[3.6].Boxplot → association with invasiveness CREBBP (muscle-invasive vs. non-invasive)
# GSE13507_data 
se=as.factor(info$invasiveness)
table(se)

plot_violin<-ggplot(info, aes(as.factor(invasiveness), CREBBP_mut))+
  geom_violin(aes(fill=as.factor(invasiveness)), scale = "width")+
  geom_boxplot(width=0.3, lwd=0.5)+
  labs(title="GSE13507", x=NULL,y="Signature CREBBP_mut")+
  geom_jitter(position = position_dodge(width = 0.1), alpha = 0.2, color = "blue")+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  coord_cartesian(ylim = c(min(info$CREBBP_mut), max(info$CREBBP_mut)+5))+
  theme_bw()+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))+
  scale_x_discrete(labels = c("Non Muscle-Inv (n=103)", "Muscle-Inv (n=62)"))

print(plot_violin)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_6_violin_plot.pdf", width = 4, height =5)
print(plot_violin)
dev.off()

3.7.Boxplot → CREBBP in different stages
# with "GSE13507 data with 165 patients
'''plot_box<-ggplot(info, aes(T_stage, CREBBP_mut, fill=factor(T_stage)))+
  geom_boxplot()+
  labs(title ="GSE13507", x="Tumor Stages", y="Signature CREBBP_mut")+
  coord_cartesian(ylim = c(min(info$CREBBP_mut), max(info$CREBBP_mut)+8))+
  geom_hline(yintercept = mean(info1$CREBBP_mut), linetype = 2) +
  stat_compare_means(method = "anova", label.y = 64) +        # Add global ANOVA p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.") +
  theme_classic()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 10))
print(plot_box)'''

library(ggplot2)
table(info$T_stage, info$invasive)
# Create a factor with custom levels for T_stage
info$T_stage1 <- factor(info$T_stage, levels = c("Ta","T1", "T2", "T3", "T4"))

# Calculate the number of samples for each T_stage
sample_counts <- info %>%
  group_by(T_stage1) %>%
  summarise(N = n())

# Create the reordered box plot
plot_box<- ggplot(info, aes(T_stage1, CREBBP_mut, fill = factor(T_stage1))) +
  geom_boxplot() +
  labs(title = "GSE13507", x = "Tumor Stages", y = "Signature CREBBP_mut") +
  coord_cartesian(ylim = c(min(info$CREBBP_mut), max(info$CREBBP_mut) + 8)) +
  geom_hline(yintercept = mean(info$CREBBP_mut), linetype = 2) +
  stat_compare_means(method = "anova", label.y = 64) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.",
                     label.y = c(27, 50, 58, 58, 50)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 10))+
  scale_x_discrete(labels = c("Ta (n=24)","T1 (n=80)", "T2 (n=31)", "T3 (n=19)", "T4 (n=11)"))
                   
print(plot_box)
                   
                   
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_7_box_plot.pdf", width = 4, height =4)
print(plot_box)
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# With Sjodahl_GSE32894_data with 224 patient 
table(info1$stage)
info1$stage<-ifelse(info1$stage == "T4", "T3", info1$stage)
# Create a factor with custom levels for T_stage
info1$stage1 <- factor(info1$stage, levels = c("Ta","T1", "T2", "T3"))

# Calculate the number of samples for each T_stage
sample_counts <- info1 %>%
  group_by(stage1) %>%
  summarise(N = n())

# T1  T2  T3  T4  Ta 
# 63  43   7   1 110 
# Create the reordered box plot
plot_box1<- ggplot(info1, aes(stage1, CREBBP_mut, fill = factor(stage1))) +
  geom_boxplot() +
  labs(title = "GSE32894", x = "Tumor Stages", y = "Signature CREBBP_mut") +
  coord_cartesian(ylim = c(min(info1$CREBBP_mut), max(info1$CREBBP_mut) + 8)) +
  geom_hline(yintercept = mean(info1$CREBBP_mut), linetype = 2) +
  stat_compare_means(method = "anova", label.y = max(info1$CREBBP_mut)) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.",
                     label.y = c(60, 80, 104,104)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 10))+
  scale_x_discrete(labels = c("Ta (n=110)","T1 (n=63)", "T2 (n=43)", "T3 (n=8)"))

print(plot_box1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_7_box_plot_GSE32894.pdf", width = 4, height =4)
print(plot_box1)
dev.off()

3.8.Boxplot → CREBBP in different molecular subtypes
library(ggpubr)
xx = as.character(info1$molecular_subtype)
se = grep("^MS1a", xx)
xx[se] = "urothelial_A"
se = grep("^MS1b", xx)
xx[se] = "urothelial_A"

se = grep("^MS2a.1", xx)
xx[se] = "genomically unstable"
se = grep("^MS2a.2", xx)
xx[se] = "genomically unstable"

se = grep("^MS2b2.1", xx)
xx[se] = "urobasal B"

se = grep("^MS2b2.2", xx)
xx[se] = "SCC-like"

se = grep("^MS2b.1", xx)
xx[se] = "infiltrated"

info1$Mol_subtypes = xx
info1$Mol_subtypes <- as.character(info1$Mol_subtypes)
# Sjodahl_GSE32894_data with 224 patient clinical and signature information 
inp1="/home/u251079/BLCA_code/BLCA_ms_plots_data/GSE32894_for_BOX_Forest.csv"
info1=read.csv(inp1, header = T, row.names = 1)
colnames(info1)<-gsub("uni.noj__", "", colnames(info1))

p3<-ggboxplot(data = info1,
              x = "Mol_subtypes",
              y = "CREBBP_mut",
              color = "Mol_subtypes",
              palette = "jco",
              legend="none")+
  ylim(min(info1$CREBBP_mut), max(info1$CREBBP_mut) + 25) +
  geom_hline(yintercept = mean(info1$CREBBP_mut), linetype = 2) + # Add horizontal line at base mean
  stat_compare_means(method = "anova", label.y = 120) +        # Add global ANOVA p-value
  stat_compare_means(label = "p.signif", method = "wilcox.test",
                     ref.group = ".all.",
                     label.y = c(50,105, 90, 110, 50)) +  # Pairwise comparison against all
  labs(title="GSE32894_data", x = "Molecular Subtype", y = "Signature_CREBBP_mut")+ 
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.1, hjust = 0.5))
print(p3)

# Calculate mean scores for each molecular subtype
mean_scores <- aggregate(CREBBP_mut ~ Mol_subtypes, data = info1, FUN = mean)

# Reorder the factor levels of Mol_subtypes based on mean scores
info1$Mol_subtypes <- factor(info1$Mol_subtypes, levels = mean_scores[order(mean_scores$CREBBP_mut), "Mol_subtypes"])

p3 <- ggboxplot(data = info1,
                x = "Mol_subtypes",
                y = "CREBBP_mut",
                color = "black",  # Set border color to black
                fill = "Mol_subtypes",   # Fill color based on molecular subtypes
                palette = "jco",
                legend = "none") +
  ylim(min(info1$CREBBP_mut), max(info1$CREBBP_mut) + 5) +
  geom_hline(yintercept = mean(info1$CREBBP_mut), linetype = 2, color = "black") + # Horizontal line for overall mean in black
  stat_summary(fun = "mean", geom = "line", color = "black", size = 1) +  # Add line indicating mean for each group
  stat_compare_means(method = "anova", label.y = 80) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", label.y = c(50, 50, 83, 105, 105)) +
  labs(title = "Wilcox.test between sig_score vs mol_subtyp", x = "Molecular Subtype", y = "Signature_CREBBP_mut") + 
  theme(plot.title = element_text(hjust = 0.5))

print(p3)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_8_Mol_sub_plot.pdf", width = 7, height =4)
print(p3)
dev.off()

3.9.Forest plot ⇒ results from stratified analysis (MI vs. non-MS)
# Fit a Cox proportional hazards model using forest plot
# Sjodahl_GSE32894_data with 224 patient clinical and signature information 
inp1="/home/u251079/BLCA_code/BLCA_ms_plots_data/GSE32894_for_BOX_Forest.csv"
info1=read.csv(inp1, header = T, row.names = 1)
colnames(info1)<-gsub("uni.noj__", "", colnames(info1))

colnames(info1)
fit.coxph1 <- coxph(Surv(t.dfs, e.dfs) ~CREBBP_mut+gender+ Mol_subtypes, info1)
n<-ggforest(fit.coxph1, info1, main="Molecular Subtypes-GSE32894",cpositions = c(0, 0.15,0.35),
            fontsize = 1)
print(n)

#@@@@@@@@@@@ ADD REF
# Ensure Mol_subtypes is a factor and set "Urothelial-A" as the reference
info1$Mol_subtypes1 <- factor(info1$Mol_subtypes, levels = c("urothelial_A", levels(info1$Mol_subtypes)[!levels(info1$Mol_subtypes) == "urothelial_A"]))

# Fit the Cox model
fit.coxph1 <- coxph(Surv(t.dfs, e.dfs) ~ CREBBP_mut+gender+Mol_subtypes1, data = info1)

# Generate the forest plot
n <- ggforest(fit.coxph1, data = info1, main = "Molecular Subtypes-GSE32894", 
              cpositions = c(0, 0.15, 0.35), fontsize = 1)
print(n)

# # Ensure Mol_subtypes is a factor and set "Urothelial-A" as the reference
#info1$Mol_subtypes <- factor(info1$Mol_subtypes, levels = c("urothelial_A", levels(info1$Mol_subtypes)[!levels(info1$Mol_subtypes) == "urothelial_A"]))
# # Ensure Mol_subtypes is a factor
# #@@@@@@@@@@@ REMOVE REF
# info1$Mol_subtypes <- factor(info1$Mol_subtypes)
# # Remove "urothelial_A" from the levels
# info1$Mol_subtypes <- factor(info1$Mol_subtypes, levels = levels(info1$Mol_subtypes)[levels(info1$Mol_subtypes) != "urothelial_A"])
# 

# 
# # Fit the Cox model
# fit.coxph1 <- coxph(Surv(t.dfs, e.dfs) ~ CREBBP_mut+age+Mol_subtypes, data = info1)
# 
# # Generate the forest plot
# n <- ggforest(fit.coxph1, data = info1, main = "Molecular Subtypes-GSE32894", 
#               cpositions = c(0, 0.15, 0.35), fontsize = 1)
# print(n)


?ggforest

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_9_Forestplot.pdf", width = 10, height =8)
print(n)
dev.off() 


# library(patchwork)
# all_plots<-p1+plot_violin+plot_box+p3+n+plot_layout(ncol = 3, nrow = 2,guides = "collect")+plot_annotation(tag_levels = "A")
# print(all_plots)
# 
# pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_3/3_5to3_9_allplots.pdf", width = 15, height =10, units = "in", res = 1200)
# print(all_plots)
# dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-----------
Section-4:
#------------------------
[4.1] Correlation Heatmap signature scores verses Thorsson scores  
rm(list = ls())
myinf1 <- "/mount/ictr1/chenglab/cc59/PubDat/Dataset/Firehose/done/Thorsson_2018_TCGA_immunelandscape.csv"
myinf2 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/BladderCancer/IntGen/data/TCGA_BLCA__epiGene_iRAS.txt"

#------------------------
# Epi gene signatures with -log10(p-value) in all patients of BLCA up and down according to the BETA coefficients up (beta>0) down (beta<0)
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="", check.names=F)
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
# genes are weighted based on their association with ER sig mutation status
# up score- down score 
xx = dat1-dat2 
# High score indicates highly association with mutation status, uncorrelated genes received small weights
colnames(xx) = gsub("\\.up", "", colnames(dat1))
score =xx 

dim(score)

data1 = read.table(myinf1, sep=",", header=T, row.names=1, quote="")
se = which(data1$TCGA.Study=="BLCA")
data1 = data1[se,]
names(data1)[1:3]
names(data1)[32:35]
se = c(4:31, 36:63)
data1 = data1[,se]
se2<-c("Proliferation","Leukocyte.Fraction","Stromal.Fraction","Lymphocyte.Infiltration.Signature.Score","TGF.beta.Response", "TCR.Richness", "BCR.Richness","Macrophage.Regulation")
1. **Proliferation**:
  - Proliferation refers to the process of cell division and reproduction, leading to an increase in the number of cells. In the context of cancer, high proliferation rates can indicate aggressive tumor behavior and poor prognosis.

2. **Leukocyte Fraction**:
  - Leukocytes, also known as white blood cells, are a crucial component of the immune system. The leukocyte fraction typically refers to the proportion or percentage of leukocytes in a sample, such as blood or tumor tissue. Changes in leukocyte fraction can indicate immune system activation or infiltration into tissues, which is relevant in the context of immune responses to cancer.

3. **Stromal Fraction**:
  - The stroma is the supportive connective tissue framework within an organ or tissue, providing structural support and facilitating interactions between cells. The stromal fraction refers to the proportion of stromal cells within a tissue sample. Changes in stromal fraction can be associated with alterations in tissue architecture and microenvironment, which can influence tumor growth and behavior.

4. **Macrophage Regulation**:
  - Macrophages are a type of white blood cell that plays a key role in immune responses and tissue homeostasis. Macrophage regulation likely refers to the modulation or control of macrophage activity within a biological system. Macrophages can have both pro-tumor and anti-tumor functions depending on their activation state and the context of the tumor microenvironment.

5. **Lymphocyte Infiltration Signature Score**:
  - Lymphocytes are a type of white blood cell involved in adaptive immune responses. The lymphocyte infiltration signature score is a quantitative measure of the presence and activity of lymphocytes within a tissue sample, particularly in the context of tumor-infiltrating lymphocytes (TILs). High lymphocyte infiltration is often associated with better prognosis in cancer patients due to enhanced anti-tumor immune responses.

6. **TGF-beta Response**:
  - Transforming growth factor-beta (TGF-beta) is a cytokine involved in various cellular processes, including cell growth, differentiation, and immune regulation. TGF-beta response refers to the cellular or tissue-level response to TGF-beta signaling. Dysregulation of TGF-beta signaling pathway is implicated in cancer progression and immune evasion.

7. **TCR Richness**:
  - TCR (T-cell receptor) richness refers to the diversity of T-cell receptor sequences within a population of T-cells. TCR diversity is crucial for recognizing a wide range of antigens, including those derived from tumor cells. Higher TCR richness is often associated with better immune surveillance and response to tumors.

8. **BCR Richness**:
  - BCR (B-cell receptor) richness is analogous to TCR richness but refers to the diversity of B-cell receptor sequences within a population of B-cells. BCR diversity is important for recognizing and responding to a wide array of antigens, including those associated with tumor cells. Similar to TCR richness, higher BCR richness may indicate better immune responsiveness against tumors.

data1<-data1[,se2]


comxx = intersect(row.names(score), row.names(data1))
data1 = data1[comxx,]
score = score[comxx,]
dim(data1)		## 408  56

#@@@@@@@@@@@@ Corrlation analysis only using uniadj of 13 aberrations
score1=score[,1:13]
names(score1)
colnames(score1)<-gsub("uni.noj__","",colnames(score1))


comxx = intersect(row.names(score1), row.names(data1))
data1 = data1[comxx,]
score1 = score1[comxx,]

# calculation of z-scores for signatures 
class(score1)
# Calculate the mean and standard deviation for each column
sig_mean = apply(score1, 2, mean, na.rm=T)
sig_sd = apply(score1, 2, sd, na.rm=T)

# Initialize a matrix for Z-data1s
z_scores = matrix(0, nrow = nrow(score1), ncol = ncol(score1))
row.names(z_scores)<- rownames(score1)
colnames(z_scores)<-colnames(score1)
# Calculate Z-data1s for each column
for (k in 1:ncol(score1)) {
  cat("\r", k)
  z_scores[, k] = (score1[, k] - sig_mean[k]) / sig_sd[k]
}
class(z_scores)
z_scores<-as.data.frame(z_scores)
score1=z_scores

dim(score1)	

# COnverting data1 in to Z-scores
# calculation of z-scores for signatures 
class(data1)
# Calculate the mean and standard deviation for each column
sig_mean = apply(data1, 2, mean, na.rm=T)
sig_sd = apply(data1, 2, sd, na.rm=T)

# Initialize a matrix for Z-data1s
z_scores = matrix(0, nrow = nrow(data1), ncol = ncol(data1))
row.names(z_scores)<- rownames(data1)
colnames(z_scores)<-colnames(data1)
# Calculate Z-data1s for each column
for (k in 1:ncol(data1)) {
  cat("\r", k)
  z_scores[, k] = (data1[, k] - sig_mean[k]) / sig_sd[k]
}
class(z_scores)
z_scores<-as.data.frame(z_scores)
data1=z_scores

dim(data1)	


xx = cor(data1, score1, method="s", use="pair")
heatmap(xx, Colv  = NA)
# Assuming correlation_matrix is your correlation matrix
library(pheatmap)

# Create a square heatmap
p<-pheatmap(
  xx,
  cluster_cols = F,
  cluster_rows = F,
  cellwidth = 15,  # Adjust the cell width to make it square
  cellheight = 15, # Adjust the cell height to make it square
  main = "Signature score of epiRG \n correlation with immune system",
  fontsize_row = 8, # Adjust font size for row labels
  fontsize_col = 8,  # Adjust font size for column labels
  breaks = seq(-1, 1, length.out = 101)
)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/Corr_imm_subtypes.pdf", width = 6, height = 5)
print(p)
dev.off()

library(pheatmap)

# Assuming 'dat' is your data matrix
max_val <- max(abs(xx), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 101)  # Create breaks from -max_val to max_val

# Define custom color palette with navy blue, white, and orange-red
custom_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

p1 <- pheatmap(xx,
               border_color = NA,
               cluster_cols = TRUE,
               cluster_rows = TRUE,
               treeheight_row = 0, 
               treeheight_col = 0,
               breaks = breaks,  # Use the dynamically created breaks
               color = custom_palette)  # Use the custom color palette

print(p1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/Corr_imm_subtypes.pdf", width = 6, height = 5)
print(p1)
dev.off()

#@@@@@@@@@@@@@@@@@@@@
  
[4.a]. Volcano plot (G): Association with Prolifiration
rm(list = ls())
library(ggplot2)
library(ggrepel)
inf1= "/home/u251079/BLCA_code/BLCA_ms_plots_data/4_cell_pro_score.txt"
dat1= read.table(inf1, sep = "\t", header = T, row.names = 1, quote = "")
rownames(dat1)
se= c("ARID1A_mut", "CHD6_mut","CHD7_mut","CREBBP_mut","EP300_mut","KDM6A_mut","CHD6_amp","CHD7_amp","PRDM9_amp","CHD3_del","CREBBP_del","HDAC4_del","PHF23_del")
dat1=dat1[se,]
names(dat1)
dat1$gene<-rownames(dat1)
#dat1 <- dat1[!rownames(dat1) %in% c("TP53_del", "TP53_delmut"),]
dat1$meandiff<-c(dat1$avg.MU-dat1$avg.WT)
dat1$sig_gene<-ifelse(dat1$negLogP > 1.30103, rownames(dat1), NA)

a <- ggplot(dat1, aes(log2FC, negLogP)) +
  geom_hline(yintercept = 1.30103, color="blue", linetype="dashed")+
  geom_point(size = 1.5, aes(color=sig_gene)) + # Label points with 'gene'
  theme_get()+
  theme(legend.position ="none", plot.title = element_text(hjust = 0.5, size = 8)) +
  labs(title = "Cell proliferation score", x = "Log2FC", y ="-log10(p-value)")+
  geom_label_repel(aes(label=rownames(dat1)), size=2)

# Print the plot
print(a)

a <- ggplot(dat1, aes(meandiff, negLogP, label=gene)) +
  geom_hline(yintercept = 1.30103, color="blue", linetype="dashed")+
  geom_vline(xintercept = 0, color="darkred", linetype="dashed")+
  geom_point(size = 1.5, aes(col=sig_gene)) + # Label points with 'gene'
  coord_cartesian(xlim = c(-0.15,0.15))+
  theme_classic()+
  geom_text_repel(size=4) +
  theme(legend.position ="none", plot.title = element_text(hjust = 0.5, size = 14, colour = "black")) +
  labs(title = "Cell proliferation score", x = "Average mean difference (mut-wild)", y ="-log10(p-value)")
#geom_label_repel(aes(label=sig_gene), size=2)

# Print the plot
print(a)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/Volcano_Proliferationimm.pdf", width = 4, height = 3.5)
print(a)
dev.off()




#----------
rm(list = ls())
library(ggplot2)
library(ggrepel)
library(ggpubr)
inf1= "/home/u251079/BLCA_code/BLCA_ms_plots_data/4_1_Leuco_frac_score.txt"
dat1= read.table(inf1, sep = "\t", header = T, row.names = 1, quote = "")
rownames(dat1)
se=c("TP53_mut","TP53_del" ,"TP53_delmut","CHD6_ampmut","CHD7_ampmut", "CREBBP_delmut", "EP300CREBBP_mut","EP300CREBBP_delmut")
dat1$gene<-rownames(dat1)
dat1 <- dat1[!rownames(dat1) %in% se,]
dat1$meandiff<-c(dat1$avg.MU-dat1$avg.WT)
dat1$sig_gene<-ifelse(dat1$negLogP > 2, rownames(dat1), NA)
dat1$val<-ifelse(dat1$gene %in% dat1$sig_gene, 1, 0)

a <- ggplot(dat1, aes(meandiff, negLogP, label = rownames(dat1))) +
  geom_hline(yintercept = 2, color = "blue", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "darkred", linetype = "dashed") +
  geom_point(size = 1.5, aes(col = factor(val))) + # Label points with 'gene'
  coord_cartesian(xlim = c(-0.15, 0.15)) +
  theme_classic() +
  geom_text_repel(size = 2) +
  scale_color_manual(values = c("0" = "black", "1" = "blue")) + # Define colors for 0 and 1
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 8, colour = "black")) +
  labs(title = "Leukocyte score", x = "Average mean difference (mut-wild)", y = "-log10(p-value)")#geom_label_repel(aes(label=sig_gene), size=2)

# Print the plot
print(a)

[4.2]. Volcano plot (G): Association with lymphocyte level
inf2= "/home/u251079/BLCA_code/BLCA_ms_plots_data/4_2_Lymph_inf_Sig_score.txt"
dat2= read.table(inf2, sep = "\t", header = T, row.names = 1, quote = "")
names(dat2)
dat2 <- dat2[!rownames(dat2) %in% se,]
dat2$meandiff<-c(dat2$avg.MU-dat2$avg.WT)
rownames(dat2)
dat2$sig_neg_log10p<- ifelse(dat2$negLogP > 1.5, rownames(dat2), NA)
dat2$val<-ifelse(dat2$sig_neg_log10p%in% rownames(dat2), 1, 0)

b <- ggplot(dat2, aes(meandiff, negLogP, label=rownames(dat2))) +
  geom_hline(yintercept = 1.5, color="blue", linetype="dashed")+
  geom_vline(xintercept = 0, color="darkred", linetype="dashed")+
  geom_point(size = 1.5, aes(color=as.factor(val))) + # Label points with 'gene'
  coord_cartesian(xlim = c(-1,1))+
  theme_classic()+
  geom_text_repel(size=2) +
  scale_color_manual(values = c("0" = "black", "1" = "blue"))+
  theme(legend.position ="none", plot.title = element_text(hjust = 0.5, size = 8)) +
  labs(title = "lymphocyte level", x = "Average mean difference (mut-wild)", y ="-log10(p-value)")
  #geom_label_repel(size=2)

# Print the plot
print(b)

c<-ggarrange(a, b, ncol=2, nrow=1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/4_1and2_volcano_plots_1013.pdf", width = 6, height = 3.5)
print(c)
dev.off()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.3]. Heatmap (G) → aberration x Imm 
# Mut vs Wild type
inp= "/home/u251079/BLCA_code/BLCA_ms_plots_data/4_3_abber_imm.csv"
data=read.csv(inp, header = T)
unique(data$name)

se= -log10(0.01)
rownames(data)
# Removed cells data
se1 <- grep("\\Cells", data$name)
data=data[-se1,]
# New jan 31 2024
se2<-c("Proliferation","Leukocyte.Fraction","Stromal.Fraction","Macrophage.Regulation","Lymphocyte.Infiltration.Signature.Score","TGF.beta.Response", "TCR.Richness", "BCR.Richness" )
data=data[data$name %in% se2,]
unique(data$name)
se1<- unique(data$name)[26:35]
data=data[!(data$name%in% se1),]

p <- ggplot(data, aes(x = reorder(info_colname, -negLogP), y = reorder(name, +negLogP))) +
  geom_tile(aes(fill = negLogP)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = se) +
  theme_classic() + xlab("")+ylab("")+
  labs(fill ="-log10(P-value)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.title = element_text(size=8, vjust = 0.1))

print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/Abbr_imm.pdf", width = 6, height = 3.5)
print(p)
dev.off()



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.4]. Boxplot (example)
rm(list = ls())
library(reshape2)
library(ggpubr)
myinf1 = "/home/u251079/BLCA_code/BLCA_ms_plots_data/4_3_Box.csv"

data_box<-read.csv(myinf1, row.names = 1, header = T)
names(data_box)
box_leuco<-data_box[,-2]
# remove non significant abbr from plots
se=c("TP53_mut","TP53_del" ,"TP53_delmut","ARID1A_mut", "CHD6_mut", "CHD7_mut", 
     "CREBBP_mut", "EP300_mut", "PIK3CA_mut", "RB1_mut" , "CREBBP_delmut", 
     "EP300CREBBP_mut", "EP300CREBBP_delmut","CHD6_ampmut","CHD7_ampmut")
box_leuco=box_leuco[,!(colnames(box_leuco) %in% se)]

box_leuco1<-melt(box_leuco, id.vars=c("Leukocyte.Fraction"))
box_leuco1<-box_leuco1[,c(3,2,1)]
names(box_leuco1)[1]<-"Mutation Status"
levels(box_leuco1$variable)

Leuco_box1<-ggplot(box_leuco1, aes(variable, Leukocyte.Fraction, fill=`Mutation Status`))+
  geom_boxplot(outlier.shape = NA,width=0.5, position=position_dodge(width = 0.7))+
  labs(title = NULL, x = NULL, y = "Leukocyte fraction") +
  coord_cartesian(ylim = c(min(box_leuco1$Leukocyte.Fraction),max(box_leuco1$Leukocyte.Fraction)+0.20))+
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     label.y = 0.85) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, color = "black"),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8, color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.text = element_text(angle = 0),
        legend.background = element_rect(size=0.3, 
                                         linetype="solid",
                                         colour ="black"))+
  scale_fill_discrete(labels=c('aberr', 'wild'),
                      guide = guide_legend(label.position = "top"))+
  coord_flip()

print(Leuco_box1)


#@@@@@@@@@@@@@ Lympho
box_lymp<-data_box[,-1]
se=c("TP53_mut","TP53_del" ,"TP53_delmut","CHD6_ampmut","CHD7_ampmut", "CREBBP_delmut", "EP300CREBBP_mut","EP300CREBBP_delmut")
box_lymp=box_lymp[,!colnames(box_lymp) %in% se]
se= c("ARID1A_mut", "CHD6_mut", "CHD7_mut", "CREBBP_mut","EP300_mut","PIK3CA_mut","CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut", "EP300CREBBP_delmut")
box_lymp<-box_lymp[, !(colnames(box_lymp) %in% se)]

box_lymp1<-melt(box_lymp, id.vars=c("Lymphocyte.Infiltration.Signature.Score"))
box_lymp1<-box_lymp1[,c(3,2,1)]
names(box_lymp1)[1]<-"Mutation Status"
names(box_lymp1)
levels(box_lymp1$variable)

lymp_box1<-ggplot(box_lymp1, aes(variable,Lymphocyte.Infiltration.Signature.Score,fill=`Mutation Status`))+
  geom_boxplot(outlier.shape = NA,width=0.5, position=position_dodge(width = 0.7))+
  labs(title = NULL, x = NULL, y = "Lymphocyte infiltration") +
  coord_cartesian(xlim = c(min(box_lymp1$Lymphocyte.Infiltration.Signature.Score),max(box_lymp1$Lymphocyte.Infiltration.Signature.Score)+2))+
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     label.y = 3.80) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, color = "black"),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none")+
  coord_flip()

print(lymp_box1)

library(patchwork)
p<-Leuco_box1+lymp_box1

print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/4_4_Boxplots.pdf", width = 8, height = 5)
print(p)
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.5]. Correlation matrix → signature x Imm
rm(list = ls())
# Using TIMER data set
mydat= "/home/u251079/BLCA_code/BLCA_ms_plots_data/Corr_Timer_TCGA_TIL_score.csv"
dat<-read.csv(mydat, row.names = 1, header = T)
dat<-t(dat)
p <- pheatmap(dat,cluster_cols = T, main ="Driver EpiRG aberration correlation \n with TIMER_representative")


library(pheatmap)

# Assuming 'dat' is your data matrix
max_val <- max(abs(dat), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 101)  # Create breaks from -max_val to max_val

# Define custom color palette with navy blue, white, and orange-red
custom_palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

p1 <- pheatmap(dat,
               border_color = NA,
               cluster_cols = TRUE,
               cluster_rows = TRUE,
               treeheight_row = 0, 
               treeheight_col = 0,
               breaks = breaks,  # Use the dynamically created breaks
               color = custom_palette)  # Use the custom color palette

print(p1)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/Cor_TCGA_BLCA_TIMER_cells.pdf", width = 6, height = 6)
print(p1)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.6]. ??? → with immune pathways
rm(list = ls())
inf="/home/u251079/BLCA_code/BLCA_ms_plots_data/4_3_immpath.csv"
corr_Imm_pathways=read.csv(inf, row.names = 1, header = T)
names(corr_Imm_pathways)
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
corr_Imm_pathways= corr_Imm_pathways[, !(colnames(corr_Imm_pathways)%in% se)]
dim(corr_Imm_pathways)
[1] 836  13

corr_Imm_pathways1<-corr_Imm_pathways
corr_Imm_pathways1$path<-rownames(corr_Imm_pathways1)
rownames(corr_Imm_pathways1)<-NULL
corr_data <- melt(corr_Imm_pathways1, id.vars = "path")

names(corr_data)
[1] "path"     "variable" "value"  
# Filter highly correlated values (greater than or equal to 0.8 or less than or equal to -0.8)
library(dplyr)
corr_data_filtered <- corr_data %>%
  filter(value >= 0.7 | value <= -0.7)
# threshold |0.8| only 47 pathways correlated
# threshold |0.7| only 157 pathways correlated

# After dim this will continue
se=unique(corr_data_filtered$path)
corr_Imm_pathways<-corr_Imm_pathways[se,]

corr_Imm_pathways$path<-rownames(corr_Imm_pathways)
rownames(corr_Imm_pathways)<-NULL
corr_data1 <- melt(corr_Imm_pathways, id.vars = "path")
se=unique(corr_data1$path)

# Find indices of names containing "CELL" in corr_data1$path
se1_indices <- grep("CELL_", corr_data1$path)
isolated_names <- corr_data1$path[se1_indices]
print(isolated_names)
corr_data2=corr_data1[c(corr_data1$path %in% isolated_names),]

# Create the heatmap plot with the filtered data and customized fill scale
heatmap_plot1 <- ggplot(corr_data2, aes(variable, path, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(-1, 1), breaks = c(-1, -0.7, 0.7, 1),
                       labels = c("-1.0", "-0.7", "0.7", "1.0"), name = "Correlation") +
  labs(title = "Highly Correlated Immune pathways (Threshold |0.7|)",
       x = NULL, y = NULL) +
  theme_pubclean() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
        plot.title = element_text(size = 10, color = "black"),
        legend.background = element_rect(size = 0.3, linetype = "solid", colour = "black"),
        legend.position = c(-1, 0.05))  # Set legend text to be vertical

# Print the heatmap with the legend scale and vertical legend labels
print(heatmap_plot1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/4_6_Imm_pathways_2.pdf", width = 12, height = 4.5)
print(heatmap_plot1)
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.7].with immune genes
rm(list = ls())
inf="/home/u251079/BLCA_code/BLCA_ms_plots_data/4_3_immgene.csv"
corr_Imm_gene= read.csv(inf, row.names = 1, header = T)
names(corr_Imm_gene)
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
corr_Imm_gene= corr_Imm_gene[, !(colnames(corr_Imm_gene)%in% se)]

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(dplyr)

corr_Imm_gene1<-corr_Imm_gene
corr_Imm_gene1$Imm_gene<-rownames(corr_Imm_gene1)
rownames(corr_Imm_gene1)<-NULL

corr_data <- melt(corr_Imm_gene1, id.vars = "Imm_gene")
names(corr_data)

corr_data_filtered <- corr_data %>% 
  filter(value >= 0.7 | value <= -0.7)
# Threshold |0.7| we have 207 Immune genes Highly correlated with ER signatures.
corr_data_filtered <- corr_data %>% 
  dplyr::filter(value >= 0.7 | value <= -0.7)

se=unique(corr_data_filtered$Imm_gene)

corr_Imm_gene<-corr_Imm_gene[se,]
corr_Imm_gene$Imm_gene<-rownames(corr_Imm_gene)
rownames(corr_Imm_gene)<-NULL

corr_data1 <- melt(corr_Imm_gene, id.vars = "Imm_gene")
names(corr_data)
unique(corr_data_filtered$Imm_gene)

# 171 immune genes highley correlated 

# Create the heatmap plot with the filtered data and customized fill scale
heatmap_plot1 <- ggplot(corr_data1, aes(variable, Imm_gene, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1),
                       breaks = c(-1, -0.7, 0, 0.7, 1),
                       labels = c("-1.0", "-0.7", "0.0", "0.7", "1.0"),
                       name = "Correlation") +
  labs(title = "Correlated with Immune genes (Threshold |0.7|)",
       x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.title = element_text(size = 10, color = "black", hjust = 0.5),
        legend.background = element_rect(size = 0.3, linetype = "solid", color = "black"))
#legend.position = c(0.9, 0.9))  # Set legend text to be vertical

# Print the heatmap with the updated color scale
print(heatmap_plot1)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/4_7_Imm_genes.pdf", width = 7, height = 25)
print(heatmap_plot1)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.8]. Volcano plot (G) → with immunotherapy response (Res vs Non-Res)
# https://sdgamboa.github.io/post/2020_volcano/
rm(list = ls())
in1<-"/home/u251079/BLCA_code/BLCA_ms_plots_data/Immune_res vs Non_res.csv"
inp<-read.csv(in1, row.names = 1, header = T)

inp$Glo_mean_diff<-inp$avg.MU-inp$avg.WT

# inp$Glo_abs_mean_log2FC<-log2(abs(inp$avg.MU))-log2(abs(inp$avg.WT))
inp$neg_log10p<- c(-log10(inp$P.t))
p_value<- -log10(0.05)
inp$sig<-ifelse(inp$Glo_mean_diff > 0 & inp$Glo_mean_diff < 0 | inp$neg_log10p >=p_value, rownames(inp), NA)
inp$val<-ifelse(inp$sig %in% rownames(inp), 1, 0)
  
names(inp)
rownames(inp)
se=c("TP53_mut","TP53_del" ,"TP53_delmut","CHD6_ampmut","CHD7_ampmut", "CREBBP_delmut", "EP300CREBBP_mut","EP300CREBBP_delmut")
inp <- inp[!rownames(inp) %in% se,]
  
library(ggrepel)

p1 <- ggplot(inp, aes(Glo_mean_diff, neg_log10p, label=sig)) + # -log10 conversion  
  geom_hline(yintercept = p_value, color="blue", linetype="dashed")+
  geom_vline(xintercept = 0, color="#E00A10", linetype="dashed")+
  geom_point(size = 1.5, aes(color=as.factor(val))) +
  theme_classic()+
  coord_cartesian(xlim = c(-12,14))+
  scale_color_manual(values = c("0" = "black", "1" = "red"))+
  labs(title = "ER genes with Immunotherapy response", 
       x = "Mean diffrence between Responders and Non-Responders", 
       y=expression("-log"[10]*"(p-value)"))+
  theme(legend.position ="none", 
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title.x = element_text(size = 8, color = "black")) +
  geom_text_repel(aes(label=sig), size=2)



print(p1)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/4_8_volcano_plot.pdf", width = 4, height =3)
print(p1)
dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[4.9]. Boxplot (E) → examples
rm(list=ls())
library(ggpubr)
library(reshape2)
data="/home/u251079/BLCA_code/BLCA_ms_plots_data/Box_Immune_res vs Non_res.csv"
data=read.csv(data, row.names = 1, header = T)
data1<-data
names(data1)
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
data1=data1[, !(colnames(data1) %in% se)]


# remove non significant 
se= c("CHD6_mut", "CHD7_mut", "CREBBP_mut", "EP300_mut",  "KDM6A_mut","PRDM9_amp")
data1=data[, ! (colnames(data)%in%se)]
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
data1=data1[, !(colnames(data1) %in% se)]

rownames(data1)<-NULL
data2<-melt(data1)
unique(data2$variable)

names(data2)
# add p-value
# https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/
p <- ggplot(data2, aes(variable, value, fill = response)) +
  geom_boxplot(outlier.shape = NA, colour = "black") +
  # scale_fill_manual(values = c("#DE676A", "#DDF68F")) +
  coord_cartesian(ylim = c(min(data2$value), max(data2$value) + 1)) +
  labs(title = "ER genes with Immunotherapy response", x = NULL, y = "ER_gene signature scores") +
  stat_compare_means(method = "wilcox.test", label = "p", 
                     label.y = max(data2$value) + 0.1) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color="black"),
        plot.title = element_text(hjust = 0.5, size = 10))+
  theme(legend.position = c(0.22,0.07),
        legend.direction = "horizontal")
  
# flip the plot
p <- ggplot(data2, aes(variable, value, fill = response)) +
  geom_boxplot(outlier.shape = NA, colour = "black", position = position_dodge(width = 0.9)) +
  coord_cartesian(xlim = c(min(data2$value), max(data2$value) + 1)) +
  labs(title = "ER genes with Immunotherapy response", x = NULL, y = "ER_gene signature scores") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     label.x = max(data2$value) + 0.1) +
  theme_classic() +
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  theme(legend.position = c(0.22, 0.07),
        legend.direction = "horizontal") +
  coord_flip()

print(p)



pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_4/Immune_responder_ER_genes_Boxplots.pdf", width = 3.2, height =4)
print(p)
dev.off() 



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-----------
Section-5:Mutations of some ER genes are associated with a global DNA-methylation change
#----------
[5.1] Volcano plot (G) → changes in median methylation levels
# t-test log10_p value vs median methylation levels.
rm(list = ls())
library(ggrepel)

res1="/home/u251079/BLCA_code/BLCA_ms_plots_data/Global_median_mean_CpG.csv"
res1=read.csv(res1, row.names = 1, header = T)
names(res1)
res1$neg_log10p<-c(-log10(res1$P.t))
res1$ER_genes<-row.names(res1)
res1$labs2_sig<-ifelse(res1$neg_log10p>=2, rownames(res1), NA)
res1$val<-ifelse(res1$labs2_sig %in% rownames(res1), 1, 0)
rownames(res1)
se=c("TP53_mut","TP53_del" ,"TP53_delmut","CHD6_ampmut","CHD7_ampmut", "CREBBP_delmut", "EP300CREBBP_mut","EP300CREBBP_delmut")
res1 <- res1[!rownames(res1) %in% se,]


Met_cp<-ggplot(res1, aes(Mean_Log2FC, neg_log10p, label=labs2_sig)) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 2, color = "#DC3A17", linetype = "dashed") +
  theme_classic()+
  geom_point(size = 1.5, aes(color=as.factor(val)))+
  scale_color_manual(values = c("0" = "black", "1" = "blue"))+
  coord_cartesian(xlim = c(-0.25, +0.25))+
  labs(title = "ER genes global median DNA methylation changes", x="Log2(Fold Change)", y="-Log10(p-value)")+
  theme(plot.title = element_text(size =12, hjust = 0.5, color = "black"),
        legend.position ="none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))+
  geom_text_repel(aes(label=labs2_sig), size=3, color="black") 

print(Met_cp)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/5_1_volcano.pdf", width = 5, height =4)
print(Met_cp)
dev.off() 


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[5.2] Boxplot (E) → KDM6A_mut  
rm(list = ls())
myinf1 = "/home/u251079/BLCA_code/BLCA_ms_plots_data/KDM6A_mut_CpG_mean_box.csv"
my_dat<-read.csv(myinf1, row.names = 1, header = T)
names(my_dat)

library(ggpubr)
p<-ggplot(my_dat, aes(x = KMD6A, y = KMD6A_mean_CpG, colour = KMD6A, shape=KMD6A)) +
  geom_boxplot(width=0.5, lwd=1,col="black")+ geom_jitter(width = 0.15)+ 
  labs(x = "KMD6A", y = "KMD6A_mean_CpG")+
  stat_compare_means(method="wilcox.test", label = "p.signif", label.x = 1.5)+
  theme_classic()+
  theme(legend.position="none",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"))

print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/5_2_KDM6A_mean_CpG.pdf", width = 4, height =4)
print(p)
dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[5.3] Boxplot (E) → KDM6A_mut  (CpG-4 categories)
#---------------------------------------------
inf="/home/u251079/BLCA_code/BLCA_ms_plots_data/5_3_mean_CpGs_types.csv"
data=read.csv(inf, row.names = 1, header = T)
data$KDM6A_mut<-ifelse(data$KM_mut>0, "Mut", "Wild")
data<- data[, c(1:4,6)]
library(reshape2)
molted=melt(data,id.vars=c("KDM6A_mut"))
library(ggpubr)
names(molted)
plot_CpGs<-ggplot(molted, aes(variable, value, color=KDM6A_mut))+
  geom_boxplot()+
  labs(x=NULL, y="mean CpG levels")+
  coord_cartesian(ylim =c(min(molted$value), max(molted$value)+0.01))+
  theme_classic2()+
  stat_compare_means(aes(group = KDM6A_mut),method = "wilcox.test",label = "p",
                     label.y = c(0.30, 0.75,0.52,0.75))+
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))+
  theme(legend.position = c(0.15, 0.85),legend.text = element_text(size=8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

print(plot_CpGs)
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/5_3_met_CpG.pdf", width = 4, height =4)
print(plot_CpGs)
dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
[5.4] Venn diagram (E) → KDM6A_mut (differential methylation CpG sites) (Total CpGs, Up CpG’s, Down CpG’s with KDM6A patient with Mut vs Wild): Ref: Fig 2A-PMC5839302
rm(list = ls())
mylist1="/home/u251079/BLCA_code/BLCA_ms_plots_data/Venn_diagram_KDM6A_met.csv"
mydat1<-read.csv(mylist1, row.names = 1, header = T)

library("ggVennDiagram")
class(mydat1)
mydat2<- as.list.data.frame(mydat1)
# CpG_lists<-mydat2[1:3]
which(mydat2$All_CpGs != mydat2$up)
names(mydat2)[3]<-"All_CpGs"
p<-ggVennDiagram(mydat2[1:3],
                 set_size = 4,
                 label= "count",edge_lty = "dashed",
                 edge_size = 1,label_color = "black",
                 label_alpha = 0)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  theme(legend.position = "NONE", plot.title = element_text(size =10, hjust = 0.5))

p1<-p + scale_x_continuous(expand = expansion(mult = .2))+ labs(title ="differential methylation CpGs sites of patients \n with KDM6A_mut vs wild")
print(p1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/5_4_Venn_KMD6A.pdf", width = 5, height =5, units = "in", res = 900)
print(p1)
dev.off() 
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
[5.5] Same as above → (CpG categories)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
mydat1<-read.csv(mylist1, row.names = 1, header = T)
se<-complete.cases(mydat1$up)
my_cpg<-mydat1[se,]
se2<-complete.cases(mydat1$down)
my_cpg2<-mydat1[se2,]
my_cat_cpg<-rbind(my_cpg, my_cpg2)
my_cat_cpg$HC<-ifelse(my_cat_cpg$HIL_CpG_class == "HC", my_cat_cpg$ID, NA)
my_cat_cpg$ IC<-ifelse(my_cat_cpg$HIL_CpG_class == " IC", my_cat_cpg$ID, NA)
my_cat_cpg$ICshore<-ifelse(my_cat_cpg$HIL_CpG_class == "ICshore", my_cat_cpg$ID, NA)
my_cat_cpg$LC<-ifelse(my_cat_cpg$HIL_CpG_class == "LC", my_cat_cpg$ID, NA)

my_cat_cpg1<-my_cat_cpg[, c(1:2,5:8)]
my_dat2<- as.list.data.frame(my_cat_cpg1)
class(my_dat2)

my_cat_cpg1 <- lapply(my_cat_cpg1, as.character)
names(my_cat_cpg1)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/5_5_Venn_KMD6A.pdf", width = 5, height =5, units = "in", res = 900)



venn<-ggVennDiagram(my_cat_cpg1, label= "count", 
                    edge_lty = "dashed",
                    edge_size = 1,label_color = "black",
                    label_alpha = 0)+
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  theme(legend.position = "NONE", plot.title = element_text(size =10, hjust = 0.5))

venn + scale_x_continuous(expand = expansion(mult = .2))+ labs(title ="Up_Down CpGs with categories of patients \n KDM6A_mut vs wild")

dev.off() 


#--------------------------------------------------------
[5.6] ?? → drug sensitivity a KMD6A mutation in CCLE/GDSC?
  #--------------------------------------------------------
rm(list = ls())
# Data
# /mount/ictr1/chenglab/cc59/PubDat/Dataset/CCLE
# CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf
# CCLE_Cell_Lines.txt
# CCLE_sample_info_file_2012-10-18.txt
# Cell_lines_annotations_20181226.txt
# isolate bladder
# CCLE_Drugs.txt
# CCLE_Expression_2012-09-29.res

# Drug sensitivity can be found at 
# /mount/ictr1/chenglab/cc59/PubDat/Dataset/CCLE/
in1= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/CCLE/CCLE_Cell_Lines.txt"
x= read.table(in1, sep = "\t", quote = "", stringsAsFactors = F)
in2= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/CCLE/CCLE_Drugs.txt"
x1= read.table(in2, sep = "\t", quote = "", stringsAsFactors = F)
in3= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/CCLE/CCLE_NP24.2009_Drug_data_2012.02.20.csv"
x2= read.csv(in3, header = T)

# More drug IC50 
# /mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS
done                                    gdsc_cell_lines_w2.csv
expU133A.txt                            gdsc_compounds_conc_w2.csv
GDSC1_fitted_dose_response_25Feb20.txt  GDSC_Drugs.txt
# permission denied GDSC2_2022Dec/mutations_summary_20221018.csv                          gdsc_en_input_w2_blood.txt
GDSC2_fitted_dose_response_25Feb20.txt  gdsc_en_input_w2.csv
GDSC_BreastCancerCellLine_profile.txt   gdsc_manova_input_w2.csv
GDSC_Cell_Lines.txt 

inp="/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/gdsc_en_input_w2.csv"
# inp1="/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/gdsc_cell_lines_w2.csv"
# inp2="/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/gdsc_compounds_conc_w2.csv"
# inp3="/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/gdsc_manova_input_w2.csv"


xx=read.csv(inp, row.names = 1, header = T)
xx_blader=xx[rownames(xx)=="bladder",]
se = which(xx_blader[1,]%in% 1)
xx1=xx[,se]

xx_blader1= xx[rownames(xx) == "KMD6A MUT",]


# xx1=read.csv(inp1, header = T) # COSMIC ID and the cell line name
# xx2=read.csv(inp2, row.names = 1, header = T) # min and max concentration micro molar of drugs
# xx3=read.csv(inp3, row.names = 1, header = T)

inp4= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/expU133A.txt"
inp5= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/GDSC1_fitted_dose_response_25Feb20.txt"
inp6= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/GDSC_Drugs.txt"
inp7= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/GDSC2_fitted_dose_response_25Feb20.txt" 
inp8= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/GDSC_BreastCancerCellLine_profile.txt"
inp9= "/mount/ictr1/chenglab/cc59/PubDat/Dataset/GDS/GDSC_Cell_Lines.txt"

xx4= read.table(inp4, sep = "\t", row.names = 1, header = T, quote = "", stringsAsFactors = F)
xx5= read.table(inp5, sep = "\t", header = T, quote = "", stringsAsFactors = F)
xx6= read.table(inp6, sep = "\t",  quote = "", stringsAsFactors = F)
xx7= read.table(inp7, sep = "\t",  header = T, quote = "", stringsAsFactors = F)
xx8= read.table(inp8, sep = "\t", row.names = 1, header = T, quote = "", stringsAsFactors = F)
xx9= read.table(inp9, sep = "\t",  quote = "", stringsAsFactors = F)


rm(x2)


# Added on January 16 2024
[6] association with MS (Mutational Signatures)
## SBS1	SBS2	SBS5	SBS13
Finally CREBBP-mut status of TCGA bladder patients revealed as highly association with MS (mutational signatures) from COSMIC.
## SBS13
https://cancer.sanger.ac.uk/signatures/sbs/sbs13/
  n1  n2  avg1  avg2 tscore         pt         pw
CREBBP_del          30 279 0.217 0.282 -2.271 0.02945401 0.02912755
CREBBP_delmut       61 248 0.233 0.286 -2.626 0.01013228 0.01048390

## SBS1 --> an endogenous mutational process initiated by spontaneous or enzymatic deamination of 5-methylcytosine to thymine which generates G
https://cancer.sanger.ac.uk/signatures/sbs/sbs1/
  n1  n2  avg1  avg2 tscore          pt          pw
CHD6_mut            28 281 0.044 0.072 -3.321 0.001950262 0.003584518
CHD7_mut            27 282 0.047 0.071 -2.317 0.026961134 0.004053723
EP300_mut           50 259 0.051 0.073 -2.700 0.008500867 0.002242213
KDM6A_mut           82 227 0.060 0.073 -1.851 0.065925508 0.060712493
EP300CREBBP_mut     80 229 0.059 0.073 -2.094 0.037853674 0.029267455
EP300CREBBP_delmut 100 209 0.060 0.074 -2.128 0.034304008 0.055406852
## --> mutation in these genes protect the spontaneous or enzymatic deamination of 5-methylcytosine to thymine

## SBS2
https://cancer.sanger.ac.uk/signatures/sbs/sbs2/
  n1  n2  avg1  avg2 tscore         pt         pw
CREBBP_mut          37 272 0.224 0.266 -2.185 0.03361776 0.03937055
CREBBP_delmut       61 248 0.233 0.268 -2.111 0.03733973 0.03913539

## SBS5
https://cancer.sanger.ac.uk/signatures/sbs/sbs5/
  n1  n2  avg1  avg2 tscore          pt          pw
CREBBP_del          30 279 0.453 0.321  3.148 0.003390694 0.001653681
CREBBP_delmut       61 248 0.406 0.316  3.035 0.003146955 0.001267301

# SBS2 + SBS13 (APOBEC)
n1  n2  avg1  avg2 tscore         pt          pw
CREBBP_delmut       61 248 0.466 0.554 -2.681 0.00866345 0.007183715


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ROC of epiRG-aber status to signature scores
# 1 Choi_data
rm(list = ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/BladderCancer/IntGen/data/Choi_GSE48277_GPL6947__GenomicEvent_iRAS.txt"
myinf2 = "/mount/ictr1//chenglab/cc59/PubDat/Cancer/Bladder/Choi_GSE48277/Clinical_info_GPL6947.txt" 

data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx


info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
info = info[!is.na(info$fgfr3.mutation), ]

comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
names(info)[2:5]


info<-info[,2:5]

names(data)
se<-grep("uni.noj__", names(data))
data<-data[,se]
colnames(data)<-gsub("uni.noj__","", colnames(data))

names(data)
se=c("FGFR3__MUT","RB1__MUT" , "TP53__MUT")
data1<-data[,se]

library(ROC)
library(survival)
# install.packages("pROC")
library(pROC)

# Create ROC curves 
roc_mut1 <- roc(response = info$p53.mutation, predictor = data1$TP53__MUT)
auc1<-auc(roc_mut1); print(auc1)

roc_mut2 <- roc(response = info$fgfr3.mutation, predictor = data1$FGFR3__MUT)
auc2<-auc(roc_mut2); print(auc2)

roc_mut3 <- roc(response = info$rb1.mutation, predictor = data1$RB1__MUT)
auc3<-auc(roc_mut3); print(auc3)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/ROC_Chao_data_test.pdf", width = 5, height =5)
par(pty="s")
plot.roc(roc_mut1, col = "black", lwd = 2, main = "Choi_data", legacy.axes = TRUE, 
         xlab = "1-Specificity (false positive)",
         ylab = "Sensitivity (true positive)")  
plot.roc(roc_mut2, add = TRUE, col = "blue", lwd = 2)
plot.roc(roc_mut3, add = TRUE, col = "red", lwd = 2)

legend("bottomright", legend = c("TP53_MUT (0.73)", "FGFR3_mut (0.67)", "RB1_mut (0.77)"),
       col = c("black","blue","red"), lwd = 2, box.col ="black", cex = 1)

dev.off() 


# 2. TCGA-BLCA
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/BladderCancer/IntGen/data/TCGA_BLCA__epiGene_iRAS.txt"
myinf3 = "/home/u251079/r_program/TCGA_BLCA_Freq_SomaticMutation_epiGenes2.txt"

data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx [,1:13]
names(data)
se=gsub("uni.noj__", "", names(data))
colnames(data)=se







info = read.table(myinf3, sep="\t", header=T, row.names=1, quote="")

comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]

info1=ifelse(info ==1, "mut", "wild")
colnames(info1)
se=names(data)
info1=info1[,se]
se=gsub("_mut", "_mutstatus", colnames(info1))
colnames(info1)=se
se=gsub("_amp", "_ampstatus", colnames(info1))
colnames(info1)=se
se=gsub("_del", "_delstatus", colnames(info1))
colnames(info1)=se
colnames(info1)

merged_data=cbind(data, info1)
colnames(merged_data)
# install.packages("ROC")
library(ROC)
library(survival)
# install.packages("pROC")
library(pROC)

# Create ROC curves "ARID1A_mut"       "CHD6_mut"         "CHD7_mut"         "CREBBP_mut"       "EP300_mut"        "KDM6A_mut" 

roc_mut1 <- roc(response = merged_data$ARID1A_mutstatus , predictor = merged_data$ARID1A_mut)
auc1<-auc(roc_mut1)
roc_mut2<- roc(response = merged_data$CHD6_mutstatus , predictor = merged_data$CHD6_mut)
auc2<-auc(roc_mut2)
roc_mut3<- roc(response = merged_data$CHD7_mutstatus , predictor = merged_data$CHD7_mut)
auc3<-auc(roc_mut3)

roc_mut4<- roc(response = merged_data$CREBBP_mutstatus , predictor = merged_data$CREBBP_mut)
auc4<-auc(roc_mut4)

roc_mut5<- roc(response = merged_data$EP300_mutstatus , predictor = merged_data$EP300_mut)
auc5<-auc(roc_mut5)

roc_mut6<- roc(response = merged_data$KDM6A_mutstatus , predictor = merged_data$KDM6A_mut)
auc6<-auc(roc_mut6)

# Add the ROC curve for CHD6_mut and display its AUC value
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/roc_TCGA-BLCA-mut.pdf", width = 6, height =6)
par(pty="s")
plot.roc(roc_mut1, col = "violet", lwd = 2, main = "TCGA-BLCA: on mutation status", legacy.axes=TRUE, 
         xlab = "1-Specificity (false positive)",
         ylab = "Sensitivity (true positive)")  # Adjust the print.auc.y value for the second AUC text
plot.roc(roc_mut2, add = TRUE, col = "blue", lwd = 2)
plot.roc(roc_mut3, add = TRUE, col = "green", lwd = 2)
plot.roc(roc_mut4, add = TRUE, col = "black", lwd = 2)
plot.roc(roc_mut5, add = TRUE, col = "orange", lwd = 2)
plot.roc(roc_mut6, add = TRUE, col = "red", lwd = 2)

# Add the legend
legend("bottomright", legend = c("ARID1A_mut (AUC = 0.78)", "CHD6_mut (AUC =0.90)", "CHD7_mut (AUC = 0.88)","CREBBP_mut (AUC = 0.74)","EP300_mut (AUC = 0.74)","KDM6A_mut (AUC = 0.71)" ),
       col = c( "violet", "blue","green","black", "orange","red"), lwd = 2, box.col ="black", cex = 0.9)
dev.off() 



# ROC curves for "CHD6_amp"         "CHD7_amp"         "PRDM9_amp"  
roc_mut1 <- roc(response = merged_data$CHD6_ampstatus , predictor = merged_data$CHD6_amp)
auc1<-auc(roc_mut1)
roc_mut2<- roc(response = merged_data$CHD7_ampstatus , predictor = merged_data$CHD7_amp)
auc2<-auc(roc_mut2)
roc_mut3<- roc(response = merged_data$PRDM9_ampstatus , predictor = merged_data$PRDM9_amp)
auc3<-auc(roc_mut3)
# ROC curves for "CHD3_del"         "CREBBP_del"       "HDAC4_del"        "PHF23_del"   
roc_mut4<- roc(response = merged_data$CHD3_delstatus , predictor = merged_data$CHD3_del)
auc4<-auc(roc_mut4)

roc_mut5<- roc(response = merged_data$CREBBP_delstatus , predictor = merged_data$CREBBP_del)
auc5<-auc(roc_mut5)

roc_mut6<- roc(response = merged_data$HDAC4_delstatus , predictor = merged_data$HDAC4_del)
auc6<-auc(roc_mut6)

roc_mut7<- roc(response = merged_data$PHF23_delstatus , predictor = merged_data$PHF23_del)
auc7<-auc(roc_mut7)

# Add the ROC curve for CHD6_mut and display its AUC value
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_2/roc_TCGA-BLCA-CNV.pdf", width = 6, height =6)
par(pty="s")
plot.roc(roc_mut1, col = "violet", lwd = 2, main = "TCGA-BLCA: on CNV status", legacy.axes=TRUE, 
         xlab = "1-Specificity (false positive)",
         ylab = "Sensitivity (true positive)")  # Adjust the print.auc.y value for the second AUC text
plot.roc(roc_mut2, add = TRUE, col = "blue", lwd = 2)
plot.roc(roc_mut3, add = TRUE, col = "green", lwd = 2)
plot.roc(roc_mut4, add = TRUE, col = "black", lwd = 2)
plot.roc(roc_mut5, add = TRUE, col = "orange", lwd = 2)
plot.roc(roc_mut6, add = TRUE, col = "red", lwd = 2)
plot.roc(roc_mut7, add = TRUE, col = "brown", lwd = 2)

# Add the legend
legend("bottomright", legend = c("CHD6_amp (AUC = 0.828)", "CHD7_amp (AUC = 0.826)", "PRDM9_amp (AUC = 0.832)" ,"CHD3_del (AUC = 0.831)", "CREBBP_del (AUC = 0.832)", "HDAC4_del (AUC = 0.883)","PHF23_del (AUC = 0.835)" ),
       col = c( "violet", "blue","green","black", "orange","red", "brown"), lwd = 2, box.col ="black", cex = 0.9)
dev.off() 




#&&&+===. KDM6a-mut methylation diffrence volcano plot
inf<-"/home/u251079/BLCA_code/BLCA_ms_plots_data/Volcano_KDM6A_met_up_down.csv"

kdm6a<-read.csv(inf, row.names = 1, header = T)
kdm6a<-na.omit(kdm6a)
library(ggplot2)
names(kdm6a)
kdm6a$log10p<-c(-log10(kdm6a$KDM6A_mut_Q_val_BH))

kdm6a <- kdm6a %>% 
  mutate(Expression = case_when(KDM6A_mut_mut_wt >= 0 & KDM6A_mut_Q_val_BH <= 0.01 ~ "Up-regulated",
                                KDM6A_mut_mut_wt <= -0 & KDM6A_mut_Q_val_BH <= 0.01 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))

legend_counts <- kdm6a %>%
  group_by(Expression) %>%
  summarise(Count = n())
View(legend_counts)


p1 <- ggplot(kdm6a, aes(x=KDM6A_mut_mut_wt, y=log10p))+
  geom_point(aes(color = Expression), size = 2/5) +
  geom_vline(xintercept = 0, color = "blue", linetype = "dashed") +
  geom_hline(yintercept =2, color = "#DC3A17", linetype = "dashed") +
  labs(title = "Diffrentially methylated CpGs of patients with KDM6A mutation to wild type",
       x="Difference between means (mut-wild)") +
  ylab(expression("-log"[10]*"(padj-BH)")) +
  theme_classic()+
  #  coord_cartesian( ylim = c(min(data$log10padj),max(data$log10padj)+5))+
  #xlim=c(-10.5,10.5),
  scale_color_manual(
    values = c("firebrick3", "gray50", "dodgerblue3"),
    labels = c(
      paste("Up-regulated (", legend_counts[legend_counts$Expression == "Up-regulated", "Count"], ")"),
      paste("NA"),
      paste("Down-regulated (", legend_counts[legend_counts$Expression == "Down-regulated", "Count"], ")")
    )
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = "top",
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6))

print(p1)


pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/Diff_meth_kdm6a.pdf", width = 6, height =5)
print(p1)
dev.off()




# Ep300 is missing in this data
# go with this data due to the EP300 is not a main resource in this session
mydat2="/mount/ictr1/chenglab/cc59/PubDat/Dataset/Firehose/Methylation/Gene_Promoter_Methyl/BLCA_Methy450KAvg_Promoter_beta.rda" 
load(mydat2)
met_data=mydata
se=grep("P30", rownames(met_data))
met_data2=met_data[se,]


met_data=as.data.frame(met_data)
se<-c("ARID1A","EP300","KDM6A","CHD6","CHD7","PRDM9","CHD3","CREBBP", "HDAC4","PHF23") 
met_data1=met_data[se,]
met_data1=na.omit(met_data1)

xx = colnames(met_data1)
se = which(substr(xx, 14,15)=="01")
met_data1 = met_data1[,se]
colnames(met_data1) = substr(colnames(met_data1), 1, 12)
met_data1<-t(met_data1)
met_data1<-as.data.frame(met_data1)
comxx = intersect(row.names(met_data1), row.names(score))
met_data1 = met_data1[comxx,]


cor_meth<-cor(met_data1, data, method = "spearman", use = "pair")

p<-pheatmap(cor_meth,cluster_rows =T, cluster_cols = T, angle_col = 45)
pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/Corr_DNA_meth.pdf", width = 5, height =3.5)
print(p)
dev.off()

[8.2] Correlation analysis using signature score to DNA methylation
rm(list=ls())
myinf1 = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/Firehose/Methylation/processed/BLCA_Methy450K.rda"
myinf2 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/BladderCancer/IntGen/data/TCGA_BLCA_Freq_SomaticMutation_epiGenes.txt"
myinf3 = "/mount/ictr1/chenglab/cc59/PubDat/organisms/human/annotation/GPL16304_illumina_HumanMethy450K_annotation.txt"


cpg = read.table(myinf3, sep="\t", header=T, quote="")
xx = as.character(cpg$HIL_CpG_class)
names(xx) = cpg$ID
cpg1 = xx

info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = which(info$KDM6A_mut>0)
sam.mu = row.names(info)[se]
se = which(info$KDM6A_mut==0)
sam.wt = row.names(info)[se]

load(myinf1)
data = mydata
xx = colnames(data)
# Assuming 'your_dataframe' is the name of your data frame
medians <- apply(data,2, median, na.rm=TRUE)
print(medians)
names(medians)
se = which(substr(names(medians), 14, 15)=="01")
medians<-medians[se]
names(medians)= substr(names(medians), 1, 12)

myinf2 = "/mount/ictr1/chenglab/cc59/WorSpa/m1_cancer/BladderCancer/IntGen/data/TCGA_BLCA__epiGene_iRAS.txt"
info <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="", check.names=F)
data<-info
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
score =xx
colnames(score)<-gsub("uni.noj__", "", colnames(score))

score=score[,1:13]


xx=intersect(rownames(score), names(medians))
score=score[xx,]
medians=medians[xx]


medians <- as.matrix(medians)
score <- as.matrix(score)

corr_meth<-cor(medians, score, method = "spearman")
plot(corr_meth)
corr_meth<-as.data.frame(corr_meth)
rownames(corr_meth)[1]<-"DNA Methylarion"
heatmap(corr_meth)

# Assuming 'corr' is your vector of correlation coefficients
corr <- c(ARID1A_mut = -0.11493086, CHD6_mut = -0.16913085, CHD7_mut = 0.08453035, 
          CREBBP_mut = 0.32817909, EP300_mut = 0.38738692, KDM6A_mut = -0.59146212, 
          CHD6_amp = -0.55269039, CHD7_amp = -0.54972258, PRDM9_amp = -0.07020391, 
          CHD3_del = -0.56200671, CREBBP_del = -0.45171140, HDAC4_del = -0.50857248, 
          PHF23_del = -0.55950554)

# Convert to a data frame for plotting
corr_df <- data.frame(Gene = names(corr), Correlation = corr)

# Order the data frame by correlation for a more intuitive visualization
corr_df <- corr_df[order(corr_df$Correlation), ]

# Plotting
library(ggplot2)
p<-ggplot(corr_df, aes(x = reorder(Gene, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() + # Flips the axes for better readability of gene names
  labs(title = "Correlation between epiRG-aber signature score and \n global DNA Methylation", 
       x = "Gene Aberration", y = "Spearman Correlation Coefficient") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", name="Correlation") +
  theme(legend.position = "right",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

# Note: Adjust 'scale_fill_gradient2' parameters as needed to fit your color preference.
print(p)

pdf("/home/u251079/BLCA_code/BLCA_ms_plots/res_sec_5/Corr_score_DNA_meth.pdf", width = 7, height =4)
print(p)
dev.off()





