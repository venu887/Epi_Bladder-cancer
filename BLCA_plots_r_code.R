#==========
Section-1: Alterations of epigenetic regulator gene aberrations in TCGA data   
#==========
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 1:Predominance of somatic mutations (SM) and copy number variations (CNV: indels) in epigenetic regulator genes (epiRG) across The Cancer Genome Atlas (TCGA) cancers
Figure 1A: plot somatic mutation rate of epigenes in all cancer types
library(ggplot2)
library(reshape)
library(reshape2)
library(patchwork)
# You can download the list of epigenetic regulator genes (epiRG) from the paper "Epigenetic modulators, modifiers and mediators in cancer aetiology and progression":Nature Reviews Genetics. 2016;17(5):284-99
# Then calculated the number of somatic mutation rate for each epiRG for all cancer types from TCGA project. 
myinf1 ="1A_Mutation_Num_All_categories_each_cancer_TCGA.csv"
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
# Save the figure in local directory



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 1B: Somatic Mutation (SM) and Copy Number Variaction (CNV) pattern pattern in BLCA
rm(list = ls())
# calculate the number of SM and CNV for each epiRG for all Bladder cancer patients 
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)
# Read and preprocess data
sm_cnv = "1B_BLCA_Epi_sm_cnv.csv" 
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
# Save the figure in local directory



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 1C: Volcano plot to show the prognostic association of them with epiRG mutation status of TCGA_BLCA patients
rm(list = ls())
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(ggplot2)
mycox_epi<-"1C_Epi_genes_survival_results.csv"

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

# Save the figure in local directory





#===========
Section-2: Differentiation and classification of driver epiRG aberration signature scores: 
# Gene signatures for driver epiRG aberrations
#===========
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 2.1: Schematic diagram
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
https://www.biorender.com

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 2B: Boxplot (Example) → The signatures can distinguish genes with mutations from without in Choi_data (GSE48075) 
rm(list = ls())
library(ggpubr)
library(rstatix)
inpt="2B_chao_box.csv"
c_box<-read.csv(inpt, header = T)
names(c_box)
table(c_box$Group~c_box$Gene)
p <- ggplot(c_box, aes(Gene, transformed_score, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.6, outlier.colour = "gray") +  # Add position_dodge to create a gap
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.4, color = "gray") +  # Add color to jitter points
  coord_cartesian(ylim = c(min(c_box$transformed_score), max(c_box$transformed_score) + 1)) +
  labs(title = "Choi_signatures", x = NULL, y = "Signature Score") +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.y = max(c_box$transformed_score) + 0.5) +
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


# Save the figure in local directory



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 2C: ROC of epiRG-aber status to signature scores
# 1 Choi_data
rm(list = ls())library(ROC)
library(survival)
library(pROC)

myinf1<-"2C_Figure.csv"
#You can download clinical information from GEO Choi_GSE48277_GPL6947
myinf2 = "Clinical_info_Choi_GSE48277_GPL6947.txt" 

data <- read.csv(myinf1,row.names=1, header = T)
info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
info = info[!is.na(info$fgfr3.mutation), ]

comSam = intersect(row.names(info), row.names(data))
data = data[comSam,]
info = info[comSam,]
names(info)[2:5]

info<-info[,2:5]
names(data)

# Create ROC curves 
roc_mut1 <- roc(response = info$p53.mutation, predictor = data$TP53__MUT)
auc1<-auc(roc_mut1); print(auc1)

roc_mut2 <- roc(response = info$fgfr3.mutation, predictor = data$FGFR3__MUT)
auc2<-auc(roc_mut2); print(auc2)

roc_mut3 <- roc(response = info$rb1.mutation, predictor = data$RB1__MUT)
auc3<-auc(roc_mut3); print(auc3)
# Save the figure in local directory
pdf("ROC_Chao_data_test.pdf", width = 5, height =5)
par(pty="s")
plot.roc(roc_mut1, col = "black", lwd = 2, main = "Choi_data", legacy.axes = TRUE, 
         xlab = "1-Specificity (false positive)",
         ylab = "Sensitivity (true positive)")  
plot.roc(roc_mut2, add = TRUE, col = "blue", lwd = 2)
plot.roc(roc_mut3, add = TRUE, col = "red", lwd = 2)

legend("bottomright", legend = c("TP53_MUT (0.73)", "FGFR3_mut (0.67)", "RB1_mut (0.77)"),
       col = c("black","blue","red"), lwd = 2, box.col ="black", cex = 1)
dev.off() 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 2D: Boxplot → CREBP and P300 signatures can cross- predict each others 
# Exp: Sig.score of CREBP significantly predict P300 Mut vs Wild
rm(list = ls())
in1="2D_Box_sec_TCGA.csv"
# Load the necessary library
library(gridExtra)
dat= read.csv(in1, row.names = 1, header = T)
colnames(dat)
dat1=as.data.frame(ifelse(dat[,7:12] == 1, "Mut", "Wild"))
dat1=cbind(dat[,1:6], dat1)
colnames(dat1)
# Sig.score of CREBBP significantly predict P300 Mut vs Wild
###transform scores If X is the current score, please use the following transformation:   X' = X/sd(abs(X))
# Assuming xx is your list of scores
xx <- dat1$Sig_CREBBP_mut
abs.xx <- abs(xx)
sd_abs_xx <- sd(abs.xx)
transformed_score <- xx / sd_abs_xx
max(transformed_score)
# Print the transformed scores
print(transformed_score)
dat1$Sig_CREBBP_mut_transformed_score<-transformed_score

p_value <- wilcox.test(Sig_CREBBP_mut_transformed_score ~ EP300_mut, data = dat1, alternative = "g")$p.value
p_value <-round(p_value,3)
p <- ggplot(dat1, aes(x = EP300_mut, y = Sig_CREBBP_mut_transformed_score, fill = EP300_mut)) +
  geom_boxplot(width = 0.7, lwd = 1, color = "black", notch = TRUE, notchwidth = 0.5) +  # Set black border here
  labs(x = "EP300", y = "Signature CREBBP") +
  ylim(min(dat1$Sig_CREBBP_mut_transformed_score), max(dat1$Sig_CREBBP_mut_transformed_score) + 0.5) +
  geom_text(aes(label = paste(signif(p_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_CREBBP_mut_transformed_score), size = 5, vjust = -1) +
  theme_classic() +
  scale_x_discrete(expand = c(0.25, 0.25)) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 10, color = "black"),
        text = element_text(family = "sans", face = "plain"))  # Set font style to plain

print(p)

# Sig.score of P300 significantly predict CREBBP Mut vs Wild
###transform scores If X is the current score, please use the following transformation:   X' = X/sd(abs(X))
# Assuming xx is your list of scores
xx <- dat1$Sig_EP300_mut
abs.xx <- abs(xx)
sd_abs_xx <- sd(abs.xx)
transformed_score <- xx / sd_abs_xx
max(transformed_score)
# Print the transformed scores
print(transformed_score)
dat1$Sig_EP300_mut_transformed_score<-transformed_score

p1_value <- wilcox.test(Sig_EP300_mut_transformed_score ~ EP300_mut, data = dat1, alternative = "g")$p.value
p1_value <-round(p1_value,9)
p1 <- ggplot(dat1, aes(x=CREBBP_mut, y=Sig_EP300_mut_transformed_score, fill=CREBBP_mut))+
  geom_boxplot(width = 0.7, lwd = 1, color = "black", notch = TRUE, notchwidth = 0.5) +  # Set black border here
  labs(x="CREBBP", y="Signature EP300")+
  ylim(min(dat1$Sig_EP300_mut_transformed_score), max(dat1$Sig_EP300_mut_transformed_score) + 0.5)+
  geom_text(aes(label = paste(signif(p1_value, digits = 4))),
            x = 1.5, y = max(dat1$Sig_EP300_mut_transformed_score), size = 5, vjust = -1)+
  theme_classic()+
  scale_x_discrete(expand = c(0.25, 0.25))+
  theme(axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        legend.position="none",
        plot.title = element_text(hjust = 0.5, size = 10))
print(p1)
library(patchwork)
p4<-p+ p1+ plot_annotation(title= "TCGA_BLCA", tag_levels = "I")+plot_layout(nrow = 1)
print(p4)
# Save the figure in local directory


#===============
Section-3:Association of driver gene signatures with prognosis using independent datasets
#===============
Figure 3A: Volcano plot (Global) X-HR vs. -log10(P) for all signatures 
rm(list = ls())
library(ggrepel)
library(tidyverse)
library(dplyr)
inf="3A_HR_with_Z_scores.csv"
my_cox_global=read.csv(inf, row.names = 1, header = T)
rownames(my_cox_global)
dim(my_cox_global)
my_cox_global<-my_cox_global[1:21,]

se= c("uni.noj__TP53_del","uni.noj__CHD6_ampmut","uni.noj__CHD7_ampmut", "uni.noj__CREBBP_delmut","uni.noj__EP300CREBBP_mut","uni.noj__EP300CREBBP_delmut", "uni.noj__TP53_delmut")
# Remove the rows with names matching the values in 'se'
my_cox_global <- my_cox_global[!(rownames(my_cox_global) %in% se), ]
rownames(my_cox_global)<-gsub("uni.noj__","", rownames(my_cox_global))
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


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 3B to 3I All KM plots
rm(list = ls())
library(survival)
library(ggplot2)
library(survminer)
library(gridExtra)
# Sjodahl_GSE32894_data_for_KM with 224 patient informatyion 
input1= "3B_3I_Sjodahl_GSE32894_data_for_KM.csv"
# GSE13507 data with 165 survival patients
input2= "3B_3I_GSE13507_data_for_KM.csv"

dat1= read.csv(input1, row.names = 1, header = T)
dat2= read.csv(input2, row.names = 1, header = T)

Figure SF2: Supplementary= KM plot (E) → EP300CREBBP_mut (in 2 datasets)
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

#Save the plot





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


#$$$$$$$$$$$$$$$$$$$$ GSE13507 data
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




#===========
Section-4: Evaluating the prognostic impact of driver epiRG aberration scores in bladder cancer using independent data sets 
#===========
rm(list = ls())
library(survival)
library(survminer)
# GSE13507_data 
inp="4_GSE13507_for_Box_Forest.csv"
info=read.csv(inp, row.names = 1, header=T)
colnames(info)<-gsub("uni.noj__", "", colnames(info))
###transform scores If X is the current score, please use the following transformation:   X' = X/sd(abs(X))
# Assuming xx is your list of scores
xx <- info$CREBBP_mut
abs.xx <- abs(xx)
sd_abs_xx <- sd(abs.xx)
transformed_score <- xx / sd_abs_xx
max(transformed_score)
# Print the transformed scores
print(transformed_score)
info$CREBBP_mut_transformed<-transformed_score

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 4A: Boxplot → association with invasiveness CREBBP (muscle-invasive vs. non-invasive)
# GSE13507_data 
se=as.factor(info$invasiveness)
table(se)

plot_violin<-ggplot(info, aes(as.factor(invasiveness), CREBBP_mut_transformed))+
  geom_violin(aes(fill=as.factor(invasiveness)), scale = "width")+
  geom_boxplot(width=0.3, lwd=0.5)+
  labs(title="GSE13507-Tumor Stages", x=NULL,y="Signature CREBBP_mut")+
  geom_jitter(position = position_dodge(width = 0.1), alpha = 0.2, color = "blue")+
  stat_summary(fun.y=median, geom="point", size=2, color="red")+
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  coord_cartesian(ylim = c(min(info$CREBBP_mut_transformed), max(info$CREBBP_mut_transformed)+1))+
  theme_bw()+
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))+
  scale_x_discrete(labels = c("Non Muscle-Inv \n (n=103)", "Muscle-Inv \n (n=62)"))

print(plot_violin)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 4B:Boxplot → association with Stages and signature scores of CREBBP
library(ggplot2)
table(info$T_stage, info$invasive)
# Create a factor with custom levels for T_stage
info$T_stage1 <- factor(info$T_stage, levels = c("Ta","T1", "T2", "T3", "T4"))
sample_counts <- info %>%
  group_by(T_stage1) %>%
  summarise(N = n())

plot_box<- ggplot(info, aes(T_stage1, CREBBP_mut_transformed, fill = factor(T_stage1))) +
  geom_boxplot() +
  labs(title = "GSE13507_Tumor Stages", x = NULL, y = "Signature CREBBP_mut") +
  coord_cartesian(ylim = c(min(info$CREBBP_mut_transformed), max(info$CREBBP_mut_transformed)+2)) +
  geom_hline(yintercept = mean(info$CREBBP_mut_transformed), linetype = 2) +
  stat_compare_means(method = "anova", label.y = 6.5) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.",
                     label.y = c(3.5, 5, 5.5, 5.7, 5)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 10))+
  scale_x_discrete(labels = c("Ta \n (n=24)","T1 \n (n=80)", "T2 \n (n=31)", "T3 \n (n=19)", "T4 \n (n=11)"))

print(plot_box)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 4C:Boxplot → CREBBP in different molecular subtypes
library(ggpubr)
inp1="4_GSE32894_for_BOX_Forest.csv"
info1=read.csv(inp1, header = T, row.names = 1)
colnames(info1)<-gsub("uni.noj__", "", colnames(info1))
xx <- info1$CREBBP_mut
abs.xx <- abs(xx)
sd_abs_xx <- sd(abs.xx)
transformed_score <- xx / sd_abs_xx
max(transformed_score)
# Print the transformed scores
print(transformed_score)
info1$CREBBP_mut_transformed<-transformed_score

# Reorder the factor levels of Mol_subtypes based on mean scores
mean_scores <- aggregate(CREBBP_mut_transformed ~ Mol_subtypes, data = info1, FUN = mean)
info1$Mol_subtypes <- factor(info1$Mol_subtypes, levels = mean_scores[order(mean_scores$CREBBP_mut_transformed), "Mol_subtypes"])

p3 <- ggboxplot(data = info1,
                x = "Mol_subtypes",
                y = "CREBBP_mut_transformed",
                color = "black",  # Set border color to black
                fill = "Mol_subtypes",   # Fill color based on molecular subtypes
                palette = "jco",
                legend = "none") +
  ylim(min(info1$CREBBP_mut_transformed), max(info1$CREBBP_mut_transformed) + 1) +
  geom_hline(yintercept = mean(info1$CREBBP_mut_transformed), linetype = 2, color = "black") + # Horizontal line for overall mean in black
  stat_summary(fun = "mean", geom = "line", color = "black", size = 1) +  # Add line indicating mean for each group
  stat_compare_means(method = "anova", label.y = 5) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.", label.y = c(3, 3, 4.2, 5.5, 5.5)) +
  labs(title = "GSE32894-Molecular Subtypes",x=NULL, y = "CREBBP_mut Sig.score") + 
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(hjust=1,vjust = 0.5, angle = 90))

print(p3)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 4D:Multivariate cox model of highly significant signature CREBBP-mut on stages and Molecular subtypes
rm(list = ls())
inp1="4_GSE32894_for_BOX_Forest.csv"
info1=read.csv(inp1, header = T, row.names = 1)
colnames(info1)<-gsub("uni.noj__", "", colnames(info1))
table(info1$tumor_stage)
info1$T_stage <- factor(info1$tumor_stage,
                        levels = c("Ta", "T1", "T2", "T2a", "T2b", "T3b", "T4a"),
                        labels = c("Pre early stage (Ta)", "Early Stage", "Early Stage","Early Stage", "Early Stage", "Late Stage", "Late Stage"))
names(info1)
se<-c("CREBBP_mut","t.dfs", "e.dfs","Mol_subtypes","stage", "age", "gender", "T_stage")
info1<-info1[se]
xx <- info1$CREBBP_mut
abs.xx <- abs(xx)
sd_abs_xx <- sd(abs.xx)
transformed_score <- xx / sd_abs_xx
max(transformed_score)
print(transformed_score)
info1$CREBBP_mut_trnsfm<-transformed_score
names(info1)
info1$Mol_subtypes <- factor(info1$Mol_subtypes)
info1$Mol_subtypes <- relevel(info1$Mol_subtypes, ref = "urothelial_A")
table(info1$stage, info1$e.dfs)
#@@@@@@@@@@@@@@@@@ Model-1
table(info1$stage)
fit.coxph2 <- coxph(Surv(t.dfs, e.dfs) ~ CREBBP_mut_trnsfm + T_stage +Mol_subtypes , data = info1)
summary(fit.coxph2)
coef_summary <- summary(fit.coxph2)$coefficients
coef_df <- as.data.frame(coef_summary)
xy <- c("CREBBP_mut (sig score)", "Early Stage (T1,T2)", "Late Stage (T3, T4)", 
        "genomically unstable", "infiltrated", "SCC-like", "urobasal B")
coef_df$Term <- xy
coef_df$Term <- factor(coef_df$Term, levels = rev(xy))
colnames(coef_df) <- c("coef", "exp(coef)", "se(coef)", "z", "p", "Term")
coef_df$lower_95 <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
coef_df$upper_95 <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
print(coef_df)
Tab_moltype_2<-coef_df
# Create the forest plot
p22<-ggplot(coef_df, aes(x = Term, y = `exp(coef)`, ymin = lower_95, ymax = upper_95)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10() +
  labs(x = "Variables",
       y = "Hazard Ratio (exp(coef))") +
  theme_minimal() +
  geom_text(aes(label = sprintf("p = %.3f", p)), 
            vjust = -0.5, 
            hjust = -0.1, 
            size = 3, 
            color = "black")

#@@@@@@@@@@@@@@@@@ Model-2
fit.coxph1 <- coxph(Surv(t.dfs, e.dfs) ~CREBBP_mut_trnsfm+Mol_subtypes, info1)
summary(fit.coxph1)
coef_summary <- summary(fit.coxph1)$coefficients
coef_df <- as.data.frame(coef_summary)
xy <- c("CREBBP_mut (sig score)", "genomically unstable", "infiltrated", "SCC-like", "urobasal B")
coef_df$Term <- xy
coef_df$Term <- factor(coef_df$Term, levels = rev(xy))
colnames(coef_df) <- c("coef", "exp(coef)", "se(coef)", "z", "p", "Term")
coef_df$lower_95 <- exp(coef_df$coef - 1.96 * coef_df$`se(coef)`)
coef_df$upper_95 <- exp(coef_df$coef + 1.96 * coef_df$`se(coef)`)
print(coef_df)
Tab_moltype_1<-coef_df
p11<-ggplot(coef_df, aes(x = Term, y = `exp(coef)`, ymin = lower_95, ymax = upper_95)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  scale_y_log10() +
  labs(x = NULL,
       y = "Hazard Ratio (exp(coef))") +
  theme_minimal() +
  geom_text(aes(label = sprintf("p = %.3f", p)), 
            vjust = -0.5, 
            hjust = -0.1, 
            size = 3, 
            color = "black")

print(p11)
ggarrange(p22+p11)




#==============
Section-5: Associations of epiRG aberration scores with cell proliferation and immune infiltration
#==============
Figure-5A: Correlation Heatmap signature scores verses Thorsson scores  
rm(list = ls())
#Download data from Thorsson et al.2018 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5982584/
myinf1 <- "Thorsson_2018_TCGA_immunelandscape.csv"
myinf2 = "TCGA_BLCA__epiGene_iRAS.txt"

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="", check.names=F)
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2 
# High score indicates highly association with mutation status, uncorrelated genes received small weights
colnames(xx) = gsub("\\.up", "", colnames(dat1))
score =xx 

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

#@@@@@@@@@@@@ Corrlation analysis only using Driver signatures 13 aberrations
score1=score[,1:13]
names(score1)
colnames(score1)<-gsub("uni.noj__","",colnames(score1))
comxx = intersect(row.names(score1), row.names(data1))
data1 = data1[comxx,]
score1 = score1[comxx,]

# calculation of z-scores for signatures score1
class(score1)
sig_mean = apply(score1, 2, mean, na.rm=T)
sig_sd = apply(score1, 2, sd, na.rm=T)
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

# calculation of z-scores for signatures data1
class(data1)
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

#Save the figure
library(pheatmap)
max_val <- max(abs(xx), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 101) 
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



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure: 5B Volcano plot (G): Association with Prolifiration
rm(list = ls())
library(ggplot2)
library(ggrepel)
inf1= "5B_cell_pro_score.txt"
dat1= read.table(inf1, sep = "\t", header = T, row.names = 1, quote = "")
rownames(dat1)
se= c("ARID1A_mut", "CHD6_mut","CHD7_mut","CREBBP_mut","EP300_mut","KDM6A_mut","CHD6_amp","CHD7_amp","PRDM9_amp","CHD3_del","CREBBP_del","HDAC4_del","PHF23_del")
dat1=dat1[se,]
names(dat1)
dat1$gene<-rownames(dat1)
dat1$meandiff<-c(dat1$avg.MU-dat1$avg.WT)
dat1$sig_gene<-ifelse(dat1$negLogP > 1.30103, rownames(dat1), NA)

a <- ggplot(dat1, aes(meandiff, negLogP, label=gene)) +
  geom_hline(yintercept = 1.30103, color="blue", linetype="dashed")+
  geom_vline(xintercept = 0, color="darkred", linetype="dashed")+
  geom_point(size = 1.5, aes(col=sig_gene)) + # Label points with 'gene'
  coord_cartesian(xlim = c(-0.15,0.15))+
  theme_classic()+
  geom_text_repel(size=4) +
  theme(legend.position ="none", plot.title = element_text(hjust = 0.5, size = 14, colour = "black")) +
  labs(title = "Cell proliferation score", x = "Average mean difference (mut-wild)", y ="-log10(p-value)")

# Print the plot
print(a)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 5B: Volcano plot (G): Association with Leukocyte score
rm(list = ls())
library(ggplot2)
library(ggrepel)
library(ggpubr)
inf1= "5B_1_Leuco_frac_score.txt"
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 5B: Volcano plot (G): Association with lymphocyte level
inf2= "5B_2_Lymph_inf_Sig_score.txt"
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
  geom_point(size = 1.5, aes(color=as.factor(val))) + 
  coord_cartesian(xlim = c(-1,1))+
  theme_classic()+
  geom_text_repel(size=2) +
  scale_color_manual(values = c("0" = "black", "1" = "blue"))+
  theme(legend.position ="none", plot.title = element_text(hjust = 0.5, size = 8)) +
  labs(title = "lymphocyte level", x = "Average mean difference (mut-wild)", y ="-log10(p-value)")
# Print the plot
print(b)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 5C: Heatmap → aberration score with Immune cells 
# Using TIMER data set
mydat= "5C_Corr_Timer_TCGA_TIL_score.csv"
dat<-read.csv(mydat, row.names = 1, header = T)
dat<-t(dat)
max_val <- max(abs(dat), na.rm = TRUE)
breaks <- seq(-max_val, max_val, length.out = 101)  
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



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 5D: Prediction of signature scores with patients with Immunotherapy responces
library(ggpubr)
library(reshape2)
data="5D_Box_Immune_res vs Non_res.csv"
data=read.csv(data, row.names = 1, header = T)
data1<-data
names(data1)
# remove non significant 
se= c("CHD6_mut", "CHD7_mut", "CREBBP_mut", "EP300_mut",  "KDM6A_mut","PRDM9_amp")
data1=data[, ! (colnames(data)%in%se)]
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
data1=data1[, !(colnames(data1) %in% se)]
rownames(data1)<-NULL
data2<-melt(data1)
unique(data2$variable)
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


#===============
Section-6:Mutations of some ER genes are associated with a global DNA-methylation change
#===============
rm(list=ls())

Figure 6A: Correlation analysis using signature score to DNA methylation
myinf1 = "BLCA_Methy450K.rda" # Download from the Firebrowse http://firebrowse.org
myinf2 = "TCGA_BLCA_Freq_SomaticMutation_epiGenes.txt" # Download from the Firebrowse http://firebrowse.org
myinf3 = "GPL16304_illumina_HumanMethy450K_annotation.txt" # Download from the Firebrowse http://firebrowse.org
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
medians <- apply(data,2, median, na.rm=TRUE)
print(medians)
names(medians)
se = which(substr(names(medians), 14, 15)=="01")
medians<-medians[se]
names(medians)= substr(names(medians), 1, 12)
myinf2 = "TCGA_BLCA__epiGene_iRAS.txt"
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
corr <- c(ARID1A_mut = -0.11493086, CHD6_mut = -0.16913085, CHD7_mut = 0.08453035, 
          CREBBP_mut = 0.32817909, EP300_mut = 0.38738692, KDM6A_mut = -0.59146212, 
          CHD6_amp = -0.55269039, CHD7_amp = -0.54972258, PRDM9_amp = -0.07020391, 
          CHD3_del = -0.56200671, CREBBP_del = -0.45171140, HDAC4_del = -0.50857248, 
          PHF23_del = -0.55950554)
corr_df <- data.frame(Gene = names(corr), Correlation = corr)
corr_df <- corr_df[order(corr_df$Correlation), ]
library(ggplot2)
p<-ggplot(corr_df, aes(x = reorder(Gene, Correlation), y = Correlation, fill = Correlation)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() + 
  labs(title = "Correlation between epiRG-aber signature score and \n global DNA Methylation", 
       x = "Gene Aberration", y = "Spearman Correlation Coefficient") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", name="Correlation") +
  theme(legend.position = "right",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))
print(p)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 6B: Volcano plot of changes in median methylation levels
rm(list = ls())
library(ggrepel)
res1="6B_Global_median_mean_CpG.csv"
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 6C: KDM6a-mut methylation Up and Down CpGs
inf<-"6C_Volcano_KDM6A_met_up_down.csv"
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 6D: Average methylation levels in patients who have mutations in the KDM6A gene (KDM6A-mut)
rm(list = ls())
myinf1 = "6D_KDM6A_mut_CpG_mean_box.csv"
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


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure 6E: Bladder cancer patients with KDM6A_mut to wild type CpG levels (CpG-4 categories)
inf="6E_Kdm6a_mean_CpGs_types.csv"
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




#==================
Supplementry 
#==================
Figure SF1: ROC curves and the AUC of aberrations (SM/CNV ≥ 30 patients) epiG signature scores on patient aberration status prediction. 
rm(list=ls())
myinf1 = "TCGA_BLCA__epiGene_iRAS.txt"
myinf3 = "TCGA_BLCA_Freq_SomaticMutation_epiGenes2.txt" # Somatic mutation data of epiRG 
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
library(ROC)
library(survival)
library(pROC)
# Create ROC curves SM
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
pdf("roc_TCGA-BLCA-mut.pdf", width = 6, height =6)
par(pty="s")
plot.roc(roc_mut1, col = "violet", lwd = 2, main = "TCGA-BLCA: on mutation status", legacy.axes=TRUE, 
         xlab = "1-Specificity (false positive)",
         ylab = "Sensitivity (true positive)")  # Adjust the print.auc.y value for the second AUC text
plot.roc(roc_mut2, add = TRUE, col = "blue", lwd = 2)
plot.roc(roc_mut3, add = TRUE, col = "green", lwd = 2)
plot.roc(roc_mut4, add = TRUE, col = "black", lwd = 2)
plot.roc(roc_mut5, add = TRUE, col = "orange", lwd = 2)
plot.roc(roc_mut6, add = TRUE, col = "red", lwd = 2)
legend("bottomright", legend = c("ARID1A_mut (AUC = 0.78)", "CHD6_mut (AUC =0.90)", "CHD7_mut (AUC = 0.88)","CREBBP_mut (AUC = 0.74)","EP300_mut (AUC = 0.74)","KDM6A_mut (AUC = 0.71)" ),
       col = c( "violet", "blue","green","black", "orange","red"), lwd = 2, box.col ="black", cex = 0.9)
dev.off() 

# ROC curves for CNV
roc_mut1 <- roc(response = merged_data$CHD6_ampstatus , predictor = merged_data$CHD6_amp)
auc1<-auc(roc_mut1)
roc_mut2<- roc(response = merged_data$CHD7_ampstatus , predictor = merged_data$CHD7_amp)
auc2<-auc(roc_mut2)
roc_mut3<- roc(response = merged_data$PRDM9_ampstatus , predictor = merged_data$PRDM9_amp)
auc3<-auc(roc_mut3)
roc_mut4<- roc(response = merged_data$CHD3_delstatus , predictor = merged_data$CHD3_del)
auc4<-auc(roc_mut4)
roc_mut5<- roc(response = merged_data$CREBBP_delstatus , predictor = merged_data$CREBBP_del)
auc5<-auc(roc_mut5)
roc_mut6<- roc(response = merged_data$HDAC4_delstatus , predictor = merged_data$HDAC4_del)
auc6<-auc(roc_mut6)
roc_mut7<- roc(response = merged_data$PHF23_delstatus , predictor = merged_data$PHF23_del)
auc7<-auc(roc_mut7)
# Add the ROC curve for CHD6_mut and display its AUC value
pdf("roc_TCGA-BLCA-CNV.pdf", width = 6, height =6)
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
legend("bottomright", legend = c("CHD6_amp (AUC = 0.828)", "CHD7_amp (AUC = 0.826)", "PRDM9_amp (AUC = 0.832)" ,"CHD3_del (AUC = 0.831)", "CREBBP_del (AUC = 0.832)", "HDAC4_del (AUC = 0.883)","PHF23_del (AUC = 0.835)" ),
       col = c( "violet", "blue","green","black", "orange","red", "brown"), lwd = 2, box.col ="black", cex = 0.9)
dev.off() 



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure SF2: Prognostic prediction of combination of CREBBP mutation status with EP300 mutation and developed signature scores of EP300CREBBP-mut.
1. Volcano Plot Refer to Figure 3A 
2. KM plots # GSE32894 data
Supplementary= KM plot (E) → EP300CREBBP_mut (in 2 datasets)
# HR and p-value was calculated from Cox model
                 GSE13507_HR GSE32894_HR
EP300CREBBP_mut   1.8707262   1.8943059
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
#@@@@@@@@@@@@@@@@@@@
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure SF3:  Boxplot → CREBBP in different stages With Sjodahl_GSE32894_data of 224 patient 
# data can used where we used Forest plots Figure 4D
table(info1$stage)
info1$stage<-ifelse(info1$stage == "T4", "T3", info1$stage)
info1$stage1 <- factor(info1$stage, levels = c("Ta","T1", "T2", "T3"))
sample_counts <- info1 %>%
  group_by(stage1) %>%
  summarise(N = n())
plot_box1<- ggplot(info1, aes(stage1, CREBBP_mut_transformed, fill = factor(stage1))) +
  geom_boxplot() +
  labs(title = "GSE32894", x = NULL, y = "Signature CREBBP_mut") +
  coord_cartesian(ylim = c(min(info1$CREBBP_mut_transformed), max(info1$CREBBP_mut_transformed) + 0.5)) +
  geom_hline(yintercept = mean(info1$CREBBP_mut_transformed), linetype = 2) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", ref.group = ".all.",
                     label.y = c(4, 4, 5.5, 5.5)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 10))+
  scale_x_discrete(labels = c("Ta (n=110)","T1 (n=63)", "T2 (n=43)", "T3 (n=8)"))

print(plot_box1)
# Save the plot

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure SF4: correlation signature scores with immune pathways
rm(list = ls())
library(dplyr)
inf="SF4_immpath.csv"
corr_Imm_pathways=read.csv(inf, row.names = 1, header = T)
names(corr_Imm_pathways)
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
corr_Imm_pathways= corr_Imm_pathways[, !(colnames(corr_Imm_pathways)%in% se)]
corr_Imm_pathways1<-corr_Imm_pathways
corr_Imm_pathways1$path<-rownames(corr_Imm_pathways1)
rownames(corr_Imm_pathways1)<-NULL
corr_data <- melt(corr_Imm_pathways1, id.vars = "path")
corr_data_filtered <- corr_data %>%
  filter(value >= 0.7 | value <= -0.7)

se=unique(corr_data_filtered$path)
corr_Imm_pathways<-corr_Imm_pathways[se,]
corr_Imm_pathways$path<-rownames(corr_Imm_pathways)
rownames(corr_Imm_pathways)<-NULL
corr_data1 <- melt(corr_Imm_pathways, id.vars = "path")
se=unique(corr_data1$path)

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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Figure SF5: correlation of signature scores with immune genes
rm(list = ls())
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(reshape2)
library(dplyr)
inf="SF5_immgene.csv"
corr_Imm_gene= read.csv(inf, row.names = 1, header = T)
names(corr_Imm_gene)
se=c("CHD6_ampmut","CHD7_ampmut","CREBBP_delmut","EP300CREBBP_mut","EP300CREBBP_delmut")
corr_Imm_gene= corr_Imm_gene[, !(colnames(corr_Imm_gene)%in% se)]
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

print(heatmap_plot1)

