library(tidyverse)
library(ggplot2)
library(ggrepel)
options(scipen=999)
source('0_data/ggplot_theme_Publication-2.R')


# Load data ####
# GO-based metrics
aspects <- c("F","P","C")
for (identifier in c("All","HTP","LTP")){
  combined <- readRDS(paste0("1_metrics/GO_exp_aspect_",identifier,".RDS"))
  # Change name of dataframe into prot_df_{identifier}_{aspect}
  for (aspect in aspects){
    assign(paste0("prot_df_",identifier,"_",aspect), combined[[aspect]])
  }
}
wider_format <- function(df, column){
  df_byyear <- df%>%
    select(1,2, !!column)%>%
    pivot_wider(id_cols = Genes, names_from = Date, values_from = !!column)
  return(df_byyear)
}
icF <- wider_format(prot_df_All_F, "IA")
icP <- wider_format(prot_df_All_P, "IA")
icC <- wider_format(prot_df_All_C, "IA")
icF = icF%>%mutate(gain = `2022` - `2013`)%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0) | gain < 0, 0, gain))%>%
  filter(`2022` > 0) # remove proteins that are currently fake
icP = icP%>%mutate(gain = `2022` - `2013`)%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0) | gain < 0, 0, gain))%>%
  filter(`2022` > 0) # remove proteins that are currently fake
icC = icC%>%mutate(gain = `2022` - `2013`)%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0) | gain < 0, 0, gain))%>%
  filter(`2022` > 0) # remove proteins that are currently fake
# Interest metric
fc = readRDS(file="1_metrics/interest_FPE.RDS")
# Combine
fc_ic <- full_join(fc[,c(1:4,13:15)],icF[,c(1,2,11,12)],by="Genes")%>%
  full_join(.,icP[,c(1,2,11,12)],by="Genes")%>%
  full_join(.,icC[,c(1,2,11,12)],by="Genes")%>%
  rename_with(~c("STRING","Genes","Common","FPE_13","FPE_22","relgain_FPE","gain_FPE","MFO_13","MFO_22","gain_MFO", "BPO_13","BPO_22", "gain_BPO","CCO_13","CCO_22", "gain_CCO"))

# Helper functions
cor_summary <- function(x,y){
  test <- cor.test(x,y, method="spearman",use="pairwise.complete.obs",quiet=TRUE)
  valid_pairs <- complete.cases(x,y)
  tibble(
    rho = test$estimate,
    pval = test$p.value,
    num = sum(valid_pairs)
  )
}


# Disease-gene associations ####
dis <- read.csv(file="0_data/ext_data/human_disease_knowledge_filtered.tsv",header=F,sep="\t")
colnames(dis) <- c("STRING","Common","DO_ID","Disease_name","Database","Curated","Confidence_score")
dis2 <- dis[,c("STRING","Common","Disease_name")]%>%
  distinct(.)
dis_fc_ic <- full_join(dis2, fc_ic, by = "STRING")%>%
  rowwise()%>%
  select(-c(Common.x,Common.y))
dis_fc_ic$Disease_name[is.na(dis_fc_ic$Disease_name)] <- "No known associations"
filtered_dis_fc_ic <- dis_fc_ic%>%
  filter(!is.na(Genes))%>%
  distinct(Disease_name,Genes, .keep_all = T)

# Preprocess ####
# Filter to >= 20-gene diseases and then sort by correlation coefficient with adjusted p-value cutoff
# MFO
bh_MFO <- filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%
  select(-c(6,7,11,12,13))%>%
  filter(!is.na(FPE_13) & !is.na(gain_MFO))%>%
  group_by(Disease_name)%>%
  mutate(num=n(), .keep = "all")%>%
  filter(num >= 20)%>%
  group_by(Disease_name,num)%>%
  summarise(
    cor_FPE13_gainMFO_rho = cor.test(FPE_13,gain_MFO,method="spearman")$estimate,
    cor_FPE13_gainMFO_pval = cor.test(FPE_13,gain_MFO,method="spearman")$p.value,
  )
bh_MFO$p.adjust <- p.adjust(bh_MFO$cor_FPE13_gainMFO_pval,method="BH")
bh_MFO$p.adjust2 <- p.adjust(bh_MFO$cor_FPE13_gainMFO_pval,method="bonferroni")
bh_MFO <- bh_MFO %>%
  mutate(Disease_name = ifelse(p.adjust < 0.05, paste0(Disease_name, "(*)"), Disease_name))%>%
  arrange(cor_FPE13_gainMFO_rho)
# BPO
bh_BPO <- filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%
  select(-c(6,7,8,9,10))%>%
  filter(!is.na(FPE_13) & !is.na(gain_BPO))%>%
  group_by(Disease_name)%>%
  mutate(num=n(), .keep = "all")%>%
  filter(num >= 20)%>%
  group_by(Disease_name,num)%>%
  summarise(
    cor_FPE13_gainBPO_rho = cor.test(FPE_13,gain_BPO,method="spearman")$estimate,
    cor_FPE13_gainBPO_pval = cor.test(FPE_13,gain_BPO,method="spearman")$p.value,
  )
bh_BPO$p.adjust <- p.adjust(bh_BPO$cor_FPE13_gainBPO_pval,method="BH")
bh_BPO$p.adjust2 <- p.adjust(bh_BPO$cor_FPE13_gainBPO_pval,method="bonferroni")
bh_BPO <- bh_BPO %>%
  mutate(Disease_name = ifelse(p.adjust < 0.05, paste0(Disease_name, "(*)"), Disease_name))%>%
  arrange(cor_FPE13_gainBPO_rho)

# Option 1: Geom_point ####
# MFO
bg_all <- cor_summary(filtered_dis_fc_ic$FPE_13,filtered_dis_fc_ic$gain_MFO)
bg_nodisease <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(FPE_13)%>%unlist(),
                            filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(gain_MFO)%>%unlist())
bg_alldisease <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(gain_MFO)%>%unlist())
bg_bigdisease <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name %in% bh_MFO$Disease_name)%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name %in% bh_MFO$Disease_name)%>%select(gain_MFO)%>%unlist())
ggplot(data=bh_MFO,aes(y=cor_FPE13_gainMFO_rho, x=fct_inorder(Disease_name))) + geom_point(shape = ifelse(bh_MFO$p.adjust < 0.05,"triangle","circle"),size=ifelse(bh_MFO$p.adjust < 0.05,5,3)) + 
  theme_Publication() + 
  theme(axis.text.x = element_text(angle = 75,hjust = 1)) + 
  xlab("Disease name") + ylab("Spearman coefficient") + 
  geom_hline(yintercept = bg_all$rho,colour = "black",show.legend = TRUE) + 
  geom_hline(yintercept = bg_nodisease$rho, colour="grey", show.legend = TRUE) + 
  annotate("text", x = Inf, y = bg_all$rho, label = "All proteins", vjust = -0.5, hjust = 1.1, colour = "black", size = 5) + 
  annotate("text", x = Inf, y = bg_nodisease$rho, label = "Proteins without disease associations", vjust = -0.5, hjust = 1.1, colour = "grey",size=5)

# BPO 
bg_all <- cor_summary(filtered_dis_fc_ic$FPE_13,filtered_dis_fc_ic$gain_BPO)
bg_nodisease <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(FPE_13)%>%unlist(),
                            filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(gain_BPO)%>%unlist())
bg_alldisease <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(gain_BPO)%>%unlist())
bg_bigdisease <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name %in% bh_BPO$Disease_name)%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name %in% bh_BPO$Disease_name)%>%select(gain_BPO)%>%unlist())
ggplot(data=bh_BPO,aes(y=cor_FPE13_gainBPO_rho, x=fct_inorder(Disease_name))) + geom_point(shape = ifelse(bh_BPO$p.adjust < 0.05,"triangle","circle"),size=ifelse(bh_BPO$p.adjust < 0.05,5,3)) + 
  theme_Publication() + 
  theme(axis.text.x = element_text(angle = 75,hjust = 1)) + 
  xlab("Disease name") + ylab("Spearman coefficient") + 
  geom_hline(yintercept = bg_all$rho,colour = "black",show.legend = TRUE) + 
  geom_hline(yintercept = bg_nodisease$rho, colour="grey", show.legend = TRUE) + 
  annotate("text", x = Inf, y = bg_all$rho, label = "All proteins", vjust = -0.5, hjust = 1.1, colour = "black", size = 5) + 
  annotate("text", x = Inf, y = bg_nodisease$rho, label = "Proteins without disease associations", vjust = -0.5, hjust = 1.1, colour = "grey", size = 5)
# Option 2: Boxplot ####
# MFO
bg_all_MFO <- cor_summary(filtered_dis_fc_ic$FPE_13,filtered_dis_fc_ic$gain_MFO)
bg_nodisease_MFO <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(FPE_13)%>%unlist(),
                            filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(gain_MFO)%>%unlist())
bg_alldisease_MFO <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(gain_MFO)%>%unlist())
bg_bigdisease_MFO <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name %in% bh_MFO$Disease_name)%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name %in% bh_MFO$Disease_name)%>%select(gain_MFO)%>%unlist())
pMFO <- ggplot(data = bh_MFO, aes(x = "", y = cor_FPE13_gainMFO_rho)) + 
  geom_boxplot(outlier.shape = NA, width=0.3) + 
  geom_jitter(shape = ifelse(bh_MFO$p.adjust < 0.05, 24, 21), size = ifelse(bh_MFO$p.adjust < 0.05, 5, 3),width = 0.03,fill="black") +
  geom_text_repel(data = bh_MFO %>% filter(cor_FPE13_gainMFO_rho > 0.4),
                  aes(label = Disease_name), size = 5,nudge_x = -0.25) +
  xlab(NULL) + ylab("Spearman coefficient") + theme_Publication() + 
  geom_hline(yintercept = bg_all$rho, colour = "black", show.legend = TRUE) + 
  geom_hline(yintercept = bg_nodisease$rho, colour = "grey", show.legend = TRUE) + 
  annotate("text", x = Inf, y = bg_all$rho, label = paste0("All proteins: ",round(bg_all$rho,3)), vjust = -0.5, hjust = 1.1, colour = "black",size=5) + 
  annotate("text", x = Inf, y = bg_nodisease$rho, label = paste0("Proteins without disease\nassociations: ",round(bg_nodisease$rho,3)), vjust = 1.1, hjust = 1.1, colour = "grey",size=5) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), text = element_text(size=20,family = "sans"), axis.text = element_text(size=20,family = "sans")) +
  coord_cartesian(ylim = c(-0.35,0.6))
# BPO
bg_all_BPO <- cor_summary(filtered_dis_fc_ic$FPE_13,filtered_dis_fc_ic$gain_BPO)
bg_nodisease_BPO <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(FPE_13)%>%unlist(),
                            filtered_dis_fc_ic%>%filter(Disease_name == "No known associations")%>%select(gain_BPO)%>%unlist())
bg_alldisease_BPO <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name != "No known associations")%>%select(gain_BPO)%>%unlist())
bg_bigdisease_BPO <- cor_summary(filtered_dis_fc_ic%>%filter(Disease_name %in% bh_BPO$Disease_name)%>%select(FPE_13)%>%unlist(),
                             filtered_dis_fc_ic%>%filter(Disease_name %in% bh_BPO$Disease_name)%>%select(gain_BPO)%>%unlist())
pBPO <- ggplot(data = bh_BPO, aes(x = "", y = cor_FPE13_gainBPO_rho)) + 
  geom_boxplot(outlier.shape = NA, width=0.3) + 
  geom_jitter(shape = ifelse(bh_BPO$p.adjust < 0.05, 24, 21), size = ifelse(bh_BPO$p.adjust < 0.05, 5, 3),width = 0.03,fill="black") +
  geom_text_repel(data = bh_BPO %>% filter(cor_FPE13_gainBPO_rho > 0.4),
                  aes(label = Disease_name), size = 5,nudge_x = 0.25) +
  
  xlab(NULL) + ylab(NULL) + theme_Publication() + 
  geom_hline(yintercept = bg_all$rho, colour = "black", show.legend = TRUE) + 
  geom_hline(yintercept = bg_nodisease$rho, colour = "grey", show.legend = TRUE) + 
  annotate("text", x = Inf, y = bg_all$rho, label = paste0("All proteins: ",round(bg_all$rho,3)), vjust = -0.5, hjust = 1.1, colour = "black",size=5) + 
  annotate("text", x = Inf, y = bg_nodisease$rho, label = paste0("Proteins without disease\nassociations: ",round(bg_nodisease$rho,3)), vjust = 1.1, hjust = 1.1, colour = "grey",size=5) + 
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),text = element_text(size=20,family = "sans"), axis.text = element_text(size=20,family = "sans")) + 
  coord_cartesian(ylim=c(-0.35,0.6))
# plot_grid(pMFO,pBPO,nrow=1)
# ggsave(file="disease_MFO_BPO.png",dpi=300,width=16,height=10)

# Option 3: boxplots in the same plot
both_n <- nrow(bh_MFO) + nrow(bh_BPO)
both_df <- data.frame(Disease_name = character(both_n), rho = double(both_n), padj = double(both_n))
both_df[1:nrow(bh_MFO),] = bh_MFO[,c(1,3,5)]
both_df[((nrow(bh_MFO)+ 1) : both_n),] = bh_BPO[,c(1,3,5)]
both_df$Aspect = c(rep("Molecular function",times=nrow(bh_MFO)),rep("Biological process",times=nrow(bh_BPO)))
both_df$nudge_x_val <- c(rep(0.15,times=nrow(bh_MFO)),rep(0.25,times=nrow(bh_BPO)))
ggplot(data = both_df, aes(x = fct_inorder(Aspect), y =rho)) + 
  geom_boxplot(outlier.shape = NA, width=0.2) + 
  geom_jitter(width = 0, shape = ifelse(both_df$padj < 0.05, 24, 21), size = ifelse(both_df$padj < 0.05, 5, 3),fill="black") +
  geom_label_repel(data = both_df %>% filter(rho > 0.4),
                  aes(label = Disease_name), size = 5, 
                  nudge_x = both_df%>%filter(rho > 0.4)%>%select(nudge_x_val)%>%unlist()) +
  xlab(NULL) + ylab("Spearman's coefficient") + theme_Publication() + 
  annotate("segment",x=0.7,xend=0.75,y= bg_all_MFO$rho,yend=bg_all_MFO$rho,linewidth=1)+
  annotate("segment",x=0.7,xend=0.75,y= bg_nodisease_MFO$rho,yend=bg_nodisease_MFO$rho,linewidth=1,color="grey")+
  annotate("text", x = 0.65, y = bg_all_MFO$rho, label = paste0("All proteins: ",round(bg_all_MFO$rho,3)), vjust = -0.5, colour = "black",size=5) + 
  annotate("text", x = 0.65, y = bg_nodisease_MFO$rho, label = paste0("Proteins without disease\nassociations: ",round(bg_nodisease_MFO$rho,3)), vjust = 1.1, colour = "grey",size=5) + 
  annotate("segment",x=2.25,xend=2.3,y= bg_all_BPO$rho,yend=bg_all_BPO$rho,linewidth=1)+
  annotate("segment",x=2.25,xend=2.3,y= bg_nodisease_BPO$rho,yend=bg_nodisease_BPO$rho,linewidth=1,color="grey")+
  annotate("text", x=2.35, y = bg_all_BPO$rho, label = paste0("All proteins: ",round(bg_all_BPO$rho,3)), vjust = -0.5, colour = "black",size=5) + 
  annotate("text", x =2.35, y = bg_nodisease_BPO$rho, label = paste0("Proteins without disease\nassociations: ",round(bg_nodisease_BPO$rho,3)), vjust = 1.1, colour = "grey",size=5) + 
  theme(axis.ticks.x = element_blank(),text = element_text(size=20,family = "sans"), axis.text = element_text(size=20,family = "sans")) + 
  coord_cartesian(ylim=c(-0.35,0.6))
# ggsave("scatter-box-diseases.png",dpi=300,width=13,height=11)

# Option 4: boxplots in the same plot, only show sig. diseases
both_n <- nrow(bh_MFO) + nrow(bh_BPO)
both_df <- data.frame(Disease_name = character(both_n), rho = double(both_n), padj = double(both_n), padj2 = double(both_n))
both_df[1:nrow(bh_MFO),] = bh_MFO[,c(1,3,5,6)]
both_df[((nrow(bh_MFO)+ 1) : both_n),] = bh_BPO[,c(1,3,5,6)]
both_df$Aspect = c(rep("Molecular function",times=nrow(bh_MFO)),rep("Biological process",times=nrow(bh_BPO)))
both_df$nudge_x_val <- c(rep(0.2,times=nrow(bh_MFO)),rep(0.25,times=nrow(bh_BPO)))
ggplot(data = both_df, aes(x = fct_inorder(Aspect), y =rho)) + 
  geom_boxplot(outlier.shape = NA, width=0.2) + 
  geom_jitter(width = 0.05, shape = ifelse(both_df$padj2 < 0.05, 24, 21), size = ifelse(both_df$padj2 < 0.05, 5, 2),fill="black") +
  geom_label_repel(data = both_df %>% filter(padj2 < 0.05),
                   aes(label = Disease_name), size = 5, 
                   nudge_x = both_df%>%filter(padj2 < 0.05)%>%select(nudge_x_val)%>%unlist()) +
  xlab(NULL) + ylab("Spearman's coefficient") + theme_Publication() + 
  annotate("segment",x=0.7,xend=0.75,y= bg_all_MFO$rho,yend=bg_all_MFO$rho,linewidth=1)+
  annotate("segment",x=0.7,xend=0.75,y= bg_nodisease_MFO$rho,yend=bg_nodisease_MFO$rho,linewidth=1,color="grey")+
  annotate("text", x = 0.65, y = bg_all_MFO$rho, label = paste0("All proteins: ",round(bg_all_MFO$rho,3)), vjust = -0.5, colour = "black",size=5) + 
  annotate("text", x = 0.65, y = bg_nodisease_MFO$rho, label = paste0("Proteins without disease\nassociations: ",round(bg_nodisease_MFO$rho,3)), vjust = 1.1, colour = "grey",size=5) + 
  annotate("segment",x=2.25,xend=2.3,y= bg_all_BPO$rho,yend=bg_all_BPO$rho,linewidth=1)+
  annotate("segment",x=2.25,xend=2.3,y= bg_nodisease_BPO$rho,yend=bg_nodisease_BPO$rho,linewidth=1,color="grey")+
  annotate("text", x=2.35, y = bg_all_BPO$rho, label = paste0("All proteins: ",round(bg_all_BPO$rho,3)), vjust = -0.5, colour = "black",size=5) + 
  annotate("text", x =2.35, y = bg_nodisease_BPO$rho, label = paste0("Proteins without disease\nassociations: ",round(bg_nodisease_BPO$rho,3)), vjust = 1.1, colour = "grey",size=5) + 
  theme(axis.ticks.x = element_blank(),text = element_text(size=20,family = "sans"), axis.text = element_text(size=20,family = "sans")) + 
  coord_cartesian(ylim=c(-0.35,0.6))
ggsave("4_figures/scatter-box-diseases-bonferroni.png",dpi=300,width=14,height=9)


sig_bhMFO <- bh_MFO%>%filter(p.adjust2 < 0.05)
genes_sig_bhMFO <- filtered_dis_fc_ic%>%
  mutate(Disease_name = paste0(Disease_name, "(*)"))%>%
  filter(Disease_name %in% sig_bhMFO$Disease_name)%>%
  select(Genes)%>%unlist()%>%unique()
print(length(genes_sig_bhMFO))