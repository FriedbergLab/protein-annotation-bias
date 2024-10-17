library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggnewscale)
library(ggrepel)
library(cowplot)
options(scipen=999)
source('0_data/ggplot_theme_Publication-2.R')

# Helper functions ####
wider_format <- function(df, column){
  df_byyear <- df%>%
    select(1,2, !!column)%>%
    pivot_wider(id_cols = Genes, names_from = Date, values_from = !!column)
  return(df_byyear)
}

cor_summary <- function(x,y){
  test <- cor.test(x,y, method="spearman",use="pairwise.complete.obs",quiet=TRUE)
  valid_pairs <- complete.cases(x,y)
  tibble(
    rho = test$estimate,
    pval = test$p.value,
    num = sum(valid_pairs)
  )
}

plot_heatmap <- function(df,method){
  plotdata = df
  cc = cor(plotdata, method=method,use="pairwise.complete.obs")
  cc_df = as.data.frame(cc)
  cc_df$metrics = row.names(cc_df)
  ccm = melt(cc_df, id = "metrics")
  ccm$metrics <- factor(ccm$metrics,levels=rev(unique(fct_inorder(ccm$metrics))))
  ccm$variable <- factor(ccm$variable, levels = unique(fct_inorder(ccm$variable)))
  
  p1 = ggplot(ccm, aes(x = (variable), y = (metrics)))
  p2 <- p1 + geom_tile(aes(fill = value), colour = "lightgrey") + 
    geom_text(aes(label = round(value, 2)),colour = ifelse(ccm$value > 0.8, "white","black")) + 
    coord_equal() + labs(colour = "Spearman's correlation coefficient") + 
    scale_fill_gradient2(low = "red", mid="white", high = "blue", midpoint = 0, limits=c(-1,1)) + 
    theme(axis.text.y = element_text(size = 12, face = "bold", colour = "black"), legend.title = element_text(size = 10, face = "bold"),
          legend.position = "bottom", axis.text.x = element_text(size = 12, angle = 90, face = "bold",colour = "black", vjust = 0.5, hjust = 0), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = NA), axis.ticks = element_blank()) + 
          labs(x= "", y = "", fill = "Spearman's correlation") +
    scale_x_discrete(position = "top") + scale_y_discrete(limits = (levels(ccm$metrics)))
  return(p2)  
}

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
# Interest metric
fc = readRDS(file="1_metrics/interest_FPE.RDS")
fc_clinical = readRDS(file="1_metrics/interest_FPE_clinical.RDS")

# ID mapping
proteinstring <- read.csv(file="0_data/ext_data/sw-human-2024.txt",sep="\t",header=F)%>%
  select(2,3,5)%>%
  rename_with(~c("Genes","Common","STRING"))


# Format data ####
prot_df_All_FC <- fc[,c(2,4:13)]%>%
  reshape2::melt()%>%
  rename_with(~c("Genes","Date","FractionalCount"))%>%
  mutate_at(2,as.character)%>%
  distinct(Genes,Date,FractionalCount,.keep_all = T)
prot_df_All_clinical <- fc_clinical[,c(2,4:13)]%>%
  reshape2::melt()%>%
  rename_with(~c("Genes","Date","FractionalCount"))%>%
  mutate_at(2,as.character)%>%
  distinct(Genes,Date,FractionalCount,.keep_all = T)

icF <- wider_format(prot_df_All_F, "IA")
icP <- wider_format(prot_df_All_P, "IA")
icC <- wider_format(prot_df_All_C, "IA")
# Calculate knowledge gain
icF = icF%>%mutate(gain = `2022` - `2013`)%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0) | gain < 0, 0, gain))%>%
  filter(`2022` > 0) # remove proteins that are currently fake
icP = icP%>%mutate(gain = `2022` - `2013`)%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0) | gain < 0, 0, gain))%>%
  filter(`2022` > 0) # remove proteins that are currently fake
icC = icC%>%mutate(gain = `2022` - `2013`)%>%
  mutate(gain = ifelse((is.na(gain) & `2013` == 0 & `2022` == 0) | gain < 0, 0, gain))%>%
  filter(`2022` > 0) # remove proteins that are currently fake

# Combine knowledge and interest 
fc_ic <- full_join(fc[,c(1:4,13:15)],icF[,c(1,2,11,12)],by="Genes")%>%
  full_join(.,icP[,c(1,2,11,12)],by="Genes")%>%
  full_join(.,icC[,c(1,2,11,12)],by="Genes")%>%
  rename_with(~c("STRING","Genes","Common","FPE_13","FPE_22","relgain_FPE","gain_FPE","MFO_13","MFO_22","gain_MFO", "BPO_13","BPO_22", "gain_BPO","CCO_13","CCO_22", "gain_CCO"))


# (Fig S7) Heatmap correlation ####
fc_plot <- fc_ic%>%
  rename_with(~c("STRING","Genes","Common","Number of FPEs in 2013","Number of FPEs in 2022",
                 "Relative gain in FPEs 2013-2022", "Gain in FPEs 2013-2022", 
                 "IC in MFO in 2013", "IC in MFO in 2022", "Gain IC in MFO 2013-2022",
                 "IC in BPO in 2013", "IC in BPO in 2022", "Gain IC in BPO 2013-2022",
                 "IC in CCO in 2013", "IC in CCO in 2022", "Gain IC in CCO 2013-2022"))
plot_heatmap(fc_plot[,c(4,7,8,10,11,13,14,16)], method="spearman") # high cor between FPE-based metrics, low cor among the rest
ggsave("4_figures/heatmap_correlation_nomarker.png",dpi=300,width=8,height=8)

# More correlation
cor_summary(fc_ic$FPE_22,fc_ic$MFO_22)
cor_summary(fc_ic$FPE_22,fc_ic$BPO_22)
cor_summary(fc_ic$FPE_22,fc_ic$CCO_22)
fc_ic_all <- full_join(fc[,c(1:13)],icF[,c(1:11)],by="Genes")%>%
  full_join(.,icP[,c(1:11)],by="Genes")%>%
  full_join(.,icC[,c(1:11)],by="Genes")%>%
  rename_with(~c("STRING","Genes","Common","FPE_13","FPE_14","FPE_15","FPE_16","FPE_17","FPE_18","FPE_19","FPE_20","FPE_21","FPE_22",
                 "MFO_13","MFO_14","MFO_15","MFO_16","MFO_17","MFO_18","MFO_19","MFO_20","MFO_21","MFO_22",
                 "BPO_13","BPO_14","BPO_15","BPO_16","BPO_17","BPO_18","BPO_19","BPO_20","BPO_21","BPO_22",
                 "CCO_13","CCO_14","CCO_15","CCO_16","CCO_17","CCO_18","CCO_19","CCO_20","CCO_21","CCO_22"))
cor_fc_ic <- function(year, aspect, fc_ic_df){
  cor_summary(fc_ic_df[[paste0("FPE_",year)]],fc_ic_df[[paste0(aspect,"_",year)]])
}
loop_cor_fc_ic <- function(yy_range, aspect, fc_ic_df){
  cor_df <- lapply(yy_range, function(x) cor_fc_ic(x, aspect, fc_ic_df))
  cor_df <- data.frame(do.call(rbind, cor_df))%>%
    mutate(year = 2013:2022)%>%
    select(year, rho, pval, num)%>%
    rename_with(~c("Year",paste0("rho_",aspect),paste0("pval_",aspect),paste0("num_",aspect)))
  return(cor_df)
}
cor_F <- loop_cor_fc_ic(13:22, "MFO", fc_ic_all)
cor_P <- loop_cor_fc_ic(13:22, "BPO", fc_ic_all)
cor_C <- loop_cor_fc_ic(13:22, "CCO", fc_ic_all)

cor_FPC <- left_join(left_join(cor_F[,c(1,2,4)],cor_P[,c(1,2,4)],by="Year"),cor_C[,c(1,2,4)],by="Year")
print(cor_FPC)

# (Fig 3) Scatterplot IC vs FPE ####
scatterplot_ic_fc <- function(prot_df, clinical_df) {
  prot_clinical <- left_join(prot_df, clinical_df, by=c("Genes","Date"))%>%
    filter(Count > 0)
  top4_FC <- prot_clinical%>%filter(Date == "2022") %>% top_n(n=4, FractionalCount) %>% select(Genes) %>% unlist()
  top7_IC <- prot_clinical%>%filter(Date == "2022") %>% top_n(n=7, IA) %>% select(Genes) %>% unlist()
  data_ends <- prot_clinical %>%
    filter(Genes %in% top4_FC | Genes %in% top7_IC) %>%
    arrange(desc(ifelse(Date == "2022", IA, 0)), Date)%>%
    left_join(proteinstring, by="Genes")

  # Threads of top genes plotted separately with other data points (same scale)
  pt <- ggplot(mapping = aes(x = FractionalCount, y = IA))+ theme_Publication() + scale_colour_Publication() +
    geom_point(data = prot_clinical%>%filter(Date == "2022"), size = 0.7)+labs(colour="")+
    new_scale_colour() +
    geom_point(data = data_ends, aes(colour = Genes), size = 1.3, show.legend = F) + 
    geom_line(data = data_ends, aes(colour = Genes), linewidth = 1,show.legend = F) +  
    geom_label_repel(aes(label = Common, color = Genes), max.overlaps = 20, data = data_ends%>%filter(Date == "2022"), show.legend = F, size = 4.5) + 
    scale_colour_manual(values = c("cyan","salmon","green","brown","#8931EF","#877731","#FF00BD","#001A46","#66a326","#E11845","#0000cd")) +
    xlab("Number of FPEs") + ylab("Information Content") +
    scale_y_continuous(limits = c(NA,NA)) + scale_x_log10(limits = c(NA, NA)) +
    theme(legend.position = "none", text = element_text(size=20), axis.text = element_text(size=20))
  return(pt)
}

ptF <- scatterplot_ic_fc(prot_df_All_F, prot_df_All_clinical)
ptP <- scatterplot_ic_fc(prot_df_All_P, prot_df_All_clinical)
ptC <- scatterplot_ic_fc(prot_df_All_C, prot_df_All_clinical)
plot_grid(ptF,ptP,ptC,align = 'h', labels = c("MFO", "BPO","CCO"), label_size = 20, hjust = -1, nrow = 1)
ggsave("4_figures/scatter-ic-fc_EXP.png",dpi=300,width=19,height=7)
