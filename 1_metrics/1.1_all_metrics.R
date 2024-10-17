library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggrepel)
source("0_data/ggplot_theme_Publication-2.R")

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
prot_df_All_FC <- fc[,c(2,4:13)]%>%
  reshape2::melt()%>%
  rename_with(~c("Genes","Date","FractionalCount"))%>%
  mutate_at(2,as.character)%>%
  distinct(Genes,Date,FractionalCount,.keep_all = T)

# ID mapping
proteinstring <- read.csv(file="0_data/ext_data/sw-human-2024.txt",sep="\t",header=F)%>%
  select(2,3,5)%>%
  rename_with(~c("Genes","Common","STRING"))


# Combine HTP and LTP annotations ####
combine_TP_dataframes <- function(LTP_df, HTP_df) {
  combined_df <- rbind(LTP_df, HTP_df) %>%
    mutate(TP = rep(c("LTP", "HTP"), times = c(nrow(LTP_df), nrow(HTP_df))))
  return(combined_df)
}
for (aspect in aspects){
  assign(paste0("prot_df_HLTP_",aspect), combine_TP_dataframes(get(paste0("prot_df_LTP_",aspect)), get(paste0("prot_df_HTP_",aspect))))
}


# (Fig S1) With and without protein-binding-only proteins ####
df1 <- prot_df_All_F%>%
  arrange(Date)%>%
  filter(complete.cases(.) & Count > 0)
df <- prot_df_All_F%>%
  arrange(Date)%>%
  filter(complete.cases(.) & Count > 0 & IA > 0.16)
combined <- rbind(df,df1)
combined$type <- c(rep("Without protein-binding-only proteins", times = nrow(df)),rep("All proteins",times=nrow(df1)))
mean_ic = aggregate(IA ~ Date + type, combined, mean)
ggplot(combined, aes(x=factor(Date),y=IA,color=type)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(data=mean_ic, aes(x=fct_inorder(as.character(Date)),y=IA),size = 3,position=position_dodge(0.75)) + 
  xlab("Year") + ylab("Information Content") + 
  theme_Publication() + scale_colour_Publication() + 
  scale_y_continuous() + coord_cartesian(ylim=c(0,50)) + scale_x_discrete(labels = c("2013", "", "2015", "", "2017", "", "2019", "", "2021", "")) + 
  scale_color_manual(name="",values=c("black","grey")) + 
  theme(legend.position = "none", text = element_text(size=20), axis.text = element_text(size=20))
ggsave("4_figures/ovt-boxplot-protein-binding.png",dpi=300,width=12,height=8)

# (Fig 1A) Boxplot IC proteins by year by thruput ####
boxplot_byTP_IC <- function(prot_df_HLTP, ontology, y_limit) {
  df <- prot_df_HLTP %>%
    arrange(Date) %>%
    group_by(Genes, Date) %>%
    mutate(All = sum(IA)) %>%
    filter(All > 0)
  
  if (ontology == "CCO") {
    df <- df %>%
      mutate(TP = replace(TP, TP == "HTP", "High-throughput experiment, > 100 proteins per article")) %>%
      mutate(TP = replace(TP, TP == "LTP", "Low-throughput experiment, \u2264 100 proteins per article"))
  }
  
  mean_ic = aggregate(IA ~ Date + TP, df, mean)
  
  p <- ggplot(df, aes(x = factor(as.character(Date)), y = IA, color = TP)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylab("Information Content") + xlab("") +
    scale_colour_Publication() + theme_Publication() + scale_y_continuous() + 
    geom_point(data = mean_ic, aes(x = fct_inorder(as.character(Date)), y = IA), size = 3, position = position_dodge(0.75)) +
    coord_cartesian(ylim = c(0, y_limit)) + 
    scale_x_discrete(labels = c("2013", "", "2015", "", "2017", "", "2019", "", "2021", ""))
  
  if (ontology == "CCO") {
    p <- p + labs(colour = "") + 
      theme(legend.text = element_text(size = 20), legend.direction = "vertical")
  }
  
  return(p)
}

pF <- boxplot_byTP_IC(prot_df_HLTP_F, "MFO", 40)
pP <- boxplot_byTP_IC(prot_df_HLTP_P, "BPO", 80)
pC <- boxplot_byTP_IC(prot_df_HLTP_C, "CCO", 40)
# Combine:
pcol <- plot_grid(pF+theme(legend.position = "none", text = element_text(size=20), 
                           axis.text = element_text(size=20)),NULL,
                  pP+theme(legend.position = "none", text = element_text(size=20), 
                           axis.text = element_text(size=20)),NULL,
                  pC+theme(legend.position = "bottom", text = element_text(size=20), 
                           axis.text = element_text(size=20), legend.margin=margin(0,0,0,0), legend.box.margin = margin(-40,0,0,0)),
                  labels = c("MFO","","BPO","","CCO"), 
                  hjust = -2.5, rel_heights = c(1,-0.06,1,-0.06,1.10),ncol=1)

# (Fig 1B) Boxplot IC proteins by year no throughput, rich genes highlighted ####
boxplot_rich_IC <- function(prot_df, ontology, sc, y_limit) {
  df <- prot_df %>%
    arrange(Date) %>%
    filter(IA > 0) %>%
    filter(complete.cases(.))

  mean_ic = aggregate(IA ~ Date, df, mean)
  
  top7_IA <- df %>% filter(Date == 2022) %>% top_n(n=7, IA) %>% select(Genes) %>% unlist()
  data_ends <- df %>%
    filter(Genes %in% top7_IA) %>%
    mutate(IA = IA/sc) %>%
    arrange(desc(IA))%>%
    left_join(proteinstring, by="Genes")

  p <- ggplot(df, aes(x = factor(as.character(Date)), y = IA)) + 
    geom_boxplot(outlier.shape = NA) + xlab("") + ylab("Information Content") +
    scale_colour_Publication() + theme_Publication() + scale_y_continuous() +
    geom_point(data = mean_ic, aes(x = fct_inorder(as.character(Date)), y = IA), size = 3, position = position_dodge(0.75)) +
    geom_line(data = data_ends, aes(color = Genes, group = Genes, x = fct_inorder(as.character(Date)), y = IA), linewidth = 1, alpha = 0.5, show.legend = F) + 
    scale_colour_manual(values = c("#8931EF","#877731","#FF00BD","#001A46","#66a326","#E11845","#0000cd")) + 
    geom_label_repel(aes(label = Common, color = Genes), max.overlaps = 20, data = data_ends %>% filter(Date == 2022), show.legend = F, size = 4.5) + 
    scale_y_continuous(name = "Information Content",
                       sec.axis = sec_axis(~.*sc, name = "Information Content (high)")) +
    scale_x_discrete(labels = c("2013", "", "2015", "", "2017", "", "2019", "", "2021", "")) + 
    coord_cartesian(ylim = c(0, y_limit))

  return(p)
}

pFa <- boxplot_rich_IC(prot_df_All_F, "F", sc = 5, y_limit=40)
pPa <- boxplot_rich_IC(prot_df_All_P, "P", sc = 7.5, y_limit=80)
pCa <- boxplot_rich_IC(prot_df_All_C, "C", sc = 4, y_limit=40)

# Combine:
pcol2 <- plot_grid(pFa+theme(legend.position = "none", text = element_text(size=20), 
                           axis.text = element_text(size=20)),NULL,
                  pPa+theme(legend.position = "none", text = element_text(size=20), 
                           axis.text = element_text(size=20)),NULL,
                  pCa+theme(legend.position = "bottom", text = element_text(size=20), 
                           axis.text = element_text(size=20)),NULL,
                  labels = c("MFO","","BPO","","CCO",""),
                  hjust = -2.5, rel_heights = c(1,-0.06,1,-0.06,1,0.105),ncol=1)
plot_grid(pcol,pcol2,nrow=1,labels="AUTO",hjust= -0.25, label_size = 20)
ggsave("4_figures/overtime-box-ic-alltp_EXP_rich-combined-new.png",width=15,height=20,dpi=300)


# (Fig S6) Boxplot of number of FPEs ####
df <- prot_df_All_FC%>%
  arrange(Date)%>%filter(FractionalCount > 0)
mean_fc = aggregate(FractionalCount ~ Date, df, mean)
print(mean_fc)
median_fc = aggregate(FractionalCount ~ Date, df, median)
print(median_fc)
sc = 150
top5_FC <- df %>% filter(Date == 2022) %>% top_n(n=5, FractionalCount) %>% select(Genes) %>% unlist()
data_ends <- df %>%
  filter(Genes %in% top5_FC) %>%
  mutate(FractionalCount = FractionalCount/sc) %>%
  left_join(proteinstring, by="Genes")
pFPE <- ggplot(df, aes(x = factor(as.character(Date)), y = FractionalCount)) + 
  geom_boxplot(outlier.shape = NA) + xlab("Year") + ylab("Number of FPEs") +
  scale_colour_Publication()+ theme_Publication() + scale_y_continuous()+
  geom_point(data=mean_fc, aes(x=fct_inorder(as.character(Date)),y=FractionalCount),size = 3,position=position_dodge(0.75)) +
  geom_line(data=data_ends, aes(color=Genes,group=Genes,x=fct_inorder(as.character(Date)),y=FractionalCount),linewidth = 1,alpha = 0.5, show.legend = FALSE) + 
  scale_colour_manual(values = c("#8931EF","#877731","#001A46","#66a326","#E11845")) + 
  geom_label_repel(aes(label = Common, color = Genes), max.overlaps = 20, data = data_ends%>% filter(Date == 2022), show.legend = F, size = 4.5) + 
  scale_y_continuous(name = "Number of FPEs",
                     sec.axis = sec_axis(~.*sc, name = "Number of FPEs (high)")) +
  scale_x_discrete(labels = c("2013", "", "2015", "", "2017", "", "2019", "", "2021", "")) + 
  coord_cartesian(ylim=c(0,400))
ggsave("4_figures/overtime-box-fc-alltp_EXP_rich.png",width=11,height=8,dpi=300)