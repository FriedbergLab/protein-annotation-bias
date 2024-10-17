library(ineq)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpubr)
source('0_data/ggplot_theme_Publication-2.R')

# Helper functions ####
not_na <- function(x) {sum(!is.na(x))}

gini_vector <- function(x,column_name){
  data=x
  prot_go <- data%>%
    select(1,2,!!column_name)%>%
    pivot_wider(id_cols=Genes, names_from=Date, values_from = !!column_name)
  prot_go_dict = prot_go
  prot_go_dict[prot_go_dict == 0] <- NA
  gini_wt0 <- sapply(prot_go_dict[,2:ncol(prot_go_dict)], ineq, na.rm = T, simplify = T)
  n_wt0 <- sapply(prot_go_dict[,2:ncol(prot_go_dict)], not_na, simplify = T)
  val = gini_wt0
  n = n_wt0
  return(list(val,n))
}

# gini_vector including zero
gini_vector_0 <- function(x,column_name){
  data=x
  prot_go <- data%>%
    select(1,2,!!column_name)%>%
    pivot_wider(id_cols=Genes, names_from=Date, values_from = !!column_name)
  prot_go_dict = prot_go
  # prot_go_dict[prot_go_dict == 0] <- NA
  gini_wt0 <- sapply(prot_go_dict[,2:ncol(prot_go_dict)], ineq, na.rm = T, simplify = T)
  n_wt0 <- sapply(prot_go_dict[,2:ncol(prot_go_dict)], not_na, simplify = T)
  val = gini_wt0
  n = n_wt0
  return(list(val,n))
}


# Load data ####
# GO-based metrics
aspects <- c("F","P","C")
for (identifier in c("All")){
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


# (Fig 2A) Bar chart for Gini over time, all metrics ####
gini1 = list()
# year = "2013"
for (aspect in aspects) {
  aspect_data = list()
  for (metric in c("Count", "UniCount", "IA")) {
    df_name = paste0("prot_df_All_", aspect)
    df = get(df_name)
    g = gini_vector(df, metric)[[1]]%>%unlist()
    aspect_data[[paste0(metric, "_mean1")]] = mean(g)
    aspect_data[[paste0(metric, "_sd1")]] = sd(g)
  }
  gini1[[aspect]] = aspect_data
}
gini1_df = as.data.frame(do.call(rbind, gini1))
gini1_df$mean1 <- unlist(gini1_df$mean1)
gini1_df$sd1 <- unlist(gini1_df$sd1)
gini1_df$Aspect = names(gini1)
rownames(gini1_df) <- NULL
gini1_df <- gini1_df %>%
  pivot_longer(cols = -Aspect, 
               names_to = c("Metric","Type"), 
               names_pattern = "(.*)_(.*)",
               values_to = "Value")
gini1_df$Value <- unlist(gini1_df$Value)
gini1_df <- gini1_df %>%
  pivot_wider(names_from = "Type", values_from = "Value") %>%
  mutate(
    Aspect = case_when(
      Aspect == "F" ~ "Molecular function",
      Aspect == "P" ~ "Biological process",
      Aspect == "C" ~ "Cellular component",
      TRUE ~ Aspect
    ),
    Metric = case_when(
      Metric == "Count" ~ "Term count",
      Metric == "UniCount" ~ "Unique term count",
      Metric == "IA" ~ "Information content",
      TRUE ~ Metric
    )
  )%>%
  mutate(ymin_val = mean1-sd1, ymax_val = mean1+sd1)
gini1_df$Metric = factor(gini1_df$Metric, levels = c("Term count","Unique term count", "Information content"))
gini1_df$Aspect = factor(gini1_df$Aspect, levels = c("Molecular function","Biological process","Cellular component"))
custom_colors = c("Term count" = "darkgreen","Unique term count"= "#7fc97f","Information content" = "black")
pA = ggplot(gini1_df, aes(x = Aspect, y =mean1, fill = Metric)) + 
  geom_bar(stat="identity",position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin = ymin_val, ymax=ymax_val),position = position_dodge(0.9),width=0.25, color="#959595",linewidth = 1) + 
  theme_Publication() + scale_fill_manual(values=custom_colors)+
  labs(x = "Aspect", y ="Average Gini Coefficient\n2013-2022", fill = "Knowledge metric") + xlab("GO aspect") +  
  theme(legend.position = "right", legend.direction = "vertical", legend.text = element_text(size=13),text = element_text(size=15), axis.text = element_text(size=15))

# (Fig 2B) Gini overtime IC_MFO of all proteins vs proteins before 2013 ####
gini_MFO = gini_vector(prot_df_All_F, "IA")[[1]]
f13 <- prot_df_All_F%>%filter(Date == "2013" & IA > 0)
gini_MFO_13 = gini_vector(prot_df_All_F%>%filter(Genes %in% f13$Genes),"IA")[[1]]
gini2 = data.frame(gini_MFO,gini_MFO_13)
gini2$Date <- rownames(gini2)
gini2 <- gini2%>%
  rename_with(~c("All proteins","Proteins annotated\nbefore 2013","Date"))%>%
  pivot_longer(cols = 1:2, names_to = "Type",values_to = "Value")%>%
  mutate(Date = as.character(Date))
pB = ggplot(gini2, aes(x=Date, y=Value,group=Type, color = fct_inorder(Type))) + xlab("Year") + ylab("Gini coefficient") +
  geom_point(size=3) + geom_line(linewidth = 2) + theme_Publication() + scale_colour_Publication() + 
  scale_x_discrete(labels = c("2013", "", "2015", "", "2017", "", "2019", "", "2021", "")) + 
  coord_cartesian(ylim=c(0.45,0.65)) + scale_colour_manual(name = "", values = c("black","darkgrey")) + 
  theme(legend.position = "right", legend.direction = "vertical", legend.text = element_text(size=13), text = element_text(size=15), axis.text = element_text(size=15))+xlab("Year")

plot_grid(pA,pB,ncol=1,labels=c("A","B"),label_size = 20,hjust=-0.05,vjust=0.88,rel_heights = c(0.55,0.45))
ggsave("4_figures/gini-comb_EXP.png",dpi=300,width=11,height=8)


# (Fig 3) GINI LINE FPEs ####
gini_FC = gini_vector(prot_df_All_FC, "FractionalCount")[[1]]
gini_FC = data.frame(gini_FC)
gini_FC$Date = as.numeric(rownames(gini_FC))
colnames(gini_FC) <- c("Value","Date")
ggplot(gini_FC, aes(x=Date, y=Value)) + xlab("Year") + ylab("Gini Coefficient") +
  geom_point(size=3) + geom_line(linewidth = 2) + theme_Publication() + scale_colour_Publication() + 
  coord_cartesian(ylim=c(0.65,0.9)) + scale_x_continuous(breaks = seq(2013, 2021, by = 2)) +
  theme(text = element_text(size=15), axis.text = element_text(size=15))
ggsave("4_figures/gini-ovt-fpe.png",dpi=300,width=7,height=4)


# (Fig S4) GINI LINE ####
gini_all_metrics <- function(df_F, df_P, df_C, leg=FALSE){
  plots <- list()
  
  for (metric in c("Count", "UniCount", "IA")) {
    gini_F = gini_vector(df_F, metric)[[1]]
    gini_P = gini_vector(df_P, metric)[[1]]
    gini_C = gini_vector(df_C, metric)[[1]]
    gini = data.frame(gini_F, gini_P, gini_C)
    gini$Date = as.character(rownames(gini))
    colnames(gini) <- c("Molecular Function","Biological Process","Cellular Component","Date")
    gini <- gini%>%
        pivot_longer(cols = 1:3, names_to = "Aspect",values_to = "Value")
    p = ggplot(gini, aes(x=Date, y=Value, color = fct_inorder(Aspect), group = Aspect)) + 
        xlab(NULL) + ylab("Gini Coefficient") +
        geom_point() + geom_line() + theme_Publication() + scale_colour_Publication() + 
        coord_cartesian(ylim=c(0.2,0.72)) + scale_colour_manual(name = "GO Aspect", values = c("red","purple","#7fc97f")) + 
        scale_x_discrete(labels = c("2013", "", "2015", "", "2017", "", "2019", "", "2021", "")) + 
        theme(text = element_text(size=15), axis.text.x=element_blank())
    if (metric == "IA"){
      p1 <- p + theme(legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size=15))
      legend_p1 <- ggpubr::get_legend(p1)
      p <- p +  theme(legend.position = "none",axis.text.x=element_text(size=13))
    } else {
      p <- p + theme(legend.position = "none")
    }
    plots[[metric]] <- p
  }
  # Combine plots using plot_grid
  combined_plot <- plot_grid(plotlist = plots, ncol = 1, align = "v",labels=c("Term Count", "Unique Term Count", "Information Content"),hjust=-0.3,label_size=12)
  if (leg){
    return(list(combined_plot, legend_p1))
  } else {
    return(combined_plot)
  }
}

p_all = gini_all_metrics(prot_df_All_F, prot_df_All_P, prot_df_All_C)

# Distinguish between before 2013 and after 2013
filter_year <- function(prot_df, year) {
  # Filter genes for the specified year and aspect
  df_year <- prot_df %>%
    filter(Date == as.character(year) & IA > 0) %>%
    select(Genes) %>%
    unlist()
  # Filter genes not in the specified year but with IA > 0
  df_not_year <- prot_df %>%
    filter(IA > 0 & !(Genes %in% df_year)) %>%
    select(Genes) %>%
    unlist()

  # Create dataframes for the specified year and non-specified year
  prot_df_year <- prot_df %>% filter(Genes %in% df_year)
  prot_df_not_year <- prot_df %>% filter(Genes %in% df_not_year)
  # Return the filtered dataframes as a list
  return(list(prot_df_year, prot_df_not_year))
}
F_13_not <- filter_year(prot_df_All_F, 2013)
P_13_not <- filter_year(prot_df_All_P, 2013)
C_13_not <- filter_year(prot_df_All_C, 2013)

p_13_leg <- gini_all_metrics(F_13_not[[1]], P_13_not[[1]], C_13_not[[1]], leg=TRUE)
p_not_13 <- gini_all_metrics(F_13_not[[2]], P_13_not[[2]], C_13_not[[2]])

# Combine plots using plot_grid
combined_plot <- plot_grid(p_all, p_13_leg[[1]], p_not_13, nrow=1, align = "h",labels = c("A","B","C"), hjust = -0.05)
legend_plot <- plot_grid(combined_plot, p_13_leg[[2]], nrow=2, rel_heights = c(1, 0.1))
ggsave("4_figures/gini-ovt-ic-13-not13_EXP.png",dpi=300,width=17,height=12)

# (Fig S5) More proof on Gini decrease-increase contrast ####
f13 <- prot_df_All_F%>%filter(Date == "2013" & IA > 0)
f22 <- prot_df_All_F%>%filter(Date == "2022" & IA > 0)
f13in22 <- prot_df_All_F%>%filter(Date == "2022" & IA > 0 & Genes %in% f13$Genes)
# f13in22_2 = prot_df_All_F%>%filter(Genes %in% f13$Genes)%>%filter(Date == "2022")
f13not <- prot_df_All_F%>%filter(!(Genes %in% f13$Genes) & IA > 0)
f13vsnot13in22 = rbind(f13in22[,c("Genes","IA")],f13not[f13not$Date == "2022",c("Genes","IA")])
f13vsnot13in22$`Year of first experimental annotation` = c(rep("In or before 2013",length=nrow(f13in22)), rep("After 2013",length=nrow(f13vsnot13in22)-nrow(f13in22)))
ggplot(f13vsnot13in22,aes(x=IA,fill=`Year of first experimental annotation`)) + geom_density(alpha=0.5) +
  ylab("Density") + xlab("IC_MFO of proteins in 2022") +
  theme_Publication() + theme(legend.position = c(0.6,0.8),legend.direction = "vertical")
ggsave("4_figures/gini_evidence.png",dpi=300,width=7,height=4)

# # EXTRA: Gini coefficient including or excluding zero ####
# gini_no0 = list()
# for (aspect in c("F","P","C")){
#   for (metric in c("Count","UniCount","IA")){
#     df_name=paste0("prot_df_All_",aspect)
#     df = get(df_name)
#     gini_no0[[aspect]][[metric]] = gini_vector(df, metric)[[1]]
#   }
# }
# flatten_list <- function(nested_list) {
#   flat_list <- list()
#   for (aspect in names(nested_list)) {
#     for (metric in names(nested_list[[aspect]])) {
#       values <- nested_list[[aspect]][[metric]]
#       for (i in seq_along(values)) {
#         flat_list[[length(flat_list) + 1]] <- data.table(
#           aspect = aspect,
#           metric = metric,
#           value = values[i]
#         )
#       }
#     }
#   }
#   return(flat_list)
# }
# flattened_list <- flatten_list(gini_no0)
# gini_no0_df <- rbindlist(flattened_list, fill = TRUE)

# gini_0 = list()
# for (aspect in c("F","P","C")){
#   for (metric in c("Count","UniCount","IA")){
#     df_name=paste0("prot_df_All_",aspect)
#     df = get(df_name)
#     gini_0[[aspect]][[metric]] = gini_vector_0(df, metric)[[1]]
#   }
# }
# flattened_list <- flatten_list(gini_0)
# gini_0_df <- rbindlist(flattened_list, fill = TRUE)
