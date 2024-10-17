library(tidyverse)
library(cowplot)
library(ggplot2)
source('0_data/ggplot_theme_Publication-2.R')

# Publication-to-annotation timeline ####
pmid_year <- read.csv(file="3_delay_curation/annotation_pub_years.csv",sep="\t",header=T)%>%
  rename_with(~c("Genes","Common","GO_ID","aspect","reference","evidence","year","PMID","pub_year"))
pmid_year_F <- pmid_year%>%filter(aspect == "F")
pmid_year_P <- pmid_year%>%filter(aspect == "P")
pmid_year_C <- pmid_year%>%filter(aspect == "C")

# (Fig 6) Delay in curation ####
# # Distinguish same protein-GO, different reference
# annt_pub2 <- pmid_year%>%
#   filter(pub_year > 0 &!is.na(pub_year))%>%
#   filter(GO_ID != "GO:0005515" & GO_ID != "GO:0042802")%>% # IntAct annotations
#   select(c(Genes, GO_ID, reference, evidence, year, pub_year))%>%
#   arrange(year)%>%
#   distinct(Genes,GO_ID, reference, .keep_all = T)%>%
#   mutate(pub_year = ifelse(year < pub_year, year, pub_year))%>%
#   mutate(delay = as.integer(year) - as.integer(pub_year))

# Not distinguishing same protein-GO, different reference
annt_pub3 <- pmid_year%>%
  filter(pub_year > 0 &!is.na(pub_year))%>%
  filter(GO_ID != "GO:0005515" & GO_ID != "GO:0042802")%>% # IntAct annotations
  select(c(Genes, GO_ID, reference, evidence, year, pub_year))%>%
  arrange(year)%>%
  distinct(Genes,GO_ID,.keep_all = T)%>%
  mutate(pub_year = ifelse(year < pub_year, year, pub_year))%>%
  mutate(delay = as.integer(year) - as.integer(pub_year))

year_pubyear = annt_pub3%>%
  select(c(Genes,GO_ID,year, pub_year))%>%
  pivot_longer(cols=3:4, values_to = "y",names_to = "time")
pB <- ggplot(data=year_pubyear, aes(x=y,fill=time)) + geom_bar(alpha=0.5,position="dodge",width = 0.8) + 
  xlab("Year") + ylab("Number of annotations") + theme_Publication() +
  scale_fill_manual(name="",values = c("salmon","cyan"),labels = c("pub_year"="Publication year of PubMed reference","year"= "Annotation year in GO")) + 
  coord_cartesian(xlim=c(1995,NA), ylim=c(0,15000)) + scale_x_continuous(breaks = c(1995,2001,2007,2012, 2017,2022)) + 
  theme(legend.position = "inside", legend.position.inside = c(0.38,0.9), legend.direction = "vertical", legend.text = element_text(size=18,family="sans"), text = element_text(size=18,family = "sans"), axis.text = element_text(size=18,family = "sans")) 
pA <- ggplot(data=annt_pub3, aes(x=delay)) + geom_histogram(binwidth=1) + xlab("Years between publication and annotation") + ylab("Number of annotations") + theme_Publication() + geom_vline(xintercept = median(annt_pub3$delay)) + scale_x_continuous(breaks = c(0,median(annt_pub3$delay),10,20,30,40,50)) + 
  coord_cartesian(ylim=c(0,15000)) + 
  theme(text = element_text(size=18,family = "sans"), axis.text = element_text(size=18,family = "sans")) 
plot_grid(pA,pB,nrow=1,align="h",labels = c("A","B"),label_size = 20,rel_widths = c(0.45,0.55))
ggsave("4_figures/delay_curation_comb.png",width=16,height=8,dpi=300)

