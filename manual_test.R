# manual test
library(shiny)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(plotly)
library(bslib)
library(tibble)
library(ape)
library(vegan)
library(usedist)
library(scales)


# Load helper functions
source("utils.R")

# Color palette for discrete groups
col21 <- rev(c("tomato1","darkblue","turquoise1","lightblue","darkred","mediumblue","purple","bisque",
               "greenyellow","yellow","violetred2","darkgreen","darkgoldenrod1","deeppink3","cadetblue4",
               "orchid2","seagreen3","purple4","dodgerblue2","red","gray27"))



md = read.csv('test_data/test1/metadata.csv')
#dm = read.table('test_data/nhlbi-167_diversity_shiny_test/weighted-unifrac.tsv',check.names = FALSE, sep = "\t")

#df <- read.table('test_data/nhlbi-167_diversity_shiny_test/weighted-unifrac.tsv', header = TRUE, sep = "\t", check.names = FALSE)
#plot_beta_pcoa(dm,md,'Stage','x')


df<- read.csv('test_data/test1/NCPHControl_merged.NIDDK26.metaphlan_bugs_list.amended.tsv', header = TRUE, check.names = FALSE)
df.species = tax_collapse(df,'genus',top_n = 9)[[3]]


plot_tax_stacked_bar(exps = df.species,pheno = md,split_by_2 = 'group',split_by_1 = 'type', show_names = F)
