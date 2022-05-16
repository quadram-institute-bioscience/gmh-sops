#!/usr/bin/env Rscript
#for: humann.nf

#loading libraries
library('ggplot2')
library('tidyr')

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
  stop("Missing inputs and parameters", call.=FALSE)
} else if (length(args) == 2) {
  genus_species = args[1] # Genus relative abundance
  outdir = args[2] # Output dir
}

#Creating directory to store plots
#If directory doesn't exist create one.
if (!file.exists(outdir)){
  dir.create(file.path(outdir))
}

#checking the end of the filename
if (grepl("Species-rel.csv", genus_species)=="TRUE"){
  datatype="Species" 
} else if (grepl("Genus-rel.csv", genus_species)=="TRUE") {
  datatype="Genus"
} else {
  datatype="Genus/Species"
}

#Filtering top N genus/species
df=read.csv(genus_species, header=TRUE) #Reading file.
colnames(df)[1]<-"Sample" 
df=subset(df, select=-c(Others))
N=ncol(df)-1

df_long=df %>% 
  gather(key="id", value="Sample_id",Sample) %>%
  gather(key = "species", value = "Abundance",colnames(df)[2]:colnames(df)[ncol(df)],factor_key = TRUE)
  
heatmap=ggplot(df_long, aes(x=Sample_id, y=species, fill= Abundance)) + 
  ggtitle(paste0("Top ",N," ",datatype, c(" in the experiment"))) +
  geom_tile() +
  #scale_fill_gradient(low="green", high="black") +
  scale_fill_distiller(palette = "Blues") +
  labs(x="", y="",size="",color="",fill="Abundance (%)") +
  theme(
    axis.text.x = element_text(colour = "black", size = 8, angle=90, hjust = 1, vjust = 0.5) ,
    axis.text.y = element_text(colour = "black", size = 8)
    )

ggsave(file.path(outdir,paste0(basename(datatype),c("_relabundance.svg"))),plot=heatmap,width=0.8*N,height=10,device="svg")
