---
title: "Wheat Microbiome paper - R script for figures"
author: "Marie Simonin"
output: html_document
---


# Code implemented in R to generate the figures of the Simonin et al. paper on the wheat microbiome

##Figure 1
###Load the diversity and structure data for the 8 wheat genotypes in one soil (FR2) - For detailed script on how to generate the diversity indices and NMDS, see Figure 3
```{r}
Div <- read.table("Diversity&Structure_data_8genotypes_FR2soil.txt", header=TRUE, check.names = FALSE)
```

###Graph ESV Richness - 8 genotypes (Figure 1A)
```{r warning=FALSE, message=FALSE}
library(Rmisc)
library(ggplot2)
color= c("#4CB5F5",  "#D70026", "#EDB83D",  "navy",  "#F77604", "#B3C100", "#EC96A4", "#5BC8AC")
Total <- summarySE(Div, measurevar="S", groupvars=c("Genotype"), na.rm = TRUE)
Bac <- summarySE(Div, measurevar="S_16S", groupvars=c("Genotype"), na.rm = TRUE)
Euk <- summarySE(Div, measurevar="S_18S", groupvars=c("Genotype"), na.rm = TRUE)
p2=ggplot() + geom_point(data=Total, aes(x=reorder(Genotype, S), y=S, color=Genotype), size=5, shape=15) + geom_point(data=Bac, aes(x=Genotype, y=S_16S, color=Genotype), size=4, shape=5) + geom_point(data=Euk, aes(x=Genotype, y=S_18S, color=Genotype), size=4, shape=1) + geom_errorbar(data=Total, aes(x=Genotype, ymin=S-se, ymax=S+se, color=Genotype), width=.2)+ geom_errorbar(data=Bac, aes(x=Genotype, ymin=S_16S-se, ymax=S_16S+se, color=Genotype), width=.2) + geom_errorbar(data=Euk, aes(x=Genotype, ymin=S_18S-se, ymax=S_18S+se, color=Genotype), width=.2)+ theme_classic()+xlab("Genotypes")+ylab("Observed ESVs Richness")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```

###Graph Phylogenetic diversity - 8 genotypes (Figure 1B)
```{r warning=FALSE}
color= c("#4CB5F5",  "#D70026", "#EDB83D",  "navy",  "#F77604", "#B3C100", "#EC96A4", "#5BC8AC")
Total <- summarySE(Div, measurevar="PD", groupvars=c("Genotype"), na.rm = TRUE)
Bac <- summarySE(Div, measurevar="PD_16S", groupvars=c("Genotype"), na.rm = TRUE)
Euk <- summarySE(Div, measurevar="PD_18S", groupvars=c("Genotype"), na.rm = TRUE)
p2=ggplot() + geom_point(data=Total, aes(x=reorder(Genotype, PD), y=PD, color=Genotype), size=5, shape=15) + geom_point(data=Bac, aes(x=Genotype, y=PD_16S, color=Genotype), size=4, shape=5) + geom_point(data=Euk, aes(x=Genotype, y=PD_18S, color=Genotype), size=4, shape=1) + geom_errorbar(data=Total, aes(x=Genotype, ymin=PD-se, ymax=PD+se, color=Genotype), width=.2)+ geom_errorbar(data=Bac, aes(x=Genotype, ymin=PD_16S-se, ymax=PD_16S+se, color=Genotype), width=.2) + geom_errorbar(data=Euk, aes(x=Genotype, ymin=PD_18S-se, ymax=PD_18S+se, color=Genotype), width=.2)+ theme_classic()+xlab("Genotypes")+ylab("Faith's Phylogenetic Diversity")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```

###Graph correlation between Prokaryotic and Eukaryotic richness - 8 genotypes (Figure 1C)
```{r warning=FALSE}
color= c("#4CB5F5",  "#D70026", "#EDB83D",  "navy",  "#F77604", "#B3C100", "#EC96A4", "#5BC8AC")
p2=ggplot(data=Div) + geom_point(aes(x=S_16S, y=S_18S, color=Genotype), alpha=0.9, size=4) + theme_classic()+xlab("Prokaryotes ESV Richness")+ylab("Eukaryotes  ESV Richness")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+geom_smooth(aes(x=S_16S, y=S_18S),method = "lm", se=FALSE, color="black")+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```

###NMDS genotype effect on total microbiome structure - 8 genotypes (Figure 1D)
```{r warning=FALSE}
color= c("#4CB5F5",  "#D70026", "#EDB83D",  "navy",  "#F77604", "#B3C100", "#EC96A4", "#5BC8AC")
p1=ggplot(data=Div, aes(x=NMDS1, y=NMDS2,color=Genotype))+geom_point(shape=5, size=3.5)+theme_classic(base_size = 15)+xlab("NMDS1")+ylab("NMDS2")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p1

```


##Figure 2
###Load the diversity and structure data for the 8 soils (3 wheat genotypes studied) - For detailed script on how to generate the diversity indices and NMDS, see Figure 3
```{r}
Div <- read.table("Diversity&Structure_data_8soils.txt", header=TRUE, check.names = FALSE)
```


###Graph ESV Richness - 8 soils (Figure 2A)
```{r warning=FALSE}
library(Rmisc)
library(ggplot2)
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
Total <- summarySE(Div, measurevar="S", groupvars=c("Soils"), na.rm = TRUE)
Bac <- summarySE(Div, measurevar="S_16S", groupvars=c("Soils"), na.rm = TRUE)
Euk <- summarySE(Div, measurevar="S_18S", groupvars=c("Soils"), na.rm = TRUE)
p2=ggplot() + geom_point(data=Total, aes(x=reorder(Soils, S), y=S, color=Soils), size=5, shape=15) + geom_point(data=Bac, aes(x=Soils, y=S_16S, color=Soils), size=4, shape=5) + geom_point(data=Euk, aes(x=Soils, y=S_18S, color=Soils), size=4, shape=1) + geom_errorbar(data=Total, aes(x=Soils, ymin=S-se, ymax=S+se, color=Soils), width=.2)+ geom_errorbar(data=Bac, aes(x=Soils, ymin=S_16S-se, ymax=S_16S+se, color=Soils), width=.2) + geom_errorbar(data=Euk, aes(x=Soils, ymin=S_18S-se, ymax=S_18S+se, color=Soils), width=.2)+ theme_classic()+xlab("Soils")+ylab("Observed ESVs Richness")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```

###Graph phylogenetic diversity - 8 soils (Figure 2B)
```{r warning=FALSE}
library(Rmisc)
library(ggplot2)
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
Total <- summarySE(Div, measurevar="PD", groupvars=c("Soils"), na.rm = TRUE)
Bac <- summarySE(Div, measurevar="PD_16S", groupvars=c("Soils"), na.rm = TRUE)
Euk <- summarySE(Div, measurevar="PD_18S", groupvars=c("Soils"), na.rm = TRUE)
p2=ggplot() + geom_point(data=Total, aes(x=reorder(Soils, PD), y=PD, color=Soils), size=5, shape=15) + geom_point(data=Bac, aes(x=Soils, y=PD_16S, color=Soils), size=4, shape=5) + geom_point(data=Euk, aes(x=Soils, y=PD_18S, color=Soils), size=4, shape=1) + geom_errorbar(data=Total, aes(x=Soils, ymin=PD-se, ymax=PD+se, color=Soils), width=.2)+ geom_errorbar(data=Bac, aes(x=Soils, ymin=PD_16S-se, ymax=PD_16S+se, color=Soils), width=.2) + geom_errorbar(data=Euk, aes(x=Soils, ymin=PD_18S-se, ymax=PD_18S+se, color=Soils), width=.2)+ theme_classic()+xlab("Soils")+ylab("Faith's Phylogenetic Diversity")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```


###Graph correlation between Prokaryotic and Eukaryotic richness - 8 soils (Figure 2C)
```{r warning=FALSE}
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
p2=ggplot(data=Div) + geom_point(aes(x=S_16S, y=S_18S, color=Soils), alpha=0.9, size=4) + theme_classic()+xlab("Prokaryotes ESV Richness")+ylab("Eukaryotes ESV Richness")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+geom_smooth(aes(x=S_16S, y=S_18S),method = "lm", se=FALSE, color="black")+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```


###NMDS soil effect on total microbiome structure - 8 soils (Figure 2D)
```{r warning=FALSE}
library(ggplot2)
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
p1=ggplot(data=Div, aes(x=NMDS1, y=NMDS2,color=Soils, shape=Genotype))+geom_point(size=3)+theme_classic(base_size = 15)+xlab("NMDS1")+ylab("NMDS2")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p1

```


##Figure 3

#Load merged 16S and 18S SV table
```{r}
SV_use1 <- read.table("ESV_table_16S&18Smerged.txt", header=TRUE, check.names = FALSE, sep="\t")
```

#Subset just european soils
```{r}
Sols=subset(SV_use1, Huit_sols=="8_sols")
SV_use1=subset(Sols, Type=="DNA")
SV_use1=subset(SV_use1, Continent=="Europe")
dim(SV_use1)
```


###Create filtered matrix for NMDS and diversity indices calculations 
```{r}
matrix<-SV_use1[c(14:6278)]
dim(matrix)
matrix_use<-matrix[,colSums(matrix)>=1]
dim(matrix_use)
```

###Ordination: NMDS with Bray-Curtis distances
```{r warning=FALSE, results='hide', message=FALSE}
library(vegan)
NMDS <- metaMDS(matrix_use, distance = "bray", trymax = 100)
NMDS
stressplot(NMDS)

##Extract scores for sites (samples)
NMDSsites=scores(NMDS, display="sites")
SV_use1=cbind(SV_use1,NMDSsites)
```


###Graph effect agricultural pratices in European soils on microbiome structure - NMDS (Figure 3B)
```{r warning=FALSE}
library(ggplot2)
color= c("#4CB5F5", "#D70026")
p1=ggplot(data=SV_use1, aes(x=NMDS1, y=NMDS2,color=Practices, shape=Soils))+geom_point(size=3)+theme_classic(base_size = 15)+xlab("NMDS1")+ylab("NMDS2")+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+scale_color_manual(values=color)+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p1
```


###Calculate Diversity index (S, Shannon, Simpson, Evenness)
```{r}
Shannon <- diversity(matrix_use)
simp <- diversity(matrix_use, "simpson")
# Species richness (S) and Pielou's evenness (J):
S <- specnumber(matrix_use) 
J <- Shannon/log(S)

divtable= cbind(S,Shannon,J,simp)
SV_use1=cbind(SV_use1,divtable)
```

###Graph effect agricultural pratices in European soils on microbiome ESV richness (Figure 3A)
```{r warning=FALSE}
library(Rmisc)
library(ggplot2)
Diversity_stat <- summarySE(SV_use1, measurevar="S", groupvars=c("Type", "Practices", "Soils"), na.rm = TRUE)
color= c("#4CB5F5", "#D70026")
p2=ggplot() + geom_point(data=Diversity_stat, aes(x=Soils, y=S, color=Practices), size=4)+geom_point(data=SV_use1, aes(x=Soils, y=S, color=Practices), alpha=0.8) + geom_errorbar(data=Diversity_stat, aes(x=Soils, ymin=S-se, ymax=S+se, color=Practices), width=.2) + theme_classic()+xlab("Soils")+ylab("Observed ESVs Richness")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```

##Figure 4 - Bubble plot "Most abundant clades"

###Load the relative abundance data of the most abundant clades (Prokaryotes and Eukaryotes)
```{r}
Fam <- read.table("Family_bubbleplot.txt", header=TRUE, check.names = FALSE, sep="\t")
```

###Bubble plot of most abundant clades (Figure 4)
```{r warning=FALSE}
library(ggplot2)
mycolors <- scale_color_manual(values = c("black","pink", "orange", "yellow", "purple"))
p2=ggplot(data=Fam) + geom_point(aes(x=Rel_abund, y=(reorder (Clade, Rel_abund)), size=number_SV, color=Group), alpha=0.8) + theme_classic()+xlab("Relative Abundance")+ylab("Clades")+ theme(axis.title = element_text(color="black", size=10, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(colour="black", size = 10, face = "bold"))+ theme(legend.title = element_text(colour="black", size=10, face="bold"))+
    scale_size_continuous(name = "Number of ESVs", breaks=c(25, 50, 100, 200, 300),
                          limits = c(4, 300),
                          range = c(0, 10) )+ scale_x_continuous(labels = scales::percent_format(accuracy = 1))+mycolors
p2
```

##Figure 5
###Load the diversity and structure data for the 8 soils (3 wheat genotypes studied)
```{r}
Div <- read.table("Diversity&Structure_data_8soils.txt", header=TRUE, check.names = FALSE)
```

### Number of core taxa per sample (Figure 5C)
```{r warning=FALSE}
library(Rmisc)
library(ggplot2)
Diversity_stat <- summarySE(Div, measurevar="core_count", groupvars=c( "Soils"), na.rm = TRUE)
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
p2=ggplot() + geom_point(data=Diversity_stat, aes(x=reorder(Soils, core_count),, y=core_count, color=Soils), size=4)+geom_point(data=Div, aes(x=Soils, y=core_count, color=Soils), alpha=0.8) + geom_errorbar(data=Diversity_stat, aes(x=Soils, ymin=core_count-se, ymax=core_count+se, color=Soils), width=.2) + theme_classic()+xlab("Soils")+ylab("Number of Core Taxa Observed per Sample")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))
p2
```

### Graph relative abundance of core taxa per sample (Figure 5D)
```{r warning=FALSE}
library(Rmisc)
library(ggplot2)
Diversity_stat <- summarySE(Div, measurevar="Core_rel_abund", groupvars=c("Soils"), na.rm = TRUE)
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
p2=ggplot() + geom_point(data=Diversity_stat, aes(x=reorder(Soils, Core_rel_abund), y=Core_rel_abund, color=Soils), size=4)+geom_point(data=Div, aes(x=Soils, y=Core_rel_abund, color=Soils), alpha=0.8) + geom_errorbar(data=Diversity_stat, aes(x=Soils, ymin=Core_rel_abund-se, ymax=Core_rel_abund+se, color=Soils), width=.2) + theme_classic()+xlab("Soils")+ylab("Relative abundance of Core Taxa ")+scale_color_manual(values=color)+ theme(axis.title = element_text(color="black", size=14, face="bold"))+ theme(axis.text = element_text(color="black", size=12, face="bold"))+ theme(legend.text = element_text(colour="black", size = 12, face = "bold"))+ theme(legend.title = element_text(colour="black", size=14, face="bold"))+ scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.8), breaks = c(0, 0.25, 0.5, 0.75))+ annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.4786, ymax = 0.5369, fill = "#AEBD38", alpha = .2, color = NA)+geom_hline(yintercept=0.5077, color="#AEBD38") 
p2
```



##Figure 6
###Heatmap relative abundance of core taxa (Figure 6)
```{r, warning=FALSE, message=FALSE}
library(pheatmap)
library(RColorBrewer)
library(viridis)
metacore<-read.table("Core_taxa_meta_core.txt", header=TRUE, sep = "\t")
matrix<-read.table("18S&16S_SVtable_heatmap_relabund.txt", header=TRUE, sep = "\t", row.names = 1)
meta2<-read.table("metadata_core.txt", header=TRUE)
```

###Scale the rows
```{r}
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
 
data_norm <- t(apply(matrix, 1, cal_z_score))
```

###Heatmap core taxa
```{r, warning=FALSE}
my_sample_col <- meta2$Soils
mat_col <- data.frame(Soils = my_sample_col)

row.names(mat_col) <- colnames(data_norm)
#faire couleur par groupe taxo (fungi...)
my_gene_col=metacore$Group2
mat_row <- data.frame(Group = my_gene_col)
row.names(mat_row) <- row.names(data_norm)
color= c("#8B5742", "#FFA07A", "#CAE2FF", "#1632AF", "#CCC0DA", "#60497A", "#C4D79B", "#6E8B3D")
ann_colors = list(
Soils = c(CAM1="#8B5742",CAM2="#FFA07A", FR1="#CAE2FF",FR2="#1632AF", IT1="#CCC0DA", IT2="#60497A", SEN1="#C4D79B", SEN2="#6E8B3D"),
Group = c( Amoebozoa="grey77", Archaea= "#1f78b4", Archaeplastida="#b2df8a", Bacteria="#1E1E1E", Hacrobia="#ff97a2",Fungi="#ff7f00", Rhizaria="#ffff00", Stramenopiles="darkmagenta", Unassigned="#cab2d6"))

pheatmap(data_norm,  annotation_col = mat_col, annotation_row = mat_row, annotation_colors = ann_colors, fontsize_row = 4,cutree_rows = 2, col=brewer.pal(9,"OrRd"), show_rownames=FALSE, show_colnames = FALSE)
```


##Figure 7

###Load packages
```{r, warning=FALSE, message=FALSE}
library(phyloseq)
library(RColorBrewer) 
library(SpiecEasi) # Network analysis for sparse compositional data  
library(network)
library(intergraph)
library(ggnet)
library(igraph)
library(microbiome)
library(ggpubr)
```

##Import SV table and taxonomy of the core taxa (n=179 ESVs) 
```{r, warning=FALSE}
otu.core <- read.table(file="Core_Taxa_SVtable_for_network.txt", sep='\t', header=TRUE,check.names=FALSE,row.names=1)
taxo.core <- read.table(file="Core_Taxa_taxo_for_network.txt", sep='\t', header=TRUE,check.names=FALSE,row.names=1)
otuall3=as.matrix(t(otu.core))
otuall4=as.matrix(otu.core)
taxonomy=as.matrix(taxo.core)
TAXall=tax_table(taxonomy)
OTUall=otu_table(otuall3,taxa_are_rows=TRUE)
physeq_all = phyloseq(OTUall, TAXall)
physeq_all
```

###Calculate network on the 179 Core taxa using the SpiecEasi package
```{r}
#net.c <- spiec.easi(otuall4, method='mb', lambda.min.ratio=5e-4, nlambda=60,icov.select.params=list(rep.num=99, ncores=7))

#export
#saveRDS(net.c, "network_core16&18S_final.rds")

```

###The output of spiec.easi is stored as network_core16&18S_final.rds
```{r}
net.c <- readRDS("network_core16&18S_final.rds")
class(net.c)
n.c <- symBeta(getOptBeta(net.c))

```

####Check stabilty of network (need to be close to 0.05). If not This problem we can be fixed by lowering lambda.min.ratio to explore denser networks.To get closer to the mark, we should bump up nlambda to more finely sample the lambda path, which gives a denser network.
```{r}
getStability(net.c)
getOptLambda(net.c)
sum(getRefit(net.c))/2
```

###Prepare data for plotting using igraph package
```{r}
#Add names to IDs
#We also add abundance values to vertex (=nodes=SVs).

colnames(n.c) <- rownames(n.c) <- colnames(otuall4)
# add log abundance as properties of vertex/nodes
vsize <- log2(apply(otuall4, 2, mean))

ig <- graph.adjacency(n.c, mode='undirected', add.rownames = TRUE, weighted = TRUE)
ig # we can see all the attributes and weights
```

##Calculate centrality index for each node (ESV)
```{r, warning=FALSE, message=FALSE}
library(qgraph)
library(erer)
allindex=centrality_auto(ig, weighted=TRUE)

#print(allindex)
```


##Graph relationship Node degree with Betweeness centrality or with Closeness centrality to identify hub/keystone species
```{r}
index=read.table(file="Network_index_core_16S&18Sfinal.txt", header=TRUE)
taxoall3 <- read.table(file="Core_Taxa_taxo_for_network.txt", sep='\t', header=TRUE,check.names=FALSE)

index_tax=merge(index, taxoall3, by="OTUID")
```


###Graph Node degree vs Betweeness centrality (Figure 7C)
```{r warning=FALSE}
library(ggplot2)
mycolors <- scale_color_manual(values = c("#a6cee3",  "#b2df8a","#1f78b4","#1E1E1E","#000075","#33a02c","#fdbf6f","#ff7f00","#cab2d6","#ffff99","#b15928",'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000',"#800000","#8E388E","#7171C6","#7D9EC0","#388E8E","#71C671","#8E8E38","#C5C1AA", "#C67171","#555555", "orange"))
myshape <- scale_shape_manual(values=c(17, 16))
p1=ggplot(data=index_tax, aes(x=Degree, y=Betweenness,color=Class, shape=Hub.x))+geom_point(size=3)+theme_classic(base_size = 8)+xlab("Node Degree (nb Correlations)")+ylab("Betweenness Centrality ")+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(colour="black", size = 6, face = "bold"))+ theme(legend.title = element_text(colour="black", size=8, face="bold"))+mycolors+myshape
p1

```

###Graph Node degree vs Closeness centrality (Figure 7D)
```{r warning=FALSE}
library(ggplot2)
mycolors <- scale_color_manual(values = c("#a6cee3",  "#b2df8a","#1f78b4","#1E1E1E","#000075","#33a02c","#fdbf6f","#ff7f00","#cab2d6","#ffff99","#b15928",'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000',"#800000","#8E388E","#7171C6","#7D9EC0","#388E8E","#71C671","#8E8E38","#C5C1AA", "#C67171","#555555", "orange"))

p1=ggplot(data=index_tax, aes(x=Degree, y=Closeness,color=Class, shape=Hub.x))+geom_point(size=3)+theme_classic(base_size = 8)+xlab("Node Degree (nb Correlations)")+ylab("Closeness Centrality ")+ theme(axis.title = element_text(color="black", size=12, face="bold"))+ theme(axis.text = element_text(color="black", size=10, face="bold"))+ theme(legend.text = element_text(colour="black", size = 6, face = "bold"))+ theme(legend.title = element_text(colour="black", size=8, face="bold"))+mycolors+myshape
p1

```

###Improved network visualization with ggnet package, add colors to edges
```{r}
net <- asNetwork(ig)
network::set.edge.attribute(net, "color", ifelse(net %e% "weight" > 0, "lightblue", "red"))
```


###add information on taxonomy, hub, nodesize for each core taxon in the network
```{r}
class <- map_levels(colnames(otuall4), from = "OTUID2", to = "Supergroup", tax_table(physeq_all))
net %v% "Supergroup" <- class
cluster <- map_levels(colnames(otuall4), from = "OTUID2", to = "Cluster", tax_table(physeq_all))
net %v% "Cluster" <- cluster
net %v% "nodesize" <- vsize

hub2 <- index_tax$Hub.x
hub2 <-  as.character(hub2)


net %v% "Hub" <- hub2
```


###Network graph - nodes colored by Taxonomic Group (Figure 7A)
```{r, warning=FALSE, message=FALSE}
mycolors <- scale_color_manual(values = c("grey77",  "#1f78b4","#b2df8a","#1E1E1E","#C67171","#ff7f00","pink","goldenrod1","darkmagenta","#cab2d6","#ffff99","#b15928",'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000',"#800000","#8E388E","#7171C6","#7D9EC0","#388E8E","#71C671","#8E8E38","#C5C1AA", "#C67171","#555555", "orange"))
myshape <- scale_shape_manual(values=c(17, 16))

p <- ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), node.color = "Supergroup", label = FALSE, node.size = "nodesize", edge.color = "color", node.shape = "Hub", max_size = 4, legend.size = 12) + guides(color=guide_legend(title="Group"), size = FALSE, mode = c("x", "y")) + mycolors +myshape
p 
```

###Network graph - nodes colored by core taxa cluster (Figure 7B)
```{r, warning=FALSE, message=FALSE}
mycolors <- scale_color_manual(values = c("darkblue",  "darkred","gray"))
myshape <- scale_shape_manual(values=c(17, 16))

p <- ggnet2(net, mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75), node.color = "Cluster", label = FALSE, node.size = "nodesize", edge.color = "color", node.shape = "Hub", max_size = 4, legend.size = 12) + guides(color=guide_legend(title="Cluster"), size = FALSE, mode = c("x", "y")) + mycolors+myshape
p 
```














