#!/usr/bin/env R
library(ape)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(microshades)
library(speedyseq)
library(cowplot)
min_read_count=1000

###STEP2: Import Mapping file (metadate file)
#1.Check mapping file before import to R, R doesn’t seem to like sample name to start with number or contain “-” in sample name. If you get error in this step, you should check file name first.

#2.First column of first row should not start with #, R will not read the first row that starts with #

metatable = file.path("sample_info.tsv")
meta = read.table(metatable,header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)
rownames(meta) <- gsub('-','.',rownames(meta))

meta <- meta[which(meta$DemuxReads > min_read_count),]

###STEP3: Check if your metadata file has been import successfully and correctly

# The output will show a table of your metadata file (mapping file).

#*If you do not have header, you might have started your first row with #

head(meta)

###STEP4: Construct sample_data-class using imported metadata

sampleData <- sample_data(meta)

###STEP5: Import OTU table

#OTU table from  ITS data is “MC2017FC.otu_table.txt”.

otutable <- file.path("ECDRE_LR_Y0_1_ITS.otu_table.txt")
otus <- read.table(otutable,header=T,sep="\t",row.names=1)
otumat <- as(as.matrix(otus), "matrix")
OTU = otu_table(otumat, taxa_are_rows = TRUE)

#Check imported OTU table

head(OTU)

###STEP6: Import taxonomy table
#Taxonmy table generated from AMPtk need to be rearranged using following script.

#perl rdp_taxonmy2mat.pl<Input_taxonmy.txt>Output_taxonomy.txt
#rdp_taxonomy2mat.pl was created by Professor Jason E. Stajich

taxin = file.path("ECDRE_LR_Y0_1_ITS.cluster.taxonomy.fix.txt")
taxmat <- read.table(taxin, header=T,sep="\t",row.names=1)
taxmat <- as(as.matrix(taxmat),"matrix")
TAX = tax_table(taxmat)

###STEP7: Import phylogenetic tree
#Phylogenetic tree can also be include for further phylogenetic analysis.

treefile = file.path("ECDRE_LR_Y0_1_ITS.cluster.tree.phy")
tree = read.tree(treefile)

###STEP8: Construct Phyloseq object
#To construct phyloseq object, otu table, taxonomy table, and sampleData are required. Phylogenetic tree can be included, but it is not necessary for constructing phyloseq object.
#Construct Phyloseq object called "Physeq"

physeq = phyloseq(OTU,TAX,sampleData,tree)

#Check phyloseq object
#This should indicate that your physeq is a "phyloseq-class experiment-level object""

physeq

###STEP9: Remove singletons
#Remove any OTUs that present only one time.

physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)
physeq.prune = subset_samples(physeq.prune, Tissue_type != "Control")

physeq.prune

###STEP10: Plot read counts to check dataset
#Check read counts: any samples that have very low reads should be removed.
#[Ref](http://evomics.org/wp-content/uploads/2016/01/phyloseq-Lab-01-Answers.html)

readcount = data.table(as(sample_data(physeq.prune), "data.frame"),
                 TotalReads = sample_sums(physeq.prune), 
                 keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")
#For plotting, use command below.
#SeqDepth = ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

#TotalReads of all the samples can be in this table (select only SampleID and TotalReads columns).
#In order to check samples with low number of reads, "order()" can be used to sort "TotalReads" column.
#In this dataset, N55.Rhizo has very low number of reads, so will will filter this sample out using the next minimum number of reads.
readcount = readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")]

head(readcount)

set.seed(1)
physeq.prune.rarefy = rarefy_even_depth(physeq.prune, sample.size = 2000, replace = FALSE, trimOTUs = TRUE)
physeq.prune.rarefy


###STEP11.1: Plot Alpha diversity by Crust_type
#Alpha diversity can be Chao1, Observed, Shannon, Simpson
#This plot include statistical analysis using "stat_compare_means" with "method = anova"

###STEP12.1 Taxonomic composition

# remove phyla we don't want to see
m = subset_taxa(physeq.prune.rarefy, Phylum != "Cercozoa")
m = subset_taxa(m, Phylum != "Choanoflagellata")
m = subset_taxa(m, Phylum != "Anthophyta")
physeq.prune.rarefy = m
FPhylum = as.character(get_taxa_unique(physeq.prune.rarefy,"Phylum"))
FPhylum = FPhylum[complete.cases(FPhylum)]
FPhylum


p<-plot_richness(physeq.prune.rarefy, x = "Tissue_type", color = "Tissue_type", measures = c("Shannon")) + 
  geom_boxplot() +
  theme_bw() + ggtitle("Alpha diversity plot by Tissue (Shannon)") + 
  theme(plot.title = element_text(hjust = 0.5),
         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_point(size=5, alpha=0.8)

ggsave("Figures/Shannon_diversity_by_Tissue_type.pdf",p,width=8,height=8)

p<-plot_richness(physeq.prune.rarefy, x = "Field", color = "Field", measures = c("Shannon")) + 
  geom_boxplot() +
  theme_bw() + ggtitle("Alpha diversity plot by Field (Shannon)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_point(size=5, alpha=0.8)

ggsave("Figures/Shannon_diversity_by_Field.pdf",p,width=8,height=8)

p<-plot_richness(physeq.prune.rarefy, x = "Location", color = "Location", measures = c("Shannon")) + 
  geom_boxplot() +
  theme_bw() + ggtitle("Alpha diversity plot by Location (Shannon)") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_point(size=5, alpha=0.8)

ggsave("Figures/Shannon_diversity_by_Location.pdf",p,width=8,height=8)

###STEP12:Taxonomy Barplot

#Taxonomy barplot by sample (Phylum level)

p<-plot_bar(physeq.prune.rarefy, x = "Sample", y = "Abundance", fill ="Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  ggtitle("Taxonomy barplot by sample (Phylum)") + scale_colour_brewer(palette="Set1") +
  theme(plot.title = element_text(hjust = 0.5))

pdf("Figures/Abundance_by_phylum.pdf",width=24)
p
dev.off()
#ggsave("Figures/Abundance_by_phylum.pdf",p,width=24,height=18)

#### make some phyloseq taxon barplot plots prettier wth microshades

mdf_prep <- prep_mdf(physeq.prune.rarefy)
# see https://karstenslab.github.io/microshades/articles/microshades-GP.html

color_objs_TR <- create_color_dfs(mdf_prep,
                                  selected_groups = c("Glomeromycota","Mucoromycota","Ascomycota", "Basidiomycota") , cvd = TRUE)

# Extract
mdf_TR <- color_objs_TR$mdf
cdf_TR <- color_objs_TR$cdf

plot <- plot_microshades(mdf_TR, cdf_TR)

plot_1 <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 4)) 

plot_1 
ggsave("Figures/Taxbarplot_All.pdf",plot_2,width=32,height=18)

plot_2 <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~Tissue_type, scales = "free_x", nrow = 2) +
  theme (strip.text.x = element_text(size = 6))

ggsave("Figures/Taxbarplot_Tissue_type.pdf",plot_2,width=24,height=18)

new_groups <- extend_group(mdf_GP, cdf_GP, "Phylum", "Genus", "Ascomycota", existing_palette = "micro_cvd_orange", new_palette = "micro_orange", n_add = 5)

TR_legend <-custom_legend(new_groups$mdf, new_groups$cdf)

plot_diff <- plot_microshades(new_groups$mdf, new_groups$cdf) + 
 scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~Tissue_type + Location, scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_text(size= 6)) + 
  theme(plot.margin = margin(6,20,6,6))


p<-plot_grid(plot_diff, TR_legend,  rel_widths = c(1, .25))
ggsave("Figures/ITS_abundances_by_Tissue_Location.pdf",p,width=32,height=18)

plot_diff <- plot_microshades(new_groups$mdf, new_groups$cdf) + 
 scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 8)) +
  facet_wrap(~Tissue_type + Field, scales = "free_x", nrow = 2) +
  theme(axis.text.x = element_text(size= 8)) + 
  theme(plot.margin = margin(6,20,6,6))

p<-plot_grid(plot_diff, TR_legend,  rel_widths = c(1, .25))
ggsave("Figures/ITS_abundances_by_Tissue_Field.pdf",p,width=32,height=18)


# beta diversity plots
physeq.prune.rarefy.ps.ord <- ordinate(physeq.prune.rarefy, "PCoA", "unifrac")
p<-plot_ordination(physeq.prune.rarefy, physeq.prune.rarefy.ps.ord, type = "samples", shape = "Tissue_type",
                color="Field")  + ggtitle("Fungi Beta Diversity (PCoA) by Tissue Type / Field") + 
  scale_colour_brewer(palette="Set1") + theme_cowplot(12)

ggsave("Figures/Betadiversity_PCoA_Tissue_Field.pdf",p,width=8,height=8)

p2 <- p + facet_wrap(~Field, 3)
ggsave("Figures/Betadiversity_PCoA_Tissue_Field_facet.pdf",p2,width=8,height=8)

taxcom_treatment = ggplot(data = psmelt(physeq.prune.rarefy), mapping = aes_string(x = "Tissue_type",
                                                                                   y = "Abundance", 
                                                                                   fill = "Phylum" )) + 
  geom_bar(stat="identity", position="fill") + 
  ggtitle("Fungi from Tissue type in Year 0-1")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
                     ) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Figures/Taxcom_Tissue.pdf",taxcom_treatment,plot_amplicon_simple.R)
#print(taxcom_treatment)

###Plot Alpha Diversity (Species richness)

physeq.prune.rarefy.plot.richness = plot_richness(physeq.prune.rarefy, x="Field", color=("Tissue_type"),
                                                  measures=c("Observed")) + geom_boxplot(lwd=0.5) +
  ggtitle("Alpha Diversity by Tissue_type") + stat_compare_means(method = "t.test", label.x.npc = c("right")) + 
  scale_colour_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  theme(plot.title = element_text(hjust = 0.5)) + theme_cowplot(12)
#


ggsave("Figures/Alphadiversity_Field_Treatment.pdf",physeq.prune.rarefy.plot.richness,width=12,height=12)


