##########################################################################################################################################
##########################################################################################################################################

## Life-stage specificity and cross generational climate effects on the microbiome of a tropical sea urchin (Echinodermata: Echinoidea) ##

##########################################################################################################################################
##########################################################################################################################################

##EMMA MARANGON

##########################################################
#################### PRE-PROCESSING ####################
##########################################################

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library (decontam) #identify contaminants 
library(scales)   #to change y axis tick marks to percentage
library(forcats) #to change order groups axix x (mutate function)
library(RColorBrewer) #to change colours graphs
library(RVAideMemoire) #post hoc tests beta diversity
library (DESeq2) #differentially abundant ASVs
library(glmmTMB) #linear mixed models for alpha diversity
library(DHARMa) #check assumptions linear mixed models
library(emmeans) ##post hoc tests alpha diversity
library("MicEco") #venn diagram
library(fct_relevel)
library(mixOmics)
library (btools) #calculate Faith's phylogenetic diversity


############################################################
### load data ####

##1 ASV TABLE
otu_table = read.csv("urchins-feature-table_R.txt", sep = '\t',  dec = ".", check.names = FALSE, row.names=1)

##2 TAXA TABLE
otu_matrix = read.csv("taxonomy_R.tsv", sep = '\t', header=T, row.names=1)
TAXA_TABLE_split<-otu_matrix %>% separate(Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=TRUE)
TAXA_TABLE_matrix<-as.matrix(TAXA_TABLE_split)

##3 METADATA
metadata9 = read.csv("urchins-metadata9.txt", sep = '\t') 


############################################################
### import ###
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table (TAXA_TABLE_matrix)
meta = sample_data(metadata9)
phy_tree = read_tree("tree.nwk") #rooted tree

phyloseq_merged = phyloseq(OTU, TAX, meta, phy_tree) #merging


############################################################
### FILTERING ###

### filtering eukaryota, mitochondria, chloroplasts ###
phyloseq_merged_filtered_yesGut <- phyloseq_merged %>% subset_taxa(Domain != "D_0__Eukaryota" & Family  != "D_4__Mitochondria" & Order   != "D_3__Chloroplast")

### excluding samples not belonging to this study ####
phyloseq_merged_filtered <- phyloseq_merged_filtered_yesGut %>% subset_samples(sampleType != "gut")
phyloseq_merged_filtered_NoBlanks <- subset_samples(phyloseq_merged_filtered, sampleType!="blank")
phyloseq_merged_filtered_NoBlanksNonegatives <- subset_samples(phyloseq_merged_filtered_NoBlanks, sampleType!="negative")


############################################################
### removing contaminants using decontam ###

##### 1. inspect library size
plot_richness(phyloseq_merged_filtered, x = "treatment", color = "sampleType") + geom_boxplot()
df <- as.data.frame(sample_data(phyloseq_merged_filtered))
df$LibrarySize <- sample_sums(phyloseq_merged_filtered)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sampleType)) + geom_point()

##### 2. identify contaminants - prevalence
sample_data(phyloseq_merged_filtered)$is.neg <- sample_data(phyloseq_merged_filtered)$sampleType == "blank"
contamdf.prev <- isContaminant(phyloseq_merged_filtered, method="prevalence", neg="is.neg") #defaut threshold = 0.1, and $contaminant=TRUE if $p < 0.1
contamdf.prev_0.5 <- isContaminant(phyloseq_merged_filtered, method="prevalence", neg="is.neg", threshold=0.5) #more aggressive classification threshold: threshold = 0.5
table(contamdf.prev$contaminant) #FALSE 45118 #TRUE 17 -> 17 ASVs have been identified as contaminants
table(contamdf.prev_0.5$contaminant) #FALSE 45103 #TRUE 32 -> 32 ASVs have been identified as contaminants 

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_filtered, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sampleType == "blank", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sampleType != "blank", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev_0.5$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Blanks)") + ylab("Prevalence (All Samples but blanks)")


### removing contaminants ###
taxa_are_rows(phyloseq_merged_filtered) #TRUE
contaminants <- subset(contamdf.prev_0.5, contaminant == "TRUE") #I am using the strict threshold
contaminants_only <- c("420e7973dcc614c06e823acd64399046", "e88e9cbc6d53354d2de24dffc7a5e93c", "47528e332c77ba4cef8d067390f8b3b1", "bbe968da0a6da2ce9f66c7ab5ea15425",
                       "8424ac265142fe201c78bc93949a1ed8", "baf7e5d44d4a61f89c1bdc6ab7a446c4", "b47dcf513fb6be7c020be3216ec4a639", "de8144e92fbe0d32a5a90b8c1133f23b",
                       "7a71f57deebcdbc6796cbab3a587bd5a", "fd98346394a5c79e554003012cb33826", "4668ae4998e1edc9996f75e7c09f7bac", "71cef5b1174477c4ffc4ad09a21f4c63",
                       "ebd0562fa73b889281c0ded93518b7fa", "0b80a2392f611b736d856d177d560bf8", "118735e135c14c9c69b787423d51ac40", "0c0044e3d1acdbb71c7cf4bc581bd27e",
                       "afb1de4000e2688b434dacab00c9ae4c", "5f56ed10a6d42d666669971daba052da", "c8b083c2a3aef63ed75421f998aaa72f", "529ff5c6b52a153757f4437c5722d155",
                       "06a6e5477a6e3b608494ee0adccb4a94", "75c904abd5443e40e28c688641317b7c", "ecb39484d5cfb68cd99bd33debef4fd4", "550a0f69768245c3debe5262c3c6d668",
                       "511a1dc820a661caeb441b3e5a55fb2f", "39385ee32866666be718160857c5bcb6", "a166af5a482fb08f9ffbdcc3cc588e83", "b99be30a570f9fb32fff706db4ae4609",
                       "8c44c6e0e7f7629878fe6666c139ef3d", "a93e5ddc25eca21fbfbdb99d40634ca7", "be066e42c37b89db6680195009e0f55b", "7d217addfd80cb5393d81dae5f160367")
AllTaxa = taxa_names(phyloseq_merged_filtered_NoBlanksNonegatives)
AllTaxa <- AllTaxa[!(AllTaxa %in% contaminants_only)]
phyloseq_merged_filtered_NoBNC = prune_taxa(AllTaxa, phyloseq_merged_filtered_NoBlanksNonegatives)
summary(sample_sums(phyloseq_merged_filtered_NoBNC)) # Between 0 and 104419 reads per sample, mean = 33980
dfnoBCN <- as.data.frame(sample_data(phyloseq_merged_filtered_NoBNC))
dfnoBCN$LibrarySize <- sample_sums(phyloseq_merged_filtered_NoBNC)
dfnoBCN <- dfnoBCN[order(dfnoBCN$LibrarySize),]
dfnoBCN$Index <- seq(nrow(dfnoBCN))
ggplot(data=dfnoBCN, aes(x=Index, y=LibrarySize, color=sampleType)) + geom_point() 


############################################################
### MORE FILTERING ###

### exclude EMBRYOS and EGGS from this study ###
phyloseq_merged_FINAL3 <- phyloseq_merged_filtered_NoBNC %>% subset_samples(sampleType != "embryos" & sampleType !="eggs") 

### Exclude samples < 8400 reads ###
phyloseq_merged_FINAL2 <- prune_samples(sample_sums(phyloseq_merged_FINAL3) > 8400, phyloseq_merged_FINAL3) 

### Exlude transplant gonads F1 and larvae F2 from this study ###
phyloseq_merged_FINAL2.2 <- phyloseq_merged_FINAL2 %>% subset_samples(sample_treatment!= "gonad_MMH" & sample_treatment !="seawater_MMH" & sample_treatment !="larvae_MMHH_first" & sample_treatment != "larvae_MMHH_second" )  #to remove samples I don't want

### keep only taxa that were observed at least twice (=remove singletons) ###
phyloseq_merged_FINAL <- prune_taxa(taxa_sums(phyloseq_merged_FINAL2.2) >=2, phyloseq_merged_FINAL2.2)

### excluding samples not relevant here ###
phyloseq_merged_FINAL_filtered <- subset_samples(phyloseq_merged_FINAL, sample_date != "gonad_29.11.17" & sample_sex != "gonad_male" & sample_generation != "larvae_F2_dead")
phyloseq_merged_FINAL_filtered = prune_taxa(taxa_sums(phyloseq_merged_FINAL_filtered) > 0, phyloseq_merged_FINAL_filtered) #remove unobserved ASVs (sum 0 across all samples)
#I checked for singletons in the database after this last filtering - and no singletons found - all good


############################################################
### RAREFACTION CURVES BY SAMPLE TYPE ###

phyloseq_merged_FINAL_sw <- subset_samples(phyloseq_merged_FINAL_filtered, sampleType == "seawater")
rarecurve(t(otu_table(phyloseq_merged_FINAL_sw)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#C77CFE")

phyloseq_merged_FINAL_gonad <- subset_samples(phyloseq_merged_FINAL_filtered, sampleType == "gonad")
rarecurve(t(otu_table(phyloseq_merged_FINAL_gonad)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#F8766C")

phyloseq_merged_FINAL_larvae <- subset_samples(phyloseq_merged_FINAL_filtered, sampleType == "larvae") 
rarecurve(t(otu_table(phyloseq_merged_FINAL_larvae)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#00BFC3")

phyloseq_merged_FINAL_juv <- subset_samples(phyloseq_merged_FINAL_filtered, sampleType == "juvenile")
rarecurve(t(otu_table(phyloseq_merged_FINAL_juv)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#7CAD00")


############################################################
### INFO FINAL DATASET ###
sum(sample_sums(phyloseq_merged_FINAL_filtered)) # 4767687 reads
summary(sample_sums(phyloseq_merged_FINAL_filtered)) # Between 8495 and 104419 reads per sample, mean = 35580
phyloseq_merged_FINAL_filtered # 40599 ASVs for 134 samples


############################################################
### NORMALIZATION ###

####### PROPORTIONS
phyloseq_merged_FINALabundances = transform_sample_counts(phyloseq_merged_FINAL_filtered, function(x){x / sum(x)}) 
FINALabundances_cutoff <- filter_taxa(phyloseq_merged_FINALabundances, function(x) mean(x) > 1e-5, TRUE) #keeping only ASvs with mean > 0.00001; from 40599 to 6187 taxa
FINALabundances_cutoff = prune_taxa(taxa_sums(FINALabundances_cutoff) > 0, FINALabundances_cutoff) #to make sure no ASVs with sum 0 across all samples 

####### RAREFACTION
rarefied = rarefy_even_depth(phyloseq_merged_FINAL, rngseed=1, sample.size = min(sample_sums(phyloseq_merged_FINAL)))
FINALrarefied <- subset_samples(rarefied, sample_date != "gonad_29.11.17" & sample_sex != "gonad_male" & sample_generation != "larvae_F2_dead")
FINAL_rarefied = prune_taxa(taxa_sums(FINALrarefied) > 0, FINALrarefied)  #remove unobserved ASVs (sum 0 across all samples) 


########################################################################################################################
########################################################################################################################






####################################################################################################################
################## Microbial dynamics across life stages and generations under ambient conditions ##################
####################################################################################################################

############################################################
### FILTERING ###

FINAL_abundances_Q1_cutoff <- subset_samples(FINALabundances_cutoff, climateFINAL == "C") #only ambient
FINAL_abundances_Q1_cutoff = prune_taxa(taxa_sums(FINAL_abundances_Q1_cutoff) > 0, FINAL_abundances_Q1_cutoff) #remove unobserved ASVs (sum 0 across all samples)

rarefied_Q1_NOcutoff <- subset_samples(FINAL_rarefied, climateFINAL == "C") #only ambient
rarefied_Q1_NOcutoff = prune_taxa(taxa_sums(rarefied_Q1_NOcutoff) > 0, rarefied_Q1_NOcutoff) #remove unobserved ASVs (sum 0 across all samples)


############################################################
### SQRT TRANSFORMATION ###
FINAL_abundances_Q1_sqrt <- transform_sample_counts(FINAL_abundances_Q1_cutoff, function (x) sqrt(x))
FINAL_abundances_Q1_sqrt = prune_taxa(taxa_sums(FINAL_abundances_Q1_sqrt) > 0, FINAL_abundances_Q1_sqrt) #remove unobserved ASVs (sum 0 across all samples)


############################################################
### NMDS ###
nmds_Q1 <- ordinate(FINAL_abundances_Q1_sqrt, "NMDS", "bray") %>%
  plot_ordination(FINAL_abundances_Q1_sqrt, ., color = "sampleTypeAll", shape = "generation", title = "Fig2_nmds_brays_Q1") + 
  geom_point(size = 4.5) +
  scale_color_discrete(name="sample type",
                                 breaks=c("gonad", "juvenile", "larvae_first", "larvae_second" , "seawater"),
                                 labels=c("Adult gonad", "Juvenile", "Larvae 1d", "Larvae > 1d", "Seawater")) +
								 theme_bw() + theme( 
								 plot.background = element_blank(),
								 panel.grid.major = element_blank(),
								 panel.grid.minor = element_blank()) + 
								 theme(legend.text = element_text(size = 15),legend.title = element_text(size = 16)) 
nmds_Q1

###stress nmds
nmds_Q1_stress <- ordinate(FINAL_abundances_Q1_sqrt, "NMDS", "bray")
cat("Stress:", nmds_Q1_stress$stress, fill=TRUE) 


############################################################
### BAR PLOT AT CALSS LEVEL ###
Q1_Cl <- tax_glom(FINAL_abundances_Q1_cutoff, taxrank = 'Class') #agglomerate data (merging taxa of the same Class)
q1_cl<- prune_taxa(taxa_sums(Q1_Cl) > 0, Q1_Cl)
Q1_Gcl_sample <- merge_samples(q1_cl, "sample_generation", fun=mean)  #I create a new group named Sample with my sample_generation categories, and I associated mean class per group 
Q1_Gcl_sample_abund = transform_sample_counts(Q1_Gcl_sample, function(x){x / sum(x)})
Q1_Gcl_sample_abund_melt <- psmelt(Q1_Gcl_sample_abund)
Q1_Gcl_sample_abund_melt$Class <- as.character(Q1_Gcl_sample_abund_melt$Class) # convert Class to a character vector from a factor
Q1_Gcl_sample_abund_melt$Class[Q1_Gcl_sample_abund_melt$Abundance < 0.01] <- "Other" #rename classes with < 1% abundance
Q1_Gcl_sample_abund_melt <- Q1_Gcl_sample_abund_melt %>% mutate(Sample = fct_relevel(Sample, "gonad_F0", "gonad_F1", "larvae_F1_first", "larvae_F2_first", "larvae_F1_second", "larvae_F2_second","juvenile_F1", "seawater")) #reorder
p <- ggplot(data=Q1_Gcl_sample_abund_melt, aes(x=Sample, y=Abundance, fill=Class))
nb.cols <- 46
mycolors <- colorRampPalette(brewer.pal(46, "Set3"))(nb.cols)
p + geom_bar(aes(), stat="identity", position="stack")+ theme_bw() +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + guides(fill=guide_legend(nrow=5)) + labs(title = "Q1_C_Gcl") +scale_fill_manual(values = mycolors) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) + 
  labs (y = "Mean relative abundance (%)") + 
  scale_y_continuous(labels=percent) 


############################################################
### BUBBLE PLOT DOMINANT FAMILIES ###
Q1_Fam <- tax_glom(FINAL_abundances_Q1_cutoff, taxrank = 'Family') #agglomerate data (merging taxa of the same Family)
q1_Fam<- prune_taxa(taxa_sums(Q1_Fam) > 0, Q1_Fam) #to delete taxa which abundance sum is 0 
Q1_Gfam_sample <- merge_samples(q1_Fam, "sample_generation", fun=mean) 
Q1_Gfam_sample_abund = transform_sample_counts(Q1_Gfam_sample, function(x){x / sum(x)}) 
Q1_Gfam_sample_abund_melt <- psmelt(Q1_Gfam_sample_abund)
Q1_Gfam_sample_abund_melt$Family <- as.character(Q1_Gfam_sample_abund_melt$Family)
#selecting families >5% relative abundance:
Q1_Fam5perc <- subset (Q1_Gfam_sample_abund_melt, Family == "D_4__Alteromonadaceae" | Family == "D_4__Vibrionaceae"| Family ==  "D_4__Saprospiraceae" | Family == "D_4__Saccharospirillaceae" | Family == "D_4__Rhodobacteraceae" | Family ==  "D_4__Prolixibacteraceae"| Family == "D_4__P3OB-42"| Family == "D_4__Oleiphilaceae" | Family == "D_4__Nitrosopumilaceae" | Family == "D_4__Nitrincolaceae" | Family == "D_4__Methyloligellaceae"| Family == "D_4__Kiritimatiellaceae" | Family == "D_4__Halomonadaceae" | Family == "D_4__Fusobacteriaceae"| Family == "D_4__Flavobacteriaceae"| Family == "D_4__Desulfobulbaceae"| Family == "D_4__Cryomorphaceae"| Family == "D_4__Cellvibrionaceae") 
p <- Q1_Fam5perc %>% mutate(name = fct_relevel(Family, "D_4__Vibrionaceae", "D_4__Saccharospirillaceae", "D_4__Oleiphilaceae", "D_4__Nitrincolaceae", "D_4__Halomonadaceae","D_4__Cellvibrionaceae", "D_4__Alteromonadaceae", "D_4__Saprospiraceae","D_4__Prolixibacteraceae","D_4__Flavobacteriaceae","D_4__Cryomorphaceae", "D_4__Rhodobacteraceae", "D_4__Methyloligellaceae",  "D_4__P3OB-42", "D_4__Desulfobulbaceae", "D_4__Fusobacteriaceae",  "D_4__Nitrosopumilaceae", "D_4__Kiritimatiellaceae")) 
p2 <- p %>% mutate(name2 = fct_relevel(Sample, "gonad_F0", "gonad_F1", "larvae_F1_first", "larvae_F2_first", "larvae_F1_second","larvae_F2_second", "juvenile_F1", "seawater")) 
Q1_bubbleFamily <- ggplot(p2, aes(x = name2, y = name, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=11), legend.title=element_text(size=13), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=30,  face="bold"), axis.text.y=element_text(size=30)) +
  labs(title = "Q1_C_plot_fam") +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank())
Q1_bubbleFamily


############################################################
### BAR PLOT DOMINANT FAMILIES ###

#D_4__Rhodobacteraceae - D_2__Alphaproteobacteria
Rhodobacteraceae <- subset(Q1_Gfam_sample_abund_melt, Family == "D_4__Rhodobacteraceae")
Rhodobacteraceae2 <- Rhodobacteraceae %>% mutate(name2 = fct_relevel(Sample, "gonad_F0", "gonad_F1", "larvae_F1_first", "larvae_F2_first","larvae_F1_second", "larvae_F2_second", "juvenile_F1", "seawater")) #to reoder x axis
ONEfamily_plot <- ggplot(data=Rhodobacteraceae2, aes(x=name2, y=Abundance))
ONEfamily_plot + geom_bar(aes(), stat="identity", position="stack", fill="#A8DDC1") + theme_bw() +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + guides(fill=guide_legend(nrow=5)) + labs(title = "Q1_Rhodobacteraceae_C ") + ylim(0, 0.32)

#D_4__Vibrionaceae - D_2__Gammaproteobacteria
Vibrionaceae <- subset(Q1_Gfam_sample_abund_melt, Family == "D_4__Vibrionaceae")
Vibrionaceae2 <- Vibrionaceae %>% mutate(name2 = fct_relevel(Sample, "gonad_F0", "gonad_F1", "larvae_F1_first", "larvae_F2_first","larvae_F1_second", "larvae_F2_second", "juvenile_F1", "seawater")) 
ONEfamily_plot <- ggplot(data=Vibrionaceae2, aes(x=name2, y=Abundance))
ONEfamily_plot + geom_bar(aes(), stat="identity", position="stack", fill="#C0BCD7") + theme_bw() +
    theme(legend.position="bottom", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + guides(fill=guide_legend(nrow=5)) + labs(title = "Q1_Vibrionaceae_C ") + ylim(0, 0.32)

#D_4__Alteromonadaceae - D_2__Gammaproteobacteria 
Alteromonadaceae <- subset(Q1_Gfam_sample_abund_melt, Family == "D_4__Alteromonadaceae")
Alteromonadaceae2 <- Alteromonadaceae %>% mutate(name2 = fct_relevel(Sample, "gonad_F0", "gonad_F1", "larvae_F1_first", "larvae_F2_first","larvae_F1_second", "larvae_F2_second", "juvenile_F1", "seawater")) 
ONEfamily_plot <- ggplot(data=Alteromonadaceae2, aes(x=name2, y=Abundance))
ONEfamily_plot + geom_bar(aes(), stat="identity", position="stack", fill="#C0BCD7") + theme_bw() +
	  theme(legend.position="bottom", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + guides(fill=guide_legend(nrow=5)) + labs(title = "Q1_Alteromonadaceae_C ") + ylim(0, 0.32)




############################################################
### STATS BETA DIVERSITY ###

###sampleType*generation ###

Q1_adonis = as (sample_data(FINAL_abundances_Q1_sqrt), "data.frame")
Q1_d = phyloseq::distance(FINAL_abundances_Q1_sqrt,'bray') 
Adonis_Q1 <-adonis2(Q1_d ~ tank + sampleType*generation, data=Q1_adonis,  permutations = 10000, method = "bray") 
Adonis_Q1 

## dispersion ##
beta <- betadisper(Q1_d, sample_data(FINAL_abundances_Q1_sqrt)$sampleType)
plot(beta)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # p<0.05 - NOT OK
TukeyHSD (beta)

beta <- betadisper(Q1_d, sample_data(FINAL_abundances_Q1_sqrt)$generation)
plot(beta)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # p<0.05 - NOT OK
TukeyHSD (beta)


## post hoc tests ##
testing_sampleType = pairwise.perm.manova(Q1_d, sample_data(FINAL_abundances_Q1_sqrt)$sampleType,
                                          nperm=10000, p.method = "BH")
testing_sampleType$p.value

testing_generation = pairwise.perm.manova(Q1_d, sample_data(FINAL_abundances_Q1_sqrt)$generation,
                                         nperm=10000, p.method = "BH")
testing_generation$p.value

sample_data(FINAL_abundances_Q1_sqrt)$SampleGen <-interaction(sample_data(FINAL_abundances_Q1_sqrt)$sampleType, sample_data(FINAL_abundances_Q1_sqrt)$generation)
testing_sampleType = pairwise.perm.manova(Q1_d, sample_data(FINAL_abundances_Q1_sqrt)$SampleGen,
                                          nperm=10000, p.method = "BH")
testing_sampleType$p.value


### larval age ###
#larvaeTimePoint is larval age. first = 1d; second = >1d
larvae_Q1 <- subset_samples (FINAL_abundances_Q1_sqrt, sampleType == "larvae")
larvae_Q1 = prune_taxa(taxa_sums(larvae_Q1) > 0, larvae_Q1) 
Q1_adonis_larvae = as (sample_data(larvae_Q1), "data.frame")
Q1_d_larvae = phyloseq::distance(larvae_Q1,'bray')
Adonis_Q1 <-adonis2(Q1_d_larvae ~ tank + larvaeTimePoint, data=Q1_adonis_larvae,  permutations = 10000, method = "bray")
Adonis_Q1 



############################################################
### DESEQ2 - differentially abundant ASVs ###

phyloseq_merged_FINAL_filtered # I need to provide count data for DESeq (no proportions)
phyloseq_merged_FINAL_deseq_Q1 <- subset_samples(phyloseq_merged_FINAL_filtered, climateFINAL == "C") #only ambient
phyloseq_merged_FINAL_deseq_Q1_noZeroZero <- prune_taxa(taxa_sums(phyloseq_merged_FINAL_deseq_Q1) > 0, phyloseq_merged_FINAL_deseq_Q1)
sample_data(phyloseq_merged_FINAL_deseq_Q1_noZeroZero)$sampleType <- as.factor(sample_data(phyloseq_merged_FINAL_deseq_Q1_noZeroZero)$sampleType) #converting sampleType into factor


## larvae F1 1d vs larvae F1 5d ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q1_noZeroZero, sampleType == "larvae" & generation == "F1")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~larvaeTimePoint) 
Q1_dds$larvaeTimePoint<-relevel(Q1_dds$larvaeTimePoint,ref="first") #reference
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("larvaeTimePoint","second", "first"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q1)
res_Q1_order <- res_Q1[order(res_Q1$padj),] #order according to padj 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ] #only significant
res_Q1_order_sig 
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix")) #create dataframe with taxonomy
res_Q1_sign_OTU <-tibble::rownames_to_column(as.data.frame(res_Q1_table_sign), var="OTU")
Abund <- subset_samples(FINAL_abundances_Q1_cutoff, sampleType == "larvae" & generation == "F1") 
Abund_cutoff <- filter_taxa(Abund, function(x) mean(x) > 0.01, TRUE) #only ASVs >1% mean relative abundance
Abund_cutoff_melted <- psmelt (Abund_cutoff)
larvaeF1_1dVS5d <- inner_join(res_Q1_sign_OTU, Abund_cutoff_melted, by = "OTU", copy = FALSE) 
larvaeF1_1dVS5d <- larvaeF1_1dVS5d[!duplicated(larvaeF1_1dVS5d$OTU),] 
larvaeF1_1dVS5d2 <- larvaeF1_1dVS5d %>% unite("OTU_fam", Class.x, Family.x,  Genus.x, Species.x, OTU, sep=",", remove=FALSE, na.rm = TRUE) #ASV list
#bubble plot
filtering <- subset_samples(FINAL_abundances_Q1_cutoff, sampleType == "larvae" & generation == "F1")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("d9bbccd9b043e43055cb13ca37c3680b", "98d4a54c43f07e8a85e574879dfb6714", "29279cf5fedb344b6a00a63947682318", "444e63ef145a3152188b3420f0343ce2", 
                                                                                         "33a5560897d81b1854da2a7b636c3ab7", "a65862dc95ce35a2f6a841e27cf89724", "24846d1a4a92a8cdb6bc3ed2c24d9307", "7f0ce4f17ba4f5a4ecccf1816c6f599f"))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_generation, scales = "free_x", space = "free_x")
bubbleplot


## larvae F1 (1d+15d) vs juvenile F1 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q1_noZeroZero, sampleType == "larvae" & generation == "F1" | sampleType == "juvenile")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~sampleType)
Q1_dds$sampleType<-relevel(Q1_dds$sampleType,ref="larvae") 
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("sampleType","juvenile", "larvae"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
res_Q1_order <- res_Q1[order(res_Q1$padj),] 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ] 
res_Q1_order_sig
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix")) 
res_Q1_sign_OTU <-tibble::rownames_to_column(as.data.frame(res_Q1_table_sign), var="OTU") 
Abund <- subset_samples(FINAL_abundances_Q1_cutoff, sampleType == "larvae" & generation == "F1" | sampleType == "juvenile") 
Abund_cutoff <- filter_taxa(Abund, function(x) mean(x) > 0.01, TRUE) #only mean rel ab >1%
Abund_cutoff_melted <- psmelt (Abund_cutoff)
larvaeF1vsJuvenile <- inner_join(res_Q1_sign_OTU, Abund_cutoff_melted, by = "OTU", copy = FALSE)
larvaeF1vsJuvenile <- larvaeF1vsJuvenile[!duplicated(larvaeF1vsJuvenile$OTU),]
larvaeF1vsJuvenile2 <- larvaeF1vsJuvenile %>% unite("OTU_fam", Class.x, Family.x, Genus.x, Species.x, OTU, sep=",", remove=FALSE, na.rm = TRUE)  #ASV list
#bubbleplot
filtering <- subset_samples(FINAL_abundances_Q1_cutoff, sampleType == "larvae" & generation == "F1" | sampleType == "juvenile")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("48a6a0e363dfb3823afdb092a0bb834d", "cd3c68e64285e960b6c46d8bf4cd75e6", "901bdd384be1d151e57df47fc0f66110", "12c4cf4a9b6a4541ef111f02191bae96", 
                                                                              "4ee16f924a58f3d50a57e9b5c692f1c8", "3b5e5025d3eb814d0abf098474ba41c7", "d9bbccd9b043e43055cb13ca37c3680b", "7f0ce4f17ba4f5a4ecccf1816c6f599f",
                                                                              "0677d9f732493f832161376bda44672f"))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_generation, scales = "free_x", space = "free_x")
bubbleplot



############################################################
### ALPHA DIVERSITY ###

## pre-processing ##
richQ1 = estimate_richness(rarefied_Q1_NOcutoff, split=TRUE) #I have to process data to run richness stats as 'richQ1' dataset doesn't include metadata
richQ1<-as.data.frame(richQ1)
RICHNESSVER<-tibble::rownames_to_column(as.data.frame(richQ1), var="sampleID")
RICHNESSVER$sampleID = gsub("X", "", RICHNESSVER$sampleID)
metadata10 <- tibble::rownames_to_column(as.data.frame(metadata9), var="sampleID")
metadata11 <- metadata10 %>% unite ("sample_gen", sampleType, generation, sep="_", remove=FALSE, na.rm = TRUE) 
metadata12 <- metadata11 %>% mutate(sample_gen = replace(sample_gen, sample_gen== "seawater_F0", "seawater")) 
metadata13 <- metadata12 %>% mutate(sample_gen = replace(sample_gen, sample_gen== "seawater_F1", "seawater"))
richness_table_ver <- inner_join(RICHNESSVER, metadata13, by="sampleID")

##calculate evenness
H<- richQ1$Shannon
S1<- richQ1$Observed
s<-log(S1)
evenness <- H/s
richness_table_ver$Evenness = evenness	

##calculate Faith's phylogenetic diversity			    
PD1 <- estimate_pd(rarefied_Q1_NOcutoff)
richness_table_ver$PhyloDiv = PD1$PD	
			    
#more preprocessing
table(sample_data(richness_table_ver)$sampleType)
table(sample_data(richness_table_ver)$generation)
richness_table_ver$sampleType <- as.factor(richness_table_ver$sampleType)
richness_table_ver$generation <- as.factor(richness_table_ver$generation)
richness_table_ver %>% dplyr::group_by(sampleType, generation)%>% count()
richness_table_ver %>% pull(sampleType) %>% levels()
richness_table_ver %>% pull(generation) %>% levels()
richness_table_ver2 <- richness_table_ver %>% mutate(Group =  factor(paste(sampleType,generation))) #better to have group instead of interaction between factors due to unbalanced design 
richness_table_ver2 %>% pull(Group) %>% levels()			    

## stats Shannon ##
mod1 <- glmmTMB(Shannon ~ Group + (1|tank), data = richness_table_ver2, family = 'gaussian', REML = TRUE)
summary(mod1)
simulateResiduals(fittedModel = mod1, plot = T) #no problems detected

## post hoc tests ##
emmeans (mod1, ~Group, type = "response") %>% pairs %>%rbind(adjust='bh')
			    
#plot Shannon is in next section
			    
#plot Evenness
p <- ggplot(for_plot, aes(x=sample_GEN, y=Evenness)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "Evenness") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)
p
			    
#plot tot ASVs
p <- ggplot(for_plot, aes(x=sample_GEN, y=Observed)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "Observed") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)
p
			    
			    
#plot Simpson
p <- ggplot(for_plot, aes(x=sample_GEN, y=Simpson)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "Simpson") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)
p

			    
#plot Faith's Phylogenetic Dviersity
p <- ggplot(for_plot, aes(x=sample_GEN, y=PhyloDiv)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "PhyloDiv") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)
p

			    
			    
############################################################
### VENN DIAGRAM ###

venn <- subset_samples(rarefied_Q1_NOcutoff, sample_generation =="gonad_F0" | sample_generation =="larvae_F1_first"| sample_generation =="larvae_F1_second" | sample_generation =="juvenile_F1")			    
ps_venn(venn,
          "sample_generation",
          fraction = 0.5, #The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count
          weight = FALSE, #If TRUE, the overlaps are weighted by abundance
          relative = TRUE, #Should abundances be made relative
          plot = TRUE, fill = c("#D6CCCC","#767C80","#D3DEDC", "#92A9BD") #figure	
          )
						
ps_venn(venn,
          "sample_generation",
          fraction = 0.5, 
          weight = FALSE, 
          relative = TRUE, 
          plot = FALSE #list shared ASVs
          )

venn <- subset_samples(rarefied_Q1_NOcutoff, sample_generation =="gonad_F0" | sample_generation =="gonad_F1" | sample_generation =="larvae_F1_first" | sample_generation =="larvae_F2_first")
ps_venn(venn,
        "sample_generation",
        fraction = 0.5, 
        weight = FALSE,
        relative = TRUE, 
        plot = TRUE, fill = c("#D6CCCC","#B5A4A4","#D3DEDC", "#9FB7B3")
)
			    
venn <- subset_samples(rarefied_Q1_NOcutoff, sample_generation =="larvae_F1_first" | sample_generation =="larvae_F1_second" |seawaterTimePoint =="larvaeF1")
ps_venn(venn,
        "sample_generation",
        fraction = 0.5, 
        weight = FALSE, 
        relative = TRUE,
        plot = TRUE, fill = c("#D3DEDC","#92A9BD","#4CB0DF")
)
			    
venn <- subset_samples(rarefied_Q1_NOcutoff, sample_generation =="juvenile_F1" | seawaterTimePoint =="juvenileF1")
ps_venn(venn,
        "sample_generation",
        fraction = 0.5,
        weight = FALSE,
        relative = TRUE,
        plot = TRUE, fill = c("#D6C298","#4CB0DF")
)			    

########################################################################################################################
########################################################################################################################







####################################################################################################################
################## Effects of OA and OW predicted for years 2050 and 2100 on the urchin microbiome ##################
####################################################################################################################

############################################################
### FILTERING ###
#climate1 = climate treatment (C=ambient, M= 2050, H=2100)

FINAL_abundances_Q2_cutoff <- subset_samples(FINALabundances_cutoff, treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant
FINAL_abundances_Q2_cutoff = prune_taxa(taxa_sums(FINAL_abundances_Q2_cutoff) > 0, FINAL_abundances_Q2_cutoff) #remove unobserved ASVs (sum 0 across all samples)

rarefied_Q2_NOcutoff <- subset_samples(FINAL_rarefied, treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant
rarefied_Q2_NOcutoff = prune_taxa(taxa_sums(rarefied_Q2_NOcutoff) > 0, rarefied_Q2_NOcutoff) #remove unobserved ASVs (sum 0 across all samples)


############################################################
### SQRT TRANSFORMATION ###

FINAL_abundances_Q2_sqrt <- transform_sample_counts(FINAL_abundances_Q2_cutoff, function (x) sqrt(x))


############################################################
### NMDS ###

Colors <- c(
  "#7A7C84", "#B05320", "#207DB0")
ordinate(FINAL_abundances_Q2_sqrt, "NMDS", "bray") %>% 
  plot_ordination(FINAL_abundances_Q2_sqrt, ., color = "climate1", shape = "sample_generation", title = "nmds_brays_Q2") + theme_bw() + geom_point(size = 4.5) + scale_color_manual(values = Colors) + scale_shape_manual(values = 0:8) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  theme(legend.text = element_text(size = 15),legend.title = element_text(size = 16)) + guides(col=guide_legend("Climate treatment"), shape = guide_legend("Sample Type and generation")) #change title 

###stress nmds
nmds_Q2_stress <- ordinate(FINAL_abundances_Q2_sqrt, "NMDS", "bray") 
cat("Stress:", nmds_Q2_stress$stress, fill=TRUE) 


############################################################
### BUBBLE PLOT DOMINANT FAMILIES ###
Q2_Fam <- tax_glom(FINAL_abundances_Q2_cutoff, taxrank = 'Family') #merging taxa of the same Family
q2_fam<- prune_taxa(taxa_sums(Q2_Fam) > 0, Q2_Fam) 
Q2_Gfam_sample <- merge_samples(q2_fam, "sample_treatment", fun=mean)  
Q2_Gfam_sample_abund = transform_sample_counts(Q2_Gfam_sample, function(x){x / sum(x)}) 
Q2_Gfam_sample_abund_melt <- psmelt(Q2_Gfam_sample_abund)
Q2_Gfam_sample_abund_melt$Family <- as.character(Q2_Gfam_sample_abund_melt$Family)
#dominant microbial families (>5%)
Q2_Fam5perc <- subset (Q2_Gfam_sample_abund_melt, Family == "D_4__Alcanivoracaceae" | Family == "D_4__Alteromonadaceae" | Family == "D_4__Cellvibrionaceae" | Family == "D_4__Cryomorphaceae" | 
                       Family == "D_4__Cycloclasticaceae" | Family == "D_4__Desulfobulbaceae" | Family == "D_4__Endozoicomonadaceae" | Family == "D_4__Flavobacteriaceae" | Family == "D_4__Fusobacteriaceae" | 
                       Family == "D_4__Halomonadaceae" | Family == "D_4__Hyphomonadaceae" | Family == "D_4__Kangiellaceae" | Family == "D_4__Kiritimatiellaceae" | 
                       Family == "D_4__Methyloligellaceae" |  Family == "D_4__Nannocystaceae" |  Family == "D_4__Nitrincolaceae"  |  Family == "D_4__Nitrosopumilaceae" |  Family == "D_4__Oleiphilaceae" | Family == "D_4__P13-46" | 
                       Family == "D_4__P3OB-42" | Family == "D_4__Pirellulaceae" |  Family == "D_4__Porticoccaceae" | Family == "D_4__Prolixibacteraceae" | Family == "D_4__Pseudomonadaceae" | 
                       Family == "D_4__Rhizobiaceae" | Family == "D_4__Rhodobacteraceae" | Family == "D_4__Saccharospirillaceae" | Family == "D_4__Saprospiraceae" | Family == "D_4__Sphingomonadaceae" | 
                       Family ==  "D_4__Spongiibacteraceae" | Family == "D_4__uncultured gamma proteobacterium"| Family == "D_4__Vibrionaceae")
p <- Q2_Fam5perc %>% mutate(name2 = fct_relevel(Sample, "gonad_C", "gonad_M", "gonad_H", "gonad_CCC", "gonad_HHH", "larvae_CC_first", "larvae_MM_first", "larvae_HH_first", "larvae_CCCC_first","larvae_HHHH_first", "larvae_CC_second", "larvae_MM_second", "larvae_HH_second", "larvae_CCCC_second", "larvae_HHHH_second","juvenile_CCC", "juvenile_MMM",  "juvenile_HHH", "seawater_C", "seawater_M", "seawater_H", "seawater_CCC", "seawater_HHH")) #to reoder x axis
bubbleFamily_Q2 <- ggplot(p, aes(x = name2, y = Family, color = Class)) + 
						     geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
						     scale_size_area(max_size = 20) + 
						     theme_bw() +
						     theme(axis.title.y=element_blank(),  axis.text.y=element_text(size=20), axis.ticks.y=element_blank(),
						           legend.text=element_text(size=11), legend.title=element_text(size=13), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=30)) +labs(title = "Q2_bubbleplot_fam")
						   bubbleFamily_Q2 +
						     scale_y_discrete(limits = c("D_4__Fusobacteriaceae", "D_4__Kiritimatiellaceae", "D_4__Nitrosopumilaceae", "D_4__Pirellulaceae",
						                                 "D_4__Desulfobulbaceae", "D_4__Nannocystaceae","D_4__P3OB-42",
						                                 "D_4__Cryomorphaceae", "D_4__Flavobacteriaceae", "D_4__Prolixibacteraceae", "D_4__Saprospiraceae",
						                                 "D_4__Hyphomonadaceae", "D_4__Methyloligellaceae", "D_4__Rhizobiaceae", "D_4__Rhodobacteraceae", "D_4__Sphingomonadaceae",
						                                 "D_4__Alcanivoracaceae", "D_4__Alteromonadaceae", "D_4__Cellvibrionaceae", "D_4__Endozoicomonadaceae",
						                                 "D_4__Halomonadaceae", "D_4__Kangiellaceae", "D_4__Nitrincolaceae", "D_4__Oleiphilaceae", 
						                                 "D_4__Cycloclasticaceae", "D_4__P13-46", "D_4__Porticoccaceae", "D_4__Pseudomonadaceae", "D_4__Saccharospirillaceae",
						                                 "D_4__Spongiibacteraceae", "D_4__uncultured gamma proteobacterium",   "D_4__Vibrionaceae")) +
						     scale_color_manual(values = c("D_2__Alphaproteobacteria" = "#A8DDC1", "D_2__Bacteroidia" = "#C3E7BC",
						                                   "D_2__Deltaproteobacteria" = "#DFDEC5", "D_2__Fusobacteriia" = "#CFCDCE",
						                                   "D_2__Gammaproteobacteria" = "#C0BCD7", "D_2__Kiritimatiellae" = "#D89FAB",
						                                   "D_2__Nitrososphaeria" = "#E58782", "D_2__Planctomycetacia" = "#CFB289"
						                                  )) +
						     theme( 
						       plot.background = element_blank(),
						       panel.grid.major = element_blank(),#to eliminate gridlines
						       panel.grid.minor = element_blank()) 




############################################################
### STATS BETA DIVERSITY ###

### sampleType*climate1 ; F0 adults, F1 larvae, F1 juvenile only ###
						    
FINAL_abundances_Q2_sqrt_part1 <- subset_samples(FINAL_abundances_Q2_sqrt, sample_generation != "gonad_F1" & sample_generation != "larvae_F2_first" & sample_generation != "larvae_F2_second" & sample_generation != "seawater")
table(sample_data(FINAL_abundances_Q2_sqrt_part1)$sample_generation)
Q2_adonis_part1 = as (sample_data(FINAL_abundances_Q2_sqrt_part1), "data.frame")
Q2_d_part1 = phyloseq::distance(FINAL_abundances_Q2_sqrt_part1,'bray') 						    
Adonis_Q2_part1 <-adonis2(Q2_d_part1 ~ sampleType*climate1 + tank, data=Q2_adonis_part1,  permutations = 10000, method = "bray") 
Adonis_Q2_part1
						    
## dispersion ##
beta <- betadisper(Q2_d_part1, sample_data(FINAL_abundances_Q2_sqrt_part1)$climate1)
disper.test = permutest(beta, permutations =10000)
disper.test #OK
						    
beta <- betadisper(Q2_d_part1, sample_data(FINAL_abundances_Q2_sqrt_part1)$sampleType)
plot(beta)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # p < 0.05, NOT OK
TukeyHSD (beta)
						    
beta <- betadisper(Q2_d_part1, sample_data(FINAL_abundances_Q2_sqrt_part1)$tank)
disper.test = permutest(beta, permutations =10000)
disper.test # p < 0.05, NOT OK

## post hoc tests ##
sample_data(FINAL_abundances_Q2_sqrt_part1)$climateSample <-interaction(sample_data(FINAL_abundances_Q2_sqrt_part1)$climate1, sample_data(FINAL_abundances_Q2_sqrt_part1)$sampleType)
testing_climateSample = pairwise.perm.manova(Q2_d_part1, sample_data(FINAL_abundances_Q2_sqrt_part1)$climateSample,
                                             nperm=10000, p.method = "BH")
testing_climateSample


### climate1 ; F1 adults only ###
						    
FINAL_abundances_Q2_sqrt_part2 <- subset_samples(FINAL_abundances_Q2_sqrt, sample_generation == "gonad_F1")
table(sample_data(FINAL_abundances_Q2_sqrt_part2)$sample_generation)
Q2_adonis_part2 = as (sample_data(FINAL_abundances_Q2_sqrt_part2), "data.frame")
Q2_d_part2 = phyloseq::distance(FINAL_abundances_Q2_sqrt_part2,'bray') 						    
Adonis_Q2_part2 <-adonis2(Q2_d_part2 ~ climate1 + tank, data=Q2_adonis_part2,  permutations = 10000, method = "bray") 
Adonis_Q2_part2 #significant
			
beta <- betadisper(Q2_d_part2, sample_data(FINAL_abundances_Q2_sqrt_part2)$climate1)
disper.test = permutest(beta, permutations =10000)
disper.test # OK
						    
						    
### climate1 ; F1 larvae only ###
						    
#1-day						    
FINAL_abundances_Q2_sqrt_part3 <-  subset_samples(FINAL_abundances_Q2_sqrt, sample_generation == "larvae_F1_first")
table(sample_data(FINAL_abundances_Q2_sqrt_part3)$sample_generation)
table(sample_data(FINAL_abundances_Q2_sqrt_part3)$climate1)	
Q2_adonis_part3 = as (sample_data(FINAL_abundances_Q2_sqrt_part3), "data.frame")
Q2_d_part3 = phyloseq::distance(FINAL_abundances_Q2_sqrt_part3,'bray') 
Adonis_Q2_part3 <-adonis2(Q2_d_part3 ~ climate1, data=Q2_adonis_part3,  permutations = 10000, method = "bray") 
Adonis_Q2_part3	
						    
beta <- betadisper(Q2_d_part3, sample_data(FINAL_abundances_Q2_sqrt_part3)$climate1)
disper.test = permutest(beta, permutations =10000)
disper.test # OK
						    
#5-day						    
FINAL_abundances_Q2_sqrt_part3 <-  subset_samples(FINAL_abundances_Q2_sqrt, sample_generation == "larvae_F1_second")
table(sample_data(FINAL_abundances_Q2_sqrt_part3)$sample_generation)
table(sample_data(FINAL_abundances_Q2_sqrt_part3)$climate1)
Q2_adonis_part3 = as (sample_data(FINAL_abundances_Q2_sqrt_part3), "data.frame")
Q2_d_part3 = phyloseq::distance(FINAL_abundances_Q2_sqrt_part3,'bray') 
Adonis_Q2_part3 <-adonis2(Q2_d_part3 ~ climate1, data=Q2_adonis_part3,  permutations = 10000, method = "bray") 
Adonis_Q2_part3	

beta <- betadisper(Q2_d_part3, sample_data(FINAL_abundances_Q2_sqrt_part3)$climate1)
disper.test = permutest(beta, permutations =10000)
disper.test # OK
						    
						    
### seawaterTimePoint*climate1  ###
FINAL_abundances_Q2_sqrt_part4 <-  subset_samples(FINAL_abundances_Q2_sqrt, sample_generation == "seawater")
table(sample_data(FINAL_abundances_Q2_sqrt_part4)$seawaterTimePoint)
table(sample_data(FINAL_abundances_Q2_sqrt_part4)$climate1)
Q2_adonis_part4 = as (sample_data(FINAL_abundances_Q2_sqrt_part4), "data.frame")
Q2_d_part4 = phyloseq::distance(FINAL_abundances_Q2_sqrt_part4,'bray') 
Adonis_Q2_part4 <-adonis2(Q2_d_part4 ~ seawaterTimePoint*climate1 + tank, data=Q2_adonis_part4,  permutations = 10000, method = "bray") 
Adonis_Q2_part4						    

beta <- betadisper(Q2_d_part4, sample_data(FINAL_abundances_Q2_sqrt_part4)$climate1)
disper.test = permutest(beta, permutations =10000)
disper.test # OK

beta <- betadisper(Q2_d_part4, sample_data(FINAL_abundances_Q2_sqrt_part4)$seawaterTimePoint)
disper.test = permutest(beta, permutations =10000)
disper.test # p < 0.05, NOT OK
						    
						    
############################################################
### DESEQ2 - differentially abundant ASVs ###

phyloseq_merged_FINAL_filtered # I need to provide count data for DESeq (no proportions)
phyloseq_merged_FINAL_deseq_Q2 <- subset_samples(phyloseq_merged_FINAL_filtered, treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant samples
phyloseq_merged_FINAL_deseq_Q2_noZeroZero <- prune_taxa(taxa_sums(phyloseq_merged_FINAL_deseq_Q2) > 0, phyloseq_merged_FINAL_deseq_Q2)
sample_data(phyloseq_merged_FINAL_deseq_Q2_noZeroZero)$climate1 <- as.factor(sample_data(phyloseq_merged_FINAL_deseq_Q2_noZeroZero)$climate1) #converting climate1 into factor

## JUVENILE AMB vs 2050 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sampleType == "juvenile" & climate1 !="H")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1)
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C")
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","M", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #22 ASVs -> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("9775dd8bbc76e65781081c268471ca3a", "e51567da19b784d647b39e5ba717e54a","408c648b8883625ecc3e098ff5f53698")) 

## JUVENILE AMB vs 2100 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sampleType == "juvenile" & climate1 !="M")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1) 
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C") 
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","H", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #22 ASVs -> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("ce4efbce621af5f2670efae86439c6ba")) #differentially abundant ASVs


## JUVENILE BUBBLE PLOT ##		
filtering <- subset_samples(FINAL_abundances_Q2_cutoff, sampleType == "juvenile")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("9775dd8bbc76e65781081c268471ca3a", "e51567da19b784d647b39e5ba717e54a","408c648b8883625ecc3e098ff5f53698",
                                                                              "ce4efbce621af5f2670efae86439c6ba")) #differentially abundant ASVs

variable1 = as.character(get_variable(physeq.subset, "sample"))
variable2 = as.character(get_variable(physeq.subset, "treatment"))
sample_data(physeq.subset)$rep_treat <- mapply(paste0, variable2, variable1,
                                               collapse = "_") #create new vairable in phyloseq
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Family,  Genus, OTU, sep=";", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
physeq.subset_m_ord <- physeq.subset_m %>% mutate(OTU_fam_2 = fct_relevel(OTU_fam, "D_4__Haliangiaceae;D_5__Haliangium;ce4efbce621af5f2670efae86439c6ba", "D_4__Flavobacteriaceae;D_5__Muricauda;9775dd8bbc76e65781081c268471ca3a", "D_4__Flavobacteriaceae;D_5__Lutibacter;408c648b8883625ecc3e098ff5f53698", "D_4__uncultured bacterium;D_5__;e51567da19b784d647b39e5ba717e54a"))
physeq.subset_m_ord$sample_treatment_g = factor (physeq.subset_m_ord$sample_treatment, levels=c('juvenile_CCC', 'juvenile_MMM', 'juvenile_HHH'))
bubbleplot <- ggplot(physeq.subset_m_ord, aes(x = rep_treat, y = OTU_fam_2, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,25), name="Relative abundance (%)") +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_treatment_g, scales = "free_x", space = "free_x") +
  scale_y_discrete(labels=c("D_4__Flavobacteriaceae;D_5__Muricauda;9775dd8bbc76e65781081c268471ca3a"="Flavobacteriaceae; Muricauda; 9775dd8bbc76e65781081c268471ca3a","D_4__uncultured bacterium;D_5__;e51567da19b784d647b39e5ba717e54a"="uncultured bacterium; e51567da19b784d647b39e5ba717e54a", "D_4__Flavobacteriaceae;D_5__Lutibacter;408c648b8883625ecc3e098ff5f53698"="Flavobacteriaceae; Lutibacter; 408c648b8883625ecc3e098ff5f53698", "D_4__Haliangiaceae;D_5__Haliangium;ce4efbce621af5f2670efae86439c6ba"="Haliangiaceae; Haliangium; ce4efbce621af5f2670efae86439c6ba")) 
bubbleplot +
  scale_color_manual(values = c("D_2__Bacteroidia" = "#C3E7BC", "D_2__Gammaproteobacteria" = "#C0BCD7", "D_2__Deltaproteobacteria" = "#DFDEC5"))

###########

## ADULT GONAD F0 AMB vs 2050 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sample_generation == "gonad_F0" & climate1 !="H")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1)
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C") 
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts")
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","M", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #11 ASVs-> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("09b78f73acb048b791cb7552a85a8e17", "99216daf0d0a03df19329d040f366982")) #differentially abundant ASVs

## ADULT GONAD F0 AMB vs 2100 ##	
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sample_generation == "gonad_F0" & climate1 !="M")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1) 
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C") 
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","H", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) 
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #15 ASVs -> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("9fb7e54678258f18659cd687049f8a55", "8a948cfcbfa0c36c6eea18e93ac6a275", "84ee7c41da835709760bcdd971f7cd56",
                                                                              "09b78f73acb048b791cb7552a85a8e17", "808949fe66333b1ceb99e1f11195e786")) #differentially abundant ASVs
																			  
## ADULT GONAD F0 BUBBLE PLOT ##	
filtering <- subset_samples(FINAL_abundances_Q2_cutoff, sample_generation == "gonad_F0")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("9fb7e54678258f18659cd687049f8a55", "8a948cfcbfa0c36c6eea18e93ac6a275", "84ee7c41da835709760bcdd971f7cd56",
                                                                              "09b78f73acb048b791cb7552a85a8e17", "808949fe66333b1ceb99e1f11195e786", "99216daf0d0a03df19329d040f366982")) 

variable1 = as.character(get_variable(physeq.subset, "sample"))
variable2 = as.character(get_variable(physeq.subset, "treatment"))
sample_data(physeq.subset)$rep_treat <- mapply(paste0, variable2, variable1,
                                               collapse = "_") #create new vairable in phyloseq
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Family,  Genus, OTU, sep=";", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
physeq.subset_m_ord <- physeq.subset_m %>% mutate(OTU_fam_2 = fct_relevel(OTU_fam, 
                                                                          "D_4__Psychromonadaceae;D_5__Psychromonas;9fb7e54678258f18659cd687049f8a55",
                                                                          "D_4__Prolixibacteraceae;D_5__Roseimarinus;8a948cfcbfa0c36c6eea18e93ac6a275",
                                                                          "D_4__Desulfobulbaceae;D_5__Desulfocapsa;09b78f73acb048b791cb7552a85a8e17",
                                                                          "D_4__Desulfobulbaceae;D_5__Desulfocapsa;84ee7c41da835709760bcdd971f7cd56",
                                                                          "D_4__Desulfobulbaceae;808949fe66333b1ceb99e1f11195e786",
                                                                          "D_4__Xenococcaceae;D_5__Pleurocapsa PCC-7319;99216daf0d0a03df19329d040f366982"))
physeq.subset_m_ord$sample_treatment_g = factor (physeq.subset_m_ord$sample_treatment, levels=c('gonad_C', 'gonad_M', 'gonad_H'))
bubbleplot <- ggplot(physeq.subset_m_ord, aes(x = rep_treat, y = OTU_fam_2, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,25), name="Relative abundance (%)" ) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_treatment_g, scales = "free_x", space = "free_x")
bubbleplot +
  scale_color_manual(values = c("D_2__Deltaproteobacteria" = "#DFDEC5",  "D_2__Bacteroidia" = "#C3E7BC","D_2__Gammaproteobacteria" = "#C0BCD7", "D_2__Oxyphotobacteria" = "#89ABC9")) 

###########

## LARVAE F1 1-day AMBIENT VS 2050 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sample_generation == "larvae_F1_first" & climate1 !="H")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
table(sample_data(Filtered)$sample_generation)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1) 
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C") 
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts")
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","M", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #18 ASVs -> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("0429a8a999c3238e12bbfaa1714d385e", "78db7fa6b16b73b101fbb8b282b5330e", "892402b1fcd3e237cd6b53d333fc16d4",
                                                                              "68ccaba761bb2635fceb85f3f9556d97","cf2eea2a38138614b3121fc8cfae6b3a", "eaae4c740095f734738d1e0f9ff01663",
                                                                              "27f8b975f616fc83f4918afb2f83a40f", "10afda2baef44de4c584a6641de399b1", "c3c32c0479bf4d183832a2e3cf235481", 
                                                                              "33a5560897d81b1854da2a7b636c3ab7","14b2d26bb9ddf9abfd5744a502fcc6d7")) #differentially abundant ASVs

## LARVAE F1 1-DAY AMBIENTT VS 2100 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sample_generation == "larvae_F1_first" & climate1 !="M")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1) 
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C")
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01 
res_Q2 <- results(Q2_dds, contrast=c("climate1","H", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE)
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #13 ASVs -> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("0429a8a999c3238e12bbfaa1714d385e", "24846d1a4a92a8cdb6bc3ed2c24d9307", "cf2eea2a38138614b3121fc8cfae6b3a",
                                                                              "eaae4c740095f734738d1e0f9ff01663",  "27f8b975f616fc83f4918afb2f83a40f", "b5b9214245fe31349cb3541b6ea9e9ef",
                                                                              "f420b5340d48941cc0e22ffef62e3df1")) #differentially abundant ASVs

## LARVAE F1 1-DAY BUBBLE PLOT ##	
filtering1 <- subset_samples(FINAL_abundances_Q2_cutoff, generation == "F1")
filtering <- subset_samples(filtering1, sampleTypeAll == "larvae_first")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("0429a8a999c3238e12bbfaa1714d385e", "24846d1a4a92a8cdb6bc3ed2c24d9307", "cf2eea2a38138614b3121fc8cfae6b3a",
                                                                              "eaae4c740095f734738d1e0f9ff01663",  "27f8b975f616fc83f4918afb2f83a40f", "b5b9214245fe31349cb3541b6ea9e9ef",
                                                                              "f420b5340d48941cc0e22ffef62e3df1", "78db7fa6b16b73b101fbb8b282b5330e", "892402b1fcd3e237cd6b53d333fc16d4",
                                                                              "68ccaba761bb2635fceb85f3f9556d97", "10afda2baef44de4c584a6641de399b1", "c3c32c0479bf4d183832a2e3cf235481", 
                                                                              "33a5560897d81b1854da2a7b636c3ab7","14b2d26bb9ddf9abfd5744a502fcc6d7"))

variable1 = as.character(get_variable(physeq.subset, "sample"))
variable2 = as.character(get_variable(physeq.subset, "treatment"))
sample_data(physeq.subset)$rep_treat <- mapply(paste0, variable2, variable1,
                                               collapse = "_") #create new vairable in phyloseq
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Family,  Genus, OTU, sep=";", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
physeq.subset_m_ord <- physeq.subset_m %>% mutate(OTU_fam_2 = fct_relevel(OTU_fam, 
                                                                          "D_4__Sphingomonadaceae;D_5__Sphingomonas;10afda2baef44de4c584a6641de399b1",
                                                                          "D_4__Hyphomonadaceae;D_5__Hyphomonas;f420b5340d48941cc0e22ffef62e3df1",
                                                                          "D_4__Rhodobacteraceae;b5b9214245fe31349cb3541b6ea9e9ef",
                                                                          "D_4__Rhodobacteraceae;D_5__Shimia;c3c32c0479bf4d183832a2e3cf235481",
                                                                          "D_4__Rhodobacteraceae;14b2d26bb9ddf9abfd5744a502fcc6d7",
                                                                          "D_4__Rhodobacteraceae;33a5560897d81b1854da2a7b636c3ab7",
                                                                          "D_4__Pseudomonadaceae;D_5__Pseudomonas;78db7fa6b16b73b101fbb8b282b5330e",
                                                                          "D_4__Nitrincolaceae;D_5__Marinobacterium;68ccaba761bb2635fceb85f3f9556d97",
                                                                          "D_4__Pseudoalteromonadaceae;D_5__Pseudoalteromonas;24846d1a4a92a8cdb6bc3ed2c24d9307",
                                                                          "D_4__Oleiphilaceae;D_5__Oleiphilus;892402b1fcd3e237cd6b53d333fc16d4",
                                                                          "D_4__Vibrionaceae;D_5__Vibrio;0429a8a999c3238e12bbfaa1714d385e",
                                                                          "D_4__Moraxellaceae;D_5__Psychrobacter;cf2eea2a38138614b3121fc8cfae6b3a",
                                                                          "D_4__Moraxellaceae;D_5__Acinetobacter;eaae4c740095f734738d1e0f9ff01663",
                                                                          "D_4__Alteromonadaceae;D_5__Aestuariibacter;27f8b975f616fc83f4918afb2f83a40f"
                                                                          ))
physeq.subset_m_ord$sample_treatment_g = factor (physeq.subset_m_ord$sample_treatment, levels=c('larvae_CC_first', 'larvae_MM_first', 'larvae_HH_first'))
bubbleplot <- ggplot(physeq.subset_m_ord, aes(x = rep_treat, y = OTU_fam_2, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,25), name="Relative abundance (%)") +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_treatment_g, scales = "free_x", space = "free_x")
bubbleplot +
  scale_color_manual(values = c("D_2__Alphaproteobacteria" = "#A8DDC1", "D_2__Gammaproteobacteria" = "#C0BCD7")) 
  
###########

## LARVAE F1 5-day AMBIENT VS 2050 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sample_generation == "larvae_F1_second" & climate1 !="H")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
head(sample_data(Filtered)$climate1, 6) 
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1)
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C") 
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts")
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","M", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #15 ASVs -> among them, I keep only ASVs present in at least half of the samples in one treatment group
KEEP <- c("3df2134cb53c720c05a6a35a83ef6d15", "c595de1845ebd5bb497bed4ad0a56431", "bed4df5dd89852e1e36acc676cc57bd6",
                                                                              "d291a36801ef1e17761e727ffcb814a5", "68b9ce06fc98f29759c01b03263a5d67", "5cd1ce661f4a05cd6c3732999d7d6053",
                                                                              "e838116200a8436ae7fdd7004643cfe5", "c41aa9c9ce9467cf47962a420beaaa94", "c855ad414d3bbdde6e3eb4f68955924c",
                                                                              "9d10a746c450af4fab9ee831e0b5b4bb", "62032461f24ee5870f288214f17d5ab2", "b454cbe5493189c498f863f3712a6272",
                                                                              "2261bf5c940e3b1722702a281dfbbbb9", "73d8983eb5ded77a4133af2d535baf7c", "217b8f3a9468e8701af9871d73440b4d"))
																				  #differentially abundant ASVs

## LARVAE F1 5-day AMBIENT VS 2100 ##
Filtered <- subset_samples(phyloseq_merged_FINAL_deseq_Q2_noZeroZero, sample_generation == "larvae_F1_second" & climate1 !="M")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q2_dds = phyloseq_to_deseq2(Filtered, ~climate1) 
Q2_dds$climate1<-relevel(Q2_dds$climate1,ref="C") 
Q2_dds = DESeq(Q2_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q2 <- results(Q2_dds, contrast=c("climate1","H", "C"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q2) #summary of results
res_Q2_order <- res_Q2[order(res_Q2$padj),] #order according to padj 
res_Q2_order = res_Q2_order[order(res_Q2_order$padj, na.last=NA), ]
res_Q2_order_sig = res_Q2_order[(res_Q2_order$padj < alpha), ] #only significant
res_Q2_order_sig #15 ASVs -> among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <-  c("47a6c3d4b77e22e2363b6f567f6b6ff5", "c41aa9c9ce9467cf47962a420beaaa94", "e8093da980543c5406b1b2a8d14ae1ac",
                                                                              "5cdf326d1392b9053408d8b7e8135342", "0a94faf00b60b6c7b2f6f3d7e8f32680", "08a38227ea0a4fa20ac48786ceecedf4",
                                                                              "9d10a746c450af4fab9ee831e0b5b4bb", "7c96d8954560aee0155b635945a904ca", "fd8dd59bca314bac1716641f4520cd01",
                                                                              "bed4df5dd89852e1e36acc676cc57bd6", "217b8f3a9468e8701af9871d73440b4d", "e03c201ace81dec799a397b73f962d09", 
                                                                              "fe6d797d623713681a107859038d7010", "651083dd15565eead4fdf9cd1344a6e5", "aa01366a9e121f9962359babd78685ac"))

## LARVAE F1 5-DAY BUBBLE PLOT ##	
filtering1 <- subset_samples(FINAL_abundances_Q2_cutoff, generation == "F1")
filtering <- subset_samples(filtering1, sampleTypeAll == "larvae_second")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("c41aa9c9ce9467cf47962a420beaaa94", "9d10a746c450af4fab9ee831e0b5b4bb", "bed4df5dd89852e1e36acc676cc57bd6",
                                                                               "217b8f3a9468e8701af9871d73440b4d",
                                                                               "47a6c3d4b77e22e2363b6f567f6b6ff5", "e8093da980543c5406b1b2a8d14ae1ac",  "08a38227ea0a4fa20ac48786ceecedf4",
                                                                               "5cdf326d1392b9053408d8b7e8135342", "0a94faf00b60b6c7b2f6f3d7e8f32680", "aa01366a9e121f9962359babd78685ac",
                                                                               "7c96d8954560aee0155b635945a904ca", "fd8dd59bca314bac1716641f4520cd01", "651083dd15565eead4fdf9cd1344a6e5",
                                                                               "e03c201ace81dec799a397b73f962d09","fe6d797d623713681a107859038d7010",
                                                                               "3df2134cb53c720c05a6a35a83ef6d15", "c595de1845ebd5bb497bed4ad0a56431", "73d8983eb5ded77a4133af2d535baf7c",
                                                                               "d291a36801ef1e17761e727ffcb814a5", "68b9ce06fc98f29759c01b03263a5d67", "5cd1ce661f4a05cd6c3732999d7d6053",
                                                                               "e838116200a8436ae7fdd7004643cfe5", "c855ad414d3bbdde6e3eb4f68955924c", "2261bf5c940e3b1722702a281dfbbbb9",
                                                                               "62032461f24ee5870f288214f17d5ab2", "b454cbe5493189c498f863f3712a6272"))

variable1 = as.character(get_variable(physeq.subset, "sample"))
variable2 = as.character(get_variable(physeq.subset, "treatment"))
sample_data(physeq.subset)$rep_treat <- mapply(paste0, variable2, variable1,
                                               collapse = "_") #create new vairable in phyloseq
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Family,  Genus, OTU, sep=";", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
physeq.subset_m_ord <- physeq.subset_m %>% mutate(OTU_fam_2 = fct_relevel(OTU_fam, 
                                                                          "D_4__uncultured archaeon;D_5__;62032461f24ee5870f288214f17d5ab2",
                                                                          "D_4__Puniceicoccaceae;217b8f3a9468e8701af9871d73440b4d",
                                                                          "D_4__Nannocystaceae;D_5__uncultured bacterium;fe6d797d623713681a107859038d7010",
                                                                          "D_4__Bdellovibrionaceae;D_5__OM27 clade;73d8983eb5ded77a4133af2d535baf7c",
                                                                          "D_4__uncultured bacterium;D_5__;fd8dd59bca314bac1716641f4520cd01",
                                                                          "D_4__Bacteriovoracaceae;5cd1ce661f4a05cd6c3732999d7d6053",
                                                                          "D_4__Bacteriovoracaceae;D_5__Peredibacter;e838116200a8436ae7fdd7004643cfe5",
                                                                          "D_4__Terasakiellaceae;D_5__Aestuariispira;0a94faf00b60b6c7b2f6f3d7e8f32680",
                                                                          "D_4__Rhodobacteraceae;c595de1845ebd5bb497bed4ad0a56431",
                                                                          "D_4__Rhodobacteraceae;7c96d8954560aee0155b635945a904ca",
                                                                          "D_4__Flavobacteriaceae;651083dd15565eead4fdf9cd1344a6e5",
                                                                          "D_4__NS9 marine group;D_5__uncultured Bacteroidetes bacterium;e03c201ace81dec799a397b73f962d09",
                                                                          "D_4__Saprospiraceae;D_5__Lewinella;c855ad414d3bbdde6e3eb4f68955924c",
                                                                          "D_4__Flavobacteriaceae;D_5__Tenacibaculum;3df2134cb53c720c05a6a35a83ef6d15",
                                                                          "D_4__Flavobacteriaceae;bed4df5dd89852e1e36acc676cc57bd6",
                                                                          "D_4__Cryomorphaceae;c41aa9c9ce9467cf47962a420beaaa94",
                                                                          "D_4__Cryomorphaceae;9d10a746c450af4fab9ee831e0b5b4bb",
                                                                          "D_4__Salinisphaeraceae;D_5__Salinisphaera;aa01366a9e121f9962359babd78685ac",
                                                                          "D_4__Vibrionaceae;D_5__Photobacterium;08a38227ea0a4fa20ac48786ceecedf4",
                                                                          "D_4__Vibrionaceae;D_5__Photobacterium;e8093da980543c5406b1b2a8d14ae1ac",
                                                                          "D_4__Vibrionaceae;D_5__Photobacterium;47a6c3d4b77e22e2363b6f567f6b6ff5",
                                                                          "D_4__Colwelliaceae;D_5__Thalassotalea;2261bf5c940e3b1722702a281dfbbbb9",
                                                                          "D_4__Kangiellaceae;D_5__Kangiella;b454cbe5493189c498f863f3712a6272",
                                                                          "D_4__Vibrionaceae;D_5__Vibrio;68b9ce06fc98f29759c01b03263a5d67",
                                                                          "D_4__Spongiibacteraceae;D_5__BD1-7 clade;d291a36801ef1e17761e727ffcb814a5",
                                                                          "D_4__Alteromonadaceae;D_5__Rheinheimera;5cdf326d1392b9053408d8b7e8135342"))
physeq.subset_m_ord$sample_treatment_g = factor (physeq.subset_m_ord$sample_treatment, levels=c('larvae_CC_second', 'larvae_MM_second', 'larvae_HH_second'))
bubbleplot <- ggplot(physeq.subset_m_ord, aes(x = rep_treat, y = OTU_fam_2, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,25), name="Relative abundance (%)") +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_treatment_g, scales = "free_x", space = "free_x")
bubbleplot +
  scale_color_manual(values = c("D_2__Alphaproteobacteria" = "#A8DDC1", "D_2__Gammaproteobacteria" = "#C0BCD7", "D_2__Bacteroidia" = "#C3E7BC", "D_2__Deltaproteobacteria" = "#DFDEC5", 
                                 "D_2__Verrucomicrobiae" = "#F3B962", "D_2__Thermoplasmata" = "#CECE66")) 




############################################################
### PCA ###
low_ab <- phyloseq::genefilter_sample(phyloseq_merged_FINAL_filtered, filterfun_sample(function(x) {x / sum(x)} > 1e-5))
for_mixomics <- phyloseq::prune_taxa(low_ab, phyloseq_merged_FINAL_filtered)
for_mixomics = prune_taxa(taxa_sums(for_mixomics) > 0, for_mixomics)
	
	
## adult gonad F0 ##	
for_mixomics_selection <- subset_samples(for_mixomics, sample_generation == "gonad_F0" & treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant
taxo <- tax_table(for_mixomics_selection) # extract the taxonomy
meta.data <- for_mixomics_selection@sam_data # extract the metadata
data.raw <- t(otu_table(for_mixomics_selection)) # extract OTU table from phyloseq object # samples should be in row and variables in column
data.offset <- data.raw+1 ## STEP 1: OFFSET
sum(which(data.offset == 0)) # double check there are no zeros
dim(data.offset) # check dimensions
pca.result <- pca(data.offset, logratio = 'CLR') # undergo PCA after CLR transformation
Colors <- c(
  "#7A7C84", "#B05320", "#207DB0")
plotIndiv(pca.result, 
          group = meta.data$climateFINAL, 
          title = 'gonad - no transplant, PCA Comps 1&2',
          legend=TRUE, ind.names = FALSE,
          star=TRUE, col.per.group = Colors) 
	
## larvae F1 1-day ##
for_mixomics_selection <- subset_samples(for_mixomics, sample_generation == "larvae_F1_first" & treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant
taxo <- tax_table(for_mixomics_selection) 
meta.data <- for_mixomics_selection@sam_data 
data.raw <- t(otu_table(for_mixomics_selection)) 
data.offset <- data.raw+1 
sum(which(data.offset == 0)) 
dim(data.offset)
pca.result <- pca(data.offset, logratio = 'CLR') 
Colors <- c(
  "#7A7C84", "#B05320", "#207DB0")
plotIndiv(pca.result, 
          group = meta.data$climateFINAL, 
          title = 'larvae F1 1d - no transplant, PCA Comps 1&2',
          legend=TRUE, ind.names = FALSE,
          star=TRUE, ellipse = FALSE, col.per.group = Colors) 
	

## larvae F1 5-day ##	
for_mixomics_selection <- subset_samples(for_mixomics, sample_generation == "larvae_F1_second" & treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant
taxo <- tax_table(for_mixomics_selection) 
meta.data <- for_mixomics_selection@sam_data 
data.raw <- t(otu_table(for_mixomics_selection)) 
data.offset <- data.raw+1 
sum(which(data.offset == 0)) 
dim(data.offset) 
pca.result <- pca(data.offset, logratio = 'CLR') 
Colors <- c(
  "#7A7C84", "#B05320", "#207DB0")
plotIndiv(pca.result, 
          group = meta.data$climateFINAL, 
          title = 'larvae F1 5d - no transplant, PCA Comps 1&2',
          legend=TRUE, ind.names = FALSE,
          star=TRUE, ellipse = FALSE, col.per.group = Colors) 

						
## juvenile F1 ##							 
for_mixomics_selection <- subset_samples(for_mixomics, sampleType == "juvenile" & treatment != "CCM" & treatment != "CCH" & treatment != "MMH" & treatment != "MMHH") #no transplant
taxo <- tax_table(for_mixomics_selection) 
meta.data <- for_mixomics_selection@sam_data 
data.raw <- t(otu_table(for_mixomics_selection)) 
data.offset <- data.raw+1 
sum(which(data.offset == 0)) 
dim(data.offset) 
pca.result <- pca(data.offset, logratio = 'CLR') 
Colors <- c(
  "#7A7C84", "#B05320", "#207DB0")
plotIndiv(pca.result, 
          group = meta.data$climateFINAL, 
          title = 'juvenile - no transplant, PCA Comps 1&2',
          legend=TRUE, ind.names = FALSE,
          star=TRUE, col.per.group = Colors) 
		  
	  				 


############################################################
### ALPHA DIVERSITY ###

## pre-processing ##
richQ2 = estimate_richness(rarefied_Q2_NOcutoff, split=TRUE) #to calculate diversity indexes: shannon, simpson, etc.
richQ2<-as.data.frame(richQ2)
RICHNESS_Q2<-tibble::rownames_to_column(as.data.frame(richQ2), var="sampleID") #I need to create a new column named 'sampleID' so I can use right_join in the next command (otherwise sampleiDs would be only on row_names and right_join would not work)
RICHNESS_Q2$sampleID = gsub("X", "", RICHNESS_Q2$sampleID) # I modify SampleIDs as there was an extra X at the beginning of each ID. I therefore substitute X in RICHNESS$SAMPLE with ""
metadata10 <- tibble::rownames_to_column(as.data.frame(metadata9), var="sampleID")
metadata11 <- metadata10 %>% unite ("sample_gen", sampleType, generation, sep="_", remove=FALSE, na.rm = TRUE) # I add a column as combination of sampleType and generation to the metadat10
metadata12 <- metadata11 %>% mutate(sample_gen = replace(sample_gen, sample_gen== "seawater_F0", "seawater")) #I want to have 'seawater' instead of 'seawater_F0'
metadata13 <- metadata12 %>% mutate(sample_gen = replace(sample_gen, sample_gen== "seawater_F1", "seawater")) #I want to have 'seawater' instead of 'seawater_F1'
richness_table_Q2 <- inner_join(RICHNESS_Q2, metadata13, by="sampleID")


## plots ##
#AMBIENT
alpha_plot_Q2_C <- subset(richness_table_Q2, climate1 =="C")
alpha_plot_Q2_C_order <- alpha_plot_Q2_C %>% mutate(name2 = fct_relevel(sample_gen, "gonad_F0", "gonad_F1", "larvae_F1", "larvae_F2", "juvenile_F1", "seawater")) #to reorder x axis
nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(4, "Set3"))(nb.cols)
p <- ggplot(alpha_plot_Q2_C_order, aes(x=name2, y=Shannon)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)+ ylim(2.1, 6.7) + labs(title = "Shannon C") + theme_classic(base_size = 15) #theme classic no background 
p


#2050
alpha_plot_Q2_M <- subset(richness_table_Q2, climate1 =="M")
alpha_plot_Q2_M_order <- alpha_plot_Q2_M %>% mutate(name2 = fct_relevel(sample_gen, "gonad_F0", "larvae_F1", "juvenile_F1", "seawater")) #to reorder x axis
nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(4, "Set3"))(nb.cols)
p <- ggplot(alpha_plot_Q2_M_order, aes(x=name2, y=Shannon)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "Shannon") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)+ ylim(2.1, 6.7) + labs(title = "Shannon M") + theme_classic(base_size = 15)
p

#2100	
alpha_plot_Q2_H <- subset(richness_table_Q2, climate1 =="H")
alpha_plot_Q2_H_order <- alpha_plot_Q2_H %>% mutate(name2 = fct_relevel(sample_gen, "gonad_F0", "gonad_F1", "larvae_F1", "larvae_F2", "juvenile_F1", "seawater")) #to reorder x axis
nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(4, "Set3"))(nb.cols)
p <- ggplot(alpha_plot_Q2_H_order, aes(x=name2, y=Shannon)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "Shannon") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=sampleType) +
  geom_point(colour = "black", size = 1)+ ylim(2.1, 6.7) + labs(title = "Shannon H") + theme_classic(base_size = 15)
p
	
								 
## more pre-processing ##
richness_table_Q2$sampleType <- as.factor(richness_table_Q2$sampleType)
richness_table_Q2$generation <- as.factor(richness_table_Q2$generation)
richness_table_Q2$climate1 <- as.factor(richness_table_Q2$climate1)
richness_table_Q2 %>% dplyr::group_by(climate1, generation, sampleType)%>% count()
richness_table_Q2 %>% pull(sampleType) %>% levels()
richness_table_Q2 %>% pull(generation) %>% levels()
richness_table_Q2 %>% pull(climate1) %>% levels()
richness_table_Q2_final <- richness_table_Q2 %>% mutate(Group =  factor(paste(sampleType,generation, climate1)))
richness_table_Q2_final %>% pull(Group) %>% levels()

## stats Shannon diversity
mod1 <- glmmTMB(Shannon ~ Group + (1|tank), data = richness_table_Q2_final, family = 'gaussian', REML = TRUE)
summary(mod1)
simulateResiduals(fittedModel = mod1, plot = T) #no problems detected
emmeans (mod1, ~Group, type = "response") %>% pairs %>%rbind(adjust='bh')


############################################################
### VENN DIAGRAM ###

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 != "M" & sample_generation =="larvae_F1_first" | climate1 != "M" & seawaterTimePoint =="larvaeF1")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5,
        weight = FALSE,
        relative = TRUE,
        plot = TRUE, fill = c("#D3DEDC","#D3DEDC", "#4CB0DF", "#4CB0DF")
)

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 != "H" & sample_generation =="larvae_F1_first" | climate1 != "H" & seawaterTimePoint =="larvaeF1")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5,
        weight = FALSE,
        relative = TRUE,
        plot = TRUE, fill = c("#D3DEDC","#D3DEDC", "#4CB0DF", "#4CB0DF")
)

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 != "M" & sample_generation =="larvae_F1_first" | climate1 != "M" & seawaterTimePoint =="larvaeF1")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5, #The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count
        weight = FALSE, #If TRUE, the overlaps are weighted by abundance
        relative = TRUE, #Should abundances be made relative
        plot = TRUE, fill = c("#D3DEDC","#D3DEDC", "#4CB0DF", "#4CB0DF")
)

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 != "H" & sample_generation =="larvae_F1_first" | climate1 != "H" & seawaterTimePoint =="larvaeF1")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5, #The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count
        weight = FALSE, #If TRUE, the overlaps are weighted by abundance
        relative = TRUE, #Should abundances be made relative
        plot = TRUE, fill = c("#D3DEDC","#D3DEDC", "#4CB0DF", "#4CB0DF")
)

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 == "M" & sample_generation =="gonad_F0" | climate1 == "M" & sample_generation =="larvae_F1_first"| climate1 == "M" & sample_generation =="larvae_F1_second" | climate1 == "M" & sample_generation =="juvenile_F1")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5, #The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count
        weight = FALSE, #If TRUE, the overlaps are weighted by abundance
        relative = TRUE, #Should abundances be made relative
        plot = TRUE
)

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 == "H" & sample_generation =="gonad_F0" | climate1 == "H" & sample_generation =="larvae_F1_first"| climate1 == "H" & sample_generation =="larvae_F1_second" | climate1 == "H" & sample_generation =="juvenile_F1")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5, #The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count
        weight = FALSE, #If TRUE, the overlaps are weighted by abundance
        relative = TRUE, #Should abundances be made relative
        plot = TRUE
)

venn <- subset_samples(rarefied_Q2_NOcutoff, climate1 == "H" & sample_generation =="gonad_F0" | climate1 == "H" & sample_generation =="gonad_F1" | climate1 == "H" & sample_generation =="larvae_F1_first" | climate1 == "H" & sample_generation =="larvae_F2_first")
ps_venn(venn,
        "sample_treatment",
        fraction = 0.5, #The fraction (0 to 1) of samples in a group in which the taxa should be present to be included in the count
        weight = FALSE, #If TRUE, the overlaps are weighted by abundance
        relative = TRUE, #Should abundances be made relative
        plot = TRUE
)


########################################################################################################################
########################################################################################################################








####################################################################################################################
########################## Parental effect on the microbiome responses to climate change ##########################
####################################################################################################################


############################################################
### FILTERING ###

FINAL_abundances_Q3_cutoff <- subset_samples(FINALabundances_cutoff, sampleType != "larvae" & sample_generation != "gonad_F1") 
FINAL_abundances_Q3_cutoff = prune_taxa(taxa_sums(FINAL_abundances_Q3_cutoff) > 0, FINAL_abundances_Q3_cutoff) 

rarefied_Q3_NOcutoff <- subset_samples(FINAL_rarefied, sampleType != "larvae" & sample_generation != "gonad_F1") 
rarefied_Q3_NOcutoff = prune_taxa(taxa_sums(rarefied_Q3_NOcutoff) > 0, rarefied_Q3_NOcutoff) 


############################################################
### SQRT TRANSFORMATION ###

FINAL_abundances_Q3_sqrt <- transform_sample_counts(FINAL_abundances_Q3_cutoff, function (x) sqrt(x))


############################################################
### multilevel sPLS-DA ###

for_mixomics_selection <- subset_samples(for_mixomics, sampleType == "juvenile" & treatment != "CCC") #excludiing ambient-ambient samples
table(sample_data(for_mixomics_selection)$sample_treatment)
taxo <- tax_table(for_mixomics_selection) # extraction of the taxonomy
meta.data <- for_mixomics_selection@sam_data # extraction of the metadata
data.raw <- t(otu_table(for_mixomics_selection)) # extract OTU table from phyloseq object # samples should be in row and variables in column
data.offset <- data.raw+1 ## STEP 1: OFFSET
sum(which(data.offset == 0)) # check how many zeroes there are
dim(data.offset) # check dimensions
pca.result <- pca(data.offset, logratio = 'CLR') # undergo PCA after CLR transformation
design <- data.frame(sample = meta.data$tank) #tank as repeated measure
X <- logratio.transfo(data.offset,logratio = "CLR") #data transformation
Y <- as.factor(meta.data$treatment)
q3.splsda <- splsda(X, Y, ncomp = 10)
perf.splsda.q3 <- perf(q3.splsda, validation = "Mfold", 
  folds = 3, nrepeat = 50,
  progressBar = FALSE, auc = TRUE) 
plot(perf.splsda.q3, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")
perf.splsda.q3$error.rate.class$mahalanobis.dist
perf.splsda.q3$choice.ncomp 
list.keepX <- c(1:10,  seq(20, 300, 10)) 
tune.splsda.q3 <- tune.splsda(X, Y, ncomp = 7,
   validation = 'Mfold',
  folds = 3, nrepeat = 50, # use repeated cross-validation
  dist = 'max.dist', 
 measure = "BER", # use  balanced error rate of dist measure instead of overall (suggested in mixomics book and forum for dataset with unbalanced design) 
 test.keepX = list.keepX, #number of variables to select on each component
  cpus = 2) # allow for paralleliation to decrease runtime
plot(tune.splsda.q3, col = color.jet(7)) # plot output of variable number tuning
#the optimal number of components accordding to our one-sided t-tests:
tune.splsda.q3$choice.ncomp$ncomp 
#the optimal keppX parameter according to minimal error rate:
tune.splsda.q3$choice.keepX 
optimal.ncomp <- tune.splsda.q3$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.q3$choice.keepX[1:optimal.ncomp]
final.multilevel.splsda.q3 <- splsda(X, Y, ncomp = optimal.ncomp, 
                                              keepX = optimal.keepX,
                                              multilevel = design)
Colors <- c("#7A7C84", "#B05320", "#207DB0")
plotIndiv(final.multilevel.splsda.q3, group = meta.data$climate1, 
          ind.names = meta.data$climateFINAL, #
          legend = TRUE, legend.title = 'Parental climate treatment',
          ellipse = FALSE, col.per.group = Colors,
          pch = as.factor(meta.data$climateFINAL), legend.title.pch = 'Offspring climate treatment',
          title = 'Sample Plot of multilevel sPLS-DA on urchin juvenile data')



############################################################
### BUBBLE PLOT DOMINANT FAMILIES JUVENILE F1 ###

juv_Q3 <- subset_samples(FINAL_abundances_Q3_cutoff, sampleType == "juvenile")
Q3_Fam_juv <- tax_glom(juv_Q3, taxrank = 'Family') 
q3_fam_juv <- prune_taxa(taxa_sums(Q3_Fam_juv) > 0, Q3_Fam_juv) 
Q3_Gfam_sample_juv <- merge_samples(q3_fam_juv, "sample_treatment", fun=mean)  
Q3_Gfam_sample_juv_abund = transform_sample_counts(Q3_Gfam_sample_juv, function(x){x / sum(x)}) 
Q3_Gfam_sample_juv_abund_melt <- psmelt(Q3_Gfam_sample_juv_abund)
Q3_Gfam_sample_juv_abund_melt$Family <- as.character(Q3_Gfam_sample_juv_abund_melt$Family) # convert Family to a character vector from a factor
Q3_Fam2perc <- subset (Q3_Gfam_sample_juv_abund_melt, Family == "D_4__Amoebophilaceae" | Family == "D_4__Bdellovibrionaceae" | Family == "D_4__Colwelliaceae" | Family == "D_4__Cryomorphaceae" | 
                         Family == "D_4__Cyanobacteriaceae" | Family == "D_4__Endozoicomonadaceae" | Family == "D_4__Flavobacteriaceae" | Family == "D_4__Halieaceae" | 
                         Family == "D_4__Marinilabiliaceae" | Family == "D_4__P13-46" | Family == "D_4__Pirellulaceae" | Family == "D_4__Prolixibacteraceae" | 
                         Family == "D_4__Rhizobiaceae" | Family == "D_4__Rhodobacteraceae" | Family == "D_4__Saccharospirillaceae" | Family == "D_4__Saprospiraceae" |
                         Family == "D_4__Vibrionaceae" | Family == "D_4__Xenococcaceae")
p <- Q3_Fam2perc %>% mutate(name2 = fct_relevel(Sample, "juvenile_CCC", "juvenile_CCM", "juvenile_CCH", "juvenile_MMM", "juvenile_MMH", "juvenile_HHH")) #to reoder x axis
bubbleFamily_Q3 <- ggplot(p, aes(x = name2, y = Family, color = Class)) + 
							   geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
							   scale_size_area(max_size = 20, name="Relative abundance (%)" ) + 
							   theme_bw() +
							   theme(axis.title.y=element_blank(),  axis.text.y=element_text(size=20), axis.ticks.y=element_blank(),
							         legend.text=element_text(size=11), legend.title=element_text(size=13), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=30)) +labs(title = "Q3_bubbleplot_fam") +
							   theme( 
							     plot.background = element_blank(),
							     panel.grid.major = element_blank(),#to eliminate gridlines
							     panel.grid.minor = element_blank())
bubbleFamily_Q3 + scale_y_discrete(limits = c( "D_4__Pirellulaceae",
							                               "D_4__Cyanobacteriaceae",   "D_4__Xenococcaceae",
							                               "D_4__Bdellovibrionaceae",
							                               "D_4__Amoebophilaceae", "D_4__Marinilabiliaceae",  "D_4__Prolixibacteraceae", "D_4__Cryomorphaceae", "D_4__Flavobacteriaceae", "D_4__Saprospiraceae",
							                               "D_4__P13-46",  "D_4__Saccharospirillaceae","D_4__Colwelliaceae","D_4__Halieaceae",  "D_4__Endozoicomonadaceae","D_4__Vibrionaceae",
							                               "D_4__Rhizobiaceae","D_4__Rhodobacteraceae")) +
							   scale_color_manual(values = c("D_2__Alphaproteobacteria" = "#A8DDC1", "D_2__Gammaproteobacteria" = "#C0BCD7", "D_2__Bacteroidia" = "#C3E7BC", "D_2__Deltaproteobacteria" = "#DFDEC5", 
							                                 "D_2__Oxyphotobacteria" = "#89ABC9",  "D_2__Planctomycetacia" = "#CFB289")) +
							   scale_x_discrete(labels = c("Ambient-Ambient", "Ambient-2050", "Ambient-2100", "2050-2050", "2050-2100", "2100-2100")) +
							   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
  


############################################################
## STATS BETA DIVERSITY JUVENILE F1	###

### climate1*sampleType ###
FINAL_abundances_Q3_sqrt <- transform_sample_counts(FINAL_abundances_Q3_cutoff, function (x) sqrt(x)) #sqrt transformation
juv_Q3_sqrt <- subset_samples(FINAL_abundances_Q3_sqrt, sampleType == "juvenile") #only juveniles
juv_Q3_adonis = as (sample_data(juv_Q3_sqrt), "data.frame")
juv_Q3_d = phyloseq::distance(juv_Q3_sqrt,'bray')
Adonis_Q3_juv <-adonis2(juv_Q3_d ~ climate1*climateFINAL + tank, data=juv_Q3_adonis,  permutations = 10000, method = "bray")
Adonis_Q3_juv

##dispersion ##
beta <- betadisper(juv_Q3_d, sample_data(juv_Q3_sqrt)$climate1)
disper.test = permutest(beta, permutations =10000)
disper.test # OK

beta <- betadisper(juv_Q3_d, sample_data(juv_Q3_sqrt)$climateFINAL)
disper.test = permutest(beta, permutations =10000)
disper.test # OK

## post hoc tests ##
sample_data(juv_Q3_sqrt)$climate1F <-interaction(sample_data(juv_Q3_sqrt)$climate1, sample_data(juv_Q3_sqrt)$climateFINAL)
testing_temp = pairwise.perm.manova(juv_Q3_d, sample_data(juv_Q3_sqrt)$climate1F,
                                    nperm=10000, p.method = "BH")
testing_temp$p.value 



############################################################
### DESEQ2 - differentially abundant ASVs  - JUVENILE F1 ###

phyloseq_merged_FINAL_filtered #count data
phyloseq_merged_FINAL_deseq_Q3 <- subset_samples(phyloseq_merged_FINAL_filtered, sample_generation != "larvae_F1_first" & sample_generation != "larvae_F1_second" & sample_generation != "gonad_F1") #filtering
phyloseq_merged_FINAL_deseq_Q3_noZeroZero <- prune_taxa(taxa_sums(phyloseq_merged_FINAL_deseq_Q3) > 0, phyloseq_merged_FINAL_deseq_Q3)
sample_data(phyloseq_merged_FINAL_deseq_Q3_noZeroZero)$treatment <- as.factor(sample_data(phyloseq_merged_FINAL_deseq_Q3_noZeroZero)$treatment) #converting treatment into factor
Filtered1 <- subset_samples(phyloseq_merged_FINAL_deseq_Q3_noZeroZero, sampleType == "juvenile")
filtering1 <- subset_samples(FINAL_abundances_Q3_cutoff,  sampleType == "juvenile")


### ambient-2050 VS 2050-2050 ###

Filtered <- subset_samples(Filtered1, treatment =="CCM" | treatment =="MMM")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q3_dds = phyloseq_to_deseq2(Filtered, ~treatment) 
Q3_dds$treatment<-relevel(Q3_dds$treatment,ref="MMM") 
Q3_dds = DESeq(Q3_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q3 <- results(Q3_dds, contrast=c("treatment","CCM", "MMM"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q3) #summary of results
res_Q3_order <- res_Q3[order(res_Q3$padj),] #order according to padj 
res_Q3_order = res_Q3_order[order(res_Q3_order$padj, na.last=NA), ]
res_Q3_order_sig = res_Q3_order[(res_Q3_order$padj < alpha), ] #only significant
res_Q3_order_sig #19 ASVs ->  among them, I keep only ASvs present in at least half of the samples in one treatment group
#bubbleplot 
filtering <- subset_samples(filtering1, treatment =="CCM" | treatment =="MMM")
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("8e213ca9798cc21f92f8566a165ef1aa", "4e84e68433f5335a1bb057f24a0687ee", "7e903dc05dd350ab74ffe7f7523c5abc",
                                                                              "408c648b8883625ecc3e098ff5f53698")) #differentially abundant ASVs

variable1 = as.character(get_variable(physeq.subset, "sample"))
variable2 = as.character(get_variable(physeq.subset, "treatment"))
sample_data(physeq.subset)$rep_treat <- mapply(paste0, variable2, variable1,
                                               collapse = "_") #create new vairable in phyloseq
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Family,  Genus, OTU, sep=";", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = rep_treat, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,25), name="Relative abundance (%)" ) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_treatment, scales = "free_x", space = "free_x")
bubbleplot +
  scale_color_manual(values = c( "D_2__Bacteroidia" = "#C3E7BC", "D_2__Verrucomicrobiae" = "#F3B962", "D_2__Gammaproteobacteria" = "#C0BCD7")) 


### ambient-2100 VS 2050-2100 ###

Filtered <- subset_samples(Filtered1, treatment =="CCH" | treatment =="MMH")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
head(sample_data(Filtered)$treatment, 20) 
Q3_dds = phyloseq_to_deseq2(Filtered, ~treatment) 
Q3_dds$treatment<-relevel(Q3_dds$treatment,ref="CCH") 
Q3_dds = DESeq(Q3_dds, test="Wald", fitType="parametric", sfType="poscounts") 
res_Q3 <- results(Q3_dds, contrast=c("treatment","MMH", "CCH"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q3)
res_Q3_order <- res_Q3[order(res_Q3$padj),] #order according to padj 
res_Q3_order = res_Q3_order[order(res_Q3_order$padj, na.last=NA), ]
res_Q3_order_sig = res_Q3_order[(res_Q3_order$padj < alpha), ] #only significant
res_Q3_order_sig #11 ASVs ->  among them, I keep only ASvs present in at least half of the samples in one treatment group
KEEP <- c("549482182dc3ab964cef7ea4ca4f82a6", "c604642bc71763fb6886501c2e0a9f77", "0bb16706a59a3944ef14883539a8e8af",
                                                                              "d5ea657a986a9ec747004ac42316b7fb")) #differentially abundant ASVs


### ambient-2100 VS 2100-2100 ###

Filtered <- subset_samples(Filtered1, treatment =="CCH" | treatment =="HHH")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
table(sample_data(Filtered)$sampleType)
head(sample_data(Filtered)$treatment, 21) 
Q3_dds = phyloseq_to_deseq2(Filtered, ~treatment) 
Q3_dds$treatment<-relevel(Q3_dds$treatment,ref="CCH") 
Q3_dds = DESeq(Q3_dds, test="Wald", fitType="parametric", sfType="poscounts") 
alpha = 0.01
res_Q3 <- results(Q3_dds, contrast=c("treatment","HHH", "CCH"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
summary (res_Q3) #summary of results
res_Q3_order <- res_Q3[order(res_Q3$padj),] #order according to padj 
res_Q3_order = res_Q3_order[order(res_Q3_order$padj, na.last=NA), ]
res_Q3_order_sig = res_Q3_order[(res_Q3_order$padj < alpha), ] #only significant
res_Q3_order_sig #16 ASVs ->  among them, I keep only ASvs present in at least half of the samples in one treatment group:
KEEP <- c("549482182dc3ab964cef7ea4ca4f82a6", "696e2d293f8815faa541ed20318549fe", "ce4efbce621af5f2670efae86439c6ba",
                                                                              "0bb16706a59a3944ef14883539a8e8af")) #differenttially abundant ASVs

  
## JUVENILE F1 AMBIENT-2100 VS 2050-2100 VS 2100-2100 BUBBLE PLOT ##	

filtering <- subset_samples(filtering1,  sampleType == "juvenile" & treatment =="CCH" | treatment =="HHH" | treatment =="MMH")
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("549482182dc3ab964cef7ea4ca4f82a6", "c604642bc71763fb6886501c2e0a9f77", "0bb16706a59a3944ef14883539a8e8af",
                                                                              "d5ea657a986a9ec747004ac42316b7fb", "696e2d293f8815faa541ed20318549fe", "ce4efbce621af5f2670efae86439c6ba"))

variable1 = as.character(get_variable(physeq.subset, "sample"))
variable2 = as.character(get_variable(physeq.subset, "treatment"))
sample_data(physeq.subset)$rep_treat <- mapply(paste0, variable2, variable1,
                                               collapse = "_") #create new vairable in phyloseq
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Family,  Genus, OTU, sep=";", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
physeq.subset_m$sample_treatment_g = factor (physeq.subset_m$sample_treatment, levels=c('juvenile_CCH', 'juvenile_MMH', 'juvenile_HHH'))
bubbleplot <- ggplot(physeq.subset_m, aes(x = rep_treat, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,25), name="Relative abundance (%)" ) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),#to eliminate gridlines
    panel.grid.minor = element_blank()) +
  facet_grid(~ sample_treatment_g, scales = "free_x", space = "free_x")
bubbleplot +
  scale_color_manual(values = c("D_2__Gammaproteobacteria" = "#C0BCD7", "D_2__Deltaproteobacteria" = "#DFDEC5", "D_2__Oxyphotobacteria" = "#89ABC9")) 




############################################################
### ALPHA DIVERSITY ###

## pre-processing ##
juve_rarefied_Q3 <- subset_samples(rarefied_Q3_NOcutoff, sampleType == "juvenile")
richQ3 = estimate_richness(juve_rarefied_Q3, split=TRUE) 
richQ3<-as.data.frame(richQ3)
RICHNESS_Q3<-tibble::rownames_to_column(as.data.frame(richQ3), var="sampleID")
RICHNESS_Q3$sampleID = gsub("X", "", RICHNESS_Q3$sampleID) 
metadata14 <- tibble::rownames_to_column(as.data.frame(metadata8), var="sampleID")
richness_table_Q3 <- inner_join(RICHNESS_Q3, metadata14, by="sampleID")
richness_table_Q3$climate1 <- as.factor(richness_table_Q3$climate1)
richness_table_Q3$climateFINAL <- as.factor(richness_table_Q3$climateFINAL)
richness_table_Q3 %>% dplyr::group_by(climate1, climateFINAL)%>% count()
richness_table_Q3 %>% pull(climate1) %>% levels()
richness_table_Q3 %>% pull(climateFINAL) %>% levels()
richness_Q3_juv <- richness_table_Q3 %>% mutate(Group =  factor(paste(climate1,climateFINAL)))
richness_Q3_juv %>% pull(Group) %>% levels()

## stats Shannon diversity
mod1 <- glmmTMB(Shannon ~ Group + (1|tank), data = richness_Q3_juv, family = 'gaussian', REML = TRUE)
summary(mod1)
simulateResiduals(fittedModel = mod1, plot = T) #no problems detected
emmeans (mod1, ~Group, type = "response") %>% pairs %>%rbind(adjust='bh')


## plot ##
alpha_plot_Q3_order <- richness_table_Q3 %>% mutate(name2 = fct_relevel(sample_treatment, "juvenile_CCC", "juvenile_CCM", "juvenile_CCH", "juvenile_MMM", "juvenile_MMH", "juvenile_HHH")) #to reorder x axis
p <- ggplot(alpha_plot_Q3_order, aes(x=name2, y=Shannon)) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + theme_bw() + labs(title = "Shannon") + theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=climate1) +
  geom_point(colour = "black", size = 1)+ ylim(2.1, 6.7) + labs(title = "Shannon transplant") + theme_classic(base_size = 15)
p + scale_x_discrete(breaks=c("juvenile_CCC", "juvenile_CCM", "juvenile_CCH", "juvenile_MMM", "juvenile_MMH", "juvenile_HHH"),
                     labels=c("Ambient - Ambient", "Ambient - 2050", "Ambient - 2100", "2050 - 2050", "2050 - 2100", "2100 - 2100")) + theme(axis.text.x=element_text(angle = 90, hjust = 0))

