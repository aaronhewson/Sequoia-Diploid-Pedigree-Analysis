# Pedigree analysis using Sequoia for Jim Dunckley Heritage Orchard and Plant & Food Research samples SNP genotyped in 2025 for Aaron Hewson's MSc project.
# The majority of samples were identified as known cultivars by duplicate matching, this is an attempt to identify parents of remaining samples.

# Prior to running this script, all samples and reference accessions were combined into one PLINK .ped/.map file pair, with ACTG allele codes, 
# using a list from Nick Howard's to identify the high-quality SNPs common across the apple 20K, 50K, and 480K SNP array platforms.


# Load Packages -----------------------------------------------------------

#install.packages("sequoia")
#install.packages("igraph")
#install.packages("tibble")

library(sequoia)
library(igraph)
library(tibble)


# Set Working Directory ---------------------------------------------------

setwd("C:/Users/curly/Desktop/Apple Genotyping/Methods/Sequoia Diploid Pedigree Analysis")


# Load Input Data for Sequoia ---------------------------------------------
# All diacritical marks and apostrophes were removed from sample names before loading!

#Load genotypes (PLINK .raw file)
geno <- GenoConvert( InFile = "Inputs/All_Samples.raw", InFormat = "raw", OutFormat = "seq", sep = "\t", Missing = "NA")

#Load life history (Sex and BirthYear given as "NA" = unknown)
LH <- read.table("Inputs/Life_History.txt", header = F, sep = "\t")

#Load pedigree (Mum and Dad given as "NA" = unknown)
pedi <- read.table("Inputs/Pedigree.txt", header=F, sep = "\t")


# Check samples for duplicates --------------------------------------------

DupCheck <- sequoia(GenoM = geno, Err = 0.025, LifeHistData = LH, Module = "dup", quiet = FALSE, Plot = TRUE)

Dup <- DupCheck$DupGenotype

#With a high error rate, Sequoia can assign erroneous duplicates - evaluate in the .csv file using "Mismatch" and "LLR".
write.csv(Dup, "Inputs/SeqDup.csv")

#Duplicate pairs involving "HW_21_MN1627" are not actual duplicates. These erroneous duplicates need to be removed later from the list of duplicate samples.


# Group duplicates -------------------------------------------------------

#Load in duplicate pair names
DupList <- read.csv("Inputs/SeqDup.csv", header = TRUE)
DupList <- DupList[!colnames(DupList) %in% c("X", "row1", "row2", "Mismatch", "SnpdBoth", "LLR")]

##Combine duplicate pairs into duplicate groups
#Group duplicates with igraph
graph <- graph_from_data_frame(DupList, directed = FALSE)
components <- components(graph)

#Sort groups by # of duplicates
group_sizes <- table(components$membership)
sorted_group_ids <- order(group_sizes)
new_ids <- match(components$membership, sorted_group_ids)
V(graph)$group <- new_ids
grouped_samples <- split(names(components$membership), new_ids)

#Pad group with length less than max length with NA's
max_len <- max(sapply(grouped_samples, length))
padded_list <- lapply(grouped_samples, function(x) {c(x, rep(" ", max_len - length(x)))})

#Add groups to dataframe, assign group number, add count of duplicates in each group
dd <- as.data.frame(do.call(rbind, padded_list))
dd <- cbind(Group = seq_len(nrow(dd)), dd)
sample_counts <- rowSums(dd[, -1] != " ")
dd <- add_column(dd, SampleCount = sample_counts, .after = "Group")

##Save .csv of duplicate groupings
write.csv(dd, "Inputs/Grouped_Duplicates.csv", row.names = FALSE)


# Remove Duplicates -------------------------------------------------------
#From the Grouped Duplicates, determine one sample of each duplicate genotype to keep. Create a .txt file that lists ALL other duplicate sample names to be removed.
#Remove false duplicates from this list, so they are not wrongly removed. Names to remove: "7_10_T145R", "Red_Baron", "DN_2302", "HW_21_MN1627"
#Triploids also need to be removed from the dataset.
#Triploids to remove from this dataset: DN_482 , DN_630 , DN_869 , DN_901 , DN_907 , DN_1023 , DN_2219 , DN_2231 , DN_2334 , DN_2338 , DN_2458 , DN_5234

RemoveDup <- readLines("Inputs/RemoveDup.txt")

geno_Clean <- geno[!(row.names(geno) %in% c(RemoveDup,"DN_482", "DN_630", "DN_869", "DN_901", "DN_907", "DN_1023", "DN_2219",
                                          "DN_2231", "DN_2334", "DN_2338", "DN_2458", "DN_5234")),]

pedi_Clean <- pedi[!(pedi$V1 %in% c(RemoveDup,"DN_482", "DN_630", "DN_869", "DN_901", "DN_907", "DN_1023", "DN_2219",
                                          "DN_2231", "DN_2334", "DN_2338", "DN_2458", "DN_5234", "7_10_T145R", 
                                          "Red_Baron", "DN_2302", "HW_21_MN1627")),]

LH_Clean <- LH[!(LH$V1 %in% c(RemoveDup,"DN_482", "DN_630", "DN_869", "DN_901", "DN_907", "DN_1023", "DN_2219",
                                          "DN_2231", "DN_2334", "DN_2338", "DN_2458", "DN_5234", "7_10_T145R", 
                                          "Red_Baron", "DN_2302", "HW_21_MN1627")),]

#Recheck all duplicates are removed
CleanCheck <- sequoia(GenoM = geno_Clean, Err = 0.025, LifeHistData = LH_Clean, Module = "dup", quiet = FALSE)

CleanCheck$DupGenotype


# Identify potential parent-offspring relationships -----------------------

getMay <- GetMaybeRel (GenoM = geno_Clean,  LifeHistData = LH_Clean, Err = 0.025, Module = "par" , quiet = FALSE, Herm="no", MaxPairs = 100000000)

#List of potential parent-offspring trios 

trios <-getMay$MaybeTrio

write.csv(trios, "Inputs/maybetrios.csv", row.names = FALSE)

#List of potential parent-offspring pairs (where one parent is missing)

duos <- getMay$MaybePar

write.csv(duos, "Inputs/maybeduos.csv", row.names = FALSE)
